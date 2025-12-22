#!/usr/bin/env python3
# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import re
import time
import click
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from importlib.metadata import version, PackageNotFoundError
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Semaphore
from tqdm import tqdm


LIST_MODULES = "modules_list.txt"
MODULES_TABLE = "modules_table.tsv"
CHANGED_TABLE = "changed.tsv"
LIST_SEPARATED_IN_DEFINITION = "mediators_list.txt"

KEGG_API_LIST_MODULES = "http://rest.kegg.jp/list/module"
KEGG_API_GET = "http://rest.kegg.jp/get"


def get_version():
    """Get package version from installed metadata"""
    try:
        return version("kegg-pathways-completeness")
    except PackageNotFoundError:
        return "unknown"


class ModulesDataFetchTool:
    def __init__(
        self,
        outdir: str,
        max_workers: int = 10,
        max_retries: int = 3,
        delay_between_requests: float = 0.5,
    ):
        """
        This script uses requests to get list of KEGG modules via KEGG API with parallel fetching.
        It generates a tab-separated table with columns: module, definition, name, class.
        It is useful to run whole fetch and update because KEGG changes modules structure for already existing modules.
        :param outdir: path to output directory
        :param max_workers: maximum number of concurrent requests (default: 10)
        :param max_retries: maximum number of retries for failed requests (default: 3)
        :param delay_between_requests: delay in seconds between requests (default: 0.5)
        :return:
        """
        self.outdir = outdir
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        self.modules_data = []
        self.max_workers = max_workers
        self.max_retries = max_retries
        self.delay_between_requests = delay_between_requests
        self.session = self._create_session_with_retries()
        # Rate limiting semaphore
        self.rate_limiter = Semaphore(max_workers)
        # Track failed requests
        self.failed_modules = []

    def _create_session_with_retries(self):
        """Create a requests session with automatic retry logic"""
        session = requests.Session()

        # Configure retry strategy
        retry_strategy = Retry(
            total=self.max_retries,
            backoff_factor=2,  # Wait 2s, 4s, 8s between retries
            status_forcelist=[
                403,
                429,
                500,
                502,
                503,
                504,
            ],  # Retry on these HTTP status codes
            allowed_methods=["GET"],  # Only retry GET requests
            raise_on_status=False,  # Don't raise exception on retry
        )

        # Mount the retry adapter to both http and https
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)

        return session

    def fetch_and_save_kegg_modules_list(self):
        # Ensure the output directory exists
        os.makedirs(self.outdir, exist_ok=True)
        output_file = os.path.join(self.outdir, LIST_MODULES)

        try:
            print(f"Fetching modules list from {KEGG_API_LIST_MODULES}")
            # Fetch the data from the KEGG API using requests
            response = self.session.get(KEGG_API_LIST_MODULES, timeout=30)
            response.raise_for_status()

            # Process the output to extract the module codes
            module_codes = [line.split("\t")[0] for line in response.text.splitlines()]

            # Write the results to the output file
            with open(output_file, "w") as f:
                for code in module_codes:
                    f.write(code + "\n")

            print(
                f"Successfully saved {len(module_codes)} KEGG module codes to {output_file}"
            )
            return module_codes
        except requests.RequestException as e:
            print(f"An error occurred while fetching modules list: {e}")
            exit(1)

    def fetch_module_info(self, module, output_dir):
        """Fetch information for a single module"""
        with self.rate_limiter:
            try:
                # Add delay between requests to avoid rate limiting
                time.sleep(self.delay_between_requests)

                line_separator_in_kegg_names = " "  # space=AND, comma=OR
                url = f"{KEGG_API_GET}/{module}"

                # Fetch the data from the KEGG API using requests
                response = self.session.get(url, timeout=30)
                response.raise_for_status()

                lines = response.text.split("\n")
                definition_lines = []
                capture_definition = False
                name = ""
                class_module = ""

                for line in lines:
                    if line.startswith("NAME"):
                        match = re.search(r"^NAME\s+(.*)", line)
                        if match:
                            name = match.group(1)
                    elif line.startswith("CLASS"):
                        match = re.search(r"^CLASS\s+(.*)", line)
                        if match:
                            class_module = match.group(1)
                    elif line.startswith("DEFINITION"):
                        # can be in multiple lines
                        match = re.search(r"^DEFINITION\s+(.*)", line)
                        if match:
                            definition_lines.append(match.group(1))
                        capture_definition = True
                    elif capture_definition and line.startswith(" "):
                        definition_lines.append(line.strip())
                    else:
                        capture_definition = False

                if len(definition_lines) > 1:
                    for i in range(len(definition_lines)):
                        if len(definition_lines[i]) > 6:
                            # not a single KOXXXX
                            definition_lines[i] = f"({definition_lines[i]})"
                    with open(os.path.join(output_dir, LIST_SEPARATED_IN_DEFINITION), 'a') as file_out:
                        file_out.write(f'{module}:{line_separator_in_kegg_names.join(definition_lines)}' + '\n')

                def_module = line_separator_in_kegg_names.join(definition_lines)

                return {
                    "module": module,
                    "definition": def_module,
                    "name": name,
                    "class": class_module,
                }

            except requests.RequestException as e:
                self.failed_modules.append(module)
                print(
                    f"\nError fetching module {module} after {self.max_retries} retries: {e}"
                )
                return None

    def fetch_all_modules_parallel(self, modules, outdir):
        """Fetch all modules in parallel with progress bar"""
        print(
            f"Fetching {len(modules)} modules with {self.max_workers} concurrent workers..."
        )

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_module = {
                executor.submit(self.fetch_module_info, module, outdir): module
                for module in modules
            }

            # Process results as they complete with progress bar
            with tqdm(
                total=len(modules), desc="Fetching modules", unit="module"
            ) as pbar:
                for future in as_completed(future_to_module):
                    result = future.result()
                    if result:
                        self.modules_data.append(result)
                    pbar.update(1)

        print(f"Successfully fetched {len(self.modules_data)} modules")

    def write_modules_table(self):
        """Write modules data to a tab-separated table"""
        output_file = os.path.join(self.outdir, MODULES_TABLE)
        try:
            with open(output_file, "w") as f:
                # Write header
                f.write("module\tdefinition\tname\tclass\n")
                # Write data rows
                for module_data in self.modules_data:
                    f.write(
                        f"{module_data['module']}\t{module_data['definition']}\t{module_data['name']}\t{module_data['class']}\n"
                    )
            print(f"Successfully saved modules table to {output_file}")
        except Exception as e:
            print(f"An error occurred while writing table: {e}")
            exit(1)

    @staticmethod
    def parse_old_file(file_path):
        """
        Parse old format file (module:value) into a dictionary
        :param file_path: path to the old file
        :return: dict with module as key and value as value
        """
        old_data = {}
        try:
            with open(file_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if line and ":" in line:
                        parts = line.split(":", 1)
                        if len(parts) == 2:
                            module, value = parts
                            old_data[module] = value
            return old_data
        except FileNotFoundError:
            print(f"Warning: File {file_path} not found")
            return {}
        except Exception as e:
            print(f"Error parsing file {file_path}: {e}")
            return {}

    def detect_changes(self, old_definitions=None, old_names=None, old_classes=None):
        """
        Detect changes between old and new module data
        :param old_definitions: dict of old definitions {module: definition}
        :param old_names: dict of old names {module: name}
        :param old_classes: dict of old classes {module: class}
        :return: list of changed module data
        """
        changed_modules = []

        for module_data in self.modules_data:
            module = module_data["module"]
            is_changed = False

            # Check if definition changed
            if old_definitions is not None:
                if module in old_definitions:
                    if old_definitions[module] != module_data["definition"]:
                        is_changed = True
                else:
                    # New module not in old data
                    is_changed = True

            # Check if name changed
            if old_names is not None:
                if module in old_names:
                    if old_names[module] != module_data["name"]:
                        is_changed = True
                else:
                    # New module not in old data
                    is_changed = True

            # Check if class changed
            if old_classes is not None:
                if module in old_classes:
                    if old_classes[module] != module_data["class"]:
                        is_changed = True
                else:
                    # New module not in old data
                    is_changed = True

            if is_changed:
                changed_modules.append(module_data)

        return changed_modules

    def write_changed_table(self, changed_modules):
        """Write changed modules data to a tab-separated table"""
        output_file = os.path.join(self.outdir, CHANGED_TABLE)
        try:
            with open(output_file, "w") as f:
                # Write header
                f.write("module\tdefinition\tname\tclass\n")
                # Write data rows
                for module_data in changed_modules:
                    f.write(
                        f"{module_data['module']}\t{module_data['definition']}\t{module_data['name']}\t{module_data['class']}\n"
                    )
            print(f"Successfully saved changed modules table to {output_file}")
            print(f"Total changed modules: {len(changed_modules)}")
        except Exception as e:
            print(f"An error occurred while writing changed table: {e}")
            exit(1)


@click.command()
@click.option(
    "-o",
    "--output-dir",
    default="fetched_data",
    help="Output directory",
    show_default=True,
)
@click.option(
    "--list-modules",
    type=click.Path(exists=True),
    help="Path to existing list_modules.txt file (if not provided, will fetch from API)",
)
@click.option(
    "--old-definitions",
    type=click.Path(exists=True),
    help="Path to old definitions file (format: module:definition)",
)
@click.option(
    "--old-names",
    type=click.Path(exists=True),
    help="Path to old names file (format: module:name)",
)
@click.option(
    "--old-classes",
    type=click.Path(exists=True),
    help="Path to old classes file (format: module:class)",
)
@click.option(
    "--max-workers",
    default=5,
    type=int,
    help="Maximum number of concurrent API requests",
    show_default=True,
)
@click.option(
    "--max-retries",
    default=3,
    type=int,
    help="Maximum number of retries for failed requests",
    show_default=True,
)
@click.option(
    "--delay",
    default=0.5,
    type=float,
    help="Delay in seconds between requests to avoid rate limiting",
    show_default=True,
)
@click.version_option(version=get_version(), prog_name="fetch_modules_data")
def main(
    output_dir,
    list_modules,
    old_definitions,
    old_names,
    old_classes,
    max_workers,
    max_retries,
    delay,
):
    """
    Script fetches KEGG API for list of modules with NAME, DEFINITION and CLASS.

    Generates a tab-separated table (modules_table.tsv) with all current module data.
    If old data files are provided, also generates a changed.tsv file with modules
    that have been modified.

    Uses parallel requests to speed up data fetching significantly.
    """
    module_fetch_tool = ModulesDataFetchTool(
        outdir=output_dir,
        max_workers=max_workers,
        max_retries=max_retries,
        delay_between_requests=delay,
    )

    # get a list of all modules - either from file or API
    if list_modules:
        print(f"Loading module list from {list_modules}")
        try:
            with open(list_modules, "r") as f:
                modules_all = [line.strip() for line in f if line.strip()]
            print(f"Loaded {len(modules_all)} modules from file")
        except Exception as e:
            print(f"Error loading module list from file: {e}")
            exit(1)
    else:
        modules_all = module_fetch_tool.fetch_and_save_kegg_modules_list()

    # get info for each module in parallel
    module_fetch_tool.fetch_all_modules_parallel(modules_all, output_dir)

    # write the tab-separated table
    module_fetch_tool.write_modules_table()

    # Check if any old files were provided for change detection
    old_definitions_data = None
    old_names_data = None
    old_classes_data = None

    if old_definitions:
        print(f"Loading old definitions from {old_definitions}")
        old_definitions_data = module_fetch_tool.parse_old_file(old_definitions)

    if old_names:
        print(f"Loading old names from {old_names}")
        old_names_data = module_fetch_tool.parse_old_file(old_names)

    if old_classes:
        print(f"Loading old classes from {old_classes}")
        old_classes_data = module_fetch_tool.parse_old_file(old_classes)

    # If any old data was provided, detect changes
    if (
        old_definitions_data is not None
        or old_names_data is not None
        or old_classes_data is not None
    ):
        print("Detecting changes...")
        changed_modules = module_fetch_tool.detect_changes(
            old_definitions=old_definitions_data,
            old_names=old_names_data,
            old_classes=old_classes_data,
        )
        if changed_modules:
            module_fetch_tool.write_changed_table(changed_modules)
        else:
            print("No changes detected.")


if __name__ == "__main__":
    main()
