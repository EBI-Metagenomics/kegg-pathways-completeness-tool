# Update KEGG Modules Database

This guide explains how to fetch the latest KEGG modules data, generate graphs, and create visualizations.

## Overview

The workflow consists of three main steps:
1. **Fetch modules data** from KEGG API → generates `modules_table.tsv`
2. **Generate graphs** from modules data → generates `graphs.pkl`
3. **Create visualizations** (optional) → generates PNG plots

Additional notes: 
- Troubleshooting
- Migration from Old Format

## Step 1: Fetch KEGG Modules Data

### Quick Start

```bash
# Fetch all modules data (recommended settings)
fetch_modules_data \
  -o fetched_data \
  --max-workers 5 \
  --delay 0.5
```

This creates `fetched_data/modules_table.tsv` with all current KEGG modules information (ID, name, class, definition).  Typical runtime ~3-5 minutes for ~500 modules with default settings

### Command Help

```
Usage: fetch_modules_data [OPTIONS]

  Script fetches KEGG API for list of modules with NAME, DEFINITION and
  CLASS.

  Generates a tab-separated table (modules_table.tsv) with all current
  module data. If old data files are provided, also generates a changed.tsv
  file with modules that have been modified.

  Uses parallel requests to speed up data fetching significantly.

Options:
  -o, --output-dir TEXT           Output directory  [default: fetched_data]
  --list-modules PATH             Path to existing list_modules.txt file (if
                                  not provided, will fetch from API)
  --old-definitions PATH          Path to old definitions file (format:
                                  module:definition)
  --old-names PATH                Path to old names file (format:
                                  module:name)
  --old-classes PATH              Path to old classes file (format:
                                  module:class)
  --max-workers INTEGER           Maximum number of concurrent API requests
                                  [default: 5]
  --max-retries INTEGER           Maximum number of retries for failed
                                  requests  [default: 3]
  --delay FLOAT                   Delay in seconds between requests to avoid
                                  rate limiting  [default: 0.5]
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```

### Advanced Usage

#### Detect Changes from Previous Version

Compare new data with old format files to identify what changed:

```bash
fetch_modules_data \
  -o fetched_data \
  --old-definitions pathways_data/all_pathways.txt \
  --old-names pathways_data/all_pathways_names.txt \
  --old-classes pathways_data/all_pathways_class.txt \
  --max-workers 5 \
  --delay 0.5
```

This generates both `modules_table.tsv` and `changed.tsv` (with only modified modules).

#### Resume from Existing List

If you already have a modules list, skip fetching it:

```bash
fetch_modules_data \
  -o fetched_data \
  --list-modules pathways_data/modules_list.txt \
  --max-workers 5 \
  --delay 0.5
```

#### Adjust for Network Conditions

For unstable connections, use more conservative settings:

```bash
# Conservative (safer, slower)
fetch_modules_data \
  -o fetched_data \
  --max-workers 3 \
  --delay 1.0 \
  --max-retries 5

# Aggressive (faster, may hit rate limits)
fetch_modules_data \
  -o fetched_data \
  --max-workers 10 \
  --delay 0.2 \
  --max-retries 3
```

### Output Files

- **`modules_table.tsv`** - Complete modules data with columns:
  - `module` - Module accession (e.g., M00001)
  - `definition` - KO definition/pathway
  - `name` - Module name
  - `class` - Module classification

- **`modules_list.txt`** - Simple list of module IDs (one per line)

- **`changed.tsv`** (optional) - Modules that changed since previous version


## Step 2: Generate Graphs

### Quick Start

```bash
# Generate graphs from TSV format
make_graphs \
  -i fetched_data/modules_table.tsv \
  -o graphs_output
```

This creates `graphs_output/graphs.pkl` containing networkx graph structures for all modules.

### Command Help

```
Usage: make_graphs [OPTIONS]

  Generates graph structures for KEGG modules and saves them to graphs.pkl.

  Supports both formats:
  - New TSV format: modules_table.tsv with columns (module, definition, name, class)
  - Old format: one module per line as module:definition

  The script automatically detects the input format.

Options:
  -i, --input PATH      Input file: TSV format (modules_table.tsv) or old
                        format (module:definition per line)  [required]
  -o, --outdir TEXT     Output directory where graphs.pkl will be stored
                        [default: outdir]
  -v, --verbose         Enable verbose logging
  --version             Show the version and exit.
  --help                Show this message and exit.
```

### Usage Examples

#### With New TSV Format (Recommended)

```bash
make_graphs \
  -i fetched_data/modules_table.tsv \
  -o graphs_output \
  -v
```

#### With Old Format (Backward Compatible)

```bash
make_graphs \
  -i pathways_data/all_pathways.txt \
  -o graphs_output
```

### Output

- **`graphs.pkl`** - Pickled dictionary containing:
  - Graph structure (networkx MultiDiGraph)
  - Edge definitions
  - Unnecessary nodes list

## Step 3: Generate Visualizations (Optional)

Create visual plots showing module pathway structures.

### Command Help

```
Usage: plot_modules_graphs [OPTIONS]

  Generate plots for KEGG module pathways.

  Can visualize modules with or without completeness information. Requires
  at least one of: --input-completeness, --modules, or --modules-file.

  Examples:
  # Plot with completeness data
  plot_modules_graphs -i completeness.tsv -g graphs.pkl -d definitions.txt

  # Plot specific modules
  plot_modules_graphs -m M00001 -m M00002 -g graphs.pkl -d definitions.txt

  # Plot modules from file
  plot_modules_graphs -l modules.txt -g graphs.pkl -d definitions.txt

Options:
  -i, --input-completeness PATH   Output table from give_completeness.py
  -m, --modules TEXT              Module accessions (can be specified
                                  multiple times)
  -l, --modules-file PATH         File containing module accessions
  -s, --file-separator TEXT       Modules separator in file  [default:
                                  ]
  -g, --graphs PATH               Graphs in pickle format  [default:
                                  pathways_data/graphs.pkl]
  -d, --definitions PATH          Pathways definitions file  [default:
                                  pathways_data/all_pathways.txt]
  -o, --outdir TEXT               Output directory for plots  [default:
                                  pathways_plots]
  --use-pydot                     Use pydot instead of graphviz
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```

### Usage Examples

#### Plot Specific Modules

```bash
# Single or multiple modules
plot_modules_graphs \
  -m M00001 -m M00002 -m M00003 \
  -g graphs_output/graphs.pkl \
  -d fetched_data/modules_table.tsv \
  -o my_plots
```

#### Plot from Modules List File

```bash
# Create file with module IDs
echo -e "M00001\nM00002\nM00003" > modules_to_plot.txt

plot_modules_graphs \
  -l modules_to_plot.txt \
  -g graphs_output/graphs.pkl \
  -d fetched_data/modules_table.tsv \
  -o my_plots
```

#### Plot with Completeness Data

After running `give_completeness` (see main documentation):

```bash
plot_modules_graphs \
  -i completeness_results/summary.kegg_pathways.tsv \
  -g graphs_output/graphs.pkl \
  -d fetched_data/modules_table.tsv \
  -o completeness_plots
```

### Output Structure

```
pathways_plots/
├── dot/           # GraphViz .dot files
│   ├── M00001.dot
│   ├── M00002.dot
│   └── ...
└── png/           # PNG images
    ├── M00001.png
    ├── M00002.png
    └── ...
```


## Troubleshooting

### 403 Forbidden Errors

If you encounter `403 Forbidden` errors from KEGG API:

1. **Reduce concurrent workers**: `--max-workers 3`
2. **Increase delay**: `--delay 1.0`
3. **Add more retries**: `--max-retries 5`

## Migration from Old Format

If you have existing data in old format (separate files), the new scripts are fully backward compatible:

```bash
# make_graphs accepts old format
make_graphs \
  -i pathways_data/all_pathways.txt \
  -o graphs_output

# give_completeness can use old format
give_completeness \
  -i input.tsv \
  -a pathways_data/all_pathways.txt \
  -n pathways_data/all_pathways_names.txt \
  -c pathways_data/all_pathways_class.txt \
  -o results

# Or use new TSV format (recommended)
give_completeness \
  -i input.tsv \
  -t fetched_data/modules_table.tsv \
  -o results
```