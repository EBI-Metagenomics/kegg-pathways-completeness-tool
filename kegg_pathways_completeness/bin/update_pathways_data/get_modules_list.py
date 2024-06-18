import subprocess
import os

def fetch_and_save_kegg_modules():
    # Ensure the output directory exists
    output_dir = 'pathways'
    os.makedirs(output_dir, exist_ok=True)
    temp_file = os.path.join(output_dir, 'temp_list_module.txt')
    output_file = os.path.join(output_dir, 'list_pathways.txt')

    try:
        print(f'Fetching modules from http://rest.kegg.jp/list/module')
        # Fetch the data from the KEGG API using wget
        wget_command = f'wget -qO- http://rest.kegg.jp/list/module'
        wget_result = subprocess.run(wget_command, shell=True, capture_output=True, text=True)
        wget_output = wget_result.stdout

        # Process the output to extract the module codes
        module_codes = [line.split('\t')[0]for line in wget_output.splitlines()]
        # Write the results to the output file
        with open(output_file, 'w') as f:
            for code in module_codes:
                f.write(code + '\n')

        print(f'Successfully saved KEGG module codes to {output_file}')
        return module_codes
    except subprocess.CalledProcessError as e:
        print(f'An error occurred: {e}')
        exit(1)


def compare_with_existing(new_modules):
    existing = []
    with open('kegg_pathways_completeness/pathways_data/all_pathways.txt', 'r') as file_in:
        for line in file_in:
            existing.append(line.strip().split(':')[0])
    result = list(set(new_modules).difference(set(existing)))
    if result:
        print(f'Found {len(result)} new modules')
        with open('new_modules.txt', 'w') as file_out:
            file_out.write('\n'.join(result))
    else:
        print('No new modules')

if __name__ == '__main__':
    new_modules = fetch_and_save_kegg_modules()
    compare_with_existing(new_modules)
