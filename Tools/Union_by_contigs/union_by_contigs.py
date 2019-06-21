import argparse

def union_KOs(filename):
    dict_by_contigs = {}
    with open(filename, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            name = line[0].split('.')[0]
            orthologs = line[1:]
            if name not in dict_by_contigs:
                dict_by_contigs[name] = []
            dict_by_contigs[name] += orthologs
    return dict_by_contigs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unites KO relateted to one contig")
    parser.add_argument("-i", "--input", dest="input_file", help="Each line = protein", required=True)
    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="Relative path to directory where you want the output file to be stored (default: cwd)",
                        default=".")
    args = parser.parse_args()
    dict_by_contigs = union_KOs(args.input_file)
    output_name_list = args.input_file.split('/')
    output_name = output_name_list[len(output_name_list)-1].split('.')[0]

    with open('union_ko_contigs.txt', 'w') as file_out:
        for contig in dict_by_contigs:
            line_out = '\t'.join([contig]+dict_by_contigs[contig]) + '\n'
            file_out.write(line_out)
