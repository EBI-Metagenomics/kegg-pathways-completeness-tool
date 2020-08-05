#!/usr/bin/python3

import subprocess
import argparse
import sys


list_class = "all_pathways_class.txt"
list_name = "all_pathways_name.txt"
list_paths = "all_pathways.txt"
dop_file = "file.txt"


def check_brackets(line):
    positions1 = [1 for m in line if m == '(']
    positions2 = [1 for m in line if m == ')']
    if len(positions1) != len(positions2):
        return False
    else:
        return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generates Graphs for each contig")
    parser.add_argument("-m", "--mode", dest="mode", help="class, name or definition", required=True)
    parser.add_argument("-l", "--list", dest="list", help="list_pathways.txt", required=True)

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        mode = args.mode
        list_pathways = args.list

        if mode == 'name':
            with open(list_pathways, 'r') as file_in:
                for line in file_in:
                    line = line.strip()
                    print(line)
                    bashCommand = "echo {mo} >> {f}".format(mo=line, f=dop_file)
                    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()

                    bashCommand = "curl -s http://rest.kegg.jp/get/{mo} | " \
                              "grep '^NAME' | cut -c 13- >> {f}".format(mo=line, f=dop_file)
                    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()
            file_in.close()
            #  create final file
            with open(dop_file, 'r') as file_in, open(list_name, 'w') as file_name, open(list_pathways, 'r') as file_list:
                for line, mo in zip(file_in, file_list):
                    mo = mo.strip()
                    file_name.write(':'.join([mo, line]))
            file_name.close()
            file_list.close()
            file_in.close()

        elif mode == 'class':
            with open(list_pathways, 'r') as file_in:
                for line in file_in:
                    line = line.strip()
                    print(line)
                    bashCommand = "echo {mo} >> {f}".format(mo=line, f=dop_file)
                    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()

                    bashCommand = "curl -s http://rest.kegg.jp/get/{mo} | " \
                                  "grep '^CLASS' | cut -c 13- >> {f}".format(mo=line, f=dop_file)
                    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()
            file_in.close()
            # create final file
            with open(dop_file, 'r') as file_in, open(list_name, 'w') as file_class, open(list_pathways, 'r') as file_list:
                for line, mo in zip(file_in, file_list):
                    mo = mo.strip()
                    file_class.write(':'.join([mo, line]))
            file_class.close()
            file_list.close()
            file_in.close()

        elif mode == 'definition':
            list_MOs = {}
            with open(list_pathways, 'r') as file_in:
                for line in file_in:
                    line = line.strip()
                    list_MOs[line] = ""

                    print(line)
                    bashCommand = "echo {mo} >> {f}".format(mo=line, f=dop_file)
                    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()

                    bashCommand = "curl -s http://rest.kegg.jp/get/{mo} > {f}".format(mo=line, f=dop_file+'1')
                    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()

                    bashCommand = "less {f1} | grep -A 5 '^DEF' >> {f}".format(mo=line, f=dop_file, f1=dop_file+'1')
                    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()

            file_in.close()

            weird_mo = []
            collection = []
            with open(dop_file, 'r') as file_def:
                for line in file_def:
                    line = line.strip()
                    if line in list_MOs:
                        cur_mo = line
                        collection = []
                    elif 'DEFINITION' in line:
                        cur_def = line.split('DEFINITION')[1].lstrip().rstrip()
                        collection.append('(' + cur_def + ')')
                    elif 'ORTHOLOGY' in line:
                        if len(collection) > 1:
                            weird_mo.append(cur_mo)
                            print(cur_mo, len(collection))
                            cur_line = ','.join(collection).rstrip()
                            print(cur_line)
                            list_MOs[cur_mo] = cur_line
                        else:
                            list_MOs[cur_mo] = collection[0]
                    else:
                        cur_def = '(' + line.lstrip().rstrip() + ')'
                        collection.append(cur_def)
            print(len(weird_mo))

            with open(list_paths, 'w') as list_path_file:
                for mo in list_MOs:
                    list_path_file.write(mo + ':' + list_MOs[mo] + '\n')
