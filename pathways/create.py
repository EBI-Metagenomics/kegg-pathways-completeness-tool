import subprocess


list_pathways = "list_pathways.txt"
list_class = "all_pathways_class.txt"
list_name = "all_pathways_name.txt"
dop_file = "file.txt"
with open(list_pathways,'r') as file_in:
    for line in file_in:
        line = line.strip()
        bashCommand = "echo {mo} >> {f}".format(mo=line, f=dop_file)
        subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True)

        bashCommand = "curl -s http://rest.kegg.jp/get/{mo} | " \
                  "grep '^NAME' | cut -c 13- >> {f}".format(mo=line, f=dop_file)
        subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True)
file_in.close()
"""
with open(dop_file, 'r') as file_in, open(list_name, 'w') as file_name, open(list_pathways,'r') as file_list:
    for line, mo in zip(file_in, file_list):
        mo = mo.strip()
        file_name.write(':'.join([mo, line]))
"""