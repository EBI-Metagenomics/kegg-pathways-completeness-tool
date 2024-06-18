#!/usr/bin/python3

import subprocess

dop_file = 'check_list.txt'
count = 1
list_all_modules = []
with open('../../pathways/list_pathways.txt', 'r') as list_modules:
    for mo in list_modules:
        number = int(mo.strip().split('M')[1])
        if count != number:
            if count < 10:
                zeros = '0000'
            elif count < 100:
                zeros = '000'
            else:
                zeros = '00'
            name = 'M' + zeros + str(count)
            list_all_modules.append(name)
            print(name)
            count = number
        else:
            list_all_modules.append(mo.strip())
        count += 1

print(len(list_all_modules))
for mo in list_all_modules:
    print(mo)
    bashCommand = "echo {mo} >> {f}".format(mo=mo, f=dop_file)
    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()

    bashCommand = "curl -s http://rest.kegg.jp/get/{mo} | grep '^NAME' | cut -c 13- >> {f}".format(mo=mo, f=dop_file)
    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True).wait()

prev = ''
with open(dop_file, 'r') as check_file:
    for line in check_file:
        line = line.strip()
        if len(line.split('M0')) > 1 and len(prev.split('M0')) > 1:
            print(prev)
        prev = line