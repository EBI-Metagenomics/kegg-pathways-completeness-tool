import pandas as pd

def get_GK(filename_koala):
    # GoastKoala reading
    koala_dict = {}
    kol_empty = 0
    with open(filename_koala, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            if line[0] not in koala_dict:
                koala_dict[line[0]] = []
            if len(line) > 1:
                koala_dict[line[0]].append(line[1])
            else:
                kol_empty += 1
    # save
    basename_list = filename_koala.split('/')
    basename = basename_list[len(basename_list) - 1].split('.')[0]
    print(basename)
    with open(basename + '_koala.txt', 'w') as file_out:
        for key in koala_dict:
            file_out.write('\t'.join([key] + koala_dict[key]) + '\n')
    print('empty koala annotations', kol_empty)
    return koala_dict


def get_hmm(filename_hmm):
    # reading hmmscan
    hmm_dict = {}
    with open(filename_hmm, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            key = line[0][:-1]
            if key not in hmm_dict:
                hmm_dict[key] = []
            if len(line) > 1:
                hmm_dict[key] += line[1:]
    return hmm_dict


def comparison_all_KO(hmm_dict, koala_dict):
    compare_tbl = [['key', 'in_common', 'in_GK', 'in_hmm']]
    kol_intersections, kol_empty = [0, 0]
    for key in all_labels:
        new_line = [key]
        common, in_hmm, in_gk = ["" for _ in range(3)]
        if key in hmm_dict and key in koala_dict:  # common
            common = ','.join(list(set(hmm_dict[key]).intersection(set(koala_dict[key]))))
            in_hmm = ','.join(list(set(hmm_dict[key]).difference(set(common))))
            in_gk = ','.join(list(set(koala_dict[key]).difference(set(common))))
            if common != "":
                kol_intersections += 1
            else:
                kol_empty += 1
        elif key in hmm_dict:
            in_hmm = ','.join(hmm_dict[key])
        elif key in koala_dict:
            in_gk = ','.join(koala_dict[key])
        new_line += [common, in_gk, in_hmm]
        compare_tbl.append(new_line)

    table = pd.DataFrame(data=compare_tbl)
    table.to_csv('compare_table.csv', sep='\t')
    print(100. * kol_intersections / len(all_labels))
    print(100. * kol_empty / len(all_labels))


def check_empty(common, in_gk, in_hmm):
    if common == '': common = '-'
    if in_gk == '': in_gk = '-gk-'
    if in_hmm == '': in_hmm = '-hmm-'
    return [common, in_gk, in_hmm]


# ANALYSIS
prefix = '_koala'

# filenames
filename_koala = '/Users/kates/Desktop/prot/UP000000625_ko.txt'  # KOALA
filename_hmm = '/Users/kates/Desktop/kegg-cwl/Workflow/merged_complete_0.00001_625_parsed.txt'
    #'/Users/kates/Desktop/kegg-cwl/Workflow/merged_complete_625_tab_parsed.txt'

koala_dict = get_GK(filename_koala)
hmm_dict = get_hmm(filename_hmm)

print('koala_dict:', len(koala_dict), 'hmm_dict:', len(hmm_dict))
all_labels = set(koala_dict).union(set(hmm_dict))
print(len(all_labels))
print('intersection of labels annotated by GK and HMM', len(set(koala_dict).intersection(set(hmm_dict))))

# compare
#comparison_all_KO(hmm_dict, koala_dict)

# comparison with restrictions
filename_ko_hmm = '/Users/kates/Desktop/prot/all_hmm_ko.txt'
filename_ko_kegg = '/Users/kates/Desktop/prot/all_kegg_ko.txt'

list_ko_hmm, list_ko_kegg = [[], []]
with open(filename_ko_hmm, 'r') as file_hmm_ko:
    for line in file_hmm_ko:
        list_ko_hmm.append(line.strip())

with open(filename_ko_kegg, 'r') as file_kegg_ko:
    for line in file_kegg_ko:
        list_ko_kegg.append(line.strip().split('\t')[0].split(':')[1])

common_ko = set(list_ko_hmm).intersection(set(list_ko_kegg))
print('common_ko', len(common_ko))

compare_tbl = [['key', 'in_common', 'in_GK', 'in_hmm', 'new_kegg']]

kol_intersections, kol_empty_non_annotated, kol_empty_non_intersected = [0 for _ in range(3)]
annotated_both, annotated_only_gk, annotated_only_hmm, kol_empty_non_annotated = [0 for _ in range(4)]
hmm_overperform, gk_overperform, annotated_by_new, annotated_dif, ann_dif_new, ann_dif_no_new = [0, 0, 0, 0, 0, 0]

labels = ['key', 'in_common', 'in_GK', 'in_hmm', 'new_kegg']
table_non_ann_both, table_hmm_better, table_gk_better, table_ann_sim, table_ann_diff = [[labels] for _ in range(5)]

for key in all_labels:
    #new_line = [key]
    new_line_common = [key]
    common, in_hmm, in_gk, new_kegg = ["" for _ in range(4)]
    ko_included_kegg = set(koala_dict[key]).intersection(set(common_ko))  # annotated by GK using old db
    new_kegg = list(set(koala_dict[key]).difference(set(ko_included_kegg)))  # annotated by GK using new db

    if key in hmm_dict and key in koala_dict:  # common (annotated by both)
        annotated_both += 1

        common = ','.join(list(set(hmm_dict[key]).intersection(set(koala_dict[key]))))
        in_hmm = ','.join(list(set(hmm_dict[key]).difference(set(common))))
        in_gk = ','.join(list(set(koala_dict[key]).difference(set(common))))

        if koala_dict[key] == []:  # annotateted by GK as ", but annotated good by HMM
            hmm_overperform += 1
            table_hmm_better.append([key]+check_empty(common, in_gk, in_hmm))
        if koala_dict[key] != [] and common == '':  # annotated differently
            annotated_dif += 1
            table_ann_diff.append([key] + check_empty(common, in_gk, in_hmm))
            if list(set(koala_dict[key]).intersection(set(common_ko))) == []:  # differences only because new KEGG db
                ann_dif_new += 1
            else:
                print(key)  # annotated differently on the ko included in HMM db and new
                ann_dif_no_new += 1
        elif koala_dict[key] != []:
            table_ann_sim.append([key] + check_empty(common, in_gk, in_hmm))

    elif key in hmm_dict:  # annotated by HMM, not GK
        in_hmm = ','.join(hmm_dict[key])
        annotated_only_hmm += 1

    elif key in koala_dict:  # annotated by GK, not annotated by HMM
        in_gk = ','.join(koala_dict[key])
        if koala_dict[key] == []:  # was annotated by GK as empty
            kol_empty_non_annotated += 1
            table_non_ann_both.append([key]+check_empty(common, in_gk, in_hmm))
        else:  # was annotated by GK only
            if koala_dict[key] in list(common_ko):
                gk_overperform += 1
            else:
                annotated_by_new += 1
                table_gk_better.append([key] + check_empty(common, in_gk, in_hmm))
        annotated_only_gk += 1

print('annotated by both tools', annotated_both)
print('annotated by GK', annotated_only_gk, 'empty', kol_empty_non_annotated)
print('GK overperform on existing db', gk_overperform)
print('HMM overperform on existing db', hmm_overperform)
print('Annonated by both and differently', annotated_dif, 'new', ann_dif_new)
print('GK overperform on new db', annotated_by_new)
print('annotated by HMM', annotated_only_hmm)

pd.DataFrame(data=table_non_ann_both).to_csv('table_non_ann_both.csv'+prefix, sep='\t', columns=None, index=False, header=False)
pd.DataFrame(data=table_hmm_better).to_csv('table_hmm_better.csv'+prefix, sep='\t', columns=None, index=False, header=False)
pd.DataFrame(data=table_gk_better).to_csv('table_gk_better.csv'+prefix, sep='\t', columns=None, index=False, header=False)
pd.DataFrame(data=table_ann_sim).to_csv('table_ann_similar.csv'+prefix, sep='\t', columns=None, index=False, header=False)
pd.DataFrame(data=table_ann_diff).to_csv('table_ann_diff.csv'+prefix, sep='\t', columns=None, index=False, header=False)

