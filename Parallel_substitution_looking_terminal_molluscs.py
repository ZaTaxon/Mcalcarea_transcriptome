#!/usr/bin/python3
"look for the parallel amino acid substitutions in given aligment fasta files (need to separate two compared conditions to two fasta)"
#Usage: ./Parallel_substitutions_looking_terminal_molluscs.py fasta1 fasta2 output_name

import sys
import Bio.SeqIO as SeqIO
import re


def write_in_dict(name_dict, count, letter):
    for key, value in name_dict.items():
        temp_val = name_dict[count]
        if letter not in temp_val:
            temp_val += letter
        name_dict[count] = temp_val

def retriever_fasta(sequence_seq, name_dict):
    count = 1
    for letter in sequence_seq:
        if count in name_dict.keys():
            write_in_dict(name_dict, count, letter)
        else:
            name_dict[count] = [letter]
        count += 1

def read_protein_file(file, name_dict):
    seq_records = SeqIO.parse(file, "fasta")
    name_dict = dict()
    for fasta in seq_records:
        name_seq, sequence_seq = fasta.id, str(fasta.seq)
        #print(sequence_seq)
        retriever_fasta(sequence_seq, name_dict)
    print(name_dict)
    return name_dict

def write_file(count, v1, v2, var1, var2, file, file_name):
    with open(file, 'a') as answer:
        answer.write(file_name + '\t')
        answer.write(str(count) + '\t' + ','.join(var1) + '\t' + ','.join(var2) + '\t' + ','.join(v1) + '\t' +
                     ','.join(v2) + '\n')
    print(str(count) + '\t' + ','.join(var1) + '\t' + ','.join(var2) + '\t' + ','.join(v1) + '\t' + ','.join(v2))


def compare_ones(v1, v2, count, file, file_name):
    if v1 == v2:
        #print(count)
        #print(' It is equal')
        var1 = var2 = ['--']
        #write_file(count,v1, v2, var1, var2, file)
    else:
        #print(count)
        #print('Not equal')
        var1 = v1
        var2 = v2
        write_file(count, v1, v2, var1, var2, file, file_name)

def compare_manies(v1, v2, count, file, file_name):
    res = list(set(v1) - set(v2))
    if len(res) > 0 and set(res).issubset(v1):
        #print('In list one')
        var1 = res
        var2 = ['--']
        write_file(count, v1, v2, var1, var2, file, file_name)
    elif len(res) > 0 and set(res).issubset(v2):
        #print('In list two')
        var1 = ['--']
        var2 = res
        write_file(count, v1, v2, var1, var2, file, file_name)
    elif len(res) == 0:
        var1 = var2 = ['--']
        #write_file(count, v1, v2, var1, var2, file)

def compare_dict(dict1, dict2, file, var1, var2, v1, v2, file_name):
    with open(file, 'a') as answer:
        answer.write('File name: '+ str(file_name) + '\n')
        answer.write('File name' + '\t'+'Position' + '\t' + var1 + '\t' + var2 + '\t' + v1 + '\t' + v2 + '\n')
    count = 1
    while count != len(dict1) +1:
        value1 = dict1[count]
        value2 = dict2[count]
        if len(value1) == 1 and len(value2) == 1:
            compare_ones(value1, value2, count, file, file_name)
        elif len(value1) == 1 and len(value2) != 1:
            var1 = value1
            var2 = value2
            write_file(count, value1, value2, var1, var2, file, file_name)
        elif len(value1) != 1 and len(value2) == 1:
            var1 = value1
            var2 = value2
            write_file(count, value1, value2, var1, var2, file, file_name)
        else:
            compare_manies(value1, value2, count, file, file_name)
        count += 1

under = read_protein_file(sys.argv[1], 'un1_aa')
terrest = read_protein_file(sys.argv[2], 'n2_aa')
file_name = re.split(r'\.', sys.argv[1])[0]
table_name = sys.argv[3]
compare_dict(under, terrest, table_name, 'Niche1', 'Niche2', 'Native_N1', 'Native_N2', file_name)


#print(under)
#print(terrest)

