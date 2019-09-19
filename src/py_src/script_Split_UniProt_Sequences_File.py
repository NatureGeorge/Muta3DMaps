# @Date:   2019-09-06T10:41:46+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: script_Split_UniProt_Sequences_File.py
# @Last modified time: 2019-09-06T11:06:11+08:00
import os

def split_fasta(file, outputFormat):
    tage = False
    seq = ''
    for line in file:
        if line[0] == '>':
            if tage:
                with open(outputFormat % tage, 'w+') as outputFile:
                    outputFile.write(seq)
                seq = ''
            tage = line.split('|')[1]
            print(tage)
        seq += line
    with open(outputFormat % tage, 'w+') as outputFile:
        outputFile.write(seq)


if __name__ == '__main__':
    dir = '/data/zzf/UniProt_files/'
    outRoute = '/data/zzf/UniProt_files/sub_fasta_files/%s.fasta'
    for file in os.listdir(dir):
        if file[-6:] == '.fasta':
            with open(dir+file, 'rt') as inputFile:
                split_fasta(inputFile, outRoute)
