# @Date:   2019-09-08T22:38:57+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: apply_template.py
# @Last modified time: 2019-09-08T22:49:47+08:00

template_path = r'C:\Users\Nature\Desktop\LiGroup\Work\SIFTS_Plus_Muta_Maps\md\script\pdb_table_template.html'
out_path = r'C:\Users\Nature\Desktop\LiGroup\Work\SIFTS_Plus_Muta_Maps\md\script\output_0908.html'
pdb_list = ['2zqt', '6iwg', '1hho', '2ake', '2xyn', '3a6p', '1dfv', '3hl2']

with open(template_path, 'rt') as inFile:
    with open(out_path, 'a+') as outFile:
        template = inFile.read()
        for pdbId in pdb_list:
            temp_tp = [pdbId]*8
            temp_tp.insert(4, pdbId[1:3])
            temp_tp.insert(7, pdbId[1:3])
            tp = tuple(temp_tp)
            # print(template,tp)
            outFile.write(template%tp)
