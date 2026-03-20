import pandas as pd 
import os
import subprocess
from sys import argv

Pre = 'Men-IGT'
inf2 = f'{Pre}_assem.txt'
inf = f'{Pre}_assem.kraken2.txt'
afile = pd.read_table(inf,header=None)
afile = afile.loc[afile[5].str.contains('Neisseria meningitidis'),]
taxidlist = afile[4].to_list()
bfile = pd.read_table(inf2,header=None)
bfile = bfile.loc[bfile[2].isin(taxidlist),]
bfile[[1]].to_csv('NMList.txt',sep='\t',index=False,header=False)
subprocess.run(f'seqkit grep -f NMList.txt {Pre}.final.fasta > {Pre}.NM.fasta',shell=True)
if not os.path.isdir('new'):
    os.makedirs('new')
subprocess.run(f'mv {Pre}.NM.fasta new',shell=True)
for i in [80,85,90,95]:
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output  -n cdhit cd-hit-est -i new/{Pre}.NM.fasta -o new/{Pre}.NMfilter{i}.fasta -c 0.{i} -T 10',shell=True)
subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output  -n cm2  checkm2 predict -i new -x fasta -o checkm2_out1 -t 10 --force --database_path /data1/shanghai_pip/meta_genome/uniref100.KO.1.dmnd',shell=True)
