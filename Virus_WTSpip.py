#!/home/wusihao/miniconda3/envs/Vlib/bin/python
# -*- coding: utf-8 -*-
#--------env Immu_pip------------
import sys
import os
import argparse
import datetime
import subprocess
import json
import pandas as pd
import re
import time
from collections import Counter
import dask.dataframe as dd
__author__='wsh'
__version__='1.0.0'
__date__='20250427'
pd.options.mode.chained_assignment = None
#---------更新档案--------------
##########v1.0############
#1.数据质控
#2.去宿主
#3.kraken2鉴定物种 & 比对鉴定正负链比例
#4.宏基因组组装
#5.组装后结果物种鉴定（blast virsort checkV）
#---------parameters ------------------
parser = argparse.ArgumentParser(description='Immu pipline')
parser.add_argument('--list','-l',type=str,default=False,help='文件列表')
parser.add_argument('--step','-s',type=str,default='All',help='分析内容')
parser.add_argument('--queue','-q',type=int,default=2,help='队列数')
parser.add_argument('--threads','-t',type=int,default=10,help='线程数量')
parser.add_argument('--output','-o',type=str,default=False,help='输出文件')
parser.add_argument('--downsample','-d',type=str,default='20000000',help='下采样')
argv = parser.parse_args()
Tims = datetime.datetime.now().strftime('%m%d_%H%M')
stime=time.time()
def format_seconds(seconds):
    # 将秒数转换为整数
    total_seconds = int(seconds) 
    # 计算小时、分钟和秒
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    # 格式化输出
    return f"{hours}小时{minutes}分钟{seconds}秒"

listn = argv.list
threads = argv.threads
ofn = argv.output
step = argv.step 
queue = argv.queue
downsample = argv.downsample
print(f''' 
input:\t{listn}
threads:\t{threads}
output:\t{ofn}
step:\t{step}
queue:\t{queue}
downsample:\t{downsample}
''')
if step == 'All':
    workstep == 'QC,RmC,Kk2&mapping,Assemble,IdenS'
else:
    workstep = step


def Preprocess(Pre,fq1,fq2,threads,downsample):
    #1.table1 原始数据过滤表 2.table2 过滤后数据库质量表 3.flash merge表格
    with open('Preprocess.log','w') as f:
        cfq1 = f'{Pre}.clean.R1.fastq.gz'
        cfq2 = f'{Pre}.clean.R2.fastq.gz'
        if not os.path.isfile(f'{Pre}.sam.1.fq') or os.path.isfile(f'{Pre}.sam.1.fq.gz'):
            if downsample != 'nosample':
                if '0.' not in downsample:
                    #subprocess.run(f'''conda run --no-capture-output -n tNGS rasusa reads -n {downsample} -s 13 {fq1} {fq2} -o {Pre}.sam.1.fq -o {Pre}.sam.2.fq ''',shell=True,stdout=f,stderr=f)
                    subprocess.run(f'''conda run --no-capture-output -n tNGS rasusa reads -b {downsample} -s 13 {fq1} {fq2} -o {Pre}.sam.1.fq -o {Pre}.sam.2.fq ''',shell=True,stdout=f,stderr=f)
                else:
                    subprocess.run(f'''conda run --no-capture-output -n tNGS rasusa reads -f {downsample} -s 13 {fq1} {fq2} -o {Pre}.sam.1.fq -o {Pre}.sam.2.fq ''',shell=True,stdout=f,stderr=f)

            else:
                if fq1.endswith('gz'):
                    subprocess.run(f'''ln -s {fq1} {Pre}.sam.1.fq.gz''',shell=True,stdout=f,stderr=f)
                    subprocess.run(f'''ln -s {fq2} {Pre}.sam.2.fq.gz''',shell=True,stdout=f,stderr=f)
                else:
                    subprocess.run(f'''ln -s {fq1} {Pre}.sam.1.fq''',shell=True,stdout=f,stderr=f)
                    subprocess.run(f'''ln -s {fq2} {Pre}.sam.2.fq ''',shell=True,stdout=f,stderr=f)

        if os.path.isfile(f'{Pre}.sam.1.fq.gz'):
            subprocess.run(f'fastp -i {Pre}.sam.1.fq.gz -I {Pre}.sam.2.fq.gz -o {cfq1} -O {cfq2} --cut_front --cut_tail --cut_front_mean_quality 20 -e 20 -n 30 --length_required 50 -w {threads} --json {Pre}.json',shell=True,stdout=f,stderr=f)
        else:
            subprocess.run(f'fastp -i {Pre}.sam.1.fq -I {Pre}.sam.2.fq -o {cfq1} -O {cfq2} --cut_front --cut_tail --cut_front_mean_quality 20 -e 20 -n 30 --length_required 50 -w {threads} --json {Pre}.json',shell=True,stdout=f,stderr=f)



        tmpdict =json.load(open(f'{Pre}.json','r'))
        QCdb = pd.concat([pd.DataFrame(tmpdict['summary']['before_filtering'],index=['before']),pd.DataFrame(tmpdict['summary']['after_filtering'],index=['after'])]).reset_index()
        QCdb['SamName'] = Pre
        QCdb.to_csv('QC_summary.tsv',sep='\t',index=False)


def removeHost(Pre,fq1,fq2,threads,Hostref):
    with open('rmhost.log','w') as f:
        if Hostref != 'noref':
            if not os.path.isfile(f'{Pre}.rmhost.1.fq'):
                subprocess.run(f'/home/wusihao/anaconda3/condabin/conda run --no-capture-output -n Vlib kneaddata --input1 {fq1} --input2 {fq2} -db {Hostref} -db /Data/RawData/Database/KneaDB/silva/  -t {threads} --output {Pre}_knea_output --bypass-trf --bypass-trim ',shell=True)
                rawfq1 = os.popen(f'''find {Pre}_knea_output -name "*_kneaddata_paired_1.fastq"''').read().strip()
                rawfq2 = os.popen(f'''find {Pre}_knea_output -name "*_kneaddata_paired_2.fastq"''').read().strip()
                subprocess.run(f'''ln -s {rawfq1} {Pre}.rmhost.1.fq''',shell=True,stdout=f,stderr=f)
                subprocess.run(f'''ln -s {rawfq2} {Pre}.rmhost.2.fq''',shell=True,stdout=f,stderr=f)
        else:
            #if not os.path.isfile(f'{Pre}.rmhost.1.fq'):
             #   subprocess.run(f'/home/wusihao/anaconda3/condabin/conda run --no-capture-output -n Vlib kneaddata --input1 {fq1} --input2 {fq2}  -db /Data/RawData/Database/KneaDB/silva/  -t {threads} --output {Pre}_knea_output --bypass-trf --bypass-trim ',shell=True)
             #   rawfq1 = os.popen(f'''find {Pre}_knea_output -name "*_kneaddata_paired_1.fastq"''').read().strip()
              #  rawfq2 = os.popen(f'''find {Pre}_knea_output -name "*_kneaddata_paired_2.fastq"''').read().strip()
               # subprocess.run(f'''ln -s {rawfq1} {Pre}.rmhost.1.fq''',shell=True,stdout=f,stderr=f)
               # subprocess.run(f'''ln -s {rawfq2} {Pre}.rmhost.2.fq''',shell=True,stdout=f,stderr=f)

            subprocess.run(f'''ln -s {fq1} {Pre}.rmhost.1.fq''',shell=True,stdout=f,stderr=f)
            subprocess.run(f'''ln -s {fq2} {Pre}.rmhost.2.fq''',shell=True,stdout=f,stderr=f)

import pandas as pd

def kk2(pre, dblist):
    tmpdict = {}
    tmpdict1 = {}

    # 分类级别初始化
    levels = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
    current = {lvl: '-' for lvl in levels}

    with open(f'{pre}.report.txt', encoding='utf-8') as f:
        for line in f:
            Abun, ReadsC, ReadsO, level, taxid, tmpName = line.strip().split('\t')
            tmpName = tmpName.strip()

            # 更新分类级别
            if level in levels:
                idx = levels.index(level)
                for l in levels[idx+1:]:
                    current[l] = '-'
                current[level] = tmpName

            if level == 'S':
                tmpdict[tmpName] = {
                    '序列数量': ReadsC,
                    '比例': Abun,
                    'NCBI物种号': taxid,
                    '分类': current['D'],
                    '属': current['G']
                }
            elif level == 'U':
                tmpdict['unclassified'] = {
                    '序列数量': ReadsC,
                    '比例': Abun,
                    'NCBI物种号': 0,
                    '分类': 'unclassified',
                    '属': 'unclassified'
                }
                tmpdict1['unclassified'] = {
                    '序列数量': ReadsC,
                    '比例': Abun,
                    'NCBI物种号': 0,
                    '分类': 'unclassified',
                    '属': 'unclassified',
                    '种': 'unclassified'
                }
            elif level in ['S1', 'S2', 'S3']:
                tmpdict1[tmpName] = {
                    '序列数量': ReadsC,
                    '比例': Abun,
                    'NCBI物种号': taxid,
                    '分类': current['D'],
                    '属': current['G'],
                    '种': current['S']
                }

    # DataFrame 构建
    kkdb = pd.DataFrame(tmpdict).T
    kkdb1 = pd.DataFrame(tmpdict1).T
    kkdb['物种'] = kkdb.index
    kkdb1['亚种'] = kkdb1.index

    # 合并元数据
    metadb = pd.read_table(dblist, encoding='utf-8')
    metadb['taxid'] = metadb['taxid'].astype(str)
    kkdb = kkdb.merge(metadb, left_on='NCBI物种号', right_on='taxid', how='left')
    kkdb1 = kkdb1.merge(metadb, left_on='NCBI物种号', right_on='taxid', how='left')

    kkdb['覆盖度%'] = ''
    kkdb.fillna('-', inplace=True)
    kkdb1.fillna('-', inplace=True)

    # 列顺序
    kkdb = kkdb[['物种', '序列数量', '比例', '覆盖度%', '中文名', '致病性', '可能引起的疾病', '危害程度等级', '分类', '属', 'NCBI物种号']]
    kkdb1 = kkdb1[['亚种', '序列数量', '比例', '中文名', '致病性', '可能引起的疾病', '危害程度等级', '分类', '属', '种', 'NCBI物种号']]

    # 类型转换与排序
    for df, key in [(kkdb, '物种'), (kkdb1, '亚种')]:
        try:
            df['序列数量'] = df['序列数量'].astype(int)
        except Exception:
            df['序列数量'] = pd.to_numeric(df['序列数量'], errors='coerce').fillna(0).astype(int)
        df.sort_values('序列数量', ascending=False, inplace=True)

    kkdb.to_csv(f'{pre}.list.txt', sep='\t', index=False, encoding='utf-8')
    kkdb1.to_csv(f'{pre}.list2.txt', sep='\t', index=False, encoding='utf-8')

    # 英文表头
    kkdb.rename(columns={
        '物种': 'Species', '序列数量': 'Reads Num', '比例': 'Abundance', 'NCBI物种号': 'TaxID',
        '中文名': 'ChineseName', '致病性': 'Pathogenicity', '可能引起的疾病': 'Cause Disease',
        '危害程度等级': 'Hazard Ranking', '分类': 'Kingdom', '属': 'Genus'
    }, inplace=True)
    kkdb1.rename(columns={
        '亚种': 'Subspecies', '序列数量': 'Reads Num', '比例': 'Abundance', 'NCBI物种号': 'TaxID',
        '中文名': 'ChineseName', '致病性': 'Pathogenicity', '可能引起的疾病': 'Cause Disease',
        '危害程度等级': 'Hazard Ranking', '分类': 'Kingdom', '属': 'Genus', '种': 'Species'
    }, inplace=True)

    kkdb.to_csv(f'{pre}.list_EN.txt', sep='\t', index=False, encoding='utf-8')
    kkdb1.to_csv(f'{pre}.list2_EN.txt', sep='\t', index=False, encoding='utf-8')

def compareid(level1,level2,rawlist):
    print(level1)
    rlevel1 = [i for i in ['R','D','K','P','C','O','F','G','S'] if i in level1][0]
    rlevel2 = [i for i in ['R','D','K','P','C','O','F','G','S'] if i in level2][0]

    if rawlist.index(rlevel1) == rawlist.index(rlevel2):
        if level1 > level2:
            return 1 
        elif level1 == level2:
            return 0
        else:
            return -1
    elif rawlist.index(rlevel1) > rawlist.index(rlevel2):
        return 1
    else:
        return -1


def proc_kra(kraken,tax,lel='S'):
    tmplist = [tax]
    if [i for i in ['R','D','K','P','C','O','F','G','S'] if i in lel]:    
        rawlist = ['R','D','K','P','C','O','F','G','S']
    else:
        rawlist = ['S1','S2','S3','S4','S5','S6']
    if tax != 0:
        kradb = pd.read_table(kraken,header=None)
        #kradb = kradb[(kradb[3]==lel)&(kradb[4]==tax)]
        kradb[4] = kradb[4].astype('str')
        tmpindex = kradb[(kradb[3]==lel)&(kradb[4]==str(tax))].index.tolist()[0]+1
        if tmpindex <=  kradb.shape[0]-1:
            def getlindex(tmpindex):
                tmpl = kradb.iloc[tmpindex,3]
                for tl in rawlist:
                    if tl in tmpl:
                        if tl == tmpl:
                            tmpl = tl
                            tmlindex = rawlist.index(tmpl)
                        else:
                            tmpl = tl
                            tmlindex = rawlist.index(tmpl)+1
                        return tmlindex
            while compareid(kradb.iloc[tmpindex,3],lel,rawlist) == 1 and tmpindex <= kradb.shape[0]-2:
                tmplist.append(kradb.iloc[tmpindex,4])
                tmpindex+=1
    return tmplist

def exreadsID(taxlist,kraresult,fq1,fq2=0):
    Maintax = taxlist[0]
    #kraredb = pd.read_table(kraresult,header=None)
    kraredb = dd.read_csv(kraresult, header=None, usecols=[1, 2], dtype={1:'str',2:'int32'},sep='\t')
    tmp2db = kraredb[kraredb[2].isin(taxlist)]
    tmp1db = pd.DataFrame(tmp2db[1].unique())
    tmp2db.to_csv(f'{Maintax}.id.tsv',sep='\t',index=False)
    #tmp1db.to_csv(f'{Maintax}_fqID.txt',index=False,header=False)
    pd.DataFrame(tmp1db).to_csv(f'{Maintax}_fqID.txt', index=False, header=False)
    subprocess.run(f'head -n 1 {Maintax}_fqID.txt > tt.txt',shell=True)
    subprocess.run(f'''cut -d '/' -f1 {Maintax}_fqID.txt|sort -u > {Maintax}.listID.txt''',shell=True)
    if os.popen(f'''head -n 1 tt.txt''').read().strip().endswith('/1') or os.popen(f'''head -n 1 tt.txt''').read().strip().endswith('/2'):
        subprocess.run(f'''sed 's/$/\/1/' {Maintax}.listID.txt > {Maintax}.listID1.txt''',shell=True)
        subprocess.run(f'''sed 's/$/\/2/' {Maintax}.listID.txt > {Maintax}.listID2.txt''',shell=True)
        subprocess.run(f'seqkit grep -f {Maintax}.listID1.txt {fq1} > {Maintax}.1.fastq',shell=True)
        if fq2:
            subprocess.run(f'seqkit grep -f {Maintax}.listID2.txt {fq2} > {Maintax}.2.fastq',shell=True)
    else:
        subprocess.run(f'seqkit grep -f {Maintax}.listID.txt {fq1} > {Maintax}.1.fastq',shell=True)
        if fq2:
            subprocess.run(f'seqkit grep -f {Maintax}.listID.txt {fq2} > {Maintax}.2.fastq',shell=True)


def exreadsID2(taxlist,kraresult,fq1,fq2=0):
    Maintax = taxlist[0]
    #kraredb = pd.read_table(kraresult,header=None)
    kraredb = dd.read_csv(kraresult, header=None, usecols=[1, 2], dtype={1:'str',2:'int32'},sep='\t')
    tmp1db = pd.DataFrame(kraredb[kraredb[2].isin(taxlist)][1].unique())
    #tmp1db.to_csv(f'{Maintax}_fqID.txt',index=False,header=False)
    pd.DataFrame(tmp1db).to_csv(f'{Maintax}_fqID.txt', index=False, header=False)
    subprocess.run(f'head -n 1 {Maintax}_fqID.txt > tt.txt',shell=True)
    subprocess.run(f'''cut -d '/' -f1 {Maintax}_fqID.txt|sort -u > {Maintax}.listID.txt''',shell=True)
    if os.popen(f'''head -n 1 tt.txt''').read().strip().endswith('/1') or os.popen(f'''head -n 1 tt.txt''').read().strip().endswith('/2'):
        subprocess.run(f'''sed 's/$/\/1/' {Maintax}.listID.txt > {Maintax}.listID1.txt''',shell=True)
        subprocess.run(f'seqkit grep -f {Maintax}.listID1.txt {fq1} --delete-matched | seqkit fq2fa | seqkit sample -n 100 -w0 > {Maintax}.sam.fasta',shell=True)
    else:
        subprocess.run(f'seqkit grep -f {Maintax}.listID.txt {fq1} --delete-matched | seqkit fq2fa | seqkit sample -n 100 -w0 > {Maintax}.sam.fasta',shell=True)

def asmvirus(Pre,threads):
    #1.提取病毒序列 2.去冗余 3.无参组装
    level = 'D'
    krakenfile = f'{Pre}.report.txt'
    tkid = 10239
    tkidl = os.popen(f'grep {tkid} {krakenfile}|cut -f4').read().strip()
    taxlist1 = proc_kra(krakenfile,tkid,level)
    taxlist1 = [int(i) for i in taxlist1]
    exreadsID(taxlist1,f'{Pre}.txt',f'{Pre}.1.fq',f'{Pre}.2.fq')
    newPre = Pre.replace('.rmhost','')
    with open('megahit.log','w') as f:
        if os.path.isdir(f'{newPre}_megahit_out'):
            subprocess.run(f'rm -r {newPre}_megahit_out',shell=True)
        subprocess.run(f'megahit -1 10239.1.fastq -2 10239.2.fastq -t {threads} -o {newPre}_megahit_out',shell=True,stdout=f,stderr=f)

def get_bestr(tmpdb,cov=0.9,iden=90):
    tmpdb1 = tmpdb.loc[(tmpdb['cov']>=cov) & (tmpdb[2]>=iden)]
    if tmpdb1.shape[0] > 0:
        tmpdb1 = tmpdb1.sort_values(11,ascending=False)
        tmpdb1['cer'] = 'HQ'
    else:
        tmpdb1 = tmpdb.loc[(tmpdb['cov']>=0.5) & (tmpdb[2]>=50)]
        if tmpdb1.shape[0] >0:
            tmpdb1 = tmpdb1.sort_values(11,ascending=False)
            tmpdb1['cer'] = 'MQ'
        else:
            tmpdb1 = tmpdb
            tmpdb1 = tmpdb1.sort_values(11,ascending=False)
            tmpdb1['cer'] = 'LQ'

    return tmpdb1.head(1)

def get_nt(cname,ntdb):
    tmpdb = ntdb.loc[ntdb['Contig']==cname]
    if tmpdb.shape[0] == 1:
        cer = tmpdb['cer'].tolist()[0]
        Spe = tmpdb['Spe'].tolist()[0]
        tmpr = f'{Spe}_{cer}'
    else:
        tmpr = 'NoIden'
    return tmpr

def getSpe(x):
    if x['RVDB'].startswith('nan_'):
        Spe = x['NT'].split('|')[0]
    else:
        Spe = ' '.join(x['RVDB'].split('_')[:-1])
    Spe=Spe.split('strain')[0].strip()
    return Spe
def compare_rvdb(tmpdb,protdb):
    Qlist = ['LQ','MQ','HQ']
    ctg = tmpdb['Contig']
    cer = tmpdb['cer']
    Spe = tmpdb['Spe']
    if ctg in protdb[0]:
        protSpe = protdb.loc[protdb[0]==ctg,'species'].tolist()[0]
        if protSpe == Spe:
            try:
                ncer = Qlist[Qlist.index(cer)+1]
            except:
                ncer = cer
    else:
        ncer = cer
    return ncer
def ContigIden(Pre,threads):
    if os.path.isfile(f'{Pre}_megahit_out/final.contigs.fa'):
        #1.比对RVDB数据库 2.将能比对到RVDB数据库的contig比对到nt_core 3. 跑checkV和 vitsort2 验证
        #20250626 优化增加RVDB prot数据库 合并RVDB和RVDB_prot的结果，哪个有算哪个，都有就按prot的来
        if not os.path.isfile(f'{Pre}.RVDB.out.tsv'):
            if not os.path.isfile(f'{Pre}.contigs.fasta'):
                subprocess.run(f'ln -s {Pre}_megahit_out/final.contigs.fa ./{Pre}.contigs.fasta',shell=True)
            if not os.path.isfile(f'{Pre}.out.txt') or os.path.getsize(f'{Pre}.out.txt') == 0:
                subprocess.run(f'blastn -db /Data/RawData/Database/RVDB/RVDB -num_threads {threads} -out {Pre}.out.txt -query {Pre}.contigs.fasta -max_target_seqs 10 -evalue 0.05 -word_size 28 -outfmt 6',shell=True)
            subprocess.run(f'seqkit fx2tab -n -l -i  {Pre}.contigs.fasta > {Pre}.contig.tsv',shell=True)
        if not os.path.isfile(f'{Pre}.rvdbprot.txt'):
            subprocess.run(f'diamond blastx -d /Data/RawData/Database/RVDB/RVDB_prot.dmnd -q {Pre}.contigs.fasta --evalue 0.05 --outfmt 6 --out {Pre}.rvdbprot.txt',shell=True)

        blastdb = pd.read_table(f'{Pre}.out.txt',header=None)
        #--add rvdbprot step
        if os.path.getsize(f'{Pre}.rvdbprot.txt'):
            blastdb2 = pd.read_table(f'{Pre}.rvdbprot.txt',header=None)
            protmetadb = pd.read_table('/Data/RawData/Database/RVDB/RVDB_prot.meta.tsv')
            blastdb2 = blastdb2.merge(protmetadb,left_on=1,right_on='name')

        contigsdb = pd.read_table(f'{Pre}.contig.tsv',header=None)
        metadb = pd.read_table('/Data/RawData/Database/RVDB/meta.tsv')
        blastdb['contig'] = blastdb[1].str.split('|').str[2]
        blastdb = blastdb.merge(metadb,left_on='contig',right_on='序列名称')
        blastdb = blastdb.merge(contigsdb,on=0)
        blastdb['cov'] = blastdb[3]/blastdb['1_y']*100
        blastdb[[0,'1_x',2,'cov']].to_csv(f'{Pre}.RVDB.raw.tsv',sep='\t',index=False)
        blastdb = blastdb[(blastdb[2]>=85) & (blastdb['cov']>=85)]
        blastout = pd.DataFrame(blastdb.groupby(0).apply(lambda x:x.head(1))['种'].value_counts()).reset_index()
        blastout.to_csv(f'{Pre}.RVDB.out.tsv',sep='\t',index=False)
        #2
        blastdb[[0]].drop_duplicates().to_csv('RVDB.contig.txt',sep='\t',index=False,header=False)
        if not os.path.isfile(f'{Pre}.nt_core.out.tsv'):
            if blastdb[[0]].drop_duplicates().shape[0] > 0:
                subprocess.run(f'seqkit grep -f RVDB.contig.txt {Pre}.contigs.fasta > {Pre}.RVDB.fasta',shell=True)
            if not os.path.isfile(f'{Pre}.ntout.txt') or os.path.getsize(f'{Pre}.ntout.txt') == 0:
                subprocess.run(f'''blastn -db /dev/shm/core_nt/core_nt -num_threads {threads} -out {Pre}.ntout.txt -query {Pre}.RVDB.fasta -max_target_seqs 10 -evalue 0.05 -word_size 28 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname' ''',shell=True)
        #3
        subprocess.run(f'/home/wusihao/anaconda3/condabin/conda run --no-capture-output -n VFind virsorter run -w {Pre}_virsort_out -i {Pre}.contigs.fasta -j {threads} all  --include-groups dsDNAphage,ssDNA,RNA',shell=True)
        subprocess.run(f'/home/wusihao/anaconda3/condabin/conda run --no-capture-output -n VFind checkv end_to_end {Pre}.contigs.fasta {Pre}_checkv_output -d /Data/RawData/Database/checkvDB/checkv-db-v1.5 -t 10',shell=True)
        if not os.path.isfile(f'{Pre}.nt.taxonlist.tsv'):
            subprocess.run(f'''cut -f13 {Pre}.ntout.txt |sort -u|taxonkit reformat -I 1 -f "{{s}}|{{t}}" -F > {Pre}.nt.taxonlist.tsv''',shell=True)
        taxiddb = pd.read_table(f'{Pre}.nt.taxonlist.tsv',header=None,names=['Taxid','Spe'])
        Rawdb = pd.read_table(f'{Pre}.contig.tsv',header=None,names=['Contig','length'])
        RVdb = pd.read_table(f'{Pre}.out.txt',header=None)
        if RVdb.shape[0] > 0:
            RVdb1 = Rawdb.merge(RVdb,left_on='Contig',right_on=0)
            RVdb1['cov'] = RVdb1[3]/RVdb1['length']
            RVdbn = RVdb1.groupby('Contig').apply(lambda x:get_bestr(x))
            RVdbn = RVdbn.reset_index(drop=True)
            RVdbn['Spe'] = RVdbn[1].str.split('|').str[4]
            if os.path.getsize(f'{Pre}.rvdbprot.txt'): 
                RVdbn['cer'] = RVdbn.apply(lambda x:compare_rvdb(x,blastdb2),axis=1)
            #print(RVdbn[['Contig','length','Spe',2,'cov','cer']])
            if os.path.isfile(f'{Pre}.ntout.txt') and os.path.getsize(f'{Pre}.ntout.txt') != 0:
                NTdb = pd.read_table(f'{Pre}.ntout.txt',header=None)
                NTdb1 = Rawdb.merge(NTdb,left_on='Contig',right_on=0)
                NTdb1['cov'] = NTdb1[3]/NTdb1['length']
                NTdbn = NTdb1.groupby('Contig').apply(lambda x:get_bestr(x))
                NTdbn = NTdbn.reset_index(drop=True)
                NTdbn = NTdbn.merge(taxiddb,left_on=12,right_on='Taxid')
                Rawdb['NT'] = Rawdb.apply(lambda x:get_nt(x['Contig'],NTdbn),axis=1)
            else:
                Rawdb['NT'] = 'NoIden'
            Rawdb['RVDB'] = Rawdb.apply(lambda x:get_nt(x['Contig'],RVdbn),axis=1)
            
            #3.checkV
            if os.path.isfile(f'{Pre}_checkv_output/contamination.tsv'):
                checkvdb = pd.read_table(f'{Pre}_checkv_output/contamination.tsv')
                Rawdb = Rawdb.merge(checkvdb,left_on='Contig',right_on='contig_id',how='left')
            else:
                checkvheaderlist = ['total_genes'   ,  'viral_genes'   ,  'host_genes'   ,   'provirus'     ,   'proviral_length', 'host_length'   ,  'region_types'   , 'region_lengths',  'region_coords_bp'     ,   'region_coords_genes'  ,   'region_viral_genes'  ,    'region_host_genes']
                for tmpheader in checkvheaderlist:
                    Rawdb[tmpheader] = 0
            #4.virsort
            if os.path.isfile(f'{Pre}_virsort_out/final-viral-score.tsv') and os.path.getsize(f'{Pre}_virsort_out/final-viral-score.tsv') != 0:
                virsortdb = pd.read_table(f'{Pre}_virsort_out/final-viral-score.tsv',usecols=['seqname','hallmark','max_score','max_score_group'])
                virsortdb['Contig'] = virsortdb['seqname'].str.split('|').str[0]
                Rawdb = Rawdb.merge(virsortdb,on='Contig',how='left')
            else:
                virsorttmph = [ 'seqname',   'max_score'  ,     'max_score_group',  'hallmark' ]
                for vs2h in virsorttmph:
                    Rawdb[virsorttmph] = 0
            Rawdb = Rawdb[['Contig',    'length'    ,'NT',  'RVDB', 'total_genes'   ,'viral_genes', 'host_genes',   'provirus', 'seqname',  'max_score' ,'max_score_group','hallmark']]
            #filter standard https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v2?version_warning=no&step=4
            Rawdb.to_csv(f'{Pre}.contigiden.tsv',sep='\t',index=False)
            Rawdb = Rawdb.loc[(Rawdb['hallmark']>2) | (Rawdb['viral_genes']>0) | (Rawdb['viral_genes']>0) | ((Rawdb['viral_genes']==0) & (Rawdb['host_genes']==0))]
            Rawdb = Rawdb.loc[~((Rawdb['viral_genes']==0)&(Rawdb['host_genes']>1))]
            Rawdb = Rawdb.loc[~((Rawdb['viral_genes']==0)&(Rawdb['host_genes']==1)&(Rawdb['length']<10000))]
            #---下一步有可能过滤掉新病毒------
            Rawdb1 = Rawdb.copy()
            Rawdb = Rawdb.loc[~((Rawdb['NT']=='NoIden')&(Rawdb['RVDB']=='NoIden'))]
            #---判断新病毒----
            #-give taxid for contigs 
            #1.如果只有一个有那个用哪个 2如果有两个判断是否一致，一致直接用，不一致判断RVDB和NT哪个质量高，如果质量一致用RVDB,
            Rawdb['FinalSpe'] = Rawdb.apply(lambda x:getSpe(x),axis=1)
            HQSpelist = list(set(Rawdb.loc[((Rawdb1['NT'].str.contains('_HQ'))|(Rawdb['RVDB'].str.contains('_HQ')))]['FinalSpe'].tolist()))
            MNvirusdb = Rawdb1.loc[~((Rawdb1['NT'].str.contains('_HQ'))|(Rawdb1['RVDB'].str.contains('_HQ')))]
            if MNvirusdb.shape[0] >0:
                MNvirusdb['FinalSpe'] = MNvirusdb.apply(lambda x:getSpe(x),axis=1)
                MNvirusdb = MNvirusdb.loc[~MNvirusdb['FinalSpe'].isin(HQSpelist)]
                MNvirusdb.to_csv(f'{Pre}.maynewvirus.tsv',sep='\t',index=False)
            Rawdb[['FinalSpe']].drop_duplicates().to_csv(f'{Pre}.finalSpe.list.txt',sep='\t',index=False,header=False)
            subprocess.run(f'''cat {Pre}.finalSpe.list.txt|taxonkit name2taxid > {Pre}.finalSpeid.list.txt ''',shell=True)
            finalSpedb = pd.read_table(f'{Pre}.finalSpeid.list.txt',header=None,names=['FinalSpe','FinalSpeTaxid'])
            finalSpedb.fillna('noid',inplace=True)
            Rawdb = Rawdb.merge(finalSpedb,on='FinalSpe')
            Rawdb.to_csv(f'{Pre}.contigiden_filter.tsv',sep='\t',index=False)

def kk2cmbl(Pre,threads):
    QCdb = pd.read_table('QC_summary.tsv')
    Allreads = QCdb['total_reads'].tolist()[1]
    Spekk2db = pd.read_table(f'{Pre}.rmhost.list.txt')
    Spekk2db['RPM'] = Spekk2db['序列数量']/Allreads*1000000
    Spekk2db = Spekk2db.loc[Spekk2db['RPM']>=1]
    Spekk2db['NCBI物种号'] = Spekk2db['NCBI物种号'].astype('str')
    Rawdb = pd.read_table(f'{Pre}.contigiden_filter.tsv')
    subprocess.run(f'''cut -f2 {Pre}.finalSpeid.list.txt|taxonkit lineage -t|cut -f3 |sort -u > Allidlineage.txt''',shell=True)
    Alllist = []
    with open('Allidlineage.txt') as f:
        for line in f:
            line = line.strip()
            Alllist.extend(line.split(';'))
    Alllist = list(set(Alllist))
    Spekk2db = Spekk2db.loc[~Spekk2db['NCBI物种号'].isin(Alllist)]
    Spekk2db = Spekk2db[Spekk2db['分类']=='Viruses']
    Spekk2db.to_csv(f'{Pre}.readsonly.tsv',sep='\t',index=False)
    for tmpid in Spekk2db['NCBI物种号'].tolist():
        level = 'S'
        krakenfile = f'{Pre}.rmhost.report.txt'
        tkid = tmpid
        tkidl = os.popen(f'grep {tkid} {krakenfile}|cut -f4').read().strip()
        taxlist1 = proc_kra(krakenfile,tkid,level)
        taxlist1 = [int(i) for i in taxlist1]
        exreadsID2(taxlist1,f'{Pre}.rmhost.txt','10239.1.fastq',f'10239.2.fastq')
        subprocess.run(f'blastn -db /Data/RawData/Database/RVDB/RVDB -num_threads {threads} -out {tkid}.out.txt -query {tkid}.sam.fasta -max_target_seqs 10 -evalue 0.05 -word_size 28 -outfmt 6',shell=True)
        subprocess.run(f'seqkit fx2tab -n -l -i  {tkid}.sam.fasta > {tkid}.sam.tsv',shell=True)
        if os.path.isfile(f'{tkid}.out.txt') and os.path.getsize(f'{tkid}.out.txt') != 0:
            blastdb = pd.read_table(f'{tkid}.out.txt',header=None)
            contigsdb = pd.read_table(f'{tkid}.sam.tsv',header=None)
            metadb = pd.read_table('/Data/RawData/Database/RVDB/meta.tsv')
            blastdb['contig'] = blastdb[1].str.split('|').str[2]
            blastdb = blastdb.merge(metadb,left_on='contig',right_on='序列名称')
            blastdb = blastdb.merge(contigsdb,on=0)
            blastdb['cov'] = blastdb[3]/blastdb['1_y']*100
            blastdb[[0,'1_x',2,'cov']].to_csv(f'{tkid}.RVDB.raw.tsv',sep='\t',index=False)
            blastdb = blastdb[(blastdb[2]>=85) & (blastdb['cov']>=85)]
            blastout = pd.DataFrame(blastdb.groupby(0).apply(lambda x:x.head(1))['种'].value_counts()).reset_index()
            blastout.to_csv(f'{tkid}.RVDB.out.tsv',sep='\t',index=False)
            #2
            blastdb[[0]].drop_duplicates().to_csv(f'{tkid}.RVDB.contig.txt',sep='\t',index=False,header=False)
            if os.path.getsize(f'{tkid}.RVDB.contig.txt') > 0:
                if not os.path.isfile(f'{tkid}.nt_core.out.tsv'):
                    if blastdb[[0]].drop_duplicates().shape[0] > 0:
                        subprocess.run(f'seqkit grep -f {tkid}.RVDB.contig.txt {tkid}.sam.fasta > {tkid}.RVDB.fasta',shell=True)
                    if not os.path.isfile(f'{tkid}.ntout.txt') or os.path.getsize(f'{tkid}.ntout.txt') == 0:
                        subprocess.run(f'''blastn -db /dev/shm/core_nt/core_nt -num_threads {threads} -out {tkid}.ntout.txt -query {tkid}.RVDB.fasta -max_target_seqs 10 -evalue 0.05 -word_size 28 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname' ''',shell=True)
                if not os.path.isfile(f'{tkid}.nt.taxonlist.tsv'):
                    subprocess.run(f'''cut -f13 {tkid}.ntout.txt |sort -u|taxonkit reformat -I 1 -f "{{s}}|{{t}}" -F > {tkid}.nt.taxonlist.tsv''',shell=True)
                taxiddb = pd.read_table(f'{tkid}.nt.taxonlist.tsv',header=None,names=['Taxid','Spe'])
                subprocess.run(f'seqkit fx2tab -n -l {tkid}.sam.fasta > {tkid}.sam.tsv',shell=True)
                Rawdb = pd.read_table(f'{tkid}.sam.tsv',header=None,names=['Contig','length'])
                RVdb = pd.read_table(f'{tkid}.out.txt',header=None)
                RVdb1 = Rawdb.merge(RVdb,left_on='Contig',right_on=0)
                RVdb1['cov'] = RVdb1[3]/RVdb1['length']
                RVdbn = RVdb1.groupby('Contig').apply(lambda x:get_bestr(x))
                RVdbn = RVdbn.reset_index(drop=True)
                RVdbn['Spe'] = RVdbn[1].str.split('|').str[4]
                #print(RVdbn[['Contig','length','Spe',2,'cov','cer']])
                NTdb = pd.read_table(f'{tkid}.ntout.txt',header=None)
                NTdb1 = Rawdb.merge(NTdb,left_on='Contig',right_on=0)
                NTdb1['cov'] = NTdb1[3]/NTdb1['length']
                NTdbn = NTdb1.groupby('Contig').apply(lambda x:get_bestr(x))
                NTdbn = NTdbn.reset_index(drop=True)
                NTdbn = NTdbn.merge(taxiddb,left_on=12,right_on='Taxid')
                Rawdb['NT'] = Rawdb.apply(lambda x:get_nt(x['Contig'],NTdbn),axis=1)
                Rawdb['RVDB'] = Rawdb.apply(lambda x:get_nt(x['Contig'],RVdbn),axis=1)
                Rawdb['FinalSpe'] = Rawdb.apply(lambda x:getSpe(x),axis=1)
                Rawdb[['FinalSpe']].drop_duplicates().to_csv(f'{Pre}.finalSpe.list.txt',sep='\t',index=False,header=False)
                subprocess.run(f'''cat {Pre}.finalSpe.list.txt|taxonkit name2taxid > {Pre}.finalSpeid.list.txt ''',shell=True)
                finalSpedb = pd.read_table(f'{Pre}.finalSpeid.list.txt',header=None,names=['FinalSpe','FinalSpeTaxid'])
                finalSpedb.fillna('noid',inplace=True)
                Rawdb = Rawdb.merge(finalSpedb,on='FinalSpe')
                Rawdb.to_csv(f'{tkid}.reads.tsv',sep='\t',index=False)

def check_and_generate_bwa_index(ref_path):

    bwa_index_files = [
        f"{ref_path}.amb",
        f"{ref_path}.ann",
        f"{ref_path}.bwt",
        f"{ref_path}.pac",
        f"{ref_path}.sa"
    ]

    # Check if all BWA index files exist
    if all(os.path.exists(file) for file in bwa_index_files):
        print("BWA index already exists. No action needed.")
        return  # Pass if index files already exist
    else:
        print("BWA index not found. Generating index...")
        # Run BWA index command to generate the index files
        try:
            subprocess.run(["bwa", "index", ref_path], check=True)
            print("BWA index generated successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error generating BWA index: {e}")

def mappingRate(Pre,threads,ref):
    if ref !='noref':
        with open('mapping.log','w') as f:
            check_and_generate_bwa_index(ref)
            subprocess.run(f'bwa mem {ref} 10239.1.fastq -t {threads}|samtools sort -o {Pre}.R1.sorted.bam',shell=True)
            subprocess.run(f'bwa mem {ref} 10239.2.fastq -t {threads}|samtools sort -o {Pre}.R2.sorted.bam',shell=True)
            subprocess.run(f'samtools index {Pre}.R1.sorted.bam',shell=True)
            subprocess.run(f'samtools index {Pre}.R2.sorted.bam',shell=True)
            if not os.path.isfile(f'{Pre}.regions.bed'):
                subprocess.run(f'minimap2 -ax sr {ref} 10239.1.fastq 10239.1.fastq -t 10 |samtools sort -o {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                subprocess.run(f'samtools index {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                subprocess.run(f'mosdepth -b 1 -n {Pre} {Pre}.sorted.bam -t 10',shell=True,stdout=f,stderr=f)
                subprocess.run(f'gunzip -f  {Pre}.regions.bed.gz',shell=True,stdout=f,stderr=f)
            tmpdb = pd.read_table(f'{Pre}.regions.bed',header=None)
            genomesize = tmpdb.shape[0]
            Cov1 = tmpdb.loc[tmpdb[3]>=1].shape[0]/genomesize*100
            Cov10 = tmpdb.loc[tmpdb[3]>=10].shape[0]/genomesize*100
            Cov50 = tmpdb.loc[tmpdb[3]>=50].shape[0]/genomesize*100
            meandepth = tmpdb[3].mean()
            open(f'{Pre}.cov_summary.tsv','w').write(f'Sample\tGenomesize\tcov(1x)\tcov(10x)\tcov(50x)\tMeanDepth\n')
            open(f'{Pre}.cov_summary.tsv','a').write(f'{Pre}\t{genomesize}\t{Cov1}\t{Cov10}\t{Cov50}\t{meandepth}')
            R1Rreadsnum = int(os.popen(f'samtools view -F 2308 -f 16 {Pre}.R1.sorted.bam|wc -l').read())
            R1Freadsnum = int(os.popen(f'samtools view -F 2308 -F 16 {Pre}.R1.sorted.bam|wc -l').read())
            R2Rreadsnum = int(os.popen(f'samtools view -F 2308 -f 16 {Pre}.R2.sorted.bam|wc -l').read())
            R2Freadsnum = int(os.popen(f'samtools view -F 2308 -F 16 {Pre}.R2.sorted.bam|wc -l').read())
            subprocess.run(f'mosdepth -b1 {Pre}_R1 {Pre}.R1.sorted.bam',shell=True)
            subprocess.run(f'mosdepth -b1 {Pre}_R2 {Pre}.R2.sorted.bam',shell=True)
            subprocess.run(f'gunzip -f {Pre}_R1.regions.bed.gz',shell=True)
            subprocess.run(f'gunzip -f {Pre}_R2.regions.bed.gz',shell=True)
            tmpR1db = pd.read_table(f'{Pre}_R1.regions.bed',header=None)
            tmpR2db = pd.read_table(f'{Pre}_R2.regions.bed',header=None)
            R1cov = tmpR1db[tmpR1db[3]>=10].shape[0]/tmpR1db.shape[0]*100
            R2cov = tmpR2db[tmpR1db[3]>=10].shape[0]/tmpR2db.shape[0]*100
            tmpdict = {'R1':{'+':R1Freadsnum,'-':R1Rreadsnum,'cov':R1cov},'R2':{'+':R2Freadsnum,'-':R2Rreadsnum,'cov':R2cov}}
            tmpdb = pd.DataFrame(tmpdict).T
            tmpdb['Ratio+'] = tmpdb['+']/(tmpdb['+']+tmpdb['-']) 
            tmpdb['Ratio-'] = tmpdb['-']/(tmpdb['+']+tmpdb['-'])
            tmpdb.reset_index().to_csv(f'{Pre}_mapping_summary.tsv',sep='\t',index=False)





def blast_result(Pre,threads):
    subprocess.run(f'''blastn -db ~/Database/nt_core/core_nt/core_nt -query contigs.filter1000.fa -out {Pre}.blast.out -num_threads {threads} -evalue 1e-5 -outfmt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname' -word_size 28 -max_target_seqs 10''',shell=True)
    afile = pd.read_table(f'{Pre}.blast.out',header=None,names=['序列ID','参考ID','一致性','长度','错配碱基','空缺碱基','序列起始位置','序列终止位置','参考起始位置','参考终止位置','evalue','比对得分','TaxID','SamName'])
    
def main_process(Pre,fq1,fq2,threads,workstep,downsample,Hostref,ref):
    timeS = time.time()
    if 'QC' in workstep:
        print(f'{Pre}\t数据质控已开始')
        Preprocess(Pre,fq1,fq2,threads,downsample)
        QCend = time.time()
        QCtime = format_seconds(QCend-timeS)
        print(f'{Pre}\t数据质控已结束.用时：{QCtime}')
        timeS = QCend
    if 'RmC' in workstep:
        print(f'{Pre}\t去宿主已开始')
        removeHost(Pre, f'{Pre}.clean.R1.fastq.gz', f'{Pre}.clean.R2.fastq.gz',threads,Hostref)
        flashend = time.time()
        flashtime = format_seconds(flashend-timeS)
        print(f'{Pre}\t去宿主已结束.用时：{flashtime}')
        timeS = flashend
    if 'Kk2' in workstep:
        print(f'{Pre}\tkmer物种鉴定已开始')
        #if not os.path.isfile(f'{Pre}.raw.report.txt') or os.path.getsize(f'{Pre}.raw.report.txt') == 0:
        #    subprocess.run(f'kraken2 --db /Data/RawData/Database/kk2db/new_kk2/ --threads  {threads} --report {Pre}.raw.report.txt --output {Pre}.raw.txt {Pre}.clean.R1.fastq.gz {Pre}.clean.R2.fastq.gz',shell=True)
        if not os.path.isfile(f'{Pre}.rmhost.report.txt') or os.path.getsize(f'{Pre}.rmhost.report.txt') == 0:
            subprocess.run(f'kraken2 --db /Data/RawData/Database/kk2db/new_kk2/ --threads  {threads} --report {Pre}.rmhost.report.txt --output {Pre}.rmhost.txt {Pre}.rmhost.1.fq {Pre}.rmhost.2.fq --confidence 0.1 --threads 20 ',shell=True)
        #临时
        subprocess.run(f'ln -s {Pre}.rmhost.report.txt {Pre}.raw.report.txt',shell=True)
        kk2(f'{Pre}.raw','/Data/Analysis/Virues_WTS/taxa_info_20210508.txt')
        kk2(f'{Pre}.rmhost','/Data/Analysis/Virues_WTS/taxa_info_20210508.txt')
        Mapend = time.time()
        Maptime = format_seconds(Mapend-timeS)
        print(f'{Pre}\tkmer物种鉴定已结束.用时：{Maptime}')
        timeS = Mapend
    if 'Asm' in workstep:
        print(f'{Pre}\t无参组装已开始')
        asmvirus(f'{Pre}.rmhost',threads)
        Asmend = time.time()
        Asmtime = format_seconds(Asmend-timeS)
        print(f'{Pre}\t无参组装已结束.用时：{Asmtime}')
        timeS = Asmend
    if 'Ctg' in workstep:
        print(f'{Pre}\tContig识别已开始')
        ContigIden(Pre,threads)
        kk2cmbl(Pre,threads)
        Conend = time.time()
        Contime = format_seconds(Conend-timeS)
        print(f'{Pre}\tContig识别已结束.用时：{Contime}')
        timeS = Conend
    if 'Rfa' in workstep:
        print(f'{Pre}\t有参组装已开始')
        mappingRate(Pre,threads,ref)
        Rfaend = time.time()
        Rfatime = format_seconds(Rfaend-timeS)
        print(f'{Pre}\t有参组装已结束.用时：{Rfatime}')
        timeS = Rfaend

import multiprocessing

if not os.path.isdir(ofn):
    os.makedirs(ofn)
os.chdir(ofn)
wkdir = os.getcwd()

def process_line(line, threads, workstep,ofn,downsample):
    colnum = len(line.strip().split('\t'))
    if colnum == 3:
        Pre, fq1, fq2  = line.strip().split('\t')
        Hostref = 'noref'
        ref = 'noref'
    elif colnum == 4:
        Pre, fq1, fq2, Hostref = line.strip().split('\t')
        ref = 'noref'
    elif colnum == 5:
        Pre, fq1, fq2, Hostref,ref = line.strip().split('\t')
    if not os.path.isdir(f'{wkdir}/{Pre}'):
        os.makedirs(f'{wkdir}/{Pre}')
    os.chdir(f'{wkdir}/{Pre}')

    main_process(Pre, fq1, fq2, threads, workstep,downsample,Hostref,ref)

def read_and_process_file(listn, threads, workstep,queue,ofn,downsample):
    with open(listn) as f1:
        # 创建一个进程池
        with multiprocessing.Pool(processes=queue) as pool:
            pool.starmap(process_line, [(line, threads, workstep,ofn,downsample) for line in f1])

if __name__ == '__main__':
    read_and_process_file(listn, threads, workstep,queue,ofn,downsample) 
    print('分析结束')
    #compare_NC(ofn)
