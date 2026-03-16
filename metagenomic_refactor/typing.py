from __future__ import annotations

import json
import math
import os
import re
import subprocess
from typing import Dict, Iterable, List, Set, Tuple

import pandas as pd
import pytaxonkit

from metagenomic_refactor.context import get_runtime_context


def serotype_B(pre):   # 大肠+致贺——ectyper
    runtime = get_runtime_context()
    runtime_nt = runtime.nt
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n cm210 ectyper -i ./{pre}.final.fasta -c {runtime_nt} -o {pre}_serotype --verify --pathotype',shell=True)
    s1 = pd.read_table(f'{pre}_serotype/output.tsv')
    s1 = s1[['Name', 'O-type', 'H-type', 'Serotype', 'Species', 'Pathotype', 'StxSubtypes']]
    s1.columns = ['样本名', 'O抗原', 'H抗原', '血清型', '物种', '分型', 'Stx亚型']
    s1['样本名'] = pre
    s1['志贺分型'] = '-'
    if os.path.isfile('flex_out.tsv'):
        subprocess.run('rm flex_out.tsv', shell=True)
    if s1['物种'].tolist()[0] == 'Shigella flexneri':
        subprocess.run(f'cat {pre}.final.fasta |seqkit amplicon -p /data/test/cptyper_test/flex_primer.tsv --bed -o flex_out.tsv', shell=True)
        flexdb = pd.read_table('flex_out.tsv', header=None)
        flexplist = set(flexdb[3].tolist())
        flextype = '-'
        if flexplist:
            if 'wzx1' in flexplist:
                if flexplist == set(['wzx1', 'gtrI']):
                    flextype = '1a'
                elif flexplist == set(['wzx1', 'gtrI', 'oac']):
                    flextype = '1b'
                elif flexplist == set(['wzx1', 'gtrI', 'oac', 'gtrIC']):
                    flextype = '1c'
                elif flexplist == set(['wzx1', 'gtrII']):
                    flextype = '2a'
                elif flexplist == set(['wzx1', 'gtrII', 'gtrX']):
                    flextype = '2b'
                elif flexplist == set(['wzx1', 'oac', 'gtrX']):
                    flextype = '3a'
                elif flexplist == set(['wzx1', 'oac']):
                    flextype = '3b'
                elif flexplist == set(['wzx1', 'gtrIV']):
                    flextype = '4a'
                elif flexplist == set(['wzx1', 'gtrIV', 'oac']):
                    flextype = '4b'
                elif flexplist == set(['wzx1', 'gtrV']):
                    flextype = '5a'
                elif flexplist == set(['wzx1', 'gtrX']):
                    flextype = 'X或Xv'
                elif flexplist == set(['wzx1']):
                    flextype = 'Y'
            elif flexplist == set(['wzx6']):
                flextype = 'F6'
        s1['志贺分型'] = flextype
    s1.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=0)
    return s1['血清型'].tolist()[0]


def serotype_D(pre):  # 副溶血弧菌血清型——VPsero
    runtime = get_runtime_context()
    runtime_nt = runtime.nt
    subprocess.run(f'mkdir VPsero;cp ./{pre}.final.fasta VPsero', shell=True)
    subprocess.run(f'python /home/dell/biosoft/VPsero-master/program.py -i VPsero -o my_out_put_2  -n {runtime_nt}', shell=True)
    s1 = pd.read_excel('my_out_put_2/serotype_predict/04.predict_result/all_strain_predict_result.xlsx')
    s1 = s1[['O_Spec_Gene', 'K_Spec_Gene', 'Predict_O_sero', 'Predict_K_sero', 'New_serotype']]
    s1.columns = ['O血清型基因', 'K血清型基因', 'O血清型', 'K血清型', '血清型类型']
    s1['血清型'] = s1['O血清型'] + ':' + s1['K血清型']
    s1.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=0)
    subprocess.run('cp my_out_put_i/serotype_predict/04.predict_result/all_strain_predict_result.xlsx {pre}_strain_predict_result.xlsx', shell=True)
    return s1['血清型'].tolist()[0]


def PathoNet(Pre, species):
    PathoSamdict = {'样本名称': Pre, '物种': species, '血清型': '-', '毒力基因': '-'}
    PathoNetdict = {
        'vcholerae': {'serotype': ['O1', 'O139'], 'vfgene': ['ctxA', 'ctxB']},
        'senterica': {'serotype': ['S.Typhi', 'S.Paratyphi A', 'S.Paratyphi B', 'S.Paratyphi C', 'S.Enteritidis', 'S.Typhimurium', 'S.Choleracsuis', 'S.Derby', 'S.London', 'S.Stanley', 'S.Calabar', 'S.Agona', 'S.Thompson', 'S.Rissen', 'S.enterica subsp. enterica serovar Typhimurium monophasic variant'], 'vfgene': []},
        'campylobacter': {'serotype': ['HS:1', 'HS:2', 'HS:4', 'HS:19', 'HS:23', 'HS:41', 'HS:44'], 'vfgene': ['hcp', 'virB', 'ciaB', 'ggt', 'cdtA', 'cdtB', 'ctdC', 'cgtA', 'cgtB', 'wlaN', 'cstII']},
        'klebsiella': {'serotype': ['K1', 'K2', 'K5', 'K20', 'K54', 'K57'], 'vfgene': []},
        'ecoli': {'serotype': ['O2', 'O45', 'O103', 'O111', 'O121', 'O145', 'O157'], 'vfgene': ['stx1A', 'stx1B', 'stx2A', 'stx2B', 'stxA']},
        'Shigella': {'serotype': ['1a', '1b', '1c', '2a', '2b', '3a', '3b', '4a', '4b', '5a', '5b', 'X', 'Xv', 'F6', 'Y'], 'vfgene': ['stx1A', 'stx1B', 'stx2A', 'stx2B', 'stxA']},
    }
    if species == 'campylobacter':
        pathodict = PathoNetdict[species]
        if os.path.isfile(f'{Pre}.vfdb.tsv') and os.path.getsize(f'{Pre}.vfdb.tsv') != 0:
            cpvfdb = pd.read_table(f'{Pre}.vfdb.tsv')
            tarvflist = [i for i in pathodict['vfgene'] if i in cpvfdb['基因名称'].tolist()]
            if tarvflist:
                PathoSamdict['毒力基因'] = ';'.join(tarvflist)
        if os.path.isfile(f'{Pre}_serotype_result.tsv') and os.path.getsize(f'{Pre}_serotype_result.tsv') != 0:
            cpserodb = pd.read_table(f'{Pre}_serotype_result.tsv')
            if cpserodb['血清型'].tolist()[0] in pathodict['serotype']:
                PathoSamdict['血清型'] = f"{cpserodb['血清型'].tolist()[0]}(重点关注)"
            else:
                PathoSamdict['血清型'] = cpserodb['血清型'].tolist()[0]
    if species == 'klebsiella':
        pathodict = PathoNetdict[species]
        if os.path.isfile(f'{Pre}.vfdb.tsv') and os.path.getsize(f'{Pre}.vfdb.tsv') != 0:
            klvfdb = pd.read_table(f'{Pre}.vfdb.tsv')
            tarvflist = [i for i in pathodict['vfgene'] if i in klvfdb['基因名称'].tolist()]
            if tarvflist:
                PathoSamdict['毒力基因'] = ';'.join(tarvflist)
        if os.path.isfile(f'{Pre}_serotype_result.tsv') and os.path.getsize(f'{Pre}_serotype_result.tsv') != 0:
            cpserodb = pd.read_table(f'{Pre}_serotype_result.tsv')
            stype = cpserodb['KO血清型'].tolist()[0].split('|')[0].replace('KL', 'K')
            if stype in pathodict['serotype']:
                PathoSamdict['血清型'] = f'{stype}(重点关注)'
            else:
                PathoSamdict['血清型'] = stype
    if species == 'senterica':
        pathodict = PathoNetdict[species]
        if os.path.isfile(f'{Pre}.vfdb.tsv') and os.path.getsize(f'{Pre}.vfdb.tsv') != 0:
            salvfdb = pd.read_table(f'{Pre}.vfdb.tsv')
            tarvflist = [i for i in pathodict['vfgene'] if i in salvfdb['基因名称'].tolist()]
            if tarvflist:
                PathoSamdict['毒力基因'] = ';'.join(tarvflist)
        if os.path.isfile(f'{Pre}_serotype_result.tsv') and os.path.getsize(f'{Pre}_serotype_result.tsv') != 0:
            cpserodb = pd.read_table(f'{Pre}_serotype_result.tsv')
            stype = cpserodb['亚型全称'].tolist()[0]
            if stype in pathodict['serotype']:
                PathoSamdict['血清型'] = f'{stype}(重点关注)'
            else:
                PathoSamdict['血清型'] = stype
    if species == 'vcholerae':
        pathodict = PathoNetdict[species]
        if os.path.isfile(f'{Pre}.vfdb.tsv') and os.path.getsize(f'{Pre}.vfdb.tsv') != 0:
            vchovfdb = pd.read_table(f'{Pre}.vfdb.tsv')
            tarvflist = [i for i in pathodict['vfgene'] if i in vchovfdb['基因名称'].tolist()]
            if tarvflist:
                PathoSamdict['毒力基因'] = ';'.join(tarvflist)
        if os.path.isfile(f'{Pre}_serotype_result.tsv') and os.path.getsize(f'{Pre}_serotype_result.tsv') != 0:
            cpserodb = pd.read_table(f'{Pre}_serotype_result.tsv').fillna('-')
            stype = cpserodb['血清型'].tolist()[0]
            if stype in pathodict['serotype']:
                PathoSamdict['血清型'] = f'{stype}(重点关注)'
            else:
                PathoSamdict['血清型'] = stype
    if species == 'ecoli_achtman_4':
        ecodb = pd.read_table(f'{Pre}_serotype_result.tsv')
        if ecodb['物种'].tolist()[0] == 'Escherichia coli':
            pathodict = PathoNetdict['ecoli']
            stype = ecodb['O抗原'].tolist()[0]
        else:
            pathodict = PathoNetdict['Shigella']
            stype = ecodb['志贺分型'].tolist()[0]
        PathoSamdict['血清型'] = f'{stype}(重点关注)' if stype in pathodict['serotype'] else stype
        PathoSamdict['物种'] = ecodb['物种'].tolist()[0]
        ecovfdb = pd.read_table(f'{Pre}.vfdb.tsv')
        tarvflist = [i for i in pathodict['vfgene'] if i in ecovfdb['基因名称'].tolist()]
        if tarvflist:
            PathoSamdict['毒力基因'] = ';'.join(tarvflist)
    pd.DataFrame(PathoSamdict, index=[0]).to_csv(f'{Pre}.pathonet_result.tsv', sep='\t', index=False)


def bp_vaccine(Pre):  # 百日咳疫苗基因型——BLASTn
    tmpdir = '/data1/shanghai_pip/meta_genome/database/BIGsdb/bordetella'
    tmpdict = {}
    genethrdict = {'23S_rRNA': 2800, 'fhaB-2400_5550': 3100, 'fim2': 600, 'fim3': 600, 'prn': 2600, 'ptxA': 800, 'ptxB': 350, 'ptxC': 650, 'ptxD': 350, 'ptxE': 350, 'ptxP': 150, 'tcfA': 1900}
    for i in os.listdir(tmpdir):
        if i.endswith('.fas'):
            tgene = i.replace('.fas', '')
            bitlen = genethrdict.get(tgene, 100)
            subprocess.run(f'''blastn -db {tmpdir}/{tgene} -query {Pre}.final.fasta -out {tgene}.blast.out -num_threads 10 -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 10  -perc_identity 90 -max_hsps 1''', shell=True)
            if os.path.getsize(f'{tgene}.blast.out') != 0:
                tmpdb = pd.read_table(f'{tgene}.blast.out', header=None)
                if tmpdb.shape[0] > 0:
                    tmpdb = tmpdb.loc[tmpdb[3] > bitlen]
                    tmpdict[tgene] = {}
                    if tmpdb.shape[0] > 0:
                        tmpdict[tgene] = {'基因名称': tgene, 'Contig名称': tmpdb.iloc[0, 0], '起始位置': tmpdb.iloc[0, 6], '终止位置': tmpdb.iloc[0, 7], '分型': tmpdb.iloc[0, 1], '一致性': tmpdb.iloc[0, 2], '差异碱基数量': tmpdb.iloc[0, 4]}
                    else:
                        tmpdict[tgene] = {'基因名称': tgene, 'Contig名称': '-', '起始位置': '-', '终止位置': '-', '分型': '-', '一致性': '-', '差异碱基数量': '-'}
            else:
                tmpdict[tgene] = {'基因名称': tgene, 'Contig名称': '-', '起始位置': '-', '终止位置': '-', '分型': '-', '一致性': '-', '差异碱基数量': '-'}
    newdb = pd.DataFrame(tmpdict).T.reset_index(drop=True)
    newdb.to_csv(f'{Pre}_scheme.tsv', sep='\t', index=False)


def serotype_HI(pre):  # 流感嗜血亚型预测——HICAP
    subprocess.run(f'''cut -d '' -f 1 {pre}.final.fasta > {pre}.fasta''', shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n GTDBtk python /home/dell/biosoft/hicap-master/hicap-runner.py -q  {pre}.final.fasta -o ./', shell=True)
    if os.path.isfile(f'{pre}.tsv'):
        hidb = pd.read_table(f'{pre}.tsv')
        hidb['样本名称'] = pre
        hidb.rename(columns={'predicted_serotype': '亚型预测', 'genes_identified': '检测基因', 'IS1016_hits': 'IS1016数量'}, inplace=True)
        hidb = hidb[['样本名称', '亚型预测', '检测基因', 'IS1016数量']]
        hidb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
        sers = hidb['亚型预测'].tolist()[0]
    else:
        open(f'{pre}_serotype_result.tsv', 'w').write('样本名称\t亚型预测\t检测基因\tIS1016数量\n')
        open(f'{pre}_serotype_result.tsv', 'a').write(f'{pre}\t-\t-\t-')
        sers = '-'
    return sers


def serotype_ys(pre):   # 耶尔森菌
    YS_dict = {'样本名': pre, 'PLA': '-', 'YPO2088': '-', 'inv': '-', 'opgG': '-', '97_predict': 'Unknown', '物种预测': '未知耶尔森'}
    subprocess.run(f'cat {pre}.final.fasta|seqkit amplicon -p /data/deploy/meta_genome/YS_primer.tsv --bed > YS_primer.bed', shell=True)
    if os.path.isfile('YS_primer.bed') and os.path.getsize('YS_primer.bed') != 0:
        primerdb = pd.read_table('YS_primer.bed', header=None)
        genelist = primerdb[3].tolist()
        for targene in genelist:
            YS_dict[targene] = '+'
    subprocess.run(f'python /data/deploy/meta_genome/Identify_Y.pestis-main/Identify_Y.pestis/Identify_Y.pestis_from_dir.py {pre}.final.fasta YS_97.txt', shell=True)
    if os.path.isfile('YS_97.txt') and os.path.getsize('YS_97.txt') != 0:
        YS97db = pd.read_table('YS_97.txt')
        YS_dict['97_predict'] = YS97db['Is_Ypestis'].tolist()[0]
    if YS_dict['PLA'] == '+' and YS_dict['YPO2088'] == '+' and YS_dict['97_predict'] == 'Yes':
        YS_dict['物种预测'] = '鼠疫耶尔森'
    elif YS_dict['opgG'] == '+' and YS_dict['inv'] == '+' and YS_dict['97_predict'] == 'No':
        YS_dict['物种预测'] = '假结核耶尔森'
    s1 = pd.DataFrame(YS_dict, index=[0])
    s1.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=0)
    return s1['物种预测'].tolist()[0]


def serotype_A(pre):  # 沙门氏菌血清型——sistr
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n GTDBtk sistr -i {pre}.final.fasta {pre} -f tab -o {pre}_result', shell=True)
    serotype_s = '/data/deploy/bio-elite/bio/load_file/pathogenic/salmonella_52seotype_v2.txt'
    sero_info = pd.read_table(serotype_s)
    sero = pd.read_table(f'{pre}_result.tab')
    sero = sero.rename(columns={'o_antigen': 'O抗原', 'h1': 'H1相抗原(fliC)', 'h2': 'H2相抗原(fljB)'})
    sero['O抗原'] = sero['O抗原'].astype('str')
    serox = sero.merge(sero_info, on=['O抗原', 'H1相抗原(fliC)'], how='left')
    serox.fillna('-', inplace=True)
    sero_result = pd.DataFrame()
    sero_result['样本'] = pre
    sero_result[['O抗原', 'H1相抗原(fliC)', 'H2相抗原(fljB)']] = serox[['O抗原', 'H1相抗原(fliC)', 'H2相抗原(fljB)_y']]
    sero_result[['菌种', '血清型', '血清型注释信息(simple)', '血清型注释信息(details)']] = serox[['serovar', 'serogroup', 'simple_description', 'details_description']]
    sero_result['抗原组成'] = sero_result['O抗原'] + ':' + sero_result['H1相抗原(fliC)'] + ":" + sero_result['H2相抗原(fljB)']
    sero_result['样本'] = pre
    sero_result['亚型全称'] = serox['name']
    sero_result.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=0)
    return sero_result['血清型'].tolist()[0]


def serotype_kb(pre):  # 克雷伯菌血清型分型——kleborate
    subprocess.run(f'kleborate --all -o results.txt -a {pre}.final.fasta > {pre}.keblo.tsv', shell=True)
    kledb = pd.read_table(f'{pre}.keblo.tsv')
    kledb['样本名称'] = pre
    kledb.rename(columns={'virulence_score': '毒力得分', 'resistance_score': '耐药得分', 'Yersiniabactin': '耶尔森菌素', 'Colibactin': '大肠菌素', 'Bla_chr': '氨苄类耐药SHV等位基因', 'SHV_mutations': 'SHV耐药突变', 'wzi': 'wzi荚膜预测'}, inplace=True)
    kledb['KO血清型'] = kledb['K_locus'].tolist()[0] + '|' + kledb['O_locus'].tolist()[0]
    kledb = kledb[['样本名称', 'ST', '毒力得分', '耐药得分', '耶尔森菌素', '大肠菌素', '氨苄类耐药SHV等位基因', 'SHV耐药突变', 'wzi荚膜预测', 'KO血清型']]
    kledb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
    return kledb['KO血清型'].tolist()[0]


def serotype_nm(pre):  # 奈瑟氏菌亚型预测——PMGA
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n RGI pmga {pre}.final.fasta  -t 10 --force --blastdir /data1/shanghai_pip/meta_genome/database/pmga/', shell=True)
    nmdb = pd.read_table(f'pmga/{pre}.finalsta.txt')
    nmdb['样本名称'] = pre
    nmdb.rename(columns={'prediction': '亚型预测', 'genes_present': '验证基因集', 'notes': '注释'}, inplace=True)
    nmdb = nmdb[['样本名称', '亚型预测', '验证基因集', '注释']]
    nmdb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
    return nmdb['亚型预测'].tolist()[0]


def serotype_groupA(pre):  # GAS，A链emm分型——emm_typing.py
    subprocess.run(f'python /home/dell/biosoft/emm_typing/emm_typing/emm_typing.py -f {pre}.final.fasta', shell=True)
    gasdb = pd.read_table('emm_results.tab')
    gasdb.rename(columns={'emm-type': '亚型预测', 'pident': '一致性', 'Isolate': '样本名称', 'length': '比对长度'}, inplace=True)
    gasdb = gasdb[['样本名称', '亚型预测', '一致性', '比对长度']]
    gasdb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
    return gasdb['亚型预测'].tolist()[0]


def serotype_E(pre):  # 霍乱亚型预测
    serodict = {}
    subprocess.run(f'''/home/dell/miniconda3/bin/conda run -n  choleraefinder python /data/test/test_VP/choleraefinder/choleraefinder.py --input {pre}.final.fasta -o ./ -p /data/test/test_VP/choleraefinder_db/  -t 0.95 -l 0.95 -q''', shell=True)
    serojson = json.load(open('data_CholeraeFinder.json'))
    serodict['样本名称'] = pre
    serodict['血清型'] = serojson['choleraefinder']['typing_cholerae']['serogroup']
    serodict['生物型'] = serojson['choleraefinder']['typing_cholerae']['biotype']
    serodb = pd.DataFrame(serodict, index=[0])
    serodb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
    return serojson['choleraefinder']['typing_cholerae']['serogroup']


def serotype_MLVA(pre):  # 布鲁氏菌MLVA分析
    if not os.path.isdir(f'{pre}_mlvafafile'):
        os.makedirs(f'{pre}_mlvafafile')
    subprocess.run(f'cp {pre}.final.fasta {pre}_mlvafafile', shell=True)
    subprocess.run(f'python /data1/shanghai_pip/meta_genome/database/MLVA_finder/MLVA_finder.py -i {pre}_mlvafafile -o ./ -p /data1/shanghai_pip/meta_genome/database/MLVA_finder/data_test/primers/Brucella_primers.txt', shell=True)
    mlvadb = pd.read_table(f'{pre}_mlvafafile_output.csv', sep=',')
    mlvadb['样本名称'] = pre
    mlvadb.rename(columns={'primer': '引物', 'position1': '起始位置', 'position2': '终止位置', 'size': '扩增片段大小', 'allele': '重复基因数量'}, inplace=True)
    mlvadb = mlvadb[['样本名称', '引物', '起始位置', '终止位置', '扩增片段大小', '重复基因数量']]
    mlvadb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
    return ';'.join(mlvadb['allele'].tolist())


def serotype_st(pre):  # 金葡家系分析——SALTY
    if not os.path.isdir(f'{pre}_st_fafile'):
        os.makedirs(f'{pre}_st_fafile')
    subprocess.run(f'cp {pre}.final.fasta {pre}_st_fafile', shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n RGI salty -i {pre}_st_fafile -o {pre}_st_fafile -t 10', shell=True)
    stdb = pd.read_table(f'{pre}_st_fafile/summaryReport.txt')
    stdb['样本名称'] = pre
    stdb.rename(columns={'Lineage': '家系', 'SACOL1908': 'SACOL1908基因座等位基因', 'SACOL0451': 'SACOL0451基因座等位基因', 'SACOL2725': 'SACOL2725基因座等位基因'}, inplace=True)
    stdb = stdb[['样本名称', '家系', 'SACOL1908基因座等位基因', 'SACOL0451基因座等位基因', 'SACOL2725基因座等位基因']]
    stdb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
    return stdb['家系'].tolist()[0]


def serotype_lm(pre):  # 单增李斯特血清型——lissero
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n RGI lissero {pre}.final.fasta > Lm_sero.tsv ', shell=True)
    stdb = pd.read_table('Lm_sero.tsv')
    stdb['样本名称'] = pre
    stdb.rename(columns={'SEROTYPE': '血清型'}, inplace=True)
    stdb = stdb[['样本名称', '血清型', 'PRS', 'LMO0737', 'LMO1118', 'ORF2110', 'ORF2819']]
    stdb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
    return stdb['血清型'].tolist()[0]


def serotype_bt(pre):  # 蜡样芽孢杆菌毒力因子——btyper3
    subprocess.run(f'btyper3 -i {pre}.final.fasta -o Bcere_sero', shell=True)
    btdb = pd.read_table(f'Bcere_sero/btyper3_final_results/{pre}.final_final_results.txt')
    btdb['样本名称'] = pre
    btdb.rename(columns={'final_taxon_names': '物种名称', 'anthrax_toxin(genes)': 'nthrax基因集', 'emetic_toxin_cereulide(genes)': 'cereulide基因集', 'diarrheal_toxin_Nhe(genes)': 'Nhe基因集', 'diarrheal_toxin_Hbl(genes)': 'Hbl基因集', 'diarrheal_toxin_CytK(top_hit)': 'CytK基因集', 'sphingomyelinase_Sph(gene)': 'Sph基因集', 'capsule_Cap(genes)': 'Cap基因集'}, inplace=True)
    btdb = btdb[['样本名称', '物种名称', 'nthrax基因集', 'cereulide基因集', 'Nhe基因集', 'Hbl基因集', 'CytK基因集', 'Sph基因集', 'Cap基因集']]
    btdb.to_csv(f'{pre}_serotype_result.tsv', sep='\t', index=False)
    return btdb['物种名称'].tolist()[0]


def serotype_SS(Pre):  # 猪链血清型
    database = '/data/deploy/meta_genome/database/SS_sero/'
    blastref = f'{database}/Ssuis_Serotyping'
    cpsref = f'{database}/Ssuis_cps2K.fasta'
    with open('sero.log', 'w') as serof:
        subprocess.run(f'blastn -query {Pre}.final.fasta -db {blastref} -out {Pre}.cps.out.tsv -outfmt 6 -perc_identity 90 ', shell=True)
        Serodb = pd.read_table(f'{Pre}.cps.out.tsv', header=None)
        lengdb = pd.read_table(f'{database}/Ssuis_Serotyping.tsv')
        Serodb = Serodb.merge(lengdb, left_on=1, right_on='Serotype')
        Serodb['Perc'] = Serodb[3] / Serodb['Length']
        Serodb = Serodb.loc[Serodb['Perc'] > 0.9]
        if Serodb.shape[0] > 0:
            Serotype = Serodb['Serotype'].tolist()[0]
        else:
            Serotype = 'notype'
        Serotype = Serotype.replace('cps-', '')
        if Serotype == '1' or Serotype == '2':
            subprocess.run(f'nucmer --maxmatch -b 200 -c 65 -d 0.12 -g 90 -l 20 {cpsref} {Pre}.final.fasta -p {Pre}', shell=True, stdout=serof, stderr=serof)
            subprocess.run(f'show-snps {Pre}.delta -T > {Pre}.snps.out', shell=True)
            snpdb = pd.read_table(f'{Pre}.snps.out', skiprows=2)
            snplist = snpdb['[P2]'].tolist()
            if Serotype == '1':
                if '483' not in snplist:
                    Serotype = '14'
            else:
                if '483' in snplist:
                    Serotype = '1/2'

        if os.path.isfile(f'{Pre}.SsuisChara.tsv'):
            subprocess.run(f'rm -r {Pre}.SsuisChara.tsv', shell=True)
        subprocess.run(f'python /data/deploy/meta_genome/database/SsuisChara/SsuisChara.py -i {Pre}.final.fasta -o {Pre}.SsuisChara.tsv', shell=True)
        SsuisCdb = pd.read_table(f'{Pre}.SsuisChara.tsv')
        SsuisCdb = SsuisCdb[['human infection potential', 'AMRG_level', 'aminoglycoside', 'macrolide', 'tetracycline']]
        SsuisCdb['样本名称'] = Pre
        SsuisCdb['血清型'] = Serotype
        SsuisCdb = SsuisCdb.rename(columns={'human infection potential': '感染等级', 'aminoglycoside': '氨基糖苷类', 'macrolide': '大环内酯类', 'tetracycline': '四环素类', 'AMRG_level': '耐药数量'})
        SsuisCdb = SsuisCdb[['样本名称', '血清型', '感染等级', '耐药数量', '氨基糖苷类', '大环内酯类', '四环素类']]
        SsuisCdb.to_csv(f'{Pre}_serotype_result.tsv', sep='\t', index=False)
        return Serotype


def len2mlvacopy(PrimerN, PrimerL):
    primerdict = {'VNTR1': [383, 15, 9], 'VNTR3a': [135, 5, 7], 'VNTR3b': [0, 0, 0], 'VNTR4': [232, 12, 9], 'VNTR5': [143, 6, 7], 'VNTR6': [234, 9, 11]}
    modlen, modcopylen, modcopynum = primerdict[PrimerN]
    if PrimerN != 'VNTR3b':
        PrimerCopy = round(modcopynum - (modlen - PrimerL) / modcopylen)
    else:
        PrimerCopy = 0
    return PrimerCopy


def bp_mlva(Pre, primer='/data/test/mlva/mlva_primer.tsv', typetable='/data1/shanghai_pip/meta_genome/BPMLVA.table.tsv'):
    subprocess.run(f'cat {Pre}.final.fasta |seqkit amplicon -p {primer} --bed > {Pre}.raw.mlva.tsv', shell=True)
    rawdb = pd.read_table(f'{Pre}.raw.mlva.tsv', header=None, names=['Chrom', 'Startpos', 'Endpos', 'PrimerName', 'Mismatch', 'Strand', 'Sequence'])
    finaldict = {}
    for PrimerN in rawdb['PrimerName'].tolist():
        tmpdb = rawdb.loc[rawdb['PrimerName'] == PrimerN, ]
        if PrimerN == 'VNTR3':
            if tmpdb.shape[0] == 1:
                PrimerN = 'VNTR3a'
                finaldict[PrimerN] = len2mlvacopy(PrimerN, tmpdb['Endpos'].tolist()[0] - tmpdb['Startpos'].tolist()[0])
                finaldict['VNTR3b'] = len2mlvacopy('VNTR3b', tmpdb['Endpos'].tolist()[0] - tmpdb['Startpos'].tolist()[0])
            else:
                VNTR3list = []
                tmpdb = tmpdb.reset_index()
                for tmpi in tmpdb.index:
                    tmptype = len2mlvacopy('VNTR3a', tmpdb.loc[tmpdb.index == tmpi, 'Endpos'].tolist()[0] - tmpdb.loc[tmpdb.index == tmpi, 'Startpos'].tolist()[0])
                    VNTR3list.append(tmptype)
                if len(set(VNTR3list)) == 1:
                    finaldict['VNTR3a'] = VNTR3list[0]
                else:
                    finaldict['VNTR3a'] = ';'.join(list(set(VNTR3list)))
                finaldict['VNTR3b'] = 0
        elif PrimerN == 'VNTR4':
            if tmpdb.shape[0] == 1:
                finaldict['VNTR4'] = len2mlvacopy(PrimerN, tmpdb['Endpos'].tolist()[0] - tmpdb['Startpos'].tolist()[0])
            else:
                finaldict['VNTR4'] = '未知'
        else:
            finaldict[PrimerN] = len2mlvacopy(PrimerN, tmpdb['Endpos'].tolist()[0] - tmpdb['Startpos'].tolist()[0])
    typedb = pd.read_table(typetable)
    typedb = typedb.astype('str')
    typedb['combine'] = typedb.apply(lambda x: x['VNTR1'] + '_' + x['VNTR3a'] + '_' + x['VNTR3b'] + '_' + x['VNTR4'] + '_' + x['VNTR5'] + '_' + x['VNTR6'], axis=1)
    typedb.to_csv('tt1.tsv', sep='\t', index=False)
    finaldb = pd.DataFrame(finaldict, index=[0]).reset_index(drop=True)
    finaldb = finaldb.astype('str')
    if finaldb.shape[1] == 6:
        finaldb['combine'] = finaldb.apply(lambda x: x['VNTR1'] + '_' + x['VNTR3a'] + '_' + x['VNTR3b'] + '_' + x['VNTR4'] + '_' + x['VNTR5'] + '_' + x['VNTR6'], axis=1)
        finaldb = finaldb.merge(typedb, on='combine', how='left')
        finaldb['样本名称'] = Pre
        finaldb.rename(columns={'VNTR1_x': 'VNTR1', 'VNTR3a_x': 'VNTR3a', 'VNTR3b_x': 'VNTR3b', 'VNTR4_x': 'VNTR4', 'VNTR5_x': 'VNTR5', 'VNTR6_x': 'VNTR6'}, inplace=True)
        finaldb = finaldb[['样本名称', 'MT', 'VNTR1', 'VNTR3a', 'VNTR3b', 'VNTR4', 'VNTR5', 'VNTR6']]
        finaldb.to_csv(f'{Pre}.mlva.tsv', sep='\t', index=False)
    else:
        print(f'{Pre} {finaldb}')


def len2mlvacopy_mp(PrimerN, PrimerL):
    primerdict = {'Mpn13': [428, 16, 4], 'Mpn14': [378, 21, 4], 'Mpn15': [192, 21, 4], 'Mpn16': [447, 47, 4]}
    modlen, modcopylen, modcopynum = primerdict[PrimerN]
    PrimerCopy = math.ceil(modcopynum - (modlen - PrimerL) / modcopylen)
    return PrimerCopy


def mp_mlva(Pre, primer='/data/test/mlva/mlva_primer_mp.tsv'):
    subprocess.run(f'cat {Pre}.final.fasta |seqkit amplicon  -p {primer} --bed > {Pre}.raw.mlva.tsv', shell=True)
    rawdb = pd.read_table(f'{Pre}.raw.mlva.tsv', header=None, names=['Chrom', 'Startpos', 'Endpos', 'PrimerName', 'Mismatch', 'Strand', 'Sequence'])
    finaldict = {}
    for PrimerN in rawdb['PrimerName'].tolist():
        tmpdb = rawdb.loc[rawdb['PrimerName'] == PrimerN, ]
        finaldict[PrimerN] = len2mlvacopy_mp(PrimerN, tmpdb['Endpos'].tolist()[0] - tmpdb['Startpos'].tolist()[0])
    finaldb = pd.DataFrame(finaldict, index=[0]).reset_index(drop=True)
    finaldb = finaldb.astype('str')
    print(finaldb)
    if finaldb.shape[1] == 4:
        finaldb['combine'] = finaldb.apply(lambda x: x['Mpn13'] + '_' + x['Mpn14'] + '_' + x['Mpn15'] + '_' + x['Mpn16'], axis=1)
        finaldb['样本名称'] = Pre
        finaldb = finaldb[['样本名称', 'combine', 'Mpn13', 'Mpn14', 'Mpn15', 'Mpn16']]
        finaldb.to_csv(f'{Pre}.mlva.tsv', sep='\t', index=False)
    else:
        print(f'{Pre} {finaldb}')


def bp_2037(Pre):
    with open('bp2037.log', 'w') as f:
        if os.path.isfile(f'{Pre}.R2.fastq.gz'):
            subprocess.run(f'bwa mem /data1/shanghai_pip/meta_genome/rrn.fasta {Pre}.R1.fastq.gz {Pre}.R2.fastq.gz -t 10 |samtools sort -o {Pre}.rrn.sorted.bam', shell=True, stdout=f, stderr=f)
        else:
            subprocess.run(f'bwa mem /data1/shanghai_pip/meta_genome/rrn.fasta {Pre}.R1.fastq.gz -t 10 |samtools sort -o {Pre}.rrn.sorted.bam', shell=True, stdout=f, stderr=f)
        subprocess.run(f'samtools index {Pre}.rrn.sorted.bam', stdout=f, stderr=f, shell=True)
        subprocess.run(f'freebayes -f /data1/shanghai_pip/meta_genome/rrn.fasta {Pre}.rrn.sorted.bam > {Pre}.rrn.vcf', stdout=f, stderr=f, shell=True)
        subprocess.run(f'/home/dell/miniconda3/envs/PathoSource/bin/vt normalize {Pre}.rrn.vcf -r /data1/shanghai_pip/meta_genome/rrn.fasta -o {Pre}.rrn.filt1.vcf', shell=True, stdout=f, stderr=f)
        skip_rows = int(os.popen(f'''grep '##' {Pre}.rrn.filt1.vcf|wc -l ''').read())
        rrndb = pd.read_table(f'{Pre}.rrn.filt1.vcf', skiprows=skip_rows)
        if rrndb.shape[0] > 0:
            rrndb['GT'] = rrndb['unknown'].str.split(':').str[0]
            rrndb = rrndb[rrndb['GT'] != '0/0']
            if rrndb.shape[0] > 0:
                rrndb['POS'] = rrndb['POS'].astype('str')
                rrndb['symbol'] = rrndb['REF'] + rrndb['POS'] + rrndb['ALT']
                if any(rrndb['symbol'] == 'A2037G'):
                    rrndb.to_csv(f'{Pre}.2037.tsv', sep='\t', index=False)


def fa_2037(Pre):
    with open('2037.log', 'w') as f:
        subprocess.run(f'nucmer /data1/shanghai_pip/meta_genome/rrn.fasta {Pre}.final.fasta', shell=True, stdout=f, stderr=f)
        if int(os.popen('show-snps out.delta|grep -w 2037|wc -l').read()) >= 1:
            subprocess.run(f'show-snps out.delta > {Pre}.2037.tsv', shell=True)


def level2Spe(Pre, taxid):
    id2dict = {'Vibrio parahaemolyticus': 'vparahaemolyticus', 'Vibrio cholerae': 'vcholerae'}
    taxid = str(taxid)
    Spedb = pd.read_table('/data/test/level2Spe/tmp.txt')
    Spedb['idlist'] = Spedb['idlist'].str.split(';')
    nSpedb = Spedb[Spedb['idlist'].apply(lambda x: taxid in x)]
    if nSpedb.shape[0] == 1:
        nSpe = nSpedb['Spe'].tolist()[0]
    else:
        if nSpedb.shape[0] > 1:
            if not os.path.isfile(f'{Pre}.kraken2.report.txt'):
                subprocess.run(f'kraken2 --db /home/dell/kraken2_custom_202101_24G {Pre}.final.fasta --report {Pre}.kraken2.report.txt --output {Pre}.kraken2.txt -t 10', shell=True)
            fak2db = pd.read_table(f'{Pre}.kraken2.report.txt', header=None)
            nSpe = id2dict.get(fak2db.loc[fak2db[3] == 'S'][5].str.strip().tolist()[0], 0)
        else:
            nSpe = 0
    return nSpe


Prim = Tuple[str, str]
A = "alpha"
B = "beta"
G = "gamma"
D = "delta"
PRIMER_TO_MIX: Dict[str, str] = {
    "Mu_HS2": A, "Mu_HS3": A, "Mu_HS4A": A, "Mu_HS6": A, "Mu_HS10": A, "Mu_HS15": A,
    "Mu_HS41": A, "Mu_HS53": A, "Mu_HS19": A, "Mu_HS63": A, "Mu_HS33": A,
    "Mu_HS1": B, "Mu_HS4B": B, "Mu_HS8": B, "Mu_HS23/36": B, "Mu_HS42": B, "Mu_HS57": B,
    "Mu_HS12": B, "Mu_HS27": B, "Mu_HS21": B, "Mu_HS31": B,
    "Mu_HS44": G, "Mu_HS45": G, "Mu_HS29": G, "Mu_HS22": G, "Mu_HS9": G, "Mu_HS37": G,
    "Mu_HS18": G, "lpxA": G,
    "Mu_HS58": D, "Mu_HS52": D, "Mu_HS60": D, "Mu_HS55": D, "Mu_HS32": D, "Mu_HS11": D,
    "Mu_HS40": D, "Mu_HS38": D,
}
ALIASES: Dict[str, str] = {"Mu_HS5": "Mu_HS31"}


def norm_primer(p: str) -> str:
    return ALIASES.get(p, p)


class Rule:
    def __init__(self, name: str, must_have: Set[Prim], must_not: Set[Prim] = None, allowed_extra: Set[Prim] = None):
        self.name = name
        self.must_have = {(m, norm_primer(p)) for m, p in (must_have or set())}
        self.must_not = {(m, norm_primer(p)) for m, p in (must_not or set())}
        self.allowed_extra = {(m, norm_primer(p)) for m, p in (allowed_extra or set())}


RULES: List[Rule] = [
    Rule("HS1", {(B, "Mu_HS1")}), Rule("HS2", {(A, "Mu_HS2")}), Rule("HS3", {(A, "Mu_HS3")}), Rule("HS4 complex", {(A, "Mu_HS4A")}),
    Rule("CG8486 (HS4 complex member)", {(B, "Mu_HS4B")}), Rule("HS5", {(B, "Mu_HS31"), (G, "Mu_HS45")}), Rule("HS6", {(A, "Mu_HS6")}),
    Rule("HS7", {(A, "Mu_HS6")}), Rule("HS8", {(B, "Mu_HS8")}), Rule("HS9", {(G, "Mu_HS9")}), Rule("HS10", {(A, "Mu_HS10")}),
    Rule("HS11", {(D, "Mu_HS11")}), Rule("HS12", {(B, "Mu_HS12")}), Rule("HS13", {(A, "Mu_HS4A")}),
    Rule("HS15", {(A, "Mu_HS15")}, allowed_extra={(D, "Mu_HS58")}), Rule("HS16", {(A, "Mu_HS4A"), (B, "Mu_HS4B")}, allowed_extra={(D, "Mu_HS52")}),
    Rule("HS17", {(B, "Mu_HS8")}), Rule("HS18", {(G, "Mu_HS18")}), Rule("HS19", {(A, "Mu_HS19")}), Rule("HS21", {(B, "Mu_HS21")}),
    Rule("HS22", {(G, "Mu_HS22")}), Rule("HS23", {(B, "Mu_HS23/36")}), Rule("HS27", {(B, "Mu_HS27")}), Rule("HS29", {(G, "Mu_HS29")}),
    Rule("HS31", {(B, "Mu_HS31")}, allowed_extra={(A, "Mu_HS15")}), Rule("HS32", {(D, "Mu_HS32"), (G, "Mu_HS45")}, allowed_extra={(B, "Mu_HS8")}),
    Rule("HS33", {(A, "Mu_HS33")}), Rule("HS35", {(A, "Mu_HS33")}), Rule("HS36", {(B, "Mu_HS23/36")}), Rule("HS37", {(G, "Mu_HS37")}),
    Rule("HS38", {(D, "Mu_HS38")}), Rule("HS40", {(D, "Mu_HS40")}), Rule("HS41", {(A, "Mu_HS41")}), Rule("HS42", {(B, "Mu_HS42")}),
    Rule("HS43", {(A, "Mu_HS4A")}), Rule("HS44", {(G, "Mu_HS44")}), Rule("HS45", {(G, "Mu_HS45")}, must_not={(B, "Mu_HS31"), (D, "Mu_HS32"), (D, "Mu_HS60")}),
    Rule("HS50", {(A, "Mu_HS4A")}), Rule("HS52", {(D, "Mu_HS52")}), Rule("HS53", {(A, "Mu_HS53")}), Rule("HS55", {(D, "Mu_HS55")}),
    Rule("HS57", {(B, "Mu_HS57")}), Rule("HS58", {(D, "Mu_HS58")}, allowed_extra={(A, "Mu_HS15")}), Rule("HS60", {(D, "Mu_HS60"), (G, "Mu_HS45")}),
    Rule("HS62", {(A, "Mu_HS4A")}), Rule("HS63", {(A, "Mu_HS63")}), Rule("HS64", {(A, "Mu_HS4A"), (B, "Mu_HS4B")}), Rule("HS65", {(A, "Mu_HS4A")}),
]
LPXA = (G, "lpxA")


def primers_to_pairs(primer_names: Iterable[str]) -> Set[Prim]:
    pairs: Set[Prim] = set()
    for raw in primer_names:
        p = norm_primer(str(raw).strip())
        mix = PRIMER_TO_MIX.get(p)
        if mix:
            pairs.add((mix, p))
    return pairs


def call_capsule_types(detected: Iterable[Prim], require_lpxA: bool = False) -> List[str]:
    det: Set[Prim] = {(m, norm_primer(p)) for m, p in detected}
    if require_lpxA and LPXA not in det:
        return ["Uninterpretable (lpxA negative)"]
    calls: List[str] = []
    for r in RULES:
        if not r.must_have.issubset(det):
            continue
        if any((m, p) in det for (m, p) in r.must_not):
            continue
        calls.append(r.name)
    if "HS4 complex" in calls:
        calls = ['HS4']
    return sorted(set(calls)) or ["Untypeable"]


def call_capsule_types_from_primernames(primer_names: Iterable[str], require_lpxA: bool = False) -> List[str]:
    return call_capsule_types(primers_to_pairs(primer_names), require_lpxA=require_lpxA)


def serotype_Cb(Pre):  # 空肠弯曲菌
    sers = '-'
    subprocess.run(f' cat {Pre}.final.fasta |seqkit amplicon -p /data/test/CB_type/CB_primer.tsv --bed > {Pre}.CB_primer.tsv', shell=True)
    serodict = {'样本名称': Pre, '血清型': '-', '物种可靠': '否'}
    if os.path.isfile(f'{Pre}.CB_primer.tsv') and os.path.getsize(f'{Pre}.CB_primer.tsv') != 0:
        CBdb = pd.read_table(f'{Pre}.CB_primer.tsv', header=None)
        if CBdb.shape[0] > 0:
            if 'lpxA' in CBdb[3].tolist():
                serodict['物种可靠'] = '是'
            sers = call_capsule_types_from_primernames(CBdb[3].tolist())[0]
            sers = sers.replace('HS', 'HS:')
            serodict['血清型'] = sers
    serodb = pd.DataFrame(serodict, index=[0])
    serodb.to_csv(f'{Pre}_serotype_result.tsv', sep='\t', index=False)
    return sers


def extract_SpeID(row):
    if pd.isna(row['FullLineageRanks']) or pd.isna(row['FullLineage']):
        return 'noSpe'
    ranks = row['FullLineageRanks'].split(';')
    if 'species' not in ranks:
        return 'noSpe'
    idx = ranks.index('species')
    lineage = row['FullLineageTaxIDs'].split(';')
    return lineage[idx] if idx < len(lineage) else 'noSpe'


def extract_Spe(row):
    if pd.isna(row['FullLineageRanks']) or pd.isna(row['FullLineage']):
        return 'noSpe'
    ranks = row['FullLineageRanks'].split(';')
    if 'species' not in ranks:
        return 'noSpe'
    idx = ranks.index('species')
    lineage = row['FullLineage'].split(';')
    return lineage[idx] if idx < len(lineage) else 'noSpe'


def is_non_numeric_in_bracket(x):
    if pd.isna(x):
        return False
    m = re.search(r'\((.*?)\)', str(x))
    return bool(m) and (not m.group(1).isdigit())


def mlst_serotype(Pre, tSpe):
    runtime = get_runtime_context()
    requested_scheme = tSpe or runtime.species
    cmdb = pd.read_table(f'{Pre}.checkm.tsv')
    Asdb = pd.read_table('Assem_info1.tsv')
    cfile = pytaxonkit.lineage(Asdb['taxid'].tolist())
    cfile['Species'] = cfile.apply(extract_SpeID, axis=1)
    cfile['SpeciesName'] = cfile.apply(extract_Spe, axis=1)
    tcfile = cfile.merge(Asdb, left_on='TaxID', right_on='taxid').drop_duplicates().groupby('Species').sum('序列长度').reset_index().sort_values('序列长度', ascending=False)
    krspeid = tcfile.loc[tcfile['Species'] != 'noSpe']['Species'].tolist()[0]
    krspe = cfile.loc[cfile['Species'] == krspeid, 'SpeciesName'].tolist()[0]
    krspeidlist = cfile.loc[cfile['Species'] == str(krspeid), 'TaxID'].tolist()
    mainper = int(Asdb.loc[Asdb['taxid'].isin(krspeidlist), '序列长度'].sum() / Asdb.序列长度.sum() * 100)
    if 'noSpe' in tcfile['Species'].tolist():
        noSpeper = round(tcfile.loc[tcfile['Species'] == 'noSpe', '序列长度'].tolist()[0] / tcfile.序列长度.sum() * 100, 2)
        cmdb['物种名称'] = f'{krspe}({mainper}%) noSpe({noSpeper}%)'
    else:
        cmdb['物种名称'] = f'{krspe}({mainper}%)'
    cmdb.to_csv(f'{Pre}.checkm.tsv', sep='\t', index=False)

    subprocess.run(f"mlst --quiet --csv {Pre}.final.fasta > {Pre}_mlst.csv", shell=True)
    mlst_gene = pd.read_table(f"{Pre}_mlst.csv", sep=",", header=None)
    if mlst_gene.shape[1] <= 4:
        if requested_scheme:
            subprocess.run(f"mlst --scheme {requested_scheme} --quiet --csv {Pre}.final.fasta > {Pre}_mlst.csv", shell=True)
        else:
            subprocess.run(f"mlst --scheme ecoli --quiet --csv {Pre}.final.fasta > {Pre}_mlst.csv", shell=True)
    mlst_gene = pd.read_table(f"{Pre}_mlst.csv", sep=",", header=None)
    sub = mlst_gene.iloc[:, 3:]
    if sub.applymap(is_non_numeric_in_bracket).sum(axis=1).tolist()[0] < 3:
        mlst_B = mlst_gene.iloc[0, 1]
    else:
        mlst_B = 'UnKnown'
    mlst_st = mlst_gene.iloc[0, 2]
    for x in mlst_gene.iloc[0, :].tolist():
        x = str(x)
        if '(' in x:
            gene = x.split('(')[0]
            gene_num = x.split('(')[1].split(')')[0]
            try:
                os.system(f"seqkit grep -p {gene}_{gene_num} /home/dell/miniconda3/envs/TB_ONT/db/pubmlst/{mlst_B}/{gene}.tfa >> mlst_all.fa")
            except Exception:
                pass
    with open('mlst.log', 'w') as mlstg:
        subprocess.run(f"dnadiff mlst_all.fa {Pre}.final.fasta", shell=True, stdout=mlstg, stderr=mlstg)
        subprocess.run("show-coords -lTH out.delta|sort -k7nr > mlst.coords", shell=True, stdout=mlstg, stderr=mlstg)
    if os.path.isfile('mlst.coords') and os.path.getsize('mlst.coords') != 0:
        mlst_gene = pd.read_table('mlst.coords', header=None)
        mlst_gene.columns = ["起始位置", "终止位置", "比对起始位置", "比对终止位置", "基因长度", "比对长度", "一致性%", "基因长度", "序列长度", "管家基因", "序列名称"]
        sss = mlst_gene.pop("管家基因")
        mlst_gene.insert(0, "管家基因", sss)
        sss = mlst_gene.pop("序列名称")
        mlst_gene.insert(3, "序列名称", sss)
        mlst_gene["序列分型(ST)"] = mlst_st
        mlst_gene["物种信息"] = mlst_B
        mlst_gene.pop("序列长度")
        mlst_gene['管家基因序号'] = mlst_gene['管家基因'].str.split('_').str[1]
        mlst_gene = mlst_gene[['管家基因', '管家基因序号', '起始位置', '终止位置', '序列名称', '比对起始位置', '比对终止位置', '基因长度', '比对长度', '一致性%', '序列分型(ST)', '物种信息']]
        mlst_gene.to_csv(f"{Pre}.mlst_Stat.txt", index=0, sep='\t')
        for x in mlst_gene["管家基因"].tolist():
            y = mlst_gene[mlst_gene["管家基因"] == x]["序列名称"].tolist()[0]
            os.system(f"show-aligns out.delta {x} {y}|sed -n '/-- Alignments/,/--   END/p' > {x}_gene_show.txt")

    cdb = pd.read_table(f'{Pre}.checkm.tsv')
    cdb['mlst 物种名称'] = mlst_B
    print(mlst_B, Pre)
    cdb[['样本名称', '物种名称', 'mlst 物种名称', '污染率', '完整性']].to_csv(f'{Pre}.checkm.tsv', sep='\t', index=False)

    if mlst_B == 'bordetella_3':
        bp_vaccine(Pre)
        if os.path.isfile(f'{Pre}.R1.fastq.gz'):
            bp_2037(Pre)
        else:
            fa_2037(Pre)
        bp_mlva(Pre)
    elif mlst_B == 'klebsiella':
        serotype_kb(Pre)
    elif 'ecoli' in mlst_B:
        serotype_B(Pre)
    elif 'senterica' in mlst_B:
        serotype_A(Pre)
    elif 'hinfluenzae' in mlst_B:
        serotype_HI(Pre)
    elif 'vparahaemolyticus' in mlst_B:
        serotype_D(Pre)
    elif 'spyogenes' in mlst_B:
        serotype_groupA(Pre)
    elif 'vcholerae' in mlst_B:
        serotype_E(Pre)
    elif 'saureus' in mlst_B:
        serotype_st(Pre)
    elif 'listeria' in mlst_B:
        serotype_lm(Pre)
    elif 'bcereus' in mlst_B:
        serotype_bt(Pre)
    elif 'mpneumoniae' in mlst_B:
        mp_mlva(Pre)
    elif 'ssuis' in mlst_B:
        serotype_SS(Pre)
    elif 'campylobacter' in mlst_B:
        serotype_Cb(Pre)

    if os.path.isfile(f'{Pre}_2.report.txt'):
        kradb = pd.read_table(f'{Pre}_2.report.txt', header=None)
    else:
        kradb = pd.read_table(f'{Pre}_assem.kraken2.txt', header=None)
    if kradb.loc[kradb[3] == 'G', 5].tolist()[0].strip() == 'Yersinia':
        serotype_ys(Pre)
    if requested_scheme:
        PathoNet(Pre, requested_scheme)
    else:
        PathoNet(Pre, mlst_B)
