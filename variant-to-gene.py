from sklearn.preprocessing import quantile_transform
import pandas as pd
import numpy as np
import os
import gzip
import chardet
import re

def transform_data(filename):
    rawdata = pd.read_table(filename,sep='\t')
    data = np.array(rawdata['-log10(Pvalue)'], dtype=np.float64)
    transformed_data = quantile_transform(data.reshape(-1, 1), n_quantiles=10, random_state=0)
    rawdata['transform']=transformed_data
    return rawdata

def transform(path):
    for filename in os.listdir(path):
      if filename.endswith('_ch.txt'):
          name=os.path.splitext(filename)[0]
          d=transform_data(os.path.join(path, filename))
          d.to_csv(os.path.join(path, name+'_transform.txt'),index=False,sep='\t')

#transform('./eQTL/DICE')
#transform('./eQTL/ImmuNexUT/E-GEAD-398')
#transform('./eQTL/BLUEPRINT/QTL/EGAD00001005200-significant_association')
#transform('./eQTL/Gtex/GTEx_Analysis_v8_eQTL')
#transform('./eQTL/eQTLGen')
#transform('./eQTL/scRNA-seq/Science2022-Single-cell_autoimmune_disease')

def score(path,i,r,a,g,s):
    output = open(os.path.join(path, 'all_score.txt'), "w")
    data_dict = {}
    for filename in os.listdir(path):
        if filename.endswith('_transform.txt'):
            input = open(os.path.join(path, filename), "r")
            #CHROM  POS     ID      REF     ALT     GeneSymbol      -log10(Pvalue)  transform
            #chr1    752566  rs3094315       G       A       LINC00115       4.89229 0.7732054068910376
            for line in input:
                list = line.rstrip('\n').split('\t')
                if line.startswith('#'):
                    continue
                #id, ref, alt,gene,pval,score = list[2:8]
                id, ref, alt,gene,score = list[i],list[r],list[a],list[g],list[s]
                key = id + '\t' + ref + '\t' + alt
                if key in data_dict and gene in data_dict[key]:
                    if data_dict[key][gene] < score:
                        data_dict[key][gene] = score
                else:
                    data_dict.setdefault(key,{})[gene]=score  #二维
                  #  data_dict.setdefault(key,score)   #一维
                  #  data_dict.get(key,score)   #一维
            input.close()
    
    output.write('#ID\tREF\tALT\tGENE_SCORE\n')
    #ds = sorted(data_dict.items(), key=lambda x: x[1], reverse=True)  #一维排序
    #ds= {k: {i: dict(sorted(j.items(), key=lambda x: x[1], reverse=True))} for k, v in data_dict.items() for i, j in v.items()} #三维排序
    ds= {k: dict(sorted(v.items(), key=lambda x: x[1], reverse=True)) for k, v in data_dict.items()} #二维排序
    for k, v in ds.items():
        if type(v) is dict:
            info=""
            for nk, nv in v.items():
                info += nk + ':' + nv + '|'
            output.write(k + '\t' + info + '\n')
    
    output.close()
#score('./eQTL/DICE',2,3,4,5,7)
#score('./eQTL/ImmuNexUT/E-GEAD-398',0,1,2,4,6)
#score('./eQTL/BLUEPRINT/QTL/EGAD00001005200-significant_association',1,0,0,2,4)
#score('./eQTL/Gtex/GTEx_Analysis_v8_eQTL',0,0,0,1,3)
#score('./eQTL/eQTLGen',0,3,4,5,7)
#score('./eQTL/scRNA-seq/Science2022-Single-cell_autoimmune_disease',2,3,4,1,7)

#rawdata = pd.read_table('./eQTL/Capture_Hi-C/hg19/OSF/combine_allcell.txt',sep='\t',header=None)
def transform_hic(filename):
    rawdata = pd.read_table(filename,sep='\t',header=None)
    data = np.array(rawdata[4], dtype=np.float64)
    transformed_data = quantile_transform(data.reshape(-1, 1), n_quantiles=10, random_state=0)
    rawdata['transform']=transformed_data
    rawdata.to_csv('./eQTL/Capture_Hi-C/hg19/OSF/combine_allcell_transform.txt',index=False,sep='\t')
#transform_hic('./eQTL/Capture_Hi-C/hg19/OSF/combine_allcell.txt')

def pos_file(path):
    input = open(path, "r",encoding="ISO-8859-1")
    POS_grch37={}
    POS_grch38={}
    ID={}
    for line in input:
        list = line.rstrip('\n').split('\t')
        if line.startswith('#LOCI'):
            continue
        id1, id2, grch37, grch38 = list[9],list[12],list[11],list[10]
        POS_grch37.setdefault(grch37,1)
        POS_grch38.setdefault(grch38,1)
        ID.setdefault(id1,1)
        ID.setdefault(id2,1)
    return POS_grch37,POS_grch38,ID
POS_grch37,POS_grch38,ID_RS = pos_file('./all_0.8_expand_current')

def qtl_ano(path):
    input = open(path, "r")
    dict={}
    for line in input:
        list = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            continue
        #id, ref, alt,gene,pval,score = list[2:8]
        id, ref, alt, gene = list[0:4]
        if id in ID_RS:
            dict.setdefault(id,gene)
    return dict

def gtex_ano(path):
    input = open(path, "r")
    dict={}
    #chr1_64764_C_T_b38      chr1_64764_C_T_b38      chr1_64764_C_T_b38      WASH7P:0.6285187949198925|
    for line in input:
        list = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            continue
        #id, ref, alt,gene,pval,score = list[2:8]
        id, ref, alt, gene = list[0:4]
        id = id.replace("chr", "")
        id1 = id.split('_')[0]+':'+id.split('_')[1]
        if id1 in POS_grch38:
            dict.setdefault(id1,gene)
    return dict

DICE = qtl_ano('./eQTL/DICE/all_score.txt')
IMU = qtl_ano('./eQTL/ImmuNexUT/E-GEAD-398/all_score.txt')
BLUE = qtl_ano('./eQTL/BLUEPRINT/QTL/EGAD00001005200-significant_association/all_score.txt')
GTEX = gtex_ano('./eQTL/Gtex/GTEx_Analysis_v8_eQTL/all_score.txt')
QTLGEN = qtl_ano('./eQTL/eQTLGen/all_score.txt')
SCRNA = qtl_ano('./eQTL/scRNA-seq/Science2022-Single-cell_autoimmune_disease/all_score.txt')

def ABC_CELL_TYPE(path):
    input = open(path, "r")
    #B_cell-ENCODE
    dict={}
    for line in input:
        list = line.rstrip('\n').split('\t')
        id = list[0]
        dict.setdefault(id,1)
    return dict
ABC_CELL = ABC_CELL_TYPE('./eQTL/Activity-By-Contact/immune_cell')

def anno_ABC(path):
    input = gzip.open(path, 'rb')
    abc={}
    ABC={}
    input.readline()  # 跳过第一行
    #less AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz|awk -F'\t' '{print $1,$2,$3,$7,$21}'
    #chr1 710010 710210 LOC100288069 0.025955
    #2:24856901
    for line in input:
        list = line.decode().rstrip('\n').split('\t')
        #id, ref, alt,gene,pval,score = list[2:8]
        chr, start, end, gene,score,cell_type = list[0],int(list[1]),int(list[2]),list[6],list[20],list[23]
        chr = chr.replace("chr", "")
        if cell_type in ABC_CELL:
          for num in range(start,end+1):
            id = chr+':'+str(num)
            if id in POS_grch37:
                if id in abc and gene in abc[id] and abc[id][gene]<score:
                    abc[id][gene] = score
                else:
                    abc.setdefault(id,{})[gene]=score
    ds= {k: dict(sorted(v.items(), key=lambda x: x[1], reverse=True)) for k, v in abc.items()} #二维排序
    for k, v in ds.items():
        if type(v) is dict:
            info=""
            for nk, nv in v.items():
                info += nk + ':' + nv + '|'
            ABC.setdefault(k,info)
    return ABC
ABC_MODEL=anno_ABC('./eQTL/Activity-By-Contact/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz')

def GENE_ano(path):
    input = open(path, "r")
    #ENSG00000000003 TSPAN6
    dict={}
    for line in input:
        list = line.rstrip('\n').split('\t')
        id, gene = list[0:2]
        dict.setdefault(id,gene)
    return dict
GENE = GENE_ano('./eQTL/EpiMap/hg19/gene_ch')

def anno_EpiMap(path):
    epi={}
    map={}
    for filename in os.listdir(path):
      if filename.endswith('.tsv.gz'):
         input = gzip.open(os.path.join(path, filename), 'rb')
         #chr1    1366120 1366360 ENSG00000179403 0.9632734656333923      E10
         for line in input:
            list = line.decode().rstrip('\n').split('\t')
            chr, start, end, gene,score = list[0],int(list[1]),int(list[2]),list[3],list[4]
            chr = chr.replace("chr", "")
            if gene in GENE:
                g=GENE[gene]
            else:
                next
            for num in range(start,end+1):
                id = chr+':'+str(num)
                if id in POS_grch37:
                    if id in epi and g in epi[id] and epi[id][g]<score:
                        epi[id][g] = score
                    else:
                        epi.setdefault(id,{})[g]=score
    ds= {k: dict(sorted(v.items(), key=lambda x: x[1], reverse=True)) for k, v in epi.items()} #二维排序
    for k, v in ds.items():
        if type(v) is dict:
            info=""
            for nk, nv in v.items():
                info += nk + ':' + nv + '|'
            map.setdefault(k,info)
    return map
EpiMap=anno_EpiMap('./eQTL/EpiMap/hg19')

def anno_HiC(path):
    input = open(path, "r")
    hic={}
    HIC={}
    input.readline()  # 跳过第一行
    #8       124700922       124703675       MIR548D1;RNU6-875P      0.5129132512114258
    for line in input:
        list = line.rstrip('\n').split('\t')
        chr, start, end,gene,score = list[0],int(list[1]),int(list[2]),list[3],list[4]
        for num in range(start,end+1):
            id = chr+':'+str(num)
            if id in POS_grch37:
                if id in hic and gene in hic[id] and hic[id][gene]<score:
                    hic[id][gene] = score
                else:
                    hic.setdefault(id,{})[gene]=score
    ds= {k: dict(sorted(v.items(), key=lambda x: x[1], reverse=True)) for k, v in hic.items()} #二维排序
    for k, v in ds.items():
        if type(v) is dict:
            info=""
            for nk, nv in v.items():
                info += nk + ':' + nv + '|'
            HIC.setdefault(k,info)
    return HIC
Captur_HIC=anno_HiC('./eQTL/Capture_Hi-C/hg19/OSF/combine_allcell_transform_max.txt')

def VLS2G_ano(path):
    input = open(path, "r")
    #rs10000518      4_11502867_G_A  HS3ST1:0.05975855130784708|
    dict={}
    for line in input:
        list = line.rstrip('\n').split('\t')
        id, gene = list[0],list[-1]
        dict.setdefault(id,gene)
    return dict
V2G = VLS2G_ano('./eQTL/Open_Targets/rs_expand_pos_v2g_new_use.txt')
L2G= VLS2G_ano('./eQTL/Open_Targets/rs_lead_ch_l2g_ch_new_use.txt')
cS2G= VLS2G_ano('./eQTL/cS2G/rs_expand_new_SGscore_ch.txt')

def consequence_ano(path):
    input = open(path, "r")
    input.readline()  # 跳过第一行
    #awk -F'\t' '{print $1,$7,$20}' rs_hg38_all_vep_filter_rmHLA_rmPseudogene.txt
    dict={}
    for line in input:
        list = line.rstrip('\n').split('\t')
        id, impact,gene,biotype = list[0],list[6],list[19],list[22]
        g=impact+'\t'+gene+"|"+biotype
        dict.setdefault(id,g)
    return dict
consequence = consequence_ano('./consequence_annotation/rs_hg38_all_vep_filter_rmHLA_rmPseudogene.txt')

def decide(a,b,c):
    d=''
    if a in c:
        d=c[a]
    elif b in c:
        d=c[b]
    else:
        d= 0
    return d

def decide_pos(a,b):
    d=''
    if a in b:
        d=b[a]
    else:
        d= 0
    return d

def decide_con(a,b,c):
    d=''
    if a in c:
        d=c[a]
    elif b in c:
        d=c[b]
    else:
        d= '0\t0'
    return d

weight = {'dice_w': 0.66, 
          'immu_w': 0.66, 
          'blue_w': 0.66,
          'gtex_w': 0.66,
          'eqtl_w': 0.66,
          'scrna_w': 0.66,
          
          'abc_w': 0.33,
          'epi_w': 0.33,
          'hic_w': 0.33,
          
          'stop_gained': 1,
          'stop_lost': 1,
          'frameshift_variant': 1,
          'missense_variant': 1,
          'synonymous_variant': 1,
          'splice_acceptor_variant': 1,
          'splice_donor_variant': 1,
          'splice_region_variant': 1,
          '3_prime_UTR_variant': 1,
          '5_prime_UTR_variant': 1,
          'intron_variant': 0.66,
          
          'v2g_w': 0.66,
          'l2g_w': 0.66,
          'cs2g_w': 0.66}


def annotation_file(path):
    output = open('./consequence_annotation/all_0.8_expand_current_summary.txt', "w",encoding="utf-8")
    input = open(path, "r",encoding="ISO-8859-1")
    for line in input:
        line=line.rstrip('\n')
        list = line.split('\t')
        if line.startswith('#LOCI'):
            head='DICE\tImmuNexUT\tBLUEPRINT\tGTEx\teQTLGen\tscRNA-seq_of_autoimmune_disease\tABC_model\tEpiMap\tPCHiC\tConsequence\tSYMBOL\tV2G\tL2G\tcS2G\teqtl_total_score\tall_total_score'
            output.write(line+'\t'+str(head)+'\n')
            continue
        id1, id2, grch37, grch38 = list[9],list[12],list[11],list[10]
        dice=decide(id1,id2,DICE)
        immu=decide(id1,id2,IMU)
        blue=decide(id1,id2,BLUE)
        gtex=decide_pos(grch38,GTEX)
        eqtlgen=decide(id1,id2,QTLGEN)
        scrna=decide(id1,id2,SCRNA)
        # medium value for each gene in above QTL studies
        qtl_medium={}
        for k in dice,immu,blue,gtex,eqtlgen,scrna:
            #ADCY3:0.8969166907137364|DNAJC27:0.2137321178599474|DTNB:0.10991217373095351|	0
            if k==0:
                continue
            else:
                for i in k.split('|'):
                    j=i.split(':')
                    if len(j) < 2:
                        continue
                    qtl_medium.setdefault(j[0],[]).append(float(j[1]))
        tmp={}
        qtl_total = 0
        for key, value in qtl_medium.items():
            if type(value) is __builtins__.list:     #"DON'T" use list keywords which are reserved by python!
            #if isinstance(value, __builtins__.list):
            #   median_value = np.median(value)    #中位值
               median_value = max(value)    #最大值
               tmp[key] = median_value
        if len(tmp) > 0:
            tmp_ds = sorted(tmp.items(), key=lambda x: x[1], reverse=True)
            qtl_total = '|'.join([i[0]+':'+str(i[1]) for i in tmp_ds])      #ADCY3:0.8969166907137364|POMC:0.7035668608782154|DNAJC27-AS1:0.4986696602536807|TP53I3:0.3779836541869014

        abc_model=decide_pos(grch37,ABC_MODEL)   #DNAJC27:0.028813|DNAJC27-AS1:0.028813|CENPO:0.025434|PTRHD1:0.025433|EFR3B:0.0221|ITSN2:0.021639|ADCY3:0.019751|DNMT3A:0.017732|PFN4:0.016092|UBXN2A:0.015056|
        epimap=decide_pos(grch37,EpiMap)   #ADCY3:0.7161606550216675|
        cap_hic=decide_pos(grch37,Captur_HIC)  #AAMP;PNKD:0.5207998732344372|CXCR1:0.3381709127515197|ARPC2:0.2336194872883838|
        
        conse=decide_con(id1,id2,consequence)  #3_prime_UTR_variant	FOSL2|protein_coding
        conse_tmp=conse.split('|')[0]
        
        v2g=decide(id1,id2,V2G)   #FOSL2:0.30744466800804826|MRPL33:0.05975855130784708|PLB1:0.05975855130784708|RBKS:0.04647887323943661
        l2g=decide(id1,id2,L2G)   #FOSL2:0.7952556014060974|PLB1:0.036290884017944336|BABAM2:0.016975289210677147
        cs2g=decide(id1,id2,cS2G)  #FOSL2:1|

        total = {}
        def process_data(data_list, weight_factor):
            for k in data_list:
                if k == 0:
                    continue
                else:
                    for i in k.split('|'):
                        j = i.split(':')
                        if len(j) < 2:
                            continue
                        t = j[0].split(';')
                        u = j[1] + '|' + str(weight_factor)
                        for q in t:
                            total.setdefault(q, []).append(u)

        data_list_1 = [qtl_total, v2g, l2g, cs2g]
        weight_factor_1 = 0.66
        process_data(data_list_1, weight_factor_1)

        data_list_2 = [abc_model, epimap, cap_hic]
        weight_factor_2 = 0.33
        process_data(data_list_2, weight_factor_2)
        if conse != 0:
            i=conse.split('\t')
            j = i[0].split(',')[0]
            if j in weight:
                num = weight[j]
                u = '1|' + str(num)
                ty,ts=i[1].split('|')[0:2]
                if ts != 'lncRNA' and ts != 'retained_intron':    #lncRNA and retained_intron are not included in the final output
                    total.setdefault(ty, []).append(u)

        #{'ADCY3': ['0.8969166907137364|0.66', '0.019751|0.33', '0.7161606550216675|0.33'], 'POMC': ['0.7035668608782154|0.66'], 'DNAJC27-AS1': ['0.4986696602536807|0.66', '0.028813|0.33'], 'TP53I3': ['0.3779836541869014|0.66'], 'FOSL2': ['0.30744466800804826|0.66', '0.7952556014060974|0.66', '1|0.66'], 'MRPL33': ['0.05975855130784708|0.66'], 'PLB1': ['0.05975855130784708|0.66', '0.036290884017944336|0.66'], 'RBKS': ['0.04647887323943661|0.66'], 'BABAM2': ['0.016975289210677147|0.66'], 'DNAJC27': ['0.028813|0.33'], 'CENPO': ['0.025434|0.33'], 'PTRHD1': ['0.025433|0.33'], 'EFR3B': ['0.0221|0.33'], 'ITSN2': ['0.021639|0.33'], 'DNMT3A': ['0.017732|0.33'], 'PFN4': ['0.016092|0.33'], 'UBXN2A': ['0.015056|0.33'], 'AAMP': ['0.5207998732344372|0.33'], 'PNKD': ['0.5207998732344372|0.33'], 'CXCR1': ['0.3381709127515197|0.33'], 'ARPC2': ['0.2336194872883838|0.33']}
        tmp={}
        final_score=0
        gene_num={}
        for key, value in total.items():
            tmp_num=0
            tmp_num2=0
            if len(value) > 1:
               for j in value:
                   m,n=j.split('|')[0:2]
                   tmp_num += float(m)*float(n)
                   tmp_num2 += float(n)
                   gene_num[key] += 1
               tmp[key] = tmp_num/tmp_num2
            else:
               tmp[key] = float(value[0].split('|')[0]) * float(value[0].split('|')[1])
               gene_num[key] += 1
        if len(tmp) > 0:
            tmp_ds = sorted(tmp.items(), key=lambda x: x[1], reverse=True)
            final_score = '|'.join([i[0]+':'+str(i[1])+':'+str(gene_num[i[0]]) for i in tmp_ds])

        output.write(line+'\t'+str(dice)+'\t'+str(immu)+'\t'+str(blue)+'\t'+str(gtex)+'\t'+str(eqtlgen)+'\t'+str(scrna)+'\t'+str(abc_model)+'\t'+str(epimap)+'\t'+str(cap_hic)+'\t'+str(conse_tmp)+'\t'+str(v2g)+'\t'+str(l2g)+'\t'+str(cs2g)+'\t'+str(qtl_total)+'\t'+str(final_score)+'\n')

annotation_file('./consequence_annotation/all_0.8_expand_current')        
