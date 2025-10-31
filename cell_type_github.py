import sys, getopt
import os
import gzip
import re

try:
    opts, args = getopt.getopt(sys.argv[1:], "", ["inputfile=", "cell_type=", "gene=", "level=", "outputfile="])
except getopt.GetoptError:
    print("Example: python", sys.argv[0], "--inputfile inputfile --cell_type cell_type --gene gene_yes/gene_no --level L1/L2/L3 --outputfile outputfile")
    sys.exit(2)

inputfile, cell_type, gene, level, outputfile = None, None, None, None, None

for opt, arg in opts:
    if opt == "--inputfile":
        inputfile = arg
    elif opt == "--cell_type":
        cell_type = arg
    elif opt == "--gene":
        gene = arg
    elif opt == "--level":
        level = arg
    elif opt == "--outputfile":
        outputfile = arg

if None in (inputfile, cell_type, gene, level, outputfile):
    print("Missing required arguments. Usage:")
    print("Example: python", sys.argv[0], "--inputfile inputfile --cell_type cell_type --gene gene_yes/gene_no --level L1/L2/L3 --outputfile outputfile")
    sys.exit(2)

def CELL_TYPE(path):
    global level
    input_file = open(path, "r")
    #CD4_NAIVE	T_cell	CD4T_cell	CD4T_naive	DICE    eQTL
    input_file.readline()  # 跳过第一行
    new_dict={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        id = li[0]
        va = li[1] if level == "L1" else li[2] if level == "L2" else li[3] if level == "L3" else "NA"
        if va != "NA" and len(va) > 0:
            new_dict.setdefault(id,va)
    input_file.close()
    return new_dict
CELL = CELL_TYPE(cell_type)

def pos_file(path):
    input_file = open(path, "r",encoding="ISO-8859-1")
    POS_grch37={}
    POS_grch38={}
    ID={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        if line.startswith('#LOCI'):
            continue
        id1, id2, grch37, grch38 = li[3],li[6],li[5],li[4]
        POS_grch37.setdefault(grch37,1)
        POS_grch38.setdefault(grch38,1)
        ID.setdefault(id1,1)
        ID.setdefault(id2,1)
    input_file.close()
    return POS_grch37,POS_grch38,ID
POS_grch37,POS_grch38,ID_RS = pos_file(inputfile)


def qtl_ano(path):
    input_file = open(path, "r")
    new_dict={}
    rs={}
    rs_new={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        for num in range(1,len(li)):
            if line.startswith('Variants'):
                new_dict.setdefault(num,li[num])
            else:
                if li[0] in ID_RS and li[num] != "0" and new_dict[num] in CELL:
                    rs.setdefault(li[0], []).append(CELL[new_dict[num]])
    for key, value in rs.items():
        rs_new[key] = list(set(value))
    input_file.close()
    return rs_new

def gtex_ano(path):
    input_file = open(path, "r")
    new_dict={}
    rs={}
    rs_new={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        for num in range(1,len(li)):
            if line.startswith('Variants'):
                new_dict.setdefault(num,li[num])
            else:
               id = li[0].replace("chr", "")
               id1 = id.split('_')[0]+':'+id.split('_')[1]
               if id1 in POS_grch38 and li[num] != "0" and new_dict[num] in CELL:
                    rs.setdefault(id1, []).append(CELL[new_dict[num]])
    for key, value in rs.items():
        rs_new[key] = list(set(value))
    input_file.close()
    return rs_new

def qtl_ano_gene(path):
    input_file = open(path, "r")
    new_dict={}
    rs={}
    rs_new={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        for num in range(1,len(li)):
            if line.startswith('Variants'):
                new_dict.setdefault(num,li[num])
            else:
                if li[0] in ID_RS and li[num] != "0" and new_dict[num] in CELL:
                  for i in li[num].split('|'):
                    j=i.split('[')[0]
                    rs.setdefault(li[0], {}).setdefault(j, []).append(CELL[new_dict[num]])
    for k, v in rs.items():
        if type(v) is dict:
            for nk, nv in v.items():
                rs_new.setdefault(k,{})[nk]=list(set(rs[k][nk]))
    input_file.close()
    return rs_new

def gtex_ano_gene(path):
    input_file = open(path, "r")
    new_dict={}
    rs={}
    rs_new={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        for num in range(1,len(li)):
            if line.startswith('Variants'):
                new_dict.setdefault(num,li[num])
            else:
               id = li[0].replace("chr", "")
               id1 = id.split('_')[0]+':'+id.split('_')[1]
               if id1 in POS_grch38 and li[num] != "0" and new_dict[num] in CELL:
                    for i in li[num].split('|'):
                      j=i.split('[')[0]
                      rs.setdefault(id1, {}).setdefault(j, []).append(CELL[new_dict[num]])
    for k, v in rs.items():
        if type(v) is dict:
            for nk, nv in v.items():
                rs_new.setdefault(k,{})[nk]=list(set(rs[k][nk]))
    input_file.close()
    return rs_new


def anno_ABC_gen(path):
    input_file = open(path, 'r')
    abc={}
    ABC={}
    input_file.readline()  # 跳过第一行
    #less AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz|awk -F'\t' '{print $1,$2,$3,$7,$21}'
    #chr1 710010 710210 LOC100288069 0.025955
    #2:24856901
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        #id, ref, alt,gene,pval,score = li[2:8]
        chr, start, end, gene,score,cell_type = li[0],int(li[1]),int(li[2]),li[6],li[20],li[23]
        chr = chr.replace("chr", "")
        if cell_type in CELL:
          for num in range(start,end+1):
            id = chr+':'+str(num)
            if id in POS_grch37 and gene != "" and cell_type in CELL:
#                if id in abc and gene in abc[id]:
                    abc.setdefault(id, {}).setdefault(gene, []).append(CELL[cell_type])
    for k, v in abc.items():
        if type(v) is dict:
            for nk, nv in v.items():
                ABC.setdefault(k,{})[nk]=list(set(abc[k][nk]))
    return ABC
def anno_ABC(path):
    input_file = open(path, 'r')
    abc={}
    ABC={}
    input_file.readline()  # 跳过第一行
    #less AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz|awk -F'\t' '{print $1,$2,$3,$7,$21}'
    #chr1 710010 710210 LOC100288069 0.025955
    #2:24856901
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        #id, ref, alt,gene,pval,score = li[2:8]
        chr, start, end, gene,score,cell_type = li[0],int(li[1]),int(li[2]),li[6],li[20],li[23]
        chr = chr.replace("chr", "")
        if cell_type in CELL:
          for num in range(start,end+1):
            id = chr+':'+str(num)
            if id in POS_grch37 and cell_type in CELL:
                abc.setdefault(id, []).append(CELL[cell_type])
    for key, value in abc.items():
        ABC[key] = list(set(value))
    return ABC

def GENE_ano(path):
    input_file = open(path, "r")
    #ENSG00000000003 TSPAN6
    new_dict={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        id, gene = li[0:2]
        new_dict.setdefault(id,gene)
    return new_dict
GENE = GENE_ano('./eQTL/EpiMap/hg19/gene_ch')

def anno_EpiMap_gen(path):
    epi={}
    map={}
    for filename in os.listdir(path):
      if filename.endswith('.tsv.gz'):
         n=filename.replace("_collated_pred.tsv.gz", "")
         input_file = gzip.open(os.path.join(path, filename), 'rb')
         #chr1    1366120 1366360 ENSG00000179403 0.9632734656333923      E10
         for line in input_file:
            li = line.decode().rstrip('\n').split('\t')
            chr, start, end, gene,score = li[0],int(li[1]),int(li[2]),li[3],li[4]
            chr = chr.replace("chr", "")
            if gene in GENE:
                g=GENE[gene]
            else:
                next
            for num in range(start,end+1):
                id = chr+':'+str(num)
                if id in POS_grch37 and n in CELL and g != "" and n in CELL:
#                    if id in epi and g in epi[id]:
                        epi.setdefault(id, {}).setdefault(g, []).append(CELL[n])

    for k, v in epi.items():
        if type(v) is dict:
            for nk, nv in v.items():
                map.setdefault(k,{})[nk]=list(set(epi[k][nk]))
    return map
def anno_EpiMap(path):
    epi={}
    map={}
    for filename in os.listdir(path):
      if filename.endswith('.tsv.gz'):
         n=filename.replace("_collated_pred.tsv.gz", "")
         input_file = gzip.open(os.path.join(path, filename), 'rb')
         #chr1    1366120 1366360 ENSG00000179403 0.9632734656333923      E10
         for line in input_file:
            li = line.decode().rstrip('\n').split('\t')
            chr, start, end, gene,score = li[0],int(li[1]),int(li[2]),li[3],li[4]
            chr = chr.replace("chr", "")
            if gene in GENE:
                g=GENE[gene]
            else:
                next
            for num in range(start,end+1):
                id = chr+':'+str(num)
                if id in POS_grch37 and n in CELL:
                    epi.setdefault(id,[]).append(CELL[n])

    for key, value in epi.items():
        map[key] = list(set(value))
    return map

def anno_HiC_gen(path):
    input_file = open(path, "r")
    new_dict={}
    hic={}
    HIC={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        chr, start, end = li[0],li[1],li[2]
        for num in range(3,len(li)):
            if line.startswith('Hic_chr'):
                new_dict.setdefault(num,li[num])
            else:
                for nu in range(int(start),int(end)+1):
                   id = chr+':'+str(nu)
                   if id in POS_grch37 and li[num] != "0" and new_dict[num] in CELL:
                          for i in li[num].split('|'):
                              j=i.split('[')[0]
                              for k in j.split(';'):
                                      hic.setdefault(id, {}).setdefault(k, []).append(CELL[new_dict[num]])
    for k, v in hic.items():
        if type(v) is dict:
            for nk, nv in v.items():
                HIC.setdefault(k,{})[nk]=list(set(hic[k][nk]))
    return HIC        

def anno_HiC(path):
    input_file = open(path, "r")
    new_dict={}
    hic={}
    HIC={}
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        chr, start, end = li[0],li[1],li[2]
        for num in range(3,len(li)):
            if line.startswith('Hic_chr'):
                new_dict.setdefault(num,li[num])
            else:
                for nu in range(int(start),int(end)+1):
                   id = chr+':'+str(nu)
                   if id in POS_grch37 and li[num] != "0" and new_dict[num] in CELL:
                        hic.setdefault(id,[]).append(CELL[new_dict[num]])
    for key, value in hic.items():
        HIC[key] = list(set(value))
    return HIC


if gene == "gene_no":
            DICE=qtl_ano("./eQTL/DICE/DICE_all_cell.txt")
            IMMU=qtl_ano("./eQTL/ImmuNexUT/E-GEAD-398/ImmuNexUT_all_cell_conditional_eQTL_FDR0.05.tsv")
            BLUE=qtl_ano("./eQTL/BLUEPRINT/QTL/EGAD00001005200-significant_association/blueprint_eqtl.txt")
            GTEX=gtex_ano("./eQTL/Gtex/GTEx_Analysis_v8_eQTL/GTEx_select_cell.txt")
            QTLGEN=qtl_ano("./eQTL/eQTLGen/all_cis-eQTLs.txt")
            SCRNA=qtl_ano("./eQTL/scRNA-seq/Science2022-Single-cell_autoimmune_disease/scRNA_autoimmune_disease_eqtl.txt")

            ABC_MODEL=anno_ABC('./eQTL/Activity-By-Contact/ABC0.015_immune_cell')
            EpiMap=anno_EpiMap('./eQTL/EpiMap/hg19')
#            Captur_HIC=anno_HiC('./eQTL/Capture_Hi-C/hg19/OSF/PCHiC_peak_matrix_cutoff5_combine_lg5.txt')
elif gene == "gene_yes":
            DICE=qtl_ano_gene("./eQTL/DICE/DICE_all_cell.txt")
            IMMU=qtl_ano_gene("./eQTL/ImmuNexUT/E-GEAD-398/ImmuNexUT_all_cell_conditional_eQTL_FDR0.05.tsv")
            BLUE=qtl_ano_gene("./eQTL/BLUEPRINT/QTL/EGAD00001005200-significant_association/blueprint_eqtl.txt")
            GTEX=gtex_ano_gene("./eQTL/Gtex/GTEx_Analysis_v8_eQTL/GTEx_select_cell.txt")
            QTLGEN=qtl_ano_gene("./eQTL/eQTLGen/all_cis-eQTLs.txt")
            SCRNA=qtl_ano_gene("./eQTL/scRNA-seq/Science2022-Single-cell_autoimmune_disease/scRNA_autoimmune_disease_eqtl.txt")
      
            ABC_MODEL=anno_ABC_gen('./eQTL/Activity-By-Contact/ABC0.015_immune_cell')
            EpiMap=anno_EpiMap_gen('./eQTL/EpiMap/hg19')
#            Captur_HIC=anno_HiC_gen('./eQTL/Capture_Hi-C/hg19/OSF/PCHiC_peak_matrix_cutoff5_combine_lg5.txt')
def decide(a,b,c):
    d=''
    if a in c:
        d=';'.join(c[a])
    elif b in c:
        d=';'.join(c[b])
    else:
        d= 0
    return d

def decide_pos(a,b):
    d=''
    if a in b:
        d=';'.join(b[a])
    else:
        d= 0
    return d

def decide_gen(a,b,c,d):
    e=''
    if a in d and c in d[a]:
        e=';'.join(d[a][c])
    elif b in d and c in d[b]:
        e=';'.join(d[b][c])
    else:
        e= 0
    return e
def decide_pos_gen(a,b,c):
    d=''
    if a in c and b in c[a]:
        d=';'.join(c[a][b])
    else:
        d= 0
    return d

def annotation_file():
    global inputfile
    global outputfile
    global gene
    output = open(outputfile, "w",encoding="utf-8")
    input_file = open(inputfile, "r")
    for line in input_file:
        li = line.rstrip('\n').split('\t')
        if line.startswith('#LOCI'):
            head='DICE\tImmuNexUT\tBLUEPRINT\tGTEx\teQTLGen\tscRNA-seq_of_autoimmune_disease\teqtl_total_cell\teqtl_Cell_count\tABC_model\tEpiMap\tAll_Cell_unique\tAll_cell_count\tSupport_evidence'
            output.write(line.rstrip('\n')+'\t'+str(head)+'\n')
            continue
        id1, id2, grch37, grch38, gene_name = li[3],li[6],li[5],li[4],li[7]
        if gene == "gene_no":
            dice=decide(id1,id2,DICE)
            immu=decide(id1,id2,IMMU)
            blue=decide(id1,id2,BLUE)
            gtex=decide_pos(grch38,GTEX)
            eqtlgen=decide(id1,id2,QTLGEN)
            scrna=decide(id1,id2,SCRNA)

            abc_model=decide_pos(grch37,ABC_MODEL)
            epimap=decide_pos(grch37,EpiMap)
#            cap_hic=decide_pos(grch37,Captur_HIC)
        elif gene == "gene_yes":
            dice=decide_gen(id1,id2,gene_name,DICE)
            immu=decide_gen(id1,id2,gene_name,IMMU)
            blue=decide_gen(id1,id2,gene_name,BLUE)
            gtex=decide_pos_gen(grch38,gene_name,GTEX)
            eqtlgen=decide_gen(id1,id2,gene_name,QTLGEN)
            scrna=decide_gen(id1,id2,gene_name,SCRNA)

            abc_model=decide_pos_gen(grch37,gene_name,ABC_MODEL)
            epimap=decide_pos_gen(grch37,gene_name,EpiMap)
#            cap_hic=decide_pos_gen(grch37,gene_name,Captur_HIC)
            
        # unique cell types in above QTL studies
        qtl_all=[]
        for k in dice,immu,blue,gtex,eqtlgen,scrna:
            if k==0:
                continue
            else:
                for i in k.split(';'):
                    qtl_all.append(i)
        l=list(set(qtl_all))
        if len(l) == 0:
           qtl_cell=0
           qtl_count=0
        else:
           qtl_cell=';'.join(l)
           qtl_count=len(l)
        
        cell_all=[]
        flag=0
#        for g in qtl_cell,abc_model,epimap,cap_hic:
        for g in qtl_cell,abc_model,epimap:
            if g==0:
                continue
            else:
                flag+=1
                for i in g.split(';'):
                    cell_all.append(i)
        ll=list(set(cell_all))
        if len(ll) == 0:
           all_cell=0
           all_count=0
        else:
           all_cell=';'.join(ll)
           all_count=len(ll)
        output.write(line.rstrip('\n')+'\t'+str(dice)+'\t'+str(immu)+'\t'+str(blue)+'\t'+str(gtex)+'\t'+str(eqtlgen)+'\t'+str(scrna)+'\t'+str(qtl_cell)+'\t'+str(qtl_count)+'\t'+str(abc_model)+'\t'+str(epimap)+'\t'+str(all_cell)+'\t'+str(all_count)+'\t'+str(flag)+'\n')

annotation_file()
