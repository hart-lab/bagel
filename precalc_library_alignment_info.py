#!/usr/bin/python
'''
## Precaluclation of library alignment info ##

Usage) python3 thisscript.py [READID2SEQ] [GENOME (bowtie1 indexed)] [GENECODE gtf]
python3 precalc_library_alignment_info.py Avana_readid2seq.tsv hg19 gencode.v28lift37.annotation.gtf


'''
import sys

outputfastq="temp_bowtie.fq"
outputbowtie="temp_bowtie_output"
outputalignment="Alignment_info.txt"
bowtiepath="bowtie" # "bowtie" = default
nomenclature="gene_id" # "gene_name" for symbol "gene_id" for ensembl geneid
threads = 6 #multicores

helptext=F"""## Precaluclation of library alignment info ##

Usage) python3 precalc_library_alignment_info.py [READID2SEQ] [GENOME (bowtie1 indexed)] [GENECODE gtf] (options)
python3 precalc_library_alignment_info.py Avana_readid2seq.tsv hg19 gencode.v28lift37.annotation.gtf

  - [READID2SEQ] format
    READID [tab] GUIDESEQ
    READ1  [tab] ATAGATGTCCTGTGGCCCCG
    ...

    If readid is same with guideseq,
    READID [tab] GUIDESEQ
    ATAGATGTCCTGTGGCCCCG [tab] ATAGATGTCCTGTGGCCCCG
    ...    

  - Options
    --output [path]       : Specify an alignment-info file (Default:{outputalignment})
    --outputfastq [path]  : Specify a temporary fastq output file (Default:{outputfastq})
    --outputbowtie [path] : Specify a temporary bowtie output file (Default:{outputbowtie})
    --custompam [String]  : Specify a PAM sequence (Default: NGG)
    --pam-loc [number]    : Loctaion of PAM (0: End (Default, Cas9), 1: Front (for Cpf1), 2: No PAM)
    
"""

if len(sys.argv) < 4:
    print(helptext)
    sys.exit(1)

import numpy as np
import scipy.stats as stats
import pandas as pd
import subprocess

readid2seq_file = sys.argv[1]
genomefile = sys.argv[2]
gencodefile = sys.argv[3]

try:
    readid2seq = pd.read_csv(readid2seq_file,sep="\t",index_col=0,header=0).iloc[:,0]
except:
    print("Error - Check READID2SEQ file format")
    sys.exit(1)

# check duplicates of readid
if readid2seq.index.is_unique != True:
    print("Error - Read ID is not unique")
    sys.exit(1)

pam="NGG"
pamloc=0

import getopt
try:
    opts, args = getopt.getopt(sys.argv[4:], "", ["outputfastq=","outputbowtie=","output=","pam-loc=","custompam="])
except getopt.GetoptError:
    print(helptext)
    sys.exit(1)
for opt, arg in opts:
    if opt == '--output':
        outputalignment = arg
    elif opt == '--outputfastq':
        outputfastq = arg
    elif opt == '--outputbowtie':
        outputbowtie = arg

    elif opt == '--custompam':
        pam = str(arg)
    elif opt == '--pam-loc':
        pamloc=int(arg)

if pamloc == 2: # no pam
    pam = ""



def gencode_parsing(genecodefile,chrlist=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY",'chrM']):
    ## Human GENCODE parsing
    
    chr2exon = {chrom:list() for chrom in chrlist}
    with open(genecodefile,'r') as fp:
        for line in fp:
            if line[0] == '#':
                continue
            linearray = line.rstrip().split("\t")
            #print linearray  # ['chr1', 'HAVANA', 'gene', '3073253', '3074322', '.', '+', '.', 'gene_id "ENSMUSG00000102693.1"; gene_type "TEC"; gene_name "RP23-271O17.1"; level 2; havana_gene "OTTMUSG00000049935.1";']
            chrom = linearray[0]
            if chrom not in chrlist:
                continue
            featuretype = linearray[2]
            start = int(linearray[3])
            end = int(linearray[4])
            strand = linearray[6]
            tags = [x.strip().replace('"','').split(' ') for x in linearray[8].split(";")]
            
            if featuretype=='exon':
                gene=""
                transcript=""
                support_lv=6
                for tag in tags:
                    if tag[0] == nomenclature:#'gene_id':
                        gene=tag[1].split(".")[0]
                    elif tag[0] == 'transcript_id':
                        transcript = tag[1]
                    elif tag[0] == 'exon_number':
                        exon = int(tag[1])
                    elif tag[0] == 'transcript_support_level':
                        if tag[1]=='NA':
                            tag[1] = 6  # treat as max to ignore
                        support_lv = int(tag[1])
                chr2exon[chrom].append({"start":start,"end":end,"exon_number":exon,"transcript":transcript,"gene":gene})

    for chrom in chr2exon:
        chr2exon[chrom] = pd.DataFrame().from_dict(chr2exon[chrom])
        print(chr2exon[chrom])
    return chr2exon

def find_genes(chr2exon,chrom,location,direction):
    
    if direction == "+":
        start = location
        end = location+20
    elif direction == "-":
        start = location+3
        end = location+23
    else:
        start = location
        end = location
    
    return chr2exon[chrom][(end >= chr2exon[chrom]['start']) & (start <= chr2exon[chrom]['end'])]['gene'].unique()


def prepare_library_alignment_info(guideid2seq_series,genome,gencode,outputfastq="guideseq_bowtie.fq",outputbowtie="bowtie_output.txt",outputalignment="Alignment_info.txt"
    ,chrlist=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY",'chrM'],pam="NGG",pamloc=0):
    
    
    chr2exon = gencode_parsing(gencode)
    # write guide_sgrna file (New version. Using sequence instead of readid)
    print("Gencode Parsing Completed!")

    ambiguous_nt = pam.count("N")
    check_duplicate = set()
    with open(outputfastq,'w') as fout:
        for readid in guideid2seq_series.index:

            rna = guideid2seq_series[readid]
            if pamloc==1: # cpf1
                newrna = pam + rna
            else: # cas9 or nopam
                newrna = rna + pam
            
            fout.write("@%s\n"%readid)
            fout.write("%s\n"%(newrna))
            fout.write("+\n")
            fout.write("%s\n"%("E"*len(newrna)))
            



    
    # run bowtie
    try:
        
        #subprocess.run(["bowtie","-p",str(threads),"-v",str(2+ambiguous_nt),"-l",str(5),"-a",genome,outputfastq,outputbowtie])
        subprocess.run([bowtiepath,"-p",str(threads),"-v",str(2+ambiguous_nt),"-l",str(5),"-a",genome,outputfastq,outputbowtie])

        
    except:
        print("Error - fail to run bowtie")
        sys.exit(1)
    

    
    # check alignment
    mismatch = dict()
    library_align = pd.read_table(outputbowtie,header=None,index_col=None)
    print("Load bowtie result.")
    count=0
    for i in range(len(library_align.index)):
        count+=1
        if count%100000 == 0:
            print(count)
        if pd.isnull(library_align.iloc[i,7]) == False:
            try:
                mismatch[i] = len(library_align.iloc[i,7].split(",")) - ambiguous_nt  # substract the mismatch at NGG
            except:
                print(f"'{library_align.iloc[i,7]}'")
        else:
            mismatch[i] = 0
        
    library_align[8] = pd.Series(mismatch)

    library_align.head(10)

    print("Start mapping")
    read2align = dict()
    count=0
    for i in library_align.index:
        count+=1
        if count%1000 == 0:
            print (count)
        readid = library_align[0][i]
        direction = library_align[1][i]
        chrom = library_align[2][i]
        if chrom not in chrlist:
            continue
        location = library_align[3][i]
        mismatch = library_align[8][i]
        if readid not in read2align:
            read2align[readid] = {0:{},1:{},2:{}}  # key = mismatch
        locationtag = "%s_%d_%s"%(chrom,location,direction)
        read2align[readid][mismatch][locationtag] = list()

        # find genes overlap with sgRNAs
        results = find_genes(chr2exon,chrom,location,direction)

        read2align[readid][mismatch][locationtag] = list(results)
        

    print("Exporting results.")
        
    with open(outputalignment,'w') as fout:
        for readid in read2align:
            printstr = readid
            for mismatch in [0,1,2]:
                genes = set()
                for position in read2align[readid][mismatch]:
                    genes |=  set(read2align[readid][mismatch][position])

                printstr += "\t%s\t%s"%(",".join(read2align[readid][mismatch].keys()),",".join(list(genes)))
            fout.write(printstr+"\n")

            






prepare_library_alignment_info(readid2seq,genomefile,gencodefile,outputfastq=outputfastq,outputbowtie=outputbowtie,outputalignment=outputalignment,pam=pam,pamloc=pamloc)

print(f"Job completed - {outputalignment}")






prepare_library_alignment_info(readid2seq,genomefile,gencodefile,outputfastq=outputfastq,outputbowtie=outputbowtie,outputalignment=outputalignment,pam=pam,pamloc=pamloc)

print(f"Job completed - {outputalignment}")
