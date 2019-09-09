#!/usr/bin/python

from pylab import *

import scipy.stats as stats
import pandas as pd
import subprocess


def gencode_parsing(genecodefile,chrlist=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY",'chrM']):
	## Human GENCODE parsing
	
	chr2gene = {x:set() for x in chrlist}
	gene2transcript = dict()
	transcript2gene = dict()
	transcript2exon = dict()
	transcript2type = dict()
	transcript2lv = dict()
	#transcript2CDS = dict()
	transcriptlength = dict()
	transcript2location = dict()
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
				for tag in tags:
					if tag[0] == 'gene_name':#'gene_id':
						gene=tag[1].split(".")[0]
					elif tag[0] == 'transcript_id':
						transcript = tag[1]
					elif tag[0] == 'exon_number':
						exon = int(tag[1])
					elif tag[0] == 'transcript_support_level':
						if tag[1]=='NA':
							tag[1] = 6  # treat as max to ignore
						support_lv = int(tag[1])
				if support_lv == 6:
					continue
				if transcript not in transcript2exon:
					transcript2exon[transcript] = dict()
				transcript2exon[transcript][exon] = (chrom,start,end)
				
			elif featuretype == 'transcript':
				gene=""
				transcript=""
				for tag in tags:
					if tag[0] == 'gene_name':#'gene_id':
						gene=tag[1].split(".")[0]
					elif tag[0] == 'transcript_id':
						transcript = tag[1]
					elif tag[0] == 'exon_number':
						exon = int(tag[1])
					elif tag[0] == 'transcript_support_level':
						if tag[1]=='NA':
							tag[1] = 6  # treat as max to ignore
						support_lv = int(tag[1])
					if tag[0] == 'tag' and 'appris_principal' in tag[1]:
						transcript2type[transcript] = tag[1]
				if support_lv == 6:
					continue
				transcript2lv[transcript] = support_lv
				transcript2location[transcript] = (chrom,start,end)
				transcriptlength[transcript] = abs(end-start)+1
				if gene not in gene2transcript:
					gene2transcript[gene]=set()
				gene2transcript[gene].add(transcript)
				transcript2gene[transcript] = gene
				chr2gene[chrom].add(gene)
	return gene2transcript,transcript2gene,transcript2exon,transcript2location,chr2gene

def check_exon_num_from_transcript(transcript2exon,transcript,chrom,location,direction):
	result = 0

	if direction == "+":
		start = location
		end = location+20
	elif direction == "-":
		start = location+3
		end = location+23
	else:
		start = location
		end = location
	for exon_num in transcript2exon[transcript]:
		if (transcript2exon[transcript][exon_num][1] <= start and start < transcript2exon[transcript][exon_num][2]) or (transcript2exon[transcript][exon_num][1] <= end and end < transcript2exon[transcript][exon_num][2]):
			result = exon_num

	return result


def find_genes(chr2gene,gene2transcript,transcript2location,chrom,location,direction):
	results = list()
	if direction == "+":
		start = location
		end = location+20
	elif direction == "-":
		start = location+3
		end = location+23
	else:
		start = location
		end = location
	for g in chr2gene[chrom]:
		for t in gene2transcript[g]:
			if (transcript2location[t][1] <=start and start < transcript2location[t][2]) or (transcript2location[t][1] <=end and end < transcript2location[t][2]):
				results.append((g,t))
	return results


def prepare_library_alignment_info(guideid2seq_series,genome,gencode,outputfastq="guideseq_bowtie.fq",outputbowtie="bowtie_output.txt",outputalignment="Alignment_info.txt"
	,chrlist=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY",'chrM']):
	
	
	gene2transcript,transcript2gene,transcript2exon,transcript2location,chr2gene = gencode_parsing(gencode)
	# write guide_sgrna file (New version. Using sequence instead of readid)

	with open(outputfastq,'w') as fout:
		for readid in guideid2seq_series:  #for i in range(len(librarydata.index)):
			#rna = librarydata[seq_col][i]
			#gene = librarydata[gene_col][i]

			rna = guideid2seq_series[readid]
			#if gene not in gene2seq:
			#	gene2seq[gene] = list()
			#gene2seq[gene].append(rna)
			fout.write("@%s\n"%readid)
			fout.write("%s\n"%(rna+"NGG"))
			fout.write("+\n")
			fout.write("%s\n"%("E"*len(rna+"NGG")))
			



	
	# run bowtie
	try:
		
		subprocess.run(["bowtie","-v",str(3),"-l",str(5),"-a",genome,outputfastq,outputbowtie])

		
	except:
		print("Error - fail to run bowtie")
		sys.exit(2)
	

	
	# check alignment
	mismatch = dict()
	library_align = pd.read_table(outputbowtie,header=None)


	for i in range(len(library_align.index)):
		mismatch[i] = len(library_align.ix[i][7].split(",")) - 1  # substract the mismatch at NGG
		
	library_align[8] = pd.Series(mismatch)

	library_align.head(10)


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
		results = find_genes(chr2gene,gene2transcript,transcript2location,chrom,location,direction)
		for (g,t) in results:
			exon = check_exon_num_from_transcript(transcript2exon,t,chrom,location,direction)
			read2align[readid][mismatch][locationtag].append((g,t,exon))
		

		
		
	with open(outputalignment,'w') as fout:
		for readid in read2align:
			printstr = readid
			for mismatch in [0,1,2]:
				genes = set()
				for position in read2align[readid][mismatch]:
					for (g,t,e) in read2align[readid][mismatch][position]:
						genes.add(g)
				printstr += "\t%s\t%s"%(",".join(read2align[readid][mismatch].keys()),",".join(list(genes)))
			fout.write(printstr+"\n")

			



### Prepare sgRNA librarys
LIBRARY_FILE = '/home/ekim8/export/project4/data/library_comparison/avana/guide_gene_map.csv'# TKO v3
librarydata = pd.read_csv(LIBRARY_FILE)


### check sequence uniqueness

print("Is unique?",librarydata.sgrna.is_unique)

### not unique -> make unique list
unique_sgrnas = librarydata.sgrna.unique()

readid2seq = pd.Series(index=unique_sgrnas,data=unique_sgrnas)


genome="/home/ekim8/export/bagelv2/otheralgorithms/ceres/bowtie_indexes/hg19"

gencode="/home/ekim8/export/project5/gencode.v28lift37.annotation.gtf"


prepare_library_alignment_info(readid2seq,genome,gencode,outputfastq="temp_bowtie.fq",outputbowtie="temp_bowtie_output",outputalignment="Alignment_info.txt")
