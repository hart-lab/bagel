#!/usr/bin/env python 

VERSION = 2.0
BUILD = 112

'''
Update history

Build 112
1. Add sgRNA filtering options

Build 111
1. Add an option to equalize # of sgRNA per gene

Build 110
1. Enable multi control for fold-change calculation
2. Now user can input column names
3. Fix Threshold function

'''

#---------------------------------
# BAGEL:  Bayesian Analysis of Gene EssentaLity
# (c) Traver Hart <traver@hart-lab.org>, Eiru Kim <rooeikim@gmail.com> 2017.
# modified 10/2017
# Free to modify and redistribute with attribtuion
#---------------------------------

# ------------------------------------
# constants
import sys, getopt

def helptext(arg):
	if (arg=='main'):
		return ('\n'
		'BAGEL.py' 
		'\n'
		'  from the Bayesian Analysis of Gene EssentiaLity (BAGEL) suite\n'
		'  Version ' + str(VERSION) + ' Build ' + str(BUILD) + '\n' 
		'\n'
		'Calculate fold changes from readcount data\n'
		'  BAGEL.py fc -i [read count file] -o [output label] -c [control column]\n'
		'\n'
		'Calculate Bayes Factors from foldchange data\n'
		'  BAGEL.py bf -i [fold change file] -o [output file] -e [reference essentials] -n [reference nonessentials] -c [columns to test]\n'
		'\n'
		'Calculate precision-recall from Bayes Factors\n'
		'  BAGEL.py pr -i [Bayes Factor file] -o [output file] -e [reference essentials] -n [reference nonessentials]\n'
		'\n')
	elif (arg=='bf' or arg=='analysis'):
		return( '\n'
		   'BAGEL.py bf -i [fold change file] -o [output file] -e [reference essentials] -n [reference nonessentials] -c [columns to test]\n' 
		   '\n'
		   '  from the Bayesian Analysis of Gene EssentiaLity (BAGEL) suite\n'
		   '  Version ' + str(VERSION) + ' Build ' + str(BUILD) + '\n' 
		   '\n'
		   '  required options:\n' 
		   '     -i  [fold change file]         Tab-delmited file of reagents and fold changes.  See documentation for format.\n' 
		   '     -o  [output file]              Output filename\n' 
		   '     -e  [reference essentials]     Training set of essential genes (text file)\n' 
		   '     -n  [reference nonessentials]  Training set of nonessential genes (text file)\n' 
		   '     -c  [columns to test]          comma-delimited list of columns in input file to include in analyisis\n' 
		   '\n' 
		   '  network options\n'
		   '     -w  [network file]				Enable Network boosting. Tab-delmited file of edges. [GeneA (\\t) GeneB]\n'
		   ''
		   '\n'
		   '  multi-target guides filtering options \n'
		   '     --align-info  [file]           Input precalculated align-info file\n'
		   '     -m                             Enable filtering multi-targeting guide RNAs\n'
		   '       --m0  [Number]                Filtering guide RNAs without mismatch targeting over than [N] loci (default = 10)\n'
		   '       --m1  [Number]                Filtering guide RNAs with 1-bp mismatch targeting over than [N] loci (default = 10)\n'
		   '\n'
		   '  other options:\n'
		   '     -b, --bootstrapping            Use bootstrapping instead of cross-validation (Slow)\n'
		   '     -s, --small-sample             Low-fat BAGEL, Only resampled training set (Bootstrapping, iteration = 100)\n'
		   '     -r                             Calculate sgRNA-wise Bayes Factor\n'
		   '     -f  [Number]                   Equalize the number of sgRNAs per gene\n'
		   '     --numcv=N                      Number of sections for cross validation (default 10)\n'
		   '     --numiter=N                    Number of bootstrap iterations (default 1000)\n'
		   '     --seed=N                       Define random seed\n'
		   '     -h, --help                     Show this help text\n'
		   '\n'
		   '  Example:\n' 
		   '  BAGEL.py bf -i foldchange_file -o experiment.bf -e essentials_training_set -n nonessentials_training_set -c 1,2,3\n'
		   '\n'
		   '  Calculates a log2 Bayes Factor for each gene; positive BFs indicate confidence that the gene is essential.\n'
		   '  writes to [output file]: gene name, mean Bayes Factor across all iterations, std deviation of BFs, and number of iterations\n'
		   '  in which the gene was part of the test set (and a BF was calculated[output file]\n' 
		   '\n')
	elif (arg=='fc'):
		return ('\n'
		   'BAGEL.py fc -i [read count file] -o [output label] -c [control column]\n' 
		   '\n'
		   '  from the Bayesian Analysis of Gene EssentiaLity (BAGEL) suite\n'
		   '  Version ' + str(VERSION) + '\n' 
		   '\n'
		   '  required options:\n' 
		   '     -i  [read count file]          Tab-delmited file of reagents and fold changes.  See documentation for format.\n' 
		   '     -o  [output label]             Label for all output files\n' 
		   '     -c  [control column]           comma-delimited list of columns of control (T0 or plasmid) columns (Either number or name)\n'
		   '\n' 
		   '  other options:\n'
		   '     --minreads=N                   Discard gRNA with T0 counts < N (default 0)\n'
		   '     --pseudo=N	                    Add a pseudocount of N to every readcount (default 5)\n'
		   '     -h, --help                     Show this help text\n'
		   '\n'
		   '  Example:\n' 
		   '  BAGEL.py fc -i readcount_file -o experiment_name -c 1\n' 
		   '\n'
		   '  calculates fold change, and writes [output label].foldchange and [output label].normalized_reads\n'
		   '\n')
	elif (arg=='pr'):
		return ('\n'
		   'BAGEL.py pr -i [Bayes factors] -o [output file] -e [reference essentials] -n [reference nonessentials]\n' 
		   '\n'
		   '  from the Bayesian Analysis of Gene EssentiaLity (BAGEL) suite\n'
		   '  Version ' + str(VERSION) + '\n' 
		   '  required options:\n' 
		   '     -i  [Bayes factors]            BAGEL output file.\n' 
		   '     -o  [output file]              Output filename\n' 
		   '     -e  [reference essentials]     File with list of training set of essential genes\n' 
		   '     -n  [reference nonessentials]  File with list of training set of nonessential genes\n' 
		   '\n'
		   '  other options:\n'
		   '     -k [column name]               Use other column (default \'BF\')'
		   '\n'
		   '  Example:\n' 
		   '  BAGEL.py pr -i input.bf -o output.PR -e ref_essentials -n ref_nonessentials\n' 
		   '\n')


if len(sys.argv) < 2:
	print helptext('main')
	sys.exit(2)

# ------------------------------------

import numpy as np
import pandas as pd
import scipy.stats as stats
from matplotlib.mlab import find


#-------------------------------------------#
#   SET CONTSTANTS; INITIALIZE VARIABLES    #
#-------------------------------------------#
	
# for bf:
NUMCV=10
NUM_BOOTSTRAPS = 1000
NETWORKBOOST = False
TRAINMETHOD = 1   # 0 == bootstrapping, 1 == cross-validation (default)
TESTMODE = False
SMALLSAMPLE = False
RNALEVEL = False
check=0
FLATSGRNA = False

# for fc:
MIN_READS = 0
pseudo = 5


#-------------------------------------------#
#											#
#   CALCULATE FOLD CHANGE FROM RAW READS    #
#											#
#-------------------------------------------#

if sys.argv[1] == 'fc':
	
	#----------------------------------#
	#   READ COMMAND LINE ARGUMENTS    #
	#----------------------------------#
	
	try:
		opts, args = getopt.getopt(sys.argv[2:], "hi:o:c:", ["minreads=","pseudo=","help"])
	except getopt.GetoptError:
		print helptext('fc')
		sys.exit(2)
	if len(opts) == 0:
		print helptext('fc')
		sys.exit(2)
	for opt, arg in opts:
		if opt in ( '-h', '--help'):
			print helptext('fc')
			sys.exit()
		elif opt == '-i':
			readcountfile = arg
		elif opt == '-o':
			label = arg
		elif opt == '-c':
			ctrl_columns = arg.split(",")
		elif opt == '--minreads':
			MIN_READS = int(arg)
		elif opt == '--pseudo':
			pseudo = float(arg)
		else:
			print helptext('fc')
			print "Error! Unknown arguments"
			sys.exit(2)

	#----------------------------------------------------------------#
	# Import raw read data, normalize, filter for T0 min readcounts  #
	# Output:   [output label].foldchange                            #
	#----------------------------------------------------------------#
	
	from numpy import *
	import scipy.stats as stats
	import pandas as pd
	
	reads = pd.read_table(readcountfile, sep='\t', index_col=0)
	

	#
	# missing gene name = replace
	# missing read count = zero count
	#
	reads[ reads.columns.values[0] ].fillna('NO_GENE_NAME', inplace=True)
	reads.fillna(0, inplace=True)
	
	#
	# check controls
	#

	try:
		try:
			ctrl_columns = map(int,ctrl_columns)
			ctrl_labels = reads.columns.values[ctrl_columns]
		except ValueError:
			ctrl_labels = ctrl_columns
			
		ctrl_sum = reads[ctrl_labels].sum(axis=1)
		reads.drop(ctrl_labels,axis=1,inplace=True)
		ctrl_label_new = ';'.join(ctrl_labels)
		reads[ctrl_label_new] = ctrl_sum
	except:
		print "Invalid controls"
		sys.exit(2)
	
	numClones, numColumns = reads.shape
	print "Controls: " + ", ".join(ctrl_labels)

	#
	# Add pseudo count
	#
	
	reads.ix[:,range(1,numColumns)] += pseudo
	
	#
	# normalize each sample to a fixed total readcount
	#
	sumReads = reads.ix[:,range(1,numColumns)].sum(0)
	normed   = pd.DataFrame( index=reads.index.values )
	normed['GENE'] = reads.ix[:,0]				# first column is gene name
	normed = reads.ix[:,range(1,numColumns)] / tile( sumReads, [numClones,1]) * 10000000	# normalize to 10M reads
	
	#
	# filter for minimum readcount
	#
	f = where( reads[ ctrl_label_new ] >= MIN_READS )[0]
	normed = normed.ix[f,:]
	
	#
	# calculate fold change
	#
	foldchange = pd.DataFrame( index=normed.index.values )
	foldchange.index.name = 'REAGENT_ID'
	foldchange['GENE'] = reads.ix[f,0]				# dataframe 'normed' has no GENE column
	for i in range( numColumns -1 ):			
		foldchange[ normed.columns.values[i] ] = log2( (normed.ix[:,normed.columns.values[i] ])   / normed[ctrl_label_new])
	#
	# we have calculated a foldchange for the control column.  Drop it.
	#
	foldchange.drop( ctrl_label_new, axis=1, inplace=True)
	
	#
	# write normed readcount file
	# write foldchange file
	#
	foldchange_filename = label + '.foldchange'
	foldchange.to_csv( foldchange_filename, sep='\t', float_format='%4.3f')
	
	normedreads_filename = label + '.normed_readcount'
	normed.to_csv( normedreads_filename, sep='\t', float_format='%3.2f')
	


#-------------------------------------------#
#											#
#  CALCULATE BAYES FACTOR FROM FOLDCHANGE   #
#											#
#-------------------------------------------#

elif sys.argv[1] in ['bf','analysis']:

	from numpy import *
	import scipy.stats as stats
	import pandas as pd
	
	import time
	seed = int(time.time() * 100000 % 100000)
	random.seed(seed) # set random seed

	MULTI_TARGET_FILTERING = False
	MULTI_TARGET_FILTERING_on = 10
	MULTI_TARGET_FILTERING_off = 10
	aligninfofile = ""

	try:
		opts, args = getopt.getopt(sys.argv[2:], "hti:o:c:e:n:w:f:mbsr", ["numcv=","numiter=","help","bootstrapping","small-sample","seed=","align-info=","m0=","m1="])
	except getopt.GetoptError:
		print helptext('bf')
		sys.exit(2)
	if len(opts) == 0:
		print  helptext('bf')
		sys.exit(2)
	for opt, arg in opts:
		if opt in ( '-h', '--help'):
			print helptext('bf')
			sys.exit()
		elif opt == '-i':
			foldchangefile = arg
			check+=1
		elif opt == '-o':
			outfilename = arg
		elif opt == '-e':
			check+=1
			ess_ref = arg
		elif opt == '-n':
			check+=1
			non_ref = arg
		elif opt == '-c':
			columns = arg.split(',')
		elif opt == '--numiter':
			NUM_BOOTSTRAPS = int(arg)
		elif opt == '--numcv':
			NUMCV = int(arg)
		elif opt == '--seed':
			seed = int(arg)
			random.seed(seed)
		elif opt == '-w':
			NETWORKBOOST = True
			print "Network boosting enabled"
			networkfile = arg
		elif opt in ('-b','--bootstrapping'):
			TRAINMETHOD = 0
		elif opt in ('-t'):
			TESTMODE = True
		elif opt in ('-s','--small-sample'):
			SMALLSAMPLE = True
			TRAINMETHOD = 0
			NUM_BOOTSTRAPS = 100
		elif opt == '-r':
			RNALEVEL=True
		elif opt == '-f':
			FLATSGRNA = True
			NUM_DESIRE_SGRNA = int(arg)
		elif opt == '-m':
			MULTI_TARGET_FILTERING = True
			print "Multi-target filtering enabled"
		elif opt in ["--m0"]:
			MULTI_TARGET_FILTERING_on = int(arg)
		elif opt in ["--m1"]:
			MULTI_TARGET_FILTERING_off = int(arg)
		elif opt == '--align-info':
			aligninfofile = arg



		else:
			print helptext('bf')
			print "Error! Unknown arguments"
			sys.exit(2)
	try:
		outfilename
	except:
		outfilename = foldchangefile.replace(" ","_") + ".bf"
	
	### Check arguments

	if check!=3:
		print helptext('bf')
		print "Error! Missing arguments"
		sys.exit(2)
	
	if MULTI_TARGET_FILTERING==True and aligninfofile == "":
		print helptext('bf')
		print "Error! Please indicate align-info file"
		sys.exit(2)

	if NETWORKBOOST == True and RNALEVEL==True:
		NETWORKBOOST = False
		print "# Network boosting is disabled in RNA-wise output"
	
	print 'Random seed = %d'%seed

	genes={}
	fc = {}
	gene2rna = {}
	rna2gene = {}
	
	def round_to_hundredth(x):
		return round( x*100) / 100.0
		
	def func_linear(x, a, b):
		return (a * x) + b
			
	class Training:
		def __init__(self,X, n=None, cvnum=10):
			if n==None:
				self._n = len(X)
			self._cvnum = cvnum
			self._bid = int(self._n/cvnum)
			self._bucket = arange(len(X))
			self._X = X
			self._step = 0
		def cross_validation(self):
			if self._bid < 1: #bid check
				print "The number of genes is too small! n<" + str(self._cvnum)
				sys.exit(2)
			drawing = list()
			mask = array([True]*self._n)
			for j in range(self._bid):
				#drawing.append(delete(self._bucket, random.randrange(len(self._bucket))))
				select = random.randint(len(self._bucket))
				drawing.append(self._bucket[select])
				mask[self._bucket[select]] = False
				self._bucket = delete(self._bucket, select)
			if self._step < self._n % self._cvnum: # for distribute remain..
				select = random.randint(len(self._bucket))
				drawing.append(self._bucket[select])
				mask[self._bucket[select]] = False
				self._bucket = delete(self._bucket, select)
			self._step+=1
			X_resample = self._X[mask]
			return X_resample, self._X[~mask]
		def get_cv_step(self):
			return self._step
		def bootstrap_resample(self):
			mask = array([False]*self._n)
			resample_i = floor(random.rand(self._n)*len(self._X)).astype(int)
	
			mask[resample_i] = True
			X_resample = self._X[mask]
			return X_resample, self._X[~mask]
		def get_data(self,METHOD=0):
			if METHOD == 0:
				train,test = self.bootstrap_resample()
			elif METHOD == 1:
				train,test = self.cross_validation()
			return train,test
	
	
	def fibo_weighted_sum(listofscore):
		value = p1 = p2 = 0.0
		c  = 1.0	#current value
		for v in listofscore:
			value += v / c
			p2 = p1   # go one step
			p1 = c
			c = p1 + p2
		return value

	#
	# LOAD ALIGN INFO FILE 
	#

	multi_targeting_sgrnas = dict()
	multi_targeting_sgrnas_info = dict()
	
	if MULTI_TARGET_FILTERING == True:
		from sklearn.linear_model import LinearRegression 
		try:
			aligninfo = pd.read_table(aligninfofile,header=None,index_col=0).fillna("")
			for seqid in aligninfo.index:
				perfectmatch = 0
				mismatch_1bp = 0
				perfectmatch_gene = 0
				mismatch_1bp_gene = 0
				if aligninfo[1][seqid] != "":
					perfectmatch = len(aligninfo[1][seqid].split(","))
				if aligninfo[2][seqid] != "":
					perfectmatch_gene = len(aligninfo[2][seqid].split(","))
				if aligninfo[3][seqid] != "":
					mismatch_1bp = len(aligninfo[3][seqid].split(","))
				if aligninfo[4][seqid] != "":
					mismatch_1bp_gene = len(aligninfo[4][seqid].split(","))
				if perfectmatch > MULTI_TARGET_FILTERING_on or mismatch_1bp > MULTI_TARGET_FILTERING_off:
					multi_targeting_sgrnas[seqid] = True
				elif perfectmatch > 1 or mismatch_1bp > 0:
					multi_targeting_sgrnas_info[seqid] = (perfectmatch, mismatch_1bp, perfectmatch_gene, mismatch_1bp_gene)

		except:
			print "Please check align-info file"
			sys.exit(2)

		print "Total %d multi-targeting gRNAs are discarded"%len(multi_targeting_sgrnas)		

	#
	# LOAD FOLDCHANGES
	#

	with open(foldchangefile) as fin:
		fieldname = fin.readline().rstrip().split('\t')
		#
		# DEFINE CONTROLS
		#
		try:
			try:
				column_list = map(int,columns)
				column_labels = [ fieldname[x+1] for x in column_list ]
			except ValueError:
				column_labels = columns
				column_list = [ x for x in range(len(fieldname) - 1) if fieldname[ x + 1] in column_labels]   # +1 because of First column start 2
			print "Using column:  " + ", ".join(column_labels)
			#print "Using column:  " + ", ".join(map(str,column_list))
			
		except:
			print "Invalid columns"
			sys.exit(2)

		for line in fin:
			fields = line.rstrip().split('\t')
			rnatag = fields[0]
			if MULTI_TARGET_FILTERING == True:  # multitargeting sgrna filtering
				if rnatag in multi_targeting_sgrnas:
					continue    # skip multitargeting sgrna.
			gsym = fields[1]

			genes[ gsym ]=1
			if gsym not in gene2rna:
				gene2rna[gsym]=[]
			gene2rna[gsym].append(rnatag)
			rna2gene[rnatag] = gsym
			fc[rnatag] = {}
			for i in column_list:
				fc[rnatag][i] = float(fields[i + 1])		# per user docs, GENE is column 0, first data column is col 1.


			
	genes_array = array( genes.keys() )
	gene_idx = arange( len( genes ) )
	print "Number of unique genes:  " + str( len(genes) )

	#
	# DEFINE REFERENCE SETS
	#
	coreEss = []
	
	with open(ess_ref) as fin:
		skip_header = fin.readline()
		for line in fin:
			coreEss.append( line.rstrip().split('\t')[0] )
	coreEss=array(coreEss)
	print "Number of reference essentials: " + str(len(coreEss))
	
	nonEss = []
	with open(non_ref) as fin:
		skip_header = fin.readline()
		for line in fin:
			nonEss.append( line.rstrip().split('\t')[0] )
	
	nonEss = array(nonEss)
	print "Number of reference nonessentials: " + str(len(nonEss))
	
	#
	# LOAD NETWORK 
	#
	
	if NETWORKBOOST == True:
		network = {}
		edgecount = 0
		with open(networkfile) as fin:
			for line in fin:
				 linearray = line.rstrip().split('\t') # GeneA \t GeneB format
				 if linearray[0] in genes_array and linearray[1] in genes_array:
					 for i in [0,1]:
						 if linearray[i] not in network:
							network[linearray[i]] = {}
						 network[linearray[i]][linearray[-1 * (i-1)]] = 1  # save edge information
					 edgecount += 1
				 
		print "Number of network edges: " + str(edgecount)
		
	
	#
	# INITIALIZE BFS
	#
	FC_THRESH = 2**(-1.1535*log(len(intersect1d(genes_array,nonEss))+13.324) + 0.7728)  #Define foldchange dynamic threshold. logarithm decay. parameters are defined by regression (achilles data)  2**-7 was used in previous version.
	bf = {}
	boostedbf = {}
	for g in genes_array:
		for rnatag in gene2rna[g]:
			bf[rnatag]=[]
				
		boostedbf[g] = [] # boosted bf at gene level
	
	#
	# TRAINING
	#
	if SMALLSAMPLE == True:
		#training_data = Training(setdiff1d(gene_idx,where(in1d(genes_array,coreEss))),cvnum=NUMCV)  # declare training class
		training_data = Training(where(in1d(genes_array,union1d(coreEss,nonEss)))[0],cvnum=NUMCV)  # declare training class (only for Gold-standard gene set)
		all_non_gs = where( logical_not( in1d(genes_array,union1d(coreEss,nonEss)) ) )[0] # all non-goldstandards
	else:
		training_data = Training(gene_idx,cvnum=NUMCV)  # declare training class
	
	if TRAINMETHOD == 0:
		LOOPCOUNT = NUM_BOOTSTRAPS
	elif TRAINMETHOD == 1:
		LOOPCOUNT = NUMCV  # 10-folds

	if TESTMODE == True:
		fp = open(outfilename+".traininfo","w")
		fp.write("#1: Loopcount\n#2: Training set\n#3: Testset\n")
		
	print "Iter",
	print "TrainEss",
	print "TrainNon",
	print "TestSet"
	sys.stdout.flush()
	for loop in range(LOOPCOUNT):
		currentbf = {}
		print str(loop),
		
		#
		# bootstrap resample (10-folds cross-validation) from gene list to get the training set
		# test set for this iteration is everything not selected in bootstrap resampled (10-folds cross-validation) training set
		# define essential and nonessential training sets:  arrays of indexes
		#

		gene_train_idx,gene_test_idx = training_data.get_data(TRAINMETHOD)
		if SMALLSAMPLE == True: # test set is union of rest of training set (gold-standard) and the other genes (all of non-gold-standard)
			gene_test_idx = union1d( gene_test_idx, all_non_gs  )

		if TESTMODE == True:
			fp.write("%d\n%s\n%s\n"%(loop,",".join(genes_array[gene_train_idx]),",".join(genes_array[gene_test_idx])))

		
		train_ess = where( in1d( genes_array[gene_train_idx], coreEss))[0]
		train_non = where( in1d( genes_array[gene_train_idx], nonEss))[0]
		print len(train_ess),
		print len(train_non),
		print len(gene_test_idx)
		sys.stdout.flush()
		#
		# define ess_train: vector of observed fold changes of essential genes in training set
		#
		ess_train_fc_list_of_lists = [ fc[rnatag] for g in genes_array[gene_train_idx[train_ess]] for rnatag in gene2rna[g] ]
		ess_train_fc_flat_list = [obs for sublist in ess_train_fc_list_of_lists for obs in sublist.values()]
		#
		# define non_train vector of observed fold changes of nonessential genes in training set
		#
		non_train_fc_list_of_lists = [ fc[rnatag] for g in genes_array[gene_train_idx[train_non]] for rnatag in gene2rna[g] ]
		non_train_fc_flat_list = [obs for sublist in non_train_fc_list_of_lists for obs in sublist.values()]
		#
		# calculate empirical fold change distributions for both
		#
		kess = stats.gaussian_kde( ess_train_fc_flat_list )
		knon = stats.gaussian_kde( non_train_fc_flat_list )
		#
		# define empirical upper and lower bounds within which to calculate BF = f(fold change)
		#
		x = arange(-10,2,0.01)
		nonfitx = knon.evaluate(x)
		# define lower bound empirical fold change threshold:  minimum FC where knon is above threshold
		f = where( nonfitx > FC_THRESH)
		xmin = round_to_hundredth( min(x[f]) )
		# define upper bound empirical fold change threshold:  minimum value of log2(ess/non)
		subx = arange( xmin, max(x[f]), 0.01)
		logratio_sample = log2( kess.evaluate(subx) / knon.evaluate(subx) )
		f = where( logratio_sample == logratio_sample.min() )
		xmax = round_to_hundredth( subx[f] )
		#
		# round foldchanges to nearest 0.01
		# precalculate logratios and build lookup table (for speed)
		#
		logratio_lookup = {}
		for i in arange(xmin, xmax+0.01, 0.01):
			logratio_lookup[round(i*100)] = log2( kess.evaluate(i) / knon.evaluate(i) )
		#
		# calculate BFs from lookup table for withheld test set
		#
		
		# liner interpolation
		testx=list()
		testy=list()

		for g in genes_array[gene_train_idx]:
			for rnatag in gene2rna[g]:
				for foldchange in fc[rnatag].values():
					if foldchange >= xmin and foldchange <= xmax:
						testx.append(round(foldchange*100)/100)  
						testy.append(logratio_lookup[round(foldchange*100)][0])
		try:
			slope, intercept, r_value, p_value, std_err = stats.linregress(array(testx),array(testy))
		except:
			print "Regression failed. Check quality of the screen"
			sys.exit(2)
		#
		# BF calculation
		#

		for g in genes_array[gene_test_idx]:
			for rnatag in gene2rna[g]:
				bayes_factor = []
				for rep in column_list:
					bayes_factor.append( slope * fc[rnatag][rep] + intercept )
				bf[rnatag].append(bayes_factor)


			
	if TESTMODE == True:
		fp.close()



	num_obs = dict()
	if RNALEVEL==False:
		bf_mean = dict()  
		bf_std = dict() 
		bf_norm = dict()  # sgRNA number complement
	if RNALEVEL==True or MULTI_TARGET_FILTERING == True:
		bf_mean_rna_rep = dict()
		bf_std_rna_rep = dict() 
		#bf_norm_rna_rep = dict()  
	
	
	
	for g in gene2rna:
		num_obs[g] = len(  bf[  gene2rna[g][0] ]  )
		if RNALEVEL==True or MULTI_TARGET_FILTERING == True:
			for rnatag in gene2rna[g]:
				bf_mean_rna_rep[rnatag] = dict()
				bf_std_rna_rep[rnatag] = dict()
				t = zip(*bf[rnatag])
				for rep in range(len(column_list)):
					bf_mean_rna_rep[rnatag][column_list[rep]] = mean(t[rep])
					bf_std_rna_rep[rnatag][column_list[rep]] = std(t[rep])

		if RNALEVEL==False:
			sumofbf_list = list()
			for i in range( num_obs[g] ):
				sumofbf = 0.0
				for rnatag in gene2rna[g]:
					sumofbf += sum(bf[rnatag][i])
				sumofbf_list.append(sumofbf)  # append each iter
			bf_mean[g] = mean(sumofbf_list)
			bf_std[g] = std(sumofbf_list)


	#
	# BUILD MULTIPLE REGRESSION MODEL FOR MULTI TARGETING GUIDE RNAs
	#
	if MULTI_TARGET_FILTERING == True:
		count=0
		trainset = dict()
		bf_multi_corrected_gene = dict()
		bf_multi_corrected_rna = dict()
		for gene in gene2rna:
			# multi_targeting_sgrnas_info[seqid] = (perfectmatch, mismatch_1bp, perfectmatch_gene, mismatch_1bp_gene)
			multitarget = list()
			onlytarget = list()
			for seqid in gene2rna[gene]:
				if seqid not in aligninfo.index:
					continue
				if seqid in multi_targeting_sgrnas_info:
					multitarget.append(seqid)
				else:
					onlytarget.append(seqid)

			if len(onlytarget) > 0: # comparsion between sgRNAs targeting one locus and multiple loci
				if len(multitarget) > 0:
					
					bf_only = mean([sum(bf_mean_rna_rep[seqid].values()) for seqid in onlytarget])
					for seqid in onlytarget:
						trainset[seqid] = [1,0,0]
		                
					for seqid in multitarget:
						if multi_targeting_sgrnas_info[seqid][2] > 1 or multi_targeting_sgrnas_info[seqid][3] > 0:  # train model using multi-targeting only targeting one protein coding gene
							continue
						
						count+=1
						increment = sum(bf_mean_rna_rep[seqid].values()) - bf_only

						trainset[seqid] = [multi_targeting_sgrnas_info[seqid][0], multi_targeting_sgrnas_info[seqid][1], increment]
		
		if count < 10:
			print("Not enough train set for calculating multi-targeting effect.\n")
			print("It may cause due to unmatched gRNA names between the foldchange file and the align info file.\n")
			print("Filtering is not finished\n")
			MULTI_TARGET_FILTERING = False

		else:

			trainset = pd.DataFrame().from_dict(trainset).T
			X = trainset[[0,1]]
			y = trainset[2]

			regressor = LinearRegression()  
			regressor.fit(X, y)  
			coeff_df = pd.DataFrame(regressor.coef_, X.columns, columns=['Coefficient'])  
			for i in [0,1]:
				if coeff_df['Coefficient'][i] < 0:
					print ("Regression coefficient is below than zero. Substituted to zero\n")
					coeff_df['Coefficient'][i] = 0.0
			print "Multiple effects from perfect matched loci = %.3f and 1bp mis-matched loci = %.3f"%(coeff_df['Coefficient'][0],coeff_df['Coefficient'][1])
			
			if RNALEVEL==False:
				for g in gene2rna:
					penalty = 0.0
					for seqid in gene2rna[gene]:
						if seqid in multi_targeting_sgrnas_info:
							penalty += float(multi_targeting_sgrnas_info[seqid][0] - 1) * coeff_df['Coefficient'][0] + float(multi_targeting_sgrnas_info[seqid][1]) * coeff_df['Coefficient'][1]
					bf_multi_corrected_gene[g] = bf_mean[g] - penalty
			else:
				for g in gene2rna:
					for seqid in gene2rna[g]:
						if seqid in multi_targeting_sgrnas_info:
							penalty = float(multi_targeting_sgrnas_info[seqid][0] - 1) * coeff_df['Coefficient'][0] + float(multi_targeting_sgrnas_info[seqid][1]) * coeff_df['Coefficient'][1]
						else:
							penalty = 0.0
						bf_multi_corrected_rna[seqid] = sum(bf_mean_rna_rep[seqid].values()) - penalty


	#
	#  NORMALIZE sgRNA COUNT
	#
	if RNALEVEL == False and FLATSGRNA == True:
		if MULTI_TARGET_FILTERING == True:
			targetbf = bf_multi_corrected_gene
		else:
			targetbf = bf_mean
			
		for g in gene2rna:
			multiple_factor = NUM_DESIRE_SGRNA / float(len(gene2rna[g]))
			bf_norm[g] = targetbf[g] * multiple_factor


	'''			
	if bf_std[rnatag] == 0.0:
		bf_norm[rnatag] = float('inf')
	else:
		bf_norm[g] = ( bf[rnatag] - bf_mean[rnatag] ) / bf_std[rnatag]
	'''	
	training_data = Training(gene_idx)  # set training class reset
	
	
	#
	# calculate network scores
	#

	if NETWORKBOOST == True and RNALEVEL==False:  # Network boost is only working for gene level	
		if TESTMODE == True: # TEST MODE
			fp = open(outfilename+".netscore","w")
		print "\nNetwork score calculation start\n"

		networkscores = {}
		for g in genes_array[gene_idx]:
			if g in network:
				templist = list()
				for neighbor in network[g]:
					if neighbor in bf_mean:
						templist.append( bf_mean[neighbor] )
					
				templist.sort(reverse = True)
				
				networkscores[g] = fibo_weighted_sum(templist)
		#
		# start training
		#

		for loop in range(LOOPCOUNT):	
			currentnbf = {}
			print str(loop),
			#
			# draw train, test sets
			#
			gene_train_idx,gene_test_idx = training_data.get_data(TRAINMETHOD)
			#
			# define essential and nonessential training sets:  arrays of indexes
			#
			train_ess = where( in1d( genes_array[gene_train_idx], coreEss))[0]
			train_non = where( in1d( genes_array[gene_train_idx], nonEss))[0]
			print len(train_ess),
			print len(train_non),
			print len(gene_test_idx)
			sys.stdout.flush()
			#
			# calculate Network BF for test set
			#
			ess_ns_list = [ networkscores[x] for x in genes_array[gene_train_idx[train_ess]] if x in networkscores]	
			non_ns_list = [ networkscores[x] for x in genes_array[gene_train_idx[train_non]] if x in networkscores]
		
			kess = stats.gaussian_kde( ess_ns_list )
			knon = stats.gaussian_kde( non_ns_list )
			#
			# set x boundary for liner regression
			#
			testx=list()
			testy=list()
			xmin = float(inf)
			xmax = float(-inf)
	
			for networkscore in arange(max(ess_ns_list),min(ess_ns_list),-0.01):
				density_ess = kess.evaluate(networkscore)[0]
				density_non = knon.evaluate(networkscore)[0]
				if density_ess == 0.0 or density_non == 0.0:
					continue
	
				if log2(density_ess / density_non) > -5 and networkscore < array(ess_ns_list).mean():  # reverse
					xmin = min(xmin,networkscore)
	
			for networkscore in arange(min(non_ns_list),max(non_ns_list),0.01):
				density_ess = kess.evaluate(networkscore)[0]
				density_non = knon.evaluate(networkscore)[0]
				if density_ess == 0.0 or density_non == 0.0:
					continue
				if log2(density_ess / density_non) < 5 and networkscore > array(non_ns_list).mean():  # reverse
					xmax = max(xmax,networkscore)
			#
			# liner regression
			#
			testx=list()
			testy=list()
			for g in genes_array[gene_train_idx]:
				if g in networkscores:
					if networkscores[g] >= xmin and networkscores[g] <= xmax:
						testx.append(round(networkscores[g]*100)/100)  
						testy.append(log2(kess.evaluate(networkscores[g])[0] / knon.evaluate(networkscores[g])[0]))
						
			slope, intercept, r_value, p_value, std_err = stats.linregress(array(testx),array(testy))
			
			for g in genes_array[gene_test_idx]:
				if g in networkscores:
					if TESTMODE == True:
						fp.write("%s\t%f\t%f\n"%(g,networkscores[g],  slope * networkscores[g] + intercept))
					nbf = slope * networkscores[g] + intercept
				else:
					nbf = 0.0
					
				boostedbf[g].append(bf_mean[g] + nbf)
				if FLATSGRNA == True:
					boostedbf[g].append(bf_norm[g] + nbf)
	
		if TESTMODE == True:
			fp.close()
	
	#
	# print out results
	#

	fout = open(outfilename, 'w')
	
	if RNALEVEL == True:
		fout.write('RNA\tGENE')
		for i in range(len(column_list)):
			fout.write('\t{0:s}'.format(column_labels[i]))
			if TRAINMETHOD == 0:
				fout.write('\t{0:s}'.format(column_labels[i]+"_STD"))
		fout.write('\tBF')
		if TRAINMETHOD == 0:	
			fout.write('\tNumObs')
		fout.write('\n')
		
		for rnatag in sorted( bf.keys() ):
			# RNA tag
			fout.write('{0:s}\t'.format(rnatag)) 
			# Gene
			gene = rna2gene[rnatag]
			fout.write('{0:s}\t'.format( gene) )

			# BF of replicates
			for rep in column_list:
				fout.write('{0:4.3f}\t'.format(bf_mean_rna_rep[rnatag][rep]))
				if TRAINMETHOD == 0:
					fout.write('{0:4.3f}\t'.format(bf_std_rna_rep[rnatag][rep]))
			# Sum BF of replicates
			if MULTI_TARGET_FILTERING == True:
				fout.write('{0:4.3f}'.format( float(bf_multi_corrected_rna[rnatag]) ))
			else:
				fout.write('{0:4.3f}'.format( float(sum(bf_mean_rna_rep[rnatag].values())) ))
			
			# Num obs
			if TRAINMETHOD == 0:
				fout.write('\t{0:d}'.format( num_obs[gene] ))
			fout.write('\n')
	else:
		fout.write('GENE')
		if NETWORKBOOST == True:
			fout.write('\tBoostedBF')
			if TRAINMETHOD == 0:
				fout.write('\tSTD_BoostedBF')
		fout.write('\tBF')
		if TRAINMETHOD == 0:
			fout.write('\tSTD\tNumObs')
		if FLATSGRNA == True:
			fout.write('\tNormBF')
		fout.write('\n')

		for g in sorted( genes.keys() ):
			# Gene
			fout.write('{0:s}'.format( g ))
			if NETWORKBOOST == True:
				boostedbf_mean = mean( boostedbf[g] )
				boostedbf_std  = std( boostedbf[g] )
				fout.write('\t{0:4.3f}'.format( float(boostedbf_mean) ))
				if TRAINMETHOD == 0:
					fout.write('\t{0:4.3f}'.format( float(boostedbf_std) ))

			# BF
			if MULTI_TARGET_FILTERING == True:
				fout.write('\t{0:4.3f}'.format( float(bf_multi_corrected_gene[g]) ))
			else:
				fout.write('\t{0:4.3f}'.format( float(bf_mean[g]) ))
			# STD, Count
			if TRAINMETHOD == 0:
				fout.write('\t{0:4.3f}\t{1:d}'.format( float(bf_std[g]), num_obs[g] ) )
			# Normalized BF
			if FLATSGRNA == True:
				fout.write('\t{0:4.3f}'.format( float(bf_norm[g]) ) )
			
			fout.write('\n')
	fout.close()



#-------------------------------------------#
#											#
#    CALCULATE PRECISION-RECALL CURVES      #
#											#
#-------------------------------------------#

elif sys.argv[1] == 'pr':
	BFCOL = 'BF'
	try:
		opts, args = getopt.getopt(sys.argv[2:], "htk:i:o:c:e:n:w:b", ["numiter=","help","bootstrapping"])
	except getopt.GetoptError:
		print helptext('pr')
		sys.exit(2)
	if len(opts) == 0:
		print  helptext('pr')
		sys.exit(2)
	for opt, arg in opts:
		if opt in ( '-h', '--help'):
			print helptext('pr')
			sys.exit()
		elif opt == '-i':
			bf_file = arg
			check+=1
		elif opt == '-o':
			outfilename = arg
		elif opt == '-e':
			check+=1
			ess_ref = arg
		elif opt == '-k':
			BFCOL = str(arg)
		elif opt == '-n':
			check+=1
			non_ref = arg
		else:
			print helptext('pr')
			print "Error! Unknown arguments"
			sys.exit(2)
	try:
		outfilename
	except:
		outfilename = foldchangefile.replace(" ","_") + ".pr"
	if check!=3:
		print helptext('pr')
		print "Error! Missing arguments"
		sys.exit(2)
	#
	# test for availability of all files
	#
	essentials = pd.read_table(ess_ref, index_col=0)
	nonessentials = pd.read_table(non_ref, index_col=0)
	bf = pd.read_table(bf_file, index_col=0)

	if BFCOL not in bf.dtypes.index:
		print "Error! the column name is not in the file"
		sys.exit(2)
		
	fout = open(outfilename, 'w')

	bf.sort_values(by=BFCOL, ascending=False, inplace=True)

	cumulative_tp = 0.
	cumulative_fp = 0.
	precision = 1.
	recall = 0.
	# note float formats

	ess = essentials.index.values
	non = nonessentials.index.values
	totNumEssentials = len( [x for x in bf.index.values if x in ess] )


	fout.write('Gene\t')
	fout.write(BFCOL)
	fout.write('\tRecall\tPrecision\n')

	for g in bf.index.values:
	    if ( g in ess ):
	        cumulative_tp += 1
	    elif (g in non):
	        cumulative_fp += 1
	    recall = cumulative_tp / totNumEssentials
	    if ( (cumulative_tp>0) | ( cumulative_fp > 0) ):
	        precision = cumulative_tp / (cumulative_tp + cumulative_fp)
	    fout.write('{0:s}\t{1:4.3f}\t{2:4.3f}\t{3:4.3f}\n'.format( g, bf.ix[g,BFCOL], recall, precision) )
		
	fout.close()

else:
	print helptext('main')
	sys.exit(2)

