#!/usr/bin/env python 

VERSION = 2.0
BUILD = 109

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
		   '     -w  [network file]				Enable Network boosting. Tab-delmited file of edges. [GeneA (\\t) GeneB]'
		   ''
		   '\n'
		   '  other options:\n'
		   '     -b, --bootstrapping            Use bootstrapping instead of cross-validation (Slow)\n'
		   '     --numcv=N                      Number of sections for cross validation (default 10)\n'
	       '     --numiter=N                    Number of bootstrap iterations (default 1000)\n'
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
		   '     -c  [control column]           Control (T0 or plasmid) column\n'
		   '\n' 
		   '  other options:\n'
		   '     --minreads=N                   Discard gRNA with T0 counts < N (default 0)\n'
		   '     --pseudo=N	                    Add a pseudocount of N to every readcount (default 5)\n'
		   '     -h, --help                     Show this help text\n'
		   '\n'
		   '  Example:\n' 
		   '  BAGEL.py fc -i readcount_file -o experiment_name -c 1\n' 
		   '\n'
		   '  Filters readcount_file for reagents with at least 30 reads in the control sample,\n'
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
check=0

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
			print helptext
			sys.exit()
		elif opt == '-i':
			readcountfile = arg
		elif opt == '-o':
			label = arg
		elif opt == '-c':
			ctrl_column = int(arg)
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
	
	control_label = reads.columns.values[ctrl_column]
	numClones, numColumns = reads.shape
	
	#
	# missing gene name = replace
	# missing read count = zero count
	#
	reads[ reads.columns.values[0] ].fillna('NO_GENE_NAME', inplace=True)
	reads.fillna(0, inplace=True)
	
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
	f = where( reads.ix[:,ctrl_column ] >= MIN_READS )[0]
	normed = normed.ix[f,:]
	
	#
	# calculate fold change
	#
	foldchange = pd.DataFrame( index=normed.index.values )
	foldchange.index.name = 'REAGENT_ID'
	foldchange['GENE'] = reads.ix[f,0]				# dataframe 'normed' has no GENE column
	for i in range( numColumns -1 ):			
		foldchange[ normed.columns.values[i] ] = log2( (normed.ix[:,normed.columns.values[i] ])   / normed.ix[:,control_label] )
	#
	# we have calculated a foldchange for the control column.  Drop it.
	#
	foldchange.drop( control_label, axis=1, inplace=True)
	
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

	try:
		opts, args = getopt.getopt(sys.argv[2:], "hti:o:c:e:n:w:b", ["numcv=","numiter=","help","bootstrapping"])
	except getopt.GetoptError:
		print helptext('bf')
		sys.exit(2)
	if len(opts) == 0:
		print  helptext('bf')
		sys.exit(2)
	for opt, arg in opts:
		if opt in ( '-h', '--help'):
			print helptext
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
		elif opt == '-w':
			NETWORKBOOST = True
			print "Network boosting enabled"
			networkfile = arg
		elif opt in ('-b','--bootstrapping'):
			TRAINMETHOD = 0
		elif opt in ('-t'):
			TESTMODE = True
		else:
			print helptext('bf')
			print "Error! Unknown arguments"
			sys.exit(2)
	try:
		outfilename
	except:
		outfilename = foldchangefile.replace(" ","_") + ".bf"
	if check!=3:
		print helptext('bf')
		print "Error! Missing arguments"
		sys.exit(2)
	from numpy import *
	import scipy.stats as stats
	import pandas as pd
	import random as rd
	
	
	column_list = [int(c) for c in columns]
	
	genes={}
	fc = {}
	
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
				select = rd.randrange(len(self._bucket))
				drawing.append(self._bucket[select])
				mask[self._bucket[select]] = False
				self._bucket = delete(self._bucket, select)
			if self._step < self._n % self._cvnum: # for distribute remain..
				select = rd.randrange(len(self._bucket))
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
	# LOAD FOLDCHANGES
	#

	with open(foldchangefile) as fin:
		skipfields = fin.readline().rstrip().split('\t')
		for i in column_list:
			print "Using column:  " + skipfields[i+1]
		for line in fin:
			fields = line.rstrip().split('\t')
			gsym = fields[1]
			genes[ gsym ]=1
			if ( not gsym in fc ):
				fc[gsym]=[]    # initialize dict entry as a list
			for i in column_list:
				fc[gsym].append( float(fields[i + 1]))		# per user docs, GENE is column 0, first data column is col 1.
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
	FC_THRESH = 2**(-1.1535*log(len(in1d(genes_array,nonEss))+13.324) + 0.7728)  #Define foldchange dynamic threshold. logarithm decay. parameters are defined by regression (achilles data)  2**-7 was used in previous version.
	bf = {}
	boostedbf = {}
	for g in genes_array:
		bf[g]=[]
		boostedbf[g] = []
	
	#
	# TRAINING
	#
	
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
		ess_train_fc_list_of_lists = [ fc[x] for x in genes_array[gene_train_idx[train_ess]] ]
		ess_train_fc_flat_list = [obs for sublist in ess_train_fc_list_of_lists for obs in sublist]
		#
		# define non_train vector of observed fold changes of nonessential genes in training set
		#
		non_train_fc_list_of_lists = [ fc[x] for x in genes_array[gene_train_idx[train_non]] ]
		non_train_fc_flat_list = [obs for sublist in non_train_fc_list_of_lists for obs in sublist]
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
			for foldchange in array(fc[g]):
				if foldchange >= xmin and foldchange <= xmax:
					testx.append(round(foldchange*100)/100)  
					testy.append(logratio_lookup[round(foldchange*100)][0])
					
		slope, intercept, r_value, p_value, std_err = stats.linregress(array(testx),array(testy))
		#
		# BF calculation
		#
		for g in genes_array[gene_test_idx]:
			foldchanges = array( fc[g] )
			bayes_factor = 0.0
			for x in foldchanges:
				bayes_factor += slope * x + intercept
			bf[g].append(bayes_factor)
		
	if TESTMODE == True:
		fp.close()
	bf_mean = dict()
	bf_std = dict()
	bf_norm = dict()
	
	for g in sorted( bf.keys() ):
		num_obs = len( bf[g] )
		bf_mean[g] = mean( bf[g] )
		bf_std[g]  = std( bf[g] )
		if bf_std[g] == 0.0:
			bf_norm[g] = ( bf[g] - bf_mean[g] )
		else:
			bf_norm[g] = ( bf[g] - bf_mean[g] ) / bf_std[g]
	
	training_data = Training(gene_idx)  # set training class reset
	
	
	#
	# calculate network scores
	#

	if NETWORKBOOST == True:
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
	
		if TESTMODE == True:
			fp.close()
	#
	# print out results
	#

	fout = open(outfilename, 'w')
	if NETWORKBOOST == True:
		fout.write('GENE\tBoostedBF\tSTD_BoostedBF\tBF\tSTD\tNumObs\n')
	else:
		fout.write('GENE\tBF\tSTD\tNumObs\n')
	for g in sorted( bf.keys() ):
		num_obs = len( bf[g] )
		
		if NETWORKBOOST == True:
			boostedbf_mean = mean( boostedbf[g] )
			boostedbf_std  = std( boostedbf[g] )
			fout.write('{0:s}\t{1:4.3f}\t{2:4.3f}\t{3:4.3f}\t{4:4.3f}\t{5:d}\n'.format( g, float(boostedbf_mean), float(boostedbf_std), float(bf_mean[g]), float(bf_std[g]), num_obs ) )
		else:
			fout.write('{0:s}\t{1:4.3f}\t{2:4.3f}\t{3:d}\n'.format( g, float(bf_mean[g]), float(bf_std[g]), num_obs ) )
	fout.close()



#-------------------------------------------#
#											#
#    CALCULATE PRECISION-RECALL CURVES      #
#											#
#-------------------------------------------#

elif sys.argv[1] == 'pr':
	try:
		opts, args = getopt.getopt(sys.argv[2:], "hti:o:c:e:n:w:b", ["numiter=","help","bootstrapping"])
	except getopt.GetoptError:
		print helptext('pr')
		sys.exit(2)
	if len(opts) == 0:
		print  helptext('pr')
		sys.exit(2)
	for opt, arg in opts:
		if opt in ( '-h', '--help'):
			print helptext
			sys.exit()
		elif opt == '-i':
			bf_file = arg
			check+=1
		elif opt == '-o':
			outfilename = arg
		elif opt == '-e':
			check+=1
			ess_ref = arg
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
	fout = open(outfilename, 'w')

	bf.sort_values(by='BF', ascending=False, inplace=True)

	cumulative_tp = 0.
	cumulative_fp = 0.
	precision = 1.
	recall = 0.
	# note float formats

	ess = essentials.index.values
	non = nonessentials.index.values
	totNumEssentials = len( [x for x in bf.index.values if x in ess] )


	fout.write('Gene\tBF\tRecall\tPrecision\n')

	for g in bf.index.values:
	    if ( g in ess ):
	        cumulative_tp += 1
	    elif (g in non):
	        cumulative_fp += 1
	    recall = cumulative_tp / totNumEssentials
	    if ( (cumulative_tp>0) | ( cumulative_fp > 0) ):
	        precision = cumulative_tp / (cumulative_tp + cumulative_fp)
	    fout.write('{0:s}\t{1:4.3f}\t{2:4.3f}\t{3:4.3f}\n'.format( g, bf.ix[g,'BF'], recall, precision) )
		
	fout.close()

else:
	print helptext_main
	sys.exit(2)

