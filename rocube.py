#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Joon An'
__date__ = 'September 29th, 2017'

description = '''

1. Aim
	Find a threshold of quality metrics for high-quality (HQ) rare variants from whole genome or whole exome data.

2. Input file
	The genotype matrix with GATK quality metrics is expected. The current version expects a matrix file generated from HAIL (export_genotype). 
	We will update more options to handle various file format. 

3. Positive and Negative
	The matrix should contain a column, called "TP", which holds a label of positive and negative for each variant.
	Our current method assumes that users create a positive and negative set of variants from pedigree sequencing data.
		- Positive set will be variants transmitted from parent to child. We assume these are likely true calls.
		- Negative set will be a Mendelian violation call in one child but also observed in one unrelated individual.

	For HQ rare variants, we refine our positive and negative set to AC(allele count) == 2 in your VCF.

3. Usage

	python rocube.py \
		-i input_file \
		-t number of threads (default: 1) \
		-o tag for output (default: "output") \
		-p yes if you need plotting (default: no) \
		-a minimum sensitivity for output (default: 0.90) \
		-b minimum specificity for output (default: 0.99)

'''

import os,sys,argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool 
from functools import partial
import pyximport; pyximport.install()
import rocube

def main(infile, number_threads, output_tag, min_sens, min_spec, plot):
	## Configuration
	## List of VCF quality metrics for this analysis
	metrics = ['va_QUAL',
				'va_AN',
				'va_BaseQRankSum',
				'va_DP',
				'va_FS',
				'va_MQ',
				'va_MQRankSum',
				'va_QD',
				'va_ReadPosRankSum',
				'va_SOR',
				'va_GQM',
				'g_GQ',
				'g_AB',
				'g_DP']

	## Load the data
	df = pd.read_csv(infile, compression='gzip', sep='\t', index_col=None)
	print '[Progress] Loaded an input data.'

	## Drop the parental lines from the data frame if exists
	df = df.drop(df[(df['Role'] == 'fa') | (df['Role'] == 'mo')].index)

	## Formatting some values
	## For example, AB (allelic balance), ExcessHet, FS and SOR
	## AB: An optimal AB value should be around 0.5. So a threshold will be calculated for the lowest (~0.0) or highest end (~1.0).Add a reverse direction to find out a threshold for high-depth bias
	## FS or SOR: These metrics have a reverse direction for high-quality
	df2 = df.g_AD.str.replace('[','').str.replace(']','').str.split(',', expand=True).astype(int)
	df['g_DP'] =  df2[0] + df2[1] ## Add genotype DP and AB if necessary
	df['g_AB2'] = df2[1] / df['g_DP'].astype(float)
	df['g_AB'] = df.apply(lambda row: 1 - row['g_AB2'] if row['g_AB2'] >= 0.5 else row['g_AB2'], axis=1) # We
	df = df.round({'va_DP': 0, 'va_GQM': 0, 'g_AB': 3})
	del df2 # Remove df2 
	df['va_ExcessHet'] = -df['va_ExcessHet'] 
	df['va_FS'] = -df['va_FS'] 
	df['va_SOR'] = -df['va_SOR'] 
	df = df.drop('g_AB2', 1) # Remove column
	# print df.head()

	## (Optional) Plot for the distribution of quality metrics
	if plot == 'yes':
		print '[Progress] Plot the distribution of quality metrics'
		import seaborn as sns
		import matplotlib.pyplot as plt
		from matplotlib.backends.backend_pdf import PdfPages
		# Creating a dataframe for plot after down-sampling an input
		plotMat = pd.melt(df.sample(n=len(df.index)/10), id_vars=['TP'], value_vars=metrics)
		g = sns.factorplot(x="TP", y="value", col="variable",
							data=plotMat, kind="violin", sharey=False, scale='width',
							col_wrap=4, size=4, aspect=1);
		g = g.set_axis_labels("", "Value").set_xticklabels(["TP", "FP"]).set_titles("{col_name}")
		outfile_plot = '.'.join(['plot.rocube.distQual', output_tag, 'pdf'])
		g.savefig(outfile_plot)
		print '[Progress] Plotting is done'
	else:
		print '[Progress] No plot option was given.'

	print '[Progress] DataFrame is ready for analysis.'

	## Count the number for raw positives and negatives
	set_pred = df.TP.values
	P_raw = np.count_nonzero(set_pred == 1)
	N_raw = np.count_nonzero(set_pred == 2)

	# Set the while loop
	i = 0
	sens_overal = 1
	spec_overal = 0
	roc_results = []
	while sens_overal > min_sens and spec_overal < min_spec:
		i = i + 1
		print '[Progress] Start ROC for the round %s' % (str(i))
		set_pred = df.TP.values
		P = np.count_nonzero(set_pred == 1)
		N = np.count_nonzero(set_pred == 2)
		print '[Progress] ROC Round %s training sets: Positives %s , Negatives %s' % (str(i), str(P), str(N))

		## ROC
		pool = mp.Pool(number_threads)
		# pool = ThreadPool(4) 
		# results = pool.map(partial(getMaxSpec, df=df, P=P, N=N), metrics)
		results = pool.map(partial(rocube.getMaxSpec, df=df, P=P, N=N), metrics)
		res_round = pd.DataFrame(results).rename(columns = {0:'metric', 1: 'cut', 2:'sens', 3:'spec'} )
		pool.close()
		pool.join() 

		## Print out the result from the current round
		print res_round.sort_values('spec', ascending=False)

		## Subset the dataframe for the next round
		metric_best_round = res_round['metric'][res_round['spec'].idxmax()]
		metrics.remove(metric_best_round) # update the list
		cut_best_round = res_round['cut'][res_round['spec'].idxmax()]

		if metric_best_round == 'g_AB':
			df = df[(df[metric_best_round] >= cut_best_round) & (df[metric_best_round] <= 1- cut_best_round)]
		else:
			df = df[df[metric_best_round] >= cut_best_round]

		## Calculate the overall sensitivity and specificify
		set_pred = df.TP.values
		TP = np.count_nonzero(set_pred == 1)
		FP = np.count_nonzero(set_pred == 2)
		# sensitivity or true positive rate (TPR); TP/P
		sens_overal = round(TP/float(P_raw), 5)
		# specificity (SPC) or true negative rate; (1-FP/N) TN/N
		spec_overal = round(1 - ( FP/float(N_raw) ), 5)

		if sens_overal > min_sens or spec_overal < min_spec:
			## Add to output
			out = [str(i), metric_best_round, cut_best_round, str(sens_overal), str(spec_overal), str(P), str(N)]
			roc_results.append(out)
			print '[Progress] Completed ROC for the round %s. The most probable metric from this round is %s' % (str(i), metric_best_round)
			print '[Progress] Round %s: sensitivity %s , specificity %s' % (str(i), str(sens_overal), str(spec_overal))
			print '------------------------------'
		else:
			## Skip for an output when sensitivity is lower than desired
			print '[Progress] Completed ROC for the round %s. The most probable metric from this round is %s' % (str(i), metric_best_round)
			print '[Progress] Round %s: sensitivity %s , specificity %s' % (str(i), str(sens_overal), str(spec_overal))
			print '[Progress] Quit the loop as the minimum sensitivity'
			print '------------------------------'

	## Write to an output file
	outfile = '.'.join(['table.rocube_result', output_tag, 'tsv'])
	o = pd.DataFrame(roc_results).rename(columns = {0: 'round', 1:'metric', 2: 'cut', 3:'sens', 4:'spec', 5:'train_P', 6:'train_N'} )
	o.to_csv(outfile, sep='\t', index=False, header=True)
	print '[Done] Completed ROC analysis'

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i','--infile', required=True, type=str, help='Input File')
	parser.add_argument('-t','--number_threads', required=False, type=int, help='Number of threads', default=1)
	parser.add_argument('-o','--output_tag', required=False, type=str, help='Output tag', default='output')
	parser.add_argument('-a','--min_sens', required=False, type=float, help='Minimum sensitivity', default=0.90)
	parser.add_argument('-b','--min_spec', required=False, type=float, help='Minimum specificity', default=0.99)
	parser.add_argument('-p','--plot', required=False, type=str, help='Do you need a plot? yes or no', default='no')
	args = parser.parse_args()
	main(args.infile, args.number_threads, args.output_tag, args.min_sens, args.min_spec, args.plot)





