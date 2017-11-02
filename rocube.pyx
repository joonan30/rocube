#!/usr/bin/env python

import pandas as pd
cimport cython
cimport numpy as np
import numpy as np


cpdef calRoc(np.ndarray[double] vals, np.ndarray[double] TPs, double P, double N, double cut):
	cdef double TP, FN, sens, spec

	TP = np.count_nonzero(TPs[vals > cut] == 1)
	FP = np.count_nonzero(TPs[vals > cut] == 2)

	# sensitivity 
	sens = round(TP/float(P), 5)
	# specificity (SPC) 
	spec = round(1 - ( FP/float(N) ), 5)

	return [cut, sens, spec]


def getMaxSpec(metric, df, P, N):
	## Get the list of score intervals
	qranges = sorted(df[metric].unique())
	if len(qranges) > 1000:
		## Select 1000 intervals if too many given
		qranges = qranges[0::len(qranges)/1000] 
	else:
		pass

	## Calculate ROC values
	a = pd.DataFrame(map(lambda cut: calRoc( df[metric].astype('float64').values, df['TP'].astype('float64').values, P, N, cut), qranges ), 
								dtype=float).rename(columns = {0:'cut', 1:'sens', 2:'spec'})

	## Hold sensitivity given the number of a negative set
	if N/float(P) > 0.01:
		a = a[a['sens'] >= 0.995]
	elif N/float(P) > 0.001:
		a = a[a['sens'] >= 0.99]
	else:
		a = a[a['sens'] >= 0.95]
	spec = a['spec'].max()
	sens = a['sens'][a['spec'].idxmax()]
	cut = a['cut'][a['spec'].idxmax()]

	print '[Progress] Achieved specificity %s and sensitivity %s from the quality metric %s (threshold: %s)' % (spec, sens, metric, str(cut))
	return [metric, cut, sens, spec]

