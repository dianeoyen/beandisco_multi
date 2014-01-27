BEANDiscoMulti - Bayesian Exact and Approximate Network Discovery for Multiple Related Networks

Copyright 2013 Diane Oyen <doyen(at)unm.cs.edu>

BEANDiscoMulti extends the BEANDisco software (http://www.cs.helsinki.fi/u/tzniinim/BEANDisco/)
for learning Bayesian network structures from multiple datasets with the same variables.
See [1] for more info.


Licence
=======

This software is distributed under GNU GPL. See COPYING.


Compilation
===========

Following libraries should be installed:

	gsl, boost, boost-program-options

Compilation:

	make


Usage
=====

Usage:

	./beandmulti [options] [infile [outfile]]

To see the list of available options:

	./beandmulti --help

The program will read an input data file for each task, named <infile>.taskX,
beginning with X=1. The input data file should contain one data sample per row, each sample
consisting of one integer value for each variable. Values on a row should be
separated by whitespace (tabs or spaces). For an example data file with 8
variables (columns) and 200 samples (rows) see example.task1 and example.task2.


Examples
========

Compute estimates with maximum in-degree of 3, (maximum) bucket size 4, burn-in
period of 1000 steps and 100 samples with 10 steps in between. First compute 
scores in to a file and then use the precomputed scores to estimate arc probabilities:

	./beandmulti example -t 2 -m 3 --score-file example.score
	beand --score-file example.score.task1 -b 4 -B 1000 -s 100 -S 10 
	beand --score-file example.score.task2 -b 4 -B 1000 -s 100 -S 10

Compute exact probabilities (generally with exact computation it is recommended
to use a bucket size as large as possible with the available memory):

	./beandmulti example -t 2 -m 3 --score-file example.score
	beand --score-file example.score.task1 -b 4 --exact
	beand --score-file example.score.task1 -b 4 --exact


Parameters
==========

Some recommendations for important parameter values:

maximum indegree (-m)
	If the number of variables is high, this is restricted by memory and time
	consumption. For over 100 variables 3 is probably reasonable. On the other
	hand, for about 30 variables this might be increased to 5.

order type (--order-type)
	Use bucket order (po), which is the default.

bucket size (-b)
	Good value depends on the number of variables and maximum indegree. Setting
	this to 10 is probably a reasonable choice. In general higher values are
	better, so it is a good idea to test different values and choose the largest
	one which still do not increase the time consumption per step too much.

number of samples (-s)
	The higher the better.

number of steps per sample (-S)
	Something like 5 or 10 is reasonable.

number of burn-in steps (-B)
	To ensure good convergence before the actual sampling starts, I would set
	this about equal to the number of total steps in sampling stage (B = s * S).
	For example: s = 2000, S = 10, B = 20000


References
==========
[1] D. Oyen and T. Lane. Bayesian Discovery of Multiple Bayesian Networks via 
Transfer Learning. ICDM 2013

[2] T. Niinimäki, P. Parviainen and M. Koivisto. Partial Order MCMC for
Structure Discovery in Bayesian Networks. UAI 2011

[3] P. Parviainen and M. Koivisto. Bayesian Structure Discovery in Bayesian
Networks with Less Space. AISTATS 2010


