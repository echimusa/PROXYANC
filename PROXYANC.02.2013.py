#!/usr/bin/env python
# -*- coding: latin-1 -*-
"""
%prog [options] files 

run PROXYANC of the input file(s).
eg. python PROXYANC.0.1.py paranc.txt
"""
from __future__ import with_statement
from __future__ import division

# This program is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
'''
1. PROXYANC is a program for choosing the best proxy ancestry for multi-way admixed populations given pools (groups) of reference populations.
2. PROXYANC makes use of two novel algorithms including the correlation between observed Linkage Disequilibrium in an admixed population and population genetic differentiation in ancestral populations, and an optimal quadratic programming based on the linear combination of population genetic distances. Mode details of the Methods can be found in chimusa et al. (2013) paper. 

3. PROXYANC can also select panel of AIMs based on can select AIMs based on the relationship between the observed local multi-locus linkage disequilibrium in a recently admixed population and ancestral population difference in allele frequency and based on the Kernel principal component analysis (Kernel-PCA), which is the extension of the linear PCA.

4. PROXYANC can identify possible unusual difference in allele frequency between pair-wise popualtions, as signal of natural selection.


'''

import numpy as np
from proxyOpt import *
import os,sys,exceptions, types, time,datetime
import fileinput,string,pickle
from scipy import *
from types import *
from random import*
from scipy.stats import *
from scipy import stats
from scipy.stats import norm
from scipy.stats import chi2
from scipy.stats import beta
from scipy.optimize import *
from rpy import*
from itertools import *
from operator import itemgetter
from matplotlib import pyplot as plt
from operator import itemgetter
from cvxopt import matrix, solvers, spmatrix, mul, div
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

####################################################################################
# Writing the result from admixture LD and Expected admixture LD
####################################################################################

def write_expLD(pop,out_d):
	import pylab
	for p in pop:
		rows_ld ={}
		pop = []
		admix = []
		pos = []
		ld_pop = []
		pop_anc =[]
		rows_ld1 ={}
        	for line in fileinput.input(p):
                	data =line.split()
			step = fileinput.lineno()
                        if fileinput.lineno()>1:
                       	        data = line.split()
                       	        rows_ld[fileinput.lineno()]=[int(data[0])]+data[1:]
                       	else:
				lows = data[:-1]+[p.split(".")[0]]
		chose = random_snp(1000,1.2,rows_ld)
		for d in chose:
			rows_ld1[d] = rows_ld[d]
		X = rows_ld1.values()
		X.sort()
		for d in X:
			pos.append(d[0])
			admix.append(float(d[1]))
			ld_pop.append(float(d[1]))
			pop_anc.append(float(d[3]))
		rows_ld.clear()
		rows_ld1.clear()
 		pylab.figure(1)
		pylab.xlabel("Relative Physical Distance")
		pylab.ylabel("Expected LD")
		#pylab.legend(tuple(["Admixed pop",lab[idx]]),loc='lower right')
		pylab.plot(pos,admix,'r',linewidth=5)
		pylab.plot(pos,ld_pop,'k',linewidth=5)
		pylab.plot(pos,pop_anc,'m',linewidth=5)
		pylab.grid(True)
		pylab.savefig(out_d+p.split("/")[-1].split(".")[0]+".png")
		plot_R_exp(admix,ld_pop,p.split("/")[-1].split(".")[0])
		plot_R_exp(admix,pop_anc,p.split("/")[-1].split(".")[0])
		
####################################################################################
# Plotting function
####################################################################################
def plot_R_exp(admix,axiz,pop):
	from scipy import stats
	from rpy import r
	ls_fit = r.lsfit(admix,axiz)
	gradient = ls_fit['coefficients']['X']
	yintercept = ls_fit['coefficients']['Intercept']
	nam = "scatter_"+pop+".png"
	r.png(nam, width=400, height=700)
	r.plot(x=admix, y=axiz, xlab="LD in SAC", ylab="Expected LD", xlim=(0,7), ylim=(-16,27),main=" ")
	r.abline(a=yintercept, b=gradient, col="red")
	r.dev_off()

def plot_LDexp(Y,out_dir):
	"""plot_LDexp(x,f) function writes the results of relationship between LD in admixed population and Expected admixture LD from ancestral populations.
	 x is a dictionary of LD in admixed population and expected admixture LD from ancestral population and relative distance betwee SNPs."""
	stat_LD = {}
	pop =[]
	
	ancestral = Y.keys()	
	Axis_ancestral = []
	for i in range(len(ancestral)):
		Axis_ancestral.append([])
	for anc in ancestral:
		idx = ancestral.index(anc)
		tmp = []
		tmp1 = []
		pos = []
		admix = []
		labels = []
		for des in Y[anc]:
			tmp.append(Y[anc][des][0])
			tmp1.append(Y[anc][des][1])
			pos.append(Y[anc][des][2])
			admix.append(Y[anc][des][-1]) 
		Axis_ancestral[idx].append(tmp)	
		Axis_ancestral[idx].append(tmp1)
		Axis_ancestral[idx].append(pos)
		Axis_ancestral[idx].append(admix)
	
	for axiz in Axis_ancestral:
		idx = Axis_ancestral.index(axiz)
		fin=open(out_dir+ancestral[idx]+".exp","wt")
		fin.writelines("POS"+"\t"+"AdmixLD"+"\t"+"ExpLD"+"\t"+"ancProrp"+"\n")
		for i in range(len(axiz[0])):
			fin.writelines(str(axiz[2][i])+"\t"+str(axiz[-1][i])+"\t"+str(axiz[0][i])+"\t"+str(axiz[1][i])+"\n")
		fin.close()
		pop.append(out_dir+ancestral[idx]+".exp")
		labels.append(ancestral[idx])
		value = stats.linregress(axiz[-1],axiz[0])
		stat_LD[ancestral[idx]] = [str(value[3]),str(value[2]),str(value[2]**2)]
	fi = open(out_dir+"Stat_LD.out","wt")
	fi.writelines("POP"+"\t"+"Pvalue"+"\t"+"R-sqaure"+"\t"+"R2"+"\n")
	for des in stat_LD:
		fi.write(des+"\t"+"\t".join(stat_LD[des])+"\n")
	fi.close()
	#write_expLD(pop,out_dir)



####################################################################################
# Plot the result from Unusual popualtion differeiation.
####################################################################################		
def _gen_data(fhs, columns, sep):
    """
    iterate over the files and yield chr, start, pvalue. 
    """
    for line in fhs:
        if line[0] == "#" or line[0] in ['C','H','R']: continue
	if len(line) != 0:
	    	toks = line.split() 
            	yield toks[columns[0]], int(toks[columns[1]]), float(toks[columns[2]])

def chr_cmp(a, b):
    a = a.lower().replace("_", "")
    b = b.lower().replace("_", "")
    achr = a[3:] if a.startswith("chr") else a
    bchr = b[3:] if b.startswith("chr") else b

    try:
        return cmp(int(achr), int(bchr))
    except ValueError:
        if achr.isdigit() and not bchr.isdigit(): return -1
        if bchr.isdigit() and not achr.isdigit(): return 1
        return cmp(achr, bchr)
def chr_loc_cmp(alocs, blocs):
    return chr_cmp(alocs[0], blocs[0]) or cmp(alocs[1], blocs[1])
def manhattan(fhs, columns, image_path, no_log, colors, sep, title, lines, ymax):
    """Generate a manhattan plot: inputs: file, column of pvalue or of meseasure of interest, path to out put the plot, none, a string of color see parameter file,none,none,none. """
    xs = []
    ys = []
    cs = []
    colors = cycle(colors)
    xs_by_chr = {}
    last_x = 0
    data = sorted(_gen_data(fhs, columns, sep), cmp=chr_loc_cmp)
    for seqid, rlist in groupby(data, key=itemgetter(0)):
        color = colors.next()
        rlist = list(rlist)
        region_xs = [last_x + r[1] for r in rlist]
        xs.extend(region_xs)
        ys.extend([r[2] for r in rlist])
        cs.extend([color] * len(rlist))
        xs_by_chr[seqid] = (region_xs[0] + region_xs[-1]) / 2
        # keep track so that chrs don't overlap.
        last_x = xs[-1]
    xs_by_chr = [(k, xs_by_chr[k]) for k in sorted(xs_by_chr.keys(), cmp=chr_cmp)]
    xs = np.array(xs)
    ys = np.array(ys) if no_log else -np.log10(ys)
    plt.close()
    f = plt.figure()
    ax = f.add_axes((0.1, 0.09, 0.88, 0.85))
    if title is not None:
        plt.title(title)
    ax.set_ylabel('-log10(p-value)')
    if lines:
        ax.vlines(xs, 0, ys, colors=cs, alpha=0.5)
    else:
        ax.scatter(xs, ys, s=3, c=cs, alpha=0.8, edgecolors='none')
    x = 7.2*(0.000000001)
    ax.axhline(y=-np.log10(x), color='0.5', linewidth=2)
    plt.axis('tight')
    plt.xlim(0, xs[-1])
    plt.ylim(ymin=0)
    if ymax is not None: plt.ylim(ymax=ymax)
    plt.xticks([c[1] for c in xs_by_chr], [c[0] for c in xs_by_chr], rotation=-90, size=8.5)
    print >>sys.stderr, "saving to: %s" % image_path
    plt.savefig(image_path)
   

####################################################################################
# Correlation and CI
#################################################################################### 
def get_r(x,y):
	"""get_r : inputs two vectors of the same size and return the correlation between these two vectors of values. """	
	if len(x)!=len(y):
		logger.info('Error: lengths of vectors do not match in get_corr')
	n = len(x)
	xSum = 0.0
	xxSum = 0.0
	ySum = 0.0
	yySum = 0.0
	xySum = 0.0
	for j in range(n):
		xVal = x[j]
		yVal = y[j]
		xSum +=  xVal
		xxSum += xVal * xVal
		ySum += yVal
		yySum += yVal * yVal
		xySum += xVal * yVal
		cov = xySum - (xSum * ySum) / float(n)
		xVar = xxSum - (xSum * xSum) / float(n)
		yVar = yySum - (ySum * ySum) / float(n)
		den = sqrt(xVar*yVar)
	if den>0:
		return cov/float(den)
	else:
		return 0

#*******************************This Section contains function to process bootstrap **********************************
####################################################################################
# model function and create gaussian noisy data from model
####################################################################################
def mean_confidence_interval(data, confidence=0.95):
	""" mean_confidence_interval(data, confidence=0.95), inputs: data vectors and return a interval of confidence."""
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), np.std(a)
	h = se * stats.t._ppf((1+confidence)/2., n-1)
	return m, m-h, m+h,se

####################################################################################
# Given a length, this function generates an identity matrice
####################################################################################
def identity(n):
	from cvxopt import matrix
	"""
	indentity(n), input: a size constant and return the identity matrix of size n
	"""
	I = matrix(0.0, (n, n))
	I[::n+1] = 1.0
	return I
####################################################################################
#Solving Quadratic programm.
####################################################################################

def object_mat(list_all,OPT):	
	""" cvxopt.solvers.coneqp(P, q[, G, h[, dims[, A, b[, initvals[, kktsolver]]]]]) Solves a pair of primal and dual quadratic cone programs
	inputs: list of popualtion  minor allele frequency vectors and popualtion size. Computer the Fst objective functions describes in Chimusa et al. and return possible solutions."""
	try:
		scale = float(OPT[0]*list_all[1]*(1-list_all[1]))
		pcale = float(list_all[1]*(1-list_all[1]))
		h_dem = float(list_all[0]*(list_all[0]-1))
		h_num = float(list_all[2]*(list_all[0]-list_all[2]))
		h_bar = h_num/h_dem
		h_scale = OPT[0]*h_bar
		tmp_lin = []
		Mat = array([[0.0 for row in range(OPT[0])] for col in range(OPT[0])])
		list_beta = list_all[3:]
		step = int(len(list_beta)/2)
		for i in range(step):
			if OPT[1]:
				z = (list_beta[i+1]**2)/float(h_scale) 
				n = -2*list_beta[i+1]*list_all[1] 
				m = list_all[1]**2 + n
				e = m/float(h_scale)
				Mat[i,i] = z 
				tmp_lin.append(e)
				for j in range(i+1,step):
					w = ((2*list_beta[i]*list_beta[j])/h_scale) 
					Mat[i,j] = w/2
					Mat[j,i] = w/2
			else:
				x = list_all[1]**2
				m = h_bar/float(list_all[0])
				n = -2*list_beta[i+1]*list_all[1]
				e = (x - m + n)/float(h_scale)
				if  list_beta[i] ==0.0:
					y = h_bar #pcale
				else:
					y = (h_bar/list_beta[i])
				z = ((list_beta[i+1]**2)-y)/float(h_scale)
				Mat[i,i] = z 
				tmp_lin.append(e)
				for j in range(i+1,step):
					w = ((2*list_beta[i]*list_beta[j])/h_scale) 
					Mat[i,j] = w/2
					Mat[j,i] = w/2
			w = 0
			z = 0
			e = 0
	except:
		raise exceptions.SystemError('Failed to compute the objective Matrice')
	
	min_fst,sol_opt = coneqp_solver(Mat,tmp_lin)
	fst_min = optimal_fst(list_all,scale,pcale,sol_opt)
	list_all = []
	return min_fst,fst_min
		
####################################################################################
### QP solver (FST optimal mixed strategy)
####################################################################################

def coneqp_solver(objective,lin_vect):
    """
    Input matrix M is m x m from the FST of tagged population and linear combination of reference populations.
    Return value x that is an optimal mixed strategy that minimizes
    sum of squares of x_i. Uses function qp from cvxopt library
    """
    P = matrix(objective)
    m, n = P.size
    # convert bjective matrix to cvxopt matrix object M and negate

    # make M all positive by adding large constant v
    I = identity(n)

    # set up P, q so that minimizing sum of squares of p_i is
    # equivalent to minimizing 1/2 x^T P x + q^T x
                        # P is m x m
    q = matrix(lin_vect)                    # q is m x 1
    G = matrix([-I, matrix(0.0, (1,n)), I])
    h = matrix(n*[0.0] + [1.0] + n*[0.0])        # h is 2m x 1
    dims = {'l': n, 'q': [n+1], 's': []}
    # According to the CVXOPT documentation the following requirement on G and A should also be met, 
    # (1)  rank(A) = p                  (where p = # rows of A)
    # (2)  rank(matrix([P,G,A]) = n     (where n = # columns in G and in A)
    # (this last has P stacked on top of G on top of A)
   
    # solve constrained least squares problem
    solvers.options['feastol']=1e-6 #1e-6      # slightly relaxed from default (avoids singular KKT messages)
    solvers.options['abstol']= 1e-9 #1e-9      # gives us good accuracy on final result
    solvers.options['show_progress']=True
    sol = solvers.coneqp(P.T*P, q, G, h, dims)
    x = sol['x']
    y = sol['primal objective']
    # if any were even slightly negative, round up to zero
    for i in range(m):
        x[i] = max(0.0,x[i])

    # return optimal mixed strategy that minimizes sum of squares
    # sum of x[i]'s should be 1.0/v.  Normalizing gives probability distribution.
    # This should be equivalent to, but more reliable than, simply multiplying by.        
    sumx = sum(x)
    x = [ xi / sumx for xi in x]

    return y,x
####################################################################################	
#Genome Fst from normal FST function and Quadratic program solution.
####################################################################################
def genome_fst(tagg_all,OPT):
	"""  Computing the genetic distance between two populations"""
	# tagg_all: is a list of sample size of population and allele frequencies
	fst_genome = []
	fst_genome1 = []
	if len(tagg_all[1]) == len(tagg_all[3]):
		for i in range(len(tagg_all[1])):
			tagg_genome = []
			for j in range(len(tagg_all)):
				tagg_genome.append(tagg_all[j][i])
			fst1,fst2 = object_mat(tagg_genome,OPT)
			fst_genome.append(fst1)
			fst_genome1.append(fst2)
		tagg_all = []
		G = []
		F = []
		for j in fst_genome1:
			if j == nan or j ==-inf:
				pass
			else:
				G.append(j)
		for j in fst_genome:
			if j == nan or j ==-inf: 
				pass
			else:
				F.append(j)
		return average(F),min(G),F,G
	else:
		sys.stdout.write("Inconsistency in data")

####################################################################################
#This function use Rpy to create a template for linear regression
####################################################################################
def make_R_strings(dv,iv):
# This creates two strings that RPy uses to execute the regression command
# First, it creates the actual R expression of the linear model, which is
# always of the form 'DependantVar ~ IndependentVar1 + IndependentVar2 + ... +IndependentVarK'.
# Then, Python must reference back to data stored in memory, therefore, the second 
# string is the data frame based on the user input.
    dv_label=dv.keys()  # Store variable labels
    iv_labels=iv.keys()
    
    # The dependent variable always comes first
    R_string=dv_label[0]+" ~ "
    frame=dv_label[0]+"=dv['"+dv_label[0]+"'], "
    
    last_var=iv_labels[-1]  # note the final iv
    
    # Then add independent variables as appropriate
    for v in iv_labels:
        if v==last_var:
            R_string=R_string+v
            frame=frame+v+"=iv['"+v+"']"
        else:
            R_string=R_string+v+" + "
            frame=frame+v+"=iv['"+v+"'], "

    return R_string,frame
####################################################################################
#This function use Rpy for linear regression
####################################################################################
def regress(dv,iv):
# Performs regression using R's linear model function (lm)
    if type(dv.values()) is list and all(type(x) is list for x in iv.values()):
    # First check that all of the data is in list form, otherwise RPy will throw an error
        set_default_mode(NO_CONVERSION) # Keeps values in R format until we need them
	R_string,frame=make_R_strings(dv,iv)    # Create strings used by RPy to run regression
        
        # R runs the linear regression
        OLS_model=eval('r.lm(R_string, data=r.data_frame('+frame+'))')
        set_default_mode(BASIC_CONVERSION)  # Now convert back to usable format
        
        model_summary=r.summary(OLS_model)      # Store resultss
	
        coeff=model_summary['coefficients'][:,0]    # Regression coeffecients
        std_err=model_summary['coefficients'][:,1]  # Standard Errors
        t_stat=model_summary['coefficients'][:,2]   # t-statistics
        p_val=model_summary['coefficients'][:,3]    # p-values
        r_sqr=model_summary['r.squared']            # R-squred
        asj_r_sqr=model_summary['adj.r.squared']    # Adjusted R-squared
	pvalue = p_val[1]	
	return pvalue
    else:
        raise TypeError("All variables must all be of type 'list'")

####################################################################################
# This function update and complete the information of the association between the ancestry  
#allele frequencies  the and LD in tagged population in all possible lists.
####################################################################################
def completeness(list1,pop):
	list2 = []
	k = len(pop[0])
	for i in range(k,len(list1)-1):
		leng = abs(len(list1[i-k])-len(list1[i]))
		for j in range(leng):
			list1[i].insert(j,list1[j][i-k])
		#if len(list1[i]) != len(list1[0]):
		#	raise TypeError("must have equal size!")
	return list1

####################################################################################
#This normalize le ancestry proxy score
####################################################################################    
def normalization(list1,list2,sample):
	"""Inputs: score stats, physical position and return score of the deviation from the mean of each score stats """
	lenght= len(list1)
	list3 = []
	for i in range(len(list1)):
		list3.append([])
	tmp = []
	int_conf = []
	for des in list1:
		idx = list1.index(des)
		Z =sum(des)/sqrt(len(des))
		tmp.append(Z)  
		m, C, D,se = mean_confidence_interval(des, confidence=0.95)
		c = "("+"%.5f" % round(C,6)+ ","+ "%.5f" %round(D,6)+")"
		list3[idx].append(list2[idx])
		a1 = mean(tmp)
		a2 = var(tmp)
		a3 = Z-a1
		b = a3/float(a2)
		list3[idx].append(b)
		list3[idx].append([c,str(se/float(a2))])
	list1=[]
	list2 = []
	return list3,tmp
####################################################################################  	
# This function takes x,y as dictionaries of both LD from admixed pop and difference 
# in allele frequencies in two populations and returns an associated pvalue.
####################################################################################
def get_weight_pvalue(x,y):
		tmp = []
		tmp1 = []
		try:
			ld1 ={}
			ld2 = {}
			for des in x:
				rsid = des.split("-")
				sam = y[rsid[0]] * y[rsid[1]]
				tmp.append(sam) 
				tmp1.append(float(x[des][0]))
			if len(tmp) == len(tmp1):
				ld1["LD"]=tmp
				ld2["FRQ"]=tmp1	
				#gradient, intercept, r_value, p_value, std_err = stats.linregress(tmp,tmp1)
				p_value = regress(ld1,ld2)
			
			return float(p_value)
		except:
			raise exceptions.SystemError('Failed to compute the pvalues')	

####################################################################################
# Get scaled ancstral allele frequency
####################################################################################
def ancestral_allele(x,y):
	scale = []
	if len(x) == len(y):
		for i in range(len(x)):
			scale.append(average([x[i],y[i]])*(1-average([x[i],y[i]])))
		return scale
####################################################################################
# Convert a dictionary to a list of keys and data
####################################################################################
def convert_dict_list(dict1):
	list1 = []
	list2 = []
	for d in dict1:
		list1.append(dict1[d])
		list2.append(d)
	return [list1,list2]
####################################################################################
# Convert two lists into a dictionary
####################################################################################	
def convert_list_dict(list1,list2):
	dict1 = {}
	if len(list1) == len(list2):
		for i in range(len(list1)):
			dict1[list1[i]]=list2[i]
	return dict1	
####################################################################################
# Compute difference in allele frequency.
# x, y are vectors of allele frequencies from two Reference populations over SNPs.
####################################################################################
def get_freq_diff(x,y):
	diff = []
	if len(x) != len(y):
		sys.stdout.write("\nLength of the vectors of allele frequencies should be of same (Skipping ...)")
	else:
		try:
			tmp = []	
			for i in range(len(x)):
				diff.append(abs(float(x[i]) - float(y[i])))
		except:
			raise exceptions.SystemError('Failed to compute diffence in allele frequency on SNP')
		return diff
####################################################################################
#Compute allele frequency
# sample is a string of allele values 0/1/2 number of reference alleles.
####################################################################################
def get_freq(sample):
	ref_al_count=sample.count("2") + sample.count("1")
	minor_al_count=sample.count("0")
	if "?" in sample:
		missing=sample.count("?")
		total_count = ref_al_count + minor_al_count + missing
	elif "9" in sample:
		missing = sample.count("9")
		total_count = ref_al_count + minor_al_count + missing
	else:
		total_count = ref_al_count + minor_al_count
	total_allele = ref_al_count + minor_al_count
    	if total_count != len(sample):
		logger.info('Genotype file does not match observed alleles this marker (exiting)')
		raise SystemExit
	else:
		ref_al_frq = ref_al_count/float(total_count)
		minor_al_frq = minor_al_count/float(total_count)
		return ref_al_frq,total_allele,minor_al_frq,minor_al_count

# Ancestral populations allele frequency.
def allele_frq_process(pop_dict,op):
	""" Compute ancestral allele frequencues from their genotype data and return allele frq."""
	# pop_dict: genetype data
	C1=[]
	C2=[]
	D = {}
	for rsid in pop_dict:
		"""loop over all SNP and provide genotype as a single string of 0,1,,2 sequences. """
		frq,count_all,minor,minor_count = get_freq(pop_dict[rsid])
		if op[0]:
			C1.append(minor)
			C2.append(count_all) 
		else:
			D[rsid] = minor
	return D, C1, C2

####################################################################################
#This function computes the correlation between two genotype data.
#genotypes1,genotypes2 are python lists
####################################################################################        
def calc_rsq(genotypes1,genotypes2,rsid1,rsid2):
	gen1,gen2 = remove_missing(genotypes1,genotypes2,rsid1,rsid2) 
	"""calc_rsq(genotypes1,genotypes2,rsid1,rsid2), inputs: two genotypes with related SNPs and return the person correlation."""
	snp=[]
	if len(gen1)==0 and len(gen2)==0:
		snp.append(rsid2)
		snp.append(0.0)
      		return snp # too much missing data to compute r
	elif len(gen1) != 0 and len(gen2) != 0:
   		
		if len(gen1) ==len(gen2):
			corr = get_r(gen1,gen2)
			snp.append([rsid1,rsid2])
			snp.append(corr)
   			return snp
		else:
			sys.stdout.write("\ngenotypes do not match (exiting)")
			raise SystemExit

####################################################################################
# This function test the optimal solution,
# from Qradratic program into the objective function of interest.
####################################################################################
def optimal_fst(frq_list,dem,semi_dem,varb) :
	p_t = semi_dem/frq_list[0]
	p_1 = frq_list[1]
	p_k = 0
	n_k = 0
	step = 0
	for i in range(int(len(frq_list[2:])/2)):
		if i%2 != 0:
			be = varb[i-step]
			p_k = p_k + frq_list[2:][i]*be
		else:
			be = varb[i-step]
			n_k = n_k + ((be**2)*semi_dem)/float(frq_list[2:][i])
			step = step + 1
	return (((p_k - p_1)**2) - p_t - n_k)/dem

####################################################################################	 
# Computing pairwise population genetic distance.
###################################################################################
def pwFST(tagg_pop,fst_frq_pop,n_size_pop,poplist,snp,out_dir,OPT):
	""" pwFST(tagg_pop,fst_frq_pop,n_size_pop,poplist,snp,out_dir,OPT), inputs: allele freq for pop1 and popK (k=1...K), sixe of populations genotypes, list of populations, SNps files, and option to compute Fst (unrelated /related individuals.) """
	n_size_pop.insert(0,tagg_pop[0])
	fst_frq_pop.insert(0,tagg_pop[1])
	tagg_pop
	expond =False
	if OPT[-1]:
		rows = {}
	else:
		rows ={}
		STD = {}
		for des in poplist:
			rows[des] = []
			STD[des] = []
	for subset in combinations(fst_frq_pop,OPT[0]):
		pops = []
		tag_frq = []
		n_size = []
		diff_frq = []
		for des in list(subset):
			idx = fst_frq_pop.index(des)
			if OPT[-1]:
				pops.append(poplist[idx])
				n_size.append(n_size_pop[idx])
				tag_frq.append(n_size_pop[idx])
				tag_frq.append(des)
				diff_frq.append(des)
			else:
				pops.append(poplist[idx])
				tag_frq.append(n_size_pop[idx])
				tag_frq.append(des)
			
		fst_genome = []
		if len(tag_frq[1]) == len(tag_frq[3]):
			for i in range(len(tag_frq[1])):
				tagg_genome = []
				for j in range(len(tag_frq)):
					tagg_genome.append(tag_frq[j][i])
				fst = pop_pairwisefst(tagg_genome,OPT[1])
				fst_genome.append(fst)
			if OPT[-1]:
				sam = abs(mean(fst))
				logger.info('Fst between:%s:%s' %("-".join(pops),str(sam)))
				diff = get_freq_diff(diff_frq[0],diff_frq[1])
				p = ancestral_allele(diff_frq[0],diff_frq[1])
				marker = snp.keys()
				diff_dict = convert_list_dict(marker,diff)
				rows["-".join(pops)]=n_size+[p]+[diff_dict]+[sam]
			else:
				f = round(mean(fst),5)
				sam = "%.3f" %f 
				if f == 0.0:
					expond = True
				logger.info('Fst between:%s:%s' %("-".join(pops),str(sam)))
				sam1 = str("%.3f" %(float(std(fst_genome)*10000)))
				
				rows[pops[0]].append(sam)
				STD[pops[0]].append(sam1)
	if OPT[-1] :
		return rows
	else:
		fin=open(out_dir+"PairwiseFST.out","wt")
		fin.write("****** Population Genetic Distance : Fst *****"+"\n")
		fin.write("\n")
		fin.write(" "+"\t"+"\t".join(poplist)+"\n")
		for pop in poplist[:-1]:
			ix = poplist.index(pop)		
			step = ix+1
			sam = ["-"]*step
			tp =sam+rows[pop]
			fin.write(pop+"\t"+"\t".join(tp)+"\n")
		fin.write("\n")
		fin.write("********* Standard Error: Fst*10000 *********"+"\n")
		fin.write("\n")
		fin.write(" "+"\t"+"\t".join(poplist)+"\n")
		for pop in poplist[:-1]:
			ix = poplist.index(pop)		
			step = ix+1
			sam = ["-"]*step
			tp =sam+STD[pop]
			if ix == 0:
				fin.write(pop+" "+" ".join(tp)+"\n")
			else:
				fin.write(pop+"\t"+"\t".join(tp)+"\n")
		fin.write("\n")
		fin.write("PROXYANC* Finished at:"+"\t"+str(datetime.datetime.today()) +"\n")
		fin.close
#Product Cathersian of diffirent group of populations to avoid repetition between populations closely Ethnic or group:
def product_cathersian(ad_pop,cath_list,OP):
	anc = []
	tag_frq = []
	tag_frq.append(ad_pop[0])
	tag_frq.append(ad_pop[1])
	tag_frq.append(ad_pop[2])		
	for des in list(cath_list):
		anc.append(des[0][0])
		tag_frq.append(des[0][2])#n_size_pop[idx])
		tag_frq.append(des[0][1])
	print"\nProcessing the linear combination of %s"%("-".join(anc))
	fst,fst1,gen,gen1= genome_fst(tag_frq,OP)
	c = "("+str("%.6f" %float(gen[0]))+ ","+str("%.6f" %float(gen[1]))+")"
	c1 = "("+str("%.6f" %float(gen1[0]))+","+str("%.6f" %float(gen1[1]))+")"
	SE = [str(abs(gen[2])),str(abs(gen1[2]))]
	return anc,abs(fst),abs(fst1),c,c1,SE

####################################################################################	 
# Scoring best ancstry proxy function based on FST througn solver optimal quadratic program.
###################################################################################
def FST_Score_anc(tagg_pop,fst_frq_pop,out_dir,OPT):

	if OPT[3]:
		pass
	else:
		fin=open(out_dir+"OptFST.out","wt")
		fin.writelines("POPs"+"\t"+"FST"+"\t"+"SE"+"\t"+"95% CI"+"\t"+"optFST"+"\t"+"95% CI"+"\t"+"SE"+"\n")
		step = 1
		if int(OPT[0]) == 2:
			for subset in product(fst_frq_pop[0],fst_frq_pop[1]):
				pops,fst_pop,fst_pop1,confident,confident1,sol2 = product_cathersian(tagg_pop,subset,OPT)
				fin.writelines("-".join(pops)+"\t"+str("%.6f" %fst_pop)+"\t"+sol2[0]+"\t"+confident+"\t"+str("%.6f" %fst_pop1)+"\t"+confident1+"\t"+sol2[1]+"\n")#"-".join(str(i) for i in sol2)+"\n")
			fin.close()
		elif int(OPT[0]) == 3:
			for subset in product(fst_frq_pop[0],fst_frq_pop[1],fst_frq_pop[2]):
				pops,fst_pop,fst_pop1,confident,confident1,sol2 = product_cathersian(tagg_pop,subset,OPT)
				fin.writelines("-".join(pops)+"\t"+str("%.6f" %fst_pop)+"\t"+sol2[0]+"\t"+confident+"\t"+str("%.6f" %fst_pop1)+"\t"+confident1+"\t"+sol2[1]+"\n")#"-".join(str(i) for i in sol2)+"\n")
			fin.close()
		elif int(OPT[0]) == 4:
			for subset in product(fst_frq_pop[0],fst_frq_pop[1],fst_frq_pop[2],fst_frq_pop[3]):
				pops,fst_pop,fst_pop1,confident,confident1,sol2 = product_cathersian(tagg_pop,subset,OPT)
				fin.writelines("-".join(pops)+"\t"+str("%.6f" %fst_pop)+"\t"+sol2[0]+"\t"+confident+"\t"+str("%.6f" %fst_pop1)+"\t"+confident1+"\t"+sol2[1]+"\n")#"-".join(str(i) for i in sol2)+"\n")
			fin.close()
		elif int(OPT[0]) == 5:
			for subset in product(fst_frq_pop[0],fst_frq_pop[1],fst_frq_pop[2],fst_frq_pop[3],fst_frq_pop[4]):
				pops,fst_pop,fst_pop1,confident,confident1,sol2 = product_cathersian(tagg_pop,subset,OPT)
				fin.writelines("-".join(pops)+"\t"+str("%.6f" %fst_pop)+"\t"+sol2[0]+"\t"+confident+"\t"+str("%.6f" %fst_pop1)+"\t"+confident1+"\t"+sol2[1]+"\n")#"-".join(str(i) for i in sol2)+"\n")
			fin.close()
		elif int(OPT[0]) == 6:
			for subset in product(fst_frq_pop[0],fst_frq_pop[1],fst_frq_pop[2],fst_frq_pop[3],fst_frq_pop[4],fst_frq_pop[5]):
				pops,fst_pop,fst_pop1,confident,confident1,sol2 = product_cathersian(tagg_pop,subset,OPT)
				fin.writelines("-".join(pops)+"\t"+str("%.6f" %fst_pop)+"\t"+sol2[0]+"\t"+confident+"\t"+str("%.6f" %fst_pop1)+"\t"+confident1+"\t"+sol2[1]+"\n")#+"-".join(str(i) for i in sol2)+"\n")
			fin.close()
		elif int(OPT[0]) == 7:
			for subset in product(fst_frq_pop[0],fst_frq_pop[1],fst_frq_pop[2],fst_frq_pop[3],fst_frq_pop[4],fst_frq_pop[5],fst_frq_pop[6]):
				pops,fst_pop,fst_pop1,confident,confident1,sol2 = product_cathersian(tagg_pop,subset,OPT)
				fin.writelines("-".join(pops)+"\t"+str("%.6f" %fst_pop)+"\t"+sol2[0]+"\t"+confident+"\t"+str("%.6f" %fst_pop1)+"\t"+confident1+"\t"+sol2[1]+"\n")#"-".join(str(i) for i in sol2)+"\n")
			fin.close()
		else:
			logger.info('PROXYANC cannot unhandle more than 8 proxy ancestry!')
	

#This function checks the population group in order to not estimate difference in allele frq between population of closely group :
def check_pop_group(pop_list,p1,p2):
	score = False
	for da in pop_list:
		if p1 in da and p2 in da :
			score = True
	return score

####################################################################################	 
# Scoring the best ancestry proxy function based on association between admixture 
#LD and difference in ancestral allele frequency.
####################################################################################
def FRQ_score_anc(group_frq,ld_g,pop_label,out_dir,opti):
	"""FRQ_score_anc(),inputs:alleles frequency of Group of populations, Ld in admixed population, list of popualtion labels, output path and options(see paramter files)."""
	pop_diff_frq = []
	poplist = []
	for des in group_frq:
		tmp = []
		for de in des:
			tmp += de
		pop_diff_frq += tmp

	for des in pop_label:
		poplist +=des
	if opti[3]:
		aims_list = [] # list of difference in allele frq.
	else:
		pvalue_list = [] # list of scaling pvalues based -2*log
		for k in range(len(poplist)):
			pvalue_list.append([])
	#Calculate pair-wise Difference in allele frequencies from pair of reference populations
	pop_done = []
	for l in range(len(pop_diff_frq)):
		list1 = convert_dict_list(pop_diff_frq[l])
		pvalue = []
		pops=poplist[l]
		for j in range(len(pop_diff_frq)):
			poplist[j]
			if check_pop_group(pop_label,poplist[l],poplist[j]):
				logger.info('Skipping pair of populations from the same Ethnic group :%s,%s' % (poplist[l],poplist[j]))
			else:
				pop_done.append(poplist[l]+"-"+poplist[j])
				list2 = convert_dict_list(pop_diff_frq[j])
				if len(list1[1]) == len(list2[1]):
					diff = get_freq_diff(list1[0],list2[0])
					diff_dict = convert_list_dict(list1[1],diff)
					if opti[3] :
 						aims_list.append(diff_dict)
						snp = list1[1]
					elif opti[2]:
						p=get_weight_pvalue(ld_g,diff_dict)
						pvalue.append(p)
		if opti[2]:
			if pvalue != []:
				for de in pvalue:
					pvalue_list[l].append(-2*log(de))
			else:
				pass
	if opti[2]:
		#Computing the overall score for each Population
		ld_g.clear()
		pop_diff_frq =[]	
		pvalue_list1 = completeness(pvalue_list,pop_label)
		pvalue_list = []
		if opti[1] == 0 or opti[1] == 0.0:
			final_score,tp = normalization(pvalue_list1,poplist,500)
			
		else:
			final_score,tp = normalization(pvalue_list1,poplist,opti[1])
		group_frq = []
		pop_label =[]
		return final_score,tp
	elif opti[3]:
		pop_diff_frq = []	
		aims_score(ld_g,aims_list,snp,out_dir,opti[1])

####################################################################################
# Unusual population differentiation test.
# difference in allele frequency on common SNP is used based on chisquare Static.
####################################################################################
def unusual_diff_allfrq(snp_diff,opti):
 	"""unusual_diff_allfrq(), inputs:minor allele frequency diference at each SNP, options."""
	un_diff_frq = []
	g_pop = []
	P_snp = []
	for pop in snp_diff:
		rows = {}
		for snp in snp_diff[pop][3]:
			da = snp_diff[pop][3].keys()
			idx = da.index(snp)
			N = float(snp_diff[pop][0][idx]) + float(snp_diff[pop][1][idx])
			p = float(snp_diff[pop][2][idx])
			P_snp.append(float(p*(1-p)))
		R = mean(P_snp)
		P_snp = []
		if opti :
			for snp in snp_diff[pop][3]:
				da = snp_diff[pop][3].keys()
				idx = da.index(snp)
				p = float(snp_diff[pop][2][idx])
				F = float(snp_diff[pop][-1])
				if p == 0.0:
					pass
				else:				
					n = float(p*(1-p))
					value = (snp_diff[pop][3][snp]**2)/n
					dof = 1
					pval = chi2.cdf(value/R, dof)
					rows[snp] = pval
			gc = 1+(N*F)	
			g_pop.append([pop,gc])		
			un_diff_frq.append([pop,rows])
		else:
			for snp in snp_diff[pop][3]:
				da = snp_diff[pop][3].keys()
				idx = da.index(snp)
				N = float(snp_diff[pop][0][idx]) + float(snp_diff[pop][1][idx])
				F = float(snp_diff[pop][-1])
				n1 = 1/float(snp_diff[pop][0][idx])
				n2 = 1/float(snp_diff[pop][1][idx])
				T = snp_diff[pop][0][idx]/float(N) + snp_diff[pop][1][idx]/float(N)
				s =1/float(N)
				p = float(snp_diff[pop][2][idx])
				fst = 2*float(snp_diff[pop][-1])
				if p != float(0.0):
					n = float((fst)*R) #G	
					q = (snp_diff[pop][3][snp])**2
					value = q/n	
					dof = 1
					pval = chi2.cdf(value, dof)
					rows[snp] = pval
				else:
					pass
			gc = 1+(N*F)
			g_pop.append([pop,gc])	
			un_diff_frq.append([pop,rows])
	return un_diff_frq,g_pop

####################################################################################
# Informativeness of Genetic Marker Score as describes in the paper.
# difference in allele frequency on common SNP is used through chisquare Static.
####################################################################################
####################################################################################
# Informativeness of Genetic Marker Score as describes in the paper.
# difference in allele frequency on common SNP is used based on chisquare Static.
####################################################################################
def aims_score(x,y,marker,out_folder,sample):
	aims = []
	aims_delta = {}
	import logging
    	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger()
	for i in range(len(y)):
		aims.append({})
	
	logger.info('Computing informativeness score ...')
	for pop in y:
		idx = y.index(pop)
		for des in x:
			rsid = des.split("-")
			delta = 0.25*(float(pop[rsid[0]]) * float(pop[rsid[1]]))
			if float(x[des][0])== 0.0:
				pass
			else:
				aims[idx][des]= delta/float(x[des][0])
	x.clear()	
	y = []
	for snp in aims[0]:
		tmp = []
		snp1= snp.split("-")[1]+"-"+snp.split("-")[0]
		for j in range(len(aims)):
			if snp in aims[j] :
				tmp.append(aims[j][snp])
			elif snp1 in aims[j]:
				tmp.append(aims[j][snp1])
			else:
				logger.info('Missing data at SNP %s' % snp)
		K = 1/float(sqrt(len(aims)))
		aims_delta[snp] = sum(tmp)*float(K)
	logger.info('Writing ...')
	fin = open(out_folder+"AIMS.score","wt")
	fin.writelines("SNP"+"\t"+"Informativeness_Score"+"\t"+"IC"+"\n")
	final_aims = {}
	f = []
	for m in marker:
		tmp = []
		for pair in aims_delta:
			p1 = pair.split("-")[0]
			p2 = pair.split("-")[1]
			if m == p1 or m == p2:
				tmp.append(aims_delta[pair])
		final_aims[m] = sum(tmp)	
		f.append(sum(tmp))
	marker = []
	aims_delta.clear()
	a1 = mean(f)
	se = std(f)
	
	for de in final_aims:
		element = final_aims[de]
		x = de
		del final_aims[de]
		C = element-std(final_aims.values())
		D = element+std(final_aims.values())
		c = "("+"%.5f" % round(C,6)+ ","+ "%.5f" %round(D,6)+")"
		final_aims[x] = element
		logger.info('Informativeness score for %s:%s' % (str(de),str(element)))
		fin.writelines(str(de)+"\t"+str(element)+"\t"+c+"\n")
	fin.close()
	final_aims.clear()

def Lexp(diff_dict,ad_ld,ld1,ld2):
	rows_ld = {}
	""" inputs: dictionary of difference in allele frequency, admixture LD in admixed popualtion, LD in first populations. """
	for des in ad_ld: 
		rsid = des.split("-")
		if des in ld2 or rsid[1]+"-"+rsid[0] in ld2:
			if des in ld1 or rsid[1]+"-"+rsid[0] in ld1:
				if rsid[0] in diff_dict and rsid[1] in diff_dict:
					if float(ld1[des][0])>=0.0 and float(ld2[des][0])>=0:	
						sam = float(diff_dict[rsid[0]] * diff_dict[rsid[1]])
						D = (float(ld1[des][0])-float(ld2[des][0])+sam)**2
						if sam == 0:
							F =0.0
							R=0.0
						else:
							R =  D/(4*sam)
							F=float(ld2[des][0])
							admax = D/(2*sam)
						rows_ld[des] = [F,R]+ad_ld[des]
					else:
						pass
				else:
					pass
			else:
				pass
	return rows_ld
def Admix_LD_Cherker(group_frq,ld_g,pop_label,out_dir):
	""" inputs: allele frequency of Group of populations, admixture LD in admixed popualtion, populations labels, output path."""
	
	pop_diff_frq = []
	poplist = []
	position = {}
	diff_rows = {}
	for de in ld_g:
		pos = abs(int(ld_g[de][1].split("-")[0]) - int(ld_g[de][1].split("-")[1]))
		position[de] = [pos,float(ld_g[de][0])]
	ld_g.clear()
	for des in group_frq:
		tmp = []
		for de in des:
			tmp += de
		pop_diff_frq += tmp
	for des in pop_label:
		poplist +=des
	#Calculate pair-wise Difference in allele frequencies from pair of reference populations
	pop_done = []
	for l in range(len(pop_diff_frq)):
		list1 = convert_dict_list(pop_diff_frq[l])
		pvalue = []
		pops=poplist[l]
		pkl_file = open(out_folder+pops+'_ld.pkl', 'rb')
		LD_P1 = pickle.load(pkl_file)
		for j in range(len(pop_diff_frq)):
			poplist[j]
			pkl_file = open(out_folder+poplist[j]+'_ld.pkl', 'rb')
			LD_P2 = pickle.load(pkl_file)
			if check_pop_group(pop_label,poplist[l],poplist[j]):
				logger.info('Skipping pair of population from the same Ethnic group :%s,%s' % (poplist[l],poplist[j]))
			else:
				if poplist[l]+"-"+poplist[j] in pop_done:
					pass
				elif  poplist[j]+"-"+poplist[l] in pop_done:
					pass
				else:	
					pop_done.append(poplist[l]+"-"+poplist[j])
					list2 = convert_dict_list(pop_diff_frq[j])
					if len(list1[1]) == len(list2[1]):
						diff = get_freq_diff(list1[0],list2[0])
						diff_dict = convert_list_dict(list1[1],diff)
						LDexp_list = Lexp(diff_dict,position,LD_P1,LD_P2)
						diff_rows[poplist[j]+"-"+poplist[l]] = LDexp_list
					else:
						sys.stdout.write("\nInconsistancy ....!")

	pop_diff_frq = []
	logger.info('Writing ... admixture LD!')
	plot_LDexp(diff_rows,out_dir)
	
	

#********************************* End Methods Functions ***********************

####################################################################################
# Pair-wise FST estimation.
####################################################################################
def pop_pairwisefst(list_all,mode):
	try:
		if list_all[1] != list_all[3]:	
			if mode : #unrelated samples
				scale_p1 = float(list_all[1]*(1-list_all[1]))/float(list_all[0])
				scale_p2 = float(list_all[3]*(1-list_all[3]))/float(list_all[2])
				mix_scale = (abs(list_all[1] - list_all[3]))**2
				scale = 2**mean([list_all[1],list_all[3]])*(1-mean([list_all[1],list_all[3]]))
				fst = (mix_scale - scale_p1 - scale_p2)/float(scale)
			else:#related samples
				mix_scale = (abs(list_all[1] - list_all[3]))**2
				scale = 2**mean([list_all[1],list_all[3]])*(1-mean([list_all[1],list_all[3]]))
				fst = mix_scale/float(scale)
			return fst
		else:
			return 0.0
	except:
		raise exceptions.SystemError('Failed to estimate population pairwise FST')

####################################################################################
# Pair-wise FST estimation.
####################################################################################
def pop_unbias_pairwisefst(list_all,mode):
	""" inputs: list of populations allele frequencies, options(related/unrelated individuals). Reterns the Fst"""
	try:
		if list_all[1] != list_all[3]:	
			if mode : #unrelated samples
				scale_p1 = float(list_all[1]*(1-list_all[1]))/float(list_all[0])
				scale_p2 = float(list_all[3]*(1-list_all[3]))/float(list_all[2])
				mix_scale = (abs(list_all[1] - list_all[3]))**2
				scale = 2**mean([list_all[1],list_all[3]])*(1-mean([list_all[1],list_all[3]]))
				fst = (mix_scale - scale_p1 - scale_p2)/float(scale)
			else:#related samples
				mix_scale = (abs(list_all[1] - list_all[3]))**2
				scale = 2**mean([list_all[1],list_all[3]])*(1-mean([list_all[1],list_all[3]]))
				fst = mix_scale/float(scale)
			return fst
		else:
			return 0.0
	except:
		raise exceptions.SystemError('Failed to estimate population pairwise FST')
####################################################################################
#These function compute the LD between two SNPs.
#The function LD_POP() computes LD in proxy ancestral population.
#The function PAIRWISE_SNP_LD() computes the LD in Admixed population.
####################################################################################
def LD_POP(tagget_dict,pos1,pop): #opt1 is LD Cutoff and opt2 is minor (maf) cutoff.
	"""inputs: genotype data of reference of populations based on retained SNP from admixed population,  positions, list of popualtions labels."""
	fin = open(out_folder+pop+'_ld.pkl','wb')
	pkl_ld = open(out_folder+'ld_admix.pkl', 'rb')
	ld_admix = pickle.load(pkl_ld)
	positions = ld_admix.keys()
	ld_admix.clear()
	ld_tagget = {}
	for des in positions:
		item_i = des.split("-")[0]
		item_j = des.split("-")[1]
		if item_i in tagget_dict and item_j in tagget_dict:
			LD = calc_rsq(tagget_dict[item_i],tagget_dict[item_j],item_i,item_j)
			if len(LD[0])== 2:
				pos = str(item_i)+"-"+str(item_j)
				ld_tagget[LD[0][0]+"-"+LD[0][1]] = [abs(LD[1]),pos]
			else:
				pos = str(item_i)+"-"+str(item_j)
				ld_tagget[LD[0][0]+"-"+LD[0][1]] = [0.0,pos]
	pickle.dump(ld_tagget,fin)
	logger.info('Computing LD ....Done %s'%pop)
	fin.close()
	# Computing pairwise SNP LD
	logger.info('Computing LD for reference population: %s'% pop)
	
def PAIRWISE_SNP_LD(filename,out_folder,op): 
	"""inputs: list of genotype and snp files, path of output and a list of options. Output: write Populations pair-wise LD and map file to be used later. """
	tagget_dict,pos1 = data_comb([],filename,op)
	news_map = {} 
	ld_tagget = {}
	if op[7] or op[8] or op[-1]:
		fin = open(out_folder+'new_map.pkl','wb')
		pickle.dump(pos1,fin)
		fin.close()		
	else:
		if op[4] :
			fin = open(out_folder+"Pairwise.ld","wt")
			fin.writelines("SNP1"+"\t"+"SNP2"+"\t"+"LD"+"\n")
		elif op[9] or op[-2]:
			fin = open(out_folder+'ld_admix.pkl','wb')
			fin1 = open(out_folder+'new_map.pkl','wb')
		else:
			fin = open(out_folder+'ld_admix.pkl','wb')
			fin1 = open(out_folder+'new_map.pkl','wb')
		positions = pos1.keys()
		# Computing pairwise SNP LD
		logger.info('Computing LD ....')
		print"\n"
		for i in range(len(positions)):
        		item_i = positions[i]
        		for j in range(i+1,min(i+op[0]+1,len(positions))):
            			item_j = positions[j]
			
				LD = calc_rsq(tagget_dict[item_i],tagget_dict[item_j],item_i,item_j)
				if len(LD[0])== 2:
					if op[9] :	
                                                if float(LD[1]) < op[3]:
                                                        pos = str(pos1[positions[i]])+"-"+str(pos1[positions[j]])
							ld_tagget[LD[0][0]+"-"+LD[0][1]] = [abs(LD[1]),pos]
                                                        news_map[LD[0][0]] = positions[i]
							news_map[LD[0][1]] = positions[j]
					
					else:
						if float(LD[1]) >= op[3]:
							if op[4]:
								fin.writelines(LD[0][0]+"\t"+LD[0][1]+"\t"+str(LD[1])+"\n")
							else:
								pos = str(pos1[positions[i]])+"-"+str(pos1[positions[j]])
								ld_tagget[LD[0][0]+"-"+LD[0][1]] = [abs(LD[1]),pos]
								news_map[LD[0][0]] = positions[i]
								news_map[LD[0][1]] = positions[j]
				else:
					if op[-2]:
						pos = str(pos1[positions[i]])+"-"+str(pos1[positions[j]])
						ld_tagget[LD[0][0]+"-"+LD[0][1]] = [0.0,pos]
						news_map[LD[0][0]] = positions[i] 
						news_map[LD[0][1]] = positions[j]
		if op[4]:
			logger.info('Writing LD on a file, Done!')
			fin.close()

		elif op[9]:
			pickle.dump(ld_tagget,fin)
			fin.close()
			pickle.dump(news_map,fin1)
			news_map.clear()
			fin1.close()
			logger.info('Writing LD on a file, Done!')
		else:
			tagget_dict.clear()
			logger.info('Saving objects ....')
			pickle.dump(ld_tagget,fin)
			pickle.dump(news_map,fin1)
			news_map.clear()
			ld_tagget.clear()
			fin.close()
			fin1.close()
			logger.info('Computing LD ....Done!')

#************************************** End LD Functions *******************************

####################################################################################	 
# Read Data inputs and Link genotype data to relative SNP and limit possible SNPs.
#The function readdata() process parameter files and inputs data.
#The function remove_missing() removes missing genotypes
#The function data_comb assess data, minor allele frequency and apply LD cutoff if 
#possible.
#The function logical() Converts a string to logical variable.
#The function scandirs() Retrieve files in the specified directory option 
#and append in a list.
#The Function random_snp() selects random SNP given a number from the entire dataset.
#The Function space_snp() selects SNPs spaced based on physical distance from genome.
####################################################################################		
def remove_missing(genotypes1,genotypes2,rsid1,rsid2):
	"""inputs: two genotypes and their related SNPs. Returns cleaned two genotypes and their related SNPs after removing missing entries."""
	if len(genotypes1) != len(genotypes2):
		sys.stdout.write("\ngenotypes should be of the same length (exiting)")
		raise SystemExit
	else:
    		n_indiv = len(genotypes1)
    	gen1 = []
    	gen2 = []
	if "?" in genotypes1:
		return gen1,gen2
	elif "?" in genotypes2 :
		return gen1,gen2
	else:
		gen1 = [int(i) for i in genotypes1]
		gen2 = [int(i) for i in genotypes2]
		return gen1,gen2

def data_comb(mapdict,filename,option):
	"""inputs" map and genotype file. Clean data by removing SNP that fails Hardy-Weinberg equilibrium."""
	tagget_dict = {}
	if len(mapdict) == 0:
		map_tmp = {}
		for line in fileinput.input(filename[0]):
			data = line.split()
			map_tmp[fileinput.lineno()] = data
		tagget_dict = {}
		pos = {}
		for line in fileinput.input(filename[1]):
			data=line.split()
			try:
				m_snp = map_tmp[fileinput.lineno()][0]
				if option[9] or option[10]:
					tagget_dict[m_snp]= data[0] 
					pos[m_snp] = int(map_tmp[fileinput.lineno()][3])
				else:
					frq,all_count,minor,m_count = get_freq(data[0]) 
					if minor > option[2]:
						tagget_dict[m_snp]= data[0]
						pos[m_snp] = int(map_tmp[fileinput.lineno()][3])			
			except:
				raise exceptions.SystemError('Failed to process, key not found')
		if option[9] or option[10] :
			dict_final = {}
			pos_final = {}
			logger.info('Spacing SNPS at %s cMorgan interval' % str(option[6]))
			SNP_list = space_snp(option[6],pos)
			for des in tagget_dict.keys():
				if des in SNP_list:
					dict_final[des] = tagget_dict[des]
					pos_final[des] = pos[des]
			tagget_dict.clear()
			pos.clear()
			map_tmp.clear()	
			return dict_final,pos_final
		elif option[11]:
			return tagget_dict,pos
		else:
			## Random selection of SNPs
			SNP_list = random_snp(option[1],option[6],pos)
			dict_final = {}
			pos_final = {}
			for des in tagget_dict.keys():
				if des in SNP_list:
					dict_final[des] = tagget_dict[des]
					pos_final[des] = pos[des]
			tagget_dict.clear()
			pos.clear()	
			return dict_final,pos_final
	else:
		map_tmp = {}
		if option[11]:
			pos ={}
		for line in fileinput.input(filename[0]):
			data = line.split()
			if data[0] in mapdict:
				map_tmp[fileinput.lineno()] = data
		tagget_dict = {}
		for line in fileinput.input(filename[1]):
			data=line.split()
			try:
				if fileinput.lineno() in map_tmp:
					tagget_dict[map_tmp[fileinput.lineno()][0]]=data[0]
					if option[11]:
						pos[map_tmp[fileinput.lineno()][0]] = int(map_tmp[fileinput.lineno()][3])
			except:
				raise exceptions.SystemError('Failed to process, key not found ')
		map_tmp.clear()
		if option[11]:
			return tagget_dict,pos
		else:	
			return tagget_dict

# Space SNPs based on physical distance from the entire genome, the function receive un physical interval in cM and dictionary of SNPs.
def space_snp(opt,pos):

	"""inputs: list of SNPs within their position. Return list of SNPs spaced at a defined cM."""
	sys.stdout.write("\nSpacing SNPs ......")
	chose = {}
	sort_snp ={}
	space = int(opt*1000000)
	for snp in pos:
		sort_snp[pos[snp]] = snp
	position = sort_snp.keys()
	position.sort()
	tmp_snp = []
	for i in range(len(position)):
		if i == 0:
			tmp_snp.append(position[i])
			space = position[i] + space
		else:
			if position[i] >= space:
				tmp_snp.append(position[i])
				space = position[i] + space
	for des in position:
		chose[sort_snp[des]] = des
	return chose

# Random selection of SNP from the entire dataset
def random_snp(SnpN,space,listSNP):
	"""inputs: list of SNPs within their position. Return list of a random selected SNPs ."""
	#Select random markers from the given number using random functiom
	RANDOMSEED = 0
	chosen = []
	if int(SnpN) == 0:
		logger.info('SNP will be spaced in %s SNPs .....' % str(space))
		print"\n"
		if len(listSNP.keys()) >= 15000 :
			return space_snp(space,listSNP)
		else:
			logger.info('SNPs are not sufficient to be spaced. All %d will be used ' % len(listSNP))
			print"\n"
			return listSNP.keys()
	else:
		if len(listSNP) <= SnpN:
			logger.info('All SNP will be used %d SNPs .....' % len(listSNP))
			print"\n"
			return listSNP.keys()
		else:
			logger.info('Initial Number of SNPs:%s'%len(listSNP))
			print"\n",SnpN
			ListSNPs=listSNP.keys()
        		if RANDOMSEED is None:
                		rand = Random()
        		else:
             			rand = Random(RANDOMSEED)
			while len(chosen) < SnpN:
				idx=rand.randint(0,len(ListSNPs))
				if idx in range(len(ListSNPs)):
					chosen.append(ListSNPs[idx])
                			logger.info('Including %d SNPs .....' % len(chosen))
			return chosen



# The function receive a path and return a list of files in the directory (old version)
#receive the main pops directory in which are located subdirectories of group of pops and N=Subgroup pop (New version)
def scandirs(path,N): 
	"""inputs: working directrory, admixture way that determined number of reference popualtion pools.  """
	pop_group_files = []
	pop_names = []
	for i in range(N):
		pop_group_files.append([])
		pop_names.append([])
	DIR = []
	for dirname, dirnames, filenames in os.walk(path):
		DIR.append(dirname)
	if N == len(DIR)-1:
		for paths in DIR[1:]:
			idx = DIR[1:].index(paths)
			dirList=os.listdir(paths)
			for fname in dirList:
				pop_group_files[idx].append(fname)
				pop_names[idx].append(fname.split(".")[0])
	
		return pop_group_files,pop_names,DIR[1:]
	else:
		logger.info('Way-Admixture is different to number of group ancestries! Check parameters:Subpop!')
		logger.info('Exiting at ...: %s'% str(datetime.datetime.today()))
		raise SystemExit
		
def logical(string):
	if string == "NO":
		return False
	elif string == "YES":
		return True
	else:
		return False

def readdata():
	""" Read instructions from parameter file. Returns user options."""
	popnames = []
	if len(sys.argv) == 1:
		logger.info('Starting at time:%s' % str(datetime.datetime.today()))
		logger.info('Command line usage: %s <parameterfile>  ' % sys.argv[0])
		raise SystemExit
	else:
		try:
			logger.info('Starting at time:%s' % str(datetime.datetime.today()))
			INFILE = sys.argv
			rows ={}				
			for line in fileinput.input(INFILE[1]):
				data = line.split()
				if not line.startswith("#"):
					if len(data) != 0:
						rows[data[0].split(":")[0]] = data[0].split(":")[1]
			my_pop = rows['admixD']
			popgenofiles,popnames,path=scandirs(rows['popDD'],int(rows['subpop']))#Retrieving files from Anc groups		
			markerfiles=rows['SNPd']
			opt = [int(rows['windows']),int(rows['subSNP']),int(rows['subpop']),float(rows['maf']),float(rows['LDcutoff']),logical(rows['relatedInd']),logical(rows['Fstscore']),logical(rows['Frqscore']),logical(rows['pwFst']),logical(rows['Ufrq']),logical(rows['writeLD']),logical(rows['WeightedLD']),logical(rows['aims']),logical(rows['PCAaims']),int(rows['Sampling']),float(rows['cMspace']),logical(rows['no-log']),str(rows['colors']),logical(rows['lines']),logical(rows['Lexp']),str(rows['outp'])]
			logger.info('Reading data is Done ...') 
			return my_pop,popgenofiles,markerfiles,opt,popnames,path
		except:
			raise exceptions.SystemError('Failing to process the input data ...')



#************************************** End Processing Data input Functions *******************************


####################################################################################
# Main function of PROXYANC program
####################################################################################
def proxy_anc():
	""" Main function to computer user otions."""

	tagget_pop,inputfiles,mapfile,opts,popnames,path = readdata()
	filename = [mapfile,tagget_pop] 
	global out_folder
	out_folder = opts[-1]
	
	#Create output folder if does not exist
	if not os.path.exists(out_folder):
    		os.makedirs(out_folder)
	
	tagg_frq = []
	tagg_all = []
	tagg_der = []

	opt1 = [opts[0],opts[1],opts[3],opts[4],opts[10],opts[11],opts[15]]
	opt2 = [opts[6],opts[8],opts[12],opts[13],opts[-2],opts[9]]
	"""
	# opts[0]: An integer specifies the Windows side to compute the LD.
	# opts[1]: An integer speicifies the number of random subset of SNPs to be chosen for the analysis.
	# opts[2]: Number of ancestry component in the admixed populations should equal to the number of groups (gathering in diffrent folder) of populations.
	# opts[3]: Float specifies the cutoff of Minor allele frequency to be considered.
	# opts[4]: Float specifies the cutoff of LD.
	# opts[4]: Boolean (YES/NO) specifies whether the population samples under studies are related or unrelated.
	# opts[6]: Boolean (YES/NO) specifies whether the Fst Optimal Cone programming method to select best proxy ancestry.
	# opts[7]: Boolean (YES/NO) specifies whether the proxy-ancestry score method, to select best proxy ancestry.
	# opts[8]: Boolean (YES/NO) specifies whether to write popualtion pair-wise Fst.
	# opts[9]: Boolean (YES/NO) specifies whether to investigate popualtion pair-wise unusual difference in allele frequency as signal of natual selection.
	# opts[10]: Boolean (YES/NO) specifies whether to write the LD in a file.
	# opts[11]: Boolean (YES/NO) specifies whether to consider a weighted model of LD
	# opts[12]: Boolean (YES/NO) specifies whether to select AIMS panel of the admixed population based on each single best proxy ancestral from each group.
	# opts[13]: Boolean (YES/NO) specifies whether to select AIMS panel using a Kernel PCA method.
	# opts[14]: An integer specifies bootstrap sampling.
	# opts[15]: A float specifies physical space between adjant markers to limit the background LD.
	# opts[16]: Boolean (YES/NO) specifies whether -log10(p) to be applied on the value in manhattan from unusual difference in allele frequencies (active if option Ufrq if TRUE) method.
	# opts[17]: Boolean (YES/NO) specifies the cycle through colors in manhattan from unusual difference in allele frequencies (active if option Ufrq if TRUE) method.
	# opts[18]: Boolean (YES/NO) whether manhattan to be lines mode.
	# opts[-2]: Boolean (YES/NO) specifies whether to examine Admixed LD based on the expected maxiumun admixture LD from pair-wise proxy ancestral populations.
	# ops[-1]: Path to where direct all outputs. """
	opt = opt1 + opt2

	METHOD = [opts[6],opts[7],opts[8],opts[12],opts[13],opts[-2],opts[9]]

	# Allow just one method at once
	if METHOD.count(True) == 2:
		logger.info('Starting at time:%s' % str(datetime.datetime.today()))
		logger.info('Multiple Methods choice in parameter file, exiting...! ')
		raise SystemExit

	#Compute LD in (or admixed) population understudy dataset
	if opts[6] or opts[8] or opts[10] or opts[9]:
		PAIRWISE_SNP_LD(filename,out_folder,opt)	
	else:
		if opts[12] or opts[7] or opts[-2]:
			#Compute pairwise LD in admixed pop, and reduce the SNP set based on relavant LD
			PAIRWISE_SNP_LD(filename,out_folder,opt) 
			LD ={} 
		All_frq_diff = [] # List of population groups allele frequecies
		for i in range(len(inputfiles)):
			All_frq_diff.append([])
		for i in range(len(popnames)):
			for j in range(len(popnames[i])):
				All_frq_diff[i].append([])

	"""Select genotypes in ancestral populations based on the reduced SNP sets after LD and minor allele freq cut-off from population under study dataset."""
	if opts[6] or opts[8] or opts[9]:
		fst_frq_pop1 = []
		pkl_file = open(out_folder+'new_map.pkl', 'rb')
		pos1 = pickle.load(pkl_file)
		if opts[8] or opts[9]:
			pop_sample = []
			pop_files = []
			pop_sample.append(tagget_pop.split("/")[-1].split(".")[0])
			for i in range(len(popnames)):
				for j in range(len(popnames[i])): 
					pop_sample.append(popnames[i][j])
					pop_files.append(path[i]+"/"+inputfiles[i][j])
			All_N = []
			All_P = []
			SNP = {}
			if opts[8]:
				popnames =[]
			for i in range(len(pop_sample)-1):
				All_N.append([])
				All_P.append([])
				fst_frq_pop1.append([])

			#Compute ancestral populations allele frequency from each group. 
			for filename1 in pop_files:
				idx = pop_files.index(filename1)
				popfile = [mapfile,filename1]
				pop_dict = data_comb(pos1,popfile,opt)
				logger.info('Calculating the allele frequencies for population:%s' % pop_sample[idx])
				for rsid in pop_dict:
					frq,count_all,minor,minor_count = get_freq(pop_dict[rsid])
					fst_frq_pop1[idx].append(minor)
					All_N[idx].append(count_all)
					All_P[idx].append(minor_count)
					if idx == 0 :
						SNP[rsid] = minor
		else:
			for i in range(len(inputfiles)):
				fst_frq_pop1.append([])
			for i in range(len(popnames)):
				for j in range(len(popnames[i])):
					fst_frq_pop1[i].append([])

		#Compute allele frequency in the (admixed) popualtion. 
		tagg_dict = data_comb(pos1,filename,opt)
		for rsid in tagg_dict:
			frq,all_count,minor,minor_count = get_freq(tagg_dict[rsid])
			tagg_frq.append(frq)
			tagg_all.append(all_count)
			tagg_der.append(minor_count)

	""" Select genotypes in ancestral populations based on the reduced SNP sets after LD and minor allele freq cut-off from population under study dataset.
	   using 		
	"""
	if opts[6] or opts[12] or opts[7] or opts[-2]: 
		"""Allele frequencies in each putative ancestral population, using the retained SNPs LD>threhold from the popualtion under study"""
		pkl_file1 = open(out_folder+'new_map.pkl', 'rb')
		news_map = pickle.load(pkl_file1)
		subop = [opts[6],opts[12],opts[7],opts[-2]]
		for group in inputfiles:
			gr = inputfiles.index(group)
			for filename in group:
				tmp = []
				idx = group.index(filename)
				filename1 = path[gr]+"/"+filename
				popfile = [mapfile,filename1]
				if opts[7] or opts[12] :
					pop_dict1 = data_comb(news_map,popfile,opt)
					logger.info('Calculating the allele frquencies for population:%s' % popnames[gr][idx])
					D1, c1, c2 = allele_frq_process(pop_dict1,subop)
					All_frq_diff[gr][idx].append(D1.copy())
					D1.clear()
				elif opts[-2]:
					pop_dict1,PO = data_comb(news_map,popfile,opt)
					LD_POP(pop_dict1,PO,popnames[gr][idx])
					logger.info('Calculating the allele frquencies for population:%s' % popnames[gr][idx])
					D1, c1, c2 = allele_frq_process(pop_dict1,subop)
					All_frq_diff[gr][idx].append(D1.copy())
					D1.clear()
				else:
					pop_dict1 = data_comb(news_map,popfile,opt)
					logger.info('Calculating the allele frquencies for population:%s' % popnames[gr][idx])
					if opts[6]:
						D1, c1, c2 = allele_frq_process(pop_dict1,subop)
						fst_frq_pop1[gr][idx].append([popnames[gr][idx],c1,c2])	
		pop_dict1.clear()
		inputfiles = []

#***********############## Run function of difference Methos options ########*************#
	print"\n"
	if opts[7] :##  proxy-ancestry score method
		logger.info('***** Ancestry Differentiation Frequency Score Methods. **********"' )
		OPT = [opts[2],opts[14],opts[7],opts[12]]	
		pkl_ld = open(out_folder+'ld_admix.pkl', 'rb')
		ld_tagget = pickle.load(pkl_ld)
		news_map.clear()

		# Calling  proxy-ancestry score main function 
		final_score,tp = FRQ_score_anc(All_frq_diff,ld_tagget,popnames,out_folder,OPT)

		ld_tagget.clear()
		Labbel = {}
		
		# Deleting temporary and unnecessary saved files
		os.system('rm'+' '+out_folder+'ld_admix.pkl')
		os.system('rm'+' '+out_folder+'new_map.pkl')
		All_frq_diff = []

		for p in popnames:
			idx = popnames.index(p)
			for i in p:
				Labbel[i] = idx
		
		# printing each ancestral population score, per group
		logger.info('******************Population Z score************')
		print"\n"
		print"\n","POPs","nSCORE","PSCORE"
		print"\n","Pops"," ","Norm Zscore"," ","Confidence"," ","Zscore"
		for h in range(len(tp)):
			print"\n",final_score[h][0],final_score[h][1],final_score[h][2],tp[h]

		#Writing proxy ancestry score into an output file
		fix = open(out_folder+"proxyanc.out","wt")
		fix.writelines("********** Population FRQ Score *********"+"\n")
		fix.write("\n")
		fix.writelines("POPs"+"\t"+"pScore"+"\t"+"SE"+"\t"+"Z"+"\t"+"95%CI"+"\n")
		step = len(popnames[0])
		gr = 1
		for pop in final_score:
			ifx = final_score.index(pop)
			if step == Labbel[pop[0]] and ifx == 0:
				fix.writelines("\t"+"\t"+"Population" +" " +"Group:"+" "+str(gr)+"\n")
				fix.writelines(str(pop[0])+"\t" +"%.3f" %float(pop[1])+"\t"+"%.3f" %float(pop[2][1])+"\t"+"%.3f" %float(tp[ifx])+"\t"+str(pop[2][0])+"\n")	
				gr = gr + 1
				step = 	Labbel[pop[0]]
			elif step == Labbel[pop[0]] and ifx != 0:
				fix.writelines(str(pop[0])+"\t" +"%.3f" %float(pop[1])+"\t"+"%.3f" %float(pop[2][1])+"\t"+"%.3f" %float(tp[ifx])+"\t"+str(pop[2][0])+"\n")
			else:
				if ifx !=0:
					gr = gr + 1
					step = 	Labbel[pop[0]]
				elif ifx == 0:
					step = 	Labbel[pop[0]]
				fix.writelines("\t"+"\t"+"Population" +" " +"Group:"+" "+str(gr)+"\n")
				fix.writelines(str(pop[0])+"\t" +"%.3f" %float(pop[1])+"\t"+"%.3f" %float(pop[2][1])+"\t"+"%.3f" %float(tp[ifx])+"\t"+"\t"+str(pop[2][0])+"\n")
					
		fix.write("\n")
		fix.write("PROXYANC* Finished at:"+"\t"+str(datetime.datetime.today()) +"\n")
		fix.close()

	if opts[9] : #Unsual Allele frqs differentiation
		logger.info('************** Running Unsual Allele frequency Differentiation Test. **************')
		OPT = [2,opts[5],opts[14],opts[9]]
		pos1.clear()
		os.system('rm'+' '+out_folder+'new_map.pkl')

		# call the main function to test unusual difference in allele freq for each pair of populations. 
		pop_fst = pwFST([tagg_all,tagg_frq],fst_frq_pop1,All_N,pop_sample,SNP,out_folder,OPT)
		map_snp = {}
		logger.info('starting chi2 test ...')
		unusual_frq,g_gc = unusual_diff_allfrq(pop_fst,opts[5])
		GC_inflat = []
		pair_pop = ["CHR","SNP","POS","A1","A2","P","adj_G1","adj_G2"]
		result = {}
		for line in fileinput.input(mapfile):
			data = line.split()
			map_snp[data[0]] = data

		#Computing genomic control from resulting chi2 test
		for data in unusual_frq:
			idx = unusual_frq.index(data)
			tmp =[]
			pva = {} 
			pva_c = {}
			for snp in data[1]:
				tmp.append(data[1][snp])
			if g_gc[idx][0] == data[0]:
				gc1 = median(tmp)/0.455
				gc2 = float(g_gc[idx][1])
			else:
				logger.info('NO Match ....!!!!')

			#Wrting SNP p-values of  chi2 test of unusual difference in allele freq for each pair of populations.
			fin = open(out_folder+data[0]+".Unsual.Frq.out","wt")
			fin.write("\t".join(pair_pop)+"\n")
			logger.info('Writing ....%s'%data[0])
			for snp in data[1]:
				p_c = data[1][snp]/float(gc1)
				p_c1 = data[1][snp]/float(gc2)
				fin.writelines(map_snp[snp][1]+"\t"+map_snp[snp][0]+"\t"+map_snp[snp][3]+"\t"+map_snp[snp][4]+"\t"+map_snp[snp][5]+"\t"+str(1-data[1][snp])+"\t"+str(1-p_c)+"\t"+str(1-p_c1)+"\n")
			fin.close()

			# Plotting manhattan plot for unusual difference in allele freq for each pair of populations.
			fhs = open(out_folder+data[0]+".Unsual.Frq.out")
			columns = [0,2,5]
			#Create output folder if does not exist
			if not os.path.exists(out_folder+"Manhattan_Plot/"):
    				os.makedirs(out_folder+"Manhattan_Plot/")
			fig_out = out_folder+"Manhattan_Plot/"+data[0]+".Mahanatha.png"
			manhattan(fhs, columns,fig_out, opts[16], opts[17], "\t"," ", opts[18], None)	

	if opts[6]:#FST Cone quadratic programming: Ancestry Proxy score opt[6]
		logger.info('************** Running the FST Score Methods. **************')

		if int(opts[14]) == 0:
			logger.info('The sampling for bootstrap is set to %s in parameter file' %str(opts[14]))
			logger.info('Setting defaut sampling to %s' %str(1000))
			OPT = [opts[2],opts[5],0,False]
		else:
			OPT = [opts[2],opts[5],opts[14],False]
		pos1.clear()

		# Calling Fst Optimal Cone programming method to select best proxy ancestry.	
		FST_Score_anc([tagg_all,tagg_frq,tagg_der],fst_frq_pop1,out_folder,OPT)
		os.system('rm'+' '+out_folder+'*.pkl')	

	if opts[8]: #Pop pairwise FST opt[8]
		logger.info('****** Running PairWise Populations Genetic Distance. ********')
		OPT = [2,opts[5],opts[14],opts[9]]
		logger.info('Populations included:%d' % len(pop_sample)),
		pos1.clear()

		# Calling the main function to compute population pair-wise FST.
		pwFST([tagg_all,tagg_frq],fst_frq_pop1,All_N,pop_sample,SNP,out_folder,OPT)
		os.system('rm'+' '+out_folder+'*.pkl')

	if opts[12]:
		logger.info('********** Selecting AIMS Panel.**********')
		OPT = [opts[2],opts[14],opts[7],opts[12]]
		news_map.clear()

		pkl_ld = open(out_folder+'ld_admix.pkl', 'rb')
		LD = pickle.load(pkl_ld)
		# Calling the main funtion of rpoxy ancestry, which also select AIMs panel of an admixed popualtion.
		FRQ_score_anc(All_frq_diff,LD,popnames,out_folder,OPT)
		LD.clear()
		os.system('rm'+' '+out_folder+'*.pkl')

	if opts[13]: # Kernel PCA-based method to select AIMs panel for admixed population
		logger.info('********* Running Kernel PCA for Selecting AIMS Panel. *********')
		sys.stderror.write("not yet implemented")
		sys.exit(1)

	if opts[-2]: #the expected maxiumun admixture LD from pair-wise proxy ancestral populations
		logger.info('********* Running Admixture LD Cherker . *********')
		pkl_ld = open(out_folder+'ld_admix.pkl', 'rb')
		ld_tagget = pickle.load(pkl_ld)
		news_map.clear()

		# examine Admixed LD based on the expected maxiumun admixture LD from pair-wise proxy ancestral populations.
		Admix_LD_Cherker(All_frq_diff,ld_tagget,popnames,out_folder)
		os.system('rm'+' '+out_folder+'*.pkl')	
		

if __name__ == '__main__':
	print "\n************************************************************************************"
	print "	Program for selecting best proxy ancestry for recently admixed population"
	print "			       Emile Chimusa"
	print "				Beta Version.						\n"
	print "			 2011, UNIVERSITY OF CAPE TOWN"
	print "***************************************************************************************\n"
	global out_folder
	proxy_anc()
	logger.info('Finish at time:%s' % str(datetime.datetime.today()))
	
