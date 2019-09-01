#!/usr/bin/env python
# coding: utf-8

# In[2]:


#loading required modules for metabolome script

import numpy as np
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import scipy.stats


# In[3]:


#load the metabolome data file
#this is a reduced data file, with only mass/charge values and abundance per sample
met_all = open('Metabolome_data_mz_and_abundance.csv', 'rU').read() 


# In[4]:


#here I split the text file into many text files (strings) corresponding to the rows
met_all2 = met_all.split('\n') [:-1] 


# In[5]:


#here i split the values by ',' for each 'row' of metabolome data
met_all3 = [i.split(',') for i in met_all2]


# In[6]:


#here i define a new class, which allows the samples to be split into their relevant groups (by tissue and age). 
#self.B_anv etc. runs a one way ANOVA on the defined groups
#self.B_fold_ etc. calculates the fold change between each group
#B = brain, L = liver, M = muscle, H = heart, S = serum

class a9a(object):
    
    def __init__(self, inp, j):
        
        self.inp = inp    
        self.j = j
        
        self.B6 = [float(self.inp[self.j][i]) for i in range(1,5) ]
        self.B10 = [float(self.inp[self.j][i]) for i in range(5,9) ]
        self.B16 = [float(self.inp[self.j][i]) for i in range(9,13) ]
        self.L6 = [float(self.inp[self.j][i]) for i in range(13,17) ]
        self.L10 = [float(self.inp[self.j][i]) for i in range(17,21) ]
        self.L16 = [float(self.inp[self.j][i]) for i in range(21,25) ]  
        self.M6 = [float(self.inp[self.j][i]) for i in range(25,29) ]
        self.M10 = [float(self.inp[self.j][i]) for i in range(29,33) ]
        self.M16 = [float(self.inp[self.j][i]) for i in range(33,37) ]
        self.H6 = [float(self.inp[self.j][i]) for i in range(37,41) ]
        self.H10 = [float(self.inp[self.j][i]) for i in range(41,45) ]
        self.H16 = [float(self.inp[self.j][i]) for i in range(45,49) ]
        self.S6 = [float(self.inp[self.j][i]) for i in range(49,53) ]
        self.S10 = [float(self.inp[self.j][i]) for i in range(53,57) ]
        self.S16 = [float(self.inp[self.j][i]) for i in range(57,61) ]
            
       
        self.B_anv = scipy.stats.f_oneway(self.B6, self.B10, self.B16)
        self.L_anv = scipy.stats.f_oneway(self.L6, self.L10, self.L16)
        self.M_anv = scipy.stats.f_oneway(self.M6, self.M10, self.M16)
        self.H_anv = scipy.stats.f_oneway(self.H6, self.H10, self.H16)
        self.S_anv = scipy.stats.f_oneway(self.S6, self.S10, self.S16)
        
        
        self.B_fold_6_10 = np.log2(np.average(self.B10)/np.average(self.B6))
        self.B_fold_10_16 = np.log2(np.average(self.B16)/np.average(self.B10))
        self.B_fold_6_16 = np.log2(np.average(self.B16)/np.average(self.B6))
        
        self.L_fold_6_10 = np.log2(np.average(self.L10)/np.average(self.L6))
        self.L_fold_10_16 = np.log2(np.average(self.L16)/np.average(self.L10))
        self.L_fold_6_16 = np.log2(np.average(self.L16)/np.average(self.L6))
        
        self.M_fold_6_10 = np.log2(np.average(self.M10)/np.average(self.M6))
        self.M_fold_10_16 = np.log2(np.average(self.M16)/np.average(self.M10))
        self.M_fold_6_16 = np.log2(np.average(self.M16)/np.average(self.M6))
        
        self.H_fold_6_10 = np.log2(np.average(self.H10)/np.average(self.H6))
        self.H_fold_10_16 = np.log2(np.average(self.H16)/np.average(self.H10))
        self.H_fold_6_16 = np.log2(np.average(self.H16)/np.average(self.H6))
        
        self.S_fold_6_10 = np.log2(np.average(self.S10)/np.average(self.S6))
        self.S_fold_10_16 = np.log2(np.average(self.S16)/np.average(self.S10))
        self.S_fold_6_16 = np.log2(np.average(self.S16)/np.average(self.S6))


# In[7]:


#here i create a new function which generates empty lists for ANOVA results and fold changes, and use the class
#created above (a9a) to generate pvalues and fold changes. The empty lists are then filled with these results
def lp_b(inp):
    pval = []
    log2fc6_10 = []
    log2fc10_16 =[]
    log2fc6_16 = []
    for j in range(1, len(inp)):
        ls = a9a(inp, j)
        pval.append(ls.B_anv[-1])
        log2fc6_10.append(ls.B_fold_6_10)
        log2fc10_16.append(ls.B_fold_10_16)
        log2fc6_16.append(ls.B_fold_6_16)
    return pval, log2fc6_10, log2fc10_16, log2fc6_16

def lp_l(inp):
    pval = []
    log2fc6_10 = []
    log2fc10_16 =[]
    log2fc6_16 = []
    for j in range(1, len(inp)):
        ls = a9a(inp, j)
        pval.append(ls.L_anv[-1])
        log2fc6_10.append(ls.L_fold_6_10)
        log2fc10_16.append(ls.L_fold_10_16)
        log2fc6_16.append(ls.L_fold_6_16)
    return pval, log2fc6_10, log2fc10_16, log2fc6_16

def lp_m(inp):
    pval = []
    log2fc6_10 = []
    log2fc10_16 =[]
    log2fc6_16 = []
    for j in range(1, len(inp)):
        ls = a9a(inp, j)
        pval.append(ls.M_anv[-1])
        log2fc6_10.append(ls.M_fold_6_10)
        log2fc10_16.append(ls.M_fold_10_16)
        log2fc6_16.append(ls.M_fold_6_16)
    return pval, log2fc6_10, log2fc10_16, log2fc6_16

def lp_h(inp):
    pval = []
    log2fc6_10 = []
    log2fc10_16 =[]
    log2fc6_16 = []
    for j in range(1, len(inp)):
        ls = a9a(inp, j)
        pval.append(ls.H_anv[-1])
        log2fc6_10.append(ls.H_fold_6_10)
        log2fc10_16.append(ls.H_fold_10_16)
        log2fc6_16.append(ls.H_fold_6_16)
    return pval, log2fc6_10, log2fc10_16, log2fc6_16

def lp_s(inp):
    pval = []
    log2fc6_10 = []
    log2fc10_16 =[]
    log2fc6_16 = []
    for j in range(1, len(inp)):
        ls = a9a(inp, j)
        pval.append(ls.S_anv[-1])
        log2fc6_10.append(ls.S_fold_6_10)
        log2fc10_16.append(ls.S_fold_10_16)
        log2fc6_16.append(ls.S_fold_6_16)
    return pval, log2fc6_10, log2fc10_16, log2fc6_16


# In[8]:


#this runs the newly generated functions from above (e.g lp_b, lp_l) and labels the output lists with relevant name
#e.g pvals_b for the ANOVA p-values for brain tissue
pvals_b, log2fc6_10_b, log2fc10_16_b, log2fc6_16_b = lp_b(met_all3)
pvals_l, log2fc6_10_l, log2fc10_16_l, log2fc6_16_l = lp_l(met_all3)
pvals_m, log2fc6_10_m, log2fc10_16_m, log2fc6_16_m = lp_m(met_all3)
pvals_h, log2fc6_10_h, log2fc10_16_h, log2fc6_16_h = lp_h(met_all3)
pvals_s, log2fc6_10_s, log2fc10_16_s, log2fc6_16_s = lp_s(met_all3)


# In[9]:


#here i create a new list of pvalues for each tissue, which excludes any pvalues that are 'not a number' (nan)

pvals2_b = [i for i in pvals_b if str(i) != 'nan']
pvals2_l = [i for i in pvals_l if str(i) != 'nan']
pvals2_m = [i for i in pvals_m if str(i) != 'nan']
pvals2_h = [i for i in pvals_h if str(i) != 'nan']
pvals2_s = [i for i in pvals_s if str(i) != 'nan']


# In[10]:


#here i plot the p-value distribution for each tissue using the list of p-values generated above

plt.figure(figsize=(15,15))

plt.subplot(3,2,1)
plt.subplots_adjust(top=1.5, hspace=0.8, wspace=0.5)
n, bins, patches = plt.hist(pvals2_b, 50, normed=1, facecolor='r', alpha=0.75, histtype='stepfilled') 
#plt.xlim(0,1.2)
#plt.ylim(0,10)
plt.title('p value distribution, brain only')

##

plt.subplot(3,2,2)
plt.subplots_adjust(top=1.5, hspace=0.8, wspace=0.5)
n, bins, patches = plt.hist(pvals2_l, 50, normed=1, facecolor='b', alpha=0.75, histtype='stepfilled') 
#plt.xlim(0,1.2)
#plt.ylim(0,8)
plt.title('p value distribution, liver only')

##

plt.subplot(3,2,3)
plt.subplots_adjust(top=1.5, hspace=0.8, wspace=0.5)
n, bins, patches = plt.hist(pvals2_m, 50, normed=1, facecolor='b', alpha=0.75, histtype='stepfilled') 
#plt.xlim(0,1.2)
#plt.ylim(0,8)
plt.title('p value distribution, muscle only')

##

plt.subplot(3,2,4)
plt.subplots_adjust(top=1.5, hspace=0.8, wspace=0.5)
n, bins, patches = plt.hist(pvals2_h, 50, normed=1, facecolor='b', alpha=0.75, histtype='stepfilled') 
#plt.xlim(0,1.2)
#plt.ylim(0,4)
plt.title('p value distribution, heart only')

##

plt.subplot(3,2,5)
plt.subplots_adjust(top=1.5, hspace=0.8, wspace=0.5)
n, bins, patches = plt.hist(pvals2_s, 50, normed=1, facecolor='b', alpha=0.75, histtype='stepfilled') 
#plt.xlim(0,1.2)
#plt.ylim(0,6)
plt.title('p value distribution, serum only')

plt.savefig('pvalue_distribution.pdf')


# In[11]:


# Home-made Benjamini Hochberg threshold. this returns only the threshold value
def BH(pvals, fdr):
#this defines a new function 'BH' with parameters 'pvals' and 'fdr'
    sortedp = sorted(pvals)
#this defines 'sortedp' as the sorted list of input 'pvals', which in our case is the pvals list 'pvals2_b'
    return max([sortedp[j] for j in range(len(pvals)) if float(j)/len(pvals)*fdr >= sortedp[j]])
#if float(j)/# of pvalues*fdr(input by user) is greater than the pvalue, then return the max value. do this for each line in pvals.
  


# In[12]:


# Home-made Benjamini Hochberg calculation. This returns the entire list of items included at calculated threshold.
def BH_2(pvals, fdr):
    sortedp = sorted(pvals)
    return ([sortedp[j] for j in range(len(pvals)) if float(j)/len(pvals)*fdr >= sortedp[j]])


# In[13]:


#qh_01 uses the BH function to calculate a threshold at a given FDR, for a list of pvalues 'BH(pvalue_list, FDR)'
qh_01 = BH(pvals2_b, 0.1) #FDR of 0.1
#qh01_ind generates an index list for all items included at calculated BH threshold
qh01_ind = [pvals_b.index(i) for i in pvals_b if i < qh_01]
#prova runs the a9a class function on met_all3, for items included in the threshold set by the BH function
prova = a9a(met_all3, qh01_ind[0]+1)
#here i create a new data array using only the data generated by the prova function above
data = map(np.array, [prova.B6, prova.B10, prova.B16])


# In[14]:


#Home-made q-value estimate function. Here I define a new function 'q' that estimates q-values, at a specified fdr
def q(pvals, fdr):
    fdr = float(fdr)
    len([i for i in pvals if i > fdr])
    return 1-float(len([i for i in pvals if i > fdr]))/float((len(pvals)*(1-fdr)))


# In[15]:


q(pvals_b, 0.1)


# In[16]:


#gives a conservative estimate of proportion of truly alternative tests, of those that are significant
print q(pvals_b, 0.1)
print q(pvals_l, 0.1)
print q(pvals_m, 0.1)
print q(pvals_h, 0.1)
print q(pvals_s, 0.1)


# In[17]:


#here i create a new list of p-values sorted in numerical order
sorted_pvals_b = sorted(pvals_b)
sorted_pvals_l = sorted(pvals_l)
sorted_pvals_m = sorted(pvals_m)
sorted_pvals_h = sorted(pvals_h)
sorted_pvals_s = sorted(pvals_s)


# In[18]:


#the following 11 cells are used to calculate the ideal q-value threshold for our data
#they are not necessary for producing the output data


# In[19]:


#need to calculate estimated
#here i create a new variable t_values filled with data points from 0 to 1, with interval 0.0001
t_values = np.arange(0,1,0.0001)

#calculating lambda for each tissue
lam_b = [(len([i for i in sorted_pvals_b if i > j])) / ((len(sorted_pvals_b)) * (1-j)) for j in t_values]
lam_l = [(len([i for i in sorted_pvals_l if i > j])) / ((len(sorted_pvals_l)) * (1-j)) for j in t_values]
lam_m = [(len([i for i in sorted_pvals_m if i > j])) / ((len(sorted_pvals_m)) * (1-j)) for j in t_values]
lam_h = [(len([i for i in sorted_pvals_h if i > j])) / ((len(sorted_pvals_h)) * (1-j)) for j in t_values]
lam_s = [(len([i for i in sorted_pvals_s if i > j])) / ((len(sorted_pvals_s)) * (1-j)) for j in t_values]


# In[20]:


#plotting lambda for each tissue
plt.subplot(2,3,1)
plt.plot(t_values, lam_b)
plt.ylabel("lambda(pi)")
plt.xlabel("lambda")
plt.title('Brain')

plt.subplot(2,3,2)
plt.plot(t_values, lam_l)
plt.ylabel("lambda(pi)")
plt.xlabel("lambda")
plt.title('Liver')

plt.subplot(2,3,3)
plt.plot(t_values, lam_m)
plt.ylabel("lambda(pi)")
plt.xlabel("lambda")
plt.title('Muscle')

plt.subplot(2,3,4)
plt.plot(t_values, lam_h)
plt.ylabel("lambda(pi)")
plt.xlabel("lambda")
plt.title('Heart')

plt.subplot(2,3,5)
plt.plot(t_values, lam_s)
plt.ylabel("lambda(pi)")
plt.xlabel("lambda")
plt.title('Serum')

#plt.savefig('/users/skean/Google Drive/MPI AGE/project/young_old_metabolome/analysis/MS1_analysis/lambda_pi.pdf')
plt.savefig('lambda_pi.pdf')


# In[21]:


#here i import more modules required for interpolation and plotting

import numpy as np
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from scipy import interpolate
import matplotlib.pyplot as plt


# In[22]:


#fit a spline y = spl(x) of degree k (default = cubic) to provided x,y data.
ius_b = scipy.interpolate.UnivariateSpline(t_values, lam_b)
ius_l = scipy.interpolate.UnivariateSpline(t_values, lam_l)
ius_m = scipy.interpolate.UnivariateSpline(t_values, lam_m)
ius_h = scipy.interpolate.UnivariateSpline(t_values, lam_h)
ius_s = scipy.interpolate.UnivariateSpline(t_values, lam_s)


# In[23]:


#set x to 1-d array of independent data
xs = t_values


# In[24]:


#here i plot the interpolated spline to the lambda values for each tissue

plt.subplot(2, 3, 1)
plt.plot(t_values, lam_b, 'r--')
plt.plot(xs, ius_b(xs), 'g')
plt.title('Interpolation using univariate spline - Brain')

plt.subplot(2, 3, 2)
plt.plot(t_values, lam_l, 'r--')
plt.plot(xs, ius_l(xs), 'g')
plt.title('Interpolation using univariate spline - Liver')

plt.subplot(2, 3, 3)
plt.plot(t_values, lam_m, 'r--')
plt.plot(xs, ius_m(xs), 'g')
plt.title('Interpolation using univariate spline - Muscle')

plt.subplot(2, 3, 4)
plt.plot(t_values, lam_h, 'r--')
plt.plot(xs, ius_h(xs), 'g')
plt.title('Interpolation using univariate spline - Heart')

plt.subplot(2, 3, 5)
plt.plot(t_values, lam_s, 'r--')
plt.plot(xs, ius_s(xs), 'g')
plt.title('Interpolation using univariate spline - Serum')

plt.savefig('cubic_splines.pdf')


# In[25]:


#calculate value of lambda(pi=1)
print ius_b(1)
print ius_l(1)
print ius_m(1)
print ius_h(1)
print ius_s(1)

#these values should be used for q-value calculations
#this is equal to the number that are NOT differentially expressed


# In[26]:


#truly alternative tests
print 1 - ius_b(1)
print 1 - ius_l(1)
print 1 - ius_m(1)
print 1 - ius_h(1)
print 1 - ius_s(1)


# In[27]:


#here i generate q-values for each tissue, based on the lambda values generated above

qvals_b = [ius_b(1) * float(len(pvals_b)) * i / (len([j for j in pvals_b if j <= i])) for i in pvals_b]
qvals_l = [ius_l(1) * float(len(pvals_l)) * i / (len([j for j in pvals_l if j <= i])) for i in pvals_l]
qvals_m = [ius_m(1) * float(len(pvals_m)) * i / (len([j for j in pvals_m if j <= i])) for i in pvals_m]
qvals_h = [ius_h(1) * float(len(pvals_h)) * i / (len([j for j in pvals_h if j <= i])) for i in pvals_h]
qvals_s = [ius_s(1) * float(len(pvals_s)) * i / (len([j for j in pvals_s if j <= i])) for i in pvals_s]


# In[28]:


#choosing a qvalue cut-off
number_ofsig_b = [len([i for i in qvals_b if i <= j]) for j in t_values]
number_ofsig_l = [len([i for i in qvals_l if i <= j]) for j in t_values]
number_ofsig_m = [len([i for i in qvals_m if i <= j]) for j in t_values]
number_ofsig_h = [len([i for i in qvals_h if i <= j]) for j in t_values]
number_ofsig_s = [len([i for i in qvals_s if i <= j]) for j in t_values]


# In[29]:


#plot q-value versus number of significant hits

plt.subplot(2,3,1)
plt.axis([0,0.4,0,3500])
matplotlib.rc('ytick', labelsize= 6) 
matplotlib.rc('xtick', labelsize= 5) 
plt.plot(t_values, number_ofsig_b)

plt.subplot(2,3,2)
plt.axis([0,0.3,0,3500])
matplotlib.rc('ytick', labelsize= 6) 
matplotlib.rc('xtick', labelsize= 5) 
plt.plot(t_values, number_ofsig_l)


plt.subplot(2,3,3)
plt.axis([0,0.5,0,3500])
matplotlib.rc('ytick', labelsize= 6) 
matplotlib.rc('xtick', labelsize= 5) 
plt.plot(t_values, number_ofsig_m)


plt.subplot(2,3,4)
plt.axis([0,0.4,0,3500])
matplotlib.rc('ytick', labelsize= 6) 
matplotlib.rc('xtick', labelsize= 5) 
plt.plot(t_values, number_ofsig_h)


plt.subplot(2,3,5)
plt.axis([0,0.3,0,3500])
matplotlib.rc('ytick', labelsize= 6) 
matplotlib.rc('xtick', labelsize= 5) 
plt.plot(t_values, number_ofsig_s)

plt.savefig('sig_with_qval.pdf')


# In[30]:


#this defines a new plotting function 'plotit2'. Creates barplots for individual metabolites that are significantly 
#changing with age, at a specified FDR value.

def plotit2(inp0, inp, fdr, tissue):
    q_fdr = BH(inp, fdr) #FDR of 0.1
    q_fdr_ind = [inp0.index(i) for i in inp0 if i < q_fdr]
    for k in range(len(q_fdr_ind)):
        data0 = a9a(met_all3, q_fdr_ind[k]+1) #the +1 is necessary, as the indexes of p_vals and met_all3 are shifted by 1
        if tissue == 'b':
            data = map(np.array, [data0.B6, data0.B10, data0.B16])
        elif tissue == 'l':
            data = map(np.array, [data0.L6, data0.L10, data0.L16])
        elif tissue == 'm':
            data = map(np.array, [data0.M6, data0.M10, data0.M16])
        elif tissue == 'h':
            data = map(np.array, [data0.H6, data0.H10, data0.H16])
        elif tissue == 's':
            data = map(np.array, [data0.S6, data0.S10, data0.S16])
        mass = met_all3[q_fdr_ind[k]+1][0] #this goes in the title
       
        
        plt.figure(figsize=(15,8))
        

        #plt.subplot(2,2,1)
        plt.boxplot(data) 
        plt.ylim(0,max([max(i)+0.1*max(i) for i in data ]))
        if tissue == 'b':
            plt.title('FDR = %s, BRAIN, %s '% (fdr, mass))
        elif tissue == 'l':
            plt.title('FDR = %s, LIVER, %s '% (fdr, mass))
        elif tissue == 'm':
            plt.title('FDR = %s, MUSCLE, %s '% (fdr, mass))
        elif tissue == 'h':
            plt.title('FDR = %s, HEART, %s '% (fdr, mass))
        elif tissue == 's':
            plt.title('FDR = %s, SERUM, %s '% (fdr, mass))
        plt.subplots_adjust(top=1.0, hspace=0.3, wspace=0.2)

        plt.xticks([1, 2, 3], ['6 weeks', '10 weeks', '16 weeks'])
        plt.ylabel('Abundance')

        for i in [1,2,3]:
            y = data[i-1]
            x = np.random.normal(i, 0.02, len(y))
            plt.plot(x, y, 'r.', alpha=0.5)


# In[31]:


#runs the new function plotit2 to generate barplots
#plotit2(pvals_b, pvals2_b, 0.01, 'b')


# In[32]:


#plotit2(pvals_l, pvals2_l, 0.05, 'l')


# In[33]:


#plotit2(pvals_m, pvals2_m, 0.1, 'm')


# In[34]:


#plotit2(pvals_h, pvals2_h, 0.05, 'h')


# In[35]:


#plotit2(pvals_s, pvals2_s, 0.1, 's')


# In[36]:


#define a new function 'print_it' with variables inp0, inp, fdr, tissue and outputfile name
#this generates many variables e.g 'mass' 'pvalb' which are then written to a defined output file

def print_it(inp0, inp, fdr, tissue, outputfile):
    #define the variable q_fdr as equal to function BH with variables inp, fdr
    #inp = input from user, i,e the pvals2 list
    #fdr = false discovery rate from user, e.g 0.1
    #BH: sorts the first variable 'inp' i.e pvals2 into a sorted list
    #BH: returns a q value for each line in sorted list
    output = open(outputfile, 'a')
    q_fdr = BH(inp, fdr) #FDR of 0.1
    #defines a new function q_fdr_ind as equal to the index i of inp0 for every value(i) in inp0, if i < q value in q_fdr
    q_fdr_ind = [inp0.index(i) for i in inp0 if i < q_fdr]
    #for values in the range length of q_fdr_ind
    for k in range(len(q_fdr_ind)):
       
        mass = met_all3[q_fdr_ind[k]+1][0] #this goes in the title
        pvalb = pvals_b[q_fdr_ind[k]]
        pvall = pvals_l[q_fdr_ind[k]]
        pvalm = pvals_m[q_fdr_ind[k]]
        pvalh = pvals_h[q_fdr_ind[k]]
        pvals = pvals_s[q_fdr_ind[k]]
        qvalb = qvals_b[q_fdr_ind[k]]
        qvall = qvals_l[q_fdr_ind[k]]
        qvalm = qvals_m[q_fdr_ind[k]]
        qvalh = qvals_h[q_fdr_ind[k]]
        qvals = qvals_s[q_fdr_ind[k]]
        fold6_10_b = log2fc6_10_b[q_fdr_ind[k]]
        fold10_16_b = log2fc10_16_b[q_fdr_ind[k]]
        fold6_16_b = log2fc6_16_b[q_fdr_ind[k]]
        fold6_10_l = log2fc6_10_l[q_fdr_ind[k]]
        fold10_16_l = log2fc10_16_l[q_fdr_ind[k]]
        fold6_16_l = log2fc6_16_l[q_fdr_ind[k]]
        fold6_10_m = log2fc6_10_m[q_fdr_ind[k]]
        fold10_16_m = log2fc10_16_m[q_fdr_ind[k]]
        fold6_16_m = log2fc6_16_m[q_fdr_ind[k]]
        fold6_10_h = log2fc6_10_h[q_fdr_ind[k]]
        fold10_16_h = log2fc10_16_h[q_fdr_ind[k]]
        fold6_16_h = log2fc6_16_h[q_fdr_ind[k]]
        fold6_10_s = log2fc6_10_s[q_fdr_ind[k]]
        fold10_16_s = log2fc10_16_s[q_fdr_ind[k]]
        fold6_16_s = log2fc6_16_s[q_fdr_ind[k]]
       
        if tissue == 'b':
            output.write('%s, BRAIN,%s, %s, %s, %s, %s, %s \n'% (fdr, mass, pvalb, qvalb, fold6_10_b, fold10_16_b, fold6_16_b))
        elif tissue == 'l':
            output.write('%s, LIVER, %s, %s, %s, %s, %s, %s \n'% (fdr, mass, pvall, qvall, fold6_10_l, fold10_16_l, fold6_16_l))  
        elif tissue == 'm':
            output.write('%s, MUSCLE, %s, %s, %s, %s, %s, %s \n'% (fdr, mass, pvalm, qvalm, fold6_10_m, fold10_16_m, fold6_16_m))  
        elif tissue == 'h':
            output.write('%s, HEART, %s, %s, %s, %s, %s, %s \n'% (fdr, mass, pvalh, qvalh, fold6_10_h, fold10_16_h, fold6_16_h))
        elif tissue == 's':
            output.write('%s, SERUM, %s, %s, %s, %s, %s, %s \n'% (fdr, mass, pvals, qvals, fold6_10_s, fold10_16_s, fold6_16_s))
            
    output.close()


# In[39]:


#here i am creating a new file and running the print_it function to generate information e.g pvalues, mass. These
#are then written to the new file
#!!! CHANGE OUTPUT FILE NAME BEFORE RUNNING AS TO NOT OVERWRITE PREVIOUS FILE!!!#
O = open("metabolites_all_qvalues_0.1fdr.csv", 'w')
O.write('fdr,tissue,mass,pval,qval,fc6to10,fc10to16,fc6to16\n')
O.close()

print_it(pvals_b, pvals2_b, 0.1, 'b', "metabolites_all_qvalues_0.1fdr.csv")
print_it(pvals_l, pvals2_l, 0.1, 'l', "metabolites_all_qvalues_0.1fdr.csv")
print_it(pvals_m, pvals2_m, 0.1, 'm', "metabolites_all_qvalues_0.1fdr.csv")
print_it(pvals_h, pvals2_h, 0.1, 'h', "metabolites_all_qvalues_0.1fdr.csv")
print_it(pvals_s, pvals2_s, 0.1, 's', "metabolites_all_qvalues_0.1fdr.csv")


# In[40]:


#same as cell above - creating new output files with different parameters
#!!! CHANGE OUTPUT FILE NAME BEFORE RUNNING AS TO NOT OVERWRITE PREVIOUS FILE!!!#
O = open("top_hit_list2.csv", 'w')
O.write('fdr,tissue,mass,pval,qval,fc6to10,fc10to16,fc6to16\n')
O.close()

print_it(pvals_b, pvals2_b, 0.1, 'b', "top_hit_list2.csv")
print_it(pvals_l, pvals2_l, 0.1, 'l', "top_hit_list2.csv")
print_it(pvals_m, pvals2_m, 0.1, 'm', "top_hit_list2.csv")
print_it(pvals_h, pvals2_h, 0.1, 'h', "top_hit_list2.csv")
print_it(pvals_s, pvals2_s, 0.1, 's', "top_hit_list2.csv")

#print_it(pvals_h, pvals2_h, 1, 'h', "top_hit_list2sphingheart.csv")
#print_it(pvals_s, pvals2_s, 1, 's', "top_hit_list2sphingserum.csv")

#print_it(pvals_b, pvals2_b, 0.1, 'b', "top_hit_list2_0.1.csv")
#print_it(pvals_l, pvals2_l, 0.1, 'l', "top_hit_list2_0.1.csv")
#print_it(pvals_m, pvals2_m, 0.1, 'm', "top_hit_list2_0.1.csv")
#print_it(pvals_h, pvals2_h, 0.1, 'h', "top_hit_list2_0.1.csv")
#print_it(pvals_s, pvals2_s, 0.1, 's', "top_hit_list2_0.1.csv")

#print_it(pvals_b, pvals2_b, 0.05, 'b', "top_hit_list2.csv")
#print_it(pvals_l, pvals2_l, 0.05, 'l', "top_hit_list2.csv")
#print_it(pvals_m, pvals2_m, 0.05, 'm', "top_hit_list2.csv")
#print_it(pvals_h, pvals2_h, 0.05, 'h', "top_hit_list2.csv")
#print_it(pvals_s, pvals2_s, 0.05, 's', "top_hit_list2.csv")

#print_it(pvals_b, pvals2_b, 0.01, 'b', "top_hit_list2.csv")
#print_it(pvals_l, pvals2_l, 0.01, 'l', "top_hit_list2.csv")
#print_it(pvals_m, pvals2_m, 0.01, 'm', "top_hit_list2.csv")
#print_it(pvals_h, pvals2_h, 0.01, 'h', "top_hit_list2.csv")
#print_it(pvals_s, pvals2_s, 0.01, 's', "top_hit_list2.csv")


# In[41]:


data = open('Metabolome_data_mz_and_abundance.csv').read()


# In[42]:


data = data.split('\n')


# In[43]:


data = [i.split(',') for i in data]


# In[44]:


mass_charges = [data[i][0] for i in range(len(data))]


# In[45]:


brain_data = [data[i][1:13] for i in range(len(data))]
liver_data = [data[i][13:25] for i in range(len(data))]
muscle_data = [data[i][25:37] for i in range(len(data))]
heart_data = [data[i][37:49] for i in range(len(data))]
serum_data = [data[i][49:65] for i in range(len(data))]


# In[46]:


brain_data = [[mass_charges[i]] + brain_data[i] for i in range(len(brain_data))]
liver_data = [[mass_charges[i]] + liver_data[i] for i in range(len(liver_data))]
muscle_data = [[mass_charges[i]] + muscle_data[i] for i in range(len(muscle_data))]
heart_data = [[mass_charges[i]] + heart_data[i] for i in range(len(heart_data))]
serum_data = [[mass_charges[i]] + serum_data[i] for i in range(len(serum_data))]


# In[114]:


import decimal


# In[214]:


output1 = open('brain_no_zero.csv', 'w+')
for i in range(1,(len(brain_data)-1)):
    y = brain_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = [str(j) for j in yint]
        output1.write(','.join(each_string))
        output1.write('\n')
        
output1.close()

output2 = open('liver_no_zero.csv', 'w+')
for i in range(1,(len(liver_data)-1)):
    y = liver_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = [str(j) for j in yint]
        output2.write(','.join(each_string))
        output2.write('\n')
        
output2.close()

output3 = open('muscle_no_zero.csv', 'w+')
for i in range(1,(len(muscle_data)-1)):
    y = muscle_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = [str(j) for j in yint]
        output3.write(','.join(each_string))
        output3.write('\n')
        
output3.close()

output4 = open('heart_no_zero.csv', 'w+')
for i in range(1,(len(heart_data)-1)):
    y = heart_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = [str(j) for j in yint]
        output4.write(','.join(each_string))
        output4.write('\n')
        
output4.close()

output5 = open('serum_no_zero.csv', 'w+')
for i in range(1,(len(serum_data)-1)):
    y = serum_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = [str(j) for j in yint]
        output5.write(','.join(each_string))
        output5.write('\n')
        
output5.close()
    


# In[215]:


#this gives ID column only

output1 = open('brain_index_nozero.csv', 'w+')
for i in range(1,(len(brain_data)-1)):
    y = brain_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = ['%r' %(yint[0])]
        output1.write(','.join(each_string))
        output1.write('\n')
        
output1.close()

output2 = open('liver_index_nozero.csv', 'w+')
for i in range(1,(len(liver_data)-1)):
    y = liver_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = ['%r' %(yint[0])]
        output2.write(','.join(each_string))
        output2.write('\n')
        
output2.close()

output3 = open('muscle_index_nozero.csv', 'w+')
for i in range(1,(len(muscle_data)-1)):
    y = muscle_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = ['%r' %(yint[0])]
        output3.write(','.join(each_string))
        output3.write('\n')
        
output3.close()

output4 = open('heart_index_nozero.csv', 'w+')
for i in range(1,(len(heart_data)-1)):
    y = heart_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = ['%r' %(yint[0])]
        output4.write(','.join(each_string))
        output4.write('\n')
        
output4.close()

output5 = open('serum_index_nozero.csv', 'w+')
for i in range(1,(len(serum_data)-1)):
    y = serum_data[i]
    yint = [float(k) for k in y]
    if sum(x>0 for x in yint[1:5]) >=3 and sum(x>0 for x in yint[5:9]) >=3 and sum(x>0 for x in yint[9:13]) >=3:
        each_string = ['%r' %(yint[0])]
        output5.write(','.join(each_string))
        output5.write('\n')
        
output5.close()


# In[216]:


#this opens the index files for each tissue, which are lists of metabolites that contain the minimum number of 
#read values (>= 3 out of 4 samples with read values).
#index files generated by script ''
brain_index = open('brain_index_nozero.csv').read()
liver_index = open('liver_index_nozero.csv').read()
muscle_index = open('muscle_index_nozero.csv').read()
heart_index = open('heart_index_nozero.csv').read()
serum_index = open('serum_index_nozero.csv').read()


# In[217]:


sig_data = open('top_hit_list2.csv').read()


# In[218]:


sig_data = sig_data.replace(' ', '')


# In[219]:


sig_data = sig_data.split('\n')


# In[220]:


sig_data = [i.split(',') for i in sig_data]


# In[221]:


sig_data = sig_data[:-1]


# In[222]:


brain_sig = [i for i in sig_data if i[1] == 'BRAIN']
liver_sig = [i for i in sig_data if i[1] == 'LIVER']
muscle_sig = [i for i in sig_data if i[1] == 'MUSCLE']
heart_sig = [i for i in sig_data if i[1] == 'HEART']
serum_sig = [i for i in sig_data if i[1] == 'SERUM']


# In[223]:


#formatting. make sure you are splitting by the correct characters e.g \r or \n. this may be different for each file 
#or different depending on operating system 
brain_index = brain_index.split('\n')
liver_index = liver_index.split('\n')
muscle_index = muscle_index.split('\n')
heart_index = heart_index.split('\n')
serum_index = serum_index.split('\n')


# In[236]:


#selects significant metabolites only if they are found in the index file (i.e those with enough read values)
#creates a new list of these items for each tissue
brain_new = []
for i in brain_sig:
    for j in brain_index:
        if i[2][:-1] in j:
            brain_new.append(i)
            
liver_new = []
for i in liver_sig:
    for j in liver_index:
        if i[2][:-1] in j:
            liver_new.append(i)
            
muscle_new = []
for i in muscle_sig:
    for j in muscle_index:
        if i[2][:-1] in j:
            muscle_new.append(i)
            
heart_new = []
for i in heart_sig:
    for j in heart_index:
        if i[2][:-1] in j:
            heart_new.append(i)
            
serum_new = []
for i in serum_sig:
    for j in serum_index:
        if i[2][:-1] in j:
            serum_new.append(i)


# In[241]:


len(serum_new)


# In[242]:


brain_new = ['\t'.join(i) for i in brain_new]
liver_new = ['\t'.join(i) for i in liver_new]
muscle_new = ['\t'.join(i) for i in muscle_new]
heart_new = ['\t'.join(i) for i in heart_new]
serum_new = ['\t'.join(i) for i in serum_new]


# In[243]:


brain_new = '\n'.join(brain_new)
liver_new = '\n'.join(liver_new)
muscle_new = '\n'.join(muscle_new)
heart_new = '\n'.join(heart_new)
serum_new = '\n'.join(serum_new)


# In[244]:


#writes data to files for each tissue
file = open('brain_top_hit_no_zero_values.txt', 'w')
file.write('fdr\ttissue\tmass\tpval\tqval\tfc6to10\tfc10to16\tfc6to16\n')
file.write(brain_new)
file.close()

file = open('liver_top_hit_no_zero_values.txt', 'w')
file.write('fdr\ttissue\tmass\tpval\tqval\tfc6to10\tfc10to16\tfc6to16\n')
file.write(liver_new)
file.close()

file = open('muscle_top_hit_no_zero_values.txt', 'w')
file.write('fdr\ttissue\tmass\tpval\tqval\tfc6to10\tfc10to16\tfc6to16\n')
file.write(muscle_new)
file.close()

file = open('heart_top_hit_no_zero_values.txt', 'w')
file.write('fdr\ttissue\tmass\tpval\tqval\tfc6to10\tfc10to16\tfc6to16\n')
file.write(heart_new)
file.close()

file = open('serum_top_hit_no_zero_values.txt', 'w')
file.write('fdr\ttissue\tmass\tpval\tqval\tfc6to10\tfc10to16\tfc6to16\n')
file.write(serum_new)
file.close()


# In[ ]:




