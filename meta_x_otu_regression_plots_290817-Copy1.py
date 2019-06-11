#!/usr/bin/env python
# coding: utf-8

# In[1]:


#importing modules required for regression script
import numpy as np
from scipy.stats.stats import pearsonr
from scipy import stats
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt


# In[2]:


#loading metabolome data
meta = open('metabolic_data_overlap_sorted.csv', 'rU').read()


# In[3]:


#loading OTU data
otu = open('otu_table_overlap_sorted.csv', 'rU').read()


# In[4]:


#formatting meta and otu data
metaL = [i.split(',') for i in meta.split('\n')] 
otuL = [i.split(',') for i in otu.split('\n')]

#below i transpose meta and otu data
#this creates a list of metabolite abundances. a list of lists
metat = [list(i) for i in zip(*metaL)] #transpose
#this creates a list of otu 'abundances'. a list of lists
otut = [list(i) for i in zip(*otuL)] #transpose


# In[5]:


#removing any meta or otu reads that have >=7 0 value reads
meta_value = [map(float, i[1:]) for i in metat[1:] if i.count('0.0') <7]
otu_value = [map(float, i[1:]) for i in otut[1:] if i.count('0.0') <7]


# In[8]:


#this creates a list of metabolite mass/charge 'identities'.
meta_key = [i[0] for i in metat[1:] if i.count('0.0') < 7]
#this creates a list of otu 'identities'.
otu_key = [i[0] for i in otut[1:] if i.count('0.0') < 7]

#this combines meta_key with metabolite abundance data (meta_value)
meta_dic = dict(zip(meta_key, meta_value))
#this combines otu_key with otu abundance data (otu_value)
otu_dic = dict(zip(otu_key, otu_value))


# In[13]:


#load table of linear regression data
correlshort = open('corr_short_zero.csv', 'rU').read()
#correlshort = fpkmvsotu


# In[14]:


class reg(object):
    
    """
    This class uses as input a linear regression file obtained by the loop_two function and returns 
    regression tables and plots, given specific parameteres, such as R-squared value, p-value, OTU or 
    transcript
    """
    
    version = 0.1
    
    def __init__(self, inp):
        self.inp = inp 
        self.inp2 = [ i.split('(') for i in self.inp.split('\n')[:-1]] 
    
        def loop1(inp, n):
            ls = []
            for i in inp[1:]:
                ls.append(i.split(',')[n])
            return ','.join(ls).replace(')', '')
        def loop2(inp, n):
            ls = []
            for i in inp:
                ls.append(loop1(i, n))
            return ls
       
        self.slope =  [map(float, i.split(',')) for i in loop2(self.inp2, 0)]
        self.intercept = [map(float, i.split(',')) for i in loop2(self.inp2, 1)]
        self.rvalue = [map(float, i.split(',')) for i in loop2(self.inp2, 2)]
        self.rsquare = [map(lambda x: x**2, i) for i in self.rvalue]
        self.pvalue = [map(float, i.split(',')) for i in loop2(self.inp2, 3)]
        self.stderr = [map(float, i.split(',')) for i in loop2(self.inp2, 4)  ]


# In[15]:


#needed for repr function to maintain decimal places
import decimal


# In[16]:


def filt_(inp, thr, sign): #filter p-values, r-values and so forth based on a chosen threshold
    ar = np.array(inp)
    li = []
    lp = []
    if sign == '>':
        lp = np.where(ar > thr)
    elif sign == '<':
        lp = np.where(ar < thr)
    else: 
        print "Error: check your input file"
    p1_l = [ list(i) for i in lp]
    la = list(ar[p1_l])
    p2 = p1_l+[la]
    p3 = [ list(i) for i in zip(*p2)]
    p4 = sorted(p3, key=lambda x: x[2])
    return ','.join([ ','.join(map(repr, i))+'\n' for i in p4]).replace('\n,','\n')[:-1]


# In[17]:


rf = reg(correlshort)


# In[18]:


Rpval = filt_(rf.pvalue, 0.05, '<')
P2_pval = [otu_key[int(i.split(',')[0])]+','+meta_key[int(i.split(',')[1])] +','+i.split(',')[2] for i in Rpval.split('\n')]
P3_pval = ','.join([ i+'\n' for i in P2_pval]).replace('\n,','\n')[:-1]
#len(P3_pval.split('\n'))
P3_pval = P3_pval.split('\n')


# In[19]:


rrsquare = filt_(rf.rsquare, 0.4, '>')
P2_rsquare = [otu_key[int(i.split(',')[0])]+','+meta_key[int(i.split(',')[1])] +','+i.split(',')[2] for i in rrsquare.split('\n') ]
P3_rsquare = ','.join([ i+'\n' for i in P2_rsquare]).replace('\n,','\n')[:-1]
len(P3_rsquare.split('\n'))
P3_rsquare = P3_rsquare.split('\n')
#need reverse for rsquare
P3_rsquare = P3_rsquare[::-1]


# In[20]:


X_values_pval = dict()
for i in range(len(P3_pval)):
    X_values_pval[i] = np.array(otu_dic[P3_pval[i].split(',')[0]])

Y_values_pval = dict()
for i in range(len(P3_pval)):
    Y_values_pval[i] = np.array(meta_dic[P3_pval[i].split(',')[1]])
    
X_values_rsquare = dict()
for i in range(len(P3_rsquare)):
    X_values_rsquare[i] = np.array(otu_dic[P3_rsquare[i].split(',')[0]])
                                   
Y_values_rsquare = dict()
for i in range(len(P3_rsquare)):
    Y_values_rsquare[i] = np.array(meta_dic[P3_rsquare[i].split(',')[1]])


# In[22]:


#PVALUE PLOTTING#
#sam - now repeat plotting code with iterations

# set x to number of regressions to plot
x = 4

for i in range(x):
    fig = plt.figure()
    fig.suptitle('#%s regression based on pval' % str(i+1), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(331)
    fig.subplots_adjust(top=0.85, hspace=0.8, wspace=0.8)
    #ax.set_title('axes title')
    ax.set_xlabel(str((P3_pval[i].split(',')[0]).split(';')[-1]))
    ax.set_ylabel(str([P3_pval[i].split(',')[1]]))
    plt.scatter(X_values_pval[i],Y_values_pval[i])
    fig.savefig('%s-top_regr_pval_290817.pdf' % i)


# In[23]:


#RSQUARE PLOTTING#
#sam - now repeat plotting code with iterations

# set x to number of regressions to plot
x = 5

for i in range(x):
    fig = plt.figure()
    fig.suptitle('#%s regression based on rsquare' % str(i+1), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(331)
    fig.subplots_adjust(top=0.85, hspace=0.8, wspace=0.8)
    #ax.set_title('axes title')
    ax.set_xlabel(str((P3_rsquare[i].split(',')[0]).split(';')[-1]))
    ax.set_ylabel(str([P3_rsquare[i].split(',')[1]]))
    plt.scatter(X_values_rsquare[i],Y_values_rsquare[i])
    fig.savefig('%s-top_regr_rsquare_290817.pdf' % i)


# In[24]:


P3_pval_split = [i.split(',') for i in P3_pval]
P3_rsquare_split = [i.split(',') for i in P3_rsquare]


# In[61]:


library = correlshort.replace(' ', '').split(')')
library2 = [i.split(',') for i in library]
library2 = [[i.replace('(', '') for i in j] for j in library2]
library2 = library2[:-1]


# In[64]:


#below code allows entry of a specific metabolite identifier, returning x number of top regressions based on pvalue or rsquare.
metabolite_of_interest = '149.0233'

rsquare_index_metabolite = [] 
for i in P3_rsquare_split:
    if i[1] == metabolite_of_interest:
        rsquare_index_metabolite.append(P3_rsquare_split.index(i)) 
        
x = rsquare_index_metabolite

#fig = plt.figure()
#fig.suptitle('top %s regressions based on pval' % x, fontsize=14, fontweight='bold')

for i in x:
    fig = plt.figure()
    fig.suptitle('rquare: #%s' % str(i+1), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(331)
    fig.subplots_adjust(top=0.85, hspace=0.8, wspace=0.8)
    #ax.set_title('axes title')
    plt.xticks(fontsize = 6)
    plt.yticks(fontsize = 6)
    ax.set_xlabel(str((P3_rsquare[i].split(',')[0]).split(';')[-1]))
    ax.set_ylabel(str([P3_rsquare[i].split(',')[1]]))
    #text = [j for j in library2 if (P3_rsquare[i].split(',')[2]) in j]
    #ax.text(20,0.4,'slope = %s\nintercept = %s\nrval = %s\nrsquare = %s\npval =%.15f\nstderr = %s'
            #%((text[0][-5]), (text[0][-4]) , (text[0][-3]), (float((text[0][-3])) **2),
              #(float(text[0][-2])), (text[0][-1])))
    plt.scatter(X_values_rsquare[i],Y_values_rsquare[i])
    fig.savefig('%s_%s290817.pdf' % (i, metabolite_of_interest))


# In[36]:


#sam - now repeat plotting code with iterations
metabolite_of_interest = '149.0233'

pval_index_metabolite = [] 
for i in P3_pval_split:
    if i[1] == metabolite_of_interest:
        pval_index_metabolite.append(P3_pval_split.index(i)) 

# set x to number of regressions to plot
x = pval_index_metabolite

#fig = plt.figure()
#fig.suptitle('top %s regressions based on pval' % x, fontsize=14, fontweight='bold')

for i in x:
    fig = plt.figure()
    fig.suptitle('pval: #%s' % str(i+1), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(331)
    fig.subplots_adjust(top=0.85, hspace=0.8, wspace=0.8)
    #ax.set_title('axes title')
    plt.xticks(fontsize = 6)
    plt.yticks(fontsize = 6)
    ax.set_xlabel(str((P3_pval[i].split(',')[0]).split(';')[-1]))
    ax.set_ylabel(str([P3_pval[i].split(',')[1]]))
    #text = [j for j in library2 if (P3[i].split(',')[2]) in j]
    #ax.text(20,0.4,'slope = %s\nintercept = %s\nrval = %s\nrsquare = %s\npval =%.15f\nstderr = %s'
            #%((text[0][-5]), (text[0][-4]) , (text[0][-3]), (float((text[0][-3])) **2),
              #(float(text[0][-2])), (text[0][-1])))
    plt.scatter(X_values_pval[i],Y_values_pval[i])
    fig.savefig('%s_%s_290817.pdf' % (i, metabolite_of_interest))


# In[32]:


#extract metabolite values from top regressions
#metabolites_by_regression = [str(P3_rsquare[i].split(',')[1]) for i in range(len(P3_rsquare))]


# In[33]:


#z = open('/Users/skean/Google Drive/MPI AGE/project/x-analysis/metabolites_by_regression.csv', 'w')
#z.write(str((',').join(metabolites_by_regression)))
#z.close()


# In[42]:





# In[43]:


#sam - now repeat plotting code with iterations

# set x to number of regressions to plot
x = pval_index_1490233

#fig = plt.figure()
#fig.suptitle('top %s regressions based on pval' % x, fontsize=14, fontweight='bold')

for i in x:
    fig = plt.figure()
    fig.suptitle('by pval: #%s' % str(i+1), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(331)
    fig.subplots_adjust(top=0.85, hspace=0.8, wspace=0.8)
    #ax.set_title('axes title')
    plt.xticks(fontsize = 6)
    plt.yticks(fontsize = 6)
    ax.set_xlabel(str((P3_pval[i].split(',')[0]).split(';')[-1]))
    ax.set_ylabel(str([P3_pval[i].split(',')[1]]))
    text = [j for j in library2 if (P3_pval[i].split(',')[2]) in j]
    ax.text(20,0.4,'slope = %s\nintercept = %s\nrval = %s\nrsquare = %s\npval =%.15f\nstderr = %s'
            %((text[0][-5]), (text[0][-4]) , (text[0][-3]), (float((text[0][-3])) **2),
              (float(text[0][-2])), (text[0][-1])))
    plt.scatter(X_values_pval[i],Y_values_pval[i])
    fig.savefig('/Users/skean/Google Drive/MPI AGE/project/x-analysis/%s_oxo-meth_regr_pval_withstats_220517.pdf' % i)


# In[267]:


#sam - now repeat plotting code with iterations

# set x to number of regressions to plot
x = rsquare_index_1490233

#fig = plt.figure()
#fig.suptitle('top %s regressions based on pval' % x, fontsize=14, fontweight='bold')

for i in x:
    fig = plt.figure()
    fig.suptitle('rval: #%s' % str(i+1), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(331)
    fig.subplots_adjust(top=0.85, hspace=0.8, wspace=0.8)
    #ax.set_title('axes title')
    plt.xticks(fontsize = 6)
    plt.yticks(fontsize = 6)
    ax.set_xlabel(str((P3_rsquare[i].split(',')[0]).split(';')[-1]))
    ax.set_ylabel(str([P3_rsquare[i].split(',')[1]]))
    text = [j for j in library2 if (P3_pval[i].split(',')[2]) in j]
    ax.text(20,0.4,'slope = %s\nintercept = %s\nrval = %s\nrsquare = %s\npval =%.15f\nstderr = %s'
            %((text[0][-5]), (text[0][-4]) , (text[0][-3]), (float((text[0][-3])) **2),
              (float(text[0][-2])), (text[0][-1])))
    plt.scatter(X_values_rsquare[i],Y_values_rsquare[i])
    fig.savefig('/Users/skean/Google Drive/MPI AGE/project/x-analysis/%s_oxo-meth_regr_rsquare_withstats_220517.pdf' % i)


# In[ ]:




