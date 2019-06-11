#!/usr/bin/env python
# coding: utf-8

# In[2]:


bac = open('Description_otu_table_sorted_L6n.csv', 'rU').read()


# In[3]:


bac1 = bac.replace('_','.').replace('YI','WT') #replaces _ for . and 'YI' for 'WT'
bac2 = bac1.split('\n')[1:] # splits by new lines
bac3 = [i.split(',') for i in bac2] #splits by comma
bac4 = [bac3[0]]+[i for i in bac3[1:] if sum(map(float, (i[1:]))) != 0.0] #this removes transcripts with only 0s. 80 rows were removed, checked by len(c3) and len(c4)
#c3[0] is the first row and contains all the individuals.
#The map(aFunction, aSequence) function applies a passed-in function to each item in an iterable object and returns a list containing all the function call results.
#float makes a number or a string into a floating number.
bact = zip(*bac4) #flip the matrix / transpose
bacT = [list(i) for i in bact] # makes parentheses in square brackets, meaning making the tuples into a list. tuples are not mutable!


# In[4]:


meta_try = open('serum_Nov2015_named_correct.txt', 'rU').read()
meta_try1 = meta_try.replace(';',',').replace('\t',',')
meta2 = meta_try1.split('\n')
meta_cut = meta2[1:]
meta3 = [i.split(',') for i in meta_cut] #meta3 correspons to cT from level
#not necessary to remove transcripts with only 0's as there are none
metat =zip(*meta3)
metaT = [list(i) for i in metat]


# In[5]:


from sets import Set
bacHead = [i[0] for i in bacT[1:]]
metaHead = [i[0] for i in metaT[1:]]
Head_all = list(Set(bacHead) & Set(metaHead)) #check whether the same, take overlapping ones (RNA-seq way more individuals)
bac_overlap = [bacT[0]] + [i for i in bacT if i[0] in Head_all]
meta_overlap = [metaT[0]] + [i for i in metaT if i[0] in Head_all]


# In[6]:


bac_sorted = [bac_overlap[0]] + sorted(bac_overlap[1:], key = lambda x: (x[0]))
meta_sorted = [meta_overlap[0]] + sorted(meta_overlap[1:], key = lambda x: (x[0]))


# In[7]:


[i[0] for i in bac_sorted][1:] == [i[0] for i in meta_sorted][1:] #to check whether header in same order


# In[8]:


bac_joined = ','.join([','.join(i)+'\n' for i in bac_sorted]).replace('\n,','\n')[:-1]
meta_joined = ','.join([','.join(i)+'\n' for i in meta_sorted]).replace('\n,','\n')[:-1]


# In[9]:


zu = open('otu_table_overlap_sorted.csv', 'w')
zu.write(bac_joined)
zu.close
#entspricht uj and cj or rna/otu in darios script

zt = open('metabolic_data_overlap_sorted.csv', 'w')
zt.write(meta_joined)
zt.close()


# In[10]:


[i.split(',')[0] for i in bac_sorted.split('\n')[1:-1]] == [i.split(',')[0] for i in meta_sorted.split('\n')[1:-1]]


# In[11]:


import numpy as np
from scipy.stats.stats import pearsonr
from scipy import stats
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt


# In[12]:


#metaL == meta_sorted


# In[13]:


#metaL = [i.split(',') for i in meta_joined.split('\n')] #actually is the same as meta_sorted
#otuL = [i.split(',') for i in bac_joined.split('\n')]

meta_sorted_transposed = [list(i) for i in zip(*meta_sorted)] #makes things into list
bac_sorted_transposed = [list(i) for i in zip(*bac_sorted)] #makes things into list

meta_value = [map(float, i[1:]) for i in meta_sorted_transposed[1:] if i.count('0.0') <7]
bac_value = [map(float, i[1:]) for i in bac_sorted_transposed[1:] if i.count('0.0') <7]

meta_key = [i[0] for i in meta_sorted_transposed[1:] if i.count('0.0') < 7]
bac_key = [ i[0] for i in bac_sorted_transposed[1:] if i.count('0.0') < 7]

meta_dic = dict(zip(meta_key, meta_value))
bac_dic = dict(zip(bac_key, bac_value))

#before: '0' instead of '0.0'


# In[14]:


len(meta_value)


# In[15]:


z = open('meta_value_zero.txt','w')
z.write(meta_value)
z.close()

zz = open('bac_value_zero.txt','w')
zz.write(bac_value)
zz.close()


# In[17]:


def loop_one(m1,otu):
    ls = [] #creates empty list
    for i in m1:
        ls.append(stats.linregress(i, otu)) #Calculate a linear least-squares regression for two sets of measurements (here i and otu). every line (i) from m1 will be compared to whole otu
    return ','.join(map(str, ls))+'\n' #joins the outcome of ls (as a string) with a new line, separated by comma 

def loop_two(m1,m2):
    ls = []
    for j in m2:
        ls.append(loop_one(m1, j)) #takes every line from m2(j) to loop 1 = one line from rna (i) to one line from otu (j)
    return ','.join(ls).replace('\n,','\n')


# In[31]:


#myList = ','.join(map(str, myList)) if you want strings from a list but list contains also numbers etc


# In[41]:


meta_vs_bac = loop_two(meta_value, bac_value) #compares every line in OTU to every line in RNA set (one line from rna (i) to one line from otu (j))                               


# In[18]:


meta_vs_bac_0 = loop_two(meta_value, bac_value) #compares every line in OTU to every line in RNA set (one line from rna (i) to one line from otu (j))                               


# In[43]:


#meta_vs_bac2 = meta_vs_bac_loaded.replace('LinregressResult', '').replace('slope=','').replace('intercept=', '').replace('rvalue=', '').replace('pvalue=', '').replace('stderr=', '') 


# In[19]:


meta_vs_bac_0_2 =meta_vs_bac_0.replace('LinregressResult', '').replace('slope=','').replace('intercept=', '').replace('rvalue=', '').replace('pvalue=', '').replace('stderr=', '') 


# In[20]:


z = open('corr_zero_300817.csv', 'w')
z.write(meta_vs_bac_0)
z.close()


# In[22]:


z = open('corr_short_zero_300817.csv', 'w')
z.write(meta_vs_bac_0_2)
z.close()


# In[42]:


meta_vs_bac_loaded = open('/beegfs/group_dv/home/MPopkes/Redo_Dario/Files/corr.csv', 'rU').read()


# In[46]:


meta_split = meta_vs_bac.split('LinregressResult')
meta_replace = [i.replace('slope=','') for i in meta_split]
meta_replace = [i.replace('intercept=', '') for i in meta_replace]
meta_replace = [i.replace('rvalue=', '') for i in meta_replace]
meta_replace = [i.replace('pvalue=', '') for i in meta_replace]
meta_replace = [i.replace('stderr=', '') for i in meta_replace]   


# In[49]:


meta_vs_bac


# In[2]:


meta_key[1:10000]


# In[51]:


len(meta_replace)


# In[1]:


regress = open('/beegfs/group_dv/home/MPopkes/Redo_Dario/Files/corr.csv', 'rU').read()


# In[2]:


regress_split = regress.split('LinregressResult')
regress_replace = [i.replace('slope=','') for i in regress_split]
regress_replace = [i.replace('intercept=', '') for i in regress_replace]
regress_replace = [i.replace('rvalue=', '') for i in regress_replace]
regress_replace = [i.replace('pvalue=', '') for i in regress_replace]
regress_replace = [i.replace('stderr=', '') for i in regress_replace]


# In[3]:


regress[1:100000]


# In[4]:


regress_replace[1:100000]


# In[7]:


regress_replace_str = str(regress_replace)


# In[8]:


zz = open('/beegfs/group_dv/home/MPopkes/Redo_Dario/Files/correlshort.csv', 'w')
zz.write(regress_replace_str)
zz.close()


# In[9]:


regress_replace_str[1:10000]


# In[18]:


bac_list = []
for i in bac_key:
    x = np.repeat(i, len(meta_key))
    bac_list.append(x)


# In[15]:


meta_list = np.repeat(meta_key, len(bac_key))


# In[52]:


#bac3 = [i.split(',') for i in bac2] #splits by comma
#bac_list_2 = np.array(bac_list)
for i in bac_list:
    for x in bac_list[i]:
        new = x.split(',')
        
#bac_list_split = [i.split(',') for i in bac_list]


# In[14]:


len(meta_list)


# In[51]:


bac_list[1]


# In[42]:


len(regress_replace)

