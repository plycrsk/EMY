#!/usr/bin/env python
# coding: utf-8

# In[21]:


#this script takes the output from 'metabolite_script_annotated_180817.ipynb' and formats it as required by PIUMet
#this opens the output file and sets it to the variable 'top_hits'
brain_top_hits = open('brain_top_hit_no_zero_values_280817.txt').read()
liver_top_hits = open('liver_top_hit_no_zero_values_280817.txt').read()
muscle_top_hits = open('muscle_top_hit_no_zero_values_280817.txt').read()
heart_top_hits = open('heart_top_hit_no_zero_values_280817.txt').read()
serum_top_hits = open('serum_top_hit_no_zero_values_280817.txt').read()


# In[34]:


#this opens the Mass Spec raw_data file and sets it to the variable 'raw_data'
raw_data = open('RAW_DATA111116.csv').read()


# In[24]:


#the below cells are formatting the top_hits data
brain_top_hits = brain_top_hits.replace(' ', '')
liver_top_hits = liver_top_hits.replace(' ', '')
muscle_top_hits = muscle_top_hits.replace(' ', '')
heart_top_hits = heart_top_hits.replace(' ', '')
serum_top_hits = serum_top_hits.replace(' ', '')


# In[26]:


brain_top_hits = brain_top_hits.split('\n')
liver_top_hits = liver_top_hits.split('\n')
muscle_top_hits = muscle_top_hits.split('\n')
heart_top_hits = heart_top_hits.split('\n')
serum_top_hits = serum_top_hits.split('\n')


# In[27]:


brain_top_hits = [i.split('\t') for i in brain_top_hits]
liver_top_hits = [i.split('\t') for i in liver_top_hits]
muscle_top_hits = [i.split('\t') for i in muscle_top_hits]
heart_top_hits = [i.split('\t') for i in heart_top_hits]
serum_top_hits = [i.split('\t') for i in serum_top_hits]


# In[29]:


brain_top_hits = [i[2] for i in brain_top_hits]
liver_top_hits = [i[2] for i in liver_top_hits]
muscle_top_hits = [i[2] for i in muscle_top_hits]
heart_top_hits = [i[2] for i in heart_top_hits]
serum_top_hits = [i[2] for i in serum_top_hits]


# In[35]:


raw_data = [i.split(',') for i in raw_data.split('\n')][1:-1]


# In[36]:


raw_data = [i[1:4] for i in raw_data]


# In[37]:


#this removes all blank data from the raw_data variable. 
#this code could be rewritten to remove blanks rather than removing those with a specific number
raw_data = raw_data[:-82]


# In[38]:


#should be equal to 3595 (i.e number of metabolites measured)
len(raw_data)


# In[52]:


#generates a new list 'top_hit_for_piumet'
brain_top_hit_for_piumet = []
liver_top_hit_for_piumet = []
muscle_top_hit_for_piumet = []
heart_top_hit_for_piumet = []
serum_top_hit_for_piumet = []

#searches raw_data for all items in top_hits. appends the item from top hits to matching item from raw_data
for i in brain_top_hits:
    for j in raw_data:
        if i == j[1]:
            brain_top_hit_for_piumet.append(j)
            
for i in liver_top_hits:
    for j in raw_data:
        if i == j[1]:
            liver_top_hit_for_piumet.append(j)
            
for i in muscle_top_hits:
    for j in raw_data:
        if i == j[1]:
            muscle_top_hit_for_piumet.append(j)
            
for i in heart_top_hits:
    for j in raw_data:
        if i == j[1]:
            heart_top_hit_for_piumet.append(j)
            
for i in serum_top_hits:
    for j in raw_data:
        if i == j[1]:
            serum_top_hit_for_piumet.append(j)


# In[60]:


#formatting
brain = [brain_top_hit_for_piumet[i][0:3] for i in range(len(brain_top_hit_for_piumet))]
liver = [liver_top_hit_for_piumet[i][0:3] for i in range(len(liver_top_hit_for_piumet))]
muscle = [muscle_top_hit_for_piumet[i][0:3] for i in range(len(muscle_top_hit_for_piumet))]
heart = [heart_top_hit_for_piumet[i][0:3] for i in range(len(heart_top_hit_for_piumet))]
serum = [serum_top_hit_for_piumet[i][0:3] for i in range(len(serum_top_hit_for_piumet))]


# In[61]:


#changing text 'Pos' to 'positive' and text 'Neg' to 'negative. This is as required by PIUMet
for i in brain:
    if i[0] == 'Pos':
        i[0] = 'positive'
    elif i[0] == 'Neg':
        i[0] = 'negative'
        
for i in liver:
    if i[0] == 'Pos':
        i[0] = 'positive'
    elif i[0] == 'Neg':
        i[0] = 'negative'

for i in muscle:
    if i[0] == 'Pos':
        i[0] = 'positive'
    elif i[0] == 'Neg':
        i[0] = 'negative'
        
for i in heart:
    if i[0] == 'Pos':
        i[0] = 'positive'
    elif i[0] == 'Neg':
        i[0] = 'negative'
        
for i in serum:
    if i[0] == 'Pos':
        i[0] = 'positive'
    elif i[0] == 'Neg':
        i[0] = 'negative'


# In[62]:


#defines a new list 'myorder' which is used for re-ordering the data
myorder = [1,0,2]


# In[63]:


#creates new variables for each tissue, and uses 'myorder' to re-order the data
brain = [[brain[j][i] for i in myorder] for j in range(len(brain))]
liver = [[liver[j][i] for i in myorder] for j in range(len(liver))]
muscle = [[muscle[j][i] for i in myorder] for j in range(len(muscle))]
heart = [[heart[j][i] for i in myorder] for j in range(len(heart))]
serum = [[serum[j][i] for i in myorder] for j in range(len(serum))]


# In[64]:


#create output files and reformat data ready for PIUMet

file = open('brain_top_hit_for_piumet_removed_zero_values_reorganised_280817.txt', 'w')
file.write('\r'.join(['\t'.join(i) for i in brain]))
file.close()
file = open('liver_top_hit_for_piumet_removed_zero_values_reorganised_280817.txt', 'w')
file.write('\r'.join(['\t'.join(i) for i in liver]))
file.close()
file = open('muscle_top_hit_for_piumet_removed_zero_values_reorganised_280817.txt', 'w')
file.write('\r'.join(['\t'.join(i) for i in muscle]))
file.close()
file = open('heart_top_hit_for_piumet_removed_zero_values_reorganised_280817.txt', 'w')
file.write('\r'.join(['\t'.join(i) for i in heart]))
file.close()
file = open('serum_top_hit_for_piumet_removed_zero_values_reorganised_280817.txt', 'w')
file.write('\r'.join(['\t'.join(i) for i in serum]))
file.close()

