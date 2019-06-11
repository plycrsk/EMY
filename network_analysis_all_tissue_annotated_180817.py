#!/usr/bin/env python
# coding: utf-8

# In[3]:


#this module is for my networking functions
import networkx as nx
#this module is for dataframe creation
import pandas as pd
#these modules are for plotting
import matplotlib as mpb
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
#another plotting in ipylab module
from pylab import show
import pylab


# In[4]:


#load and format data
brain_network_data = open("new_edge_file_brain_noPE_020717.txt").read()
b2 = brain_network_data.split("\r")
union_list_brain = []
for i in b2:
    a = i.split('\t')
    union_list_brain.append(a)
#create data frame from data
bdf = pd.DataFrame(union_list_brain)

liver_network_data = open("new_edge_file_liver_020717_noPE.txt").read()
l2 = liver_network_data.split("\r")
union_list_liver = []
for i in l2:
    a = i.split('\t')
    union_list_liver.append(a)
#create data frame from data
ldf = pd.DataFrame(union_list_liver)

muscle_network_data = open("new_edge_file_muscle_noPE_020717.txt").read()
m2 = muscle_network_data.split("\r")
union_list_muscle = []
for i in m2:
    a = i.split('\t')
    union_list_muscle.append(a)
#create data frame from data
mdf = pd.DataFrame(union_list_muscle)

heart_network_data = open("new_edge_file_heart_110517.txt").read()
h2 = heart_network_data.split("\n")
union_list_heart = []
for i in h2:
    a = i.split('\t')
    union_list_heart.append(a)
#create data frame from data
hdf = pd.DataFrame(union_list_heart)

serum_network_data = open("new_edge_file_serum_noPE_020717.txt").read()
s2 = serum_network_data.split("\r")
union_list_serum = []
for i in s2:
    a = i.split('\t')
    union_list_serum.append(a)
#create data frame from data
sdf = pd.DataFrame(union_list_serum)


# In[5]:


union_list_brain = union_list_brain[:-1]
union_list_liver = union_list_liver[:-1]
union_list_muscle = union_list_muscle[:-1]
union_list_heart = union_list_heart[:-1]
union_list_serum = union_list_serum[:-1]


# In[6]:


union_list_brain


# In[7]:


#create an empty graph
Gb = nx.Graph()
#add edges from data frame (for each element of data frame i)
for i in union_list_brain:
    Gb.add_edge(i[0],i[1], weight = i[2])
    
#create an empty graph
Gl = nx.Graph()
#add edges from data frame (for each element of data frame i)
for i in union_list_liver:
    Gl.add_edge(i[0],i[1], weight = i[2])
    
#create an empty graph
Gm = nx.Graph()
#add edges from data frame (for each element of data frame i)
for i in union_list_muscle:
    Gm.add_edge(i[0],i[1], weight = i[2])
    
#create an empty graph
Gh = nx.Graph()
#add edges from data frame (for each element of data frame i)
for i in union_list_heart:
    Gh.add_edge(i[0],i[1], weight = i[2])
    
#create an empty graph
Gs = nx.Graph()
#add edges from data frame (for each element of data frame i)
for i in union_list_serum:
    Gs.add_edge(i[0],i[1], weight = i[2])


# In[70]:


#generate coords to use for all networks using spring layout
b_coords = nx.spring_layout(Gb)
l_coords = nx.spring_layout(Gl)
m_coords = nx.spring_layout(Gm)
h_coords = nx.spring_layout(Gh)
s_coords = nx.spring_layout(Gs)


# In[71]:


#generate coords to use for all networks using spectral layout
b2_coords = nx.spectral_layout(Gb)
l2_coords = nx.spectral_layout(Gl)
m2_coords = nx.spectral_layout(Gm)
h2_coords = nx.spectral_layout(Gh)
s2_coords = nx.spectral_layout(Gs)


# In[72]:


#generate coords to use for all networks using fruchterman reingold layout
b3_coords = nx.fruchterman_reingold_layout(Gb)
l3_coords = nx.fruchterman_reingold_layout(Gl)
m3_coords = nx.fruchterman_reingold_layout(Gm)
h3_coords = nx.fruchterman_reingold_layout(Gh)
s3_coords = nx.fruchterman_reingold_layout(Gs)


# In[73]:


#create a dictionary with nodes and their degree
b_degrees = nx.degree(Gb).items()
l_degrees = nx.degree(Gl).items()
m_degrees = nx.degree(Gm).items()
h_degrees = nx.degree(Gh).items()
s_degrees = nx.degree(Gs).items()


# In[74]:


#this was added 03.07.17 - can print list of top hubs
brain_hubs = sorted(b_degrees, key = lambda x: x[1])
liver_hubs = sorted(l_degrees, key = lambda x: x[1])
muscle_hubs = sorted(m_degrees, key = lambda x: x[1])
heart_hubs = sorted(h_degrees, key = lambda x: x[1])
serum_hubs = sorted(s_degrees, key = lambda x: x[1])


# In[97]:


brain_hubs = [list(i) for i in brain_hubs]
liver_hubs = [list(i) for i in liver_hubs]
muscle_hubs = [list(i) for i in muscle_hubs]
heart_hubs = [list(i) for i in heart_hubs]
serum_hubs = [list(i) for i in serum_hubs]


# In[111]:





# In[112]:


#save degree information as txt file
import csv
with open('brain_degrees_030717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in brain_hubs:
        wr.writerow([i][0])
        
with open('liver_degrees_030717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in liver_hubs:
        wr.writerow([i][0])
        
with open('muscle_degrees_030717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in muscle_hubs:
        wr.writerow([i][0])
        
with open('heart_degrees_030717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in heart_hubs:
        wr.writerow([i][0])
        
with open('serum_degrees_030717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in serum_hubs:
        wr.writerow([i][0])


# In[75]:


#plot network with labels for nodes with degree => 5
Gb_degree = Gb.copy()
def b_degree_threshold(Gb):
    for k, v in b_degrees:
        if v < 4:
            Gb_degree.remove_node(k)
    return Gb_degree

Gb_degree = b_degree_threshold(Gb)

#plot network with labels for nodes with degree => 5
Gl_degree = Gl.copy()
def l_degree_threshold(Gl):
    for k, v in l_degrees:
        if v < 4:
            Gl_degree.remove_node(k)
    return Gl_degree

Gl_degree = l_degree_threshold(Gl)

#plot network with labels for nodes with degree => 5
Gm_degree = Gm.copy()
def m_degree_threshold(Gm):
    for k, v in m_degrees:
        if v < 4 :
            Gm_degree.remove_node(k)
    return Gm_degree

Gm_degree = m_degree_threshold(Gm)

#plot network with labels for nodes with degree => 5
Gh_degree = Gh.copy()
def h_degree_threshold(Gh):
    for k, v in h_degrees:
        if v < 2:
            Gh_degree.remove_node(k)
    return Gh_degree

Gh_degree = h_degree_threshold(Gh)


#plot network with labels for nodes with degree => 5
Gs_degree = Gs.copy()
def s_degree_threshold(Gs):
    for k, v in s_degrees:
        if v < 4:
            Gs_degree.remove_node(k)
    return Gs_degree

Gs_degree = s_degree_threshold(Gs)


# In[114]:


nx.draw_graphviz(Gb)


# In[76]:


###BRAIN DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gb,b_coords,node_color='blue',alpha=1,node_size=1)
nx.draw_networkx_edges(Gb,b_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gb_degree,b_coords,node_color='red',alpha=0.3,node_size=64)
# also the labels this time
nx.draw_networkx_labels(Gb_degree,b_coords,font_size=4,font_color='black')
plt.savefig("brain_FR_030717.png", dpi=1000)
plt.show()


# In[36]:


###LIVER DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gl,l_coords,node_color='blue',alpha=1,node_size=1)
nx.draw_networkx_edges(Gl,l_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gl_degree,l_coords,node_color='red',alpha=0.3,node_size=64)
# also the labels this time
nx.draw_networkx_labels(Gl_degree,l_coords,font_size=4,font_color='black')
plt.savefig("liver_030717.png", dpi=1000)
plt.show()


# In[37]:


###MUSCLE DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gm,m_coords,node_color='blue',alpha=1,node_size=1)
nx.draw_networkx_edges(Gm,m_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gm_degree,m_coords,node_color='red',alpha=0.3,node_size=64)
# also the labels this time
nx.draw_networkx_labels(Gm_degree,m_coords,font_size=4,font_color='black')
plt.savefig("muscle_030717.png", dpi=1000)
plt.show()


# In[39]:


###HEART DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gh,h_coords,node_color='blue',alpha=1,node_size=1)
nx.draw_networkx_edges(Gh,h_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gh_degree,h_coords,node_color='red',alpha=0.3,node_size=64)
# also the labels this time
nx.draw_networkx_labels(Gh_degree,h_coords,font_size=4,font_color='black')
plt.savefig("heart_030717.png", dpi=1000)
plt.show()


# In[41]:


###SERUM DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gs,s_coords,node_color='blue',alpha=1,node_size=1)
nx.draw_networkx_edges(Gs,s_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gs_degree,s_coords,node_color='red',alpha=0.3,node_size=64)
# also the labels this time
nx.draw_networkx_labels(Gs_degree,s_coords,font_size=4,font_color='black')
plt.savefig("serum_030717.png", dpi=1000)
plt.show()


# In[51]:


#report lists of degree >=5 nodes for each tissue
import csv
Brain_degree_5 = Gb_degree.nodes()
with open('brain_degree_5+_noPE_020717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Brain_degree_5:
        wr.writerow([i])
          
Liver_degree_5 = Gl_degree.nodes()
with open('liver_degree_5+_noPE_020717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Liver_degree_5:
        wr.writerow([i])
    
Muscle_degree_5 = Gm_degree.nodes()
with open('muscle_degree_5+_noPE_020717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Muscle_degree_5:
        wr.writerow([i])

    
Serum_degree_5 = Gs_degree.nodes()
with open('serum_degree_5+_noPE_020717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Serum_degree_5:
        wr.writerow([i])


# In[52]:


Muscle_degree_5


# In[53]:


set(Brain_degree_5) & set(Liver_degree_5) & set(Muscle_degree_5)


# In[54]:


Liver_degree_5


# In[55]:


#plot network with labels for nodes with degree => 4
Gb_degree = Gb.copy()
def b_degree_threshold(Gb):
    for k, v in b_degrees:
        if v < 4:
            Gb_degree.remove_node(k)
    return Gb_degree

Gb_degree = b_degree_threshold(Gb)

#plot network with labels for nodes with degree => 4
Gl_degree = Gl.copy()
def l_degree_threshold(Gl):
    for k, v in l_degrees:
        if v < 4:
            Gl_degree.remove_node(k)
    return Gl_degree

Gl_degree = l_degree_threshold(Gl)

#plot network with labels for nodes with degree => 4
Gm_degree = Gm.copy()
def m_degree_threshold(Gm):
    for k, v in m_degrees:
        if v < 4:
            Gm_degree.remove_node(k)
    return Gm_degree

Gm_degree = m_degree_threshold(Gm)


#plot network with labels for nodes with degree => 4
Gs_degree = Gs.copy()
def s_degree_threshold(Gs):
    for k, v in s_degrees:
        if v < 4:
            Gs_degree.remove_node(k)
    return Gs_degree

Gs_degree = s_degree_threshold(Gs)


# In[56]:


###BRAIN DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gb,b_coords,node_color='blue',alpha=0.2,node_size=8)
nx.draw_networkx_edges(Gb,b_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gb_degree,b_coords,node_color='red',alpha=0.5,node_size=128)
# also the labels this time
nx.draw_networkx_labels(Gb_degree,b_coords,font_size=8,font_color='blue')
plt.savefig("brain_degree4+.png", dpi=1000)
plt.show()


# In[57]:


###LIVER DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gl,l_coords,node_color='blue',alpha=0.2,node_size=8)
nx.draw_networkx_edges(Gl,l_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gl_degree,l_coords,node_color='red',alpha=0.5,node_size=254)
# also the labels this time
nx.draw_networkx_labels(Gl_degree,l_coords,font_size=8,font_color='blue')
plt.savefig("liver_degree4+.png", dpi=1000)
plt.show()


# In[58]:


###MUSCLE DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gm,m_coords,node_color='blue',alpha=0.2,node_size=8)
nx.draw_networkx_edges(Gm,m_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gm_degree,m_coords,node_color='red',alpha=0.5,node_size=254)
# also the labels this time
nx.draw_networkx_labels(Gm_degree,m_coords,font_size=8,font_color='blue')
plt.savefig("muscle_degree4+.png", dpi=1000)
plt.show()


# In[59]:


###SERUM DEGREE >=5####
# draw the nodes and the edges (all)
nx.draw_networkx_nodes(Gs,s_coords,node_color='blue',alpha=0.2,node_size=8)
nx.draw_networkx_edges(Gs,s_coords,alpha=0.1)

# draw the most important nodes with a different style
nx.draw_networkx_nodes(Gs_degree,s_coords,node_color='red',alpha=0.5,node_size=254)
# also the labels this time
nx.draw_networkx_labels(Gs_degree,s_coords,font_size=8,font_color='blue')
plt.savefig("serum_degree4+.png", dpi=1000)
plt.show()


# In[60]:


#report lists of degree >=5 nodes for each tissue
import csv
Brain_degree_5 = Gb_degree.nodes()
with open('brain_degree_4+_noPE_020717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Brain_degree_5:
        wr.writerow([i])
          
Liver_degree_5 = Gl_degree.nodes()
with open('liver_degree_4+_noPE_020717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Liver_degree_5:
        wr.writerow([i])
    
Muscle_degree_5 = Gm_degree.nodes()
with open('muscle_degree_4+_noPE_020717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Muscle_degree_5:
        wr.writerow([i])
    
Serum_degree_5 = Gs_degree.nodes()
with open('serum_degree_4+_noPE_020717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Serum_degree_5:
        wr.writerow([i])


# In[9]:


#list for all nodes in tissue networks

#report lists of degree >=5 nodes for each tissue
import csv
Brain = Gb.nodes()
with open('brain_network_all_nodes_040717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Brain:
        wr.writerow([i])
          
Liver = Gl.nodes()
with open('liver_network_all_nodes_040717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Liver:
        wr.writerow([i])
    
Muscle = Gm.nodes()
with open('muscle_network_all_nodes_040717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Muscle:
        wr.writerow([i])
          
Heart = Gh.nodes()
with open('heart_network_all_nodes_040717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Heart:
        wr.writerow([i])
    
Serum = Gs.nodes()
with open('serum_network_all_nodes_040717.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in Serum:
        wr.writerow([i])


# In[ ]:




