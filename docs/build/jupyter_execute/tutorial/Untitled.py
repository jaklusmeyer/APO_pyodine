#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Automatic reloading of imports
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import sys
import os
import matplotlib.pyplot as plt
import h5py

sys.path.append('/home/paul/pyodine/')

import pyodine


# In[2]:


a = {'1': 1, '2': 2}

b = {'3': 3, 'four': 4}

c = {'a': a, 'b': b}


# In[18]:


with h5py.File('blub.h5', 'w') as h:
    pyodine.lib.h5quick.dict_to_group(c, h, 'blub')


# In[22]:


with h5py.File('blub.h5', 'r') as h:
    print(h['blub/a/1'][()])


# In[ ]:




