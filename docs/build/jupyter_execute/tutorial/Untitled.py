#!/usr/bin/env python
# coding: utf-8

# In[1]:


import dill


# In[2]:


class A:
    
    dic = {
        'a': 1,
        'b': 2
    }
    
    @classmethod
    def change(cls, dic):
        cls.dic = dic


# In[3]:


a = A
a.dic


# In[4]:


a.change({'a': 4, 'b': 5})
a.dic


# In[5]:


b = A
b.dic


# In[6]:


with open('save.pkl', 'wb') as f:
    dill.dump(a, f)


# In[7]:


with open('save.pkl', 'rb') as f:
    c = dill.load(f)


# In[8]:


c.dic


# In[ ]:




