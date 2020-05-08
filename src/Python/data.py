# -*- coding: utf-8 -*-
"""
Created on Fri May  1 21:03:18 2020

@author: Bhaskar
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
import csv,re,sys,time


'''
Function not working as of now to read the column names in a csv file
'''
start_time=time.time()
    
emp=pd.DataFrame()
directory=os.path.join(r'C:\Users\Bhaskar\Downloads\nssac-ncov-data-country-state\nssac-ncov-data-country-state')
os.chdir(directory)
#name = input("Enter region name: ")
for root,dirs,files in os.walk(directory):
    i=0
    reg = pd.Series([])
    for file in files:
        
        if file.endswith(".csv"):
           
           df=pd.read_csv(file)
           cols=df.columns
           if i==0:
               reg=df.name
           else:
               reg=reg.append(df.name)
               
        reg_un=reg.unique() 
        i+=1
    for j in range(len(reg_un)):
        k=0
        for file in files:
            
            if file.endswith(".csv"):
                
                df=pd.read_csv(file)
           
                data=df[df.name==reg_un[j]]
                
                if k==0:
                    
                  
                    df1=emp.append(data,ignore_index=True)
                    
                else:
                    df1=df1.append(data,ignore_index=True)
                    
            k+=1
        
           
                 
        df1.to_csv(r'C:\Users\Bhaskar\Spyder\PanDem\data/'+reg_un[j]+'.csv',index=False)
        
print("Sob files done ebar kaaj kor\n")
print("--- %s seconds ---" % (time.time() - start_time))
          





   