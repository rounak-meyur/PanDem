# -*- coding: utf-8 -*-
"""
Created on Fri May  1 21:03:18 2020

@author: Bhaskar
"""
"""
This code segregates the dataset collected into specific name wise or region wise
segregation
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

def need(i):
    want={
        'name':0,
        'region':1}
    return want.get(i,"Invalid Entry")

i=input('How do you want the data name or region wise: ') #User choice input
    
emp=pd.DataFrame()
directory=os.path.join(r'C:\Users\Bhaskar\Downloads\nssac-ncov-data-country-state\nssac-ncov-data-country-state')
os.chdir(directory)

for root,dirs,files in os.walk(directory): #grabs the root, directory and file names
    out=need(i)  # 0:name, 1:region
    i=0
    reg = pd.Series([])
    
    for file in files: #For all the files in the folder
        
        if file.endswith(".csv"): #Specifically looking for .csv files
           
           df=pd.read_csv(file)
           if out==0:   #if name wise data segregation is needed
               
               if i==0:
                   reg=df.name
               else:
                   reg=reg.append(df.name)
           else:
               if i==0:
                   reg=df.Region
               else:
                   reg=reg.append(df.Region)
               
               
        reg_un=reg.unique() #Gets all the unique name or region from the data
        i+=1
    for j in range(len(reg_un)): #Iterating through all the variable names
    
        k=0
        for file in files:
            
            if file.endswith(".csv"):
                
                df=pd.read_csv(file)
           
                if out==0:
                    
                    data=df[df.name==reg_un[j]] #Get all the correspoing data with name or region match
                else:
                    
                    data=df[df.Region==reg_un[j]]
                
                
                if k==0:
                    
                  
                    df1=emp.append(data,ignore_index=True)
                    
                else:
                    df1=df1.append(data,ignore_index=True)
                    
            k+=1
        
           
        #Create a new csv file region or name specific        
        df1.to_csv(r'C:\Users\Bhaskar\Spyder\PanDem\data/'+reg_un[j]+'.csv',index=False)
        
print("Sob files done ebar kaaj kor\n`")
print("--- %s seconds ---" % (time.time() - start_time))
          





   