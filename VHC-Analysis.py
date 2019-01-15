# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 14:57:54 2018

@author: Ezra
"""

##This script should be run from within the 'Host Assigned' directory of the run you want to analyse. 
#It will go through and write each 'Counts.csv' file to a single excel file available in the Host Assigned directory called 'Total_Counts'
#Having all the counts in one place makes it easy to do comparisons and analysis of host composition.

import os 
import csv


##One / write to column one, two // write to column two, three /// write to column

#os.chdir('Run7_270618/Virus/Host Assigned')
cwd = os.getcwd()
t1 = open('Total Counts.csv','w')
tw =  csv.writer(t1,lineterminator='\n')

col1 = []
col2 = []
col3 = []

for root,dirs,files in os.walk('.'):
    rc = root.count('\\')
    if rc == 3:
      col3.append([root,' '])
    if rc == 2:
      col2.append([root,' '])
    if rc == 1:
      col1.append([root,' '])
    root1 = (cwd+'/'+ root)
    os.chdir(cwd+'/'+root)
    for f in files:
        if f == 'COUNTS.csv': 
            with open(f) as count:
                counts = csv.reader(count)
                for row in counts:
                    if rc == 3:
                        col3.append(row)
                    if rc == 2:
                        col2.append(row)
                    if rc == 1:
                        col1.append(row)
                
    os.chdir(cwd)
    
    
lcol1 = len(col1)
lcol2 = len(col2)
lcol3 = len(col3)

ml = max(lcol1,lcol2,lcol3)

while len(col1) < ml:
    col1.append([' ',' '])
while len(col2) < ml:
    col2.append(['  ',' '])
while len(col3) < ml:
    col3.append(['  ','  ']) 
    
for  i in zip(col1,col2,col3):
          tw.writerow(i[0]+i[1]+i[2])

t1.close()          


