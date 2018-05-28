# -*- coding: utf-8 -*-
"""
Created on Sun May 27 17:20:30 2018

@author: Ezra
"""
#!usr/bin/env python 3

##A python script that will filter viruses out of a list of Taxon IDs and then sort the viruses according to the host they infect.

#27/05/2018 - Script now fully functional, run time of 1min 48 secs on 436000 list with 2000 viruses.
# To do: Neaten up the code
# add functionality to specify taxonomic resolution 


##The non-indent open and close here is a nice way to deal with multiple files that need to be read/written simulataneously 
# If a taxon ID doesn't have a kingdom assignment - it won't have a name
#newline = '' used to make each entry go to row directly below

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("Taxon_IDs", help="List of Taxon IDs to classify")
parser.add_argument("vhost_db", help ="Name of vhost file")  
parser.add_argument("output_dir", help="Name of output directory")

parser.add_argument("-i","--index",type=int, help="value to index search terms from [default 0]")
args = parser.parse_args()

if args.index:
   print ("Indexing search terms from:")
   print(str(args.index))
   
    
import os
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import csv

print('Filtering for viruses...')


#Initial filter of all the data for viruses 


os.makedirs(args.output_dir,exist_ok=True)
e1 = open(args.Taxon_IDs, 'r') 
f2 = open(args.vhost_db, 'r')

##set this n to 0 to index from 0 
n = 0
if args.index:
    n = args.index
    
os.chdir(args.output_dir)
with open('Virus.csv','w',newline='') as e3, open('Non-Virus.csv','w',newline='') as e4:
    d1 = csv.reader(e1,delimiter='\t')
    d3 = csv.writer(e3)
    d4 = csv.writer(e4)
    for row in d1:
                taxid2name = ncbi.get_taxid_translator([row[0]])
                SN = (list(taxid2name.values())[0])
                row.append(str(n))
                
                if 'virus' in SN:
                    d3.writerow(row) 
                elif 'phage' in SN and not 'phagedenis' in SN:
                    d3.writerow(row)
                elif 'viridae' in SN:
                    d3.writerow(row)
                else:
                    d4.writerow(row)
                n = n+1    
Lifeform = ['Virus','NonVirus']
num_lines = ['Blank']*2
num_lines[0] = sum(1 for line in open('Virus.csv'))
num_lines[1] = sum(1 for line in open('Non-Virus.csv'))
with open ('Counts.csv','w') as counts:
     for i in zip(Lifeform,num_lines):
            co = csv.writer(counts)
            co.writerow(i)
e1.close()


os.makedirs('Virus',exist_ok=True)



f1 = open('Virus.csv', 'r')
os.chdir('Virus')
f3 = open('Eukaryota.csv', 'w',newline='')
f4 = open('Bacteria.csv','w',newline='')
f5 = open('Virus.csv','w',newline='')
f6 = open('Archea.csv','w',newline='')
f7 = open('NA.csv','w',newline='')

c1 = csv.reader(f1)
c2 = csv.reader(f2,delimiter='\t')
c3 = csv.writer(f3)
c4 = csv.writer(f4)
c5 = csv.writer(f5)
c6 = csv.writer(f6)
c7 = csv.writer(f7)



print('Analysing the vhost database...')

vhostdb = list(c2)
##get the unique groups, write them into a dictionary
#The dictionary is cool
def tax_groups (vhostdb):
    n = 0
    e = 0
    f = 0 
    
    group_one = ['']*1000
    group_two = ['']*1000
    group_three = ['']*1000
    for row in vhostdb[1:len(vhostdb)]:
        host_lineage = row[9]
        groups = host_lineage.split(';')
        try:
            gr1 = groups[1].split(' ')[1]
            gr2 = groups[2].split(' ')[1]
            gr3 = groups[3].split(' ')[1]
        except:
                pass
       
        if gr1 not in group_one:
                group_one[n]=gr1
                n = n + 1 
       
        if gr2 not in group_two:
                group_two[e]=gr2
                e = e + 1 
      
        if gr3 not in group_three:
                group_three[f]=gr3
                f = f+ 1
                
      
      


    #get rid of empty elements
    group_one = list(filter(None, group_one))
    group_one.append('NA')
    group_two = list(filter(None, group_two))
    group_two.append('NA')
    group_three = list(filter(None, group_three))
    group_three.append('NA')
  

    g1 = {}
    keys = range(len(group_one))
    values = group_one
    for i in keys:
            g1[values[i]] = [i]

    g2 = {}
    keys = range(len(group_two))
    values = group_two
    for i in keys:
            g2[values[i]] = [i]  

    g3 = {}
    keys = range(len(group_three))
    values = group_three
    for i in keys:
            g3[values[i]] = [i]
    
    return g1,g2,g3,group_one,group_two,group_three,groups

(g1,g2,g3,group_one,group_two,group_three,groups) = tax_groups(vhostdb)

##need to remove the '/' character from the group list (this prevents the .csv file writing)
([s.strip('/') for s in group_one]) # remove the / from the string borders
group_one=([s.replace('/', '') for s in group_one]) # remove all the /s
([s.strip('/') for s in group_two]) # remove the 8 from the string borders
group_two=([s.replace('/', '') for s in group_two]) # remove all the /s
([s.strip('/') for s in group_three]) # remove the / from the string borders
group_three=([s.replace('/', '') for s in group_three]) # remove all the /s

for results_row in c1:
    qtaxid = results_row[0]
    row = 0
    found = False
    for master_row in vhostdb:
        dbtaxid = master_row[0]
        host_name = master_row[8]
        host_lineage = master_row[9]
        groups = host_lineage.split(';')
        try:
            gr1 = groups[1].split(' ')[1]
        except:
            gr1 = 'NA'
        try:    
            gr2 = groups[2].split(' ')[1]
        except:
            gr2 = 'NA'
        try:    
            gr3 = groups[3].split(' ')[1]
        except:
            gr3 = 'NA'
            
        super_kingdom = host_lineage.split(';')[0]
        if dbtaxid == qtaxid:
            if super_kingdom == 'Eukaryota':
                results_row.append(host_name)
                c3.writerow(results_row)
                found = True
                ##BIN GROUP 1 EUKARYOTES

                os.makedirs('Eukaryote',exist_ok=True)
                os.chdir('Eukaryote')
                key = g1[gr1]
                family = group_one[key[0]]
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                d8.close()
                #os.chdir('..')
                ####BIN GROUP 2 EUKARYOTES
                os.makedirs(gr1,exist_ok=True)
                os.chdir(gr1)
                key = g2[gr2]
                family = group_two[key[0]]
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                d8.close()
               
                #os.chdir('../..')
                ####BIN GROUP 3 EUKARYOTES
                os.makedirs(gr2,exist_ok=True)
                os.chdir(gr2)
                key = g3[gr3]
                family = group_three[key[0]]
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                os.chdir('../../..')
                   
                
                break
            if super_kingdom == 'Bacteria':
                results_row.append(host_name)
                c4.writerow(results_row)
                found = True
                ### Bin group 1 Bacteria
                os.makedirs('Bacteria',exist_ok=True)
                os.chdir('Bacteria')
                key = g1[gr1]
                family = group_one[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                d8.close()
                #os.chdir('..')
                #break
                ####BIN GROUP 2 BACTERIA
                os.makedirs(family,exist_ok=True)
                os.chdir(family)
                key = g2[gr2]
                family = group_two[key[0]]
                # 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                d8.close()
                #os.chdir('../..')
                ####BIN GROUP 3 BACTERIA
                os.makedirs(family,exist_ok=True)
                os.chdir(family)
                key = g3[gr3]
                family = group_three[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                os.chdir('../../..')
                break
            if super_kingdom == 'Archea':
                results_row.append(host_name)
                c5.writerow(results_row)
                found = True
                #### IN GROUP 1 ARCHEA
                os.makedirs('Archea',exist_ok=True)
                os.chdir('Archea')
                key = g1[gr1]
                family = group_one[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                d8.close()
                #os.chdir('..')
                ##BIN GROUP 2 ARCHEA
                os.makedirs(gr1,exist_ok=True)
                os.chdir(gr1)
                key = g2[gr2]
                family = group_two[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                d8.close()
                #now do the count, this is the tricky part....
                #os.chdir('../..')
                #os.chdir('..')
                ####BIN GROUP 3 ARCHEA
                os.makedirs(gr2,exist_ok=True)
                os.chdir(gr2)
                key = g3[gr3]
                family = group_three[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                os.chdir('../../..')
                break
            if super_kingdom == 'Virus':
                #results_row.append(rrow)
                results_row.append(host_name)
                #results_row.append(row)
                c6.writerow(results_row)
                found = True
                #####
                os.makedirs('Virophage',exist_ok=True)
                os.chdir('Virophage')
                key = g1[gr1]
                family = group_one[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                d8.close()
                #os.chdir('..')
                #BIN GROUP 2 VIRUS
                os.makedirs(gr1,exist_ok=True)
                os.chdir(gr1)
                key = g2[gr2]
                family = group_two[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                d8.close()
                #now do the count, this is the tricky part....
                #os.chdir('../..')
                #os.chdir('..')
                ####BIN GROUP 3 VIRUS
                os.makedirs(gr2,exist_ok=True)
                os.chdir(gr2)
                key = g3[gr3]
                family = group_three[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(results_row)
                os.chdir('../../..')
                break
            break
        row = row + 1
    if not found:
        #results_row.append('NOT FOUND in master list')
        #print('not found',str(row))
        #results_row.append(rrow)
        c7.writerow(results_row)
        
    #rrow = rrow +1

f1.close()
f2.close()
#f3.close()
#f4.close()
#f5.close()
#f6.close()
f7.close()

#Figure out what the unassigned viral IDs are
#One issue - Phage - could be archea or bacteria
#For now - it will just assign phage to bacteria 
print('Analysing NCBI database...')


f8 = open('NA.csv','r')
c8 = csv.reader(f8)

f9 = open('NAf.csv','w',newline='')
c9 = csv.writer(f9)

for row in c8:
    q = row[0]
    taxid2name = ncbi.get_taxid_translator([q])
    taxid = (list(taxid2name.values())[0])
    host = (taxid.split('virus')[0])
    fhost = host.split(' ')[0]
    name2taxid = ncbi.get_name_translator([fhost])
    try:
         lid = (list(name2taxid.values())[0])
#taxid2name = ncbi.get_taxid_translator([2828])
#print(taxid2name)
    #these SN's aren't in NCBI
    #order of criteria search is important
    except IndexError:
        if 'phycodna' in host or 'mimi' in host or 'Klosneu' in host or 'Cato' in host or 'Hoko' in host or 'Indi' in host:
            #print('Eukaryota')
            row.append(host)
            c3.writerow(row)
            gr1 ='NA' 
            gr2 ='NA' 
            gr3 ='NA'
            #Bin group 1 entries
            os.makedirs('Eukaryote',exist_ok=True)
            os.chdir('Eukaryote')
            key = g1[gr1]
            family = group_one[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('..')
        elif 'sipho' in host or 'caudo' in host or 'podo' in host or 'myo' in host:
            #print('Bacteria or Archea')
            row.append(host)
            c4.writerow(row)
            gr1 = 'NA'
            gr2 = 'NA'
            gr3 = 'NA'
            os.makedirs('Bacteria',exist_ok=True)
            os.chdir('Bacteria')
            key = g1[gr1]
            family = group_one[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('..')
        elif 'virophage' in host:
            #print('Virus')
            row.append(host)
            c5.writerow(row)
            gr1 = 'NA'
            gr2 = 'NA'
            gr3 = 'NA'
            os.makedirs('Virophage',exist_ok=True)
            os.chdir('Virophage')
            key = g1[gr1]
            family = group_one[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('..')
        elif 'phage' in host and not 'virophage' in host:
            #print('Bacteria or Archea')
            row.append(host)
            c4.writerow(row)
            gr1 = 'NA'
            gr2 = 'NA'
            gr3 = 'NA'
            os.makedirs('Bacteria',exist_ok=True)
            os.chdir('Bacteria')
            key = g1[gr1]
            family = group_one[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('..')
        else:
            #print('Didnt Work')
            row.append(host)
            c9.writerow(row)
        
    else:      
        lineage = ncbi.get_lineage(lid[0])
        names = ncbi.get_taxid_translator(lineage)
        lineage_trace = ([names[taxid] for taxid in lineage])
        gr1 = 'NA'
        for i in range(0,len(lineage_trace)):
                if lineage_trace[i] in group_one:
                    gr1 = lineage_trace[i]
        gr2 = 'NA'    
        for i in range(0,len(lineage_trace)):
                if lineage_trace[i] in group_two:
                    gr2 = lineage_trace[i]
        if lineage_trace[2] == 'Eukaryota':
            row.append(host)
            c3.writerow(row)
            ##bin group 1 entries
            #gr1 = 'NA'
            #for i in range(0,len(lineage_trace)):
                #if lineage_trace[i] in group_one:
                    #gr1 = lineage_trace[i]
            
        
            os.makedirs('Eukaryote',exist_ok=True)
            os.chdir('Eukaryote')
            key = g1[gr1]
            family = group_one[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            #os.chdir('..')
            ##bin group 2 entries
            
            os.makedirs(gr1,exist_ok=True)
            os.chdir(gr1)
            key = g2[gr2]
            family = group_two[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            
            ####BIN GROUP 3 EUKARYOTES
            os.makedirs(gr2,exist_ok=True)
            os.chdir(gr2)
            key = g3[gr3]
            family = group_three[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('../../..')
                
            
        elif lineage_trace[2] == 'Bacteria':
            row.append(host)
            c4.writerow(row)      
            os.makedirs('Bacteria',exist_ok=True)
            os.chdir('Bacteria')
            key = g1[gr1]
            family = group_one[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            #os.chdir('..')
            ####BIN GROUP 2 BACTERIA
            os.makedirs(family,exist_ok=True)
            os.chdir(family)
            key = g2[gr2]
            family = group_two[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            #now do the count, this is the tricky part....
            #os.chdir('../..')
            #os.chdir('..')
            ####BIN GROUP 3 EBACTERIA
            os.makedirs(family,exist_ok=True)
            os.chdir(family)
            key = g3[gr3]
            family = group_three[key[0]]
            ### 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('../../..')
            
        elif lineage_trace[2] == 'Virus':
            row.append(host)
            c5.writerow(row)
            #bin group 1 entries        
            os.makedirs('Virophage',exist_ok=True)
            os.chdir('Virophage')
            key = g1[gr1]
            family = group_one[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            #os.chdir('..')
            ####BIN GROUP 2 VIRUS
            os.makedirs(gr1,exist_ok=True)
            os.chdir(gr1)
            key = g2[gr2]
            family = group_two[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            #now do the count, this is the tricky part....
            #os.chdir('../..')
            #os.chdir('..')
            ####BIN GROUP 3 VIRUS
            os.makedirs(gr2,exist_ok=True)
            os.chdir(gr2)
            key = g3[gr3]
            family = group_three[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('../../..')
            
        elif lineage_trace[2] == 'Archea':
            row.append(host)
            c6.writerow(row)
            os.makedirs('Archea',exist_ok=True)
            os.chdir('Archea')
            key = g1[gr1]
            family = group_one[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            #os.chdir('..')
            ####BIN GROUP 2 aRCHEA
            os.makedirs(gr1,exist_ok=True)
            os.chdir(gr1)
            key = g2[gr2]
            family = group_two[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            #now do the count, this is the tricky part....
            #os.chdir('../..')
            #os.chdir('..')
            ####BIN GROUP 3 ARCHEA
            os.makedirs(gr2,exist_ok=True)
            os.chdir(gr2)
            key = g3[gr3]
            family = group_three[key[0]]
            ## 'a' to append to file
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('../../..')
            
            
        else:   
            #these lineages aren't in the ncbi - NEED TO WRITE GROUP 1 BINNING HERE AT SOME POINT
            if 'phycodna' in host or 'mimi' in host or 'Klosneu' in host or 'Cato' in host or 'Hoko' in host or 'Indi' in host:
                #print('Eukaryota')
                row.append(host)
                c3.writerow(row)
                gr1 ='NA' 
                gr2 ='NA' 
                gr3 ='NA'
                #Bin group 1 entries
                os.makedirs('Eukaryote',exist_ok=True)
                os.chdir('Eukaryote')
                key = g1[gr1]
                family = group_one[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(row)
                d8.close()
                os.chdir('..')
            elif 'sipho' in host or 'caudo' in host or 'podo' in host or 'myo' in host:
                #print('Bacteria or Archea')
                row.append(host)
                c4.writerow(row)
                gr1 = 'NA'
                gr2 = 'NA'
                gr3 = 'NA'
                os.makedirs('Bacteria',exist_ok=True)
                os.chdir('Bacteria')
                key = g1[gr1]
                family = group_one[key[0]]
                ## 'a' to append to file
                d8 = open(family+'.csv','a',newline='')
                pr = csv.writer(d8)
                pr.writerow(row)
                d8.close()
                os.chdir('..')   
            #elif 'sipho' in host:
                #print('Bacteria or Archea')
                #row.append(host)
                #c4.writerow(row)
            #elif 'caudo' in host:
                #print('Bacteria')
             #   row.append(host)
              #  c4.writerow(row)
            #elif 'podo' in host:
                #print('Bacteria')
             #   row.append(host)
              #  c4.writerow(row)
            #elif 'myo' in host:
                #print('Bacteria')
             #   row.append(host)
              #  c4.writerow(row)
            elif 'virophage' in host:
                #print('Virus')
                 row.append(host)
                 c5.writerow(row)
            elif 'phage' in host and not 'virophage' in host:
                #print('Bacteria or Archea')
                 row.append(host)
                 c4.writerow(row)
                 gr1 = 'NA'
                 gr2 = 'NA'
                 gr3 = 'NA'
                 os.makedirs('Bacteria',exist_ok=True)
                 os.chdir('Bacteria')
                 key = g1[gr1]
                 family = group_one[key[0]]
                 ## 'a' to append to file
                 d8 = open(family+'.csv','a',newline='')
                 pr = csv.writer(d8)
                 pr.writerow(row)
                 d8.close()
                 os.chdir('..')   
            else:
                #print('Didnt Work')
                row.append(host)
                c9.writerow(row)

        
##Count the number of rows in each .csv 
#d4 = csv.reader(f4)
#row_count = sum(1 for row in d4)  
#print(row_count)

            
f3.close()
f4.close()
f5.close()
f6.close()
f8.close()
f9.close()
os.remove('NA.csv')



###This part of the script will do all the counts
super_kingdoms = ['Eukaryota','Bacteria','Virus','Archea','NA']

num_lines = ['Blank']*5
num_lines[0] = sum(1 for line in open('Eukaryota.csv'))
num_lines[1] = sum(1 for line in open('Bacteria.csv'))
num_lines[2] = sum(1 for line in open('Virus.csv'))
num_lines[3] = sum(1 for line in open('Archea.csv'))
num_lines[4] = sum(1 for line in open('NAf.csv'))

with open ('Counts.csv','w') as counts:
     for i in zip(super_kingdoms,num_lines):
            co = csv.writer(counts)
            co.writerow(i)
            

##This will go through and do all the rest of counting 

for fname in os.listdir():
    path = fname
    #print(fname)
    if os.path.isdir(path):
        #print(path)
        os.chdir(path)
        for g in group_one:
            try:  
                d8 = open('COUNTS.csv','a',newline='')
                pr = csv.writer(d8)
                num_lines = sum(1 for line in open(g+'.csv'))
                i = [g]+ [str(num_lines)]
                #i in zip(g,num_lines):
                pr.writerow(i)
                #print(num_lines)
                d8.close()
                ###group 2 counts 
                os.chdir(g)
                #print(g)
            
                for g2 in group_two:
                    try:
                            d8 = open('COUNTS.csv','a',newline='')
                            pr = csv.writer(d8)
                            num_lines = sum(1 for line in open(g2+'.csv'))
                            i = [g2]+ [str(num_lines)]
                            #i in zip(g,num_lines):
                            pr.writerow(i)
                            #print(num_lines)
                            d8.close()
                            os.chdir(g2)
                            #print(g2)
            
                            for g3 in group_three:
                                try:
                                    d8 = open('COUNTS.csv','a',newline='')
                                    pr = csv.writer(d8)
                                    num_lines = sum(1 for line in open(g3+'.csv'))
                                    i = [g3]+ [str(num_lines)]
                                    #i in zip(g,num_lines):
                                    pr.writerow(i)
                                    #print(num_lines)
                                    d8.close()
                                except FileNotFoundError:
                                       #print('ezra')
                                        pass        
                            os.chdir('..')
                    except FileNotFoundError:
                       #print('ezra')
                        pass        
                os.chdir('..')
            except FileNotFoundError:
                   #print('ezra')
                    pass
                
        os.chdir('..')            
            
print('Classification Complete!')        

os.chdir('../..')
#!dir
 
    
#list = os.listdir(dir) # dir is your directory path
                #number_files = len(list)
                #num_lines = ['']*number_files    
            

