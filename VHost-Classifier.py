# -*- coding: utf-8 -*-
"""
Created on Sun May 27 17:20:30 2018

@author: Ezra
"""
#!usr/bin/env python 3

#######################################################################################
#FUNCTIONS
#######################################################################################

def tax_groups (choice):
   
    group_one = []
    group_two = []
    group_three = []
   #make groups Phylum class order or phylum order family
    with open ('Phyla_in_NCBI.txt') as phyla, open('Order_in_NCBI.txt') as Orders, open('Class_in_NCBI.txt') as Classs, open('Families_in_NCBI.txt') as Families:
        pat = csv.reader(phyla)
        cat = csv.reader(Classs)
        oat = csv.reader(Orders)
        fat = csv.reader(Families)
        if choice == 'PCO':
            for rowe in pat:
                row = rowe[0]
                group_one.append(row)
            for rowe in cat:
                row = rowe[0]
                group_two.append(row)
            for rowe in oat:
                row = rowe[0]
                group_three.append(row)
        if choice == 'POF':
            for rowe in pat:
                row = rowe[0]
                group_one.append(row)
            for rowe in oat:
                row = rowe[0]
                group_two.append(row)
            for rowe in fat:
                row = rowe[0]
                group_three.append(row)
                
    #get rid of empty elements
    group_one = list(filter(None, group_one))
    group_one.append('NA')
    group_two = list(filter(None, group_two))
    group_two.append('NA')
    group_three = list(filter(None, group_three))
    group_three.append('NA')
    
    #sort out some discrepancy between databases
    if choice == 'PCO':
        group_three.append('Cetartiodactyla')
    if choice == 'POF':
        group_two.append('Cetartiodactyla')


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
    
    return g1,g2,g3,group_one,group_two,group_three



def host_locate(taxid,genera): 
    #don't worry about case
    rhost = 'NA'
    host = 'NA'
    #name_change- the direct viral prefix has a latin root so will misclassify. However it can be useful to assign host. 
    #If we have used it to assign host name_change = true and therefore don't delete it later in code. 
    #Ignore RHOST if we are looking at phage
    name_change = False
    taxid = taxid.lower()
    if 'virus' in taxid:
        host = taxid.split('virus')[0]
        rhost = taxid.split('virus')[1]
    if 'phage' in taxid:
        host = taxid.split('phage')[0]
    if 'viridae' in taxid:
        host = taxid.split('viridae')[0]
        rhost = taxid.split('viridae')[1]
    if 'viroid' in taxid:
        host = taxid.split('viroid')[0]
        rhost = taxid.split('viroid')[1]
    if 'virinae' in taxid:
        host = taxid.split('virinae')[0]
        rhost = taxid.split('virinae')[1]
    if 'satellite' in taxid:
        host = taxid.split('satellite')[0]
        rhost = taxid.split('satellite')[1]
    if 'virales' in taxid:
        host = taxid.split('virales')[0]
        rhost = taxid.split('virales')[1]
    if 'hiv-1' in taxid or 'hiv-2' in taxid:
        host = taxid
    lhost = host.split(' ')
    rlhost = rhost.split(' ')
    ref = len(lhost)-1        
    lhost = lhost+rlhost
    ##associated means it is not an isolate
    ##we ignore this - ver 13
    ##Common names are often used for common animals in place of genus species, this needs to be corrected
    ###Might be best to turn this into a .csv file and read through it?
    count = range(0,len(lhost))
    for c,w in zip(count,lhost):
            #if there is a genus name, make the next word genus+species providing a) there is a next word b) the next word isnt a comname or a sci name 
            change = True
            for g in genera:
                if w == g:
                    try:
                        if lhost[c+1] == '':
                            change = False
                            break
                        if change == True:
                            for nam,sci in zip(comname,sciname):
                                if lhost[c+1] == nam or lhost[c+1] == sci:
                                    change = False
                                    break
                        if change == True:
                            lhost[c+1] = lhost[c]+' '+lhost[c+1]
                            break
                        if change == False:
                            break
                    except:
                        break
                    
            ##convert common names to their scientific names so they can be used to look up NCBI database
            for nam,sci in zip(comname,sciname):
                if nam == w:
                    lhost[c] = sci
                    ##if we are dealing with influenza, find the strain - check it isn't illegitimate
                    if w == 'influenza':
                        
                        try:
                            lhost[c+4] = lhost[c+4].split('/')[1]
                            for ifs in infstrains:
                                    ifstrain = ifs[0]
                                    if lhost[c+4] == ifstrain:
                                        lhost[c+4] = 'na'
                                        break
                        except:
                            pass    
                        if lhost[c+1] == 'a':
                            lhost[c] = 'aves'

                        elif lhost[c+1] == 'b' or lhost[c+1] == 'c':
                            lhost[c] = 'mammalia'

                        else:
                            lhost[c] = 'mammalia'
                        
                        #move the influenza to the first position - spiny eel influenza virus - pop deletes and returns element
                        lhost.insert(0,lhost.pop(c))

                    #if we are dealing with norovirus
                    if w =='noro':
                     try: 
                        lhost[c+2] = lhost[c+2].split('/')[0]
                     except:
                         pass
                    ##if we have changed the virus prefix
                    if lhost[c] == lhost[ref]:
                        name_change = True
                        lhost.insert(0,lhost.pop(c))
                    break
    
    if name_change == False:
        lhost.remove(lhost[ref])
    ref = len(lhost)-1        
    if ref > -1:
        while ref > -1:
            name2taxid = ncbi.get_name_translator([lhost[ref]])
            if name2taxid:
                break
            else:
                ref = ref-1 
    else:
        name2taxid = {} 
    ##include the IMGER database
                
    return name2taxid     


def env_groups (iecos,iecoc,iecot,iecost):
    n = 0
    e = 0
    f = 0 
    r = 0 
    ecosys = ['']*1000
    ecocat = ['']*1000
    ecotyp = ['']*1000
    ecosub = ['']*1000
    ecos = iecos[1:len(iecos)]
    ecoc = iecoc[1:len(iecoc)]
    ecot = iecot[1:len(iecot)]
    ecost = iecost[1:len(iecost)]



    for item in ecos:
        if item not in ecosys:
                ecosys[n]=item
                n = n + 1 
    for item in ecoc:
        if item not in ecocat:
                ecocat[e]=item
                e = e + 1 
    for item in ecot:
        if item not in ecotyp:
                ecotyp[f]=item
                f = f+ 1         
    for item in ecost:
        if item not in ecosub:
                ecosub[r]=item
                r = r+ 1     

      
    #get rid of empty elements
    ecosys = list(filter(None, ecosys))
    ecosys.append('NA')
    ecocat = list(filter(None, ecocat))
    ecocat.append('NA')
    ecocat.append('Bioreactor')
    ecotyp = list(filter(None, ecotyp))
    ecotyp.append('NA')
    ecosub = list(filter(None, ecosub))
    ecosub.append('NA')
    


    e1 = {}
    keys = range(len(ecosys))
    values = ecosys
    for i in keys:
            e1[values[i]] = [i]
    e2 = {}
    keys = range(len(ecocat))
    values = ecocat
    for i in keys:
            e2[values[i]] = [i]  
    e3 = {}
    keys = range(len(ecotyp))
    values = ecotyp
    for i in keys:
            e3[values[i]] = [i]
    e4 = {}
    keys = range(len(ecosub))
    values = ecosub
    for i in keys:
            e4[values[i]] = [i]
    
    
    return e1,e2,e3,e4,ecosys,ecocat,ecotyp,ecosub


def environ_locate(name):
    name = name.lower()
    words = name.split(" ")
    ecoc = 'NA'
    ecos = 'NA'
    ecot = 'NA'
    ecost = 'NA'
    fd = False
    for w in words:
            if fd == True:
                break
            for nam,wa,wb,wc,wd in zip(dw,decoc,decos,decot,decost):
                if nam == w:
                    ecos = wa
                    ecoc = wb
                    ecot = wc
                    ecost = wd
                    fd = True
                    break
    return ecoc,ecos,ecot,ecost
         
                    
#################################################################################
##MAIN SCRIPT
###################################################################################


#Initial filter of all the data for viruses 
##The non-indent open and close here is a nice way to deal with multiple files that need to be read/written simulataneously 
# If a taxon ID doesn't have a kingdom assignment - it won't have a name
#newline = '' used to make each entry go to row directly below

###############################################################################
#INITALISE DATABASES
################################################################################
print('Initialising databases....')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("Taxon_IDs", help="List of Taxon IDs to classify")
parser.add_argument("vhost_db", help ="Name of vhost file")  
parser.add_argument("output_dir", help="Name of output directory")

parser.add_argument("-i","--index",type=int, help="value to index search terms from [default 0]")
parser.add_argument("-g","--groups",type=str,help="taxonomic groups to bin to [default PCO]")
parser.add_argument("-n","--names",type=str,help="parse file containing scientific names")

args = parser.parse_args()

if args.index:
   print ("Indexing search terms from:")
   print(str(args.index))
   
if args.names:  
    print("parsing names file")

if args.groups == 'PCO':
    print('Classifying virus host to Phylum Class Order')
    choice = args.groups

elif args.groups == 'POF':
    print('Classifying virus host to Phylum Order Family')
    choice = args.groups
    
else:
    print('Classifying virus host to Phylum Class Order')
    choice = 'PCO'

    
import os
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import csv

print('Filtering for viruses...')

##open com_sci
comname = []
sciname = []
cs2 = open('Common_to_Sci.csv') 
comsci = csv.reader(cs2)
for row in comsci:
    comname.append(row[0])
    sciname.append(row[1])

##create list of illegal influenza names
infstrains = []
with open('infstrains.txt') as ifs:
    infs = csv.reader(ifs)
    for row in infs:
        infstrains.append(row)


f2 = open(args.vhost_db, 'r')
c2 = csv.reader(f2,delimiter='\t')
vhostdb = list(c2)

##create taxonomic groups
(g1,g2,g3,group_one,group_two,group_three) = tax_groups(choice)

##create environmental groups
(e1,e2,e3,e4,ecosys,ecocat,ecotyp,ecosub) = env_groups(iecos,iecoc,iecot,iecost)


##need to remove the '/' character from the group list (this prevents the .csv file writing)
([s.strip('/') for s in group_one]) # remove the / from the string borders
group_one=([s.replace('/', '') for s in group_one]) # remove all the /s
([s.strip('/') for s in group_two]) # remove the 8 from the string borders
group_two=([s.replace('/', '') for s in group_two]) # remove all the /s
([s.strip('/') for s in group_three]) # remove the / from the string borders
group_three=([s.replace('/', '') for s in group_three]) # remove all the /s

#create genus list - used in binning
genera = []
import csv
with open ('Genus_in_NCBI.txt') as genus:
        sat = csv.reader(genus)
        for rowe in sat:
            row = rowe[0].lower()
            genera.append(row)


###############################################################################
##BEGIN HOST ASSIGNMENT 
##############################################################################
print('Classifying hosts')
os.makedirs(args.output_dir,exist_ok=True)
e1 = open(args.Taxon_IDs, 'r') 

##set this n to 0 to index from 0 
n = 0
if args.index:
    n = args.index
    
os.chdir(args.output_dir)
with open('Virus.csv','w',newline='') as e3, open('Non-Virus.csv','w',newline='') as e4:
    d1 = csv.reader(e1)
    d3 = csv.writer(e3)
    d4 = csv.writer(e4)
    for row in d1:
                taxid2name = ncbi.get_taxid_translator([row[0]])
                ##deal with taxonids in ncbi but not in the database
                try:
                    SN = (list(taxid2name.values())[0])
                    
               ##added below to support viruses not in NCBI
                except IndexError:
                    if args.names:
                        os.chdir('..')
                        with open(args.names) as r1:
                            #row.append(SN)
                            t1 = csv.reader(r1)
                            SN=[row for idx, row in enumerate(t1) if idx ==n]
                            SN = SN[0][0]
                            os.chdir(args.output_dir)
                    else:
                        SN = 'Not in database'
                
                SN = SN.lower()
                row.append(str(n))
                row.append(SN)
                if 'virus' in SN:
                    d3.writerow(row) 
                elif 'phage' in SN and not 'phagedenis' in SN:
                    d3.writerow(row)
                elif 'viridae' in SN:
                    d3.writerow(row)
                elif 'satellite' in SN:
                    d3.writerow(row)
                elif 'virinae' in SN:
                    d3.writerow(row)
                elif 'viroid' in SN:
                    d3.writerow(row)
                elif 'hiv-1' in SN or 'hiv-2' in SN:
                    d3.writerow(row)
                elif 'virales' in SN:
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
os.makedirs('Host Assigned',exist_ok=True)
os.chdir('Host Assigned')
f3 = open('Eukaryota.csv', 'w',newline='')
f4 = open('Bacteria.csv','w',newline='')
f5 = open('Virus.csv','w',newline='')
f6 = open('Archaea.csv','w',newline='')

c1 = csv.reader(f1)
c3 = csv.writer(f3)
c4 = csv.writer(f4)
c5 = csv.writer(f5)
c6 = csv.writer(f6)


icount = range(0,len(incbiids)) 
f9 = open('NAf.csv','w',newline='')
c9 = csv.writer(f9)

for row in c1:
    qtaxid = row[0]
    qname = row[2]
    found = False
    ##First parse the vhostdb
    for master_row in vhostdb:
        dbtaxid = master_row[0]
        host_name = master_row[8]
        if dbtaxid == qtaxid:
           #hn = [host_name.split(' ')[0]]
           name2taxid = ncbi.get_name_translator([host_name])
           if 'root' in name2taxid or 'bacteria' in name2taxid:
               name2taxid = []
           if not name2taxid:
               pass
           else:
               row.append('VHostDB')
               found = True
           break
    
    if found == False:
        ##if not present try to predict the host
        name2taxid = host_locate(qname,genera)
        
    
    try:
         lid = (list(name2taxid.values())[0])
    except IndexError:
            c9.writerow(row)
        
    else:      
        lineage = ncbi.get_lineage(lid[0])
        names = ncbi.get_taxid_translator(lineage)
        lineage_trace = ([names[taxid] for taxid in lineage])
        ##deal with things not in NCBI
        while len(lineage_trace) < 3:
            lineage_trace.append('NA')
            
        gr1 = 'NA'
        for i in range(0,len(lineage_trace)):
                if lineage_trace[i] in group_one:
                    gr1 = lineage_trace[i]
        gr2 = 'NA'    
        for i in range(0,len(lineage_trace)):
                if lineage_trace[i] in group_two:
                    gr2 = lineage_trace[i]
        gr3 = 'NA'    
        for i in range(0,len(lineage_trace)):
                if lineage_trace[i] in group_three:
                    gr3 = lineage_trace[i]   
                    
        if lineage_trace[2] == 'Eukaryota':
            c3.writerow(row)
            os.makedirs('Eukaryote',exist_ok=True)
            os.chdir('Eukaryote')
            key = g1[gr1]
            family = group_one[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
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
            ####BIN GROUP 3 
            os.makedirs(gr2,exist_ok=True)
            os.chdir(gr2)
            key = g3[gr3]
            family = group_three[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('../../..')
            
        elif lineage_trace[2] == 'Bacteria':
            c4.writerow(row)      
            os.makedirs('Bacteria',exist_ok=True)
            os.chdir('Bacteria')
            key = g1[gr1]
            family = group_one[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            ####BIN GROUP 2 BACTERIA
            os.makedirs(family,exist_ok=True)
            os.chdir(family)
            key = g2[gr2]
            family = group_two[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            ####BIN GROUP 3 EBACTERIA
            os.makedirs(family,exist_ok=True)
            os.chdir(family)
            key = g3[gr3]
            family = group_three[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('../../..')
            
        elif lineage_trace[2] == 'Viruses':
            c5.writerow(row)
            #bin group 1 entries        
            os.makedirs('Virus',exist_ok=True)
            os.chdir('Virus')
            key = g1[gr1]
            family = group_one[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            ####BIN GROUP 2 VIRUS
            os.makedirs(gr1,exist_ok=True)
            os.chdir(gr1)
            key = g2[gr2]
            family = group_two[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            ####BIN GROUP 3 VIRUS
            os.makedirs(gr2,exist_ok=True)
            os.chdir(gr2)
            key = g3[gr3]
            family = group_three[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            os.chdir('../../..')
            
        elif lineage_trace[2] == 'Archaea':
            c6.writerow(row)
            os.makedirs('Archaea',exist_ok=True)
            os.chdir('Archaea')
            key = g1[gr1]
            family = group_one[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            ####BIN GROUP 2 aRCHEA
            os.makedirs(gr1,exist_ok=True)
            os.chdir(gr1)
            key = g2[gr2]
            family = group_two[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            pr.writerow(row)
            d8.close()
            ####BIN GROUP 3 ARCHEA
            os.makedirs(gr2,exist_ok=True)
            os.chdir(gr2)
            key = g3[gr3]
            family = group_three[key[0]]
            d8 = open(family+'.csv','a',newline='')
            pr = csv.writer(d8)
            os.chdir('../../..')
            
            
            pr.writerow(row)
            d8.close()
        else:   
                c9.writerow(row)

        
cs2.close()          
f3.close()
f4.close()
f5.close()
f6.close()
f9.close()


##################################################################################
#Begin Host_Counts
#####################################################################################
print('Counting host assignments....') 
super_kingdoms = ['Eukaryota','Bacteria','Virus','Archaea','NA']

num_lines = ['Blank']*5
num_lines[0] = sum(1 for line in open('Eukaryota.csv'))
num_lines[1] = sum(1 for line in open('Bacteria.csv'))
num_lines[2] = sum(1 for line in open('Virus.csv'))
num_lines[3] = sum(1 for line in open('Archaea.csv'))
num_lines[4] = sum(1 for line in open('NAf.csv'))

with open ('Counts.csv','w') as counts:
     for i in zip(super_kingdoms,num_lines):
            co = csv.writer(counts)
            co.writerow(i)
            

##This will go through and do all the rest of counting 

for fname in os.listdir():
    path = fname
    if os.path.isdir(path):
        os.chdir(path)
        for g in group_one:
            try:  
                d8 = open('COUNTS.csv','a',newline='')
                pr = csv.writer(d8)
                num_lines = sum(1 for line in open(g+'.csv'))
                i = [g]+ [str(num_lines)]
                pr.writerow(i)
                d8.close()
                ###group 2 counts 
                os.chdir(g)
            
                for g2 in group_two:
                    try:
                            d8 = open('COUNTS.csv','a',newline='')
                            pr = csv.writer(d8)
                            num_lines = sum(1 for line in open(g2+'.csv'))
                            i = [g2]+ [str(num_lines)]
                            pr.writerow(i)
                            d8.close()
                            os.chdir(g2)
            
                            for g3 in group_three:
                                try:
                                    d8 = open('COUNTS.csv','a',newline='')
                                    pr = csv.writer(d8)
                                    num_lines = sum(1 for line in open(g3+'.csv'))
                                    i = [g3]+ [str(num_lines)]
                                    pr.writerow(i)
                                    d8.close()
                                except FileNotFoundError:
                                        pass        
                            os.chdir('..')
                    except FileNotFoundError:
                        pass        
                os.chdir('..')
            except FileNotFoundError:
                    pass
                
        os.chdir('..')            
            


#Finally do the host-assigned/unassigned counts (ugly code)    

d8 = open('COUNTS.csv','a',newline='')
pr = csv.writer(d8)
os.chdir('Host Assigned')
bum_lines = sum(1 for line in open('NAf'+'.csv'))
i = ['Host Unassigned']+ [str(bum_lines)]
os.chdir('..')
pr.writerow(i)
d8.close()
d8 = open('COUNTS.csv','a',newline='')
pr = csv.writer(d8)
os.chdir('..')
num_lines = sum(1 for line in open('Virus'+'.csv'))
os.chdir('Virus')
i = ['Host Assigned']+ [str(num_lines-bum_lines)]
pr.writerow(i)
d8.close()

          
f1.close()
f2.close()

