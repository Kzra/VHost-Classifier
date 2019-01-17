# VHost-Classifier
For a list of taxonIDs, VHost-Classifier will filter out the viruses and then sort these viruses into groups based on their host lineage.

The [VHost-Classifier algorithm][2] uses the Virus-Host DB, the NCBI Taxonomy DB and inbuilt predictive rules to achieve a high rate of virus host classification. VHost-Classifier will classify virus taxonIDs to family resolution. 

VHost-Classifier will sort viruses it could not assign a host to by the environment they were sequenced from. To do this it uses the IMG/VR database and inbuilt predictive rules. 

When [benchmarked][3] on 1000 randomly selected viral taxonids on NCBI, the software could classify 93% of vtaxids to the rank of Class, and 37% of vtaxids to the rank of Family, with an accuracy of 100%. A list of these random taxids can be found in the random_ids.csv file.

**Usage:**

Clone the directory and run from within cloned directory.

```shell
python vhost_classifier.py [TaxonID.tsv] [VirusHostDB.tsv] [Output Dir] [-i] [-g] [-n]
```

```[TaxonID.tsv]```: a .tsv list of taxonIDs to be classified (one taxon ID per row).

```[VHostDB.tsv]```: a copy of the Virus Host DB which can be downloaded [here](http://www.genome.jp/virushostdb/)</br>
    or by running :
```wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv```

```[Output Dir] ```: the name of the directory to output results to (must be unique). 

```[-i]```: optional argument, specify the value to start indexing the input taxonIDs from (default 0). 

```[-g]```: optional argument, taxonomic ranks to bin to. PCO, Phylum Class Order or POF, Phylum Order Family (default PCO). 

```[-n]```: optional argument, supply file of scientific names alongside taxon ids (use if taxonid list returns an index error).  

**Example**:

```shell
python VHost_Classifier.py TaxonIDs.tsv VirusHostDB.tsv VHC_Run_1 -i 1 -g POF -n Sci_Names.csv
``` 
Virus host classify a list of taxonIDs in ```TaxonIDs.tsv```, use the VHost-DB file supplied by ```VirusHostDB.tsv``` and output the results to ```VHC_RUN_1```. Index the input taxonIDs from 1 in the output csv files. Classify taxonIDs to Phylum Order Family. Parse the ```Sci_Names.csv``` file. 

**Dependencies:**<br/>
Python 3 <br/>
[ETE3 Toolkit for Python 3](http://etetoolkit.org/download/)  
Note:  On first run through NCBI taxonomy database will be downloaded by ETE3.  

**Output**:
VHost Classifier will [create directories][1] and in each directory write .csv files.<br/>

**Reading the .csv files**: the first column contains taxon IDs, the second column the index position (indexed from -i) of the taxon id in the input file. The final column contains the virus name. In each directory a counts.csv file is also written which contains the counts of how many taxon IDs are in each taxonomic class. 

**VHC-Analysis:** run this script from within the Host-Assigned directory of the run you want to analyse. The script will write walk the directory tree and write each Counts.csv file to a Total_Counts.csv file which will be saved in the Host-Assigned directory. This file makes it easier to compare the overall host diversity of viruses in your input.   

**References:**<br/>
Virus-Host DB: Mihara, Tomoko, et al. "Linking virus genomes with host taxonomy." Viruses 8.3 (2016): 66.

IMG/VR: Paez-Espino, David, et al. "IMG/VR: a database of cultured and uncultured DNA Viruses and retroviruses." Nucleic acids research (2016): gkw1030.

*If you use this software for publication please cite my github.*



[1]:https://github.com/Kzra/VHost-Classifier/blob/master/Directory%20Navigation%20Example.pdf
[2]:https://github.com/Kzra/VHost-Classifier/blob/master/Host%20Classification.pdf
[3]:https://github.com/Kzra/VHost-Classifier/blob/master/benchmark.png
