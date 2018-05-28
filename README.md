# VHost-Classifier
For a list of taxon IDs, VHost-Classifier will filter out the viruses and then sort these viruses into groups based on their host lineage.
The VHost-Classifier algorithm uses the Virus-Host DB , the NCBI Taxonomy DB and in-built predictive rules to acheive a high rate of virus host classification. VHost-Classifier will classify virus taxon IDs to family resolution. 

**Usage:**
```shell
!python vhost_classifier.py [TaxonID.tsv] [VirusHostDB.tsv] [Output Dir] [-i]
```

```[TaxonID.tsv]```: a .tsv list of TaxonIDs to be classified (one Taxon ID per row).

```[VHostDB.tsv]```: a copy of the Virus Host DB which can be downloaded [here](http://www.genome.jp/virushostdb/).

```[Output Dir] ```: the name of the directory to output results to (must be unique). 

```[-i]```: optional argument, specify the value to start indexing the input TaxonIDs from (default 0). 


**Dependencies:**<br/>
Python 3 <br/>
[ETE3 Toolkit for Python 3](http://etetoolkit.org/download/)  
Note: On first run through NCBI taxonomy database will be downloaded by ETE3.  

**Output**:
VHost Classifier will write create directories and in each directory write .csv files.<br/>
**Reading the .csv files**: the first column contains Taxon IDs, the second column the index position (indexed from -i) of the taxon id in the input file. The final column contains the host name, predicted host name, or virus name (if host name can't be predicted). In each directory a counts.csv file is also written which contains the counts of how many Taxon IDs are in each taxonomic class. 

**Reference:**
Virus Host DB:
Mihara, Tomoko, et al. "Linking virus genomes with host taxonomy." Viruses 8.3 (2016): 66.
