# VHost-Classifier
For a list of taxonIDs, VHost-Classifier will filter out the viruses and then sort these viruses into groups based on their host lineage.
The VHost-Classifier algorithm uses the Virus-Host DB , the NCBI Taxonomy DB and in-built predictive rules to achieve a high rate of virus host classification. VHost-Classifier will classify virus taxonIDs to family resolution. 

**Usage:**
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
python VHost-Classifier.py TaxonIDs.tsv VirusHostDB.tsv VHC_Run_1 -i 1 -g POF -n Sci_Names.csv
``` 
Virus host classify a list of taxonIDs in ```TaxonIDs.tsv```, use the VHost-DB file supplied by ```VirusHostDB.tsv``` and output the results to ```VHC_RUN_1```. Index the input taxonIDs from 1 in the output csv files. Classify taxonIDs to Phylum Order Family. Parse the ```Sci_Names.csv``` file. 

**Dependencies:**<br/>
Python 3 <br/>
[ETE3 Toolkit for Python 3](http://etetoolkit.org/download/)  
Note:  On first run through NCBI taxonomy database will be downloaded by ETE3.  

**Output**:
VHost Classifier will![create directories](https://github.com/Kzra/VHost-Classifier/blob/master/Dir%20navigation%20example.pdf) and in each directory write .csv files.<br/>
**Reading the .csv files**: the first column contains taxon IDs, the second column the index position (indexed from -i) of the taxon id in the input file. The final column contains the host name, predicted host name, or virus name (if host name can't be predicted). In each directory a counts.csv file is also written which contains the counts of how many taxon IDs are in each taxonomic class. 

**Reference:**
Virus Host DB:
Mihara, Tomoko, et al. "Linking virus genomes with host taxonomy." Viruses 8.3 (2016): 66.
