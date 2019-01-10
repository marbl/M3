# Winter 2019 Mid-Atlantic Microbiome Meetup
# Taxonomic Identification Workshop

## Logging onto the server
There is username and password at the back of your name tag.
```bash
ssh USERNAME@openclass.umiacs.umd.edu
```
Enter your password. Once you are logged in, you should be able to see m3taxworkshop folder in your directory.

**Tip:** You might want to type your password in to a text document on your computer, so you can copy and paste it if you need to log in later.

Now, download this github repository:

```bash
git clone https://github.com/shahnidhi/m3-taxonomy-workshop.git
```

We will run all our exercises in an interactive job. To start an interactive job run the following command
```bash
srun --pty --partition class --account=class --qos class --mem=8g --time=04:00:00 bash
```

We've also created an alias for this, so you can just type `int` and that should work too. As a reminder you can close your interactive job by typing `exit`.

**If you are using your UMIACS account** you need to run the following to set up your environment when you first log in:

```bash
/fs/m3taxworkshop/.interal_umiacs_people_run_this.sh
exit
```

This will set up your environement and log you out. Log back in and you should be good to go.

## Tools we will run today:

Here are the tools we'll be working with today:

1. [RDP](#running-rdp-classifier)
2. [Kraken](#running-kraken)

## Datasets 

### HMP Stool Samples 

We selected a subset of 10 stool samples sequenced by the Human Microbiome Project using V1V3 marker gene methods. <br />

Location - ```~/m3taxworkshop/data/1-datasets/hmp/```

### Halite (Salt Rock) from the Atacama Desert, Chile 

This whole metagenomic shotgun sequencing study characterized the microbial communities in halite modules in the Atacama Desert, and in particular their response to unusual rainfall in August 2015. Samples are from four time points (pre-rain in Sep 2014 and June 2015 and post rain in Feb 2016 and Feb 2017) with five replicates for each time point. <br />
Special thanks to Gherman Uritskiy and Jocelyne DiRuggiero for providing the data! <br />
The preprint including this data can be found at https://www.biorxiv.org/content/early/2018/10/13/442525 <br />

Location - ```~/m3taxworkshop/data/1-datasets/atacama_halite_timeline/```

## Machine Learning Approaches 

### Running RDP classifier

For this part, we are going to use the hmp stool dataset and the RDP classifier.

To run the RDP classifier, you input you representative set sequences and it outputs a text file with taxonomic assignments and confidence values for each sequence. 

You can run RDP using the following command:

```bash

java -jar /fs/m3taxworkshop/software/RDPTools/classifier.jar classify /fs/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_filtered_final.fna -o hmp_stool_rdp.txt

```

OR if you don't want to type the full command, you can use our wrapper script:

```bash

cd ~/m3-taxonomy-workshop/run_rdp

assign_taxonomy.sh -h

assign_taxonomy.sh /fs/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_filtered_final.fna -o hmp_stool_rdp.txt

```

RDP is also part of QIIME's script assign_taxonomy.py. If you have QIIME installed, the command would look like this to obtain classifications with a 80% confidence:

```bash
assign_taxonomy.py -m rdp -o rdp_taxonomy_stool_v1v3 -i /fs/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_filtered_final.fna -c 0.8
```

We have also provided a python script that extracts some information about our taxonomic assignments.

```bash
python format_rdp_output_to_csv.py -t hmp_stool_rdp.txt -o hmp_stool_rdp -q /fs/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_filtered_final.fna 
```


## Fast Metagenomic Profiling Methods

### Running Kraken

Kraken is a k-mer-based taxonomic classification tool. Here's a study comparing the performance of Kraken with other k-mer-based (and some non-k-mer-based) taxonomic classifiers: https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-017-1299-7.

For the exercise, make sure you're interactively logged in and then change to your home folder: `cd ~`

Let's search with Kraken1:

```bash
kraken --db m3taxworkshop/databases/kraken/minikraken1_8GB/ --threads 4 m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_filtered_final.fna > hmp.kraken
```

Next we'll create the Kraken report:

```bash
kraken-report --db m3taxworkshop/databases/kraken/minikraken1_8GB/ hmp.kraken > hmp.kreport
head -n5 hmp.kreport
```

#### Running Kraken2

Next we'll try the updated version, Kraken2:

```bash
kraken2 --db m3taxworkshop/databases/kraken/minikraken2_8GB/ --threads 4 --report hmp.kreport2 m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_filtered_final.fna > hmp.kraken2
```

Let's check the output now:

```bash
head -n5 hmp.kreport2
```

## Database searching - sequence alignment based approaches

### Running BLAST

```bash
cd ~/m3-taxonomy-workshop/run_blast/blast_toy_example
ls
makeblastdb -in database.fasta -out database.fasta -dbtype nucl
blastn -h
blastn -query query.fasta  -db database.fasta  -outfmt " 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen " -out blast.out -num_threads 4
```
Check this [site](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) to learn more about BLAST tabulated output format headers. <br /> 
### Running outlier detection pipeline

We need python3 environment and python packages such as scipy, networkx, python-louvain. We have created a virtual environment with all these installed, you just have to source it

```bash
export PYTHONPATH=''
source ~/m3taxworkshop/software/outlier_env/bin/activate
```

We are going to work on HMP stool dataset here, and use SILVA v.128 database in this example. <br />

- Query sequences: ```~/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_10_filtered.fna```
- Database sequences: ```~/m3taxworkshop/databases/silva/~/m3taxworkshop/databases/silva/SILVA_132_SSURef_Nr99_tax_silva_subset.fasta```
- Database taxonomy: ```~/m3taxworkshop/databases/silva/silva_nr_99_subset_taxonomy.tsv```

Staging BLAST output for hmp stool dataset

```bash
cd ~/m3-taxonomy-workshop/run_blast/hmp_example
```
You can copy the BLAST output to your current directory or run the blastn command (takes about 5 mins).
```bash
cp ~/m3taxworkshop/previous_run/blast/stool_blast.out.gz .
gunzip stool_blast.out.gz 
```
```bash
blastn -query ~/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_10_filtered.fna -db ~/m3taxworkshop/databases/silva/~/m3taxworkshop/databases/silva/SILVA_132_SSURef_Nr99_tax_silva_subset.fasta -outfmt " 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen " -out stool_blast.out -num_threads 4
```


Check outlier detection pipeline options and run on stool sample dataset

```bash
run_pipeline.py -h
run_pipeline.py -q ~/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_10_filtered.fna -b stool_blast.out -t ~/m3taxworkshop/databases/silva/silva_nr_99_subset_taxonomy.tsv -o stool_blast_outlier
```

Check specifically these files in the output folder:

1. results_outliers.txt - a tsv with reads and relevant DB sequences
2. results_partition_map_FINAL.txt - partition number to DB sequence mapping
3. results_read_to_partition_assignment.txt - partition assignment for reads
4. consensus_taxonomy_based_on_outliers.txt - LCA of outliers

Computing the number of reads classified at each taxonomic rank
```
python ~/m3-taxonomy-workshop/utils/format_output_to_csv.py \
    -q ~/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_10_filtered.fna \
    -t stool_blast_outlier/consensus_taxonomy_based_on_outliers.txt
    -o FINAL_stool_blast_outlier
python ~/m3-taxonomy-workshop/utils/format_output_to_csv.py \
    -q ~/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_10_filtered.fna \
    -t stool_blast_outlier/consensus_taxonomy_based_on_partition.txt
    -o FINAL_stool_blast_partition
```
and examining the read count for genus-level classification
```
cat FINAL_stool_blast_outlier_species.csv
cat FINAL_stool_blast_partition_species.csv
```
## Phylogenetic methods

### Running TIPP pipeline

If you are still in python3 environment, run ```deactivate```

Because TIPP needs python2, switch to python2 environment by running 

```bash
source ~/m3taxworkshop/test_user_profile.sh
source ~/m3taxworkshop/software/tipp_env/bin/activate
```

Inputs

- Set of query sequences, i.e., fragments/reads of unknown origin
- Reference alignment and tree or taxonomy
#### Run on small example
```bash
cd ~/m3-taxonomy-workshop/run_tipp/tipp_small_example/
run_tipp.py -a /fs/m3taxworkshop/software/tipp/refpkg/RDP_2016_Bacteria.refpkg/pasta.fasta \
            -t /fs/m3taxworkshop/software/tipp/refpkg/RDP_2016_Bacteria.refpkg/pasta.taxonomy \
            -r /fs/m3taxworkshop/software/tipp/refpkg/RDP_2016_Bacteria.refpkg/RAxML_info.taxonomy \
            -tx /fs/m3taxworkshop/software/tipp/refpkg/RDP_2016_Bacteria.refpkg/taxonomy.table \
            -txm /fs/m3taxworkshop/software/tipp/refpkg/RDP_2016_Bacteria.refpkg/species.mapping \
            -A 100 \
            -P 1000 \
            -at 0.95 \
            -pt 0.95 \
            -f small_example.fasta \
            -o TIPP_RDP_small_example \
            --tempdir tmp \
            --cpu 2
```
This will take 5-6 minutes to finish. In the meantime, let's break down the command. The first five options specify files included for the reference
+ `-a [reference multiple sequence alignment -- fasta format]`
+ `-t [reference taxonomy -- newick format]`
+ `-r [reference tree model parameters -- RAxML info file]`
+ `-tx [mapping taxonomic id to taxonomy information -- csv]`
+ `-txm [mapping sequence names to taxonomic IDs -- csv]`

The next two options specify the decomposition of the reference alignment and tree into subsets.
+ `-A [alignment subset size]`
+ `-P [placement subset size]`

TIPP was run with support thresholds of 0.95, which is the default. 

The next two options specify the input and output.
+ `-f [fragment file -- fasta]`
+ `-o [prefix of output files]`

To see all of the [TIPP options](https://github.com/ekmolloy/stamps-tutorial/blob/master/tipp-help.md), run
```
run_tipp.py -h
```
By now TIPP must have finished and written the following files
+ classification information -- csv
+ phylogenetic placement information -- json
+ alignment on both the reference and query sequences -- fasta

The classification file shows the support of classifying sequences at each taxonomic rank. Check out the support for each read classified at the species level
```
grep ",species," TIPP_RDP_small_example_classification.txt
```
Computing the number of reads classified at each taxonomic rank
```
python ../utils/restructure_tipp_classification.py \
    -i TIPP_RDP_small_example_classification.txt \
    -o FINAL_TIPP_small_example
```
and examining the read count for species-level classification
```
cat FINAL_TIPP_small_example_species.csv
```
What do read counts look like at the genus and family level? 

*Before moving on, repeat this portion of the tutorial running TIPP with a lower alignment/placement support threshold (e.g., 0.50). What do the support values look like for reads classified at the species level? How does the number of reads unclassified at the species level compare to TIPP run with an alignment/placement support threshold of 0.95?*

#### Run TIPP on HMP stool sample
```bash 
cd ~/m3-taxonomy-workshop/run_tipp/hmp_example
sh ../wrapper_tipp_script.sh ~/m3taxworkshop/data/1-datasets/hmp/stool_sample_subset_rep_set_10_filtered.fna FINAL_stool_hmp
```
If you don't want to run this, we have already provided output in ```out/``` folder. You can analyze the result files ```FINAL_stool_hmp*``` and what are the differences from the previous tools' outputs.

