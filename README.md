# Winter 2019 Mid-Atlantic Microbiome Meetup
# Taxonomic Identification Workshop

## Logging onto the server
There is username and password at the back of your name tag.
```bash
ssh USERNAME@openclass.umiacs.umd.edu
```
Enter password. Once you are logged in, you should be able to see m3taxworkshop folder in your directory. <br />
In order to load all of the pre-installed software, please run:

```bash
# Load all software modules
source m3taxworkshop/test_user_profile.sh
```
Now, download this github repository 
```bash
git clone https://github.com/shahnidhi/m3-taxonomy-workshop.git
```
## Datasets 

### HMP Stool Samples 

We selected ten stool samples sequenced by the Human Microbiome Project using both V1V3 marker gene methods and whole metagenomic shotgun sequencing. <br />
Location - /fs/m3taxworkshop/data/1-datasets/hmp/

### Halite (Salt Rock) from the Atacama Desert, Chile 

This whole metagenomic shotgun sequencing study characterized the microbial communities in halite modules in the Atacama Desert, and in particular their response to unusual rainfall in August 2015. Samples are from four time points (pre-rain in Sep 2014 and June 2015 and post rain in Feb 2016 and Feb 2017) with five replicates for each time point. <br />
Special thanks to Gherman Uritskiy and Jocelyne DiRuggiero for providing the data! <br />
The preprint including this data can be found at https://www.biorxiv.org/content/early/2018/10/13/442525
Location - /fs/m3taxworkshop/data/1-datasets/atacama_halite_timeline/

## Machine Learning Approaches 
### Running RDP classifier
For this part, we are going to use stool dataset and greengenes 16S rRNA gene database. 
```bash
cd ~/m3-taxonomy-workshop/run_rdp
assign_taxonomy.py -h

assign_taxonomy.py -m rdp -o rdp_taxonomy_stool_v1v3 -i /fs/m3taxworkshop/data/1-datasets/hmp/
```
