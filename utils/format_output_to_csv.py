import os
import argparse
import subprocess
from subprocess import Popen, PIPE
from collections import defaultdict
def main():
    parser = argparse.ArgumentParser(description="Converts tab seperated taxonomic assignment files to csv with counts at each taxonomic level")
    parser.add_argument("-q","--query_file", help="A fasta file of query sequences",required=True)
    parser.add_argument("-t","--taxa_assign_file", help="Tab seperated taxonomic assignment file",required=True)
    parser.add_argument("-o","--output_prefix", help="Output prefix", default="assignment" , required=False)
    
    args = parser.parse_args()
    total_num_reads = 0
    with open(args.query_file) as f:
        for line in f:
            if line.startswith('>'):
                total_num_reads += 1
    rank_to_name = {0:'kingdom', 1:'phylum', 2:'class', 3:'order', 4:'family', 5:'genus', 6:'species'}
    counts = {}
    total_counts = defaultdict(int)
    for i in range(0,7):
        counts[i] = defaultdict(int)
    #Read the file 
    fw = open(args.output_prefix+'_assignment.csv', 'w')
    with open(args.taxa_assign_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            val = line.strip().split('\t')
            linetoprint = [val[0]]
            taxonomy = val[1].split(';')
            counter = 0
            for counter in range(0, len(taxonomy)):
                if taxonomy[counter] =='' or taxonomy[counter].lower()=='unclassified' or taxonomy[counter].lower()=='unassigned':
                    linetoprint.append('Unclassified')
                    counts[counter]['Unclassified'] += 1
                else:
                    linetoprint.append(taxonomy[counter])
                    counts[counter][taxonomy[counter]] += 1
                    total_counts[counter] += 1
            
            while counter <= 6:
                linetoprint.append('Unclassified')
                counts[counter]['Unclassified'] += 1
                counter += 1
            fw.write(','.join(linetoprint)+'\n')

    fw.close()
    fs = open(args.output_prefix+'_total_counts.csv', 'w')
    fs.write('total # reads,'+str(total_num_reads)+'\n')
    for rank in rank_to_name:
        fs.write('#reads classified at '+rank_to_name[rank]+' level ,'+str(total_counts[rank])+'\n')
        fw = open(args.output_prefix+'_'+rank_to_name[rank]+'.csv', 'w')
        for value in counts[rank]:
            fw.write(value+','+str(counts[rank][value])+'\n')
        fw.close()
    fs.close()
    
if __name__ == '__main__':
    main()