import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        required=True, help="TIPP restructured classification .csv file")
    
    parser.add_argument("-o", "--output", type=str,
                        required=True, help="Prefix of output TSV file")

    args = parser.parse_args()
    fw = open(args.output+'_taxonomy_assignment.txt', 'w')
    fw.write('#read_name\ttaxonomy\tassigned_below_phylum?\n')
    with open(args.input) as f:
        for line in f:
            val = line.strip().split(',')
            fw.write(''.join([val[0], '\t', ';'.join(val[2:]), '\t', val[1], '\n']))
    fw.close()

if __name__ == '__main__':
    main()