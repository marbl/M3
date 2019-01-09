import argparse
from itertools import groupby
def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq
    fh.close()

def reverse_complement_fasta(ifas, ofas):

    fw = open(ofas, 'w')
    fiter = fasta_iter(ifas)
    #Inline function to get reverse complement of a DNA string
    revcompl_map = {'A':'T','C':'G','G':'C','T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a'}
    revcompl = lambda x: ''.join([revcompl_map[B] if B in revcompl_map else B for B in x][::-1])
    for ff in fiter:
        header = ff[0]
        rev_seq = revcompl(ff[1])
        fw.write(''.join(['>',header,'\n',rev_seq,'\n']))


    fw.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        required=True, help='Input fasta file')
    parser.add_argument("-o", "--output", type=str,
                        required=True, help="Output fasta file")

    args = parser.parse_args()
    reverse_complement_fasta(args.input, args.output)
