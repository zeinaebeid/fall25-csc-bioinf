from dbg_codon import DBG
from utils_codon import read_data
import sys

def usage():
    sys.stderr.write("usage: main_codon.py <dataset_dir>\n")
    sys.exit(2)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        usage()

    dataset_path = sys.argv[1]  # e.g., week1/data/data1
    short1, short2, long1 = read_data(dataset_path)

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    # Write FASTA contigs to STDOUT (so your evaluate.sh can just redirect)
    # Limit to 20 contigs exactly like your Python version.
    for i in range(20):
        c = dbg.get_longest_contig()
        if c == "":            # empty string = no more contigs
            break
        sys.stdout.write(f">contig_{i}\n{c}\n")
