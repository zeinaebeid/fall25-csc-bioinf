from dbg_codon import DBG   
from utils_codon import read_data
import sys

if __name__ == "__main__":
    argv = sys.argv
    dataset_path = "./" + argv[1]
    short1, short2, long1 = read_data(dataset_path)

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    output_file = dataset_path + "/contig.fasta"
    with open(output_file, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write(f">contig_{i}\n")
            f.write(c + '\n')
