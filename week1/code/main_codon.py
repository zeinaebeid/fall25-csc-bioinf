from typing import List
import sys
from dbg_codon import DBG
from utils_codon import read_data, join  


def main(argv: List[str]) -> int:
    if len(argv) < 2:
        print("Usage: main_codon.py <dataset_dir>\n"
              "Example: codon run -release main_codon.py week1/data/data1")
        return 2

    dataset_dir = argv[1]

    short1, short2, long1 = read_data(dataset_dir)

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    out_fp = join(dataset_dir, "contig.fasta")
    with open(out_fp, "w") as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if not c:
                break
            print(i, len(c))
            f.write(f">contig_{i}\n")
            f.write(c + "\n")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
