from dbg_codon import DBG
from utils_codon import read_data
import sys


def usage():
    sys.stderr.write("usage: main_codon.py <dataset_dir>\n")
    sys.exit(2)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        usage()

    # Codon does not support os.path.join → use simple string concat
    dataset_path = "./" + sys.argv[1]

    short1, short2, long1 = read_data(dataset_path)

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    output_file = dataset_path + "/contigs_codon.fasta"
    with open(output_file, "w") as f:
        for i in range(20):
            contig = dbg.get_longest_contig()
            if contig == "":   # use empty string instead of None
                break
            print(i, len(contig))
            f.write(f">contig_{i}\n{contig}\n")

    # write status to stderr (so stdout stays clean if redirected)
    print(f"[Codon] Wrote contigs to {output_file}", file=sys.stderr)
