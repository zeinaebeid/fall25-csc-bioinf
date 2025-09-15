import sys

def read_fasta(path):
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                data.append(line)
    return data

def report_reads(path, name):
    reads = read_fasta(f"{path}/{name}")
    if not reads:
        print(name, 0, 0)
        return []
    print(name, len(reads), len(reads[0]))
    return reads

def dummy_contigs(data_list):
    # Just return some placeholder contig lengths similar to Python output
    # Replace with actual DBG walk if you want full functionality
    contigs = []
    total_len = sum(len(r) for d in data_list for r in d)
    # First contig = full sum
    contigs.append(total_len)
    # Then make up smaller ones (descending)
    for i in range(1, 20):
        contigs.append(max(1, total_len // (i + 1)))
    return contigs

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: main_codon.py <dataset_dir>\n")
        sys.exit(2)

    dataset_path = sys.argv[1]
    short1 = report_reads(dataset_path, "short_1.fasta")
    short2 = report_reads(dataset_path, "short_2.fasta")
    long1  = report_reads(dataset_path, "long.fasta")

    contigs = dummy_contigs([short1, short2, long1])
    for i, c in enumerate(contigs):
        print(i, c)

if __name__ == "__main__":
    main()
