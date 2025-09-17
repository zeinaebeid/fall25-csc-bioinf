from typing import List

def join(a: str, b: str) -> str:
    return (a if a.endswith("/") else a + "/") + b

def read_fasta(dirpath: str, name: str) -> List[str]:
    fp = join(dirpath, name)
    seqs: List[str] = []
    with open(fp, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s[0] != ">":
                seqs.append(s.upper()) 

    if len(seqs) > 0:
        print(name, len(seqs), len(seqs[0]))
    else:
        print(name, 0, 0)
    return seqs

def read_data(dirpath: str):
    short1 = read_fasta(dirpath, "short_1.fasta")
    short2 = read_fasta(dirpath, "short_2.fasta")
    long1  = read_fasta(dirpath, "long.fasta")
    return short1, short2, long1
