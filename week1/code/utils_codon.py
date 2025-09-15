import os
from typing import List, Tuple

def read_fasta(path: str, name: str) -> List[str]:
    data: List[str] = []
    filename = os.path.join(path, name)
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] == ">":
                continue
            data.append(line)
    if data:
        print(name, len(data), len(data[0]))
    else:
        print(name, 0, 0)
    return data

def read_data(path: str) -> Tuple[List[str], List[str], List[str]]:
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1  = read_fasta(path, "long.fasta")
    return short1, short2, long1
