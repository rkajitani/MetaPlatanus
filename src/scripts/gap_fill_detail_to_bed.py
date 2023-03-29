#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print("usage:", sys.argv[0], "tgsgapcloser_gap_fill_detail.txt", file=sys.stderr)
    sys.exit(1)

with open(sys.argv[1], "r") as fin:
    seq_name = ""
    for ln in fin:
        if len(ln) == 0:
            continue
        if ln[0] == ">":
            seq_name = ln[1:].rstrip("\n")
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 3 and f[2] == "F" and len(seq_name) > 0:
            print(seq_name, int(f[0]) - 1, int(f[1]), sep="\t")
