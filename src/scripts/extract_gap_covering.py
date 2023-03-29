#!/usr/bin/env python

import sys
import my_bio

if len(sys.argv) != 3:
    print("usage:", sys.argv[0], "q_name_sorted.paf ref.fa", file=sys.stderr)
    sys.exit(1)

ref_seq = dict()
with open(sys.argv[2], "r") as fin:
    for name, seq in my_bio.read_fasta(fin):
        ref_seq[name] = seq.upper()

with open(sys.argv[1], "r") as fin:
    pos_range_dict = dict()
    cur_q_name = ""
    for line in fin:
        f = line.rstrip("\n").split("\t")
        if f[0] != cur_q_name:
            for name, pos_range in pos_range_dict.items():
                if "N" in ref_seq[name][pos_range[0]:pos_range[1]]:
                    print(cur_q_name)
            pos_range_dict = dict()
            cur_q_name = f[0]

        start = int(f[7])
        end = int(f[8])
        if not f[5] in pos_range_dict:
            pos_range_dict[f[5]] = [start, end]
        else:
            pos_range = pos_range_dict[f[5]]
            if start < pos_range[0]:
                pos_range[0] = start
            if end > pos_range[1]:
                pos_range[1] = end

    for name, pos_range in pos_range_dict.items():
        if "N" in ref_seq[name][pos_range[0]:pos_range[1]]:
            print(cur_q_name)
