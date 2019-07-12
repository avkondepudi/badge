#!/usr/bin/env python3

'''
Get differentially expressed genes.
'''

import argparse

from badge import Badge
from utils import str2bool

ap = argparse.ArgumentParser(description="Get differentially expressed genes.")

ap.add_argument("input", type=str, help="input file")
ap.add_argument("ctypes", type=str, help="cell types", nargs="*")

ap.add_argument("-g", "--graph", default=True, type=str2bool, help="whether to graph (default: True)")

args = vars(ap.parse_args())

b = Badge(args.get("input"))

cg = b.candidate_genes(args.get("ctypes"))
print(cg)

if args.get("graph"):
	for ctype in args.get("ctypes"): b.graph(cg[ctype])
