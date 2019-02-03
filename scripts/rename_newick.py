# coding: utf-8

import argparse
import pandas as pd

parser = argparse.ArgumentParser(add_help=True)

parser.add_argument('-i', '--input')
parser.add_argument('-t', '--table')
parser.add_argument('-o', '--output')

args = parser.parse_args()

with open(args.input) as f:
    tree = f.readlines()[0]

df = pd.read_csv(args.table, header=None)
df.columns = ["before", "after"]

for index, row in df.iterrows():
    tree = tree.replace(row["before"], row["after"])

with open(args.output, "w") as f:
    f.write(tree)