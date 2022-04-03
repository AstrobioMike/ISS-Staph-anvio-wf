#!/usr/bin/env python

import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('in_tab')
parser.add_argument('out_tab')

args = parser.parse_args()

tab = pd.read_csv(args.in_tab, sep = "\t")
mat = tab.pivot(index = "item", columns = "layer", values = 'value')
mat.to_csv(args.out_tab, sep = "\t")

os.remove(args.in_tab)
