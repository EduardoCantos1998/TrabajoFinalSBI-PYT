#! /usr/bin/env python

import os
import sys
import warnings

pdb_list = [f"./PDB/{i[0:5]}" for i in os.listdir("./PDB/")]

for i in pdb_list:
    os.mkdir(i)

