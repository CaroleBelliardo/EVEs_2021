#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys 

# arg 1 = input bed
# arg 2 = output bed with family column

bed = pd.read_csv(sys.argv[1], sep='\t', names=["Contig", "start", "stop","stranded","id%",'evalue', "lineage"], dtype=str)
sf = pd.read_csv('specieToFamilyTab.txt', sep='\t', names=["lineage", "Family"], dtype=str)
bed=bed.merge(sf,how='left',on='lineage')
bed.to_csv(sys.argv[2], index=False, sep='\t', header=False)
