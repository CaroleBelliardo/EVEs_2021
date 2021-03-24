#!/usr/bin/env python
# coding: utf-8

# # Taxo IMGM

import csv
import io
import os
import sys
import pandas as pd
import multiprocessing as mp
from time import process_time
import argparse

parser = argparse.ArgumentParser(description='Need PYTHON 3.8 version with packages ...')

parser.add_argument('-m', '--mg', type=str, help='correspondance Table\n',
                    default='/lerins/hub/DB/Metagenomics/3_Concat_files/0_TAXO_seqidToTaxid_v4/dataOri_adrien_out_sortUniq')
parser.add_argument('-n', '--nodes', type=str, help='nodes\n', default='/lerins/hub/DB/TAXONOMY/online/nodes.dmp')
parser.add_argument('-l', '--lineage', type=str, help='lineage\n',
                    default='/lerins/hub/DB/TAXONOMY/online/taxonomy_lineage.txt')
parser.add_argument('-L', '--logs', type=str, help='logs\n', default='/lerins/hub/projects/25_tmp/logs/IMGM_taxo2')
parser.add_argument('-o', '--out', type=str, help='output\n',
                    default='/lerins/hub/DB/Metagenomics/3_Concat_files/0_TAXO_seqidToTaxid_v4/dataOri_adrien_out_sortUniq_taxid')
args = parser.parse_args()


t1_start = process_time()
lineages = pd.read_table(args.lineage, sep='\t', names=['taxid', 'lineage', 'complTax'])
t1_stop = process_time()
print(f'time import lineage : {t1_stop - t1_start}')


t1_start = process_time()
nodes = pd.read_csv(args.nodes, sep='\t', header=None, dtype=str)[[0, 4]]
nodes = nodes.rename(columns={0: 'taxid', 4: 'rank'})
t1_stop = process_time()
print(f'time import nodes: {t1_stop - t1_start}')


def time(func):
    t1_start = process_time()
    func
    t1_stop = process_time()
    print(f'time : {t1_stop - t1_start}')


def affectTaxoDef(lineage):
    taxid = '0'
    rank = ''
    if lineage == 'Bacteria': taxid = '2'
    if lineage == 'Eukaryota': taxid = '2759'
    if lineage == 'Viruses': taxid = '10239'
    if lineage == 'Archaea': taxid = '2157'
    if taxid != '0':
        rank = 'superkingdom'
    return taxid, rank


def affectTaxoDef2(lineage):
    taxid = '0'
    rank = ''
    if 'bacteria' in lineage.lower(): taxid = '2'
    if 'eukaryota' in lineage.lower(): taxid = '2759'
    if 'virus' in lineage.lower(): taxid = '10239'
    if 'archaea' in lineage.lower(): taxid = '2157'
    if taxid != '0':
        rank = 'superkingdom'
    return taxid, rank


def progress(iteration, steps, max_value, no_limit=False):
    if int(iteration) == max_value:
        if no_limit == True:
            sys.stdout.write('\r')
            print("[x] \t%d%%" % (100), end='\r')
        else:
            sys.stdout.write('\r')
            print("[x] \t%d%%" % (100))
    elif int(iteration) % steps == 0:
        sys.stdout.write('\r')
        print("[x] \t%d%%" % (float(int(iteration) / int(max_value)) * 100), end='\r')
        sys.stdout.flush()
    else:
        pass


def parrallelize(func, jobL):
    pool = mp.Pool(os.cpu_count())
    for i, _ in enumerate(pool.imap_unordered(func, jobL), 1):
        progress(i, 1, len(jobL))


def affectTaxo(lineage):
    global lineages
    taxid = '0'
    rank = ''
    lineage=str(lineage)
    l_ = lineage

    if ':' in lineage:
        l_list = l_.split(':')
        l_ = max(l_list, key=len)

    taxid, rank = affectTaxoDef(l_)

    if taxid != '0':
        return (taxid)
    else:
        l_ = l_ + ';'
        t = lineages.taxid[lineages.lineage.str.endswith(l_)]
        if len(t) > 0: taxid = t.values[0]

        while (len(l_.rsplit(' ', 1)) > 1) and (taxid == '0'):
            l_ = l_.rsplit(' ', 1)[0]
            l_ = l_ + ';'
            t = lineages.taxid[lineages.lineage.str.endswith(l_)]
            if len(t) > 0: taxid = t.values[0]
        if taxid == '0':
            taxid, rank = affectTaxoDef2(lineage)
    return (taxid)


## FUNCTION to parrallize
## lister fichier repo
ld=os.listdir(args.mg)

def main(file):
    global lineages
    global nodes
    global args
    mgF=f'{args.mg}/{file}'
    out=f'{args.out}/{file}'
    logF=f'{args.logs}/{file}'
    taxid = '0'
    mid=file.split('.')[0]    

    t1_start = process_time()
    mg = pd.read_csv(mgF, sep='\t', names=['contig', 'species', 'taxid', 'rank'], dtype=str)
    t1_stop = process_time()
    print(f'time import metag: {t1_stop - t1_start}')

    for s in mg[mg['taxid'].isnull()].species.unique():
        taxid = affectTaxo(s)
        mg.loc[mg['species'] == s, 'taxid'] = taxid
        n = nodes[nodes.taxid == str(taxid)]
        if not n.empty:
            rank = n['rank'].values[0]
            mg.loc[mg['species'] == s, 'rank'] = rank

    t1_start = process_time()
    with open(logF, 'a') as log:
        log.write(f'{len(mg.contig)}=')
        log.write(f'{len(mg.contig.unique())}=')
    t1_stop = process_time()
    print(f'log len mg: {t1_stop - t1_start}')

    t1_start = process_time()
    mg = mg.drop_duplicates()
    t1_stop = process_time()
    print(f'dropping duplicated rows: {t1_stop - t1_start}')

    t1_start = process_time()
    with open(logF, 'a') as log:
        log.write(f'dropped tab = {len(mg.contig.unique())}=')
    t1_stop = process_time()
    print(f'log len mg: {t1_stop - t1_start}')

    t1_start = process_time()
    mg['newName'] = mg['contig'].astype(str) + '_' + str(mid) + '__mg__IMGM'
    mg['genId'] = mid
    mg.to_csv(f'{out}.csv', index=False, sep='\t', header=False)
    t1_stop = process_time()
    print(f'printing complet tab: {t1_stop - t1_start}')

    t1_start = process_time()
    mg = mg.drop(columns=['species', 'rank'])
    mgf= { 'seqId1': mg['newName'], 'seqId2':  mg['newName'], 'taxid':  mg['taxid'] }
    mg = pd.DataFrame(mgf)    
    #mg=mg[['newName', 'newName', 'taxid']
    mg.to_csv(f'{out}_taxo.txt', index=False, header=False, sep='\t')
    t1_stop = process_time()
    print(f'formating 2 col and printing: {t1_stop - t1_start}')

    with open(logF, 'a') as log:
        log.write(f'{len(mg.seqId1)}=')
        log.write(f'{len(mg.seqId1.unique())}=')
        if len(mg.seqId1) != len(mg.seqId1.unique()) : log.write(f'duplicated contig in : {file}=')


pool = mp.Pool(70)
parrallelize(main,ld)

