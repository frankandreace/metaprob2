#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 24 Feb 2021
@author: francescoandreace
"""
# Import libraries
import logging
import argparse
import numpy as np
import scipy as sp
import time
import csv
from sknetwork.clustering import Louvain
import resource

#Impor local libraries
import community_utils as cu

# Sample command
# -------------------------------------------------------------------
# python ReadGraph_SGA.py    --gfa input gfa file
#                            --reads input fasta file
#                            --out /path/to/output_folder
#                            --info /path/to/initial_infos
#                            --opt optimization parameter for modularity optimization
#                            --maxlen max length chained utgs
#                            --lout include the unassembled reads in the output fasta 
#                            --rangeutg set min and max length of unitigs that are passed to the last step
# -------------------------------------------------------------------


# Setup logger
# -----------------------

logger = logging.getLogger('community_detection_clear')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
logger.addHandler(consoleHeader)

start_time = time.time()

# Setup argument parser
# -----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--gfa", required=True, help="path to miniasm output file (gfa format)")
ap.add_argument("--reads", required=True, help="path to the reads file (fasta or fastq")
ap.add_argument("--out", required=True, help="output file path")
ap.add_argument("--info", required=True, help="path to the initial info csv")
ap.add_argument("--opt", required=False, type=float, default=0.001,
                help="path to the output fasta")
ap.add_argument("--maxlen", required=False, type=int, default=20000,
                help="maximum length of concatenated unitigs in the output fasta file. [default: 20000]")
ap.add_argument("--lout", default=False, action='store_true',
                help="call it to skip the reads left out. [default: False]")
ap.add_argument("--rangeutg", required=False, type=int, nargs=2, default=[0,-1],
                help="set to min and max length of unitigs that are passed to the last step. [default: 0,unlimited]")
args = vars(ap.parse_args())

# Getting the info from arguments

gfa_path = args["gfa"]
reads_path = args["reads"]
out_path = args["out"]
max_length_grouped = args["maxlen"]
opt_par = args["opt"]
info_p = args["info"]
skip_left_out = args["lout"]
range_utgs = args["rangeutg"]

out_csv_path = out_path+ '.csv'
out_reads_path = out_path+ '.fasta'

# Setup output path for log file
# ---------------------------------------------------

fileHandler = logging.FileHandler(out_path + "community_detection.log")
fileHandler.setLevel(logging.INFO)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)
logger.info(
    "Community detection on Unitigs graph using GFA file from Miniasm")


#max num reads
info_file= open(info_p,'r')
n_reads = int(info_file.read()) # MAX INDEX
info_file.close()
logger.info("Total number of reads: "+str(n_reads))

#### Getting infos from GFA

t1 = time.time()
utgs, grps, adj , ww = cu.read_gfa(gfa_path)
logger.info("GFA into Network info in "+str(time.time()-t1)+"seconds.")
#number of unique reads in all utgs
uniq = cu.get_unique_elements(grps)
n_unique = len(uniq)
logger.info("Unique number of reads in GFA: "+str(n_unique))
logger.info('Number of utgs: '+str(len(grps)))
logger.info('Number reads left out'+str(n_reads-n_unique))
# creo matrice sparse per clustering
### IN PIU
with open(out_csv_path, 'w', newline="") as f:
    writer = csv.writer(f)
    writer.writerows(adj)
### IN PIU

leng = len(utgs) 
network = sp.sparse.csr_matrix((adj[2],(adj[0],adj[1])),shape=(leng,leng))
#print(network.get_shape())

# modularity opt for community detection
logger.info('Louvain alg with optimization level = '+str(opt_par))
opt_lev = opt_par # 0.001
louvain = Louvain(random_state=0,tol_aggregation=opt_lev,tol_optimization=opt_lev)
out = louvain.fit_transform(network)

clusters,n_out = cu.get_clusters(out)
n_groups = len(clusters)
logger.info('Number of clusters: '+str(n_groups))

### REPRESENTATIVES CHOICE ###

#evaluating degree of each utg
deg = np.zeros(leng,dtype=np.uint32)
for i in range(len(adj[0])):
    deg[adj[0][i]] += adj[2][i]
    deg[adj[1][i]] += adj[2][i]

#create representatives based on deg
#max_length_grouped = 10000
grpd_utgs = cu.get_grpd_utgs(utgs,deg,clusters,max_length_grouped)

#select unitigs in specified range
if range_utgs[1] == -1:
    range_utgs[1] = max_length_grouped*2

grpd_utgs,clusters = cu.get_range(range_utgs,clusters,grpd_utgs)


#csv groups of utgs
csv_cl = cu.get_outputcsv(clusters,grps)

#eliminating redundancies
check = np.zeros(n_reads,dtype=np.bool)
for grp in csv_cl:
    to_rmv = []
    for el in grp:
        if check[el] == False:
            check[el] = True
        else:
            to_rmv.append(el)
    for el in to_rmv:
        grp.remove(el)

# + output
logger.info('Writing csv...')

with open(out_csv_path, 'w', newline="") as f:
    writer = csv.writer(f)
    writer.writerows(csv_cl)


### HANDLING READS LEFT OUT & OUTPUT READS ###
leftOut = cu.get_leftOut(grps,n_reads)
n_left_out = len(leftOut)

#create new fasta + utgs
logger.info('Writing fasta...')
new_fasta = ''
new_fasta += cu.add_grps(grpd_utgs)

# add left out
if skip_left_out == False:
    logger.info('Adding external reads')
    reads_fq = open(reads_path,'r').read().split('\n')
    div = cu.get_div(reads_path)
    logger.info(div)
    new_fasta += cu.add_reads(leftOut,reads_fq,div,n_reads)

# write 
out_file = open(out_reads_path ,"w")
out_file.write(new_fasta)
out_file.close()

max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

logger.info('Time elapsed: '+str(time.time()-start_time))
logger.info('Max memory used (in MB): '+str(round(max_mem / 1000.0, 3)))
logger.info('DONE')