# prune_taxon

# Removes all sequences from a certain taxon that are not present in an input list
# i.e. remove all Pseudomonas sequences that aren't in subsample_pseudomonas.txt

import argparse
import pandas as pd
from src import config
from Bio import SeqIO

def clean_id(record):
    # gets the clean id string
    # edit depending on your dirty strings
    return(record.id.split(':')[0])    

def in_taxon(record, taxon, taxon_df):
    sub_df = taxon_df[taxon_df['acc'] == clean_id(record)]
    for col in sub_df.columns:
        if len(sub_df[col].tolist()) > 0:
            if sub_df[col].tolist()[0] == taxon:
                return(True)
    else:
        return(False)

def read_keep_list(keep):
    keep_df = pd.read_table(keep, header=None)
    return(keep_df[0].tolist())
 

def main(taxon_table, fasta, taxon, keep):
    print('pruning ' + taxon + ' based on input taxon list')

    taxon_df = pd.read_csv(taxon_table)
    keep_list = read_keep_list(keep)
    to_keep = []
    dropped = 0

    for record in SeqIO.parse(fasta, "fasta"):
        if in_taxon(record, taxon, taxon_df):
            if (clean_id(record) in keep_list):
                # append seqs in taxon only if they are in the keep list
                to_keep.append(record)
            else:
                dropped +=1
        else:
            # append all other seqs not in taxon
            to_keep.append(record)

    keep_accs = [clean_id(record) for record in to_keep]
        
    
    print('Dropped ' + str(dropped) + ' sequences from ' + taxon)
    SeqIO.write(to_keep, fasta.split('.fasta')[0] + '_pruned.fasta', 'fasta')
    
    return()

if __name__ == "__main__":

    print("pruning taxa")  
    import sys
    import argparse
    import os
    import re
    import time
    from src import config
    parser = argparse.ArgumentParser(description="All")
   # parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("input_directory", nargs='?', default=config.data, type=str, help="name of directory with raw data")
    parser.add_argument("output_directory", nargs='?', default=config.data, type=str, help="name of directory with processed data")
    parser.add_argument("-k", "--keep", action = "store", default = False, help="list of sequence IDs to keep")
    parser.add_argument("-t", "--taxon", action = "store", default = False, help="taxid of taxon to trim within")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    parser.add_argument("--taxon_table", action="store", default = None, help='taxonomy table for the fasta')

    args = parser.parse_args()

    if not args.taxon_table:
        os.system('python src/data/fasta_taxonomy.py -f ' + args.input_directory + '/' + args.fasta)
        taxon_table = args.fasta.split('.fasta')[0] + '_taxonomy.csv'
    else:
        taxon_table = args.taxon_table
    
    #change dir if given
    try:
        os.chdir(args.input_directory)
    except:
        print ("didn't change dir")

    #run the thing
    main(taxon_table, args.fasta, args.taxon, args.keep)

    print("done")
