# remove_duplicates.py
# input: fasta file
# output:
#   - fasta file with no duplicate sequences
#   - table of the descriptions and protein accessions of the duplicates ("mapper")

import argparse
import pandas as pd
import subprocess
from src import config
from Bio import SeqIO
from Bio import Entrez

def get_bracketed_organism_name(description):
    # get the name in brackets from a fasta description string (i.e. [Pseudomonas cremoris])
    # returns 'NONE' for strings without names 
    try:
        return(description[description.find("[")+1:description.find("]")])
    except:
        return('NONE')

def get_sequence_dict(fasta):
    # returns a dictionary in which the key is the sequence 
    # the items are dictionaries with the accession, organism, and description of associated sequences

    sequence_dict = {}

    for index, record in enumerate(SeqIO.parse(fasta, "fasta")):
        if str(record.seq) not in sequence_dict.keys():
            sequence_dict[str(record.seq)] = []
        d = sequence_dict[str(record.seq)]
        acc = record.id.split(':')[0]
        organism = get_bracketed_organism_name(record.description)
        description = record.description
        d.append({'acc': acc, 'organism': organism, 'description':description})

    return(sequence_dict)

def find_dupes(fasta):
    # finds duplicate sequences and returns dupe IDs to drop and a dataframe of the kept sequences

    sequence_dict = get_sequence_dict(fasta)
    # two+ sequences that are identical and form the same organism

    dupe_ids_to_drop = [] # the ids of the dupe sequences to drop
    # one of the dupes is kept, a file will be output that maps kept sequences to the dupes

    kept_dupes = {'acc': [], 'organism': [], 'dupe_accs': [], 'dupe_organisms': []}

    for s in sequence_dict.values():
        # just keep the first in the list
        kept_dupes['acc'].append(s[0]['acc'])
        kept_dupes['organism'].append(s[0]['organism'])

        dropped_dupe_ids = []
        dropped_dupe_organisms = []
        for r in s[1:]: # skip the first one
            dupe_ids_to_drop.append(r['acc'])
            dropped_dupe_ids.append(r['acc'])
            dropped_dupe_organisms.append(r['organism'])
        kept_dupes['dupe_accs'].append(';'.join(dropped_dupe_ids))
        kept_dupes['dupe_organisms'].append(';'.join(dropped_dupe_organisms))

    kept_dupes_df = pd.DataFrame.from_dict(kept_dupes)
    
    return(dupe_ids_to_drop, kept_dupes_df)


def get_reduced_fasta(fasta, dupe_ids_to_drop):
    # writes a fasta file with the duplicate sequences removed

    to_write = []
    for index, record in enumerate(SeqIO.parse(fasta, "fasta")):
        if record.id.split(':')[0] not in dupe_ids_to_drop:
            to_write.append(record)

    print('filtered out ' + str(len(dupe_ids_to_drop)) + ' duplicates')
    print('remaining: ' + str(len(to_write)))
    
    return(to_write)

def write_output(reduced_fasta, fasta, dupe_ids_to_drop, kept_dupes_df):
    SeqIO.write(reduced_fasta, fasta.split('.')[0] + '_no_dupes.fasta', "fasta")
    kept_dupes_df.to_csv(fasta.split('.')[0] + '_dupes_key.csv')
    return()

def read_blast_table(table):
    table_df = pd.read_csv(table,
        usecols=[0, 1, 2, 6,7, 10, 11],
        names=['query_acc', 'subject_acc', 'pctid', 'qstart', 'qend', 'evalue', 'bit_score'])
    return(table_df)

def main(fasta, output_directory, hits_table):
    dupe_ids_to_drop, kept_dupes_df = find_dupes(fasta)
    reduced_fasta = get_reduced_fasta(fasta, dupe_ids_to_drop)
    if hits_table:
        hits_df = read_blast_table(hits_table)
        hits_df = hits_df[hits_df['subject_acc'].isin(kept_dupes_df['acc'])]

    os.chdir(output_directory)

    

    write_output(reduced_fasta, fasta, dupe_ids_to_drop, kept_dupes_df)
    if hits_table:
        hits_df.to_csv(fasta.split('.')[0] + '_no_dupes_hits.csv')

    return()

if __name__ == "__main__":

    print("removing duplicate sequences")  
    import sys
    import argparse
    import os
    import re
    import time
    from src import config
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("input_directory", nargs='?', default=config.raw, type=str, help="name of directory with raw data")
    parser.add_argument("output_directory", nargs='?', default=config.raw, type=str, help="name of directory with processed data")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    parser.add_argument("-t", "--hits_table", action="store", default=None, help='give a hits table (optional)')
    
    args = parser.parse_args()

    try:
        os.chdir(args.input_directory)
    except:
        print ("didn't change dir")

    main(args.fasta, args.output_directory, args.hits_table)