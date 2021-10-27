# filter_blast

# locally filters blast results based on quality

import argparse
import pandas as pd
import subprocess
from src import config
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "ecutts@mit.edu"
Entrez.api_key = "057f36326473f362b4123bee57854f2c3208"



def get_unique_bidirectional_hits(table_df):
    # remove non-bidirectional hits remove multiple hits against the same subject (keeping only the 1st one)
    unique_subjects = table_df.subject_acc.unique()
    table_df = table_df[table_df.subject_acc.isin(unique_subjects)]
    table_df.drop_duplicates(subset=['query_acc', 'subject_acc'], inplace=True)
    print(table_df.head())
    return(table_df)

def ensure_self_aligning_queries(table_df):

    # calculate self bit scores
    self_bit_scores = table_df.loc[table_df.query_acc==table_df.subject_acc, 
        ['query_acc', 'bit_score']].copy()
    self_bit_scores.set_index('query_acc', inplace=True)
    self_bit_scores = self_bit_scores.bit_score

    # remove queries that don't align against themselves
    table_df = table_df[table_df.query_acc.isin(self_bit_scores.index)].copy()
    return(table_df, self_bit_scores)

def filter(table_df, evalue=10**-5, pct=0, qcover=0, bit_score=0):
    qlen = table_df[table_df.query_acc == table_df.subject_acc]['qend'].tolist()[0]
    table_df['qcover'] = (table_df.qend - table_df.qstart + 1)/qlen * 100 
    table_df = table_df[table_df.evalue < evalue]
    table_df = table_df[table_df.pctid > pct]
    table_df = table_df[table_df.qcover > qcover]
    table_df = table_df[table_df.bit_score > bit_score]
    return(table_df)


def read_blast_table(table):
    table_df = pd.read_csv(table,
        usecols=[0, 1, 2, 6,7, 10, 11],
        names=['query_acc', 'subject_acc', 'pctid', 'qstart', 'qend', 'evalue', 'bit_score'])
    orig_len = len(table_df)
    table_df = get_unique_bidirectional_hits(table_df)
    table_df, self_bit_scores = ensure_self_aligning_queries(table_df)

    return(table_df, orig_len)

def get_filtered_fasta_records(fasta, table_df):
    print(table_df.subject_acc.tolist()[0:4])

    filtered_records = []

    for record in SeqIO.parse(fasta, 'fasta'):
        print(record.id)
        if record.id in table_df.subject_acc.tolist():
            filtered_records.append(record)

    return(filtered_records)

def write_filtered_data(fasta, table, filtered_records, filtered_table, name):     
    fasta_base = fasta.split('.')[0]
    table_base = table.split('.')[0]
    
    filename = str(fasta_base + '_' + name + '_filtered.fasta')
    SeqIO.write(filtered_records, filename, "fasta")

    print('wrote ' + str(len(filtered_records)) + ' sequences to ' + fasta_base + '_' + name + '_filtered.fasta')

    filename = str(table_base + '_' + name + '_filtered.csv')
    filtered_table.to_csv(filename)
    
    return
   

def get_filtered_data(fasta, table, evalue, bitscore, pctid, coverage):
    
    print('filtering blast results by quality')

    table_df, orig_len = read_blast_table(table)

    # filter by e-value etc
    table_df = filter(table_df, evalue=evalue, bit_score=bitscore, pct=pctid, qcover=coverage)


    new_len = len(table_df)

    print('filtered  out ' + str(orig_len - new_len) + ' sequences')
    print('total original: ' +  str(orig_len))
    print('total filtered: ' + str(new_len))

    print('making filtered fasta')
    filtered_records = get_filtered_fasta_records(fasta, table_df)
    filtered_table = table_df

    return(filtered_records, filtered_table)
    

# parser 
if __name__ == "__main__":
    print("running filter_blast in terminal")
    import argparse
    import os
    from src import config
    parser = argparse.ArgumentParser(description="Run filter_blast")
    parser.add_argument("raw_directory", nargs='?', default=config.raw, type=str, help="name of directory with raw data")
    parser.add_argument("processed_directory", nargs='?', default=config.processed, type=str, help="name of directory with processed data")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help=".fasta file with BLAST results")
    parser.add_argument("-t", "--table", action = "store", default = False, help=".csv file with BLAST hits table")
    parser.add_argument("-e", "--evalue", type=float, action = "store", default = 10**-5, help="e-value to filter (default: 10^-5)") 
    parser.add_argument("-b", "--bitscore", action = "store", default = 0, help="bit score to filter (default: 0)") 
    parser.add_argument("-p", "--pctid", action = "store", default = 0, help="pct identity to filter (default: 0)") 
    parser.add_argument("-c", "--coverage", action = "store", default = 0, help="query coverage pct to filter (default: 0)") 
    parser.add_argument("-n", "--name", action="store", default=None, help="name for the run")

    args = parser.parse_args()

    # get original working dir
    working_dir = os.getcwd()

    # change to directory with raw data

    try:
        os.chdir(args.raw_directory)
    except:
        print ("didn't change dir")

    #run the data processing steps

    filtered_records, filtered_table = get_filtered_data(args.fasta, args.table, args.evalue, args.bitscore, args.pctid, args.coverage)

    # change to directory with processed data

    try:
        os.chdir(working_dir)
        os.chdir(args.processed_directory)
    except:
        print ("didn't change dir")
    
    write_filtered_data(args.fasta, args.table, filtered_records, filtered_table, args.name)

    # change directory back to working dir

    os.chdir(working_dir)

    print("done")
    