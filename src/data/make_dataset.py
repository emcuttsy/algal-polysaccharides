# -*- coding: utf-8 -*-
import os
import logging
from pathlib import Path
from src import config

fasta = 'S_japonica_MC5E_234-452_BLAST_1e-5.fasta'
table = 'S_japonica_MC5E_234-452_BLAST_1e-5_hits.csv'
ssh_inst = "ssh ecutts@eofe7.mit.edu"
clus_head = "ecutts@eofe7.mit.edu:/home/ecutts/"
evalue = 0.00005 #1e-5
keep_number = 2 #number of best hits to keep when filtering per-species hits

# eventually rewrite so this can be done on multiple sequences

def main(raw_directory, processed_directory):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info('making final data set from raw data')

    # remove duplicate sequences
    os.system('python src/data/remove_duplicates.py ' + config.data + 'blast_raw/ ' + config.data + 'blast_post/ -f ' + fasta +\
        ' -t ' + table)
    print()

    # no dupes name
    fasta_no_dupes = fasta.split('.')[0] + '_no_dupes.fasta'
    hits_no_dupes = fasta.split('.')[0] + '_no_dupes_hits.csv'


    # get taxonomy information for sequences
    os.system('python src/data/fasta_taxonomy.py ' + config.data + 'blast_post/ ' + config.data + 'blast_post/ -f' + fasta_no_dupes)

    # remove all of the pseudomonas bloat
    os.system('python src/data/prune_taxon.py ' + config.data + 'blast_post/ ' + config.data + 'blast_post/ -f ' + fasta_no_dupes +\
        ' -k pseudomonas_to_keep.txt -t Pseudomonas --taxon_table ' + fasta_no_dupes.split('.fasta')[0] + '_taxonomy.csv')

    #pruned_name
    fasta_pruned = fasta_no_dupes.split('.')[0] + '_pruned.fasta'

    # retain only the best 2 sequences for each terminal taxon.
    os.system('python src/data/species_best_hits.py ' + config.data + 'blast_post/ ' + config.data + 'blast_post/ -f ' + fasta_pruned + \
        ' -t ' + hits_no_dupes + ' --tax ' + fasta_no_dupes.split('.')[0] + '_taxonomy.csv  -g Eukarya Archaea -n 2')
    
    fasta_best_hits = fasta_pruned.split('.')[0] + '_' + str(keep_number) + 'besthits.fasta'

    # run mafft G-INS-I on the cluster
    os.system('python src/data/cluster_mafft.py ' + config.data + 'blast_post/ ' + config.data + 'alignments/ -f ' + fasta_best_hits + \
        ' -s ' + fasta.split('.')[0] + 'mafftGINSI --type ginsi')

    
    os.system("scp " + clus_head + 'MAFFT/' + fasta_best_hits.split('.')[0] + '_MAFFT.aln ' + config.data + 'alignments')

if __name__ == '__main__':
    
    import argparse
    from datetime import datetime

    # datetime object containing current date and time
    now = datetime.now()

    # dd/mm/YY-H:M
    dt_string = now.strftime("%d-%m-%Y-%H:%M")

    parser = argparse.ArgumentParser(description="Make dataset")
    parser.add_argument("raw_directory", nargs='?', default=config.raw, type=str, help="name of directory with raw data")
    parser.add_argument("processed_directory", nargs='?', default=config.processed, type=str, help="name of directory with processed data")
    args = parser.parse_args()


    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    main(args.raw_directory, args.processed_directory)
