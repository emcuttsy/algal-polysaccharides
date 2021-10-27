# fasta_taxnomy

# generates a table with the taxon names and taxon ids associated with the protein seq IDs of a fasta
# by default, it outputs the table in the raw directory since this is 'queried' data
# IMPORTANT!
# Change the Entrez email and key to your own.

import argparse
import pandas as pd
from src import config
from Bio import SeqIO
from Bio import Entrez
import ete3

Entrez.email = "ecutts@mit.edu"
Entrez.api_key = "057f36326473f362b4123bee57854f2c3208"

'''
def get_bracketed_organism_name(description):
    # get the name in brackets from a fasta description string (i.e. [Pseudomonas cremoris])
    # returns 'NONE' for strings without names 
    try:
        return(description[description.find("[")+1:description.find("]")])
    except:
        return('NONE')
'''

def get_sequence_info(fasta):
    # gets the sequence protein accession
    # returns the list

    seq_accessions = []
    for index, record in enumerate(SeqIO.parse(fasta, "fasta")):
        seq_accessions.append(record.id.split(':')[0])
    return(seq_accessions)

def query_NCBI_taxonomy(accessions):
    print('fetching NCBI data')

    ncbi = ete3.NCBITaxa()
    ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    organisms = []
    taxids = []
    lineage_ids = []
    lineage_names = []
    success_accessions=[]

    
    for a in accessions:
 
        handle = Entrez.efetch(db='protein', id=a, retmode="xml")
        record = Entrez.read(handle)
        organism = record[0]['GBSeq_organism']
        print(organism)


        taxid = ncbi.get_name_translator([organism])

        try:
            taxid = list(taxid.values())[0][0]
            lineage_id = {j:i for i,j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()} #rank:id
            lineage_id = {j:i for j,i in lineage_id.items()if j in ranks}
            taxid_translator = ncbi.get_taxid_translator(lineage_id.values()) #id:name
            lineage_name = {rank:taxid_translator[idd] for rank, idd in lineage_id.items()}
            for r in ranks:
                if r not in lineage_name.keys():
                    lineage_name[r] = 'NONE'
                    lineage_id[r] = 'NONE'
    
            organisms.append(organism)
            taxids.append(taxid)
            lineage_ids.append(lineage_id)
            lineage_names.append(lineage_name)
            success_accessions.append(a)

        except:
            print('error on ' + organism)
            print(taxid)
        

    tax_df = pd.DataFrame.from_dict({'acc': success_accessions, 'organism': organisms, 'taxid': taxids})
    all_lineage_ids = {}
    all_lineage_names = {}

    for k in ranks:

        all_lineage_ids[k+'_id'] = [d[k] for d in lineage_ids]
        all_lineage_names[k] = [d[k] for d in lineage_names]

    tax_df = pd.concat([tax_df, pd.DataFrame.from_dict(all_lineage_ids), pd.DataFrame.from_dict(all_lineage_names)], axis=1)

    return(tax_df)


def write_output(fasta, tax_df, output_directory):
    os.chdir(output_directory)
    tax_df.to_csv(fasta.split('.')[0] + '_taxonomy.csv')
    return()

def main(fasta, output_directory):

    # check if there is not already a taxonomy file — this is time intensive
    seq_accessions = get_sequence_info(fasta)

    tax_df = query_NCBI_taxonomy(seq_accessions)

    write_output(fasta, tax_df, output_directory)
    return()

if __name__ == "__main__":

    print("getting taxonomy for input fasta sequences")  
    import argparse
    import os
    import warnings

    parser = argparse.ArgumentParser(description="All")
   # parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("input_directory", nargs='?', default=config.data + 'sequences', type=str, help="name of directory with raw data")
    parser.add_argument("output_directory", nargs='?', default=config.data + 'sequences', type=str, help="name of directory with processed data")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    
    args = parser.parse_args()
    #change dir if given
    try:
        os.chdir(args.input_directory)
    except:
        print ("didn't change dir")

    # check if there is not already a file:
    fasta_base = args.fasta.split('.')[0]

    if (os.path.isfile(fasta_base + '_taxonomy.csv')):

        print('found preexisting taxonomy table at ' + fasta_base + '_taxonomy.csv')
    
    elif (os.path.isfile(args.output_directory + '/' + fasta_base + 'taxonomy.csv')):
        
        print('found preexisting taxonomy table at ' + fasta_base + '_taxonomy.csv')

    # otherwise, run the thing:
    else:   

        print('no preexisting taxonomy table found. Building new table.')

        main(args.fasta, args.output_directory)

        print("done")