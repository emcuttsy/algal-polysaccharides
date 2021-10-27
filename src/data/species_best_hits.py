# subsample_terminal_taxa
# gets the best hit representative from each species
# note that it will take multiple copies of the same sequence with different ranges
# so if the sequence has two ranges that are included then it will take both potentially

import pandas as pd
from Bio import SeqIO

def clean_id(record):
    # gets the clean id string
    # edit depending on your dirty strings
    return(record.id.split(':')[0])    

def get_accessions_by_species(taxon_table, fasta):
    # creates dictionary of key: species item: list of accessions
    taxon_df = pd.read_csv(taxon_table, index_col=0)
    fasta_accs = [clean_id(record) for record in SeqIO.parse(fasta, 'fasta')]
    taxon_df = taxon_df[taxon_df['acc'].isin(fasta_accs)] # get only the sequences in the fasta

    accessions_by_species = dict.fromkeys(set(taxon_df['species']))
    original_number = len(taxon_df)
    for species in accessions_by_species.keys():
        accessions_by_species[species] = taxon_df[taxon_df['species'] == species]['acc'].tolist()
    return(accessions_by_species, original_number)

def get_best_hits(accessions_list, hits_df, keep_num ):
    hits_subdf = hits_df[hits_df['subject_acc'].isin(accessions_list)]
    hits_subdf_sorted = hits_subdf.sort_values(by='bit_score', ascending=False)
    best_hits = hits_subdf_sorted['subject_acc'].tolist()[0:keep_num]
    return(best_hits)

def ignore(species, taxon_df, ignore_list):
    # return true if species is in a taxon that was set to ignore
    # PROBLEM
    taxonomy = taxon_df[taxon_df['species'] == species]
    
    for i in ignore_list:
        if i in taxonomy.values:
            return True
    else:
        return False

def get_fasta(best_hits, fasta):

    to_write = []
    counter = 0
    for index, record in enumerate(SeqIO.parse(fasta, "fasta")):
        if record.id.split(':')[0] in best_hits:
            to_write.append(record)
            counter +=1
    to_write_ids = [i.id.split(':')[0] for i in to_write]


    return(to_write)



def main(taxon_table, fasta, hits, keep, keep_num, ignore_list, output_directory):

    accessions_by_species, original_number = get_accessions_by_species(taxon_table, fasta)
    print('Input sequences: ' + str(original_number))
    hits_df = pd.read_csv(hits)

    print('Keeping the top ' + str(keep_num) + ' hits from each species')

    if ignore_list != None:
        print('Ignored taxa will not be subsampled:')
        for i in ignore_list:
            print(i)
        taxon_df = pd.read_csv(taxon_table, index_col=0)

    if keep != None:
        print('Keeping sequences: ')
        for i in keep:
            print(i)
        print('These sequences will not be dropped')

    best_hits = []
    for species in accessions_by_species.keys():
        if not ignore(species, taxon_df, ignore_list):
            best = get_best_hits(accessions_by_species[species], hits_df, keep_num)
        else:
            best = hits_df[hits_df['subject_acc'].isin(accessions_by_species[species])]['subject_acc'].tolist()
        for i in best:
            best_hits.append(i)

    to_write = get_fasta(best_hits, fasta)

    os.chdir(output_directory)
    print()
    print('Dropped sequences: ' + str(original_number - len(best_hits)))
    print('Remaining sequences: ' + str(len(best_hits)))

    # add in some reports about how many sequences are removed!

    SeqIO.write(to_write, fasta.split('.fasta')[0] + '_' + str(keep_num) + 'besthits.fasta', 'fasta')

    return()


if __name__ == "__main__":
    print("getting species best hits") 
    import argparse
    import os
    from src import config
    import warnings

    parser = argparse.ArgumentParser(description="All")
   # parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("input_directory", nargs='?', default=config.data, type=str, help="name of directory with raw data")
    parser.add_argument("output_directory", nargs='?', default=config.data, type=str, help="name of directory with processed data")
    parser.add_argument("-k", "--keep", action = "store", default = None, help="list of sequence IDs to keep")
    parser.add_argument("-n", "--number", type=int, action="store", default = 1, help="number of sequences to keep per species (default = 1)")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    parser.add_argument("-t", "--hits", action="store", default=False, help='give a hits table')
    parser.add_argument("-g", "--ignore_list", action="store", nargs='+', default=None, help='give taxa to ignore (not to trim)')
    parser.add_argument("--tax", action="store", default = None, help='taxonomy table for the fasta')

    args = parser.parse_args()

    if not args.tax:
        print('A taxon table is needed for this script. Making or finding taxon table...')
        os.system('python src/data/fasta_taxonomy.py -f ' + args.raw_directory + '/' + args.fasta)
        taxon_table = args.fasta.split('.fasta')[0] + '_taxonomy.csv'
    else:
        taxon_table = args.tax
    
    #change dir if given
    try:
        os.chdir(args.input_directory)
    except:
        print ("didn't change dir")

    #run the thing
    main(taxon_table, args.fasta, args.hits, args.keep, args.number, args.ignore_list, args.output_directory)

    print("done")
