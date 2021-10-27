# get_species_representatives
# gets a single representative of every species 
# ignores sequences with undefined species, returns a list of these sequences.


def main(taxon_table, fasta, hits):

    return()


if __name__ == "__main__":

    print("getting species representatives")  
    import sys
    import argparse
    import os
    import re
    import time
    from src import config
    parser = argparse.ArgumentParser(description="All")
   # parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("input_directory", nargs='?', default=config.data + 'sequences', type=str, help="path to input directory")
    parser.add_argument("output_directory", nargs='?', default=config.data + 'sequences', type=str, help="path to output directory")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    parser.add_argument("-t", "--table", action = "store", default = False, help="give a hits table (.txt)")
    parser.add_argument("--taxon_table", action="store", default = None, help='taxonomy table for the fasta')
    parser.add_argument("--exclude_list", action="store", default=None, help='list of taxon names to exclude from pruning')

    args = parser.parse_args()

    if not args.taxon_table:
        os.system('python src/data/fasta_taxonomy.py -f ' + args.input_directory + '/' + args.fasta)
        taxon_table = args.fasta.split('.fasta')[0] + '_taxonomy.csv'
    else:
        taxon_table = args.taxon_table
    
    #change dir if given
    try:
        os.chdir(args.input_direftory)
    except:
        print ("didn't change dir")

    #run the thing
    main(taxon_table, args.fasta, args.table)

    print("done")
