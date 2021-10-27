#!/usr/bin/python

# adapted from version written by abigailc@Actaeon on october 19 2016
# Note!!!!!
#   There is a problem with scp-ing from the cluster within this script
#   It copies a blank file 
#   The file is actually copied in the make_dataset script
#   Unsure why it is copying a blank file, but the work-around is to just scp the file outside of the script
#   Sorry!
# by Elise Cutts

#this script will take a fasta file, send it to mafft to align on the cluster, and then bring it home.

#set these yourself
ssh_inst = "ssh ecutts@eofe7.mit.edu"
clus_head = "ecutts@eofe7.mit.edu:/home/ecutts/"

#imports
import sys
import argparse
import os
import re
import time

#makes dir if need be
def check_directory_existance(ssh_inst):
    import os
    print("checking dirs")
    os.system(ssh_inst+" \'mkdir MAFFT\'")

#removes extra files
def remove_slurm_files(ssh_inst,pattern):
    print("removing")
    os.system(ssh_inst+" \'cd MAFFT; rm "+pattern+"\'")

#does everything
def mafft_on_cluster(filename, scriptname, runtype):
    print("starting")
    #define mafft-out name
    outname = filename.split('.')[0] +"_MAFFT.aln"
    #make dir on cluster if need be
    check_directory_existance(ssh_inst)
    #artefact
    clus_path = "/MAFFT"
    #make the script
    a = gen_script(scriptname, filename, outname, runtype)
    #a = name of script
    #current dir
    direct = os.getcwd()

    os.system("scp "+clus_head[:-1]+clus_path+"/"+outname+" "+direct)
    #move files to cluster
    move_to_cluster([filename,a], clus_path)
    #submit the .sh file on the cluster

    os.system(ssh_inst+" 'cd ~/MAFFT/;echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for MAFFT output file name
    os.system("sleep 60")
    while finished is not True:
        #try and copy it home
        print('trying ' + "scp "+clus_head[:-1]+clus_path+"/"+outname+" "+direct)
        os.system("scp -v "+clus_head[:-1]+clus_path+"/"+outname+" "+direct)
        #see if it exists at home
        exists = os.path.isfile(outname)
        if exists is True:
            finished = "yes"
        #if not, wait and try again
        else:
            finished = False
            print("waiting 5 minutes and trying again")
            time.sleep(300)
        if finished == "yes":
            print("Should be done!")
            finished = True
    print("Your file should exist")
    return outname

#this creates a single-use script to run muscle
def gen_script(scriptfile, inputfile, outputfile, runtype):

## script template
    a =  """#!/bin/bash                                                                                          
#SBATCH -p sched_mit_g4nier                                                                          
#SBATCH -t 7-00:00:00   
#SBATCH --nodes=1
#SBATCH -J """+scriptfile+"""   
#SBATCH -o """+scriptfile+""".out                                                                                         

. /etc/profile.d/modules.sh
module add engaging/mafft/7.245-with-extensions
"""+runtype+""" """+inputfile+""" > """+outputfile+"""
exit"""
    
    with open(scriptfile+".sh", "w") as script:
        script.write(a)
    return scriptfile+".sh"

#moves your .fasta and script to cluster
def move_to_cluster(list_of_files, clus_path):
    #requires os.system
    #requires scp(?)
    #do the thing
    for item in list_of_files:
        os.system("scp "+item+" "+clus_head+clus_path)
    print("Finished moving files to cluster in place:"+clus_path)

# parser

if __name__ == "__main__":

    print("Running in terminal")  
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
    parser.add_argument("-s", "--script", action = "store", default = False, help="give a name for your script/jop")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    parser.add_argument("--type", action='store', default='mafft', help='ginsi, linsi, or ensi, default mafft (auto)')


    
    
    args = parser.parse_args()
    #change dir if given
    try:
        os.chdir(args.input_directory)
    except:
        print ("didn't change dir")
    #run the thing

    outname = mafft_on_cluster(args.fasta, args.script, args.type)

    os.system("mv " + outname +  " " + args.output_directory)
    os.system("mv " + args.script + '.sh ' + args.output_directory)
    print("done")
