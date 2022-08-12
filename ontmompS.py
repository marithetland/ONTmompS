#!/usr/bin/env python3

#Import modules
from argparse import ArgumentParser, FileType
import os, sys, re, collections, operator
from subprocess import call, check_output, CalledProcessError, STDOUT
from pathlib import Path
import pandas as pd
import glob
import time
import logging


#Define functions
def parse_args():
    #Version
    parser = ArgumentParser(description='Assign mompS from a long-read assembly')
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + 'v.1.0.0')

    #Arguments
    parser.add_argument('-r', '--run_folder', type=str, required=True, help='Input directory.')
    parser.add_argument('-db', '--database_folder', type=str, required=True, help='Provide a path to database location.')
    parser.add_argument('-a', '--assembly', type=Path, required=True, help='Provide a consensus assembly in fasta format.')

    parser.add_argument('--threads', type=int, default=20, required=False, help='Specify number of threads to use. Default: 20')

    return parser.parse_args()


def run_command(command, **kwargs): #yekwahs
    command_str = ''.join(command)
    #logging.info('Running shell command: {}'.format(command_str))
    try:
        exit_status = call(command_str, **kwargs)
    except OSError as e:
        message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
        raise CommandError({"Error:": message})
    if exit_status != 0:
        message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
        raise CommandError({"Error:": message})


def check_version():

    try:
        run_command(['blastn -version >> program_versions.txt'], shell=True)
    except:
        logging.exception("Some programs were not found in PATH.")
        sys.exit("Error in checking program versions. Are you in the correct conda env?")


def check_db(db_location):
    
    try:
        if not os.path.isdir(db_location):
            raise
        if not os.path.isfile(db_location + 'mompS2_ref.tfa'):
            raise
        if not os.path.isfile(db_location + 'mompS_db.fa'):
            raise
        if not os.path.isfile(db_location + '1116R.fasta'):
            raise
    except:
        logging.exception("Location of database files not found.")
        sys.exit("Unable to find location of files in " + db_location)

    print("Database files located.")


def check_assembly(assembly_file):
    
    try:
        if not os.path.isfile(assembly_file):
            raise
    except:
        logging.exception("Assembly file not located.")
        sys.exit("Unable to find assembly file. Please check path again.")


def blast_mompS2_ref(run, assembly_file, db_location):
    
    #Define mompS2_ref file (found in db)
    mompS2_ref = db_location + 'mompS2_ref.tfa'

    try:
        #Define output file
        blast_mompS2_output = run + 'mompS2_blast.tsv'
        
        run_command(['blastn -subject ', str(assembly_file), ' -query ', str(mompS2_ref), ' -outfmt 6 > ', blast_mompS2_output], shell=True)
    except:
        logging.exception("Error in blast assembly file to mompS2_ref.")
        sys.exit("Error in blast assembly file to mompS2_ref in db.")

    return blast_mompS2_output


def read_mompS2_output(blast_mompS2_output):

    try:
        mompS2_df = pd.DataFrame(pd.read_csv(blast_mompS2_output, sep="\t", names=["Query_ID", "Subject_ID", "Perc_match", "Length", "Num_mismatches", "Num_gaps", "Query_start", "Query_end", "Subject_start", "Subject_end", "E_value", "Bitscore"], engine='python'))
    except:
        logging.exception("Unable to read mompS2 blast output.")
        sys.exit("Error in reading output from mompS2 blast with assembly fasta file")

    return mompS2_df


def get_start_stop_positions(run, mompS2_df):

    #Write the results to a file
    key_blast_file = run + 'key_blast.txt'
    key_blast = open(key_blast_file, 'w')
    key_blast.write('hit\tstart\tstop\n')

    try:
        for index, row in enumerate(mompS2_df.iterrows()):
            count = index
            if row [1][8] < row[1][9]:
                start = row[1][8]
                end = row[1][9]
            else:
                start = row[1][9]
                end = row[1][8]
            key_blast.write(str(count) + '\t' + str(start) + '\t' + str(end) + '\n')
            
    except:
        logging.exception("Unable to get start/stop positions in mompS2 dataframe.")
        sys.exit("Error in retrieving start/stop positions in mompS2 dataframe.")

def mompS_workflow(run, assembly_file, db_location, threads):
    check_version()
    check_db(db_location)
    check_assembly(assembly_file)
    blast_mompS2_output = blast_mompS2_ref(run, assembly_file, db_location)
    mompS2_df = read_mompS2_output(blast_mompS2_output)
    get_start_stop_positions(run, mompS2_df)

# main function
def main():
    args = parse_args()
    threads = args.threads
    assembly_file = args.assembly
    db_location = args.database_folder
    if db_location[-1] != '/':
        db_location = db_location + '/'
    run = args.run_folder
    if run[-1] != '/':
        run = run + '/'
    
    start_time = time.time()

    #Logging
    logging.basicConfig(filename="ONTmomps.log", format='%(asctime)s %(message)s', filemode='w', level=logging.DEBUG)
    logging.info("Started ONTmomps.")

    #mompS workflow
    mompS_workflow(run, assembly_file, db_location, threads)


    #Display times
    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    print("ONTmomps finished in " + str(time_mins) + ' mins.')


# call main function
if __name__ == '__main__':
    main()
