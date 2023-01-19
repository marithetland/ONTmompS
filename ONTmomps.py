#!/usr/bin/env python3

#Import modules
from argparse import ArgumentParser, FileType
from dataclasses import replace
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
    parser = ArgumentParser(description='Assign mompS from long-read assemblies')
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + 'v.1.1.0')

    #Argsgroups
    input_args = parser.add_argument_group('Input options (required)')
    output_args = parser.add_argument_group('Output options')
    optional_args = parser.add_argument_group('Optional flags')

    #Arguments 
    #TODO: input_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=True, help='FASTA file(s) for assemblies (*.fasta)')
    input_args.add_argument('-a', '--assemblies', type=Path, required=True, help='Provide assemblies in fasta format.')
    output_args.add_argument('-o', '--outdir', type=str, required=False, default='./ONTmomps_output/', help='Output directory for all output files. Default: ./ONTmomps_output/')
    optional_args.add_argument('-d', '--database_folder', type=str, required=False, help='Provide a path to database location if different than that provided by this tool.')
    optional_args.add_argument('-t', '--threads', type=int, default=4, required=False, help='Specify number of threads to use. Default: 4')
    optional_args.add_argument('-k', '--keep_intermediate_files', type=str, choices=["on","off"], required=False, default="off", help='Keep intermediate files. Default=off.')

    return parser.parse_args()


def run_command(command, **kwargs):
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


def check_version(out):
    """Check that programs are in PATH before running and record version"""
    logging.info("Checking program versions.")

    try:
        run_command(['blastn -version >> ', out , 'program_versions.txt'], shell=True)
    except:
        logging.exception("Some programs were not found in PATH.")
        sys.exit("Error in checking program versions. Are you in the correct conda env?")


def check_dbdir(args):
    if args.database_folder:
        db_location = args.database_folder
        if db_location[-1] != '/':
            db_location = db_location + '/'
    else:
        base_dir = os.path.dirname(os.path.realpath(__file__))
        db_location=(base_dir+'/db/')

    return db_location



def check_db(db_location):
    """Check that database location exists and contains required files"""
    logging.info("Checking database for required files.")

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

def check_outdir(args):
    if args.outdir:
        out = args.outdir
        if out[-1] != '/':
            out = out + '/'
    else:
        out = "./ONTmomps_output/"

    if not os.path.isdir(out):
        os.makedirs(out)
    else:
        sys.exit('Error: Output folder already exists.')

    print("Output will be written to: " + os.path.abspath(out))
    return out


def check_assembly(assembly_file,out):
    """Check that assembly is a file and symlink"""
    logging.info("Checking assembly input.")

    try:
        if not os.path.isfile(assembly_file):
            raise
    except:
        logging.exception("Assembly file not located.")
        sys.exit("Unable to find assembly file. Please check path again.")

    try:
        source=os.path.abspath(assembly_file)
        base_name=os.path.basename(assembly_file)
        destination=out+base_name
        os.symlink(source, destination)
        assembly_file=destination

    except:
        logging.exception("Symlink failed.")
        sys.exit("Symlink failed. Please check your input.")

    print("Processing assembly file: " + source)

    return assembly_file



def blast_mompS2_ref(out, assembly_file, db_location):
    """Run blast of mompS2 ref against assembly"""
    logging.info("Running blast of mompS2 ref against assembly.")

    #Define mompS2_ref file (found in db)
    mompS2_ref = db_location + 'mompS2_ref.tfa'

    try:
        #Define output file
        blast_mompS2_output = out + 'mompS2_blast.tsv'
        
        run_command(['blastn -subject ', str(assembly_file), ' -query ', str(mompS2_ref), ' -outfmt 6 > ', blast_mompS2_output], shell=True)
    except:
        logging.exception("Error in blast assembly file to mompS2_ref.")
        sys.exit("Error in blast assembly file to mompS2_ref in db.")

    return blast_mompS2_output


def read_mompS2_output(blast_mompS2_output):
    """Read output of blastn from mompS2 ref against assembly."""
    logging.info("Reading output of blast from mompS2 ref against assembly.")

    try:
        mompS2_df = pd.DataFrame(pd.read_csv(blast_mompS2_output, sep="\t", names=["Query_ID", "Subject_ID", "Perc_match", "Length", "Num_mismatches", "Num_gaps", "Query_start", "Query_end", "Subject_start", "Subject_end", "E_value", "Bitscore"], engine='python'))
        
        extracted_sequences_list = []
        for index, row in enumerate(mompS2_df.iterrows()):
            contig = index
            extracted_sequences_list.append(contig)

    except:
        logging.exception("Unable to read mompS2 blast output.")
        sys.exit("Error in reading output from mompS2 blast with assembly fasta file")

    return mompS2_df, extracted_sequences_list


def get_start_stop_positions(out, mompS2_df, extracted_sequences_list):
    """Get start and stop positions from blast output"""
    logging.info("Getting start and stop positions of momps contigs in assembly file.")

    #Write the results to a file
    key_blast_file = out + 'key_blast.txt'
    key_blast = open(key_blast_file, 'w')
    key_blast.write('contig\tstart\tstop\n')

    try:
        for index, row in enumerate(mompS2_df.iterrows()):
            for contig in extracted_sequences_list:
                if index == contig:
                    if row [1][8] < row[1][9]:
                        start = row[1][8]
                        end = row[1][9]
                    else:
                        start = row[1][9]
                        end = row[1][8]
                    key_blast.write(str(contig) + '\t' + str(start) + '\t' + str(end) + '\n')
        
        key_blast.close()
            
    except:
        logging.exception("Unable to get start/stop positions in mompS2 dataframe.")
        sys.exit("Error in retrieving start/stop positions in mompS2 dataframe.")

    return key_blast_file


def samtools_index(assembly_file):
    """Run samtools index"""
    logging.info("Indexing assembly file.")

    try:
        run_command(['samtools faidx ', str(assembly_file)], shell=True)
    except:
        logging.exception("Unable to index assembly.")
        sys.exit("Error in indexing assembly.")


def samtools_extract(assembly_file, key_blast_file, extracted_sequences_list, out):
    """Run samtools faidx"""
    logging.info("Extracting momps contigs from assembly.")

    print(assembly_file)
    print(key_blast_file)
    print(extracted_sequences_list)

    try:
        key_blast_df = pd.DataFrame(pd.read_csv(key_blast_file, sep="\t"))
        key_blast_df = key_blast_df.set_index('contig')

        contig_start_stop_list = []
        for index, row in enumerate(key_blast_df.iterrows()):
            for contig in extracted_sequences_list:
                if index == contig:
                    if key_blast_df.loc[index, 'start']:
                        start = key_blast_df.loc[index, 'start']
                        stop = key_blast_df.loc[index, 'stop']
                        extracted_sequence = out + str(index) + '_contig.fasta'
                        to_append = [index, start, stop]
                        contig_start_stop_list.append(to_append)
                        run_command(['samtools faidx ', str(assembly_file), ' NC_002942.5:' + str(start) + '-' + str(stop), ' > ', extracted_sequence], shell=True) #TODO update samtools faidx command

    except:
        logging.exception("Unable to extract momps from assembly.")
        sys.exit("Error in extracting momps from assembly.")

    return contig_start_stop_list


def blast_mompS(db_location, extracted_sequences_list, out):
    """Run blastn of contig against momps database"""
    logging.info("Running blastn of contig against momps database.")

    momps_db_file = db_location + 'mompS_db.fa'

    try:
        for item in extracted_sequences_list:
            file = out + str(item) + '_contig.fasta'
            blast_output = out + str(item) + '_contig_mompS_blast.tsv'
            run_command(['blastn -query ', file, ' -subject ', momps_db_file, ' -perc_identity 100 -max_target_seqs 1 -outfmt 6 > ', blast_output], shell=True)

    except:
        logging.exception("Error in blastn mompS db against extracted contig.")
        sys.exit("Error in blastn mompS db against extracted contig.")


def read_blast_momps_output(out, extracted_sequences_list):
    """Read output files from blastn of contig fasta against momps database"""
    logging.info("Reading blast output from contig against momps database.")

    key_blast_momps_file = out + 'key_momps_blast.txt'
    key_blast_momps = open(key_blast_momps_file, 'w')
    key_blast_momps.write('contig\tallele\tlength\n')

    try:
        for contig in extracted_sequences_list:
            file = out + str(contig) + '_contig_mompS_blast.tsv'
            if os.stat(file).st_size == 0:
                key_blast_momps.write(str(contig) + '\t' + str('.') + '\t' + str('.') + '\n')
                continue
            df = pd.DataFrame(pd.read_csv(file, sep="\t", names=["Query_ID", "Subject_ID", "Perc_match", "Length", "Num_mismatches", "Num_gaps", "Query_start", "Query_end", "Subject_start", "Subject_end", "E_value", "Bitscore"], engine='python'))
            for index, row in enumerate(df.iterrows()):
                length = df.loc[index, 'Length']
                momps = df.loc[index, 'Subject_ID']
            key_blast_momps.write(str(contig) + '\t' + str(momps) + '\t' + str(length) + '\n')

        key_blast_momps.close()

    except:
        logging.exception("Error in reading blast output.")
        sys.exit("Error in reading momps blast output.")

    return key_blast_momps_file


def run_water_alignment(db_location, extracted_sequences_list, out):
    """Run pairwise alignment"""
    logging.info("Running water alignment on contig fasta files with 1116R primer.")

    primer_1116R = str(db_location + '1116R.fasta')

    try:
        for contig in extracted_sequences_list:
            file = out + str(contig) + '_contig.fasta'
            water_output = out + str(contig) + '_contig.water'
            run_command(['water ', file, ' ', primer_1116R, ' -gapopen 10 -gapextend 0.5 -outfile ', water_output], shell=True)

    except:
        logging.exception("Error running water pairwise alignment.")
        sys.exit("Error running water pairwise alignment.")


def read_water_alignment(out, extracted_sequences_list):
    """Read alignment files"""
    logging.info("Reading alignment files.")

    key_water_alignment_file = out + 'key_water_alignment.txt'
    key_water_alignment = open(key_water_alignment_file, 'w')
    key_water_alignment.write('contig\tmomps_copy\n')

    try:
        momps2_contig_list = []
        for contig in extracted_sequences_list:
            file = out + str(contig) + '_contig.water'
            with open(file, 'r', encoding='utf-8') as output:
                file_contents = output.readlines()
                for pos, line in enumerate(file_contents):
                    if line.startswith("# Score: 125.0"):
                        momps2_contig = contig
                        momps2_contig_list.append(momps2_contig)
        
        if momps2_contig_list:
            for i in momps2_contig_list:
                key_water_alignment.write(str(i) + '\t' + 'mompS2' + '\n')

            mompS1_list = list(set(extracted_sequences_list).difference(momps2_contig_list))
            for i in mompS1_list:
                key_water_alignment.write(str(i) + '\t' + 'mompS1' + '\n')

        else:
            print("#######NB#######")
            print("mompS2 could not be separated from mompS1 by exact match to 1116R.")
            print("Please check .water files to manually deduce.")
            print("#######NB#######")
            for i in extracted_sequences_list:
                key_water_alignment.write(str(i) + '\t' + 'unable_to_distinguish_mompS_copy' + '\n')
        
        key_water_alignment.close()

    except:
        logging.exception("Error in reading blast output.")
        sys.exit("Error in reading momps blast output.")

    return key_water_alignment_file


def make_report(out, key_blast_file, key_blast_momps_file, key_water_alignment_file):
    """Create the final report"""
    logging.info("Creating the final report.")

    blast_df = pd.DataFrame(pd.read_csv(key_blast_file, sep="\t"))
    momps_df = pd.DataFrame(pd.read_csv(key_blast_momps_file, sep="\t"))
    water_df = pd.DataFrame(pd.read_csv(key_water_alignment_file, sep="\t"))

    merged_momps_water_df = pd.merge(water_df, momps_df,  on='contig', how='left')
    final_df = pd.merge(merged_momps_water_df, blast_df, on='contig', how='left')
    pd.DataFrame.to_csv(final_df, out + 'report_momps.csv', sep=',', index=False)
    print("Report written to report_momps.csv")
    logging.info("Run report created.")


def mompS_workflow(out, assembly_input, db_location, threads):
    check_version(out)
    check_db(db_location)
    assembly_file = check_assembly(assembly_input,out)
    blast_mompS2_output = blast_mompS2_ref(out, assembly_file, db_location)
    mompS2_df, extracted_sequences_list = read_mompS2_output(blast_mompS2_output)
    key_blast_file = get_start_stop_positions(out, mompS2_df, extracted_sequences_list)
    samtools_index(assembly_file)
    contig_start_stop_list = samtools_extract(assembly_file, key_blast_file, extracted_sequences_list, out)
    blast_mompS(db_location, extracted_sequences_list,out)
    key_blast_momps_file = read_blast_momps_output(out, extracted_sequences_list)
    run_water_alignment(db_location, extracted_sequences_list, out)
    key_water_alignment_file = read_water_alignment(out, extracted_sequences_list)
    make_report(out, key_blast_file, key_blast_momps_file, key_water_alignment_file)

# main function
def main():
    args = parse_args()
    threads = args.threads
    start_time = time.time()
    out = check_outdir(args)
    db_location = check_dbdir(args)

    #Set up input (TODO: accept multiple input assemblies):
    assembly_input = args.assemblies

    #Logging 
    logging.basicConfig(filename=out+"ONTmomps.log", format='%(asctime)s %(message)s', filemode='w', level=logging.DEBUG)
    logging.info("Started ONTmomps.")

    #run mompS workflow
    mompS_workflow(out, assembly_input, db_location, threads)

    #TODO: remove intermediate files 
    #if not args.keep_intermediate_files:
    #    os.remove(


    #Display times
    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    logging.info("ONTmomps completed.")
    print("ONTmomps finished in " + str(time_mins) + ' mins.')


# call main function
if __name__ == '__main__':
    main()
