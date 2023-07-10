#!/usr/bin/env python3

"""
Copyright 2023 Marit Hetland (marit[dot]hetland[at]outlook[dot]com)
Copyright 2023 Markus Soma
https://github.com/marithetland/ONTmompS
"""

#Import modules 
from argparse import ArgumentParser, FileType
import distutils.spawn
from dataclasses import replace
import os, sys, re, collections, operator
from subprocess import call, check_output, CalledProcessError, STDOUT
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Emboss.Applications import WaterCommandline
import glob
import time
import logging


# Define functions
def parse_args():
    #Version 
    parser = ArgumentParser(description='In silico SBT of Legionella pneumophila from long-read or hybrid assemblies')
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + 'v2.0.0')

    #Argsgroups
    input_args = parser.add_argument_group('Input options (required)')
    optional_args = parser.add_argument_group('Optional flags')
    output_args = parser.add_argument_group('Output options')

    #Arguments 
    input_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=True, help='FASTA file(s) for assemblies (*.fasta)')

    optional_args.add_argument('--db', type=str, required=False, help='Provide a path to database location if different than that provided by this tool.')
    optional_args.add_argument('--store_mompS_alleles', action='store_true', required=False, help='Print mompS alleles to files named {allele}_{assembly}.fna.')
    optional_args.add_argument('--store_novel_alleles', action='store_true', required=False, help='Print novel alleles to files named {allele}_{assembly}.fna.')
    optional_args.add_argument('--store_all_alleles', action='store_true', required=False, help='Print all alleles (7 genes in SBT scheme + mompS1) to files named {allele}_{assembly}.fna. ')
    optional_args.add_argument('--verbose', action='store_true', required=False,  help='Log more details and keep intermediate files for debugging.') #Set this to store the water alignment.
    optional_args.add_argument('-l','--log', type=str, required=False, help='Write logging to specified file name instead of stdout') #Set this to store the water alignment.

    output_args.add_argument('--ST_outfile', type=str, required=False, default='./LpST_ONTmompS.tsv', help='Output filename for STs. Default: ./LpST_ONTmompS.tsv')
    output_args.add_argument('--mompS_outfile', type=str, required=False, default='./mompS_alleles_ONTmompS.tsv', help='Output filename for mompS copy allele numbers. Default: ./mompS_alleles_ONTmompS.tsv')
    output_args.add_argument('-o','--outdir', type=str, default='./ONTmompS_allele_sequences', help='Output directory to store novel alleles in. Default is current working directory')

    return parser.parse_args()

#####################################################################
### Defs for checking verions, arguments, databases
#####################################################################

def set_up_logging(args):
    """ Set up logging """
    if args.log:
        logfile = args.log #Log to provided file path
    else:
        logfile = None #Log to terminal output
    # Set up log to stdout
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG, 
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m-%d-%Y %H:%M:%S')

def run_command(command, **kwargs):
    """ Set up to run shell commands """
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

def check_dbdir(args):
    """ Check that db dir exists and return the path """
    if args.db:
        db_location = args.db
        if db_location[-1] != '/':
            db_location = db_location + '/'
    else:
        base_dir = os.path.dirname(os.path.realpath(__file__))
        db_location=(base_dir+'/db/')

    return db_location


def check_db(db_location):
    """ Check that database location exists and contains required files, and save to variable """
    #Assign variable names to files
    sts = db_location + 'lpneumophila.txt'
    flaA_db_file = db_location + 'flaA.fna'
    pilE_db_file = db_location + 'pilE.fna'
    asd_db_file = db_location + 'asd.fna'
    mip_db_file = db_location + 'mip.fna'
    mompS_db_file = db_location + 'mompS.fna'
    proA_db_file = db_location + 'proA.fna'
    neuA_db_file = db_location + 'neuA.fna'

    base_dir = os.path.dirname(os.path.realpath(__file__))
    mompS2_ref=(base_dir+'/db/mompS2_ref.fna')
    primer_1116R=(base_dir+'/db/1116R_rev.fasta')

    #Check that the files exist
    try:
        if not os.path.isdir(db_location):
            raise
        if not os.path.isfile(sts):
            raise
        if not os.path.isfile(mompS2_ref):
            raise
        if not os.path.isfile(primer_1116R):
            raise
        if not os.path.isfile(flaA_db_file):
            raise
        if not os.path.isfile(pilE_db_file):
            raise
        if not os.path.isfile(asd_db_file):
            raise
        if not os.path.isfile(mip_db_file):
            raise
        if not os.path.isfile(mompS_db_file):
            raise
        if not os.path.isfile(proA_db_file):
            raise
        if not os.path.isfile(neuA_db_file):
            raise

    except:
        logging.exception("Location of database files not found.")
        sys.exit("Unable to find location of files in " + db_location)

    logging.info("Using database located at: " + db_location)

    return sts, mompS2_ref, primer_1116R, flaA_db_file, pilE_db_file, asd_db_file, mip_db_file, mompS_db_file, proA_db_file, neuA_db_file

def load_st_db(sts):
    """ Load the lpneumophila.txt ST database, using same logic as Kleborate """
    alleles_to_st = {}
    header = []
    st_names = []
    with open(sts, 'r') as f:
        for line in f:
            fields = line.rstrip().split('\t')
            if len(header) == 0:
                header = fields
                header.pop(0)  # remove ST label
            else:
                st = fields.pop(0)
                st_names.append(st)
                alleles_to_st[','.join(fields)] = st   #dict, like: '6,4,15,28,21,14,1': '2809'

    return header, st_names, alleles_to_st

def check_outputs(args):
    """ Check that specified outputs/paths are OK """
    #Outfile for STs 
    if args.ST_outfile:
        ST_outfile = args.ST_outfile
    else:
        ST_outfile = "./LpST_ONTmompS.tsv"

   #Outfile for mompS  
    if args.mompS_outfile:
        mompS_outfile = args.mompS_outfile
    else:
        mompS_outfile = "./mompS_alleles_ONTmompS.tsv"

    logging.info("Output STs will be written to: " + os.path.abspath(ST_outfile))
    logging.info("Output mompS copy allele numbers will be written to: " + os.path.abspath(mompS_outfile))

    #Outdir for storing allele sequences or verbose files
    if args.store_mompS_alleles or args.store_novel_alleles or args.store_all_alleles or args.verbose:
        if args.outdir:
            out = args.outdir
            if out[-1] != '/':
                out = out + '/'
        else:
            out = "./ONTmompS_allele_sequences"

        if not os.path.isdir(out):
            os.makedirs(out)
        else:
            sys.exit('Error: Output folder already exists.')
        if args.store_mompS_alleles:
            store="The mompS allele sequences"
        if args.store_novel_alleles:
            store="Any novel allele sequences"
        if args.store_all_alleles:
            store="All allele sequences"
        if args.verbose:
            store="Verbose mode is on. The 987 bp mompS-region sequences and pairwise alignment files"

        logging.info(store + " will be written to files in: " + os.path.abspath(out) + "/") 
    else:
         out = "./" #Write any output (ST and mompS alleles) to current working directory

    return out, ST_outfile, mompS_outfile

def isfasta(file):
    """ Check that input file is in fasta format """
    with open(file, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

def check_assemblies(args):
    """ Check input and return list to iterate over """
    #Set up a list of all input files:
    fasta=args.assemblies
    #Create list of input sequences   
    genomes = []
    sequence_list = []

    for sequence in fasta:
        if os.path.isfile(sequence) and isfasta(sequence) and sequence not in sequence_list: #If file exists and is a fasta sequence and is not already in list
            sequence_list.append(sequence)
        else: 
            logging.info("Provided file " + sequence +" does not exist or is not in FASTA format. Please check your input. Removing from downstream analysis.")

    if sequence_list != []:
        return sequence_list
    else:
        sys.exit("No FASTA files were provided. Please check that your input/paths are correct.")


def check_assembly(assembly_file,out):
    """ Check that assembly is a file and return path and name """
    try:
        assembly_file=os.path.abspath(assembly_file)
        base_name=os.path.basename(assembly_file)
        assembly_name=os.path.splitext(base_name)[0]

    except:
        logging.exception("Symlink failed.")
        sys.exit("Symlink failed. Please check your input.")

    return assembly_file, assembly_name

#####################################################################
### Defs for BLASTn and parsing BLASTn output 
#####################################################################
def run_blastn(subject, query):
    """ Run blastn """
    #Run BLASTn to find allele hits. Uses the same BLASTn logic as Kleborate's MLSTblast.py. 
    try:
        blastn_command = NcbiblastnCommandline(cmd='blastn',
                                  task='blastn',
                                  query=query, 
                                  subject=subject, 
                                  dust='no',
                                  evalue='1E-20',
                                  word_size='32',
                                  max_target_seqs='10000',
                                  perc_identity='90',
                                  outfmt='"6 sacc pident slen qlen length score gaps mismatch bitscore sseq qseq sstrand sstart send qacc qstart qend qframe"')
        blastn_output = blastn_command()[0].strip() 

    except:
        logging.exception("Error in blast of assembly file to mompS2_ref db.")

    return blastn_output

def read_blastn_output(blastn_output):
    """ Read BLASTn output as a pandas dataframe """
    try:
        headers = ['sacc','pident','slen','qlen','length','score','gaps','mismatch','bitscore','sseq','qseq','sstrand','sstart','send','qacc','qstart','qend','qframe'] 
        rows = [line.split() for line in blastn_output.splitlines()]    
        df = pd.DataFrame(rows, columns=headers)

        #BLASTn runs an alignment. Remove any "-" characters from sequence (deletions relative to reference):
        df['sseq'] = df['sseq'].str.replace('-', '')
        df['qseq'] = df['qseq'].str.replace('-', '')

    except:
        logging.exception("Error in reading BLASTn df.")
        sys.exit("Error in reading BLASTn df.") 

    return df

def get_sseq_start_stop_positions(row,strand,seq,start,end):
    """ Get the sequence, start and stop positions from blast output. Reverse complement if on minus strand """
    #Reverse complement if on minus strand. Store sequence and start + end positions 
    if row[strand] == 'minus':
        my_seq = Seq(row[seq])
        sseq=my_seq.reverse_complement()
        sstart=row[end]
        send=row[start]
    elif row[strand] == 'plus':
        sseq = Seq(row[seq])
        sstart=row[start]
        send=row[end]

    return sseq, sstart, send 

#####################################################################
### Defs for identifying mompS and SBT 
#####################################################################

def get_mompS_blastn_hits(assembly_file, mompS2_ref,mompS_db_file,primer_1116R, out, assembly_name, args):
    """ Identify the mompS1 and mompS2 sequences """
    try:
        blastn_output = run_blastn(assembly_file, mompS2_ref)
        df = read_blastn_output(blastn_output)

        mompS_blastn_hits = []
        mompS2_copy = []
        mompS2_sequence = []
        mompS1_copy = [] 
        mompS1_sequence = []
        num_mompS_hits = 0

        for index, row in df.iterrows():
            if float(row['length'])/float(row['qlen']) >= 0.5: #Minimum half of the alignment needs to be met
                num_mompS_hits+=1
                sseq, sstart, send = get_sseq_start_stop_positions(row,strand='sstrand',seq='sseq',start='sstart',end='send')
                mompS_copy = run_water_alignment(row, primer_1116R, out, args, assembly_name, num_mompS_hits, seq='sseq')
                to_append = [mompS_copy, sseq, sstart, send]
                mompS_blastn_hits.append(to_append)

        #Assign mompS copies and store sequence to respective files:
        for hit in mompS_blastn_hits:
            if hit[0]=='mompS2':
                mompS2_copy=hit[0]
                mompS2_sequence=hit[1]
                run_command(['echo ">',mompS2_copy,'__',assembly_name,'" >>  ',out,mompS2_copy,'_seq_with_flank_', assembly_name,'.fasta ; echo ', str(mompS2_sequence),' >>  ',out,mompS2_copy,'_seq_with_flank_', assembly_name,'.fasta'  ], shell=True) 
            elif hit[0]=='mompS1':
                mompS1_copy=hit[0]
                mompS1_sequence=hit[1]
                run_command(['echo ">',mompS1_copy,'__',assembly_name,'" >>  ',out,mompS1_copy,'_seq_with_flank_', assembly_name,'.fasta ; echo ', str(mompS1_sequence),' >>  ',out,mompS1_copy,'_seq_with_flank_', assembly_name,'.fasta'  ], shell=True) 

        #Warn if num_mompS_hits != 2
        if num_mompS_hits > 2:
            logging.info("Warning, more than two mompS hits were found. Legionella pneumophila typically has two mompS copies. Please inspect your assembly quality. If you believe this is true, please let us know in a GitHub issue so that we can update the code.")
        elif num_mompS_hits < 2: 
            logging.info("Warning, only one mompS hit was found in your genome. Legionella pneumophila typically has two mompS copies. Please inspect your assembly quality - it is possible that mompS1 and mompS2 have been merged in assembly, which could lead to incorrect typing.")

        #If any of the mompS copies are not found, store them as missing "-" 
        mompS_not_identified = []
        if not mompS2_copy:
            mompS2_copy="-"
            to_print = (assembly_name + ": mompS2 not found")
            mompS_not_identified.append(to_print)
        if not mompS1_copy:
            mompS1_copy="-"
            to_print = (assembly_name + ": mompS1 not found")
            mompS_not_identified.append(to_print)


    except:
        logging.exception("Error: Did not find any hits to mompS in file: " + assembly_name)
        sys.exit("Error: Did not find any hits to mompS in file: " + assembly_name)

    return mompS2_copy, mompS2_sequence, mompS1_copy, mompS1_sequence, mompS_not_identified


def get_SBT(assembly_file, mompS2_ref,mompS_db_file, primer_1116R, asd_db_file, flaA_db_file, mip_db_file, neuA_db_file, pilE_db_file, proA_db_file, out, assembly_name, sts, ST_alleles,ST_outfile, mompS_alleles, mompS_outfile,args):
    """ Read all locus alleles and assign ST """
    mompS2_copy, mompS2_sequence, mompS1_copy, mompS1_sequence, mompS_not_identified = get_mompS_blastn_hits(assembly_file, mompS2_ref,mompS_db_file,primer_1116R, out, assembly_name, args) 

    #Get allele numbers for the mompS copies first
    #Read mompS sequences for each copy:    
    mompS2_sequence = (out + 'mompS2_seq_with_flank_' + assembly_name + '.fasta')
    mompS1_sequence = (out + 'mompS1_seq_with_flank_' + assembly_name + '.fasta')

    if os.path.isfile(mompS2_sequence):
        mompS2_allele_annotated,  mompS2_allele_number, mompS2_allele, mompS2_sequence_temp = get_SBT_allele_blastn_hits(mompS2_sequence, mompS_db_file,assembly_name,out,args) 
        if not args.verbose:
            run_command(['rm ',mompS2_sequence ], shell=True) 
    elif not os.path.isfile(mompS2_sequence):
        mompS2_allele = "-"
        mompS2_allele_annotated = "-"
    if os.path.isfile(mompS1_sequence):
        mompS1_allele_annotated, mompS1_allele_number, mompS1_allele, mompS2_sequence_temp = get_SBT_allele_blastn_hits(mompS1_sequence, mompS_db_file, assembly_name,out,args) 
        if not args.verbose:
            run_command(['rm ',mompS1_sequence ], shell=True)
    elif not os.path.isfile(mompS1_sequence):
        mompS1_allele = "-"
        mompS1_allele_annotated = "-"

    #mompS 987 bp sequences: If verbose mode is not on, remove the stored 987 bp mompS-region sequences and water alignment sequences.  
    if args.verbose:
        logging.info("Verbose mode is on. Will print the 987 bp mompS sequences of all input assemblies to files in: " + out)

    #For each of the remaining 6 SBT genes, get the allele number and check if it is *, ? or - (done for mompS1 and mompS2 above):
    flaA_allele_annotated, flaA_allele_number, flaA_allele, flaA_sequence = get_SBT_allele_blastn_hits(assembly_file, flaA_db_file,assembly_name,out,args)
    pilE_allele_annotated, pilE_allele_number, pilE_allele, pilE_sequence = get_SBT_allele_blastn_hits(assembly_file, pilE_db_file,assembly_name,out,args)
    asd_allele_annotated, asd_allele_number, asd_allele, asd_sequence  = get_SBT_allele_blastn_hits(assembly_file, asd_db_file,assembly_name,out,args)
    mip_allele_annotated, mip_allele_number, mip_allele, mip_sequence = get_SBT_allele_blastn_hits(assembly_file, mip_db_file,assembly_name,out,args)
    proA_allele_annotated, proA_allele_number, proA_allele, proA_sequence = get_SBT_allele_blastn_hits(assembly_file, proA_db_file,assembly_name,out,args)
    neuA_allele_annotated, neuA_allele_number, neuA_allele, neuA_sequence = get_SBT_allele_blastn_hits(assembly_file, neuA_db_file, assembly_name,out,args)

    #Create a combination with all of the ST and allele numbers. 
    allele_combinations = {"flaA":flaA_allele,"pilE":pilE_allele,"asd":asd_allele,"mip":mip_allele,"mompS":mompS2_allele,"proA":proA_allele,"neuA":neuA_allele}
    annotated_allele_combinations = {"flaA":flaA_allele_annotated,"pilE":pilE_allele_annotated,"asd":asd_allele_annotated,"mip":mip_allele_annotated,"mompS":mompS2_allele_annotated,"proA":proA_allele_annotated,"neuA":neuA_allele_annotated}    

    #Assign ST
    ST_annotated = get_closest_ST(sts, allele_combinations, annotated_allele_combinations)
    #Append ST to output files and print to stdout
    append_to_ST_report(assembly_name, ST_annotated, flaA_allele_annotated, pilE_allele_annotated, asd_allele_annotated, mip_allele_annotated, mompS2_allele_annotated, proA_allele_annotated, neuA_allele_annotated, mompS1_allele_annotated, ST_alleles,ST_outfile, mompS_alleles, mompS_outfile, mompS_not_identified)

def get_SBT_allele_blastn_hits(query, subject, assembly_name, out, args):
    """ Run blastn of contig against locus database to find closest matching allele """

    try:
        blastn_output = run_blastn(subject, query) 
        df = read_blastn_output(blastn_output)

        if df.empty: #If there are no BLASTn hits, the annotation for the allele should be "-"
            annotated_allele = "-"
            allele_number = "-"
            allele = "-"
            assembly_seq = "-"
        if not df.empty:
            blastn_hits = []
            hit_number=0

            for index, row in df.iterrows(): #For each blast hit, do:
                hit_number+=1
                #qseq, qstart, qend = get_sseq_start_stop_positions(row,strand='sstrand',seq='qseq',start='qstart',end='qend') 
                to_append = [hit_number, row['sacc'], row['pident'], row['score'], row['length'], row['slen'], row['qlen'], row['gaps'], row['mismatch'], row['qstart'], row['qend'], row['qseq']] 
                blastn_hits.append(to_append)

            #Keep only the best hit. 
            allele_annotated, allele_number, assembly_strand, assembly_seq, assembly_start, assembly_end = cull_redundant_hits(blastn_hits) 
    
            #Get locus, allele number and annotated allele number
            locus_number = allele_number.split('_')
            locus = locus_number[0]
            allele = locus_number[1]
            locus_number = allele_annotated.split('_')
            annotated_allele = locus_number[1]

            #Write allele seqences to files if flags are specified
            write_allele_sequences_to_files(locus, annotated_allele, allele, assembly_name, assembly_strand, assembly_start, assembly_end,assembly_seq,out, args) 


        return annotated_allele, allele_number, allele, assembly_seq

    except:
        logging.exception("Error in blastn db against sequence.")
        sys.exit("Error in blastn db against sequence.")


def get_closest_ST(sts, allele_combinations, annotated_allele_combinations): #Adapted from Kleborate
    """ Get the closest matching ST if <3 alleles are imprecise matches. """
    header, st_names, alleles_to_st = load_st_db(sts) #Load the lpneumophila.txt database

    best_allele_numbers_hits = []
    best_allele_numbers_annotated = []
    mismatch_loci, missing_loci = 0, 0
    loci=["flaA","pilE","asd","mip","mompS","proA","neuA"] #make sure only mompS2 is used here, not mompS1

    #For each locus, count if mismatch, and append the unannotated allele number to best_allele_numbers_hits
    for locus in loci:
        allele_number = allele_combinations[locus] #Get the unannotated allele number from the dictionary  
        allele = annotated_allele_combinations[locus] #Get the annotated allele number from the dictionary 

        #Count number of genes (out of the 7 in the SBT scheme) with mismatches or incomplete coverage
        if '*?' in allele:
            mismatch_loci += 1 
        elif '*' in allele:
            mismatch_loci += 1
        elif '?' in allele:
            mismatch_loci += 1
        elif '-' in allele:
            mismatch_loci += 1
        
        #Append the allele number for each locus 
        best_allele_numbers_hits.append(allele_number)
        best_allele_numbers_annotated.append(allele)

    best_allele_number_combination = ','.join(best_allele_numbers_hits) #Convert to comma-separated

    if mismatch_loci <= 3: # only report ST if at least 3 loci are precise matches 
        if best_allele_number_combination in alleles_to_st:
            closest_st = alleles_to_st[best_allele_number_combination]
        else: #If not exact match with db, determine the closest ST
            closest_st, _, mismatch_loci = get_closest_locus_variant(best_allele_numbers_hits, best_allele_numbers_annotated, alleles_to_st)
    else:
        closest_st = '0'

    if mismatch_loci > 0 and closest_st != '0':
        closest_st += '-' + str(mismatch_loci) + 'LV'

    ST_annotated = closest_st
    return ST_annotated


def append_to_ST_report(assembly_name, ST_annotated, flaA_allele_annotated, pilE_allele_annotated, asd_allele_annotated, mip_allele_annotated, mompS2_allele_annotated, proA_allele_annotated, neuA_allele_annotated, mompS1_allele_annotated, ST_alleles,ST_outfile, mompS_alleles, mompS_outfile, mompS_not_identified):
    """ Once an iteration (assembly) is complete, append the results to the outfiles and write to stdout """
    ST_alleles["Strain"].append(assembly_name)
    ST_alleles["ST"].append(ST_annotated)
    ST_alleles["flaA"].append(flaA_allele_annotated)
    ST_alleles["pilE"].append(pilE_allele_annotated)
    ST_alleles["asd"].append(asd_allele_annotated)	
    ST_alleles["mip"].append(mip_allele_annotated)	
    ST_alleles["mompS"].append(mompS2_allele_annotated)	
    ST_alleles["proA"].append(proA_allele_annotated)	
    ST_alleles["neuA"].append(neuA_allele_annotated)

    #Convert dictionary to dataframe    
    ST = pd.DataFrame.from_dict(ST_alleles, orient='index').T
    #Add "ST" prefix to ST name 
    ST['ST'] = 'ST' + ST['ST'].astype(str)
    #Write to output file
    ST.to_csv(ST_outfile,index=False, sep = "\t")

    mompS_alleles["Strain"].append(assembly_name)
    mompS_alleles["mompS1"].append(mompS1_allele_annotated)
    mompS_alleles["mompS2"].append(mompS2_allele_annotated)

    #Convert dictionary to dataframe    
    mompS = pd.DataFrame.from_dict(mompS_alleles, orient='index').T
    #Write to output file
    mompS.to_csv(mompS_outfile,index=False, sep = "\t")
    
    #Print results temp to outfile 
    print(assembly_name + "\tST" + ST_annotated + "\t" + flaA_allele_annotated + "\t" + pilE_allele_annotated + "\t" + asd_allele_annotated + "\t" + mip_allele_annotated + "\t" + mompS2_allele_annotated + "\t" + proA_allele_annotated + "\t" + neuA_allele_annotated + "\t" + "\t" + "\t" + mompS1_allele_annotated)

    if mompS_not_identified != []: 
        print(mompS_not_identified)

    return ST, mompS_not_identified


def get_closest_locus_variant(best_allele_numbers_hits, best_allele_numbers_annotated, alleles_to_st):
    """ Get the closest locus variant in case of mismatches/incomplete sequence match. Logic adapted from Kleborate. """
    best_allele_numbers_annotated = list(best_allele_numbers_annotated) 
    closest = []
    closest_alleles = {}   #key = st, value = list
    min_dist = len(best_allele_numbers_hits) #number mismatching loci, ignoring SNPs

    for index, item in enumerate(best_allele_numbers_hits):
        if item == '-':
            best_allele_numbers_hits[index] = '0'

    # get distance from closest ST, ignoring SNPs
    for st in alleles_to_st:
        if st != '': #In case there are any empty strings in db file, ignore them
            d = sum(map(lambda x, y: bool(int(x)-int(y)), st.split(','), best_allele_numbers_hits))
            if d == min_dist:
                closest.append(int(alleles_to_st[st]))
                closest_alleles[alleles_to_st[st]] = st
            elif d < min_dist:# reset 
                closest = [int(alleles_to_st[st])]
                closest_alleles[alleles_to_st[st]] = st
                min_dist = d  # distance from closest ST, ignoring SNPs

    closest_st = str(min(closest))

    for index, item in enumerate(best_allele_numbers_annotated):
        best_allele_numbers_annotated[index] = re.sub(r'-\d+%', '', item)
        if item == '-' or '*' in item:
            best_allele_numbers_annotated[index] = '0'

    # get distance from closest ST, including SNPs
    min_dist_incl_snps = sum(map(lambda x, y: bool(int(x)-int(y)),
                                 closest_alleles[closest_st].split(','), best_allele_numbers_annotated))

    return closest_st, min_dist, min_dist_incl_snps

def write_allele_sequences_to_files(locus, annotated_allele, allele, assembly_name, assembly_strand, assembly_start, assembly_end,assembly_seq,out, args):
    """ Save/delete sequences in files depending on args given """

    #SBT allele + mompS1 allele sequences - For each allele, store if it is meant to be  
    alleles_to_store = []

    if args.store_all_alleles: #Store all alleles regarless of novel or not
        alleles_to_store.append(locus)
    elif args.store_mompS_alleles and not args.store_novel_alleles and locus == "mompS":  
        alleles_to_store.append(locus)
    elif args.store_novel_alleles and annotated_allele != allele: #Store novel alleles   
        alleles_to_store.append(locus)

    #Then loop over all and print?
    for locus in alleles_to_store:
        fasta_header = ">" + "Locus:" + locus +  " Allele number:" + annotated_allele + " Assembly:" + assembly_name + " Strand:" + assembly_strand + " Pos:" + str(assembly_start) + "-" + str(assembly_end) 
        fasta_sequence = assembly_seq
        filename = out + locus + "_" + assembly_name + ".fasta"
        f=open(filename,"w+")
        f.writelines(fasta_header)
        f.writelines("\n")
        f.writelines(fasta_sequence)
        f.writelines("\n")


def cull_redundant_hits(blast_hits):
    """ Cull out redundant hits (using Kleborate's logic). """
    #Read df 
    headers = ['hit_number','sacc','pident','score', 'length','slen','qlen','gaps','mismatch','qstart','qend','qseq'] 
    df = pd.DataFrame(blast_hits, columns=headers)
    df[["pident", "score", "length", "slen", "qlen", "qstart","qend"]] = df[["pident", "score", "length", "slen", "qlen", "qstart","qend"]].apply(pd.to_numeric) #Set these cols as numeric

    #Sort the hits from best to worst. Hit quality is defined as the product of gene coverage, identity and score.        
    #Code in Kleborate: blast_hits = sorted(blast_hits, key=lambda x: (1/(x.pident[2] * x.score[3] * x.slen[4]), x.gene_id[1])) #
    df = df.assign(sorting = lambda x: (1/x['pident'] * x['score'] * x['slen'])) #Assign hit quality to new column 
    df.sort_values(by='sorting', ascending=True) #Sort the values

    allele_hit = df.head(1).values.tolist() #Take the top hit as the allele.

    for hit in allele_hit:
        #Check which strand the sequence is on
        if hit[9] > hit[10]:
            strand="minus"
            start = hit[10]
            end = hit[9]
            my_seq = Seq(hit[11])  
            seq=my_seq.reverse_complement()

        elif hit[9] < hit[10]: 
            strand="plus"
            start = hit[9]
            end = hit[10] 
            seq = Seq(hit[11])

        #Check if the best allele hit has 100% sequence identity and coverage. Add * for inexact match, ? for incomplete coverage, and - for missing or too low coverage/pident matches.
        locus_number = hit[1].split('_')
        locus = locus_number[0]

        if hit[2] < 90 or (hit[4]/hit[5]) < 0.8: #Do not call allele number if pident is below 90 or cov is below 80
            allele_annotated = locus + "_-"
            allele_number=locus + "_-"
        elif hit[2] == 100 and (hit[4]/hit[5]) == 1: #pident is 100 and the alignment is the same length as the reference
            allele_annotated = hit[1]
            allele_number= hit[1]
        elif hit[2] < 100 and (hit[4]/hit[5]) == 1: #pident is not 100 and the alignment is the same length as the reference
            allele_annotated = (hit[1] + "*") #
            allele_number= hit[1]
        elif hit[2] < 100 and (hit[4]/hit[5]) != 1: #pident is not 100 and the alignment is not the same length as the reference 
            allele_annotated = (hit[1] + "*?")
            allele_number= hit[1]
        elif hit[2] == 100 and (hit[4]/hit[5]) != 1: #pident is 100 and the alignment is not the same length as the reference
            allele_annotated = (hit[1] + "?")
            allele_number= hit[1]

    return allele_annotated, allele_number, strand, seq, start, end 


def get_1116R_hit(std_output):
    """ Get the perc id and cov from the 1116R alignment """
    for line in std_output.splitlines():
        if line.startswith("# Identity:"):
            pidentline=line
            split_1 = pidentline.split("(") #Split by ( 
            split_1 = str(split_1[1]) #Get only the second element, i.e. the percentage part  
            split_2 = split_1.split("%") #Split by %
            pident = float(split_2[0]) #Remove the characters "%)" from the percentage.

        if line.startswith("# Gaps:"):
            covline=line
            split_1 = covline.split("(") #Split by ( 
            split_1 = str(split_1[1]) #Get only the second element, i.e. the percentage part  
            split_2 = split_1.split("%") #Split by %
            cov = (100 - float(split_2[0])) #Remove the characters "%)" from the percentage.

    return pident, cov

def run_water_alignment(df, primer_1116R, out, args, assembly_name, iteration, seq):
    """ Align query sequence against the 1116R primer """
    aseq=('asis:' + df[seq])
    try: #run water alignment on df[seq]
        water_cmd = WaterCommandline(gapopen=10, gapextend=0.5,
                                      stdout=True, auto=True,
                                       asequence=aseq, bsequence=primer_1116R)
        std_output, err_output = water_cmd()

        if args.verbose: 
            filename = out + assembly_name + "_" + str(iteration) + "_" + "water.txt" 
            logging.info("Verbose mode is on. The two mompS water alignments will be stored in file: " + filename)
            f=open(filename,"w+")
            f.writelines(std_output)

        #Then check the hit scores and assign mompS copy:
        pident, cov = get_1116R_hit(std_output) 
        if pident >= 90.0 and cov == 100.0:
            mompS_copy="mompS2"
        else:
            mompS_copy="mompS1"

    except:
        logging.exception("Error running water pairwise alignment.")
        logging.info("##############")
        logging.info("mompS2 could not be separated from mompS1 by ≥90% nucleotide identity and 100% coverage match to 1116R.")
        logging.info("Please run the command again with option --verbose to print the .water alignment and the 987bp mompS sequences to files and manually inspect them.")
        logging.info("##############")

    return mompS_copy
   

def get_output_headers():
    """ Set up headers for outfiles and stdout """
    ST_alleles = {"Strain":[],"ST":[],"flaA":[],"pilE":[],"asd":[],"mip":[],"mompS":[],"proA":[],"neuA":[]} 
    mompS_alleles = {"Strain":[],"mompS1":[],"mompS2":[]} 

    print("Strain\tST\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA\t\t\tmompS1 (Not used in the SBT scheme)")

    return ST_alleles, mompS_alleles


####################################################################################################################################
# main function
####################################################################################################################################
def main():
    args = parse_args()
    start_time = time.time()

    set_up_logging(args)
    logging.info('This is ONTmompS v2.0.0')
    logging.info('Running command: {0}'.format(' '.join(sys.argv)))
    
    #Set up inputs/outputs
    verbose=args.verbose
    db_location = check_dbdir(args) 
    sts, mompS2_ref, primer_1116R, flaA_db_file, pilE_db_file, asd_db_file, mip_db_file, mompS_db_file, proA_db_file, neuA_db_file = check_db(db_location)    
    sequence_list = check_assemblies(args)
    out, ST_outfile, mompS_outfile = check_outputs(args)


    #Get mompS + SBT alleles + ST for each assemblies in a loop. 
    ST_alleles, mompS_alleles = get_output_headers()
    for assembly_input in sequence_list:
        assembly_file, assembly_name = check_assembly(assembly_input,out)
        get_SBT(assembly_file, mompS2_ref,mompS_db_file, primer_1116R, asd_db_file, flaA_db_file, mip_db_file, neuA_db_file, pilE_db_file, proA_db_file, out, assembly_name,sts, ST_alleles,ST_outfile, mompS_alleles, mompS_outfile,args)

    #Display times
    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    logging.info("ONTmompS completed in " + str(time_mins) + ' mins.')

if __name__ == '__main__':
    main()