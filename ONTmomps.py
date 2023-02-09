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
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + 'v1.2.0')

    #Argsgroups
    input_args = parser.add_argument_group('Input options (required)')
    optional_args = parser.add_argument_group('Optional flags')
    output_args = parser.add_argument_group('Output options')

    #Arguments 
    input_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=True, help='FASTA file(s) for assemblies (*.fasta)')

    optional_args.add_argument('-d', '--database_folder', type=str, required=False, help='Provide a path to database location if different than that provided by this tool.')

    #TODO: optional_args.add_argument('--store_mompS_alleles', action='store_true', required=False, help='Print mompS alleles to files named {assembly}_{allele}.fna.')

    #TODO: optional_args.add_argument('--store_novel_alleles', action='store_true', required=False, help='Print novel alleles to files named {assembly}_{allele}.fna.')

    #TODO: optional_args.add_argument('--store_all_alleles', action='store_true', required=False, help='Print all alleles (7 genes in SBT scheme + mompS1) to files named {assembly}_{allele}.fna. ')

    #TODO: optional_args.add_argument('--verbose', action='store_true', required=False,  help='Keep intermediate files for debugging.') #Set this to store the water alignment.

    output_args.add_argument('--outfilename', type=str, required=False, default='./LpST_ONTmompS.tsv', help='Output filename for STs. Default: ./LpST_ONTmompS.tsv')

    output_args.add_argument('-outdir', type=str, required=False, default='.', help='Output directory to store novel alleles in. Default is current working directory')

    return parser.parse_args()

#####################################################################
### Defs for checking verions, arguments, databases
#####################################################################

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
    """Check that database location exists and contains required files, and save to variable"""
    logging.info("Checking database for required files and storing variables.")

    #Assign variable names to files
    sts = db_location + 'lpneumophila.txt'
    mompS2_ref = db_location + 'mompS2_ref.fna'
    primer_1116R = db_location + '1116R_rev.fasta'
    flaA_db_file = db_location + 'flaA.fna'
    pilE_db_file = db_location + 'pilE.fna'
    asd_db_file = db_location + 'asd.fna'
    mip_db_file = db_location + 'mip.fna'
    mompS_db_file = db_location + 'mompS.fna'
    proA_db_file = db_location + 'proA.fna'
    neuA_db_file = db_location + 'neuA.fna'

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

    print("Database location: " + db_location)

    return sts, mompS2_ref, primer_1116R, flaA_db_file, pilE_db_file, asd_db_file, mip_db_file, mompS_db_file, proA_db_file, neuA_db_file

def load_st_db(sts): #Uses same logic as Kleborate
    """Load the lpneumophila.txt ST database"""
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
    #Outdir for temp and storing
    if args.outdir:
        out = args.outdir
        if out[-1] != '/':
            out = out + '/'
    else:
        out = "./"

    if not os.path.isdir(out):
        os.makedirs(out)
    else:
        sys.exit('Error: Output folder already exists.')

    #Outfile for STs 
    if args.outfilename:
        outfilename = args.outfilename
    else:
        outfilename = "./LpST_ONTmompS.tsv"

    print("Output STs will be written to: " + os.path.abspath(outfilename))
    print("Any extra files will be written to: " + os.path.abspath(out) + "/") #TODO: add arg for keeping extra files and make temp_dir instead of this specific name, but give option for specific name.
    return out, outfilename


def check_assemblies(args):
    #Set up a list of all input files:
    fasta=args.assemblies
    #Create list of input sequences 
    genomes = []
    sequence_list = []
    for sequence in fasta: 
        if sequence.find('.fasta') != -1 and sequence not in sequence_list:
            sequence_list.append(sequence)
    return sequence_list

def check_assembly(assembly_file,out):
    """Check that assembly is a file and symlink"""
    logging.info("Checking assembly input.")

    try:
        if not os.path.isfile(assembly_file): #TODO: Check that file is a fasta file >
            raise
    except:
        logging.exception("Assembly file not located.")
        sys.exit("Unable to find assembly file. Please check path again.")

    try:
        source=os.path.abspath(assembly_file)
        base_name=os.path.basename(assembly_file)
        destination=out+base_name
        os.symlink(source, destination) #TODO: remove symlink, or add to temp dir that is deleted in the end
        assembly_file=destination
        assembly_name=os.path.splitext(base_name)[0]

    except:
        logging.exception("Symlink failed.")
        sys.exit("Symlink failed. Please check your input.")

    return assembly_file, assembly_name

#####################################################################
### Defs for BLASTn and parsing BLASTn output 
#####################################################################
def run_blastn(subject, query):
    """Run blast of mompS2 ref against assembly"""
    logging.info("Running blast of mompS2 ref against assembly.")
    #Run BLASTn to find mompS hits. TODO: Store positions, strand and sequence.
    #Uses the same BLASTn logic as Kleborate's MLSTblast.py
    try:
        mompS_blastn_command = NcbiblastnCommandline(cmd='blastn',
                                  task='blastn',
                                  query=query, 
                                  subject=subject, 
                                  dust='no',
                                  evalue='1E-20',
                                  word_size='32',
                                  max_target_seqs='10000',
                                  perc_identity='90',
                                  outfmt='"6 sacc pident slen qlen length score gaps mismatch bitscore sseq qseq sstrand sstart send qacc qstart qend qframe"')
        blastn_output = mompS_blastn_command()[0].strip()
        #TODO: Should change to?: stdout, stderr = mompS_blastn_command()

    except:
        logging.exception("Error in blast assembly file to mompS2_ref.")
        sys.exit("Error in blast assembly file to mompS2_ref in db.")

    return blastn_output

def read_blastn_output(blastn_output):
    #Read blastn output as a pandas df
    try:
        headers = ['sacc','pident','slen','qlen','length','score','gaps','mismatch','bitscore','sseq','qseq','sstrand','sstart','send','qacc','qstart','qend','qframe'] #Some not used - might be removed.

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
    """Get the sequence, start and stop positions from blast output. Reverse complement if on minus strand."""
    logging.info("Getting sequence, start and stop positions of mompS hits in the assembly file.")
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

    return sseq, sstart, send #TODO: could add strand to this if we want to report it?

#####################################################################
### Defs for identifying mompS and SBT 
#####################################################################

def get_mompS_blastn_hits(assembly_file, mompS2_ref,mompS_db_file,primer_1116R, out, assembly_name):
    try:
        blastn_output = run_blastn(assembly_file, mompS2_ref)
        df = read_blastn_output(blastn_output)

        #TODO:If the sequence is on the reverse strand we should flip it when storing? Check for minus 
        mompS_blastn_hits = []
        hit_number=0
        for index, row in df.iterrows():
            hit_number+=1
            sseq, sstart, send = get_sseq_start_stop_positions(row,strand='sstrand',seq='sseq',start='sstart',end='send')
            mompS_copy = run_water_alignment(row, primer_1116R, seq='sseq')
            to_append = [mompS_copy, sseq, sstart, send]
            mompS_blastn_hits.append(to_append)
            #TODO: Count how many hits in the list. Should be exactly two - 1 mompS1 and 1 mompS2. Add warning if not.

        #Assign mompS copies and store sequence to respective files:
        for hit in mompS_blastn_hits:
            if hit[0]=='mompS2':
                mompS2_copy=hit[0]
                mompS2_sequence=hit[1]
                run_command(['echo ">',mompS2_copy,'__',assembly_name,'" >>  ',out,mompS2_copy,'_sequence_', assembly_name,'.fasta ; echo ', str(mompS2_sequence),' >>  ',out,mompS2_copy,'_sequence_', assembly_name,'.fasta'  ], shell=True) #TODO: Add start and stop positions to the fasta header.
                #TODO: if verbose option, print these to screen.
            elif hit[0]=='mompS1':
                mompS1_copy=hit[0]
                mompS1_sequence=hit[1]
                run_command(['echo ">',mompS1_copy,'__',assembly_name,'" >>  ',out,mompS1_copy,'_sequence_', assembly_name,'.fasta ; echo ', str(mompS1_sequence),' >>  ',out,mompS1_copy,'_sequence_', assembly_name,'.fasta'  ], shell=True) #TODO: Add start and stop positions to the fasta header. 
             
     except:
        logging.exception("Error in blast assembly file to mompS2_ref.")
        sys.exit("Error in blast assembly file to mompS2_ref in db.")

    return mompS2_copy, mompS2_sequence, mompS1_copy, mompS1_sequence


def get_SBT(assembly_file, mompS2_ref,mompS_db_file, primer_1116R, asd_db_file, flaA_db_file, mip_db_file, neuA_db_file, pilE_db_file, proA_db_file, out, assembly_name, sts, ST_alleles):
    """Read all alleles and assign ST """
    mompS2_copy, mompS2_sequence, mompS1_copy, mompS1_sequence = get_mompS_blastn_hits(assembly_file, mompS2_ref,mompS_db_file,primer_1116R, out, assembly_name) #TODO: maybe add mompS start and stop back to what is returned


    #Read mompS sequences for each copy:    
    mompS2_sequence = (out + 'mompS2_sequence_' + assembly_name + '.fasta') #write same header as for rest?
    mompS1_sequence = (out + 'mompS1_sequence_' + assembly_name + '.fasta')

    #For each of the 7 SBT genes (+ mompS1), get the allele number and check if it is *, ? or -:
    flaA_allele_annotated, flaA_allele_number, flaA_allele = get_SBT_allele_blastn_hits(assembly_file, flaA_db_file)
    pilE_allele_annotated, pilE_allele_number, pilE_allele = get_SBT_allele_blastn_hits(assembly_file, pilE_db_file)
    asd_allele_annotated, asd_allele_number, asd_allele  = get_SBT_allele_blastn_hits(assembly_file, asd_db_file)
    mip_allele_annotated, mip_allele_number, mip_allele = get_SBT_allele_blastn_hits(assembly_file, mip_db_file)
    mompS2_allele_annotated,  mompS2_allele_number, mompS2_allele = get_SBT_allele_blastn_hits(mompS2_sequence, mompS_db_file) #TODO: read variable as file? add seq to save to file
    proA_allele_annotated, proA_allele_number, proA_allele = get_SBT_allele_blastn_hits(assembly_file, proA_db_file) #could 
    pilE_allele_annotated, pilE_allele_number, pilE_allele = get_SBT_allele_blastn_hits(assembly_file, pilE_db_file)
    neuA_allele_annotated, neuA_allele_number, neuA_allele = get_SBT_allele_blastn_hits(assembly_file, neuA_db_file)

    mompS1_allele_annotated, mompS1_allele_number, mompS1_allele = get_SBT_allele_blastn_hits(mompS1_sequence, mompS_db_file) #TODO: read variable as file?

    #Create a combination with all of the ST and allele numbers
    allele_combinations = {"flaA":flaA_allele,"pilE":pilE_allele,"asd":asd_allele,"mip":mip_allele,"mompS":mompS2_allele,"proA":proA_allele,"neuA":neuA_allele}
    annotated_allele_combinations = {"flaA":flaA_allele_annotated,"pilE":pilE_allele_annotated,"asd":asd_allele_annotated,"mip":mip_allele_annotated,"mompS":mompS2_allele_annotated,"proA":proA_allele_annotated,"neuA":neuA_allele_annotated}
    #TODO: see if we can combine these dictionaries

    #Assign ST
    ST_annotated = get_closest_ST(sts, allele_combinations, annotated_allele_combinations)

    append_to_ST_report(assembly_name, ST_annotated, flaA_allele_annotated, pilE_allele_annotated, asd_allele_annotated, mip_allele_annotated, mompS2_allele_annotated, proA_allele_annotated, neuA_allele_annotated, ST_alleles)


def append_to_ST_report(assembly_name, ST_annotated, flaA_allele_annotated, pilE_allele_annotated, asd_allele_annotated, mip_allele_annotated, mompS2_allele_annotated, proA_allele_annotated, neuA_allele_annotated, ST_alleles):
 
    ST_alleles["Strain"].append(assembly_name)
    ST_alleles["ST"].append(ST_annotated)
    ST_alleles["flaA"].append(flaA_allele_annotated)
    ST_alleles["pilE"].append(pilE_allele_annotated)
    ST_alleles["asd"].append(asd_allele_annotated)	
    ST_alleles["mip"].append(mip_allele_annotated)	
    ST_alleles["mompS"].append(mompS2_allele_annotated)	
    ST_alleles["proA"].append(proA_allele_annotated)	
    ST_alleles["neuA"].append(neuA_allele_annotated)

    #Print results temp to outfile 
    print(assembly_name + "\t" + ST_annotated + "\t" + flaA_allele_annotated + "\t" + pilE_allele_annotated + "\t" + asd_allele_annotated + "\t" + mip_allele_annotated + "\t" + mompS2_allele_annotated + "\t" + proA_allele_annotated + "\t" + neuA_allele_annotated)


def get_closest_ST(sts, allele_combinations, annotated_allele_combinations): #Adapted from Kleborate
    """Get the closest matching ST if <3 alleles are imprecise matches. """
    header, st_names, alleles_to_st = load_st_db(sts) #Load the lpneumophila.txt database

        best_allele_numbers_hits = []
    best_allele_numbers_annotated = []
    mismatch_loci, missing_loci = 0, 0
    loci=["flaA","pilE","asd","mip","mompS","proA","neuA"] 

    #For each locus, count if mismatch, and append the unannotated allele number to best_allele_numbers_hits
    for locus in loci:
        allele_number = allele_combinations[locus] #Get the unannotated  allele number from the dictionary 
        allele = annotated_allele_combinations[locus] #Get the annotated allele number from the dictionary 

        if '*?' in allele:
            mismatch_loci += 1 #Count number of genes (out of the 7) with a mismatch or incomplete coverage
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

def get_closest_locus_variant(best_allele_numbers_hits, best_allele_numbers_annotated, alleles_to_st): #Adapted from Kleborate
    best_allele_numbers_annotated = list(best_allele_numbers_annotated) 
    closest = []
    closest_alleles = {}   # key = st, value = list
    min_dist = len(best_allele_numbers_hits)  # number mismatching loci, ignoring SNPs

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


def get_SBT_allele_blastn_hits(query, subject): #TODO: Must check that these are the correct way around for db and seq
    """Run blastn of contig against locus database"""
    logging.info("Running blastn of contig against locus database.")

    #Run sequence against the locus database to find closest matching allele
    try:
        blastn_output = run_blastn(subject, query) 
        df = read_blastn_output(blastn_output)

        if df.empty: #If there are no BLASTn hits, the annotation for the allele should be "-"
            annotated_allele = "-"
            allele_number = "-"
            allele = "-"
        if not df.empty:
            blastn_hits = []
            hit_number=0

            for index, row in df.iterrows(): #For each blast hit, do:
                hit_number+=1
                #qseq, qstart, qend = get_sseq_start_stop_positions(row,strand='sstrand',seq='qseq',start='qstart',end='qend') #can we do this for only the one hit? 
                to_append = [hit_number, row['sacc'], row['pident'], row['score'], row['length'], row['slen'], row['qlen'], row['gaps'], row['mismatch'], row['qstart'], row['qend'], row['qseq']] #TODO: Bitscore? #TODO Add strand, and switch start + end if minus strand
                #To remember: Query is assembly. Seq/db is mlst db. 
                blastn_hits.append(to_append)

                #Keep only the line with the best hit
                #TODO: Make sure this logic is correct - e.g. overlapping hits?
            allele_annotated, allele_number, assembly_strand, assembly_seq, assembly_start, assembly_end   = cull_redundant_hits(blastn_hits) #TODO: ALSO RETURN SEQUENCE AND START STOP POSITIONS?
    
            #Get locus, allele number and annotated allele number
            locus_number = allele_number.split('_')
            allele = locus_number[1]
            locus_number = allele_annotated.split('_')
            annotated_allele = locus_number[1]

        return annotated_allele, allele_number, allele

    except:
        logging.exception("Error in blastn mompS db against extracted contig.")
        sys.exit("Error in blastn mompS db against extracted contig.")


def cull_redundant_hits(blast_hits): #Adapted from Kleborate  
    """
    Cull out redundant hits here (essentially implementing BLAST's -culling_limit 1 feature but with Kleborate's own logic).
    """
    # Sort the hits from best to worst. Hit quality is defined as the product of gene coverage, identity and score.
    headers = ['hit_number','sacc','pident','score', 'length','slen','qlen','gaps','mismatch','qstart','qend','qseq'] 
    df = pd.DataFrame(blast_hits, columns=headers)
    df[["pident", "score", "length", "slen", "qlen", "qstart","qend"]] = df[["pident", "score", "length", "slen", "qlen", "qstart","qend"]].apply(pd.to_numeric) #Set these cols as numeric

        
    #Assign sorting score using Kleborate's lamba function for this, and then sort by the score:
    #From Kleborate: blast_hits = sorted(blast_hits, key=lambda x: (1/(x.pident[2] * x.score[3] * x.slen[4]), x.gene_id[1])) #
    df = df.assign(sorting = lambda x: (1/x['pident'] * x['score'] * x['slen']))
    df.sort_values(by='sorting', ascending=True) #Sort the values

    allele_hit = df.head(1).values.tolist() #Take the top hit as the allele.
    # NB TODO add more checks for overlap etc before doing this

    for hit in allele_hit: #should only be one hit, there is probably a neater way to look through these, but it works for now
        #Check which strand the sequence is on
        if hit[9] > hit[10]:
            strand="minus"
            start = hit[10]
            end = hit[9]
            my_seq = Seq(hit[11]) #TODO: Make sure that we are storing the correct sequence - i.e. the sequence from the assembly, not from the database.
            seq=my_seq.reverse_complement()

        elif hit[9] < hit[10]: 
            strand="plus"
            start = hit[9]
            end = hit[10] 
            seq = Seq(hit[11])

        #Check if the best allele hit has 100% sequence identity and coverage. Add * for inexact match ? for incomplete coverage, and "-" for missing or too low coverage/pident matches.

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
        else:
            print("Error") #TODO: Add error or warning if not work


    return allele_annotated, allele_number, strand, seq, start, end 

#TODO: Adapt features from Kleborate.
#filtered_blast_hits = []
#for h in blast_hits:
#    if not overlapping(h, filtered_blast_hits):
#        filtered_blast_hits.append(h)

#return filtered_blast_hits

def overlapping(hit, existing_hits):
    # Only consider hits in the same reading frame.
    existing_hits = [h for h in existing_hits if
                     h.strand == hit.strand and h.frame == hit.frame and
                     h.contig_name == hit.contig_name]

    for existing_hit in existing_hits:
        if hits_overlap(hit, existing_hit):
            return True

    return False


def hits_overlap(a, b):
    if a.contig_start <= b.contig_end and b.contig_start <= a.contig_end:  # There is some overlap
        allowed_overlap = 50
        overlap_size = len(range(max(a.contig_start, b.contig_start),
                                 min(a.contig_end, b.contig_end) + 1))
        return overlap_size > allowed_overlap
    else:
        return False
#END TODO Adapt features from Kleborate.

def get_1116R_hit(std_output):
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

def run_water_alignment(df, primer_1116R,seq):
    """Run pairwise alignment"""
    logging.info("Running water alignment on contig fasta files with 1116R primer.")

    aseq=('asis:' + df[seq])
    try: #run water alignment on df[seq]
        water_cmd = WaterCommandline(gapopen=10, gapextend=0.5,
                                      stdout=True, auto=True,
                                       asequence=aseq, bsequence=primer_1116R)
        std_output, err_output = water_cmd()

        #print(std_output) #Type string - TODO: print in verbose mode only
        #Then check the hit scores and assign mompS copy:
        pident, cov = get_1116R_hit(std_output) 
        if pident >= 90.0 and cov == 100.0: #TODO: add option to specify minimum %ID and %COV for 1116R
            mompS_copy="mompS2"
        else:
            mompS_copy="mompS1"
        #TODO: Check if 1116R is always on rev, and if we should flip?

    except: #TODO update this exception
        logging.exception("Error running water pairwise alignment.")
        sys.exit("Error running water pairwise alignment.")
        print("#######NB#######")
        print("mompS2 could not be separated from mompS1 by exact match to 1116R.")
        print("Please check .water files to manually deduce.")
        print("#######NB#######")
        logging.exception("Error in reading blast output.")
        sys.exit("Error in reading momps blast output.")

    return mompS_copy
   
def write_output(ST_alleles,outfilename):
    """Write output to files"""
    #Convert dictionary to dataframe    
    ST = pd.DataFrame.from_dict(ST_alleles, orient='index').T 
    #Write to output file
    ST.to_csv(outfilename,index=False, sep = "\t")
    print("STs written to: " + outfilename)
    logging.info("STs written to: " + outfilename)


####################################################################################################################################
# main function
####################################################################################################################################
def main():
    args = parse_args()
    threads = str(args.threads) 
    start_time = time.time()

    print("This is ONTmompS v1.2.0")

    #TODO: Set all of this in a setup def.
    db_location = check_dbdir(args) 
    sts, mompS2_ref, primer_1116R, flaA_db_file, pilE_db_file, asd_db_file, mip_db_file, mompS_db_file, proA_db_file, neuA_db_file = check_db(db_location) #TODO: move these into the necessary defs instead of here
    sequence_list = check_assemblies(args)
    out, outfilename = check_outputs(args)

    #Logging
    #TODO: only log if --log flag has been activated. Set to always log, but only to file if specified?    
    logging.basicConfig(filename=out+"ONTmompS.log", format='%(asctime)s %(message)s', filemode='w', level=logging.DEBUG) 
    logging.info("Started ONTmompS.")

    #Run SBT + mompS on all assemblies in loop.
    ST_alleles = {"Strain":[],"ST":[],"flaA":[],"pilE":[],"asd":[],"mip":[],"mompS":[],"proA":[],"neuA":[]};
    print("Strain\tST\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA")

    for assembly_input in sequence_list:
        assembly_file, assembly_name = check_assembly(assembly_input,out)
        get_SBT(assembly_file, mompS2_ref,mompS_db_file, primer_1116R, asd_db_file, flaA_db_file, mip_db_file, neuA_db_file, pilE_db_file, proA_db_file, out, assembly_name,sts, ST_alleles)

    write_output(ST_alleles, outfilename)

    #Display times
    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    logging.info("ONTmompS completed.")
    print("ONTmompS finished in " + str(time_mins) + ' mins.')


# call main function
if __name__ == '__main__':
    main()
