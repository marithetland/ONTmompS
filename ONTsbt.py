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
            df = pd.DataFrame(pd.read_csv(file, sep="\t", names=["Query_ID", "Subject_ID", "Perc_match", "Length", "Num_mismatches", 
"Num_gaps", "Query_start", "Query_end", "Subject_start", "Subject_end", "E_value", "Bitscore"], engine='python'))
            for index, row in enumerate(df.iterrows()):
                length = df.loc[index, 'Length']
                momps = df.loc[index, 'Subject_ID']
            key_blast_momps.write(str(contig) + '\t' + str(momps) + '\t' + str(length) + '\n')

        key_blast_momps.close()

    except:
        logging.exception("Error in reading blast output.")
        sys.exit("Error in reading momps blast output.")

    return key_blast_momps_file
