# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

try:
  import HTSeq
except(Exception):
  raise("[ERROR] HTSeq was not imported properly. Exiting.")
  sys.exit(-1)

try:
  from argparse import ArgumentParser, RawTextHelpFormatter, FileType
except(Exception):
  raise("[ERROR] argparse was not imported properly. Exiting.")
  sys.exit(-1)

try:
  import sys
except(Exception):
  raise("[ERROR] sys was not imported properly. Exiting.")
  sys.exit(-1)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():

    """ Main function """

    __doc__ = "Filter gtf file based on the transcript support level."

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)


    parser.add_argument("--ref-gtf",
                      dest = "gtf",
                      help = "[INPUT] Annotation file in GTF format that should be extended by novel transcripts occuring in a list of gtf files.",
                      required = True,
                      metavar = "FILE")

    parser.add_argument("--list_of_gtf_files",
                      dest = "gtfs_to_add",
                      help = "[INPUT] List of GTF files that contain novel transcripts (having intronic poly(A) sites) that will be add to the ref-gtf file. NOTE: The script will not correct for novel gene boarders, since intronic poly(A) sites will not change gene boarders. Thus, do not use this script for adding transcripts that are not within defined gene boarders.",
                      required = True,
                      metavar = "FILE")

    parser.add_argument("--out",
                      dest = "out",
                      help = "[OUTPUT] The GTF annotation file (ref-gtf) that is enriched by enriched extended GTF file output file",
                      required = True,
                      metavar="FILE")

    parser.add_argument("-v",
                      "--verbose",
                      action = "store_true",
                      dest = "verbose",
                      default = False,
                      required = False,
                      help = "Verbose")

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if options.verbose:
        sys.stdout.write("Parsing gtf file and keep transcripts with support level <= %s\n" % (str(options.support_level)))

    # dictionary of accepted transcript ids
    transcripts_dict = dict()

    # dictionary of accepted gene ids
    genes_dict = dict()

    # dictionary of min start of accepted genes ids
    genes_dict_start = dict()

    # dictionary of max end of accepted gene ids
    genes_dict_end = dict()

    # parse gtf file
    gtf_file_with_support_level = HTSeq.GFF_Reader(options.gtf)

    for gtf_line in gtf_file_with_support_level:

        if gtf_line.type == 'transcript':

            if "transcript_support_level" in gtf_line.attr:

                support_level = str(gtf_line.attr['transcript_support_level'])

                # skip NAs
                if "NA" in support_level:
                    continue
                else:

                    support_level = support_level.strip().split(" ")[0]

                    # keep the  requested support levels
                    if int(support_level) <= int(options.support_level):

                        transcript_id = gtf_line.attr['transcript_id']
                        gene_id = gtf_line.attr['gene_id']

                        transcripts_dict[transcript_id] = transcript_id
                        genes_dict[gene_id] = gene_id

                        # define new gene borders
                        if options.fix_gene_coordinates:

                            transcript_start = int(gtf_line.iv.start)
                            transcript_end = int(gtf_line.iv.end)

                            # define gene starts
                            if gene_id not in genes_dict_start:
                                genes_dict_start[gene_id] = transcript_start
                            else:
                                if transcript_start < genes_dict_start[gene_id]:
                                    genes_dict_start[gene_id] =  transcript_start

                            # define gene ends
                            if gene_id not in genes_dict_end:
                                genes_dict_end[gene_id] = transcript_end
                            else:
                                if transcript_end > genes_dict_end[gene_id]:
                                    genes_dict_end[gene_id] =  transcript_end


    if options.verbose:
        sys.stdout.write("Re-parse gtf file and write out selected genes and transcripts\n")

    w = open(options.out, 'w')

    # parse again the gtf file and keep
    # only the transcripts and genes that
    # were stored in the previous alignment files
    for gtf_line in gtf_file_with_support_level:

        if gtf_line.type == 'gene':

            if gtf_line.attr['gene_id'] in genes_dict:

                if options.fix_gene_coordinates:
                    gtf_line.iv.start = genes_dict_start[gtf_line.attr['gene_id']]
                    gtf_line.iv.end = genes_dict_end[gtf_line.attr['gene_id']]
                    w.write(gtf_line.get_gff_line())
                else:
                    w.write(gtf_line.get_gff_line())

        if gtf_line.type == 'transcript':

            if gtf_line.attr['transcript_id'] in transcripts_dict and gtf_line.attr['gene_id'] in genes_dict:
                w.write(gtf_line.get_gff_line())

        elif gtf_line.type == 'exon':

            if gtf_line.attr['transcript_id'] in transcripts_dict and gtf_line.attr['gene_id'] in genes_dict:
                w.write(gtf_line.get_gff_line())

        elif gtf_line.type == 'CDS':

            if gtf_line.attr['transcript_id'] in transcripts_dict and gtf_line.attr['gene_id'] in genes_dict:
                w.write(gtf_line.get_gff_line())

        elif gtf_line.type == 'start_codon':

            if gtf_line.attr['transcript_id'] in transcripts_dict and gtf_line.attr['gene_id'] in genes_dict:
                w.write(gtf_line.get_gff_line())

        elif gtf_line.type == 'stop_codon':

            if gtf_line.attr['transcript_id'] in transcripts_dict and gtf_line.attr['gene_id'] in genes_dict:
                w.write(gtf_line.get_gff_line())

    w.close()

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)
