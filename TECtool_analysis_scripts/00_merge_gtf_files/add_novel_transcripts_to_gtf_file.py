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

try:
  import os
except(Exception):
  raise("[ERROR] os was not imported properly. Exiting.")
  sys.exit(-1)


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():

    """ Main function """

    __doc__ = "Add novel transcripts to a given annotation (GTF file)."

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument("-t",
                        "--tectool_installation_dir",
                        dest = "tectool_installation_dir",
                        help = "The pathway to the directory to which TECtool was installed. [REQUIRED]")

    parser.add_argument("--list_of_gtf_files",
                        dest = "list_of_gtf_files",
                        help = "[INPUT] Tab separated file (TSV-format) that contains in the first column pathways to GTF files that contain novel transcripts (having intronic poly(A) sites) that will be add to the ref-gtf file. NOTE: The script will not correct for novel gene boarders, since intronic poly(A) sites will not change gene boarders. Thus, do not use this script for adding transcripts that are not within defined gene boarders.",
                        required = True,
                        metavar = "FILE")

    parser.add_argument("--out-dir",
                        dest = "out",
                        help = "[OUTPUT] The directory that will be used to write files, including the enriched annotation (GTF file).",
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


    # _____________________________________________________________________________
    # -----------------------------------------------------------------------------
    # import our own modules
    # -----------------------------------------------------------------------------
    # import tectool classes
    sys.path.insert(0, options.tectool_installation_dir)

    from gene_structure import Exon
    from gene_structure import Transcript
    from gene_structure import Gene
    from aln_analysis import DetailedAlignment
    from aln_analysis import SplitEvent
    from aln_analysis import AnalysisUnit
    from aln_analysis import FeatureCounts
    from annotation import Annotation
    #from machine_learning import MachineLearningUnit

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Create output directories
    # -------------------------------------------------------------------------
    if not os.path.exists(options.out):
        os.makedirs(options.out)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Go over all GTF files in the list_of_gtf_files and add novel transcripts
    # to the annotation.
    # -------------------------------------------------------------------------

    # create an annotation
    annotation = Annotation(annotation_id="merged_annotation",
                            tmp=options.out)

    # open the file
    list_of_gtf_files_handle = open(options.list_of_gtf_files,"r")

    # create an annotation from each GTF file and then extend the 
    # reference annotation (ref-gtf) by this annotation.
    for line in list_of_gtf_files_handle:

        # parse the current GTF file
        if options.verbose:
            sys.stdout.write((80*"_" + "\n"))
            sys.stdout.write((80*"-" + "\n"))

        gtf_file_path = str(line.strip().split("\t")[0])
        sys.stdout.write("Reading annotation file:\t%s\n" % gtf_file_path)
	curr_annotation = Annotation(annotation_id=gtf_file_path,
                                     tmp=options.out)
	curr_annotation.parse(gtf_file_path, verbose=options.verbose)

        if options.verbose:
            sys.stdout.write((80*"-" + "\n"))

        # extend the annotation by the current annotation
        annotation.extend(annotation=curr_annotation, 
                          verbose=options.verbose)

    # close the filehandle
    list_of_gtf_files_handle.close()

    # write the output GTF file
    output_gtf_file_path = os.path.join(options.out, "merged_annotation.gtf")
    output_gtf_file_handle = open(output_gtf_file_path, 'w')
    annotation.write_gtf(output_file=output_gtf_file_handle, 
                         with_CDS = True)
    output_gtf_file_handle.close()

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
