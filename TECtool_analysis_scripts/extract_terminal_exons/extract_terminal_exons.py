# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

try:
    import os
except(Exception):
    raise("[ERROR] os was not imported properly. Exiting.")
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
    from pybedtools import BedTool
except(Exception):
    raise("[ERROR] pybedtools was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from pyfasta import Fasta
except(Exception):
    raise("[ERROR] pyfasta was not imported properly. Exiting.")
    sys.exit(-1)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():

    """ Main function """

    __doc__ = "Keep annotated terminal exons from gtf file"

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument("--tectool_installation_dir",
                        dest = "tectool_installation_dir",
                        required = True,
                        help = "The pathway to the directory to which TECtool was installed.")

    parser.add_argument("--gtf",
                        dest = "gtf",
                        help = "Annotation file",
                        required = True,
                        metavar = "FILE")

    parser.add_argument("--genome",
                        dest="genome",
                        help = "Genome sequence",
                        required=True,
                        metavar = "FILE")

    parser.add_argument("--classified_as_terminal_exons_file",
                        dest="classified_as_terminal_exons_file",
                        help="Classified as terminal exons from TECtool",
                        required=True,
                        metavar="FILE")

    parser.add_argument("--out",
                        dest = "out",
                        help = "Output directory",
                        required = True)

    parser.add_argument("-v",
                        "--verbose",
                        action = "store_true",
                        dest = "verbose",
                        default = False,
                        required = False,
                        help = "Verbose")

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # DEBUG mode
    #tectool_installation_dir = "/scicore/home/zavolan/gypas/projects/TECtool/tectool"
    #gtf = "/scicore/home/zavolan/GROUP/PolyA/TECtool/resources/hg38/ENSEMBL/version87/Homo_sapiens.GRCh38.87.chr.support_level_5.correct_gene_coordinates.gtf"
    #out = "this"

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

    if options.verbose:
        sys.stdout.write("")

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Create output directories
    # -------------------------------------------------------------------------
    if not os.path.exists(options.out):
        os.makedirs(options.out)

    # Read the genome
    genome = Fasta(options.genome)

    # create annotation object
    annotation = Annotation(annotation_id=options.gtf,
                                     tmp=options.out)
    # parse the annotation
    annotation.parse(options.gtf, verbose=options.verbose)

    # find terminal exons and write them in a bed file (contains duplicates)
    annotation.write_all_terminal_exons_as_bed(os.path.join(options.out, "all_terminal_exons.bed"))

    # find all first and intermediate exons and write them in a bed file
    # (contains duplicates)
    annotation.write_all_first_and_intermediate_exons_as_bed(os.path.join(options.out, "all_first_and_intermediate_exons.bed"))

    # remove duplicate entries
    os.system("sort -u %s | sort -k 1,1 -k2,2n > %s" % (os.path.join(options.out, "all_terminal_exons.bed"),
                                                        os.path.join(options.out, "all_terminal_exons_filtered_and_sorted.bed")))

    # remove duplicate entries
    os.system("sort -u %s | sort -k 1,1 -k2,2n > %s" % (os.path.join(options.out, "all_first_and_intermediate_exons.bed"),
                                                        os.path.join(options.out, "all_first_and_intermediate_exons_filtered_and_sorted.bed")))

    # read the two files with pybedtools
    terminal_exons = BedTool(os.path.join(options.out, "all_terminal_exons_filtered_and_sorted.bed"))
    other_exons = BedTool(os.path.join(options.out, "all_first_and_intermediate_exons_filtered_and_sorted.bed"))

    # keep terminal exons that do not overlap with intermediate
    # bedtools intersect -a all_terminal_exons_filtered_and_sorted.bed -b
    # all_first_and_intermediate_exons_filtered_and_sorted.bed -s -v
    non_overlapping_terminal_exons = terminal_exons.intersect(other_exons, s=True, v=True)
    non_overlapping_terminal_exons.saveas(os.path.join(options.out, "non_overlapping_terminal_exons.bed"))

    # Extract the exon sequences
    w = open(os.path.join(options.out, "non_overlapping_terminal_exons.fa"),'w')

    for terminal_exon in non_overlapping_terminal_exons:

        terminal_exon_header=":".join([">"+str(terminal_exon.chrom),
                                       str(int(terminal_exon.start)+1),
                                       str(terminal_exon.end),
                                       str(terminal_exon.name),
                                       str(terminal_exon.strand),
                                       "annotated"])

        terminal_exon_seq = genome.sequence({'chr': terminal_exon.chrom,
                                             'start': int(terminal_exon.start),
                                             'stop': int(terminal_exon.end),
                                             'strand':terminal_exon.strand},
                                             one_based=False)
        w.write(terminal_exon_header+"\n")
        w.write(terminal_exon_seq+"\n")

    w.close()

    # Read file with novel exons and create a bed file
    w = open(os.path.join(options.out, "novel_exons.bed"), 'w')
    with open(options.classified_as_terminal_exons_file) as fp:
        next(fp) # skip header line
        for line in fp:
            line_sp = line.strip().split("\t")
            region = line_sp[0]
            gene_id = line_sp[1]
            chrom, start, end, strand = region.split(":")
            bed_line = "\t".join([chrom,
                                  str(int(start)-1),
                                  str(end),
                                  gene_id,
                                  ".",
                                  strand])
            w.write(bed_line+"\n")
    w.close()

    # extract sequences
    novel_exons = BedTool(os.path.join(options.out, "novel_exons.bed"))
    w = open(os.path.join(options.out, "novel_exons.fa"), 'w')
    for novel_exon in novel_exons:

        novel_exon_header = ":".join([">"+str(novel_exon.chrom),
                                      str(int(novel_exon.start)+1),
                                      str(novel_exon.end),
                                      str(novel_exon.name),
                                      str(novel_exon.strand),
                                      "novel"])

        novel_exon_seq = genome.sequence({'chr': novel_exon.chrom,
                                          'start': int(novel_exon.start),
                                          'stop': int(novel_exon.end),
                                          'strand': novel_exon.strand},
                                          one_based=False)

        w.write(novel_exon_header+"\n")
        w.write(novel_exon_seq+"\n")

    w.close()

    # merge the two files
    os.system("cat %s %s > %s" % (os.path.join(options.out, "non_overlapping_terminal_exons.fa"),
                                  os.path.join(options.out, "novel_exons.fa"),
                                  os.path.join(options.out, "novel_and_annotated_terminal_exons.fa")))

    # create mapping file transcript/exon id to gene id
    w = open(os.path.join(options.out, "novel_and_annotated_terminal_exons_mapping_list.tsv"), 'w')
    with open(os.path.join(options.out, "novel_and_annotated_terminal_exons.fa")) as fp:
        for line in fp:
            if line[0] == ">":
                chrom, start, end, name, strand, typ = line.strip().split(":")
                w.write("\t".join([line.strip()[1:], name+"\n"]))
    w.close()

    # Write the 5p border only once per gene
    w = open(os.path.join(options.out, "non_overlapping_terminal_exons_with_same_3p_site.bed"), 'w')
    for non_overlapping_terminal_exon in non_overlapping_terminal_exons:
        if non_overlapping_terminal_exon.strand == "+":
            w.write("\t".join([non_overlapping_terminal_exon.chrom, 
                               str(non_overlapping_terminal_exon.start), 
                               str(int(non_overlapping_terminal_exon.start)+100),
                               non_overlapping_terminal_exon.name,
                               ".",
                               non_overlapping_terminal_exon.strand+"\n"]))
        elif non_overlapping_terminal_exon.strand == "-":
            w.write("\t".join([non_overlapping_terminal_exon.chrom, 
                   str(int(non_overlapping_terminal_exon.end)-100), 
                   str(non_overlapping_terminal_exon.end),
                   non_overlapping_terminal_exon.name,
                   ".",
                   non_overlapping_terminal_exon.strand+"\n"]))
    w.close()

    # remove duplicate entries
    os.system("sort -u %s | sort -k 1,1 -k2,2n > %s" % (os.path.join(options.out, "non_overlapping_terminal_exons_with_same_3p_site.bed"),
                                                        os.path.join(options.out, "non_overlapping_terminal_exons_with_same_3p_site_unique.bed")))

    # write mapping list for novel exons
    non_overlapping_terminal_exons_with_same_3p_site_unique = BedTool(os.path.join(options.out, "non_overlapping_terminal_exons_with_same_3p_site_unique.bed"))

    w = open(os.path.join(options.out, "non_overlapping_terminal_exons_with_same_3p_site_unique_mapping_list.tsv"), 'w')
    for te in non_overlapping_terminal_exons_with_same_3p_site_unique:
        w.write("\t".join([":".join([te.chrom,str(int(te.start)+1),str(te.end),te.name,te.strand,"annotated"]), te.name+"\n"]))
    w.close()

    # write mapping list for novel exons
    novel_exons = BedTool(os.path.join(options.out, "novel_exons.bed"))

    w = open(os.path.join(options.out, "novel_exons_mapping_list.tsv"), 'w')
    for ne in novel_exons:
        w.write("\t".join([":".join([ne.chrom,str(int(ne.start)+1),str(ne.end),ne.name,ne.strand,"novel"]), ne.name+"\n"]))
    w.close()

    # merge the two lists
    os.system("cat %s %s > %s" % (os.path.join(options.out, "non_overlapping_terminal_exons_with_same_3p_site_unique_mapping_list.tsv"),
                                  os.path.join(options.out, "novel_exons_mapping_list.tsv"),
                                  os.path.join(options.out, "exons2genes_list_final.tsv")))



# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)
