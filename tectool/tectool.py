# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

try:
    import sys
except(Exception):
    raise("[ERROR] sys was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import HTSeq
except(Exception):
    raise("[ERROR] HTSeq was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import os
except(Exception):
    raise("[ERROR] os was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import math
except(Exception):
    raise("[ERROR] math was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from argparse import ArgumentParser, RawTextHelpFormatter
except(Exception):
    raise("[ERROR] argparse was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import itertools
except(Exception):
    raise("[ERROR] itertools was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import collections
    from collections import defaultdict
except(Exception):
    raise("[ERROR] collections was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import pybedtools
except(Exception):
    raise("[ERROR] pybedtools was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from pyfasta import Fasta
except(Exception):
    raise("[ERROR] pyfasta was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import pandas as pd
except(Exception):
    raise("[ERROR] pandas was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
except(Exception):
    raise "[ERROR] plt from matplotlib.pyplot was not imported properly. Exiting."
    sys.exit(-1)

try:
    from sklearn.metrics import accuracy_score
    from sklearn.metrics import confusion_matrix
    from sklearn.metrics import classification_report
    from sklearn import metrics
    from sklearn.metrics import roc_curve, auc
    from sklearn.preprocessing import label_binarize
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn import linear_model
    from sklearn.cross_validation import train_test_split
    from sklearn import neighbors
    from sklearn.model_selection import StratifiedKFold
except(Exception):
    raise("[ERROR] sklearn was not imported properly")
    sys.exit(-1)

try:
    from scipy import interp
    from scipy import stats
except(Exception):
    raise("[ERROR] scipy was not imported properly")
    sys.exit(-1)

try:
    # for feature selection
    # http://rasbt.github.io/mlxtend/user_guide/feature_selection/SequentialFeatureSelector/
    from mlxtend.feature_selection import SequentialFeatureSelector as SFS
except(Exception):
    raise("[ERROR] mlxtend was not imported properly")
    sys.exit(-1)

try:
    import random
except(Exception):
    raise("[ERROR] random was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import numpy as np
except(Exception):
    raise("[ERROR] numpy was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import copy
except(Exception):
    raise("[ERROR] copy was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import csv
except(Exception):
    raise("[ERROR] csv was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from interval import Interval, IntervalSet
except(Exception):
    raise("[ERROR] interval was not imported properly. Exiting.")
    sys.exit(-1)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import our own modules
# -----------------------------------------------------------------------------

from gene_structure import Exon
from gene_structure import Transcript
from gene_structure import Gene
from aln_analysis import DetailedAlignment
from aln_analysis import SplitEvent
from aln_analysis import AnalysisUnit
from aln_analysis import FeatureCounts
from annotation import Annotation
from machine_learning import MachineLearningUnit
from machine_learning import BayesClassifier

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------


def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--annotation",
        dest="annotation",
        help="Annotation file GTF/GFF [REQUIRED]",
        metavar="FILE"
    )

    parser.add_argument(
        "--polyasites",
        dest="polyasites",
        help="Bed file that contains polya sites [REQUIRED]",
        metavar="FILE"
    )

    parser.add_argument(
        "--bam",
        dest="bam_file",
        help="The BAM file that should be analysed. [REQUIRED] " +
        "Note that the BAM file should be shorted by coordinates." +
        "And index file should be also present in the same directory.",
        metavar="FILE"
    )

    parser.add_argument(
        "--sequencing_direction",
        dest="sequencing_direction",
        help="Are the reads annotated with 'First' in the " +
        "BAM file mapped forward (on the right strand) " +
        "or reverse (on the opposite strand)? " +
        "Please choose 'forward' or 'reverse', " +
        "respectively. If the data are unstranded " +
        "please select 'unstranded'. [REQUIRED]"
    )

    parser.add_argument(
        "--genome",
        dest="genome",
        help="Genome in fasta format [REQUIRED]",
        metavar="FILE"
    )

    parser.add_argument(
        "--minimum_spliced_reads_for_cryptic_exon_start_site",
        dest="minimum_spliced_reads_for_cryptic_exon_start_site",
        default=5,
        help="Minimum number of spliced reads required to characterize " +
        "the start site of a cryptic exon. [default=5]"
    )

    # remove in the future
    parser.add_argument(
        "--min_readthough_count",
        dest="min_readthough_count",
        default=10,
        help="Minimum number of readthough reads required to extend " +
        "an existing exon. [default=10]"
    )

    # maybe we can also remove
    parser.add_argument(
        "--min_region_overlap",
        dest="min_region_overlap",
        default=1,
        help="min_region_overlap [default=1]"
    )

    # remove in the future
    parser.add_argument(
        "--max_splice_fuzziness",
        dest="max_splice_fuzziness",
        default=3,
        help="Maximum splice fuzziness [default=3]"
    )

    # remove in the future
    parser.add_argument(
        "--feature_length_threshold",
        dest="feature_length_threshold",
        default=10,
        help="Some poly(A) sites are located a few bases away from " +
        "annotated exons. It is difficult for TECtool to predict " +
        "if the reads correspond to splicing of an annotated exon " +
        "or if there is actual terminal exons. For this reason TECtool " +
        "reports such cases as computationally predicted. The flag " +
        "feature_length_threshold specifies the distance of the poly(A) " +
        "site to the annotated exons to allow such predictions. " +
        "[default=10]"
    )

    parser.add_argument(
        "--drop_manually_selected_features",
        dest="drop_manually_selected_features",
        action="store_true",
        default=False,
        help="Flag to not use the manually selected " +
        "features (['ReadsOUTvsIN_all', " +
        "'entropy_efficiency']). When this flag is used " +
        "all the features will be selected by greedy."
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        dest="output_dir",
        help="The path to the output directory. [REQUIRED]"
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="Verbose"
    )

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

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # some values we have to hand over
    # -------------------------------------------------------------------------
    minimum_spliced_reads_for_cryptic_exon_start_site = \
        int(options.minimum_spliced_reads_for_cryptic_exon_start_site)
    min_readthough_count = int(options.min_readthough_count)
    min_region_overlap = int(options.min_region_overlap)
    max_splice_fuzziness = int(options.max_splice_fuzziness)
    feature_length_threshold = int(options.feature_length_threshold)

    if options.drop_manually_selected_features:
        manually_selected_features = list()
    else:
        manually_selected_features = ['ReadsOUTvsIN_all', 'entropy_efficiency']

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Create output directories
    # -------------------------------------------------------------------------

    # Define directory for annotation files that are generated
    annotation_dir = os.path.join(options.output_dir, 'annotation')

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    if not os.path.exists(annotation_dir):
        os.makedirs(annotation_dir)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # create the annotation object and read in the annotation file
    # -------------------------------------------------------------------------
    sys.stdout.write("Reading annotation file:\t%s\n" % options.annotation)
    annotation = Annotation(annotation_id=options.annotation,
                            tmp=options.output_dir)
    annotation.parse(options.annotation)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # read in polya site information
    # -------------------------------------------------------------------------
    sys.stdout.write("Reading poly(A) sites file:\t%s\n" % options.polyasites)

    # read bed-regions (pybedtools.bedtool.BedTool object)
    polya_regions = pybedtools.BedTool(options.polyasites)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine intronic polya sites
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining poly(A) sites located within introns " +
                     "that do not overlap with any other genes (located on " +
                     "the same strand)...\n")

    exons_per_gene_bed_file = os.path.join(
        annotation_dir,
        '1.exons_per_gene.bed')
    transcripts_per_gene_bed_file = os.path.join(
        annotation_dir,
        '2.transcripts_per_gene.bed')
    genes_bed_file = os.path.join(
        annotation_dir,
        '3.genes.bed')
    polya_sites_in_transcripts_bed_file = os.path.join(
        annotation_dir,
        '4.polya_sites_in_transcripts.bed')
    polya_sites_in_introns_bed_file = os.path.join(
        annotation_dir,
        '5.polya_sites_in_introns.bed')
    union_exons_bed_file = os.path.join(
        annotation_dir,
        '6.union_exons.bed')
    union_introns_bed_file = os.path.join(
        annotation_dir,
        '7.union_introns.bed')
    introns_that_do_not_overlap_with_exons_stranded_bed_file = os.path.join(
        annotation_dir,
        '8.introns_that_do_not_overlap_with_exons_stranded.bed')
    introns_that_do_not_overlap_with_exons_unstranded_bed_file = os.path.join(
        annotation_dir,
        '9.introns_that_do_not_overlap_with_exons_unstranded.bed')
    polyasites_in_introns_stranded_bed_file = os.path.join(
        annotation_dir,
        '10.polyasites_in_introns_stranded.bed')
    polyasites_in_introns_unstranded_bed_file = os.path.join(
        annotation_dir,
        '11.polyasites_in_introns_unstranded.bed')

    annotation.determine_intronic_sites(
        polya_regions=polya_regions,
        sequencing_direction=options.sequencing_direction,
        exons_per_gene_bed_file=exons_per_gene_bed_file,
        transcripts_per_gene_bed_file=transcripts_per_gene_bed_file,
        genes_bed_file=genes_bed_file,
        polya_sites_in_transcripts_bed_file=polya_sites_in_transcripts_bed_file,
        polya_sites_in_introns_bed_file=polya_sites_in_introns_bed_file,
        union_exons_bed_file=union_exons_bed_file,
        union_introns_bed_file=union_introns_bed_file,
        introns_that_do_not_overlap_with_exons_stranded_bed_file=introns_that_do_not_overlap_with_exons_stranded_bed_file,
        introns_that_do_not_overlap_with_exons_unstranded_bed_file=introns_that_do_not_overlap_with_exons_unstranded_bed_file,
        polyasites_in_introns_stranded_bed_file=polyasites_in_introns_stranded_bed_file,
        polyasites_in_introns_unstranded_bed_file=polyasites_in_introns_unstranded_bed_file
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine intergenic poly(A) sites
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining intergenic poly(A) sites...\n")

    intergenic_polyasites_stranded_bed_file = os.path.join(
        annotation_dir, '12.intergenic_polyasites_stranded.bed')
    intergenic_polyasites_unstranded_bed_file = os.path.join(
        annotation_dir, '13.intergenic_polyasites_unstranded.bed')

    annotation.determine_intergenic_sites(
        polya_regions=polya_regions,
        genes_bed_file=genes_bed_file,
        intergenic_polyasites_stranded_bed_file=\
            intergenic_polyasites_stranded_bed_file,
        intergenic_polyasites_unstranded_bed_file=intergenic_polyasites_unstranded_bed_file
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine dataset to train the model.
    # Here we find all the last exons of a transcript and all the intermediate
    # exons of the transcripts. They are written in bed format
    # -------------------------------------------------------------------------

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # TERMINAL EXONS
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining annotated terminal that do not" +
                     "overlap with other annotated terminal exons...\n")

    # create file path
    terminal_exons_file_core_name = "annotated_terminal_exons"
    terminal_exons_bed_file_path = \
        os.path.join(options.output_dir,
                     (terminal_exons_file_core_name + ".bed"))

    # get terminal exons
    annotation.determine_terminal_exons(
        terminal_exons_bed_file_path=terminal_exons_bed_file_path)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # INTERMEDIATE EXONS
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining annotated intermediate that do not" +
                     "overlap with other annotated intermediate exons...\n")

    # create file path
    intermediate_exons_file_core_name = "annotated_intermediate_exons"
    intermediate_exons_bed_file_path = \
        os.path.join(options.output_dir,
                     (intermediate_exons_file_core_name + ".bed"))

    # get intermediate exons
    annotation.determine_intermediate_exons(
        intermediate_exons_bed_file_path=intermediate_exons_bed_file_path)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # TERMINAL EXONS (CLEAN) that do not overlap with intermediate exons
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining terminal exons that do not" +
                     "overlap with intermediate exons...\n")

    terminal_exons_clean_file_core_name = "annotated_terminal_exons_clean"
    terminal_exons_clean_bed_file_path = \
        os.path.join(options.output_dir,
                     (terminal_exons_clean_file_core_name + ".bed"))

    annotation.determine_clean_terminal_exons(
        terminal_exons_clean_bed_file_path=terminal_exons_clean_bed_file_path,
        terminal_exons_bed_file_path=terminal_exons_bed_file_path,
        intermediate_exons_bed_file_path=intermediate_exons_bed_file_path
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # INTERMEDIATE EXONS (CLEAN) that do not overlap with terminal exons
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining intermediate exons that do not" +
                     "overlap with terminal exons...\n")

    intermediate_exons_clean_file_core_name = \
        "annotated_intermediate_exons_clean"
    intermediate_exons_clean_bed_file_path = \
        os.path.join(options.output_dir,
                     (intermediate_exons_clean_file_core_name + ".bed"))

    annotation.determine_clean_intermediate_exons(
        intermediate_exons_clean_bed_file_path=intermediate_exons_clean_bed_file_path,
        terminal_exons_bed_file_path=terminal_exons_bed_file_path,
        intermediate_exons_bed_file_path=intermediate_exons_bed_file_path
    )

    # _________________________________________________________________________
    # ----------------------------------------------------------------------
    # BACKGROUND REGIONS (Create extended exons containing intergenic regions)
    # ----------------------------------------------------------------------

    # sys.stdout.write("Determining background regions... \n")

    # final_terminal_exons_bed_file = os.path.join(annotation_dir, '14.final_terminal_exons.bed')
    # background_regions_file_core_name = "background_regions"
    # background_regions_bed_file_path_stranded = os.path.join(options.output_dir,background_regions_file_core_name+"_stranded.bed")
    # background_regions_bed_file_path_unstranded = os.path.join(options.output_dir,background_regions_file_core_name+"_unstranded.bed")

    # annotation.create_background_with_intergenic_polya_sites(
    #     genes_bed_file=genes_bed_file,
    #     final_terminal_exons_bed_file=final_terminal_exons_bed_file,
    #     intergenic_polyasites_stranded_bed_file=intergenic_polyasites_stranded_bed_file,
    #     intergenic_polyasites_unstranded_bed_file=intergenic_polyasites_unstranded_bed_file,
    #     background_regions_stranded_bed_file=background_regions_bed_file_path_stranded,
    #     background_regions_unstranded_bed_file=background_regions_bed_file_path_unstranded,
    #     bases_to_extend=5000
    # )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # INTRONIC REGIONS without poly(A) sites
    # -------------------------------------------------------------------------

    intronic_regions_file_core_name = "intronic_regions"
    intronic_regions_no_polya_bed_file_path = os.path.join(
        options.output_dir,
        (intronic_regions_file_core_name + "_without_polya.bed"))

    sys.stdout.write(
        "Counting itnronic regions that do not contain polya sites... \n")

    annotation.determine_introns_without_polyAsite(
        polya_regions=polya_regions,
        intronic_regions_no_polya_bed_file_path=intronic_regions_no_polya_bed_file_path
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Create bed file with the non overlapping genes and update information
    # in the gene objects.
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining non overlapping genes...\n")

    non_overlapping_genes_bed = os.path.join(
        annotation_dir, "14.non_overlapping_genes.bed")

    annotation.determine_non_overlapping_genes(
        all_gene_regions=genes_bed_file,
        sequencing_direction=options.sequencing_direction,
        non_overlapping_genes_bed=non_overlapping_genes_bed
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine union exons length. The length for each of the genes is stored
    # in the annotations.gene[gene_id].union_exon_length variable
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining union exons lengths...\n")

    annotation.determine_union_exon_length_per_gene(
        union_exons_bed=union_exons_bed_file)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine feature regions
    # -------------------------------------------------------------------------
    sys.stdout.write(
        "Determining feature regions (=region between an " +
        "intronic poly(A) site and the 3' end of the closest " +
        "upstream exon)...\n")
    annotation.determine_feature_regions()

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # CREATE an AnalysisUnit for each FEATURE REGION
    # -------------------------------------------------------------------------
    aunits_dict = dict()

    for idx, feat_region in enumerate(annotation.feature_regions):
        if feat_region not in aunits_dict:

            aunits_dict[feat_region] = AnalysisUnit(
                unit_id=feat_region,
                potential_5pSS_exons=annotation.feature_regions_upstream_coordinates[feat_region],
                gene_id=annotation.feature_regions_genes[idx]
            )
        else:

            sys.stderr.write("ERROR: Region " +
                             feat_region +
                             " has multiple entries in the " +
                             " feature regions. \n")
            sys.exit(-1)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Open the BAM file
    # -------------------------------------------------------------------------
    # now count the things
    bam_file_path = options.bam_file
    bam = HTSeq.BAM_Reader(bam_file_path)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Go over all AnalysisUnit objects, fetch the reads, count, and enrich the
    # annotation object
    # -------------------------------------------------------------------------

    sys.stdout.write("Counting feature regions...\n")

    unit_nr = 0
    # for region in feature regions:
    for unit_id in aunits_dict.keys():
        unit_nr += 1

        # give some feedback about the state of the script
        # (how many units have been analyzed so far?)
        if (unit_nr % 100) == 0:
            sys.stderr.write("Regions processed:\t" + str(unit_nr) + "\n")

        # get the AnalysisUnit object
        aunits_dict[unit_id].analyze_reads(
            bam=bam,
            sequencing_direction=options.sequencing_direction,
            min_region_overlap=min_region_overlap,
            splice_fuzziness=max_splice_fuzziness,
            minimum_spliced_reads_for_cryptic_exon_start_site=minimum_spliced_reads_for_cryptic_exon_start_site,
            min_readthough_count=min_readthough_count,
            feature_length_threshold=feature_length_threshold,
            count_unique_mapping_reads_only=True,
            tmp=options.output_dir,
            annotation=annotation,
            exons_per_gene_bed_file=exons_per_gene_bed_file,
            verbose=False)

        # free memory
        try:
            del(aunits_dict[unit_id])
        except(KeyError):
            pass

    # dict_of_validated_novel_exons = "write_all"

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # MACHINE LEARNING
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # CREATE AN ML UNIT
    # -------------------------------------------------------------------------
    ml_train_unit = MachineLearningUnit()

    # -------------------------------------------------------------------------
    # Find genes with potential novel terminal exons
    # -------------------------------------------------------------------------
    expressed_genes = dict()

    for gene_id in annotation.genes:
        if (annotation.genes[gene_id].has_potential_novel_terminal_exon()) and (annotation.genes[gene_id].overlaps_with_other_gene is False):

            expressed_genes[gene_id] = gene_id

    # -------------------------------------------------------------------------
    # Find non overlapping genes
    # -------------------------------------------------------------------------
    non_overlapping_genes = dict()

    for gene_id in annotation.genes:
        if annotation.genes[gene_id].overlaps_with_other_gene is False:
            non_overlapping_genes[gene_id] = gene_id

    # -------------------------------------------------------------------------
    # Merge non overlapping and expressed genes
    # -------------------------------------------------------------------------
    merged_genes = expressed_genes.copy()
    merged_genes.update(non_overlapping_genes)

    # -------------------------------------------------------------------------
    # TERMINAL EXONS - TRAINING SET CREATION
    # -------------------------------------------------------------------------
    # create the terminal exon training set
    terminal_exons_statistics_file_path = \
        os.path.join(options.output_dir,
                     (terminal_exons_clean_file_core_name + ".tsv"))
    # create the terminal exons training data set
    ml_train_unit.create_terminal_exon_training_set(
        terminal_exons_bed_file_path=terminal_exons_clean_bed_file_path,
        sequencing_direction=options.sequencing_direction,
        max_splice_fuzziness=max_splice_fuzziness,
        output_dir=options.output_dir,
        genes_to_consider_dict=merged_genes,
        bam_file_path=options.bam_file,
        annotation=annotation,
        threshold_to_filter=minimum_spliced_reads_for_cryptic_exon_start_site
    )

    # ---------------------------------------------------------------------
    # TERMINAL EXONS - FILTER TRAINING SET
    # Remove terminal exons that the 3p crossing is not 0
    # and the sum of the last 100 bases of the profile is not > 0 reads
    # ---------------------------------------------------------------------
    for gene_id in annotation.genes:
        # copy list
        list_of_terminal_exons = list(
            annotation.genes[gene_id].annotated_terminal_exons)
        # loop over the annotated terminal exons
        for terminal_exon in list_of_terminal_exons:
            # get profiles
            profile = str(terminal_exon.profile).split(",")
            # keep the last 100 bases
            if len(profile) > 100:
                profile = profile[-100:]
            # count the per base reads
            sum_profile_end = sum([int(i) for i in profile])
            # apply filters and remove from original list
            if (int(terminal_exon.unspliced_3pSS) != 0) or (sum_profile_end == 0):
                annotation.genes[gene_id].annotated_terminal_exons.remove(
                    terminal_exon)
                # annotation.genes[gene_id].background.append(terminal_exon)

    # ---------------------------------------------------------------------
    # Find genes with expressed terminal exons that do not overlap with other
    #
    # Enrich the expressed genes dictionary
    # ---------------------------------------------------------------------
    for gene_id in annotation.genes:
        if annotation.genes[gene_id].has_annotated_terminal_exon() and (annotation.genes[gene_id].overlaps_with_other_gene is False):
            if gene_id not in expressed_genes:
                expressed_genes[gene_id] = gene_id

    # ---------------------------------------------------------------------
    # INTERMEDIATE EXONS - TRAINING SET CREATION
    # ---------------------------------------------------------------------
    # create the intermediate exon training set
    intermediate_exons_statistics_file_path = \
        os.path.join(options.output_dir,
                     (intermediate_exons_clean_file_core_name + ".tsv"))
    ml_train_unit.create_intermediate_exon_training_set(
        intermediate_exons_bed_file_path=intermediate_exons_clean_bed_file_path,
        sequencing_direction=options.sequencing_direction,
        max_splice_fuzziness=max_splice_fuzziness,
        output_dir=options.output_dir,
        genes_to_consider_dict=expressed_genes,
        bam_file_path=options.bam_file,
        annotation=annotation,
        threshold_to_filter=minimum_spliced_reads_for_cryptic_exon_start_site
    )

    # # ---------------------------------------------------------------------
    # # BACKGROUND EXONS - TRAINING SET CREATION
    # # ---------------------------------------------------------------------
    # # create the background training set
    # background_regions_statistics_file_path = \
    #     os.path.join(options.output_dir,
    #                 (background_regions_file_core_name + ".tsv"))

    background_regions_file_core_name = "background_regions"
    background_regions_statistics_file_path = os.path.join(
        options.output_dir, (background_regions_file_core_name + ".tsv"))

    # background_bed_file_path = background_regions_bed_file_path_stranded

    # if options.sequencing_direction == "unstranded":
    #     background_bed_file_path = \
    #       background_regions_bed_file_path_unstranded
    # ml_train_unit.create_background_training_set(
    #     background_bed_file_path = background_bed_file_path,
    #     sequencing_direction = options.sequencing_direction,
    #     max_splice_fuzziness = max_splice_fuzziness,
    #     output_dir = options.output_dir,
    #     genes_to_consider_dict = expressed_genes,
    #     bam_file_path = options.bam_file,
    #     annotation = annotation,
    #     threshold_to_filter = \
    #       minimum_spliced_reads_for_cryptic_exon_start_site
    # )

    # # ----------------------------------------------------------------------
    # # INTRONIC REGIONS without polya site
    # # ----------------------------------------------------------------------
    # # estimate intronic expression
    # intronic_regions_without_polya_statistics_file_path = os.path.join(
    #     options.output_dir,
    #     (intronic_regions_file_core_name + "_without_polya.tsv")
    # )
    # ml_train_unit.estimate_intronic_expression(
    #     intron_regions_bed_file_path = intronic_regions_no_polya_bed_file_path,
    #     sequencing_direction = options.sequencing_direction,
    #     max_splice_fuzziness = max_splice_fuzziness,
    #     output_dir = options.output_dir,
    #     genes_to_consider_dict = annotation.genes,
    #     bam_file_path = options.bam_file,
    #     annotation = annotation,
    #     threshold_to_filter = 0
    # )

    # ---------------------------------------------------------------------
    # Estimate gene expression
    # ---------------------------------------------------------------------
    ml_train_unit.estimate_gene_expression(
        annotation=annotation,
        genes_to_consider_dict=annotation.genes,
        bam_file_path=options.bam_file,
        sequencing_direction=options.sequencing_direction,
        count_unique_mapping_reads_only=True,
        verbose=False
    )

    # ---------------------------------------------------------------------
    # Create a dataframe and supplement it with gene expression information
    # ---------------------------------------------------------------------
    novel_terminal_output_file = os.path.join(
        options.output_dir, "novel_terminal_exons.tsv")

    ml_train_unit.create_terminal_exon_candidates_dataframe(
        annotation=annotation,
        novel_terminal_output_file=novel_terminal_output_file+"_pre",
        verbose=False
    )

    ml_train_unit.filter_terminal_exon_candidates_that_overlap_with_annotated_exons(
        annotation=annotation,
        novel_terminal_output_file=novel_terminal_output_file,
        sequencing_direction=options.sequencing_direction,
        exons_per_gene_bed_file=exons_per_gene_bed_file,
        verbose=False
    )

    # # ---------------------------------------------------------------------
    # # Create dataframes and supplement it with gene expression information
    # # ---------------------------------------------------------------------
    # novel_readthough_output_file_core_name = "novel_readthough_exons"
    # novel_readthough_output_file = os.path.join(
    #     options.output_dir,
    #     novel_readthough_output_file_core_name+".tsv"
    # )
    # ml_train_unit.create_terminal_exon_readthough_candidates_dataframe(
    #   annotation = annotation,
    #   novel_readthough_output_file = novel_readthough_output_file,
    #   verbose = False)


    # Generate N=10 training and validation datasets
    for run_number in range(10):

        sys.stdout.write("Generating #%d training and validation datasets ...\n" % (run_number))

        # ---------------------------------------------------------------------
        # Create dataframes and supplement with gene expression information
        # ---------------------------------------------------------------------
        ml_train_unit.create_training_dataframes(
            annotation=annotation,
            intermediate_output_file=intermediate_exons_statistics_file_path,
            terminal_output_file=terminal_exons_statistics_file_path,
            background_output_file=background_regions_statistics_file_path
        )
    
        # ---------------------------------------------------------------------
        # Calculate features for training data
        # ---------------------------------------------------------------------
        ml_train_unit.add_features_to_training_dataframes(
            nr_to_subsample=1100,  # "all",
            output_files_dir=options.output_dir,
            verbose=True
        )
    
        # ---------------------------------------------------------------------
        # Sample training data
        # ---------------------------------------------------------------------
        ml_train_unit.load_training_data(
            training_data_set_size=1000,  # "max_equal_size",
            validation_data_fraction=0.2,
            output_files_dir=options.output_dir,
            run_number = run_number,
            verbose=True
        )

        # Append to lists of training and validation datasets
        ml_train_unit.list_of_training_dataframes.append(ml_train_unit.training_df)
        ml_train_unit.list_of_validation_dataframes.append(ml_train_unit.validation_df)

    # ---------------------------------------------------------------------
    # Calculate features
    # ---------------------------------------------------------------------
    ml_train_unit.add_features_to_terminal_exon_candidates_dataframe(
        output_files_dir=options.output_dir,
        verbose=True
    )

    # ---------------------------------------------------------------------
    # Train a classifier
    # ---------------------------------------------------------------------
    sys.stdout.write("Training and validating classifiers...\n")
    # ["Bayes", "KNeighbors", "multiclass_SVC"]
    for classifier in ["Bayes"]:

        # ml_train_unit.training_df =  ml_train_unit.list_of_training_dataframes[0]
        # ml_train_unit.validation_df = ml_train_unit.list_of_validation_dataframes[0]

        # -----------------------------------------------------------------
        # classifier = KNeighborsClassifier(id=classifier_name)
        # ml_train_unit.register_classifier(
        #     classifier=classifier,
        #     name=classifier_name
        # )
        # ml_train_unit.load_classifier(classifier_name)
        # ml_train_unit.load_training_data(ml_train_unit.training_data...)
        # ml_train_unit.greedy_feature_selection()
        # ml_train_unit.train_on_selected_features()
        # ml_train_unit.load_terminal_exon_candidates()
        # ml_train_unit.classify_terminal_exon_candidates()
        # -----------------------------------------------------------------

        # FIXME: that's a dirty fix!
        ml_train_unit.selected_features = None

        # -----------------------------------------------------------------
        # Create a results directory for the current classifier
        # -----------------------------------------------------------------
        classifier_results_dir = os.path.join(
            options.output_dir, ("RESULTS_" + classifier))

        # create the directory
        if not os.path.exists(classifier_results_dir):
            os.makedirs(classifier_results_dir)

        # give user feedback
        sys.stdout.write(
            " :: Used classifier: %s -> Results can be found in: %s\n" %
            (classifier, classifier_results_dir)
        )

        # ---------------------------------------------------------------------
        # Create dataframes and supplement with gene expression information
        # ---------------------------------------------------------------------
        ml_train_unit.create_training_dataframes(
            annotation=annotation,
            intermediate_output_file=intermediate_exons_statistics_file_path,
            terminal_output_file=terminal_exons_statistics_file_path,
            background_output_file=background_regions_statistics_file_path
        )
    
        # ---------------------------------------------------------------------
        # Calculate features for training data
        # ---------------------------------------------------------------------
        ml_train_unit.add_features_to_training_dataframes(
            nr_to_subsample="all",  # "all",1100
            output_files_dir=options.output_dir,
            verbose=True
        )
    
        # ---------------------------------------------------------------------
        # Sample training data
        # ---------------------------------------------------------------------
        ml_train_unit.load_training_data(
            training_data_set_size="max_equal_size",  # "max_equal_size",1000
            validation_data_fraction=0.2,
            output_files_dir=options.output_dir,
            run_number=-1,
            verbose=True
        )

        # -----------------------------------------------------------------
        # Select features using Greedy algorithm
        # -----------------------------------------------------------------
        # FIXME: hand over the default features
        ml_train_unit.greedy_feature_selection(
            number_of_randomization=10,
            classifier=classifier,
            manually_selected_features=manually_selected_features,
            verbose=True
        )

        # -----------------------------------------------------------------
        # Calculate features
        # -----------------------------------------------------------------
        # FIXME: remove the selected_features dependency from this function
        ml_train_unit.load_terminal_exon_candidates(verbose=True)
        
        for run_number in range(len(ml_train_unit.list_of_training_dataframes)):

            ml_train_unit.training_df =  ml_train_unit.list_of_training_dataframes[run_number]
            ml_train_unit.validation_df = ml_train_unit.list_of_validation_dataframes[run_number]

            # -----------------------------------------------------------------
            # Train a classifier on the greedy selected features
            # -----------------------------------------------------------------
            ml_train_unit.train_classifier_on_selected_features(
                results_dir_path=options.output_dir,
                classifier=classifier,
                nr_of_train_vs_test_runs=25,
                verbose=True
            )

            # -----------------------------------------------------------------
            # Calculate probabilities
            # -----------------------------------------------------------------
            ml_train_unit.fill_probabilities_dictionary(
                classifier=classifier,
                results_dir=options.output_dir,
                verbose=True
            )

        # -----------------------------------------------------------------
        # Classify novel terminal exon candidates as...
        #  - terminal exons
        #  - intermediate exons
        #  - background regions
        # -----------------------------------------------------------------
        ml_train_unit.classify_terminal_exon_candidates(
            results_dir=options.output_dir,
            classifier=classifier,
            verbose=True
        )

    # ---------------------------------------------------------------------
    # Create dictionary of classified novel exons (0 based)
    # ---------------------------------------------------------------------
    dict_of_validated_novel_exons = dict()

    for exon_id in ml_train_unit.selected_novel_terminal_exons.index.tolist():

        exon_id_sp = exon_id[0].split(":")
        exon_id_0_based = ":".join(
            [exon_id_sp[0],
             str(int(exon_id_sp[1]) - 1),
             exon_id_sp[2],
             exon_id_sp[3]]
        )

        dict_of_validated_novel_exons[exon_id_0_based] = exon_id_0_based

    # ---------------------------------------------------------------------
    # FIXME: implement me! TRAIN ML UNIT (using different approaches)
    # ---------------------------------------------------------------------
    # SVM
    # ml_train_unit.determine_weights_by_SVM(weights_key_name="SVM")
    # ml_train_unit.validate_weights_performance(weights_key_name="SVM")

    # Logistic regression
    # ml_train_unit.determine_weights_by_LogisticRegression(weights_key_name="LogisticRegression")
    # ml_train_unit.validate_weights_performance(weights_key_name="LogisticRegression")

    # compare the performance of weights determined by two different
    # machine learning approaches
    # ml_train_unit.compare_weights_performance(weights_key_name_A="SVM",
    # weights_key_name_B="LogisticRegression")

    # ml_train_unit.classify(weights_key_name="SVM"s

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Read genome file
    # -------------------------------------------------------------------------
    sys.stdout.write("Reading genome file ...\n")
    genome = Fasta(options.genome)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Annotate CDS for novel transcripts
    # -------------------------------------------------------------------------
    sys.stdout.write(
        "Annotating CDS for novel transcripts (novel splicing)...\n")
    annotation.determine_CDS_novel_transcrtips_from_spicing(genome)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Write novel GTF file
    # -------------------------------------------------------------------------
    sys.stdout.write("Writing GTF file for known and novel transcripts ...\n")
    output_file = "enriched_annotation.gtf"
    output_file_path = os.path.join(options.output_dir, output_file)
    w = open(output_file_path, 'w')
    annotation.write_gtf(
        w, with_CDS=True, accepted_exons_dict=dict_of_validated_novel_exons)
    w.close()

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Write file with transcript id and novel transcript id
    # -------------------------------------------------------------------------
    output_file = "transcript_id_2_novel_transcript_id.tsv"
    output_file_path = os.path.join(options.output_dir, output_file)
    w = open(output_file_path, 'w')
    annotation.write_transcript2novel_transcript_mapping_list(
        output_file_handle=w,
        accepted_exons_dict=dict_of_validated_novel_exons
    )
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
