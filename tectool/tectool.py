#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

__version__ = "0.2.0"

try:
    import sys
except Exception:
    raise Exception("[ERROR] sys was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import HTSeq
except Exception:
    raise Exception("[ERROR] HTSeq was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import os
except Exception:
    raise Exception("[ERROR] os was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import math
except Exception:
    raise Exception("[ERROR] math was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from argparse import ArgumentParser, RawTextHelpFormatter
except Exception:
    raise Exception("[ERROR] argparse was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import itertools
except Exception:
    raise Exception("[ERROR] itertools was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import collections
    from collections import defaultdict
except Exception:
    raise Exception("[ERROR] collections was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import pybedtools
except Exception:
    raise Exception("[ERROR] pybedtools was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from pyfasta import Fasta
except Exception:
    raise Exception("[ERROR] pyfasta was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import pandas as pd
except Exception:
    raise Exception("[ERROR] pandas was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
except Exception:
    raise Exception("[ERROR] plt from matplotlib.pyplot was not imported properly. Exiting.")
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
    # from sklearn.model_selection import StratifiedKFold
except Exception:
    raise Exception("[ERROR] sklearn was not imported properly")
    sys.exit(-1)

try:
    from scipy import interp
    from scipy import stats
except Exception:
    raise Exception("[ERROR] scipy was not imported properly")
    sys.exit(-1)

# try:
#     # for feature selection
#     # http://rasbt.github.io/mlxtend/user_guide/feature_selection/SequentialFeatureSelector/
#     from mlxtend.feature_selection import SequentialFeatureSelector as SFS
# except(Exception):
#     raise("[ERROR] mlxtend was not imported properly")
#     sys.exit(-1)

try:
    import random
except Exception:
    raise Exception("[ERROR] random was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import numpy as np
except Exception:
    raise Exception("[ERROR] numpy was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import copy
except Exception:
    raise Exception("[ERROR] copy was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import csv
except Exception:
    raise Exception("[ERROR] csv was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import functools
except Exception:
    raise Exception("[ERROR] functools was not imported properly. Exiting.")

# try:
#     import multiprocessing as mp
# except Exception:
#     raise Exception("[ERROR] multiprocessing was not imported properly. Exiting.")
#     sys.exit(-1)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import our own modules
# -----------------------------------------------------------------------------

import gene_structure.exon
import gene_structure.transcript
import gene_structure.gene
import aln_analysis.detailed_alignment
import aln_analysis.split_event
import aln_analysis.analysis_unit
import aln_analysis.feature_counts
import annotation.annotation
import machine_learning.machine_learning_unit
import machine_learning.bayes_classifier

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Classes
# -----------------------------------------------------------------------------


class FileCounter(object):

    def __init__(self, counting_start=1):

        self.file_number = counting_start - 1

    def get_new_filenumber(self):

        self.file_number += 1

        return(self.file_number)

    def get_new_filenumber_str(self):

        self.file_number += 1

        return(str(self.file_number))

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Custom functions
# -----------------------------------------------------------------------------


def estimate_expression_of_selected_genes(
    annotation,
    selected_genes_dict,
    bam_file_path,
    sequencing_direction,
    count_unique_mapping_reads_only=True,
    verbose=False
):
    # _____________________________________________________________________
    # ---------------------------------------------------------------------
    # create an AnalysisUnit object for each gene region
    # and store it in a dictionary
    # ---------------------------------------------------------------------

    # dictionary for the genes
    aunits_genes_dict = dict()

    for gene_id in annotation.genes:

        # filter genes that do not have expressed terminal exons
        if gene_id in selected_genes_dict:

            chromosome = annotation.genes[gene_id].chromosome
            start = int(annotation.genes[gene_id].start)
            end = int(annotation.genes[gene_id].end)
            strand = annotation.genes[gene_id].strand

            gene_iv = HTSeq.GenomicInterval(chromosome, start, end, strand)

            if gene_iv not in aunits_genes_dict:

                aunits_genes_dict[gene_iv] = AnalysisUnit(
                    unit_id=gene_iv,
                    potential_5pSS_exons=None,
                    gene_id=gene_id
                )

    # _____________________________________________________________________
    # ---------------------------------------------------------------------
    # Open the BAM file
    # ---------------------------------------------------------------------
    bam_file_path = bam_file_path
    bam = HTSeq.BAM_Reader(bam_file_path)

    # _____________________________________________________________________
    # ---------------------------------------------------------------------
    # Go over all AnalysisUnit objects for non overlapping genes, fetch
    # the reads and count
    # ---------------------------------------------------------------------

    sys.stdout.write("Estimating gene expression ..." + os.linesep)

    # go over each unit
    unit_nr = 0
    for unit_id in aunits_genes_dict.keys():

        unit_nr += 1

        # give some feedback about the state of the script
        # (how many units have been analyzed so far?)
        if (unit_nr % 100) == 0:
            sys.stderr.write("Regions processed:\t" +
                             str(unit_nr) + os.linesep)

        gene_reads, gene_id = aunits_genes_dict[unit_id
                                                ].estimate_gene_expression(
            bam=bam,
            unit_id=unit_id,
            sequencing_direction=sequencing_direction,
            count_unique_mapping_reads_only=True,
            annotation=annotation,
            verbose=False)

        annotation.genes[gene_id].total_reads = int(gene_reads)
        annotation.genes[gene_id].estimate_ExpressionPerKBApproximated()

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
        help="Annotation file gtf format [REQUIRED] (tested with ENSEMBL v87)",
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
        "Note that the BAM file should be shorted by coordinates. " +
        "An index file should be also present in the same directory.",
        metavar="FILE"
    )

    parser.add_argument(
        "--sequencing_direction",
        dest="sequencing_direction",
        help="Are the reads annotated with 'First' in the " +
        "BAM file mapped forward (on the right strand) " +
        "Please choose 'forward'. If the data are unstranded " +
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

    # maybe we can also remove
    parser.add_argument(
        "--min_region_overlap",
        dest="min_region_overlap",
        default=10,
        help="min_region_overlap [default=10]"
    )

    # remove in the future
    parser.add_argument(
        "--max_splice_fuzziness",
        dest="max_splice_fuzziness",
        default=0,
        help="Maximum splice fuzziness [default=3]"
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
        "--drop_intronic_polya_sites_of_overlapping_genes",
        dest="drop_intronic_polya_sites_of_overlapping_genes",
        action="store_true",
        default=False,
        help="Flag to ignore intronic polya sites shared " +
        "among overlapping genes"
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
    except Exception:
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # initialize a file counter
    fc = FileCounter()

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # some values we have to hand over
    # -------------------------------------------------------------------------
    minimum_spliced_reads_for_cryptic_exon_start_site = \
        int(options.minimum_spliced_reads_for_cryptic_exon_start_site)
    min_region_overlap = int(options.min_region_overlap)
    max_splice_fuzziness = int(options.max_splice_fuzziness)
    polyasites = options.polyasites

    # -------------------------------------------------------------------------
    sequencing_direction = options.sequencing_direction
    supported_sequencing_directions_forward = "forward"
    supported_sequencing_directions_unstranded = "unstranded"
    supported_sequencing_directions = \
        [supported_sequencing_directions_forward,
         supported_sequencing_directions_unstranded]

    if sequencing_direction not in supported_sequencing_directions:
        raise Exception("Provided sequencing direction parameter \
                         ({}) is not supported. Please choose: {} \
                         {}".format(sequencing_direction,
                                    " ".join(supported_sequencing_directions),
                                    os.linesep)
                        )

    # -------------------------------------------------------------------------
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
    sys.stdout.write("Reading annotation file:\t {} {}".format(
        options.annotation, os.linesep))
    annotation = Annotation(annotation_id=options.annotation,
                            tmp=options.output_dir)
    annotation.parse(options.annotation)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # Part to determine intronic poly(A) sites
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine exons of all genes and write them to a bed file
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting exons per gene in bed format" +
                     os.linesep)

    exons_per_gene_bed_file = os.path.join(
        annotation_dir,
        fc.get_new_filenumber_str() + '.exons_per_gene.bed')

    annotation.write_exons_bed(exons_per_gene_bed_file)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine transcripts of all genes and write them to a bed file
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting transcripts per gene in bed format" +
                     os.linesep)

    transcripts_per_gene_bed_file = os.path.join(
        annotation_dir,
        fc.get_new_filenumber_str() + '.transcripts_per_gene.bed')

    annotation.write_transcripts_bed(transcripts_per_gene_bed_file)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine genes and write them to a bed file
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting genes in bed format" +
                     os.linesep)

    genes_bed_file = os.path.join(
        annotation_dir,
        fc.get_new_filenumber_str() + '.genes.bed')

    annotation.write_genes_bed(genes_bed_file)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine the union exon per gene
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting union exons per gene in bed format" +
                     os.linesep)

    union_exons_bed_file = os.path.join(
        annotation_dir,
        fc.get_new_filenumber_str() + '.union_exons.bed')

    annotation.write_union_exon_bed(union_exons_bed_file)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine the union intron per gene
    # Important note: In case of overlapping genes it might be that an intron
    # in a gene might be part of an exon in another gene and vice versa.
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting union introns per gene in bed format" +
                     os.linesep)

    union_introns_bed_file = os.path.join(
        annotation_dir,
        fc.get_new_filenumber_str() + '.union_introns.bed')

    annotation.write_union_intron_bed(
        union_exons_bed=union_exons_bed_file,
        union_introns_bed=union_introns_bed_file
    )

    if sequencing_direction == supported_sequencing_directions_forward:

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Find introns that do not overlap with other exons (stranded)
        # ---------------------------------------------------------------------

        sys.stdout.write("Extracting union introns per gene that do not " +
                         "overlap with exons in bed format for stranded " +
                         "protocols" + os.linesep)

        introns_not_overlap_with_exons_bed_file = os.path.join(
            annotation_dir,
            fc.get_new_filenumber_str() +
            '.introns_that_do_not_overlap_with_exons_stranded.bed'
        )

        annotation.write_introns_that_do_not_overlap_with_exons_stranded_bed(
            union_exons_bed=union_exons_bed_file,
            union_introns_bed=union_introns_bed_file,
            selected_introns_bed=introns_not_overlap_with_exons_bed_file
        )

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Find intronic poly(A) sites for stranded protocol
        # -------------------------------------------------------------------------

        sys.stdout.write("Extracting intronic poly(A) sites for stranded " +
                         "protocols" + os.linesep)

        polyasites_in_introns_bed_file = os.path.join(
            annotation_dir,
            fc.get_new_filenumber_str() +
            '.polyasites_in_introns_stranded.bed')

        annotation.write_intronic_polya_sites(
            polyasites=polyasites,
            introns=introns_not_overlap_with_exons_bed_file,
            drop_intronic_polya_sites_of_overlapping_genes=options.drop_intronic_polya_sites_of_overlapping_genes,
            intronic_polya_sites_bed=polyasites_in_introns_bed_file
        )

    elif sequencing_direction == supported_sequencing_directions_unstranded:

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Find introns that do not overlap with other exons (unstranded)
        # ---------------------------------------------------------------------

        sys.stdout.write("Extracting union introns per gene that do not " +
                         "overlap with exons in bed format for unstranded " +
                         "protocols" + os.linesep)

        introns_not_overlap_with_exons_bed_file = os.path.join(
            annotation_dir,
            fc.get_new_filenumber_str() +
            '.introns_that_do_not_overlap_with_exons_unstranded.bed')

        annotation.write_introns_that_do_not_overlap_with_exons_unstranded_bed(
            union_exons_bed=union_exons_bed_file,
            union_introns_bed=union_introns_bed_file,
            selected_introns_bed=introns_not_overlap_with_exons_bed_file
        )

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Find intronic poly(A) sites for unstranded protocol
        # -------------------------------------------------------------------------

        sys.stdout.write("Extracting intronic poly(A) sites for unstranded " +
                         "protocols" + os.linesep)

        polyasites_in_introns_bed_file = os.path.join(
            annotation_dir,
            fc.get_new_filenumber_str() +
            '.polyasites_in_introns_unstranded.bed')

        annotation.write_intronic_polya_sites(
            polyasites=polyasites,
            introns=introns_not_overlap_with_exons_bed_file,
            drop_intronic_polya_sites_of_overlapping_genes=options.drop_intronic_polya_sites_of_overlapping_genes,
            intronic_polya_sites_bed=polyasites_in_introns_bed_file
        )

    else:
        raise Exception("Provided sequencing direction parameter ({}) is not supported. Please choose: {} {} ".format(
            sequencing_direction,
            " ".join(supported_sequencing_directions),
            os.linesep)
        )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # What regions to count for machine learning
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine all terminal exons from multi-exonic transcripts in bed format
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting all terminal exons from " +
                     "multi-exonic transcripts" + os.linesep)

    all_terminal_exons_bed_file = os.path.join(
        annotation_dir,
        fc.get_new_filenumber_str() +
        '.all_terminal_exons.bed')

    annotation.write_terminal_exons(bed=all_terminal_exons_bed_file)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine all intermediate exons from multi-exonic transcripts in bed
    # format
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting all intermediate exons from " +
                     "multi-exonic transcripts" + os.linesep)

    all_intermediate_exons_bed_file = os.path.join(
        annotation_dir,
        fc.get_new_filenumber_str() + '.all_intermediate_exons.bed')

    annotation.write_intermediate_exons(bed=all_intermediate_exons_bed_file)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine all start exons from multi-exonic transcripts in bed
    # format
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting all start exons from " +
                     "multi-exonic transcripts" + os.linesep)

    all_start_exons_bed_file = os.path.join(
        annotation_dir,
        fc.get_new_filenumber_str() + '.all_start_exons.bed')

    annotation.write_start_exons(bed=all_start_exons_bed_file)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine all terminal exons that do not overlap with other terminal
    # exons
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting terminal exons that do not overlap " +
                     "with other terminal exons" + os.linesep)

    terminal_exons_not_overlap_other_terminal_exons_bed_file = \
        os.path.join(
            annotation_dir,
            fc.get_new_filenumber_str() +
            '.terminal_exons_not_overlap_other_terminal_exons.bed'
        )

    annotation.write_non_overlapping_regions_to_bed(
        bed_in=all_terminal_exons_bed_file,
        bed_out=terminal_exons_not_overlap_other_terminal_exons_bed_file,
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine all intermediate exons that do not overlap with other
    # intermediate exons
    # -------------------------------------------------------------------------

    sys.stdout.write("Extracting intermediate exons that do not overlap " +
                     "with other intermediate exons" + os.linesep)

    intermediate_exons_not_overlap_other_intermediate_exons_bed_file = \
        os.path.join(
            annotation_dir,
            fc.get_new_filenumber_str() +
            '.intermediate_exons_not_overlap_other_intermediate_exons.bed'
        )

    annotation.write_non_overlapping_regions_to_bed(
        bed_in=all_intermediate_exons_bed_file,
        bed_out=intermediate_exons_not_overlap_other_intermediate_exons_bed_file,
    )

    if sequencing_direction == supported_sequencing_directions_forward:


        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Determine terminal exons than not ovelap with other exons
        # for stranded protocols
        # -------------------------------------------------------------------------

        sys.stdout.write("Extracting terminal exons that do not overlap " +
                         "with other exons (first exons, last exon or " +
                         "intermediate exons) for stranded protocols" +
                         os.linesep)

        terminal_exons_non_overlapping = \
            os.path.join(
                annotation_dir,
                fc.get_new_filenumber_str() +
                '.terminal_exons_non_overlapping_stranded.bed'
            )

        annotation.write_non_overlapping_regions(
            selected_bed=terminal_exons_not_overlap_other_terminal_exons_bed_file,
            compare_bed_1=all_intermediate_exons_bed_file,
            compare_bed_2=all_start_exons_bed_file,
            strand=True,
            bed_out=terminal_exons_non_overlapping
        )

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Determine intermediate exons than not ovelap with other exons
        # for stranded protocols
        # -------------------------------------------------------------------------

        sys.stdout.write("Extracting intermediate exons that do not overlap " +
                         "with other exons (first exons, last exon or " +
                         "intermediate exons) for stranded protocols" +
                         os.linesep)

        intermediate_exons_non_overlapping = \
            os.path.join(
                annotation_dir,
                fc.get_new_filenumber_str() +
                '.intermediate_exons_non_overlapping_stranded.bed'
            )

        annotation.write_non_overlapping_regions(
            selected_bed=intermediate_exons_not_overlap_other_intermediate_exons_bed_file,
            compare_bed_1=all_terminal_exons_bed_file,
            compare_bed_2=all_start_exons_bed_file,
            strand=True,
            bed_out=intermediate_exons_non_overlapping
        )

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Remove intermediate exons that overlap on the opposite strand.
        # This means that we keep only intermediate that are unique based on
        # the chromosome, start and end position
        # -------------------------------------------------------------------------

        sys.stdout.write("Remove intermediate exons that overlap with other " +
                         "intermediate exons on the opposite strand" +
                         os.linesep)

        intermediate_exons_non_overlapping_filtered = \
            os.path.join(
                annotation_dir,
                fc.get_new_filenumber_str() +
                '.intermediate_exons_non_overlapping_stranded_filtered.bed'
            )

        annotation.remove_regions_that_overlap_in_the_opposite_strand(
            bed_in=intermediate_exons_non_overlapping,
            bed_out=intermediate_exons_non_overlapping_filtered
        )

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Find non overlapping genes
        # -------------------------------------------------------------------------

        sys.stdout.write("Determining non overlapping genes for stranded " +
                         "protocols" + os.linesep)

        non_overlapping_genes_bed = os.path.join(
            os.path.join(
                annotation_dir,
                fc.get_new_filenumber_str() +
                '.non_overlapping_genes_stranded.bed'
            )
        )

        annotation.write_non_overlapping_genes_to_bed(
            all_gene_regions=genes_bed_file,
            strand=True,
            non_overlapping_genes_bed=non_overlapping_genes_bed
        )

    elif sequencing_direction == supported_sequencing_directions_unstranded:

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Determine terminal exons than not ovelap with other exons
        # for unstranded protocols
        # -------------------------------------------------------------------------

        sys.stdout.write("Extracting terminal exons that do not overlap " +
                         "with other exons (first exons, last exon or " +
                         "intermediate exons) for stranded protocols" +
                         os.linesep)

        terminal_exons_non_overlapping = \
            os.path.join(
                annotation_dir,
                fc.get_new_filenumber_str() +
                '.terminal_exons_non_overlapping_unstranded.bed'
            )

        annotation.write_non_overlapping_regions(
            selected_bed=terminal_exons_not_overlap_other_terminal_exons_bed_file,
            compare_bed_1=all_intermediate_exons_bed_file,
            compare_bed_2=all_start_exons_bed_file,
            strand=False,
            bed_out=terminal_exons_non_overlapping
        )

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Determine intermediate exons than not ovelap with other exons
        # for unstranded protocols
        # -------------------------------------------------------------------------

        sys.stdout.write("Extracting intermediate exons that do not overlap " +
                         "with other exons (first exons, last exon or " +
                         "intermediate exons) for unstranded protocols" +
                         os.linesep)

        intermediate_exons_non_overlapping = \
            os.path.join(
                annotation_dir,
                fc.get_new_filenumber_str() +
                '.intermediate_exons_non_overlapping_unstranded.bed'
            )

        annotation.write_non_overlapping_regions(
            selected_bed=intermediate_exons_not_overlap_other_intermediate_exons_bed_file,
            compare_bed_1=all_terminal_exons_bed_file,
            compare_bed_2=all_start_exons_bed_file,
            strand=False,
            bed_out=intermediate_exons_non_overlapping
        )

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Remove intermediate exons that overlap on the opposite strand.
        # This means that we keep only intermediate that are unique based on
        # the chromosome, start and end position
        # -------------------------------------------------------------------------

        sys.stdout.write("Remove intermediate exons that overlap with other " +
                         "intermediate exons on the opposite strand" +
                         os.linesep)

        intermediate_exons_non_overlapping_filtered = \
            os.path.join(
                annotation_dir,
                fc.get_new_filenumber_str() +
                '.intermediate_exons_non_overlapping_unstranded_filtered.bed'
            )

        annotation.remove_regions_that_overlap_in_the_opposite_strand(
            bed_in=intermediate_exons_non_overlapping,
            bed_out=intermediate_exons_non_overlapping_filtered
        )

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Find non overlapping genes
        # -------------------------------------------------------------------------

        sys.stdout.write("Determining non overlapping genes for unstranded " +
                         "protocols" + os.linesep)

        non_overlapping_genes_bed = os.path.join(
            os.path.join(
                annotation_dir,
                fc.get_new_filenumber_str() +
                '.non_overlapping_genes_unstranded.bed'
            )
        )

        annotation.write_non_overlapping_genes_to_bed(
            all_gene_regions=genes_bed_file,
            strand=False,
            non_overlapping_genes_bed=non_overlapping_genes_bed
        )

    else:
        raise Exception("Provided sequencing direction parameter ({}) is not supported. Please choose: {} {} ".format(
            sequencing_direction,
            " ".join(supported_sequencing_directions),
            os.linesep)
        )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Determine feature regions
    # The bed file provided should contain the gene name in the name field
    # If an intronic polya site is shared between two genes, then the bed file
    # should contain 2 entries. One for each of the genes that is shared with.
    # -------------------------------------------------------------------------
    sys.stdout.write("Determining feature regions (=region between an " +
                     "intronic poly(A) site and the 3' end of the closest " +
                     "upstream exon)..." + os.linesep)

    annotation.determine_feature_regions(
        polyasites=polyasites_in_introns_bed_file
    )

    # _______________________________________________________________________
    # -----------------------------------------------------------------------
    # Determine union exons length. The length for each of the genes is
    # stored in the annotations.gene[gene_id].union_exon_length variable
    # -----------------------------------------------------------------------
    sys.stdout.write("Determining union exons lengths..." + os.linesep)

    annotation.determine_union_exon_length_per_gene(
        union_exons_bed=union_exons_bed_file
    )

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
                             " feature regions. " + os.linesep)
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

    sys.stdout.write("Counting feature regions..." + os.linesep)

    unit_nr = 0
    # for region in feature regions:
    for unit_id, unit_value in list(aunits_dict.items()):
        unit_nr += 1

        # give some feedback about the state of the script
        # (how many units have been analyzed so far?)
        if (unit_nr % 100) == 0:
            sys.stdout.write("Regions processed:\t" +
                             str(unit_nr) +
                             os.linesep)

        # get the AnalysisUnit object
        aunits_dict[unit_id].analyze_reads(
            bam=bam,
            sequencing_direction=options.sequencing_direction,
            min_region_overlap=min_region_overlap,
            splice_fuzziness=max_splice_fuzziness,
            minimum_spliced_reads_for_cryptic_exon_start_site=minimum_spliced_reads_for_cryptic_exon_start_site,
            count_unique_mapping_reads_only=True,
            tmp=options.output_dir,
            annotation=annotation,
            exons_per_gene_bed_file=exons_per_gene_bed_file,
            verbose=False
        )

        # free memory
        try:
            del(aunits_dict[unit_id])
        except(KeyError):
            pass

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # MACHINE LEARNING
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # CREATE A MACHINE LEARNING UNIT
    # -------------------------------------------------------------------------
    ml_train_unit = MachineLearningUnit()

    # -------------------------------------------------------------------------
    # Create dictionary of multiexonic genes
    # (transcripts that consist of more than one exon).
    # -------------------------------------------------------------------------

    genes_with_multiexonic_transcripts_dict = \
        annotation.get_genes_with_multiexonic_transcripts(
            bed=union_exons_bed_file
        )

    if options.verbose:
        sys.stdout.write(
            "Number of genes with multiple exons : {} {}".format(
                str(len(genes_with_multiexonic_transcripts_dict)),
                os.linesep
            )
        )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # # Create dictionaty of non overlapping genes
    # -------------------------------------------------------------------------

    non_overlapping_genes_dict = \
        annotation.create_dictionary_of_non_overlapping_genes(
            non_overlapping_genes_bed=non_overlapping_genes_bed
        )

    if options.verbose:
        sys.stdout.write(
            "Number of non overlapping genes: {} {}".format(
                str(len(non_overlapping_genes_dict)),
                os.linesep
            )

        )

    genes_to_consider_for_ML = {x: non_overlapping_genes_dict[x] for x in non_overlapping_genes_dict if x in genes_with_multiexonic_transcripts_dict}

    if options.verbose:
        sys.stdout.write(
            "Number of genes to consider for training (intersection of non overlapping genes and genes with muliple exons): {} {}".format(
                str(len(genes_to_consider_for_ML)), os.linesep
            )
        )

    # -------------------------------------------------------------------------
    # TERMINAL EXONS - TRAINING SET CREATION
    # -------------------------------------------------------------------------
    # HINT: We only consider genes that have more than one exon, because
    # otherwise we might get genes that have only one exon. Such 'terminal'
    # exons would not have spliced reads that go in. Consequently, all the
    # counted features would confound our training set, since we do not
    # consider such candidates!

    ml_train_unit.create_terminal_exon_training_set(
        terminal_exons_bed_file_path=terminal_exons_non_overlapping,
        sequencing_direction=options.sequencing_direction,
        max_splice_fuzziness=max_splice_fuzziness,
        output_dir=options.output_dir,
        genes_to_consider_dict=genes_to_consider_for_ML,
        bam_file_path=options.bam_file,
        annotation=annotation,
        threshold_to_filter=minimum_spliced_reads_for_cryptic_exon_start_site,
        verbose=options.verbose
    )

    # ---------------------------------------------------------------------
    # FILTER TERMINAL EXONS TRAINING SET
    # Remove terminal exons for which the 3p crossing is not 0
    # and the sum of the last 100 bases of the profile is not > 0 reads.
    # HINT: Many terminal exons have multiple cleavage sites. In order to
    #       get a clean training set we need to filter out candidates that
    #       are in fact longer (3p crossing is not 0) or shorter
    #       (sum of the last 100 bases of the profile is not > 0 reads).
    # ---------------------------------------------------------------------
    annotation.filter_terminal_exon_training_candidates()

    # ---------------------------------------------------------------------
    # INTERMEDIATE EXONS - TRAINING SET CREATION
    # ---------------------------------------------------------------------
    # FIXME: this parts are weird...
    # 1. analyze_reads_for_annotated_regions
    #    AnalysisUnit.analyze_reads_for_annotated_regions()
    #    -> Writes a summary? where?
    # 2. directly afterwards they are deleted
    #
    ml_train_unit.create_intermediate_exon_training_set(
        intermediate_exons_bed_file_path=intermediate_exons_non_overlapping_filtered,
        sequencing_direction=options.sequencing_direction,
        max_splice_fuzziness=max_splice_fuzziness,
        output_dir=options.output_dir,
        genes_to_consider_dict=genes_to_consider_for_ML,
        bam_file_path=options.bam_file,
        annotation=annotation,
        threshold_to_filter=minimum_spliced_reads_for_cryptic_exon_start_site
    )

    # ---------------------------------------------------------------------
    # FILTER INTERMEDIATE EXONS TRAINING SET
    # Remove intermediate exons for which the profile is not > 0 reads.
    # This might happen in rare cases where the same exon is shared among
    # overlapping genes.
    # ---------------------------------------------------------------------
    annotation.filter_intermediate_exon_training_candidates()

    # ---------------------------------------------------------------------
    # Estimate gene expression
    # ---------------------------------------------------------------------
    estimate_expression_of_selected_genes(
        annotation=annotation,
        selected_genes_dict=annotation.genes,
        bam_file_path=options.bam_file,
        sequencing_direction=options.sequencing_direction,
        count_unique_mapping_reads_only=True,
        verbose=options.verbose
    )

    # ---------------------------------------------------------------------
    # Create a dataframe and supplement it with gene expression information
    # ---------------------------------------------------------------------

    pre_novel_terminal_output_file \
        = os.path.join(options.output_dir, "pre_novel_terminal_exons.tsv")

    ml_train_unit.create_terminal_exon_candidates_dataframe(
        annotation=annotation,
        novel_terminal_output_file=pre_novel_terminal_output_file,
        verbose=options.verbose
    )

    # pre_novel_terminal_filtered_output_file \
    #     = os.path.join(options.output_dir,
    #                    "pre_novel_terminal_exons_filtered.tsv")

    # # filter potential novel terminal exons that
    # # might overlap with annotated exons
    # ml_train_unit.remove_terminal_exon_candidates_that_overlap_annotated_exons(
    #     annotation=annotation,
    #     novel_terminal_output_file=pre_novel_terminal_filtered_output_file,
    #     sequencing_direction=options.sequencing_direction,
    #     exons=exons_per_gene_bed_file,
    #     verbose=False
    # )

    # create the terminal exon training set
    terminal_exons_file_core_name = "terminal_exons"
    terminal_exons_statistics_file_path = \
        os.path.join(
            options.output_dir,
            terminal_exons_file_core_name + ".tsv"
        )

    # create the intermediate exon training set
    intermediate_exons_file_core_name = "intermediate_exons"
    intermediate_exons_statistics_file_path = \
        os.path.join(
            options.output_dir,
            intermediate_exons_file_core_name + ".tsv"
        )

    # create the background training set
    background_regions_file_core_name = "background_regions"
    background_regions_statistics_file_path = \
        os.path.join(
            options.output_dir,
            background_regions_file_core_name + ".tsv"
        )

    # Generate N=10 training and validation datasets
    for run_number in range(10):

        sys.stdout.write("Generating #{} training and validation \
                          datasets ... {}".format(str(run_number),
                                                  os.linesep))

        # ---------------------------------------------------------------------
        # Create dataframes and supplement with gene expression information
        # ---------------------------------------------------------------------
        ml_train_unit.create_training_dataframes(
            annotation=annotation,
            intermediate_output_file=intermediate_exons_statistics_file_path,
            terminal_output_file=terminal_exons_statistics_file_path,
            background_output_file=background_regions_statistics_file_path,
            verbose=options.verbose
        )

        sys.stdout.write("After create_training_dataframes #{} \
                          ... {}".format(str(run_number),
                                         os.linesep))

        # ---------------------------------------------------------------------
        # Calculate features for training data
        # ---------------------------------------------------------------------
        ml_train_unit.add_features_to_training_dataframes(
            nr_to_subsample=1100,  # "all",
            output_files_dir=options.output_dir,
            verbose=options.verbose
        )

        sys.stdout.write("After add_features_to_training_dataframes #{} \
                          ...{}".format(str(run_number),
                                        os.linesep))

        # ---------------------------------------------------------------------
        # Sample training data
        # ---------------------------------------------------------------------
        ml_train_unit.load_training_data(
            training_data_set_size=1000,  # "max_equal_size",
            validation_data_fraction=0.2,
            output_files_dir=options.output_dir,
            run_number=run_number,
            verbose=options.verbose
        )

        sys.stdout.write("After load_training_data #{} \
                          ...{}".format(str(run_number),
                                        os.linesep))

        # Append to lists of training and validation datasets
        ml_train_unit.list_of_training_dataframes.append(
            ml_train_unit.training_df)
        ml_train_unit.list_of_validation_dataframes.append(
            ml_train_unit.validation_df)

        sys.stdout.write("After appending the lists #{} \
                          ... {}".format(str(run_number),
                                         os.linesep))

    # ---------------------------------------------------------------------
    # Calculate features
    # ---------------------------------------------------------------------
    ml_train_unit.add_features_to_terminal_exon_candidates_dataframe(
        output_files_dir=options.output_dir,
        verbose=options.verbose
    )

    sys.stdout.write("After add_features_to_terminal_exon_" +
                     "candidates_dataframe..." + os.linesep)

    # AJG: Go on from here!

    # ---------------------------------------------------------------------
    # Train a classifier
    # ---------------------------------------------------------------------
    sys.stdout.write("Training and validating classifiers..." + os.linesep)
    # ["Bayes", "KNeighbors", "multiclass_SVC"]
    for classifier in ["Bayes"]:

        # ml_train_unit.training_df = \
        #   ml_train_unit.list_of_training_dataframes[0]
        # ml_train_unit.validation_df = \
        #   ml_train_unit.list_of_validation_dataframes[0]

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
            " :: Used classifier: {} -> Results can be found in: {} {} \
            ".format(classifier,
                     classifier_results_dir,
                     os.linesep)
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
            verbose=options.verbose
        )

        # ---------------------------------------------------------------------
        # Sample training data
        # ---------------------------------------------------------------------
        ml_train_unit.load_training_data(
            training_data_set_size="max_equal_size",  # "max_equal_size",1000
            validation_data_fraction=0.2,
            output_files_dir=options.output_dir,
            run_number=-1,
            verbose=options.verbose
        )

        # -----------------------------------------------------------------
        # Select features using Greedy algorithm
        # -----------------------------------------------------------------
        # FIXME: hand over the default features
        ml_train_unit.greedy_feature_selection(
            number_of_randomization=10,
            classifier=classifier,
            manually_selected_features=manually_selected_features,
            verbose=options.verbose
        )

        # -----------------------------------------------------------------
        # Calculate features
        # -----------------------------------------------------------------
        # FIXME: remove the selected_features dependency from this function
        ml_train_unit.load_terminal_exon_candidates(verbose=options.verbose)

        for run_number in range(len(ml_train_unit.list_of_training_dataframes)):

            ml_train_unit.training_df = \
                ml_train_unit.list_of_training_dataframes[run_number]
            ml_train_unit.validation_df = \
                ml_train_unit.list_of_validation_dataframes[run_number]

            # -----------------------------------------------------------------
            # Train a classifier on the greedy selected features
            # -----------------------------------------------------------------
            ml_train_unit.train_classifier_on_selected_features(
                results_dir_path=options.output_dir,
                classifier=classifier,
                nr_of_train_vs_test_runs=25,
                verbose=options.verbose
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
            verbose=options.verbose
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
    sys.stdout.write("Reading genome file ..." + os.linesep)
    genome = Fasta(options.genome)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Annotate CDS for novel transcripts
    # -------------------------------------------------------------------------
    sys.stdout.write(
        "Annotating CDS for novel transcripts (novel splicing)..." +
        os.linesep)
    annotation.determine_CDS_novel_transcrtips_from_spicing(genome)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Write novel GTF file
    # -------------------------------------------------------------------------
    sys.stdout.write(
        "Writing GTF file for known and novel transcripts ..." + os.linesep)
    output_file = "enriched_annotation.gtf"
    output_file_path = os.path.join(options.output_dir, output_file)
    w = open(output_file_path, 'w')
    annotation.write_gtf(
        w,
        with_CDS=True,
        accepted_exons_dict=dict_of_validated_novel_exons
    )
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
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(-1)
    except Exception as exc:
        sys.stderr.write(exc.args[0] + os.linesep)
        sys.exit(-1)
