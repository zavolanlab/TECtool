# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import HTSeq
from collections import defaultdict
import hashlib
import numpy as np
from functools import partial
import os

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import our own modules
# -----------------------------------------------------------------------------

from aln_analysis import DetailedAlignment
from aln_analysis import SplitEvent
from aln_analysis import FeatureCounts

from gene_structure import Gene
from gene_structure import Transcript
from gene_structure import Exon

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# CLASSES
# -----------------------------------------------------------------------------

class AnalysisUnit:

    """
        Handling analysis unit class.

        :param unit_id: String of the containing feature id
        in the format: chr:start:end:strand:gene_name
        :param potential_5pSS_exons: List containing all
        potential 5'SS exons for the specific feature/unit_id
        :param gene_id: The gene_id

        :rtype: AnalysisUnit object

    *Class members*

        *annotation*
            String. Annotation of the analysis Unit.

        *norm_region*
            HTSeq Genomic Interval. The norm region of the analysis unit.

        *feat_region*
            HTSeq Genomic Interval. The feature region of the analysis unit.

        *iv*
            HTSeq Genomic Interval. The interval of the analysis unit.

        *gaos*
            GenomicArrayOfSets. The genomic array of sets for
            the norm, feature and iv region.

        *number_considered_reads*
            Int. Total number of reads that were considered in the analysis.

        *norm_region_count*
            Int. Number of reads falling in the norm region.

        *readthrough_count*
            Int. Number of reads supporting read throughs.

        *splice_into_feature_count*
            Int. Number of reads supporting novel splicing.

        *possible_exonic_start_sites*
            defaultdict(defaultdict(int)).  Dictionary that contains
            as key the possible start sites of the novel exon and the
            values is a new dictionary where the key  is the upstream 
            5'SS of an upstream exon and the values is how many times
            this spliced read is observed.

        *confident_novel_exons*
            defaultdict(list). Dictionary that contains as key the
            possible start sites of the novel exon and values is a list
            of sets where each set contains the  upstream 5'SS and the
            corresponding number of spliced reads that support the junction.

        *global_profile*
            HTSeq.GenomicArray. Per base profile of mapped reads.

        *splice_in*
            List of reads (detailed alignments) that splice
            in the feature region

        *splice_out*
            List of reads (detailed alignments) that splice out
            of the feature region

        *unspliced_reads*
            List of reads that are not spliced. Related to
            the feature region.

        *splice_in_reference*
            List of reads (detailed alignments) that splice in
            the reference region

        *splice_out_reference*
            List of reads (detailed alignments) that splice out
            of the reference region

        *unspliced_reads_reference*
            List of reads that are not spliced. Related to the
            reference region.

        *annotated_splice_out_all*
            Int. Number of splice out reads for the annotated
            region (intermediate or terminal exon)

        *annotated_splice_out_borders*
            Int. Number of splice out reads that fall in the
            borders of the annotated region (intermediate or terminal exon)

        *annotated_splice_in_all*
            Int. Number of splice in reads for
            the annotated region (intermediate or terminal exon)

        *annotated_splice_in_borders*
            Int. Number of splice in reads that fall in the
            borders of the annotated region (intermediate or terminal exon)

        *annotated_unspliced_5pSS*
            Int. Number of unspliced reads for the annotateed regions
            (intermediate or terminal exon) that cross the 5pSS

        *annotated_unspliced_3pSS*
            Int. Number of unspliced reads for the annotateed regions
            (intermediate or terminal exon) that cross the 3pSS

        *annotated_unspliced_feature*
            Int. Number of unspliced reads for the annotateed regions
            (intermediate or terminal exon) that fall within the borders
    """

    def __init__(self,
                 unit_id,
                 potential_5pSS_exons,
                 gene_id):

        # members THAT ARE CRUCIAL FOR INITIALIZATION
        self.unit_id = unit_id
        self.potential_5pSS_exons = potential_5pSS_exons
        self.gene_id = gene_id

        # Analysis unit annotation
        self.annotation = None

        # NORM REGION
        self.norm_region = None

        # FEAT REGION
        self.feat_region = None

        # Overall AnalysisUnit REGION
        self.iv = None
        self.gaos = None

        # READ COUNTS
        self.number_considered_reads = 0
        self.norm_region_count = 0
        self.readthrough_count = 0
        self.splice_into_feature_count = 0

        # Dictionary that contains as key the possible start sites of the novel
        # exon and the values is a new dictionary where the key  is the
        # upstream 5'SS of an upstream exon and the values is how many times
        # this spliced read is observed
        self.possible_exonic_start_sites = defaultdict(partial(defaultdict,
                                                               int))

        # Same as above but keep the possible exonic start sites as key
        # and values is a list of set where each set contains the
        # upstream 5'SS and the corresponding number of spliced reads
        # that support the junction
        self.confident_novel_exons = defaultdict(list)

        # Per base profile of mapped reads
        self.global_profile = HTSeq.GenomicArray("auto",
                                                 stranded=True,
                                                 typecode="i")

        # List of reads that splice in the feature region
        self.splice_in = []
        # List of reads that splice out from the feature region
        self.splice_out = []
        # List of reads that are not spliced. Related to the feature region.
        self.unspliced_reads = []

        # List of reads that splice in the reference region
        self.splice_in_reference = []
        # List of reads that splice out from the reference region
        self.splice_out_reference = []
        # List of reads that are not spliced. Related to the reference region.
        self.unspliced_reads_reference = []

        # Annotated exon stats
        self.annotated_splice_out_all = 0
        self.annotated_splice_out_borders = 0
        self.annotated_splice_in_all = 0
        self.annotated_splice_in_borders = 0
        self.annotated_unspliced_5pSS = 0
        self.annotated_unspliced_3pSS = 0
        self.annotated_unspliced_feature = 0

    def determine_norm_region(self):
        """
        Function that determines the norm region.
        It finds the exon that is closest to the feature region using the
        dictionary potential_5pSS_exons. In case there are multiple possible
        norm regions it selects the longest one.
        """

        min_start = None
        min_end = None
        exon_length = None

        if self.feat_region.strand is "+":

            for exon in self.potential_5pSS_exons:

                # This should be improved.
                # All potential_5pSS_exons should be an
                # HTSeq.GenomicInterval object.
                exon_splited = exon.split(":")

                if min_end is None:

                    min_start = int(exon_splited[1])
                    min_end = int(exon_splited[2])
                    exon_length = int(exon_splited[2]) - int(exon_splited[1])
                    continue

                if int(exon_splited[2]) >= int(min_end):

                    if int(exon_splited[2]) == int(min_end):

                        if int(exon_splited[2]) - int(exon_splited[1]) > exon_length:
                            min_end = int(exon_splited[2])
                            min_start = int(exon_splited[1])
                            exon_length = \
                                int(exon_splited[2]) - int(exon_splited[1])
                    else:

                        min_end = int(exon_splited[2])
                        min_start = int(exon_splited[1])
                        exon_length = \
                            int(exon_splited[2]) - int(exon_splited[1])

        elif self.feat_region.strand is "-":

            for exon in self.potential_5pSS_exons:

                # This should be improved.
                # All potential_5pSS_exons should be
                # an HTSeq.GenomicInterval object.
                exon_splited = exon.split(":")

                if min_start is None:

                    min_start = int(exon_splited[1])
                    min_end = int(exon_splited[2])
                    exon_length = int(exon_splited[2]) - int(exon_splited[1])
                    continue

                if int(exon_splited[1]) >= int(min_start):

                    if int(exon_splited[1]) == int(min_start):

                        if int(exon_splited[2]) - int(exon_splited[1]) > exon_length:

                            min_start = int(exon_splited[1])
                            min_end = int(exon_splited[2])
                            exon_length = \
                                int(exon_splited[2]) - int(exon_splited[1])

                else:
                    min_start = int(exon_splited[1])
                    min_end = int(exon_splited[2])
                    exon_length = \
                        int(exon_splited[2]) - int(exon_splited[1])

        else:

            sys.stderr.write("No strand info available for %s", self.unit_id)
            sys.exit(-1)

        self.norm_region = \
            HTSeq.GenomicInterval(exon_splited[0],
                                  min_start,
                                  min_end,
                                  exon_splited[3])

    def init_iv(self):
        """Initialize intervals"""

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Initialize feature region and determine norm region
        # ---------------------------------------------------------------------

        if self.unit_id is None or self.potential_5pSS_exons is None:
            sys.stderr.write("ERROR: unit_id or potential_5pSS_exons " +
                             "not initialized properly" + os.linesep)
            sys.exit(-1)

        # FEAT REGION
        feat_region_splited = str(self.unit_id).split(":")
        self.feat_region = HTSeq.GenomicInterval(feat_region_splited[0],
                                                 int(feat_region_splited[1]),
                                                 int(feat_region_splited[2]),
                                                 feat_region_splited[3])

        # NORM REGION
        self.determine_norm_region()

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Create an interval for the region of interest
        # ---------------------------------------------------------------------

        # calculate the entire region we are interested in
        aunit_start = min(self.norm_region.start, self.norm_region.end,
                          self.feat_region.start, self.feat_region.end)
        aunit_end = max(self.norm_region.start, self.norm_region.end,
                        self.feat_region.start, self.feat_region.end)

        # ensure that we talk about identical chromosomes
        if self.norm_region.chrom == self.feat_region.chrom:
            aunit_chrom = self.norm_region.chrom
        else:
            sys.stderr.write("ERROR: Chromosomes for AnalysisUnit " +
                             "object " + self.unit_id +
                             "are not identical!" + os.linesep)
            sys.exit(-1)

        # ensure that we talk about identical strands
        if self.norm_region.strand == self.feat_region.strand:
            aunit_strand = self.norm_region.strand
        else:
            sys.stderr.write("ERROR: Strands for AnalysisUnit " +
                             "object " + self.unit_id + "are not " +
                             "identical!" + os.linesep)
            sys.exit(-1)

        # Create an interval for the region of interest
        self.iv = HTSeq.GenomicInterval(aunit_chrom,
                                        aunit_start,
                                        aunit_end,
                                        aunit_strand)

    def init_gaos(self):
        """
        Initialize GenomicArrayOfSets
        """

        if self.norm_region is None or self.feat_region is None:
            sys.stderr.write("ERROR: Not initialized properly one" +
                             "of the regions in analysis unit: " +
                             self.unit_id + os.linesep)
            sys.exit(-1)

        # Create Genomic Array of Sets for the three regions
        self.gaos = HTSeq.GenomicArrayOfSets("auto",
                                             stranded=True)
        self.gaos[self.norm_region] += "norm_region"
        self.gaos[self.feat_region] += "feat_region"

    def analyze_reads(self,
                      bam,
                      sequencing_direction,
                      min_region_overlap,
                      splice_fuzziness,
                      minimum_spliced_reads_for_cryptic_exon_start_site,
                      count_unique_mapping_reads_only,
                      tmp,
                      annotation,
                      exons_per_gene_bed_file,
                      verbose=False):
        """
        Function that initializes the genomic interval and the
        genomic array of sets generates the region of interest,
        fetches the reads for the specific region and counts them.
        It then calls the functions to enrich the annotation with
        the novel exons and transcripts.
        """

        # initialize the entire genomic interval
        self.init_iv()

        # initialize the genomic array of sets
        self.init_gaos()

        # create the region string
        min_start = min(self.feat_region.start,
                        self.feat_region.end,
                        self.norm_region.start,
                        self.norm_region.end)

        max_end = max(self.feat_region.start,
                      self.feat_region.end,
                      self.norm_region.start,
                      self.norm_region.end)

        region_string = self.feat_region.chrom + \
            ":" + \
            str(min_start) + \
            "-" + \
            str(max_end)

        if verbose:
            sys.stdout.write("AnalysisUnit.iv (for fetching):\t" +
                             region_string + os.linesep)

            norm_region = self.norm_region.chrom + \
                ":" + \
                str(self.norm_region.start) + \
                "-" + \
                str(self.norm_region.end)

            sys.stdout.write("AnalysisUnit.norm_region:\t" +
                             norm_region + os.linesep)

            feat_region = self.feat_region.chrom + \
                ":" + \
                str(self.feat_region.start) + \
                "-" + \
                str(self.feat_region.end)

            sys.stdout.write("AnalysisUnit.feat_region:\t" +
                             feat_region + os.linesep)

        # go over each alignment
        for aln in bam.fetch(region=region_string):

            # In case we have unstranded data we change
            # the strand for the ALIGNMENT and the CIGAR STRING
            if "unstranded" in sequencing_direction:
                if self.feat_region.strand is not aln.iv.strand:
                    aln.iv.strand = self.feat_region.strand
                    for co in aln.cigar:
                        co.ref_iv.strand = self.feat_region.strand

            # Case reverse
            if "reverse" in sequencing_direction:
                if aln.iv.strand is "+":
                    aln.iv.strand = "-"
                    for co in aln.cigar:
                        co.ref_iv.strand = "-"
                elif aln.iv.strand is "-":
                    aln.iv.strand = "+"
                    for co in aln.cigar:
                        co.ref_iv.strand = "+"

            # count both of the reads
            self.count(
                aln,
                sequencing_direction=sequencing_direction,
                min_region_overlap=min_region_overlap,
                splice_fuzziness=splice_fuzziness,
                count_unique_mapping_reads_only=count_unique_mapping_reads_only,
                verbose=verbose
            )

        # Find novel cryptic exon coordinates and upstream 5'SS
        self.characterize_cryptic_exon_start_sites(
            bam,
            minimum_spliced_reads_for_cryptic_exon_start_site,
            tmp,
            sequencing_direction,
            exons_per_gene_bed_file
        )

        # Write novel exon statistics
        self.write_novel_exon_statistics(annotation)

        # Enrich annotation with the novel exons/transcripts
        self.enrich_annotation_novel_splicing(annotation, verbose)

    def analyze_reads_for_annotated_regions(self,
                                            bam,
                                            unit_id,
                                            sequencing_direction,
                                            splice_fuzziness,
                                            count_unique_mapping_reads_only,
                                            tmp,
                                            threshold_to_filter,
                                            feature_type,
                                            annotation,
                                            verbose=False):
        """
           Function that analyzes the reads for the
           annotated regions:
            - Annotated intermediate exons
            - Annotated terminal exons
            - Annotated background regions
        """

        # Initialize inverval
        self.iv = unit_id

        # Create Genomic Array of Sets for the regions
        self.gaos = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        self.gaos[self.iv] += "feat_region"

        # region string of the exon
        region_string = \
            unit_id.chrom + ":" + str(unit_id.start) + "-" + str(unit_id.end)

        # go over each alignment
        for aln in bam.fetch(region=region_string):

            # In case we have unstranded data we change the
            # strand for the ALIGNMENT and the CIGAR STRING
            if "unstranded" in sequencing_direction:
                if unit_id.strand is not aln.iv.strand:
                    aln.iv.strand = unit_id.strand
                    for co in aln.cigar:
                        co.ref_iv.strand = unit_id.strand

            self.count(
                aln,
                sequencing_direction=sequencing_direction,
                min_region_overlap=0,
                splice_fuzziness=splice_fuzziness,
                count_unique_mapping_reads_only=count_unique_mapping_reads_only,
                annotated=True,
                verbose=verbose
            )

        # Write annotated exon statistics
        self.write_annotated_exon_statistics(
            feature_type,
            annotation,
            threshold_to_filter
        )

    def estimate_gene_expression(
        self,
        bam,
        unit_id,
        sequencing_direction,
        count_unique_mapping_reads_only,
        annotation,
        verbose=False
    ):

        total_reads = 0

        # Initialize inverval
        self.iv = unit_id

        # region string of the gene
        region_string = \
            unit_id.chrom + ":" + str(unit_id.start) + "-" + str(unit_id.end)

        # go over each alignment
        for aln in bam.fetch(region=region_string):

            # -----------------------------------------------------------------
            # in case we want to count unique mapping reads only, we have to
            # exclude multimappers.
            # -----------------------------------------------------------------
            if (
                (count_unique_mapping_reads_only) and not
                (aln.optional_field("NH") == 1)
            ):
                continue

            # -----------------------------------------------------------------
            # Reads entirely fall in the gene coordinates
            # -----------------------------------------------------------------
            if (
                str(self.iv.chrom) == str(aln.iv.chrom) and
                int(self.iv.start) <= int(aln.iv.start) and
                int(self.iv.start) <= int(aln.iv.end) and
                int(self.iv.end) >= int(aln.iv.start) and
                int(self.iv.end) >= int(aln.iv.end)
            ):

                # -------------------------------------------------------------
                # count reads based on the orientation
                # -------------------------------------------------------------
                if sequencing_direction == "forward":
                    if self.iv.strand == aln.iv.strand:
                        total_reads += 1
                elif sequencing_direction == "unstranded":
                    total_reads += 1

        return(total_reads, self.gene_id)

    def write_novel_exon_statistics(self, annotation):
        """
        Write statistics for the novel exons
        """

        # loop over the novel exon
        for novel_exon in self.confident_novel_exons:

            # split the exon
            sp = novel_exon.split(":")

            # determine splice in reads (all and
            # the ones that fall exactly in the SS)
            splice_in_all = 0
            splice_in_borders = 0
            for splice_in in list(set(self.splice_in)):

                for event in splice_in.split_event_list:

                    if event.strand == '+':

                        # the 3'SS lies within the exon borders
                        # and the 5'SS is somewhere upstream
                        if (
                            int(event.three_prime_ss) >= int(sp[1]) and
                            int(event.three_prime_ss) <= int(sp[2]) and
                            int(event.five_prime_ss) < int(sp[1]) and
                            int(event.five_prime_ss) < int(sp[2])
                        ):

                            splice_in_all += 1

                        # the 3'SS lies exactly on the exon borders
                        # and the 5'SS is somewhere upstream
                        if int(event.three_prime_ss) == int(sp[1]):

                            splice_in_borders += 1

                    if event.strand == '-':

                        # same as above but for the minus strand

                        if (
                            int(event.three_prime_ss) >= int(sp[1]) and
                            int(event.three_prime_ss) <= int(sp[2]) and
                            int(event.five_prime_ss) > int(sp[1]) and
                            int(event.five_prime_ss) > int(sp[2])
                        ):

                            # this is a splice in reads
                            splice_in_all += 1

                        if int(event.three_prime_ss) == int(sp[2]):

                            splice_in_borders += 1

            # determine splice out reads
            splice_out_all = 0
            splice_out_borders = 0
            for splice_out in list(set(self.splice_out)):

                for event in splice_out.split_event_list:

                    # the five_prime_ss lies within the exon borders
                    # and the three_prime_ss is somewhere downstream

                    if event.strand == '+':

                        if (
                            int(event.five_prime_ss) >= int(sp[1]) and
                            int(event.five_prime_ss) <= int(sp[2]) and
                            int(event.three_prime_ss) > int(sp[2]) and
                            int(event.three_prime_ss) > int(sp[1])
                        ):

                            splice_out_all += 1

                        if int(event.five_prime_ss) == int(sp[2]):

                            splice_out_borders += 1

                    if event.strand == '-':

                        if (
                            int(event.five_prime_ss) >= int(sp[1]) and
                            int(event.five_prime_ss) <= int(sp[2]) and
                            int(event.three_prime_ss) < int(sp[1]) and
                            int(event.three_prime_ss) < int(sp[2])
                        ):

                            # This is a splice out read
                            splice_out_all += 1

                        if int(event.five_prime_ss) == int(sp[1]):

                            splice_out_borders += 1

            # count unspliced reads that either fall entirely in the
            # exon borders or they cross the 5' or 3' border

            reads_entirely_fall_in_exon = 0
            reads_fall_in_5p_border = 0
            reads_fall_in_3p_border = 0

            for det in self.unspliced_reads:

                if det.aln.iv.strand == '+':

                    # Find unspliced reads that entirely
                    # fall within the exon coordinates
                    if (
                        int(det.aln.iv.start) >= int(sp[1]) and
                        int(det.aln.iv.start) <= int(sp[2]) and
                        int(det.aln.iv.end) >= int(sp[1]) and
                        int(det.aln.iv.end) <= int(sp[2])
                    ):

                        reads_entirely_fall_in_exon += 1

                    # Find unspliced reads that crosss the 5' border
                    if (
                        int(det.aln.iv.start) <= int(sp[1]) and
                        int(det.aln.iv.end) >= int(sp[1])
                    ):

                        reads_fall_in_5p_border += 1

                    # Find unspliced reads that crosss the 3' border
                    if (
                        int(det.aln.iv.start) <= int(sp[2]) and
                        int(det.aln.iv.end) >= int(sp[2])
                    ):

                        reads_fall_in_3p_border += 1

                elif det.aln.iv.strand == '-':

                    # Find unspliced reads that entirely fall
                    # within the exon coordinates
                    if (
                        int(det.aln.iv.start) >= int(sp[1]) and
                        int(det.aln.iv.start) <= int(sp[2]) and
                        int(det.aln.iv.end) >= int(sp[1]) and
                        int(det.aln.iv.end) <= int(sp[2])
                    ):

                        reads_entirely_fall_in_exon += 1

                    # Find unspliced reads that crosss the 5' border
                    if (
                        int(det.aln.iv.start) <= int(sp[2]) and
                        int(det.aln.iv.end) >= int(sp[2])
                    ):

                        reads_fall_in_5p_border += 1

                    # Find unspliced reads that crosss the 3' border
                    if (
                        int(det.aln.iv.start) <= int(sp[1]) and
                        int(det.aln.iv.end) >= int(sp[1])
                    ):

                        reads_fall_in_3p_border += 1

            # Generate the profile for the novel exon
            exon_profile = list(str(x) for x in list(
                self.global_profile[
                    HTSeq.GenomicInterval(sp[0],
                                          int(sp[1]),
                                          int(sp[2]),
                                          sp[3])]))

            # reverse profile in case of minus strand
            if sp[3] == "-":
                exon_profile = exon_profile[::-1]

            region = ":".join([sp[0],
                               str(int(sp[1]) + 1),
                               sp[2],
                               sp[3]])
            profile = str(",".join(str(x) for x in exon_profile))

            region_feature = FeatureCounts(
                region=region,
                annotation=str("novel_splicing"),
                gene_id=str(self.gene_id),
                splice_in_all=int(splice_in_all),
                splice_in_borders=int(splice_in_borders),
                splice_out_all=int(splice_out_all),
                splice_out_borders=int(splice_out_borders),
                unspliced_feature=int(reads_entirely_fall_in_exon),
                unspliced_5pSS=int(reads_fall_in_5p_border),
                unspliced_3pSS=int(reads_fall_in_3p_border),
                profile=profile
            )

            annotation.genes[
                self.gene_id].potential_novel_exons.append(region_feature)

    def write_annotated_exon_statistics(
        self,
        feature_type,
        annotation,
        threshold_to_filter=5
    ):
        """
        Store counted features in the gene structure
        """

        chrom = str(self.iv.chrom)
        # this should be 0-based
        start = str(self.iv.start)
        end = str(self.iv.end)
        strand = str(self.iv.strand)

        exon_profile = list(str(x) for x in list(
            self.global_profile[HTSeq.GenomicInterval(
                chrom,
                int(start),
                int(end),
                strand)])
        )

        # reverse profile in case of minus strand
        if strand == "-":
            exon_profile = exon_profile[::-1]

        region = str(":".join([chrom,
                               str(int(start) + 1),
                               end,
                               strand]))
        profile = str(",".join(str(x) for x in exon_profile))

        region_feature = FeatureCounts(
            region=region,
            annotation=str(self.annotation),
            gene_id=str(self.gene_id),
            splice_in_all=int(self.annotated_splice_in_all),
            splice_in_borders=int(self.annotated_splice_in_borders),
            splice_out_all=int(self.annotated_splice_out_all),
            splice_out_borders=int(self.annotated_splice_out_borders),
            unspliced_feature=int(self.annotated_unspliced_feature),
            unspliced_5pSS=int(self.annotated_unspliced_5pSS),
            unspliced_3pSS=int(self.annotated_unspliced_3pSS),
            profile=profile
        )

        # We check if we have actual terminal exons or intermediate exons
        if int(self.annotated_splice_in_borders) >= int(threshold_to_filter):

            if feature_type == "terminal_exon":
                annotation.genes[
                    self.gene_id].annotated_terminal_exons.append(
                        region_feature)

            elif feature_type == "intermediate_exon":
                annotation.genes[
                    self.gene_id].annotated_intermediate_exons.append(
                        region_feature)

        # or if  we have background
        elif (
            int(self.annotated_splice_in_borders) > 0 and
            int(self.annotated_splice_in_borders) < int(threshold_to_filter)
        ):

            if feature_type == "terminal_exon":
                annotation.genes[
                    self.gene_id].background.append(
                        region_feature)

    # _________________________________________________________________________
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Method: Enrich annotation
    # -------------------------------------------------------------------------
    def enrich_annotation_novel_splicing(self, annotation, verbose):
        """
        Function that enriches the annotation objects with the novel
        exons/transcsripts that were detected based on novel splicing
        """

        # Loop over the confident novel exons
        for novel_exon in self.confident_novel_exons:

            # Each novel exon might have more than one upstream
            # exons that can splice to.
            # For this reason for each novel exon we loop as
            # many times as the upstream exons that splices to.
            for up in self.confident_novel_exons[novel_exon]:

                # For the novel exon and the specific upstream splice site
                # we identify all the possible upstream exons.
                upstream_exons = []

                for current_feature_region in annotation.feature_regions_upstream_coordinates:

                    # Loop over the possible exons
                    for upstream_exon in annotation.feature_regions_upstream_coordinates[current_feature_region]:

                        if upstream_exon.split(":")[3] is "+":

                            # Check if the 5'SS is the same
                            if upstream_exon.split(":")[2] == up[0].split(":")[2]:

                                if upstream_exon not in upstream_exons:

                                    upstream_exons.append(upstream_exon)

                        elif upstream_exon.split(":")[3] is "-":

                            # Check if the 5'SS is the same
                            if upstream_exon.split(":")[1] == up[0].split(":")[2]:

                                if upstream_exon not in upstream_exons:

                                    upstream_exons.append(upstream_exon)

                # Now that we have the upstream exons, we find annotated
                # transcripts that contain these exons

                # loop over the upstream exons
                for exon in upstream_exons:

                    if exon in annotation.genes[self.gene_id
                                                ].exon_coordinates_dict:

                        # find exon ids for the exon coordinates
                        for exon_id in annotation.genes[self.gene_id].exon_coordinates_dict[exon]:

                            # Get transcripts of the gene
                            for transcript in annotation.genes[self.gene_id].get_known_transctipt_ids():

                                # Split exon
                                novel_exon_sp = novel_exon.strip().split(":")

                                # make sure that the novel exon is not
                                # found after the transcript end
                                if (
                                    (
                                        (annotation.transcripts[
                                            transcript].strand == '+') and
                                        (int(novel_exon_sp[1]) <
                                            int(annotation.transcripts[
                                                transcript].end)) and
                                        (int(novel_exon_sp[2]) <
                                            int(annotation.transcripts[
                                                transcript].end))
                                    ) or
                                    (
                                        (annotation.transcripts[
                                            transcript].strand == '-') and
                                        (int(novel_exon_sp[1]) >
                                            int(annotation.transcripts[
                                                transcript].start)) and
                                        (int(novel_exon_sp[2]) >
                                            int(annotation.transcripts[
                                                transcript].start))
                                    )
                                ):

                                    # Loop over the exons of each transcript
                                    for transcript_exon in annotation.transcripts[transcript].exon_list_sorted_by_end_coord:

                                        # If exon exists in this transcript
                                        if exon_id == transcript_exon.exon_id:

                                            if annotation.transcripts[
                                                transcript].strand == "+":

                                                # Create novel transcript annotation
                                                exons_before_novel_one = annotation.transcripts[
                                                    transcript].get_existing_and_upstream_exons(
                                                        exon.split(":")[1],
                                                        exon.split(":")[2],
                                                        exon.split(":")[3]
                                                    )

                                                # In case the list is not empty
                                                if len(exons_before_novel_one) > 0:

                                                    exons_before_novel_one_name = []
                                                    for x in exons_before_novel_one:
                                                        exons_before_novel_one_name.append(":".join([str(x.start), str(x.end)]))

                                                    # create unique transcript id based on all the exons found before the novel one
                                                    upstream_coords = hashlib.md5(('_'.join(exons_before_novel_one_name).encode('utf-8'))).hexdigest()
                                                    novel_transcript_id = ""
                                                    novel_transcript_id += "novel_"+annotation.transcripts[transcript].gene_id
                                                    novel_transcript_id += "|UE_"+upstream_coords
                                                    novel_transcript_id += "|5pSS_"+str(str(novel_exon.split(":")[1]))
                                                    novel_transcript_id += "|PAS_"+str(str(novel_exon.split(":")[2]))

                                                    # create transcript object
                                                    novelTranscript = Transcript(
                                                        chromosome=annotation.transcripts[transcript].chromosome,
                                                        source="TECtool_annotated",
                                                        feature="transcript",
                                                        start=str(annotation.transcripts[transcript].start),
                                                        end=novel_exon.split(":")[2],
                                                        score=annotation.transcripts[transcript].score,
                                                        strand="+",
                                                        frame=".",
                                                        gene_id=annotation.transcripts[transcript].gene_id,
                                                        transcript_id=novel_transcript_id,
                                                        gene_name=annotation.transcripts[transcript].gene_name,
                                                        gene_biotype=annotation.transcripts[transcript].gene_biotype,
                                                        transcript_name=novel_transcript_id,
                                                        transcript_biotype="novel_splicing"
                                                    )

                                                    # Add the transcript in the list of novel transcripts
                                                    annotation.genes[self.gene_id].novel_transcripts.append(novelTranscript)

                                                    # Store to dictionary key :
                                                    # novel transcript id
                                                    # value: [transcipt id 1,
                                                    # transcipt id 2] the
                                                    # potential transcripts
                                                    # that a novel transcript
                                                    # can originate from
                                                    annotation.genes[self.gene_id].mother_transcripts_of_novel_transcripts.setdefault(novel_transcript_id, []).append(transcript)

                                                    # Create novel exon
                                                    # annotation
                                                    exon_count = 1
                                                    for gtf_exon in annotation.transcripts[transcript].exon_list_sorted_by_end_coord:

                                                        if gtf_exon.end < int(novel_exon.split(":")[1]):

                                                            novelExon = Exon(
                                                                chromosome = gtf_exon.chromosome,
                                                                source="TECtool_annotated",
                                                                feature="exon",
                                                                start=gtf_exon.start,
                                                                end=gtf_exon.end,
                                                                score=".",
                                                                strand="+",
                                                                frame=gtf_exon.frame,
                                                                gene_id=self.gene_id,
                                                                transcript_id=novel_transcript_id,
                                                                exon_number=str(exon_count),
                                                                gene_name=annotation.transcripts[transcript].gene_name,
                                                                gene_biotype=annotation.transcripts[transcript].gene_biotype,
                                                                transcript_name=novel_transcript_id,
                                                                transcript_biotype='novel_splicing',
                                                                exon_id="novel_"+exon_id+"_"+novel_transcript_id
                                                            )

                                                            novelExon.CDS = gtf_exon.CDS
                                                            novelExon.start_codon = gtf_exon.start_codon
                                                            novelExon.stop_codon = gtf_exon.stop_codon

                                                            novelTranscript.novel_exons.append(novelExon)

                                                            exon_count += 1

                                                    novelExon = Exon(
                                                        chromosome = novel_exon.split(":")[0],
                                                        source="TECtool_annotated",
                                                        feature="exon",
                                                        start=str(novel_exon.split(":")[1]),
                                                        end=str(novel_exon.split(":")[2]),
                                                        score=".",
                                                        strand="+",
                                                        frame=None,
                                                        gene_id=self.gene_id,
                                                        transcript_id=novel_transcript_id,
                                                        exon_number=str(exon_count),
                                                        gene_name=annotation.transcripts[transcript].gene_name,
                                                        gene_biotype=annotation.transcripts[transcript].gene_biotype,
                                                        transcript_name=novel_transcript_id,
                                                        transcript_biotype='novel_splicing',  # novel_splicing_exon_last
                                                        exon_id="novel_terminal_exon_"+novel_transcript_id
                                                    )

                                                    novelExon.CDS = None
                                                    novelExon.start_codon = None
                                                    novelExon.stop_codon = None

                                                    novelTranscript.novel_exons.append(novelExon)

                                            elif annotation.transcripts[transcript].strand == "-":

                                                exons_before_novel_one = annotation.transcripts[transcript].get_existing_and_upstream_exons(exon.split(":")[1], exon.split(":")[2], exon.split(":")[3])

                                                # In case the list is not empty
                                                if len(exons_before_novel_one) > 0:

                                                    exons_before_novel_one_name = []
                                                    for x in exons_before_novel_one:
                                                        exons_before_novel_one_name.append(":".join([str(x.start), str(x.end)]))

                                                    # create unique transcript id based on all the exons found before the novel one
                                                    upstream_coords = hashlib.md5(('_'.join(exons_before_novel_one_name).encode('utf-8'))).hexdigest()
                                                    novel_transcript_id =  ""
                                                    novel_transcript_id += "novel_"+annotation.transcripts[transcript].gene_id
                                                    novel_transcript_id += "|UE_"+upstream_coords
                                                    novel_transcript_id += "|5pSS_"+str(str(novel_exon.split(":")[2]))
                                                    novel_transcript_id += "|PAS_"+str(str(novel_exon.split(":")[1]))

                                                    novelTranscript = Transcript(
                                                        chromosome=annotation.transcripts[transcript].chromosome,
                                                        source="TECtool_annotated",
                                                        feature="transcript",
                                                        start=str(novel_exon.split(":")[1]),
                                                        end=str(annotation.transcripts[transcript].end),
                                                        score=annotation.transcripts[transcript].score,
                                                        strand="-",
                                                        frame=".",
                                                        gene_id=annotation.transcripts[transcript].gene_id,
                                                        transcript_id=novel_transcript_id,
                                                        gene_name=annotation.transcripts[transcript].gene_name,
                                                        gene_biotype=annotation.transcripts[transcript].gene_biotype,
                                                        transcript_name=novel_transcript_id,
                                                        transcript_biotype="novel_splicing"
                                                    )

                                                    # Add the transcript in the list of novel transcripts
                                                    annotation.genes[self.gene_id].novel_transcripts.append(novelTranscript)

                                                    # Store to dictionary key : novel transcript id value: [transcipt id 1, transcipt id 2]  
                                                    # the  potential transcripts that a novel transcript can originate from 
                                                    annotation.genes[self.gene_id].mother_transcripts_of_novel_transcripts.setdefault(novel_transcript_id, []).append(transcript)

                                                    # Create novel exon annotation
                                                    exon_count =  1
                                                    # for gtf_exon in annotation.transcripts[transcript].exon_list_sorted_by_end_coord[::-1]:
                                                    for gtf_exon in exons_before_novel_one:

                                                        if(gtf_exon.start >= int(novel_exon.split(":")[2])): # make sure that this check is enough ... Well IT IS NOT....

                                                            novelExon = Exon(
                                                                chromosome = novel_exon.split(":")[0],
                                                                source="TECtool_annotated",
                                                                feature="exon",
                                                                start=gtf_exon.start,
                                                                end=gtf_exon.end,
                                                                score=".",
                                                                strand="-",
                                                                frame=gtf_exon.frame,
                                                                gene_id=self.gene_id,
                                                                transcript_id=novel_transcript_id,
                                                                exon_number=str(exon_count),
                                                                gene_name=annotation.transcripts[transcript].gene_name,
                                                                gene_biotype=annotation.transcripts[transcript].gene_biotype,
                                                                transcript_name=novel_transcript_id,
                                                                transcript_biotype='novel_splicing',
                                                                exon_id="novel_"+exon_id+"_"+novel_transcript_id
                                                            )


                                                            novelExon.CDS = gtf_exon.CDS
                                                            novelExon.start_codon = gtf_exon.start_codon
                                                            novelExon.stop_codon = gtf_exon.stop_codon

                                                            novelTranscript.novel_exons.append(novelExon)

                                                            exon_count += 1

                                                    novelExon = Exon(
                                                        chromosome = novel_exon.split(":")[0],
                                                        source="TECtool_annotated",
                                                        feature="exon",
                                                        start=str(novel_exon.split(":")[1]),
                                                        end=str(novel_exon.split(":")[2]),
                                                        score=".",
                                                        strand="-",
                                                        frame=None,
                                                        gene_id=self.gene_id,
                                                        transcript_id=novel_transcript_id,
                                                        exon_number=str(exon_count),
                                                        gene_name=annotation.transcripts[transcript].gene_name,
                                                        gene_biotype=annotation.transcripts[transcript].gene_biotype,
                                                        transcript_name=novel_transcript_id,
                                                        transcript_biotype='novel_splicing', # novel_splicing_exon_last
                                                        exon_id="novel_terminal_exon_"+novel_transcript_id
                                                    )

                                                    novelExon.CDS = None
                                                    novelExon.start_codon = None
                                                    novelExon.stop_codon = None

                                                    novelTranscript.novel_exons.append(novelExon)

                                            else:
                                                stderr.write("[ERROR] Problem with strand info.")
                                                sys.exit(-1)

    def keep_novel_exon(self, profile_list, threshold=0.2):
        """
        Function that decides to keep or discard the novel exon based on
        coverage of the possible new exon
        """

        max_list = max(list(profile_list))
        profile_list_norm = \
            np.array(list(profile_list)) / float(max_list)
        above = (profile_list_norm >= threshold).sum()
        below = (profile_list_norm < threshold).sum()

        if above > below:
            return(True)
        else:
            return(False)

    # _________________________________________________________________________
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Method: Determine cryptic/novel exon start site
    # -------------------------------------------------------------------------
    def characterize_cryptic_exon_start_sites(
        self,
        bam,
        minimum_spliced_reads_for_cryptic_exon_start_site,
        tmp,
        sequencing_direction,
        exons_per_gene_bed_file
    ):
        """
        Find novel cryptic exon coordinates and upstream 5'SS
        """

        # Loop over the posible splice sites and start filtering
        for possible_exonic_start_site in self.possible_exonic_start_sites:

            for possible_5pSS in self.possible_exonic_start_sites[possible_exonic_start_site]:

                # if the spliced reads overpasses the user defined threshold
                if (
                    int(self.possible_exonic_start_sites[
                        possible_exonic_start_site][possible_5pSS]) >=
                        int(minimum_spliced_reads_for_cryptic_exon_start_site)
                ):

                    novel_start_site_info = \
                        possible_exonic_start_site.split(":")
                    # One more thing we have to check is that the novel
                    # exon is not falling in the region of an annotated
                    # overlapping gene. For this we have to make sure
                    # that the start site of the novel exon is not
                    # overlapping with any other exon
                    if novel_start_site_info[3] is "+":
                        novel_exon = ":".join(
                            [novel_start_site_info[0],
                             str(novel_start_site_info[1]),
                             str(self.feat_region.end),
                             novel_start_site_info[3]]
                        )

                        counts = str(self.possible_exonic_start_sites[
                            possible_exonic_start_site][possible_5pSS])
                        self.confident_novel_exons[novel_exon].append(
                            [possible_5pSS,
                             counts,
                             self.feat_region]
                        )

                    elif novel_start_site_info[3] is "-":
                        novel_exon = ":".join(
                            [novel_start_site_info[0],
                             str(self.feat_region.start),
                             str(novel_start_site_info[2]),
                             novel_start_site_info[3]]
                        )

                        counts = str(self.possible_exonic_start_sites[
                            possible_exonic_start_site][possible_5pSS])
                        self.confident_novel_exons[novel_exon].append(
                            [possible_5pSS,
                             counts,
                             self.feat_region]
                        )

                    else:
                        sys.stderr.write("[ERROR] No strand info available")
                        sys.exit(-1)

    # _________________________________________________________________________
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Method: count
    # -------------------------------------------------------------------------
    def count(
        self,
        aln,
        sequencing_direction,
        min_region_overlap,
        splice_fuzziness,
        count_unique_mapping_reads_only=True,
        annotated=False,
        verbose=False
    ):
        """
        Function that counts the number of reads that fall in each region
        of interest
        """

        # ---------------------------------------------------------------------
        # in case we want to count unique mapping reads only, we have to
        # exclude multimappers.
        if (
            count_unique_mapping_reads_only and
            not aln.optional_field("NH") == 1
        ):
            return

        # ---------------------------------------------------------------------
        # NOTE:
        # We could discard all reads that have an insertion or deletion
        # somewhere within the read. But that would also mean, that if we
        # have sequenced an individual that has a SNP compared to the
        # standard genome, this regions will be wrongly analyzed!!
        # -> Therefore: it makes sense to build up genomic intervals,
        #               at least in case the insertion/deletion/mutation
        #               is not too strong (which should already be ensured
        #               by running the aligner with the right settings).

        # ---------------------------------------------------------------------
        # we cannot treat paired-end mappings
        if aln.paired_end or not aln.pe_which == "not_paired_end":
            sys.stderr.write("ERROR: Paired-end mapping " +
                             "found, but not supported! " +
                             os.linesep)
            sys.exit(-1)

        # -----------------------------------------------------------------
        # in case we consider the read as valid, we have to recognize this
        self.number_considered_reads += 1

        # -----------------------------------------------------------------
        # give some feedback (if wanted)
        if verbose:

            sys.stdout.write("_" * 80 + os.linesep)

            sys.stdout.write("Read:\t\t" +
                             aln.read.name +
                             "\tpe_which: " +
                             aln.pe_which +
                             os.linesep)

            aln_string = \
                aln.iv.chrom + \
                ":" + \
                str(aln.iv.start) + \
                "-" + \
                str(aln.iv.end) + \
                str(aln.iv.strand)

            sys.stdout.write("Alignment:\t" + aln_string + "\n")

        # -----------------------------------------------------------------
        # check whether we have mappings that are supported by this script
        if (
            ("forward" not in sequencing_direction) and
            ("unstranded" not in sequencing_direction)
        ):
            sys.stderr.write("ERROR: Forward, reverse and unstranded" +
                             "mappings supported!" + os.linesep)
            sys.exit(-1)

        # _________________________________________________________________
        # -----------------------------------------------------------------
        # ANALYZE THE CIGAR STRING
        # -----------------------------------------------------------------
        # Now we need to loop over the cigar string and store all
        # details of the alignment (within a DetailedAlignment object).
        #
        # After we have seen everything(!!), we can decide how to count
        # this read/alignment.
        # -----------------------------------------------------------------
        deal = DetailedAlignment(aln=aln)
        current_genomic_interval = None
        current_split_event = None

        # check to which strand the alignment maps to
        if aln.iv.strand == "+":
            cigar_list = aln.cigar
        elif aln.iv.strand == "-":
            cigar_list = reversed(aln.cigar)
        else:
            sys.stderr.write("ERROR: unknown strand in alignment!\n")
            sys.exit(-1)

        # Loop over the cigar string
        for cig_op in cigar_list:

            # give some feedback if wanted
            if verbose:
                sys.stdout.write("CIGAR type: " +
                                 cig_op.type +
                                 " -> coords: " +
                                 str(cig_op.ref_iv.chrom) +
                                 ":" +
                                 str(cig_op.ref_iv.start) +
                                 "-" + str(cig_op.ref_iv.end) +
                                 os.linesep)

            # Case: we found S -> just remember how many we have seen
            if cig_op.type == "S":
                deal.number_of_S += cig_op.size

            # Case: we found M -> set as current genomic interval or add it
            #                     to the current genomic interval
            elif cig_op.type == "M":

                # Case: we have not seen matches before
                if current_genomic_interval is None:
                    current_genomic_interval = cig_op.ref_iv

                # Case: we have seen matches before,
                # that might have been interupted by insertions
                # or deletions
                else:
                    current_genomic_interval.extend_to_include(cig_op.ref_iv)

                # -------------------------------------------------------------
                # Create profile for the matching regions
                # -------------------------------------------------------------
                self.global_profile[cig_op.ref_iv] += 1

            # Case: we have a split region
            elif cig_op.type == "N":

                # -------------------------------------------------------------
                # Case: if we have not seen any split before, we open a new one
                #       (5' splice site found)
                if current_split_event is None:

                    # we have to care about the strand here
                    if aln.iv.strand == "+":
                        current_split_event = SplitEvent(
                            chrom=current_genomic_interval.chrom,
                            strand=current_genomic_interval.strand,
                            five_prime_ss=current_genomic_interval.end
                        )
                    elif aln.iv.strand == "-":
                        current_split_event = SplitEvent(
                            chrom=current_genomic_interval.chrom,
                            strand=current_genomic_interval.strand,
                            five_prime_ss=current_genomic_interval.start
                        )

                # -------------------------------------------------------------
                # Case: if we have seen a split event before, we close the
                #       one that is still open, since now we know the end
                #       of it (3' splice site found)
                else:

                    # we have to care about the strand here
                    if aln.iv.strand == "+":
                        current_split_event.three_prime_ss = \
                            current_genomic_interval.start
                    elif aln.iv.strand == "-":
                        current_split_event.three_prime_ss = \
                            current_genomic_interval.end

                    # add the SplitEvent object to the list of split events
                    deal.split_event_list.append(current_split_event)

                    # set the current split event None (not necessary, but we
                    # want to be clean here).
                    current_split_event = None

                    # We open another split event since we see another N
                    if aln.iv.strand == "+":
                        current_split_event = SplitEvent(
                            chrom=current_genomic_interval.chrom,
                            strand=current_genomic_interval.strand,
                            five_prime_ss=current_genomic_interval.end
                        )
                    elif aln.iv.strand == "-":
                        current_split_event = SplitEvent(
                            chrom=current_genomic_interval.chrom,
                            strand=current_genomic_interval.strand,
                            five_prime_ss=current_genomic_interval.start
                        )

                # -------------------------------------------------------------
                # ANALYZE
                # -------------------------------------------------------------
                # since we know that a genomic interval ends here...
                # 1. we count it
                # 2. we set the current_genomic_interval to None (so that we
                #    can open another one in case we see another match)
                self.analyze_genomic_interval(
                    genomic_interval=current_genomic_interval,
                    detailed_alignment=deal,
                    min_region_overlap=min_region_overlap,
                    verbose=verbose
                )
                current_genomic_interval = None

        # ---------------------------------------------------------------------
        # ANALYZE
        # ---------------------------------------------------------------------
        # since we know that a genomic interval ends here...
        # 1. we count it
        # 2. we set the current_genomic_interval to None (so that we
        #    can open another one in case we see another match)
        if current_genomic_interval is not None:

            # if there exists still an open split event, we need to close it
            if current_split_event is not None:

                    # we have to care about the strand here
                    if aln.iv.strand == "+":
                        current_split_event.three_prime_ss = \
                            current_genomic_interval.start
                    elif aln.iv.strand == "-":
                        current_split_event.three_prime_ss = \
                            current_genomic_interval.end

                    # add the SplitEvent object to the list of split events
                    deal.split_event_list.append(current_split_event)

                    # set the current split event None (not necessary, but we
                    # want to be clean here).
                    current_split_event = None

            # finally, count the current genomic interval
            self.analyze_genomic_interval(
                genomic_interval=current_genomic_interval,
                detailed_alignment=deal,
                min_region_overlap=min_region_overlap,
                verbose=verbose
            )
            current_genomic_interval = None

        else:
            sys.stderr.write("ERROR: Weird case, check it!\n")
            sys.exit(-1)

        # case where try to identify novel exons
        if not annotated:

            # -----------------------------------------------------------------
            # COUNT
            # -----------------------------------------------------------------
            # Finally, knowing all details about the current alignment (deal)
            # we can accurately count the alignment for each region.
            # -----------------------------------------------------------------
            self.count_detailed_alignment(
                detailed_alignment=deal,
                splice_fuzziness=splice_fuzziness,
                verbose=verbose
            )
        # case where try to count annotated exons (terminal or intermediate)
        elif annotated:

            # -----------------------------------------------------------------
            # COUNT
            # -----------------------------------------------------------------
            # Finally, knowing all details about the current alignment (deal)
            # we can accurately count the alignment for each region.
            # -----------------------------------------------------------------
            self.count_detailed_alignment_annotated(
                detailed_alignment=deal,
                verbose=verbose
            )

    # _________________________________________________________________________
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Method: analyze_genomic_interval
    # -------------------------------------------------------------------------
    def analyze_genomic_interval(
        self,
        genomic_interval,
        detailed_alignment,
        min_region_overlap,
        verbose=False
    ):
        """
        Analyze Genomic Interval
        """

        if verbose:
            sys.stdout.write("-" * 80 + os.linesep)
            sys.stdout.write(
                "Genomic interval:\t" +
                genomic_interval.chrom + ":" +
                str(genomic_interval.start) +
                "-" +
                str(genomic_interval.end) +
                os.linesep
            )

        # ---------------------------------------------------------------------
        # Create a set in which we store the regions that intersect with
        # the current genomic interval
        # ---------------------------------------------------------------------
        gi_regions_set = set()

        # ---------------------------------------------------------------------
        # loop over the regions that intersect with the genomic interval
        for region, intersection_set in self.gaos[genomic_interval].steps():

            # give some feedback if wanted
            if verbose:
                sys.stdout.write(
                    str(region) +
                    "\t" +
                    str(intersection_set) +
                    os.linesep
                )

            # Case: the genomic interval intersects with
            #       at least one of the defined regions -> everything is fine
            if len(intersection_set) == 1:
                # everything is fine, we can count
                if verbose:
                    sys.stdout.write("COUNT! -> Found NON-EMPTY SET of " +
                                     "length 1 intersecting with the " +
                                     " genomic interval!" +
                                     os.linesep)

                # we take an overlap only into account, in case it is more than
                # a minimum overlap!
                if region.length >= min_region_overlap:
                    covered_region = next(iter(intersection_set))
                    detailed_alignment.regions_set.add(covered_region)
                    gi_regions_set.add(covered_region)

            # Case: the genomic interval intersects with several regions
            #       of the gaos (that should not happen, if the regions
            #       have been defined correctly).
            elif len(intersection_set) > 1:
                sys.stderr.write("ERROR: norm_region and feat_region " +
                                 "seemto overlap!" + os.linesep)
                sys.exit(-1)

            # Case: the genomic interval intersects with the chromosome outside
            #       of the defined gaos -> do nothing, but report it if wanted.
            elif len(intersection_set) == 0:
                if verbose:
                    sys.stdout.write("Found EMPTY SET intersecting " +
                                     "with the genomic interval!" +
                                     os.linesep)

            # Case: this should never happen
            else:
                sys.stderr.write("ERROR: Length of set smaller 0!\n")
                sys.exit(-1)

        # ---------------------------------------------------------------------
        # Case: after looping over all regions that intersect with the
        #       genomic_interval we can tell wether it spans the boarder.
        if len(gi_regions_set) > 1:
            detailed_alignment.spans_regions_boarder = True

    # _________________________________________________________________________
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Method: count_detailed_alignment
    # -------------------------------------------------------------------------
    def count_detailed_alignment(self,
                                 detailed_alignment,
                                 splice_fuzziness=3,
                                 verbose=False):

        # count
        for region in detailed_alignment.regions_set:

            # # Case: count for the norm region
            # if region == "norm_region":
            #     self.norm_region_count += 1

            # Case: count for the feature region (either split or readthrough)
            if region == "feat_region":

                # -------------------------------------------------------------
                # New addition for counting features:
                # Splice in, splice out, readthroughs
                # -------------------------------------------------------------

                # Case: no spliced reads
                if len(detailed_alignment.split_event_list) == 0:

                    self.unspliced_reads.append(detailed_alignment)

                # Case: spliced reads
                elif len(detailed_alignment.split_event_list) > 0:

                    for split_event in detailed_alignment.split_event_list:

                        # -----------------------------------------------------
                        # PLUS STRAND.
                        # -----------------------------------------------------
                        if detailed_alignment.aln.iv.strand == "+":

                            # -------------------------------------------------
                            # CHECKS
                            # -------------------------------------------------

                            # The 3' splice site is in the feature region and
                            # the 5' splice site is located somewhere
                            # upstream (Splice in case)
                            if (
                                int(split_event.three_prime_ss) >=
                                int(self.feat_region.start) and
                                int(split_event.three_prime_ss) <=
                                int(self.feat_region.end) and
                                int(split_event.five_prime_ss) <=
                                int(self.feat_region.start) and
                                int(split_event.five_prime_ss) <=
                                int(self.feat_region.end)
                            ):

                                self.splice_in.append(detailed_alignment)

                            # The 5' splice site is in the feature region and
                            # the 3' splice site is located somewhere
                            # downstream (Splice out case)

                            if (
                                int(split_event.five_prime_ss) >=
                                int(self.feat_region.start) and
                                int(split_event.five_prime_ss) <=
                                int(self.feat_region.end) and
                                int(split_event.three_prime_ss) >=
                                int(self.feat_region.end) and
                                int(split_event.three_prime_ss) >=
                                int(self.feat_region.start)
                            ):

                                self.splice_out.append(detailed_alignment)

                        elif detailed_alignment.aln.iv.strand == "-":

                            # Same as above but for the minus strand

                            # splice in case
                            if (
                                int(split_event.three_prime_ss) >=
                                int(self.feat_region.start) and
                                int(split_event.three_prime_ss) <=
                                int(self.feat_region.end) and
                                int(split_event.five_prime_ss) >=
                                int(self.feat_region.start) and
                                int(split_event.five_prime_ss) >=
                                int(self.feat_region.end)
                            ):

                                self.splice_in.append(detailed_alignment)

                            # splice out case
                            if (
                                int(split_event.five_prime_ss) >=
                                int(self.feat_region.start) and
                                int(split_event.five_prime_ss) <=
                                int(self.feat_region.end) and
                                int(split_event.three_prime_ss) <=
                                int(self.feat_region.start) and
                                int(split_event.three_prime_ss) <=
                                int(self.feat_region.end)
                            ):

                                self.splice_out.append(detailed_alignment)

                # -------------------------------------------------------------
                # Case: clear readthrough alignment
                #       can also contain a split (e.g. before the readthrough)
                if detailed_alignment.spans_regions_boarder:
                    self.readthrough_count += 1

                # -------------------------------------------------------------
                # Case: the alignment does not span the regions boarder
                #       and there have been splice mappings, we have to check
                #       whether the splice event falls within the norm and
                #       the feature region.
                elif (
                    not detailed_alignment.spans_regions_boarder and
                    len(detailed_alignment.split_event_list) > 0
                ):

                    # ---------------------------------------------------------
                    # Idea: go over all split events and check whether one of
                    #       them starts in the norm region (5' splice site)
                    #       and ends (3' splice site) in the feature region.
                    #       Count such a read, in case the read ends within the
                    #       feature region (if not, the read would not be
                    #       an indication for a terminal exon.

                    for split_event in detailed_alignment.split_event_list:

                        # give some feedback to the user, if wanted
                        if verbose:
                            sys.stdout.write(
                                "SPLIT EVENT: " +
                                split_event.chrom +
                                ":" +
                                str(split_event.five_prime_ss) +
                                "-" +
                                str(split_event.three_prime_ss) +
                                os.linesep
                            )

                        # -----------------------------------------------------
                        # PLUS STRAND.
                        # -----------------------------------------------------
                        if detailed_alignment.aln.iv.strand == "+":

                            # Go over all upstream exon coordinates and check
                            # if the split reads starts in one of these exons
                            if verbose:
                                sys.stdout.write("... loop over all" +
                                                 " possible upstream" +
                                                 " exons ..." +
                                                 os.linesep)

                            observed_5pSS = []
                            # loop over all possible 5'SSs
                            for exon_coordinate in self.potential_5pSS_exons:

                                exon_coordinate_splited = \
                                    exon_coordinate.split(":")

                                # consider only 5'SS that were
                                # not counted before
                                current_5pSS = ":".join(
                                    [exon_coordinate_splited[0],
                                     exon_coordinate_splited[2],
                                     exon_coordinate_splited[2],
                                     exon_coordinate_splited[3]]
                                )

                                if current_5pSS in observed_5pSS:
                                    continue
                                else:
                                    observed_5pSS.append(current_5pSS)

                                # -------------------------------------------------
                                # CHECKS
                                # -------------------------------------------------

                                # 1. The 5' splice site is located at the 3'end
                                # of an upstream exon
                                # (+/- some tollerance = splice_fuzziness).

                                if (
                                    abs(int(exon_coordinate_splited[2]) -
                                        int(split_event.five_prime_ss)) <=
                                    splice_fuzziness
                                ):

                                    # 2. The 3' splice site is located withing
                                    #    the feature region.

                                    if (
                                        self.feat_region.start <=
                                        split_event.three_prime_ss and
                                        split_event.three_prime_ss <=
                                        self.feat_region.end
                                    ):

                                        if verbose:
                                            sys.stdout.write(
                                                "...3' splice site falls " +
                                                "within feature region." +
                                                os.linesep
                                            )

                                        # 3. The read ends within the feature
                                        #    region (if this is not the case,
                                        #    the read might go over another
                                        #    splice boarder and thus is not
                                        #    indicative for a terminal
                                        #    exon).
                                        if (
                                            detailed_alignment.aln.iv.end <=
                                            self.feat_region.end
                                        ):

                                            # 4. Create a dictionaty
                                            #    (possible_exonic_start_sites)
                                            #    that contains as key the
                                            #    possible start sites of the
                                            #    novel exon and the values is a
                                            #    new dictionary where the key
                                            #    is the upstream 5'SS of an
                                            #    upstream exon and the values
                                            #    is how many times this spliced
                                            #    read is observed.

                                            three_prime_ss_id = ":".join(
                                                [split_event.chrom,
                                                 str(split_event.three_prime_ss),
                                                 str(split_event.three_prime_ss),
                                                 split_event.strand]
                                            )

                                            self.possible_exonic_start_sites[
                                                three_prime_ss_id][
                                                    current_5pSS] += 1

                                            self.splice_into_feature_count += 1

                        elif detailed_alignment.aln.iv.strand == "-":

                            # Go over all upstream exon coordinates and check
                            # if the split reads starts in one of these exons
                            if verbose:
                                sys.stdout.write("... loop over all possible" +
                                                 "upstream exons ..." +
                                                 os.linesep)

                            observed_5pSS = []
                            # loop over all possible 5'SSs
                            for exon_coordinate in self.potential_5pSS_exons:

                                exon_coordinate_splited = \
                                    exon_coordinate.split(":")

                                # consider only 5'SS that are were
                                # not counted before
                                current_5pSS = ":".join(
                                    [exon_coordinate_splited[0],
                                     exon_coordinate_splited[1],
                                     exon_coordinate_splited[1],
                                     exon_coordinate_splited[3]]
                                )

                                if current_5pSS in observed_5pSS:
                                    continue
                                else:
                                    observed_5pSS.append(current_5pSS)

                                # -------------------------------------------------
                                # CHECKS
                                # -------------------------------------------------

                                # 1. The 5' splice site is located at the 3'end
                                #    of an upstream exon
                                #    (+/- some tollerance = splice_fuzziness).

                                if (
                                    abs(int(exon_coordinate_splited[1]) -
                                        int(split_event.five_prime_ss)) <=
                                    splice_fuzziness
                                ):

                                    # 2. The 3' splice site is located withing
                                    #    the feature region.

                                    if (
                                        self.feat_region.start <=
                                        split_event.three_prime_ss and
                                        split_event.three_prime_ss <=
                                        self.feat_region.end
                                    ):

                                        if verbose:
                                            sys.stdout.write(
                                                "...3' splice site " +
                                                "falls within feature" +
                                                "region." + os.linesep
                                            )

                                        # 3. The read ends within the
                                        #    feature region (if this is
                                        #    not the case, the read
                                        #    might go over another
                                        #    splice boarder and thus
                                        #    is not indicative for a
                                        #    terminal exon).

                                        if (
                                            detailed_alignment.aln.iv.start >=
                                            self.feat_region.start
                                        ):
                                            # 4. Create a dictionaty
                                            #    (possible_exonic_start_sites)
                                            #    that contains as key the
                                            #    possible start sites of the
                                            #    novel exon and the values is
                                            #    a new dictionary where the
                                            #    key is the upstream 5'SS of
                                            #    an upstream exon and the
                                            #    values is how many times
                                            #    this spliced read is
                                            #    observed.

                                            three_prime_ss_id = ":".join(
                                                [split_event.chrom,
                                                 str(split_event.three_prime_ss),
                                                 str(split_event.three_prime_ss),
                                                 split_event.strand])

                                            self.possible_exonic_start_sites[
                                                three_prime_ss_id][
                                                    current_5pSS] += 1

                                            self.splice_into_feature_count += 1

                        else:
                            sys.stderr.write("ERROR: unknown strand in " +
                                             "alignment!" +
                                             os.linesep)
                            sys.exit(-1)

    # _________________________________________________________________________
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Method: count_detailed_alignment
    # -------------------------------------------------------------------------
    def count_detailed_alignment_annotated(
        self,
        detailed_alignment,
        verbose=False
    ):

        # Case: no spliced reads
        if len(detailed_alignment.split_event_list) == 0:

                if detailed_alignment.aln.iv.strand == '+':

                    # Find unspliced reads that entirely fall within the
                    # exon coordinates
                    if (
                        int(detailed_alignment.aln.iv.start) >=
                        int(self.iv.start) and
                        int(detailed_alignment.aln.iv.start) <=
                        int(self.iv.end) and
                        int(detailed_alignment.aln.iv.end) >=
                        int(self.iv.start) and
                        int(detailed_alignment.aln.iv.end) <=
                        int(self.iv.end)
                    ):

                        self.annotated_unspliced_feature += 1

                    # Find unspliced reads that crosss the 5' border
                    if (
                        int(detailed_alignment.aln.iv.start) <=
                        int(self.iv.start) and
                        int(detailed_alignment.aln.iv.end) >=
                        int(self.iv.start)
                    ):

                        self.annotated_unspliced_5pSS += 1

                    # Find unspliced reads that crosss the 5' border
                    if (
                        int(detailed_alignment.aln.iv.start) <=
                        int(self.iv.end) and
                        int(detailed_alignment.aln.iv.end) >=
                        int(self.iv.end)
                    ):

                        self.annotated_unspliced_3pSS += 1

                elif detailed_alignment.aln.iv.strand == '-':

                    # Find unspliced reads that entirely fall within
                    # the exon coordinates
                    if (
                        int(detailed_alignment.aln.iv.start) >=
                        int(self.iv.start) and
                        int(detailed_alignment.aln.iv.start) <=
                        int(self.iv.end) and
                        int(detailed_alignment.aln.iv.end) >=
                        int(self.iv.start) and
                        int(detailed_alignment.aln.iv.end) <=
                        int(self.iv.end)
                    ):

                        self.annotated_unspliced_feature += 1

                    # Find unspliced reads that crosss the 5' border
                    if (
                        int(detailed_alignment.aln.iv.start) <=
                        int(self.iv.end) and
                        int(detailed_alignment.aln.iv.end) >=
                        int(self.iv.end)
                    ):

                        self.annotated_unspliced_5pSS += 1

                    # Find unspliced reads that crosss the 3' border
                    if (
                        int(detailed_alignment.aln.iv.start) <
                        int(self.iv.start) and
                        int(detailed_alignment.aln.iv.end) >
                        int(self.iv.start)
                    ):

                        self.annotated_unspliced_3pSS += 1

        # Case spliced reads
        elif len(detailed_alignment.split_event_list) > 0:

            for split_event in detailed_alignment.split_event_list:

                # -----------------------------------------------------
                # PLUS STRAND.
                # -----------------------------------------------------
                if detailed_alignment.aln.iv.strand == "+":

                    # -------------------------------------------------
                    # Find splice out reads
                    # -------------------------------------------------

                    # The 5' splice site is in the exon region and
                    # and the 3' splice site is located somwhere
                    # downstream (Splice out case)

                    if (
                        int(split_event.five_prime_ss) >=
                        int(self.iv.start) and
                        int(split_event.five_prime_ss) <=
                        int(self.iv.end) and
                        int(split_event.three_prime_ss) >=
                        int(self.iv.end) and
                        int(split_event.three_prime_ss) >=
                        int(self.iv.start)
                    ):

                        self.annotated_splice_out_all += 1

                    # The 5' splice site is located at the border
                    # of the feature region
                    if int(split_event.five_prime_ss) == int(self.iv.end):

                        self.annotated_splice_out_borders += 1

                    # -------------------------------------------------
                    # Find splice in reads
                    # -------------------------------------------------

                    # The 3' splice site is located withing the feature region.

                    if (
                        int(self.iv.start) <=
                        int(split_event.three_prime_ss) and
                        int(split_event.three_prime_ss) <=
                        int(self.iv.end)
                    ):

                        self.annotated_splice_in_all += 1

                    # The 3' splice site is located at the border
                    # of the feature region
                    if int(self.iv.start) == int(split_event.three_prime_ss):

                        self.annotated_splice_in_borders += 1

                elif detailed_alignment.aln.iv.strand == "-":

                    # -------------------------------------------------
                    # Find splice out reads
                    # -------------------------------------------------

                    # Same as above but for the minus strand

                    if (
                        int(split_event.five_prime_ss) >=
                        int(self.iv.start) and
                        int(split_event.five_prime_ss) <=
                        int(self.iv.end) and
                        int(split_event.three_prime_ss) <=
                        int(self.iv.start) and
                        int(split_event.three_prime_ss) <=
                        int(self.iv.end)
                    ):

                        self.annotated_splice_out_all += 1

                    if int(split_event.five_prime_ss) == int(self.iv.start):

                        self.annotated_splice_out_borders += 1

                    # -------------------------------------------------
                    # Find splice in reads
                    # -------------------------------------------------
                    if (
                        int(self.iv.start) <=
                        int(split_event.three_prime_ss) and
                        int(split_event.three_prime_ss) <=
                        int(self.iv.end)
                    ):

                        self.annotated_splice_in_all += 1

                    if int(split_event.three_prime_ss) == int(self.iv.end):

                        self.annotated_splice_in_borders += 1
