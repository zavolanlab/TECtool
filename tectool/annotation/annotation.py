# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import our own modules
# -----------------------------------------------------------------------------

from gene_structure import Exon
from gene_structure import Transcript
from gene_structure import Gene

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import pandas as pd
import HTSeq
import pybedtools
import os
import sys
from collections import defaultdict
from itertools import chain
flatten = chain.from_iterable
import itertools

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# CLASSES
# -----------------------------------------------------------------------------


class Annotation(object):
    """
    Handling annotation class.

        :param tmp: string that contains the path to the directory in which
            intermediate files will be written
        :rtype: Annotation object

    *Class members*

        *annotation_id*
            String. The id of the annotation.

        *exons*
            Dictionary. key: exon id, value: Exon object

        *transcripts*
            Dictionary. key: transcript id, value: Transcript object

        *genes*
            Dictionary. key: gene_id, value, Gene object

        *polyasites_in_introns*
            None pybedtool (bed) file that contains intronic polyasites

        *introns*
            None file path to bed file with intronic coordinates

        *feature_regions*
            List of feature regions

        *feature_regions_genes*
            Gene names/ids for feature regions

        *feature_regions_upstream_coordinates*
            Dictionary. key: Feature region value: list of upstream exons.

    """

    def __init__(self, annotation_id, tmp):

        self.annotation_id = annotation_id  # the annotation id
        self.tmp = tmp  # path to the directory to which file will be written

        self.exons = dict()  # key: exon id, value: Exon object
        # TODO: dictionary where key is exon id and values is a list of Exons.
        self.transcripts = dict()  # key: trx id, value: Transcript object
        self.genes = dict()  # key: gene_id, value, Gene object

        self.polyasites_in_introns = None  # pybedtool with intronic polyasites

        self.introns = None  # file path to bed file with intronic coordinates

        self.feature_regions = []  # List of feature regions
        self.feature_regions_genes = []  # Gene names/ids for feature regions
        # key: Feature region value: list of upstream exons.
        self.feature_regions_upstream_coordinates = dict()

    def extend(self, annotation, verbose=False):
        """
        Function that extends the annotation by another annotation.
        """

        if verbose:
            sys.stdout.write("Extending annotation '%s' by '%s'...%s" %
                             (self.annotation_id,
                              annotation.annotation_id,
                              os.linesep))

        # check for every gene in the annotation whether it exists already
        for gene_id in annotation.genes:

            if gene_id in self.genes:
                # if the gene exists already, extend it by the transcripts
                self.genes[gene_id].extend(annotation.genes[gene_id],
                                           verbose=verbose)
            else:
                # add the gene to the annotation
                if verbose:
                    sys.stdout.write(" :: adding gene '%s' (%s)%s" %
                                     (gene_id,
                                      annotation.genes[gene_id].gene_name,
                                      os.linesep))
                self.genes[gene_id] = annotation.genes[gene_id]

    def parse(
        self,
        gtf_file_path,
        verbose=False
    ):
        """
        Function that parses an annotation file (GTF format) and
        creates three dictionaries (genes, transcripts, exons).
        The keys of the dictionaries are the corresponding ids
         (gene_id, transcript_id, exon_id).

        *Examples for valid formated lines in the used GTF file:*

            *gene*
                3
                protein_coding
                gene
                107762145
                107809872
                .
                -
                .
                gene_id "ENSG00000196776";
                gene_name "CD47";
                gene_source "ensembl_havana";
                gene_biotype "protein_coding";

                | chromosome: 3
                | source:     protein_coding
                | feature:    gene
                | start:      107762145
                | end:        107809872
                | score:      .
                | strand:     -
                | frame:      .
                | attribute:
                |    gene_id "ENSG00000196776";
                |    gene_name "CD47";
                |    gene_source "ensembl_havana";
                |    gene_biotype "protein_coding";

            *transcript*
                3
                protein_coding
                transcript
                107762146
                107809872
                .
                -
                .
                gene_id "ENSG00000196776";
                transcript_id "ENST00000355354";
                gene_name "CD47";
                gene_source "ensembl_havana";
                gene_biotype "protein_coding";
                transcript_name "CD47-002";
                transcript_source "ensembl_havana";
                tag "CCDS";
                ccds_id "CCDS43125";

                | chromosome: 3
                | source:     protein_coding
                | feature:    transcript
                | start:      107762145
                | end:        107809872
                | score:      .
                | strand:     -
                | frame:      .
                | attribute:
                |    gene_id "ENSG00000196776";
                |    transcript_id "ENST00000355354";
                |    gene_name "CD47";
                |    gene_source "ensembl_havana";
                |    gene_biotype "protein_coding";
                |    transcript_name "CD47-002";
                |    transcript_source "ensembl_havana";
                |    tag "CCDS";
                |    ccds_id "CCDS43125";

            *exon*
                3
                protein_coding
                exon
                107809710
                107809872
                .
                -
                .
                gene_id "ENSG00000196776";
                transcript_id "ENST00000355354";
                exon_number "1";
                gene_name "CD47";
                gene_source "ensembl_havana";
                gene_biotype "protein_coding";
                transcript_name "CD47-002";
                transcript_source "ensembl_havana";
                tag "CCDS";
                ccds_id "CCDS43125";
                exon_id "ENSE00001813207";

                | chromosome: 3
                | source:     protein_coding
                | feature:    exon
                | start:      107809710
                | end:        107809872
                | score:      .
                | strand:     -
                | frame:      .
                | attribute:
                |    gene_id "ENSG00000196776";
                |    transcript_id "ENST00000355354";
                |    exon_number "1";
                |    gene_name "CD47";
                |    gene_source "ensembl_havana";
                |    gene_biotype "protein_coding";
                |    transcript_name "CD47-002";
                |    transcript_source "ensembl_havana";
                |    tag "CCDS";
                |    ccds_id "CCDS43125";
                |    exon_id "ENSE00001813207";

        :param gtf_file_path: string that contains the path to the gtf file

        """

        # get the file handle on the GTF/GFF file
        gtf_file = HTSeq.GFF_Reader(gtf_file_path)

        # loop over the file line by line
        for gtf_line in gtf_file:

            # Case gene
            if gtf_line.type == 'gene':

                gene_biotype = ''
                if "gene_biotype" in gtf_line.attr:
                    gene_biotype = gtf_line.attr['gene_biotype']
                elif "gene_type" in gtf_line.attr:
                    gene_biotype = gtf_line.attr['gene_type']

                # Create a Gene object
                gene = Gene(chromosome=gtf_line.iv.chrom,
                            source=gtf_line.source,
                            feature=gtf_line.type,
                            start=gtf_line.iv.start,
                            end=gtf_line.iv.end,
                            score=gtf_line.score,
                            strand=gtf_line.iv.strand,
                            frame=gtf_line.frame,
                            gene_id=gtf_line.attr['gene_id'],
                            gene_biotype=gene_biotype,
                            gene_name=gtf_line.attr['gene_name'])

                # Add the gene to the genes dict
                self.genes[gene.gene_id] = gene

            # Case transcript
            if gtf_line.type == 'transcript':

                gene_biotype = ''
                if "gene_biotype" in gtf_line.attr:
                    gene_biotype = gtf_line.attr['gene_biotype']
                elif "gene_type" in gtf_line.attr:
                    gene_biotype = gtf_line.attr['gene_type']

                transcript_biotype = ''
                if "transcript_biotype" in gtf_line.attr:
                    transcript_biotype = gtf_line.attr['transcript_biotype']
                elif "transcript_type" in gtf_line.attr:
                    transcript_biotype = gtf_line.attr['transcript_type']

                # Create a Transcript object
                transcript = Transcript(
                    chromosome=gtf_line.iv.chrom,
                    source=gtf_line.source,
                    feature=gtf_line.type,
                    start=gtf_line.iv.start,
                    end=gtf_line.iv.end,
                    score=gtf_line.score,
                    strand=gtf_line.iv.strand,
                    frame=gtf_line.frame,
                    gene_id=gtf_line.attr['gene_id'],
                    transcript_id=gtf_line.attr['transcript_id'],
                    gene_name=gtf_line.attr['gene_name'],
                    gene_biotype=gene_biotype,
                    transcript_name=gtf_line.attr['transcript_name'],
                    transcript_biotype=transcript_biotype
                )

                # add Transcript to the dictionary of transcripts
                self.transcripts[gtf_line.attr['transcript_id']] = transcript

            # Case exon
            if gtf_line.type == 'exon':

                gene_biotype = ''
                if "gene_biotype" in gtf_line.attr:
                    gene_biotype = gtf_line.attr['gene_biotype']
                elif "gene_type" in gtf_line.attr:
                    gene_biotype = gtf_line.attr['gene_type']

                transcript_biotype = ''
                if "transcript_biotype" in gtf_line.attr:
                    transcript_biotype = gtf_line.attr['transcript_biotype']
                elif "transcript_type" in gtf_line.attr:
                    transcript_biotype = gtf_line.attr['transcript_type']

                # Create an Exon object
                exon = Exon(chromosome=gtf_line.iv.chrom,
                            source=gtf_line.source,
                            feature=gtf_line.type,
                            start=gtf_line.iv.start,
                            end=gtf_line.iv.end,
                            score=gtf_line.score,
                            strand=gtf_line.iv.strand,
                            frame=gtf_line.frame,
                            gene_id=gtf_line.attr['gene_id'],
                            transcript_id=gtf_line.attr['transcript_id'],
                            exon_number=gtf_line.attr['exon_number'],
                            gene_name=gtf_line.attr['gene_name'],
                            gene_biotype=gene_biotype,
                            transcript_name=gtf_line.attr['transcript_name'],
                            transcript_biotype=transcript_biotype,
                            exon_id=gtf_line.attr['exon_id'])

                # Add exon to the dictionay of exons
                # (One representative exon is added in this one.)
                if gtf_line.attr['exon_id'] not in self.exons:
                    self.exons[gtf_line.attr['exon_id']] = exon

                # Add Exon to the sorted exon list of the corresponding
                # transcript
                if gtf_line.attr['transcript_id'] in self.transcripts:
                    self.transcripts[
                        gtf_line.attr['transcript_id']
                    ].insert_into_exon_list_sorted_by_end_coord(exon)
                else:
                    sys.stderr.write("ERROR: Transcript ID was not found "
                                     "for exon: " + gtf_line.attr['exon_id'])
                    sys.exit(-1)

                # Add Transcript to the transcripts dictionary of the gene
                if (gtf_line.attr['transcript_id'
                                  ] not in self.genes[transcript.gene_id
                                                      ].transcripts):

                    self.genes[gtf_line.attr['gene_id']
                               ].transcripts[gtf_line.attr['transcript_id']
                                             ] = transcript

                # Add Exon to the sorted exon list of the corresponding gene
                if gtf_line.attr['gene_id'] in self.genes:
                    self.genes[gtf_line.attr['gene_id']
                               ].insert_exon_into_exon_coordinates_dict(exon)
                else:
                    sys.stderr.write("ERROR: Gene ID was not found for exon" +
                                     gtf_line.attr['exon_id'])
                    sys.exit(-1)

            # Case start codon/CDS/stop codon
            if (
                gtf_line.type == 'start_codon' or
                gtf_line.type == 'CDS' or
                gtf_line.type == 'stop_codon'
            ):

                # check that the transcript is in the known
                # gene transcripts list
                if gtf_line.attr['transcript_id'
                                 ] in self.genes[gtf_line.attr['gene_id']
                                                 ].get_known_transctipt_ids():

                    # for all the transcripts (transcript ids)...
                    for transcript_id in self.genes[gtf_line.attr['gene_id']
                                                    ].transcripts:

                        # find the one that we are interested in
                        if transcript_id == gtf_line.attr['transcript_id']:

                            # ... and for all of the exons
                            for ex in self.genes[
                                gtf_line.attr['gene_id']].transcripts[
                                    transcript_id].exon_list_sorted_by_end_coord:

                                # ... find the one that has
                                # the specific exon number
                                if (
                                    int(ex.exon_number) ==
                                        int(gtf_line.attr['exon_number'])
                                ):

                                    # update start_codon
                                    if gtf_line.type == 'start_codon':
                                        ex.start_codon = gtf_line.iv

                                    if gtf_line.type == 'CDS':
                                        ex.CDS = gtf_line.iv
                                        ex.frame = gtf_line.frame

                                    if gtf_line.type == 'stop_codon':
                                        ex.stop_codon = gtf_line.iv

                else:
                    sys.stderr.write("ERROR: transcript id was not found the" +
                                     " list of transcripts in the Gene object")
                    sys.exit(-1)

    def write_exons_bed(self, bed):
        """
        Create bed file with exons
        Write them to a bed file where the column name contains:
        exon id $ gene id
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        all_exons = []
        for exon_id in self.exons:
            current_exon_as_bed = \
                self.exons[exon_id].get_bed_line_with_geneid_as_list()
            all_exons.append(current_exon_as_bed)
        all_exons_df = pd.DataFrame(all_exons)
        all_exons_df.columns = columns
        all_exons_df[['start', 'end']] = \
            all_exons_df[['start', 'end']].apply(pd.to_numeric)
        # sort dataframe based on chromosome and start position
        all_exons_df.sort_values(["chrom", "start"], inplace=True)

        all_exons_df.to_csv(bed,
                            sep="\t",
                            index=False,
                            header=False)

        return

    def write_transcripts_bed(self, bed):
        """
        Create bed file with transcripts
        Write them to a bed file where the column name contains gene id
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        all_transcripts = []
        for transcript_id in self.transcripts:
            current_transcript_as_bed = \
                self.transcripts[transcript_id].get_bed_with_geneid_as_list()
            all_transcripts.append(current_transcript_as_bed)
        all_transcripts_df = pd.DataFrame(all_transcripts)
        all_transcripts_df.columns = columns
        all_transcripts_df[['start', 'end']] = \
            all_transcripts_df[['start', 'end']].apply(pd.to_numeric)
        # sort dataframe based on chromosome and start position
        all_transcripts_df.sort_values(["chrom", "start"],
                                       ascending=[True, True],
                                       inplace=True)

        all_transcripts_df.to_csv(bed,
                                  sep="\t",
                                  index=False,
                                  header=False)

        return

    def write_genes_bed(self, bed):
        """
        Create bed file with genes
        Write them to a bed file
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        all_genes = []
        for gene_id in self.genes:
            current_gene_as_bed = \
                self.genes[gene_id].get_bed_with_geneid_as_list()
            all_genes.append(current_gene_as_bed)

        all_genes_df = pd.DataFrame(all_genes)
        all_genes_df.columns = columns
        all_genes_df[['start', 'end']] = \
            all_genes_df[['start', 'end']].apply(pd.to_numeric)
        # sort dataframe based on chromosome and start position
        all_genes_df.sort_values(["chrom", "start"],
                                 ascending=[True, True],
                                 inplace=True)

        all_genes_df.to_csv(bed,
                            sep="\t",
                            index=False,
                            header=False)

        return

    def write_union_exon_bed(self, bed):
        """
        Create a union exon bed file of all the genes
        If a gene overlaps with another gene, only the exons of the gene
        of interest are considered.
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        all_union_exons = []

        for gene_id in self.genes:
            current_gene_union_exons = self.genes[gene_id].get_union_exon_bed()
            all_union_exons += current_gene_union_exons
        # create dataframe
        all_union_exons_df = pd.DataFrame(all_union_exons)
        all_union_exons_df.columns = columns
        all_union_exons_df[["start", "end"]] = \
            all_union_exons_df[["start", "end"]].apply(pd.to_numeric)
        # sort dataframe based on chromosome and start position
        all_union_exons_df.sort_values(["chrom", "start"],
                                       ascending=[True, True],
                                       inplace=True)

        all_union_exons_df.to_csv(bed,
                                  sep="\t",
                                  index=False,
                                  header=False)

        return

    def write_union_intron_bed(self,
                               union_exons_bed,
                               union_introns_bed):
        """
        Create a union intron bed file of all the genes
        If a gene overlaps with another gene, only the introns of the gene
        of interest are considered.
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        columns_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "name": "object",
            "score": "object",
            "strand": "object"
        }

        # read union exons in a pandas dataframe
        union_exons_df = pd.read_csv(union_exons_bed,
                                     header=None,
                                     sep="\t",
                                     names=columns,
                                     dtype=columns_dtype)

        all_union_introns = []

        for gene_id, gene_df in union_exons_df.groupby("name"):
            union_exons = gene_df.loc[:, ["start", "end"]].values.tolist()

            union_introns = \
                self.genes[gene_id].generate_union_introns_per_gene(
                    union_exons=union_exons)

            all_union_introns += union_introns

        # free mem
        del union_exons_df

        all_union_introns_df = pd.DataFrame(all_union_introns)
        all_union_introns_df.columns = columns

        all_union_introns_df[["start", "end"]] = \
            all_union_introns_df[["start", "end"]].apply(pd.to_numeric)
        # sort dataframe based on chromosome and start position
        all_union_introns_df.sort_values(["chrom", "start"],
                                         ascending=[True, True],
                                         inplace=True)

        all_union_introns_df.to_csv(union_introns_bed,
                                    sep="\t",
                                    index=False,
                                    header=False)

    def write_introns_that_do_not_overlap_with_exons_stranded_bed(
        self,
        union_exons_bed,
        union_introns_bed,
        selected_introns_bed
    ):
        """
        Find introns that do not overlap with exons of other
        genes and write them in bed format.
        Note: Function assumes that the data are stranded.
        """

        union_exons = pybedtools.BedTool(union_exons_bed)
        union_introns = pybedtools.BedTool(union_introns_bed)

        union_introns.intersect(
            union_exons,
            s=True,
            wa=True,
            v=True).sort().saveas(selected_introns_bed)

        return

    def write_introns_that_do_not_overlap_with_exons_unstranded_bed(
        self,
        union_exons_bed,
        union_introns_bed,
        selected_introns_bed
    ):
        """
        Find introns that do not overlap with exons of other
        genes and write them in bed format.
        Note: Function assumes that the data are stranded.
        """

        union_exons = pybedtools.BedTool(union_exons_bed)
        union_introns = pybedtools.BedTool(union_introns_bed)

        union_introns.intersect(
            union_exons,
            s=False,
            wa=True,
            v=True).sort().saveas(selected_introns_bed)

        return

    def write_intronic_polya_sites(
        self,
        polyasites,
        introns,
        drop_intronic_polya_sites_of_overlapping_genes,
        intronic_polya_sites_bed
    ):
        """
        Write intronic polya sites
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        columns_custom_bed = [
            "chrom",
            "start",
            "end",
            "name_bad",
            "score",
            "strand",
            "chrom2",
            "start2",
            "end2",
            "name",
            "score2",
            "strand2"]

        columns_custom_bed_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "name_bad": "object",
            "score": "object",
            "strand": "object",
            "chrom2": "object",
            "start2": "int64",
            "end2": "int64",
            "name": "object",
            "score2": "object",
            "strand2": "object"
        }

        # read bedtools objects
        polyasites_bed = pybedtools.BedTool(polyasites)
        introns_bed = pybedtools.BedTool(introns)

        # tmp file
        tmp_intronic_polya_sites_bed = intronic_polya_sites_bed + "tmp.bed"

        # find all intronic polya sites
        # note:
        # At this point it can contain duplicated entries,
        # polya sites fount introns of overlapping genes
        polyasites_bed.intersect(
            introns_bed,
            s=True,
            wa=True,
            wb=True,
            f=1
        ).sort().saveas(tmp_intronic_polya_sites_bed)

        # read them as dataframes
        intronic_polya_sites_df = pd.read_csv(
            tmp_intronic_polya_sites_bed,
            sep="\t",
            header=None,
            names=columns_custom_bed,
            dtype=columns_custom_bed_dtype
        )

        intronic_polya_sites_df = intronic_polya_sites_df[columns]

        if drop_intronic_polya_sites_of_overlapping_genes:

            # remove all duplicates
            intronic_polya_sites_df.drop_duplicates(
                subset=["chrom", "start", "end", "strand"],
                keep=False,
                inplace=True
            )

        # write file
        intronic_polya_sites_df.to_csv(
            intronic_polya_sites_bed,
            sep="\t",
            header=False,
            index=False
        )

        # and clean up tmp file
        os.remove(tmp_intronic_polya_sites_bed)

        return

    def write_terminal_exons(self,
                             bed):
        """
        Finds terminal exons of the transcripts and write them to a bed file.
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        terminal_exons = []

        for transcript in self.transcripts:

            # make sure that the transcripts has more than 1 exon
            if (
                len(self.transcripts[transcript].exon_list_sorted_by_end_coord
                    ) > 1
            ):

                # determine terminal exon
                terminal_exon = \
                    self.transcripts[transcript
                                     ].get_bed_terminal_exon_as_list()

                terminal_exons.append(terminal_exon)

        terminal_exons_df = pd.DataFrame(terminal_exons)
        terminal_exons_df.columns = columns

        terminal_exons_df[["start", "end"]] = \
            terminal_exons_df[["start", "end"]].apply(pd.to_numeric)

        terminal_exons_df.sort_values(["chrom", "start"],
                                      ascending=[True, True],
                                      inplace=True)

        terminal_exons_df.drop_duplicates(
            subset=["chrom", "start", "end", "name", "strand"],
            keep="first",
            inplace=True
        )

        terminal_exons_df.to_csv(
            bed,
            sep="\t",
            header=False,
            index=False
        )

        return

    def write_intermediate_exons(self,
                                 bed):
        """
        Finds intermediate exons of the transcripts and write them to
        a bed file.
        """
        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        intermediate_exons = []

        for transcript in self.transcripts:

            # make sure that the transcripts has more than 1 exon
            if (
                len(self.transcripts[transcript].exon_list_sorted_by_end_coord
                    ) > 1
            ):

                # determine intermediate exons
                cutrrent_intermediate_exons = \
                    self.transcripts[transcript
                                     ].get_bed_intermediate_exons_as_list()

                intermediate_exons += cutrrent_intermediate_exons

        intermediate_exons_df = pd.DataFrame(intermediate_exons)
        intermediate_exons_df.columns = columns

        intermediate_exons_df[["start", "end"]] = \
            intermediate_exons_df[["start", "end"]].apply(pd.to_numeric)

        intermediate_exons_df.sort_values(["chrom", "start"],
                                          ascending=[True, True],
                                          inplace=True)

        intermediate_exons_df.drop_duplicates(
            subset=["chrom", "start", "end", "name", "strand"],
            keep="first",
            inplace=True
        )

        intermediate_exons_df.to_csv(
            bed,
            sep="\t",
            header=False,
            index=False
        )

        return

    def write_start_exons(self,
                          bed):
        """
        Finds intermediate exons of the transcripts and write them to
        a bed file.
        """
        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        start_exons = []

        for transcript in self.transcripts:

            # make sure that the transcripts has more than 1 exon
            if (
                len(self.transcripts[transcript].exon_list_sorted_by_end_coord
                    ) > 1
            ):

                # determine intermediate exons
                cutrrent_start_exons = \
                    self.transcripts[transcript
                                     ].get_bed_start_exon_as_list()

                start_exons.append(cutrrent_start_exons)

        start_exons_df = pd.DataFrame(start_exons)
        start_exons_df.columns = columns

        start_exons_df[["start", "end"]] = \
            start_exons_df[["start", "end"]].apply(pd.to_numeric)

        start_exons_df.sort_values(["chrom", "start"],
                                   ascending=[True, True],
                                   inplace=True)

        start_exons_df.drop_duplicates(
            subset=["chrom", "start", "end", "name", "strand"],
            keep="first",
            inplace=True
        )

        start_exons_df.to_csv(
            bed,
            sep="\t",
            header=False,
            index=False
        )

        return

    def write_non_overlapping_regions_to_bed(
        self,
        bed_in,
        bed_out
    ):
        """
        Extract region that overlap only one time form a bed file
        """

        # custom header
        columns_custom_bed = [
            "chrom",
            "start",
            "end",
            "strand",
            "count",
            "name"
        ]

        columns_custom_bed_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "strand": "object",
            "count": "int64",
            "name": "object"
        }

        # normal header
        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        # tmp file name
        tmp_file = bed_out + "_tmp"

        # find overlapping regions within the file
        pybedtools.BedTool(bed_in).merge(
            s=True,
            c=4,
            o="count,collapse"
        ).sort().saveas(tmp_file)

        # read file as a dataframe
        df = pd.read_csv(tmp_file,
                         sep="\t",
                         header=None,
                         names=columns_custom_bed,
                         dtype=columns_custom_bed_dtype)

        selected_df = df[df["count"] == 1]

        # add score columns and get only the required columns
        selected_df = \
            selected_df.assign(score="0")[columns]

        selected_df[["start", "end"]] = \
            selected_df[["start", "end"]].apply(pd.to_numeric)

        selected_df.sort_values(["chrom", "start"],
                                ascending=[True, True],
                                inplace=True)

        selected_df.to_csv(
            bed_out,
            sep="\t",
            header=False,
            index=False
        )

        os.remove(tmp_file)

        return

    def write_one_copy_for_each_region(self,
                                       bed_in,
                                       bed_out):
        """
        Given a bed file extract one copy for each of the regions
        Example: Keep one copy for each exon (since one exon
        might appear to multiple transcripts)
        """
        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        columns_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "name": "object",
            "score": "object",
            "strand": "object"
        }

        regions_df = pd.read_csv(bed_in,
                                 sep="\t",
                                 header=None,
                                 names=columns,
                                 dtype=columns_dtype)

        regions_df[["start", "end"]] = \
            regions_df[["start", "end"]].apply(pd.to_numeric)

        regions_df.drop_duplicates(
            subset=["chrom", "start", "end", "name", "strand"],
            keep="first",
            inplace=True
        )

        regions_df.to_csv(
            bed_out,
            sep="\t",
            header=False,
            index=False
        )

        return

    def write_non_overlapping_regions(
        self,
        selected_bed,
        compare_bed_1,
        compare_bed_2,
        strand,
        bed_out
    ):
        """
        Take a bed file and a list of bed files and
        find regions from A (selected_bed) that
        are not in B (bed_list)
        """

        selected_bed = pybedtools.BedTool(selected_bed)
        compare_bed_1 = pybedtools.BedTool(compare_bed_1)
        compare_bed_2 = pybedtools.BedTool(compare_bed_2)

        selected_bed.intersect(
            compare_bed_1,
            s=strand,
            v=True,
            wa=True
        ).intersect(
            compare_bed_2,
            s=strand,
            v=True,
            wa=True
        ).saveas(bed_out)

        return

    def remove_regions_that_overlap_in_the_opposite_strand(
        self,
        bed_in,
        bed_out
    ):
        """
        Function that takes as input a bed file and removes
        regions that with the same chromosome, start and 
        end positions.
        """


        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        columns_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "name": "object",
            "score": "object",
            "strand": "object"
        }

        df = pd.read_csv(bed_in,
                         sep="\t",
                         header=None,
                         names=columns,
                         dtype=columns_dtype)

        df.sort_values(["chrom", "start"],
                       ascending=[True, True],
                       inplace=True)

        df.drop_duplicates(
            subset=["chrom", "start", "end"],
            keep=False,
            inplace=True
        )

        df.to_csv(
            bed_out,
            sep="\t",
            header=False,
            index=False
        )

        return

    def write_regions_identified_only_once_from_two_beds(
        self,
        selected_regions_bed,
        all_regions_bed,
        non_ovelapping_bed,
        strand
    ):
        """
        Take two bed files and find regions from
        A (selected_regions_bed) that are not
        in B (all_regions_bed) more than once
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        columns_custom_bed = ["chrom",
                              "start",
                              "end",
                              "name",
                              "score",
                              "strand",
                              "count"]

        columns_custom_bed_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "name": "object",
            "score": "object",
            "strand": "object",
            "count": "int64"
        }

        tmp_file = non_ovelapping_bed + "_tmp"

        selected_regions = pybedtools.BedTool(selected_regions_bed)

        all_regions = pybedtools.BedTool(all_regions_bed)

        selected_regions.intersect(
            all_regions,
            s=strand,
            wa=True,
            c=True
        ).sort().saveas(tmp_file)

        non_ovelapping_bed_df = pd.read_csv(tmp_file,
                                            sep="\t",
                                            header=None,
                                            names=columns_custom_bed,
                                            dtype=columns_custom_bed_dtype)

        non_ovelapping_bed_df[["start", "end"]] = \
            non_ovelapping_bed_df[["start", "end"]].apply(pd.to_numeric)

        non_ovelapping_bed_df[non_ovelapping_bed_df[
            "count"] == 1][columns].to_csv(non_ovelapping_bed,
                                           sep="\t",
                                           header=False,
                                           index=False)

        os.remove(tmp_file)

        return

    def write_non_overlapping_genes_to_bed(
        self,
        all_gene_regions,
        strand,
        non_overlapping_genes_bed
    ):
        """
        Function that takes all gene coordinates bed file and writes a genes
        coordinates bed file that contains only the non overlapping genes.
        It considers also whether the protocol is stranded or not.
        """

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        columns_custom_bed = ["chrom",
                              "start",
                              "end",
                              "name",
                              "score",
                              "strand",
                              "count"]

        columns_custom_bed_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "name": "object",
            "score": "object",
            "strand": "object",
            "count": "int64"
        }

        tmp_file = non_overlapping_genes_bed + "_tmp"

        # load gene coordinates
        all_genes = pybedtools.BedTool(all_gene_regions)

        all_genes.intersect(all_genes,
                            s=strand,
                            c=True).sort().saveas(tmp_file)

        non_overlapping_genes_df = pd.read_csv(tmp_file,
                                               sep="\t",
                                               header=None,
                                               names=columns_custom_bed,
                                               dtype=columns_custom_bed_dtype)

        non_overlapping_genes_df[non_overlapping_genes_df["count"] == 1][
            columns].to_csv(
                non_overlapping_genes_bed,
                sep="\t",
                header=False,
                index=False
        )

        os.remove(tmp_file)

        return

    def determine_union_exon_length_per_gene(
        self,
        union_exons_bed
    ):
        """
        Find the length of the union exon for each gene
        and adds it in the gene annotation class.
        """

        union_exon_lengths = dict()

        union_exons = pybedtools.BedTool(union_exons_bed)

        for exon in union_exons:

            current_exon_length = int(exon.end) - int(exon.start)
            gene_id = str(exon.name)

            if gene_id not in union_exon_lengths:
                union_exon_lengths[gene_id] = current_exon_length
            else:
                union_exon_lengths[gene_id] += current_exon_length

        # Add the union exon length in the self.genes dictionary
        for gene_id in union_exon_lengths:

            self.genes[gene_id].union_exon_length = \
                int(union_exon_lengths[gene_id])

        return

    def determine_feature_regions(
        self,
        polyasites
    ):
        """
        Determine feature regions and upstream exonic coordinates.
        Creates a list for the features and a dictionary where key
        is the feature and the value a list of the upstream exons.
        Example:
        key: 22:29670442:29670929:+
        value: ['22:29669729:29669853:+',
                '22:29664304:29664338:+',
                '22:29670253:29670442:+',
                '22:29664279:29664338:+', ...]
        """

        polyasites_bed = pybedtools.BedTool(polyasites)

        # for site in self.polyasites_in_introns:
        for site in polyasites_bed:

            gene_id = site.name

            minimum_distance = 1000000000

            upstream_exons_coordinates = []

            feature_chr = ''
            feature_start = ''
            feature_end = ''
            feature_strand = ''
            feature_id = ''

            if site.strand == '+':

                for exon_coordinates in self.genes[gene_id].exon_coordinates_dict:

                    # select one representative exon based on the coordinates.
                    current_exon = \
                        self.exons[self.genes[
                            gene_id].exon_coordinates_dict[
                                exon_coordinates][0]]

                    diff = int(site.start) - int(current_exon.end)

                    if diff >= 0:
                        if diff < minimum_distance:

                            minimum_distance = diff

                            feature_chr = current_exon.chromosome
                            feature_start = current_exon.end
                            feature_end = site.start
                            feature_strand = '+'
                            feature_id = ":".join([feature_chr,
                                                   str(feature_start),
                                                   str(feature_end + 1),
                                                   feature_strand,
                                                   gene_id])

                        upstream_exons_coordinates.append(exon_coordinates)

                if feature_id is not '':
                    self.feature_regions.append(feature_id)
                    # might be possible to remove now
                    self.feature_regions_genes.append(gene_id)
                    self.feature_regions_upstream_coordinates[feature_id] \
                        = upstream_exons_coordinates

            elif site.strand == '-':

                for exon_coordinates in self.genes[gene_id].exon_coordinates_dict:

                    current_exon = self.exons[
                        self.genes[
                            gene_id].exon_coordinates_dict[
                                exon_coordinates][0]]

                    diff = int(current_exon.start) - int(site.end)

                    if diff >= 0:
                        if diff < minimum_distance:
                            minimum_distance = diff

                            feature_chr = current_exon.chromosome
                            feature_start = site.end
                            feature_end = current_exon.start
                            feature_strand = '-'
                            feature_id = ":".join([feature_chr,
                                                   str(feature_start - 1),
                                                   str(feature_end),
                                                   feature_strand,
                                                   gene_id])

                        upstream_exons_coordinates.append(exon_coordinates)

                if feature_id is not '':
                    self.feature_regions.append(feature_id)
                    self.feature_regions_genes.append(gene_id)
                    self.feature_regions_upstream_coordinates[feature_id] \
                        = upstream_exons_coordinates

            else:
                sys.stderr.write('[ERROR] No strand information available.')
                sys.exit(-1)

    def get_genes_with_multiexonic_transcripts(self, bed):
        """
        Function that reads a bed file with union exons of genes
        and returns a dictionary of genes with multiple exons per gene"
        """

        df_selected_dict = dict()

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        columns_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "name": "object",
            "score": "object",
            "strand": "object"
        }

        # read them as dataframes
        df = pd.read_csv(
            bed,
            sep="\t",
            header=None,
            names=columns,
            dtype=columns_dtype
        )

        df_grouped = df.groupby("name").count()["score"]

        df_selected_list = list(df_grouped[df_grouped > 1].index)

        for i in df_selected_list:
            df_selected_dict[i] = i

        return df_selected_dict

    def create_dictionary_of_non_overlapping_genes(
        self,
        non_overlapping_genes_bed
    ):
        """
        Function that takes a a bed file of non overlapping genes
        and returens a dictionary that contains as key the gene_id
        and value the gene_id. It uses the information stored in
        the column name of the bed file.
        """

        non_overlapping_genes_dict = dict()

        columns = ["chrom",
                   "start",
                   "end",
                   "name",
                   "score",
                   "strand"]

        columns_dtype = {
            "chrom": "object",
            "start": "int64",
            "end": "int64",
            "name": "object",
            "score": "object",
            "strand": "object"
        }

        # read them as dataframes
        df = pd.read_csv(
            non_overlapping_genes_bed,
            sep="\t",
            header=None,
            names=columns,
            dtype=columns_dtype
        )

        df_selected_list = list(set(df["name"].tolist()))

        for i in df_selected_list:
            non_overlapping_genes_dict[i] = i

        return non_overlapping_genes_dict

    def filter_terminal_exon_training_candidates(self):
        """
        Remove terminal exons that the 3p crossing is not 0
        and the sum of the last 100 bases of the profile is not > 0 reads
        """

        for gene_id in self.genes:

            # copy list
            list_of_terminal_exons = list(
                self.genes[gene_id].annotated_terminal_exons)

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
                if (
                    int(terminal_exon.unspliced_3pSS) != 0 or
                    sum_profile_end == 0
                ):

                    self.genes[gene_id].annotated_terminal_exons.remove(
                        terminal_exon)

    def determine_CDS_novel_transcrtips_from_spicing(
        self,
        genome
    ):
        """
        Determine CDS sequence and stop codon for novel transcripts
        """

        stop_codons = ['TAA', 'TAG', 'TGA']

        # go over all the genes
        for gene in self.genes:

            # and all the novel transcritps
            for transcript in self.genes[gene].novel_transcripts:

                start_codon_observed = False
                stop_codon_observed = False
                start_codon_observed = \
                    transcript.contains_start_codon_for_novel_transcripts()
                stop_codon_observed = \
                    transcript.contains_stop_codon_for_novel_transcripts()

                # Start and stop codon were observed. This means that the novel
                # exon falls in the 3'UTR. The existing (annotated) CDS will
                # be used without any changes.
                if start_codon_observed and stop_codon_observed:

                    transcript.write_CDS = True

                # Only transcsripts that have a start codon and not a stop
                # codon are considered here. These are the ones that might
                # have a novel CDS/stop.
                elif start_codon_observed and not stop_codon_observed:

                    # get the novel exons of the transcript
                    novel_exons = transcript.get_novel_exons_sorted_by_end()

                    # take the last exon
                    last_exon = novel_exons[-1] if transcript.strand == '+' else novel_exons[0]

                    # and the penultimate exon
                    penultimate_exon = novel_exons[-2] if transcript.strand == '+' else novel_exons[1]

                    # if the frame from the one but last exon is missing
                    # do not write the CDS for the novel transcript
                    if penultimate_exon.frame is None:
                        transcript.write_CDS = False
                        continue

                    else:

                        # find the cds sequence of the penultimate exon
                        # what we do is to change the start of the CDS based
                        # one the frame (CDS.start + frame)
                        penultimate_CDS_seq = genome.sequence(
                            {
                                'chr': penultimate_exon.chromosome,
                                'start': int(penultimate_exon.CDS.start) +
                                int(penultimate_exon.frame),
                                'stop': penultimate_exon.CDS.end,
                                'strand': penultimate_exon.strand
                            },
                            one_based=False
                        )

                        last_exon_seq = genome.sequence(
                            {
                                'chr': last_exon.chromosome,
                                'start': int(last_exon.start),
                                'stop': int(last_exon.end),
                                'strand': last_exon.strand
                            },
                            one_based=False
                        )

                        # Use the modulo to check if we get 0,1, or 2
                        # (number of bases left)
                        bases_left = len(penultimate_CDS_seq) % 3

                        possible_last_exon_frame = None

                        # determine frame for novel exon
                        if bases_left == 0:

                            possible_last_exon_frame = 0
                            # keep the bases (sequence) left from
                            # the previous exon
                            previous_bases = ''

                        elif bases_left == 1:

                            possible_last_exon_frame = 2
                            # keep the bases (sequence) left from
                            # the previous exon
                            previous_bases = penultimate_CDS_seq[-bases_left:]

                        elif bases_left == 2:

                            possible_last_exon_frame = 1
                            # keep the bases (sequence) left from
                            # the previous exon
                            previous_bases = penultimate_CDS_seq[-bases_left:]

                        # merge the bases left from the previous exon/CDS
                        # with the last exon
                        possible_cds = previous_bases + last_exon_seq

                        # Position of stop codon that will be observed
                        first_stop_codon_position = None

                        # Search for a stop codon
                        for i in range(0, len(possible_cds), 3):

                            # print(i)
                            # print(possible_cds[i:i+3])

                            # stop codon was detected
                            if possible_cds[i:i + 3].upper() in stop_codons:
                                first_stop_codon_position = i
                                break

                        # A stop codon was observed. The CDS for the novel
                        # transcript is calculated.
                        if first_stop_codon_position is not None:

                            if bases_left > first_stop_codon_position:
                                transcript.write_CDS = False
                                sys.stderr.write("-" * 80 + os.linesep)
                                sys.stderr.write("[WARNING] Manually check " +
                                                 "this exon for translation." +
                                                 os.linesep)
                                sys.stderr.write("-" * 80 + os.linesep)
                                sys.stderr.write("chromosome: " +
                                                 str(last_exon.chromosome) +
                                                 os.linesep)
                                sys.stderr.write("start: " +
                                                 str(last_exon.start) +
                                                 os.linesep)
                                sys.stderr.write("end: " +
                                                 str(last_exon.end) +
                                                 os.linesep)
                                sys.stderr.write("strand: " +
                                                 str(last_exon.strand) +
                                                 os.linesep)
                                sys.stderr.write(
                                    "first stop codon position: " +
                                    str(first_stop_codon_position) +
                                    os.linesep
                                )
                                sys.stderr.write("bases left: " +
                                                 str(bases_left) +
                                                 os.linesep)
                                sys.stderr.write(
                                    "possible last exon frame: " +
                                    str(possible_last_exon_frame) +
                                    os.linesep
                                )

                                sys.stderr.write("-" * 80 +
                                                 os.linesep +
                                                 os.linesep)
                                continue

                            if last_exon.strand == "+":

                                # CDS for the novel exon
                                CDS = HTSeq.GenomicInterval(
                                    last_exon.chromosome,
                                    last_exon.start,
                                    last_exon.start +
                                    first_stop_codon_position -
                                    bases_left,
                                    last_exon.strand
                                )

                                # stop codon for the novel exon
                                stop_codon = HTSeq.GenomicInterval(
                                    last_exon.chromosome,
                                    last_exon.start +
                                    first_stop_codon_position -
                                    bases_left,
                                    last_exon.start +
                                    first_stop_codon_position -
                                    bases_left + 2,
                                    last_exon.strand
                                )

                                last_exon.CDS = CDS
                                last_exon.frame = possible_last_exon_frame
                                last_exon.stop_codon = stop_codon

                                transcript.write_CDS = True

                            elif last_exon.strand == "-":

                                # CDS for the novel exon
                                CDS = HTSeq.GenomicInterval(
                                    last_exon.chromosome,
                                    last_exon.end -
                                    first_stop_codon_position +
                                    bases_left,
                                    last_exon.end,
                                    last_exon.strand
                                )

                                # stop codon for the novel exon
                                stop_codon = HTSeq.GenomicInterval(
                                    last_exon.chromosome,
                                    last_exon.end -
                                    first_stop_codon_position +
                                    bases_left - 2,
                                    last_exon.end -
                                    first_stop_codon_position +
                                    bases_left,
                                    last_exon.strand
                                )

                                # update information for the last exon
                                last_exon.CDS = CDS
                                last_exon.frame = possible_last_exon_frame
                                last_exon.stop_codon = stop_codon

                                transcript.write_CDS = True

                        else:

                            transcript.write_CDS = False

                # No stop codon was found, so the CDS should not be
                # written for this transcript
                else:
                    transcript.write_CDS = False

    def write_gtf(
        self,
        output_file,
        with_CDS,
        accepted_exons_dict="write_all"
    ):
        """
        Write out GFT file
        """

        for gene in self.genes:

            # write gene
            self.genes[gene].write_gtf(output_file)

            # Find which novel transcripts to use
            novel_transcripts_fitlered = \
                self.genes[gene].decide_novel_transcripts_to_use()

            novel_and_known_merged_sorted = \
                sorted(novel_transcripts_fitlered +
                       self.genes[gene].get_known_transcripts(),
                       key=lambda x: x.start, reverse=False)

            if accepted_exons_dict == "write_all":

                for transcript in novel_and_known_merged_sorted:

                    # write transcript
                    transcript.write_transcript_gtf(output_file)

                    # write exons
                    transcript.write_exons_gtf(output_file, with_CDS)

            elif isinstance(accepted_exons_dict, dict):

                for transcript in novel_and_known_merged_sorted:

                    if "TECtool" in transcript.source:

                        for exon in transcript.novel_exons:

                            exon_id = ":".join([
                                exon.chromosome,
                                str(exon.start),
                                str(exon.end),
                                exon.strand])

                            if exon_id in accepted_exons_dict:

                                # write transcript
                                transcript.write_transcript_gtf(output_file)

                                # write exons
                                transcript.write_exons_gtf(output_file,
                                                           with_CDS)

                                break
                    else:
                        # write transcript
                        transcript.write_transcript_gtf(output_file)

                        # write exons
                        transcript.write_exons_gtf(output_file, with_CDS)

    def write_transcript2novel_transcript_mapping_list(
        self,
        output_file_handle,
        accepted_exons_dict="write_all"
    ):
        """
        Function that takes as input a list of accepted novel exons
        and writes out a list of that contains transcripts ids  and
        novel transcipt ids. The transcript ids are all the potential
        "mother" transcripts that generate the novel transcript ids.
        """

        # in case of write all then we should write all novel transcripts
        if accepted_exons_dict == "write_all":
            pass
        elif isinstance(accepted_exons_dict, dict):

            # go over the genes
            for gene in self.genes:

                # Find which novel transcripts to use
                novel_transcripts_fitlered = \
                    self.genes[gene].decide_novel_transcripts_to_use()

                for novel_transcript in novel_transcripts_fitlered:

                    if "TECtool" in novel_transcript.source:

                        for exon in novel_transcript.novel_exons:

                            exon_id = ":".join([
                                exon.chromosome,
                                str(exon.start),
                                str(exon.end),
                                exon.strand]
                            )

                            if exon_id in accepted_exons_dict:

                                novel_transcript_id = \
                                    novel_transcript.transcript_id

                                for transcript_id in self.genes[gene].mother_transcripts_of_novel_transcripts[novel_transcript_id]:

                                    output_file_handle.write(
                                        "{} \t {} {}".format(
                                            transcript_id,
                                            novel_transcript_id,
                                            os.linesep)
                                    )

                                break

    def write_transcript2gene_mapping_list(
        self,
        transcript2gene_mapping_file_path
    ):
        """
        Creates a file that contains transcript to gene mappings for
        all transcripts that exist in the Annotation object.
        """

        # check whether self.transcripts exists.
        if not self.transcripts:
            sys.stderr.write(
                "ERROR: 'Annotation.transcripts' is empty! \
                 Thus, writing {} is not possible! {}".format(
                    transcript2gene_mapping_file_path,
                    os.linesep)
            )
            sys.exit(-1)

        # open the file for writing
        mapping_file_handle = open(transcript2gene_mapping_file_path, 'w')
        file_header = "transcript_id\tgene_id\tgene_name" + os.linesep
        mapping_file_handle.write(file_header)

        # write out mappings for all transcripts that
        # exist in the Annotation object.
        for transcript_id in self.transcripts:

            # create the line to write
            line_to_write = \
                self.transcripts[transcript_id].transcript_id + \
                "\t" + \
                self.transcripts[transcript_id].gene_id + \
                "\t" + \
                self.transcripts[transcript_id].gene_name + \
                os.linesep

            # write the line
            mapping_file_handle.write(line_to_write)

        # close the file
        mapping_file_handle.close()

    def get_all_terminal_exons(self):
        """
        Return list of all terminal exons.
        It can contain duplicate terminal exons.
        """

        list_of_terminal_exons = []

        # go over all genes
        for gene_id in self.genes:

            # and all transcripts
            for transcript in self.genes[gene_id].get_known_transcripts():

                # get terminal exon for each
                list_of_terminal_exons.append(transcript.get_terminal_exon())

        return(list_of_terminal_exons)

    def get_all_first_and_intermediate_exons(self):
        """
        Return list of all first and intermediate exons.
        So all exons except the terminal ones.
        It can contain duplicate exons.
        """

        list_of_first_and_intermediate_exons = []

        # go over all genes
        for gene_id in self.genes:

            # and all transcripts
            for transcript in self.genes[gene_id].get_known_transcripts():

                exons = transcript.get_all_exons_but_last()

                for exon in exons:

                    list_of_first_and_intermediate_exons.append(exon)

        return(list_of_first_and_intermediate_exons)

    def write_all_terminal_exons_as_bed(self, outfile):
        """
        Writes a file with all terminal exons as bed.
        It can contain duplicate terminal exons.
        """

        all_terminal_exons = self.get_all_terminal_exons()

        fh = open(outfile, 'w')
        for terminal_exon in all_terminal_exons:
            fh.write("\t".join([terminal_exon.chromosome,
                                str(terminal_exon.start),
                                str(terminal_exon.end),
                                terminal_exon.gene_id,
                                ".",
                                terminal_exon.strand + os.linesep]))
        fh.close()

    def write_all_first_and_intermediate_exons_as_bed(self, outfile):
        """
        Writes a file with all first and intermediate exons.
        It can contain duplicate terminal exons.
        """

        all_first_and_intermediate_exons = \
            self.get_all_first_and_intermediate_exons()

        fh = open(outfile, 'w')
        for exon in all_first_and_intermediate_exons:
            fh.write("\t".join([exon.chromosome,
                                str(exon.start),
                                str(exon.end),
                                exon.gene_id,
                                ".",
                                exon.strand + os.linesep]))
        fh.close()
