# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import HTSeq
import pybedtools
import os
import sys
from collections import defaultdict
from collections import Counter
from itertools import chain
from itertools import combinations
flatten = chain.from_iterable
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import random

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import our own modules
# -----------------------------------------------------------------------------

from gene_structure import Exon
from gene_structure import Transcript
from gene_structure import Gene

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

        self.annotation_id = annotation_id # the annotation id
        self.tmp = tmp                     # path to the directory to which file will be written
        
        self.exons = {}         # Dictionary. key: exon id, value: Exon object # TODO: dictionary where key is exon id and values is a list of Exons.
        self.transcripts = {}   # Dictionary. key: transcript id, value: Transcript object
        self.genes = {}         # Dictionary. key: gene_id, value, Gene object

        self.polyasites_in_introns = None  # pybedtool (bed) file that contains intronic polyasites

        self.introns = None # file path to bed file with intronic coordinates

        self.feature_regions = []   # List of feature regions
        self.feature_regions_genes = [] # Gene names/ids for feature regions
        self.feature_regions_upstream_coordinates = {} # Dictionary. key: Feature region value: list of upstream exons.


    def extend(self, annotation, verbose=False):

        """
        Function that extends the annotation by another annotation.
        """

        if verbose:
            sys.stdout.write("Extending annotation '%s' by annotation '%s'...\n" % (self.annotation_id, annotation.annotation_id))

        # check for every gene in the annotation whether it exists already
        for gene_id in annotation.genes:

            if gene_id in self.genes:
                # if the gene exists already, extend it by the transcripts (if necessary)
                self.genes[gene_id].extend(annotation.genes[gene_id], 
                                           verbose=verbose)
            
            else:
                # add the gene to the annotation
                if verbose:
                    sys.stdout.write(" :: adding gene '%s' (%s)\n" \
                                     % (gene_id, 
                                        annotation.genes[gene_id].gene_name))
                self.genes[gene_id] = annotation.genes[gene_id]


    def parse(self, gtf_file_path, verbose=False):

        """

        Function that parses an annotation file (GTF format) and creates three dictionaries (genes, transcripts, exons).
        The keys of the dictionaries are the corresponding ids (gene_id, transcript_id, exon_id).

        *Examples for valid formated lines in the used GTF file:*

            *gene*
                3   protein_coding  gene    107762145   107809872   .   -   .   gene_id "ENSG00000196776"; gene_name "CD47"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
                
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
                3   protein_coding  transcript  107762146   107809872   .   -   .   gene_id "ENSG00000196776"; transcript_id "ENST00000355354"; gene_name "CD47"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CD47-002"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS43125";

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
                3   protein_coding  exon    107809710   107809872   .   -   .   gene_id "ENSG00000196776"; transcript_id "ENST00000355354"; exon_number "1"; gene_name "CD47"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CD47-002"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS43125"; exon_id "ENSE00001813207";

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
                transcript = Transcript(chromosome=gtf_line.iv.chrom, 
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
                                        transcript_biotype=transcript_biotype)
                
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

                # Add exon to the dictionay of exons (One representative exon is added in this one.)
                if gtf_line.attr['exon_id'] not in self.exons:
                    self.exons[gtf_line.attr['exon_id']] = exon

                # Add Exon to the sorted exon list of the corresponding transcript
                if gtf_line.attr['transcript_id'] in self.transcripts:
                    self.transcripts[gtf_line.attr['transcript_id']].insert_into_exon_list_sorted_by_end_coord(exon)
                else:
                    sys.stderr.write("ERROR: Transcript ID was not found for exon: " + gtf_line.attr['exon_id'])
                    sys.exit(-1)

                # Add Transcript to the transcripts dictionary of the gene
                if gtf_line.attr['transcript_id'] not in self.genes[transcript.gene_id].transcripts:
                    self.genes[gtf_line.attr['gene_id']].transcripts[gtf_line.attr['transcript_id']] = transcript

                # Add Exon to the sorted exon list of the corresponding gene
                if gtf_line.attr['gene_id'] in self.genes:
                    self.genes[gtf_line.attr['gene_id']].insert_exon_into_exon_coordinates_dict(exon)
                else:
                   sys.stderr.write("ERROR: Gene ID was not found for exon " + gtf_line.attr['exon_id'])
                   sys.exit(-1)

            # Case start codon/CDS/stop codon
            if gtf_line.type == 'start_codon' or gtf_line.type == 'CDS' or gtf_line.type == 'stop_codon':

                # check that the transcript is in the known gene transcripts list
                if gtf_line.attr['transcript_id'] in (self.genes[gtf_line.attr['gene_id']].get_known_transctipt_ids()):

                    # for all the transcripts (transcript ids)...
                    for transcript_id in self.genes[gtf_line.attr['gene_id']].transcripts:

                        # find the one that we are interested in
                        if transcript_id == gtf_line.attr['transcript_id']:

                            # ... and for all of the exons
                            for ex in self.genes[gtf_line.attr['gene_id']].transcripts[transcript_id].exon_list_sorted_by_end_coord:

                                # ... find the one that has the specific exon number 
                                if int(ex.exon_number) == int(gtf_line.attr['exon_number']):

                                    # update start_codon
                                    if gtf_line.type == 'start_codon':
                                        ex.start_codon = gtf_line.iv

                                    if gtf_line.type == 'CDS':
                                        ex.CDS = gtf_line.iv
                                        ex.frame = gtf_line.frame

                                    if gtf_line.type == 'stop_codon':
                                        ex.stop_codon = gtf_line.iv

                else:
                    sys.stderr.write("ERROR: transcript id was not found the list of transcripts in the Gene object")
                    sys.exit(-1)


    def determine_intronic_sites(self,
                                polya_regions,
                                sequencing_direction,
                                exons_per_gene_bed_file,
                                transcripts_per_gene_bed_file,
                                genes_bed_file,
                                polya_sites_in_transcripts_bed_file,
                                polya_sites_in_introns_bed_file,
                                union_exons_bed_file,
                                union_introns_bed_file,
                                introns_that_do_not_overlap_with_exons_stranded_bed_file,
                                introns_that_do_not_overlap_with_exons_unstranded_bed_file,
                                polyasites_in_introns_stranded_bed_file,
                                polyasites_in_introns_unstranded_bed_file):

        """
        Function that takes as input the polya sites and finds intronic poly(A) sites.
        Returns bedtools object containing the intronic polya sites.
        The file is different depending on whether you work on stranded or 
        unstranded data. More poly(A) sites will be examined when a stranded 
        protocol is used.
        """

        # _____________________________________________________________________
        # STEP 1
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Create bedtools objects for the exons
        # ---------------------------------------------------------------------

        # create tmp file
        tmp_bed_file = exons_per_gene_bed_file+'.tmp'

        # Loop over self.exons dictionary and write exons in bed format
        w = open(tmp_bed_file,'w')
        for exon in self.exons:

            w.write("\t".join([str(self.exons[exon].chromosome), 
                                str(self.exons[exon].start), 
                                str(self.exons[exon].end), 
                                self.exons[exon].exon_id+'$'+self.exons[exon].gene_id, 
                                '.', 
                                self.exons[exon].strand+"\n"]))
        w.close()

        # sort bed file
        os.system("sort -k1,1 -k2,2n %s > %s" % (tmp_bed_file, exons_per_gene_bed_file))

        # remove tmp file
        os.system("rm %s" % (os.path.join(tmp_bed_file)))

        # Read exons bed file using pybedtools
        exons_bed = pybedtools.BedTool(exons_per_gene_bed_file)

        # _____________________________________________________________________
        # STEP 2
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Create bedtools objects for the transcripts
        # ---------------------------------------------------------------------

        # create tmp file
        tmp_bed_file = transcripts_per_gene_bed_file+'.tmp'

        # Loop over self.exons dictionary and write exons in bed format
        w = open(tmp_bed_file,'w')
        for transcript_id in self.transcripts:

            w.write("\t".join([str(self.transcripts[transcript_id].chromosome), 
                                str(self.transcripts[transcript_id].start), 
                                str(self.transcripts[transcript_id].end), 
                                self.transcripts[transcript_id].gene_id, 
                                '.', 
                                self.transcripts[transcript_id].strand+"\n"]))
        w.close()

        # sort bed file
        os.system("sort -k1,1 -k2,2n %s > %s" % (tmp_bed_file, transcripts_per_gene_bed_file))

        # remove tmp file
        os.system("rm %s" % (os.path.join(tmp_bed_file)))

        # Read exons bed file using pybedtools
        transcripts_bed = pybedtools.BedTool(transcripts_per_gene_bed_file)

        # _____________________________________________________________________
        # STEP 3
        #
        # ---------------------------------------------------------------------
        # Create bedtools objects for the genes
        # ---------------------------------------------------------------------

        # create tmp file
        tmp_bed_file = genes_bed_file+'.tmp'

        # Loop over self.genes dictionary and write genes in bed format
        w = open(tmp_bed_file,'w')

        for gene in self.genes:

            w.write(self.genes[gene].get_actual_gene_coordinates_bed()+"\n")

        w.close()

        # sort bed file
        os.system("sort -k1,1 -k2,2n %s > %s" % (tmp_bed_file, genes_bed_file))
        
        # remove tmp file
        os.system("rm %s" % (tmp_bed_file))

        # Read genes bed file using pybedtools
        genes_bed = pybedtools.BedTool(genes_bed_file)

        # _____________________________________________________________________
        # STEP 4
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Find poly(A) sites that fall within transcripts
        # (remove intergenic poly(A) sites) and store it in a file.
        # ---------------------------------------------------------------------

        # Find polya sites that fall within transcript coordinates
        polyasites_in_transcripts = polya_regions.intersect(transcripts_bed, s=True, wa=True, f=1, wb=True)

        # Reformat file
        w = open(polya_sites_in_transcripts_bed_file+"_tmp",'w')
        for site in polyasites_in_transcripts:
            tmp_sites = site.fields
            tmp_sites[3] = ":".join([tmp_sites[3], tmp_sites[9]])
            tmp_sites[5] = tmp_sites[5]+"\n"
            w.write("\t".join(tmp_sites[0:6]))
        w.close()

        # sort bed file
        os.system("sort -u %s | sort -k1,1 -k2,2n > %s" % (polya_sites_in_transcripts_bed_file+"_tmp", polya_sites_in_transcripts_bed_file))

        # remove tmp file
        os.system("rm %s" % (polya_sites_in_transcripts_bed_file+"_tmp"))

        # Find polyasites that fall within the gene loci
        polyasites_in_transcripts_enriched = pybedtools.BedTool(polya_sites_in_transcripts_bed_file)

        # _____________________________________________________________________
        # STEP 5
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Find poly(A) sites that fall within introns (remove exonic and 
        # antisense poly(A) sites) and store it in a file.
        # ---------------------------------------------------------------------

        # Keep polyasites that fall within introns
        polyasites_in_introns = polyasites_in_transcripts_enriched.intersect(exons_bed, s=True, wa=True, v=True)
        polyasites_in_introns.saveas(polya_sites_in_introns_bed_file)

        # _____________________________________________________________________
        # STEP 6
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Generate union exon bed file
        # ---------------------------------------------------------------------

        self.determine_union_exon_bed(union_exons_bed_file)

        # _____________________________________________________________________
        # STEP 7
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Generate union exon bed file
        # ---------------------------------------------------------------------

        self.determine_union_intron_bed(union_introns_bed_file)
        union_introns_bed = pybedtools.BedTool(union_introns_bed_file)

        # _____________________________________________________________________
        # STEP 8
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Find introns that do not overlap with other exons (stranded)
        # ---------------------------------------------------------------------

        introns_that_do_not_overlap_with_exons_stranded = union_introns_bed.intersect(union_exons_bed_file, s=True, wa=True, v=True)
        introns_that_do_not_overlap_with_exons_stranded.saveas(introns_that_do_not_overlap_with_exons_stranded_bed_file)

        # _____________________________________________________________________
        # STEP 9
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Find introns that do not overlap with other exons (unstranded)
        # ---------------------------------------------------------------------

        introns_that_do_not_overlap_with_exons_unstranded = union_introns_bed.intersect(union_exons_bed_file, s=False, wa=True, v=True)
        introns_that_do_not_overlap_with_exons_unstranded.saveas(introns_that_do_not_overlap_with_exons_unstranded_bed_file)

        # _____________________________________________________________________
        # STEP 10
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Find intronic poly(A) sites for stranded protocol
        # ---------------------------------------------------------------------

        polyasites_in_introns_stranded = polyasites_in_introns.intersect(introns_that_do_not_overlap_with_exons_stranded, s=True, wa=True, f=1)
        polyasites_in_introns_stranded.saveas(polyasites_in_introns_stranded_bed_file+"_tmp")

        os.system("sort -u %s > %s" % (polyasites_in_introns_stranded_bed_file+"_tmp",polyasites_in_introns_stranded_bed_file))

        os.system("rm %s" % (polyasites_in_introns_stranded_bed_file+"_tmp"))

        # _____________________________________________________________________
        # STEP 11
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Find intronic poly(A) sites for unstranded protocol
        # ---------------------------------------------------------------------
        polyasites_in_introns_unstranded = polyasites_in_introns.intersect(introns_that_do_not_overlap_with_exons_unstranded, s=True, wa=True, f=1)
        polyasites_in_introns_unstranded.saveas(polyasites_in_introns_unstranded_bed_file+"_tmp")

        os.system("sort -u %s > %s" % (polyasites_in_introns_unstranded_bed_file+"_tmp", polyasites_in_introns_unstranded_bed_file))

        os.system("rm %s" % (polyasites_in_introns_unstranded_bed_file+"_tmp"))

        # -----------------------------------------------------------------------------
        # Decide according to the data that you have (first, second, unstranded)
        # -----------------------------------------------------------------------------

        # Choose intronic poly(A) sites according to the strand
        if "forward" in sequencing_direction:
            self.polyasites_in_introns = pybedtools.BedTool(polyasites_in_introns_stranded_bed_file)
            self.introns = pybedtools.BedTool(introns_that_do_not_overlap_with_exons_stranded_bed_file)
        elif "reverse" in sequencing_direction:
            self.polyasites_in_introns = pybedtools.BedTool(polyasites_in_introns_stranded_bed_file)
            self.introns = pybedtools.BedTool(introns_that_do_not_overlap_with_exons_stranded_bed_file)
        elif "unstranded" in sequencing_direction:
            self.polyasites_in_introns = pybedtools.BedTool(polyasites_in_introns_unstranded_bed_file)
            self.introns = pybedtools.BedTool(introns_that_do_not_overlap_with_exons_stranded_bed_file)
        else:
            sys.stderr.write("[ERROR]: Strand information strand was not specified properly")


    def determine_intergenic_sites(self, polya_regions, genes_bed_file, intergenic_polyasites_stranded_bed_file, intergenic_polyasites_unstranded_bed_file):

        """Determine intergenic sites"""

        # Read genes bed file using pybedtools
        genes_bed = pybedtools.BedTool(genes_bed_file)

        # find intergenic polya sites stranded
        intergenic_polyasites_stranded = polya_regions.intersect(genes_bed, s=True, wa=True, v=True)

        intergenic_polyasites_stranded.saveas(intergenic_polyasites_stranded_bed_file)

        # find intergenic polya sites unstranded
        intergenic_polyasites_unstranded = intergenic_polyasites_stranded.intersect(genes_bed_file, S=True, wa=True, v=True)

        intergenic_polyasites_unstranded.saveas(intergenic_polyasites_unstranded_bed_file)


    def determine_intronic_coordinates(self, polya_regions, background_regions_bed_file_path, feature_length_threshold):

        """ Determine intronic coordinates """

        if self.introns is not None:

            # reads introns
            introns = pybedtools.BedTool(self.introns)

            # use bedtools to keep only the introns that do NOT contain any poly(A) site
            introns_without_polya_sites = introns.intersect(polya_regions, s=False, v=True)

            # create dictionary of the introns that do not contain poly(A) sites
            introns_without_polya_sites_dict = dict()
            for intron in introns_without_polya_sites:

                coord = HTSeq.GenomicInterval(str(intron.chrom), int(intron.start), int(intron.end), str(intron.strand))

                if coord not in introns_without_polya_sites_dict:
                    introns_without_polya_sites_dict[coord] = coord


            # write output bed file
            w = open(background_regions_bed_file_path, 'w')
            for intron in introns:

                coord = HTSeq.GenomicInterval(str(intron.chrom), int(intron.start), int(intron.end), str(intron.strand))

                chrom = ''
                if 'chr' in intron.chrom:
                    chrom = intron.chrom
                else:
                    chrom = 'chr'+intron.chrom

                # flag to keep long enough introns
                if int(intron.end) - int(intron.start) > int(feature_length_threshold):

                    intron_feature = 'intron_without_polya_site' if coord in introns_without_polya_sites_dict  else 'intron_with_polya_site'

                    w.write("\t".join([chrom, 
                            str(int(intron.start)+1), 
                            str(int(intron.end)-1),
                            str(intron.name)+","+intron_feature,
                            '.',
                            str(intron.strand)+"\n"
                            ]))
            w.close()

        else:
            sys.stderr.write("Error. self.introns was not initialized properly")
            sys.exit(-1)


    def determine_introns_without_polyAsite(self, polya_regions, intronic_regions_no_polya_bed_file_path):

        """
        Determine introns that do not contain poly(A) sites
        """

        if self.introns is not None:
            
            # reads introns
            introns = pybedtools.BedTool(self.introns)
    
            # use bedtools to keep only the introns that do NOT contain any poly(A) site
            introns_without_polya_sites = introns.intersect(polya_regions, s=False, v=True)

            # write output bed file
            w = open(intronic_regions_no_polya_bed_file_path, 'w')
            for intron in introns_without_polya_sites:

                # flag to keep long enough introns
                if int(intron.end) - int(intron.start) > 50:

                    w.write("\t".join([str(intron.chrom), 
                            str(int(intron.start)+1), 
                            str(int(intron.end)-1),
                            str(intron.name)+",intron_without_polya_site",
                            '.',
                            str(intron.strand)+"\n"
                            ]))
            w.close()

        else:
            sys.stderr.write("Error. self.introns was not initialized properly")
            sys.exit(-1)


    def determine_introns_with_polyAsite(self, polya_regions, background_regions_bed_file_path, feature_length_threshold):

        """
        Determine introns that contain poly(A) sites
        """

        if self.introns is not None:
            
            # reads introns
            introns = pybedtools.BedTool(self.introns)
    
            # use bedtools to keep only the introns that poly(A) sites
            introns_with_polya_sites = introns.intersect(polya_regions, s=True, wa=True).saveas(background_regions_bed_file_path+'_tmp.bed')

            # 
            os.system('sort -u %s > %s' % (background_regions_bed_file_path+'_tmp.bed', background_regions_bed_file_path+'_tmp.2.bed'))

            introns_with_polya_sites = pybedtools.BedTool(background_regions_bed_file_path+'_tmp.2.bed')

            w = open(background_regions_bed_file_path, 'w')

            for intron in introns_with_polya_sites:

                # flag to keep long enough introns
                if int(intron.end) - int(intron.start) > int(feature_length_threshold):
                    w.write("\t".join([intron.chrom, 
                            str(int(intron.start)+1), 
                            str(int(intron.end)-1),
                            str(intron.name)+",intron_with_polya_site",
                            '.',
                            str(intron.strand)+"\n"]))

            w.close()

            os.system("rm %s" % (background_regions_bed_file_path+'_tmp.bed'))
            os.system("rm %s" % (background_regions_bed_file_path+'_tmp.2.bed'))

        else:
            sys.stderr.write("Error. self.introns was not initialized properly")
            sys.exit(-1)


    def determine_non_overlapping_genes(self, all_gene_regions, sequencing_direction, non_overlapping_genes_bed):

        """Function that takes all gene coordinates bed file and writes a genes coordinates bed file that
            contains only the non overlapping genes. It considers also whether the protocol is stranded or
            not. It also returns a dictionary with the specific non overlapping genes."""

        non_overlapping_genes = dict()

        # load gene coordinates
        all_genes = pybedtools.BedTool(all_gene_regions)

        # select strand option for bedtools based on protocol
        if sequencing_direction == "forward" or sequencing_direction == "reverse":
            strand = True
        elif sequencing_direction == "unstranded":
            strand = False

        # determine non overlapping genes
        w = open(non_overlapping_genes_bed, 'w')

        for gene in all_genes.intersect(all_genes, s=strand, c=True):

            gene_id = gene.fields[3]

            # if the gene overlaps with itself
            if int(gene.fields[6]) == 1:

                self.genes[gene_id].overlaps_with_other_gene = False
                # gene is not in the dictionary
                if gene_id not in non_overlapping_genes:

                    non_overlapping_genes[gene_id] = gene_id
                    # write bed file
                    w.write("\t".join([str(gene.fields[0]), 
                                       str(gene.fields[1]),
                                       str(gene.fields[2]),
                                       str(gene.fields[3]),
                                       str(gene.fields[4]),
                                       str(gene.fields[5])+"\n"]))
        w.close()


        
    def determine_non_overlapping_genes_old(self, all_gene_regions, sequencing_direction, non_overlapping_genes_bed):

        """Function that takes all gene coordinates bed file and writes a genes coordinates bed file that
            contains only the non overlapping genes. It considers also whether the protocol is stranded or
            not. It also returns a dictionary with the specific non overlapping genes."""

        non_overlapping_genes = dict()
        
        
        genes_collapsed = pybedtools.BedTool(all_gene_regions).merge(s=True, c=4, o="count,collapse")
        # elif sequencing_direction == "unstranded":
            # genes_collapsed = pybedtools.BedTool(all_gene_regions).merge(s=False, c=4, o="count,collapse")

        w = open(non_overlapping_genes_bed, 'w')
        
        for gene in genes_collapsed:

            chromosome = str(gene.fields[0])
            start = str(gene.fields[1])
            end = str(gene.fields[2])
            strand = str(gene.fields[3])

            number_of_overlaps = int(str(gene.fields[4]))
            gene_id = str(gene.fields[5])

            # gene does not overlap with any other gene
            if number_of_overlaps == 1:

                self.genes[gene_id].overlaps_with_other_gene = False

                # gene is not in the dictionary
                if gene_id not in non_overlapping_genes:

                    non_overlapping_genes[gene_id] = gene_id

                    # write bed file
                    w.write("\t".join([chromosome, start, end, gene_id, ".", strand+"\n"]))

        w.close()


    def determine_union_exon_bed(self, union_exons_bed):

        """Create a union exon bed file of all the genes
           If a gene overlaps with another gene, only the exons of the gene 
           of interest are considered.
        """

        w = open(union_exons_bed+"_tmp", 'w')

        for gene_id in self.genes:

            w.write(self.genes[gene_id].get_union_exon_bed())

        w.close()

        os.system("sort -k 1,1 -k2,2n %s > %s" % (union_exons_bed+"_tmp", union_exons_bed))

        os.system("rm %s" % (union_exons_bed+"_tmp"))


    def determine_union_intron_bed(self, union_introns_bed_file):

        """Create a union intron bed file of all the genes
           If a gene overlaps with another gene, only the introns of the gene 
           of interest are considered.
        """
        
        w = open(union_introns_bed_file+"_tmp", 'w')
        for gene_id in self.genes:
            for intron in self.genes[gene_id].generate_union_introns_per_gene():
                w.write(intron+"\n")
        w.close()

        os.system("sort -k 1,1 -k2,2n %s > %s" % (union_introns_bed_file+"_tmp", union_introns_bed_file))

        os.system("rm %s" % (union_introns_bed_file+"_tmp"))
        

    def determine_union_exon_bed_for_non_overlapping_genes(self, all_exon_coordinates, genes_to_consider, union_exons_bed):

        all_exons = pybedtools.BedTool(all_exon_coordinates)

        # create temporary bed file with the exons of the non overlapping genes

        tmp_exons_file = union_exons_bed+'tmp'

        w = open(tmp_exons_file, 'w')
        for exon in all_exons:

            gene_id = str(exon.name).split("$")[1]

            if gene_id in genes_to_consider:

                w.write("\t".join([str(exon.chrom), 
                                   str(exon.start),
                                   str(exon.end), 
                                   str(gene_id), 
                                   '.', 
                                   str(exon.strand)+"\n"]))
        w.close()

        # Load bed file and generate union exons
        union_exons = pybedtools.BedTool(tmp_exons_file).merge(s=True, c=4, o="count,collapse")


        w = open(union_exons_bed, 'w')
        for exon in union_exons:

            chromosome = str(exon.fields[0])
            start = str(exon.fields[1])
            end = str(exon.fields[2])
            strand = str(exon.fields[3])
            number_of_overlaps = int(str(exon.fields[4]))
            gene_id = str(exon.fields[5]).split(",")[0]

            w.write("\t".join([chromosome, start, end, gene_id, '.', strand+"\n"]))

        w.close()

        # remove tmp file
        os.system("rm %s" % (tmp_exons_file))


    def determine_union_exon_length_per_gene(self, union_exons_bed):

        """ Find the length of the union exon for each gene
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

            self.genes[gene_id].union_exon_length = int(union_exon_lengths[gene_id])


    def determine_feature_regions(self):

        """
        Determine feature regions and upstream exonic coordinates.
        Creates a list for the features and a dictionary where key is the feature and the value a list of the upstream exons.
        Example:
        key: 22:29670442:29670929:+
        value: ['22:29669729:29669853:+', '22:29664304:29664338:+', '22:29670253:29670442:+', '22:29664279:29664338:+', ...]
        """

        for site in self.polyasites_in_introns:

            gene_name = site.name.split(':')[-1]

            minimum_distance = 1000000000

            upstream_exon_coordinates = ''
            upstream_exons_coordinates = []

            feature_chr = ''
            feature_start = ''
            feature_end = ''
            feature_strand = ''
            feature_id = ''

            if site.strand == '+':

                for exon_coordinates in self.genes[gene_name].exon_coordinates_dict:

                    current_exon = self.exons[self.genes[gene_name].exon_coordinates_dict[exon_coordinates][0]] # select one representative exon based on the coordinates. WHY ???

                    diff = int(site.start) - int(current_exon.end)

                    if diff >= 0:
                        if diff < minimum_distance:

                            minimum_distance = diff
                            upstream_exon_coordinates = exon_coordinates

                            feature_chr = current_exon.chromosome
                            feature_start = current_exon.end                    
                            feature_end = site.start
                            feature_strand = '+'
                            feature_id = ":".join([feature_chr, str(feature_start), str(feature_end+1), feature_strand, gene_name]) # added gene name here

                        upstream_exons_coordinates.append(exon_coordinates)

                if feature_id is not '':
                    self.feature_regions.append(feature_id)
                    self.feature_regions_genes.append(gene_name)          # might be possible to remove now
                    self.feature_regions_upstream_coordinates[feature_id] = upstream_exons_coordinates

            elif site.strand == '-':

                for exon_coordinates in self.genes[gene_name].exon_coordinates_dict:

                    current_exon = self.exons[self.genes[gene_name].exon_coordinates_dict[exon_coordinates][0]]
                    
                    diff = int(current_exon.start) - int(site.end)

                    if diff >= 0:
                        if diff < minimum_distance:
                            minimum_distance = diff
                            upstream_exon_coordinates = exon_coordinates

                            feature_chr = current_exon.chromosome
                            feature_start = site.end
                            feature_end = current_exon.start
                            feature_strand = '-'
                            feature_id = ":".join([feature_chr, str(feature_start-1), str(feature_end), feature_strand, gene_name])

                        upstream_exons_coordinates.append(exon_coordinates)

                if feature_id is not '':
                    self.feature_regions.append(feature_id)
                    self.feature_regions_genes.append(gene_name)
                    self.feature_regions_upstream_coordinates[feature_id] = upstream_exons_coordinates

            else:
                sys.stderr.write('[ERROR] No strand information available.')
                sys.exit(-1)



    def determine_intermediate_exons(self, intermediate_exons_bed_file_path):

        """
           Function that finds intermediate exons of the transcripts
           and writes them to a file.
        """
        
        intermediate_exons_tmp_file_path = \
            intermediate_exons_bed_file_path + "_tmp"
        w_intermediate = open(intermediate_exons_tmp_file_path, 'w')

        for transcript in self.transcripts:

            # make sure that the transcripts has more than 1 exon
            if len(self.transcripts[transcript].exon_list_sorted_by_end_coord)>1:

                # determine intermediate exons
                annotated_intermediate_exons = \
                    self.transcripts[transcript].get_intermediate_exons()
    
                for annotated_intermediate_exon in annotated_intermediate_exons:
    
                    w_intermediate.write("\t".join([annotated_intermediate_exon.chromosome,\
                        str(annotated_intermediate_exon.start),\
                        str(annotated_intermediate_exon.end),\
                        str(annotated_intermediate_exon.gene_id),\
                        str("."),\
                        str(annotated_intermediate_exon.strand+"\n") ]))

        w_intermediate.close()

        # Keep only unique entries and sort by coordinates
        os.system("sort -u %s | sort -k1,1 -k2,2n > %s" % (intermediate_exons_tmp_file_path, intermediate_exons_bed_file_path))
        
        # Remove tmp file
        os.system("rm %s" % (intermediate_exons_tmp_file_path))

        # Read the bed file with pybedtools and merge it in order to identify which 
        # are the intermediate exons that overlap with other intermediate exons and
        # which not.
        # The command that is executed is:
        # bedtools merge -i annotated_intermediate_exons.bed -s -c 4 -o count,collapse
        intermediate_exons = pybedtools.BedTool(intermediate_exons_bed_file_path).merge(s=True, c=4, o="count,collapse")


        # Write the output file in bed format and in the name field
        # we mix the gene id together with whether it overlaps with
        # some other terminal exon or not
        w_intermediate = open(intermediate_exons_bed_file_path, 'w')
        
        for intermediate_exon in intermediate_exons:
            
            # find if overlaps with other terminal exons
            overlap = 'no_overlap_with_other_intermediate_exons' if int(str(intermediate_exon.fields[4])) == 1 else 'overlap_with_other_intermediate_exons'

            if overlap == "no_overlap_with_other_intermediate_exons":

                # get gene id
                gene_id = str(str(intermediate_exon.fields[5]).split(',')[0])
    
                # write bed file
                w_intermediate.write("\t".join([str(intermediate_exon.fields[0]),
                                    str(intermediate_exon.fields[1]),
                                    str(intermediate_exon.fields[2]),
                                    ",".join([gene_id, overlap]),
                                    '.',
                                    str(intermediate_exon.fields[3])+"\n"]))
        w_intermediate.close()


    def determine_terminal_exons(self, terminal_exons_bed_file_path):

        """
           Finds terminal exons of the transcripts
           and writes them to a bed file.
        """
        
        terminal_exons_tmp_file_path = \
            terminal_exons_bed_file_path + "_tmp"
        w_terminal = open(terminal_exons_tmp_file_path, 'w')

        for transcript in self.transcripts:

            # make sure that the transcripts has more than 1 exon
            if len(self.transcripts[transcript].exon_list_sorted_by_end_coord)>1:

                # determine terminal exon
                annotated_terminal_exon = self.transcripts[transcript].get_terminal_exon()
    
                w_terminal.write("\t".join([annotated_terminal_exon.chromosome,\
                 str(annotated_terminal_exon.start),\
                 str(annotated_terminal_exon.end),\
                 str(annotated_terminal_exon.gene_id),\
                 str("."),\
                 str(annotated_terminal_exon.strand+"\n") ]))

        w_terminal.close()


        # Keep only unique entries and sort by coordinates
        os.system("sort -u %s | sort -k1,1 -k2,2n > %s" % (terminal_exons_tmp_file_path, terminal_exons_bed_file_path))

        # Remove tmp file
        os.system("rm %s" % (terminal_exons_tmp_file_path))

        # Read the bed file with pybedtools and merge it in order to identify which 
        # are terminal exons overlap with other terminal exons and which not.
        #
        # bedtools merge -i annotated_terminal_exons.bed -s -c 4 -o count,collapse
        terminal_exons = pybedtools.BedTool(terminal_exons_bed_file_path).merge(s=True, c=4, o="count,collapse")

        # Write the output file in bed format and in the name field
        # we mix the gene id together with whether it overlaps with
        # some other terminal exon or not
        w_terminal = open(terminal_exons_bed_file_path, 'w')
        
        for terminal_exon in terminal_exons:
            
            # find if overlaps with other terminal exons
            overlap = 'no_overlap_with_other_terminal_exons' if int(str(terminal_exon.fields[4])) == 1 else 'overlap_with_other_terminal_exons'

            if overlap == 'no_overlap_with_other_terminal_exons':

                # get gene id
                gene_id = str(str(terminal_exon.fields[5]).split(',')[0])
    
                # write bed file
                w_terminal.write("\t".join([str(terminal_exon.fields[0]),
                                    str(terminal_exon.fields[1]),
                                    str(terminal_exon.fields[2]),
                                    ",".join([gene_id, overlap]),
                                    '.',
                                    str(terminal_exon.fields[3])+"\n"]))
        w_terminal.close()


    def determine_clean_terminal_exons(self, terminal_exons_clean_bed_file_path, terminal_exons_bed_file_path, intermediate_exons_bed_file_path):

        """Find terminal exons that do not overlap with intermediate exons"""
        
        terminal_exons_bed = pybedtools.BedTool(terminal_exons_bed_file_path)
        intermediate_exons_bed = pybedtools.BedTool(intermediate_exons_bed_file_path)

        terminal_exons_bed_clean = terminal_exons_bed.intersect(intermediate_exons_bed, s = True, v = True).saveas(terminal_exons_clean_bed_file_path)


    def determine_clean_intermediate_exons(self, intermediate_exons_clean_bed_file_path, terminal_exons_bed_file_path, intermediate_exons_bed_file_path):

        """Find intermediate exons that do not overlap with terminal exons"""
        
        terminal_exons_bed = pybedtools.BedTool(terminal_exons_bed_file_path)
        intermediate_exons_bed = pybedtools.BedTool(intermediate_exons_bed_file_path)

        intermediate_exons_bed_clean = intermediate_exons_bed.intersect(terminal_exons_bed, s = True, v = True).saveas(intermediate_exons_clean_bed_file_path)


    def create_background_with_intergenic_polya_sites(self,
                                                      genes_bed_file,
                                                      final_terminal_exons_bed_file,
                                                      intergenic_polyasites_stranded_bed_file,
                                                      intergenic_polyasites_unstranded_bed_file,
                                                      background_regions_stranded_bed_file,
                                                      background_regions_unstranded_bed_file,
                                                      bases_to_extend = 10000):

        """Construct artificial terminal exons to be used as background by combining the 
           terminal exons and the intergenic poly(A) sites
        """

        # _____________________________________________________________________
        # STEP 1
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Identify final terminal exons of the gene.
        # Final means that the exon end overlaps with the gene end for + strand
        # or the exon start overlaps with the gen start for - strand
        # ---------------------------------------------------------------------

        w = open(final_terminal_exons_bed_file+"_tmp", 'w')
        for gene in self.genes:
            for exon in self.genes[gene].get_final_terminal_exons():
                w.write(exon.to_bed())
        w.close()
        # remove duplicates
        os.system("sort -u %s > %s" % (final_terminal_exons_bed_file+"_tmp", final_terminal_exons_bed_file))
        os.system("rm  %s" % (final_terminal_exons_bed_file+"_tmp"))

        # _____________________________________________________________________
        # STEP 2
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Construct extended regions using the final terminal exons
        # ---------------------------------------------------------------------
        final_terminal_exons = pybedtools.BedTool(final_terminal_exons_bed_file)

        extended_terminal_regions = os.path.join(self.__tmp, "extended_terminal_regions_temp.bed")
        w = open(extended_terminal_regions, 'w')
        for final_terminal_exon  in final_terminal_exons:
            if final_terminal_exon.strand == "+":
                final_terminal_exon.end = final_terminal_exon.end + int(bases_to_extend)
            elif final_terminal_exon.strand == "-":
                final_terminal_exon.start = final_terminal_exon.start - int(bases_to_extend)
            w.write(str(final_terminal_exon))
        w.close()

        # load extended part to memory
        extended_terminal_regions_mem =  pybedtools.BedTool(extended_terminal_regions)

        # _____________________________________________________________________
        # STEP 3
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Create extended terminal exons for stranded protocols
        # ---------------------------------------------------------------------

        # load stranded intergenic polya sites to memory
        intergenic_polyasites_stranded = pybedtools.BedTool(intergenic_polyasites_stranded_bed_file)

        # find overlap between the extended region and the intergenic polya sites
        overlap_extended_regions_and_intergenic_polya_sites = extended_terminal_regions_mem.intersect(intergenic_polyasites_stranded, s=True, wa=True, wb=True)

        # create files that contain only the extended exons (background regions)
        w = open(intergenic_polyasites_stranded_bed_file+"_tmp", 'w')
        for overlap in overlap_extended_regions_and_intergenic_polya_sites:

            if (str(overlap.fields[5]) == '+') and  (str(overlap.fields[11]) == '+'):

                w.write("\t".join([str(overlap.fields[0]),
                                   str(overlap.fields[1]),
                                   str(overlap.fields[7]),
                                   str(overlap.fields[3]),
                                   str(overlap.fields[4]),
                                   str(overlap.fields[5])+"\n"
                                 ]))
            elif (str(overlap.fields[5]) == '-') and  (str(overlap.fields[11]) == '-'):

                w.write("\t".join([str(overlap.fields[0]),
                                   str(overlap.fields[8]),
                                   str(overlap.fields[2]),
                                   str(overlap.fields[3]),
                                   str(overlap.fields[4]),
                                   str(overlap.fields[5])+"\n"
                                ]))
        w.close()

        # load extended exons to memory
        extended_exons = pybedtools.BedTool(intergenic_polyasites_stranded_bed_file+"_tmp")

        # load genes
        genes = pybedtools.BedTool(genes_bed_file)

        # keep the extended exons that overlap only with one gene
        w = open(background_regions_stranded_bed_file, 'w')
        for extended_exon in extended_exons.intersect(genes, s=True, c=True):

            if int(extended_exon.fields[6]) == 1:

                w.write("\t".join([str(extended_exon.fields[0]),
                                   str(extended_exon.fields[1]),
                                   str(extended_exon.fields[2]),
                                   str(extended_exon.fields[3])+",background",
                                   str(extended_exon.fields[4]),
                                   str(extended_exon.fields[5])+"\n"
                                ]))
        w.close()

        # _____________________________________________________________________
        # STEP 4
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Create extended terminal exons for unstranded protocols
        # ---------------------------------------------------------------------

        # load unstranded intergenic polya sites to memory
        intergenic_polyasites_unstranded = pybedtools.BedTool(intergenic_polyasites_unstranded_bed_file)

        # find overlap between the extended region and the intergenic polya sites
        overlap_extended_regions_and_unstranded_intergenic_polya_sites = extended_terminal_regions_mem.intersect(intergenic_polyasites_unstranded, s=True, wa=True, wb=True)

        # create files that contain only the extended exons (background regions)
        w = open(intergenic_polyasites_unstranded_bed_file+"_tmp", 'w')
        for overlap in overlap_extended_regions_and_unstranded_intergenic_polya_sites:

            if (str(overlap.fields[5]) == '+') and  (str(overlap.fields[11]) == '+'):

                w.write("\t".join([str(overlap.fields[0]),
                                   str(overlap.fields[1]),
                                   str(overlap.fields[7]),
                                   str(overlap.fields[3]),
                                   str(overlap.fields[4]),
                                   str(overlap.fields[5])+"\n"
                                 ]))
            elif (str(overlap.fields[5]) == '-') and  (str(overlap.fields[11]) == '-'):

                w.write("\t".join([str(overlap.fields[0]),
                                   str(overlap.fields[8]),
                                   str(overlap.fields[2]),
                                   str(overlap.fields[3]),
                                   str(overlap.fields[4]),
                                   str(overlap.fields[5])+"\n"
                                ]))
        w.close()

        # load extended exons to memory
        extended_exons = pybedtools.BedTool(intergenic_polyasites_unstranded_bed_file+"_tmp")

        # keep the extended exons that overlap only with one gene
        w = open(background_regions_unstranded_bed_file, 'w')
        for extended_exon in extended_exons.intersect(genes, s=True, c=True):

            if int(extended_exon.fields[6]) == 1:

                w.write("\t".join([str(extended_exon.fields[0]),
                                   str(extended_exon.fields[1]),
                                   str(extended_exon.fields[2]),
                                   str(extended_exon.fields[3])+",background",
                                   str(extended_exon.fields[4]),
                                   str(extended_exon.fields[5])+"\n"
                                ]))
        w.close()
        
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Clean up temporary files
        # ---------------------------------------------------------------------
        os.system("rm %s" % (extended_terminal_regions))
        os.system("rm %s" % (intergenic_polyasites_stranded_bed_file+"_tmp"))
        os.system("rm %s" % (intergenic_polyasites_unstranded_bed_file+"_tmp"))


    def create_background_with_intronic_distances(self,
                                                  genes_bed_file,
                                                  union_introns_bed_file,
                                                  polyasites_in_introns_stranded_bed_file,
                                                  polyasites_in_introns_unstranded_bed_file,
                                                  sequencing_direction,
                                                  final_terminal_exons_bed_file,
                                                  background_regions_bed_file):

        """ Create background regions using intronic distances """
        
        # _____________________________________________________________________
        # STEP 1
        #
        # _____________________________________________________________________
        # Create histogram of distances of intronic polya sites
        # ---------------------------------------------------------------------

        # Find to which intron each polya site falls to
        polyasites_in_introns_extended = ''

        if (sequencing_direction == "forward") or (sequencing_direction == "reverse"):

            union_introns = pybedtools.BedTool(union_introns_bed_file)
            polyasites_in_introns = pybedtools.BedTool(polyasites_in_introns_stranded_bed_file)

            polyasites_in_introns_extended = polyasites_in_introns.intersect(union_introns, s=True, wa=True, wb=True)

        elif (sequencing_direction == "unstranded"):

            union_introns = pybedtools.BedTool(union_introns_bed_file)
            polyasites_in_introns = pybedtools.BedTool(polyasites_in_introns_unstranded_bed_file)

            polyasites_in_introns_extended = polyasites_in_introns.intersect(union_introns, s=True, wa=True, wb=True)

        # dictionary of introns 
        dict_of_introns = defaultdict(list)
        # dictionary of observed intronic polya sites
        dict_of_observed_polya_sites = dict()

        for line in polyasites_in_introns_extended:

            polya_chr     = str(line.fields[0])
            polya_start   = int(line.fields[1])
            polya_end     = int(line.fields[2])
            polya_strand  = str(line.fields[5])

            intron_chr    = str(line.fields[6])
            intron_start  = int(line.fields[7])
            intron_end    = int(line.fields[8])
            intron_strand = str(line.fields[11])
            intron_gene   = line.fields[9]

            # generate intronic id
            intron_id = ":".join([intron_chr, str(intron_start), str(intron_end), intron_strand])

            # generate intronic poly(A) site id
            intronic_polya_site_id = ":".join([polya_chr, str(polya_start), str(polya_end), polya_strand])

            # if this is an intonic polya site that was not used before
            if intronic_polya_site_id not in dict_of_observed_polya_sites:

                # add it to the dicitionary
                dict_of_observed_polya_sites[intronic_polya_site_id] = intronic_polya_site_id

                # create intronic polya site object
                intronic_polya_site = HTSeq.GenomicInterval(polya_chr, polya_start, polya_end, polya_strand)

                # and append it to the dictionary of introns
                dict_of_introns[intron_id].append(intronic_polya_site)

        # calculate distances
        distances = self.calculate_distances(dict_of_introns)
        
        # create histogram of intronic polya site distances
        distances_np = np.array(distances)
        n, bins, patches = plt.hist(distances_np, 500, normed=1, facecolor='green', alpha=0.75)
        plt.savefig(os.path.join(self.__tmp, "histogram_of_intronic_polya_sites.pdf"))

        # _____________________________________________________________________
        # STEP 2
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Identify final terminal exons of the gene.
        # Final means that the exon end overlaps with the gene end for + strand
        # or the exon start overlaps with the gene start for - strand
        # ---------------------------------------------------------------------

        w = open(final_terminal_exons_bed_file+"_tmp", 'w')
        for gene in self.genes:
            for exon in self.genes[gene].get_final_terminal_exons():
                w.write(exon.to_bed())
        w.close()
        # remove duplicates
        os.system("sort -u %s > %s" % (final_terminal_exons_bed_file+"_tmp", final_terminal_exons_bed_file))
        os.system("rm  %s" % (final_terminal_exons_bed_file+"_tmp"))

        # _____________________________________________________________________
        # STEP 3
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Create extended terminal exons
        # ---------------------------------------------------------------------

        # load terminal exons in memory
        terminal_exons = pybedtools.BedTool(final_terminal_exons_bed_file)

        # 
        w = open(background_regions_bed_file+"_tmp", 'w')
        for terminal_exon in terminal_exons:

            random_extention = int(random.choice(distances))

            if (str(terminal_exon.fields[5]) == '+'):

                w.write("\t".join([str(terminal_exon.fields[0]),
                    str(terminal_exon.fields[1]),
                    str(int(terminal_exon.fields[2])+random_extention),
                    str(terminal_exon.fields[3]),
                    str(terminal_exon.fields[4]),
                    str(terminal_exon.fields[5])+"\n"
                ]))


            elif (str(terminal_exon.fields[5]) == '-'):

                minus_start = int(terminal_exon.fields[1])-random_extention
                if minus_start < 1:
                    minus_start = 1

                w.write("\t".join([str(terminal_exon.fields[0]),
                    str(minus_start),
                    str(terminal_exon.fields[2]),
                    str(terminal_exon.fields[3]),
                    str(terminal_exon.fields[4]),
                    str(terminal_exon.fields[5])+"\n"
                ]))
        w.close()

        # _____________________________________________________________________
        # STEP 4
        #
        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Remove regions that overlap with other genes
        # ---------------------------------------------------------------------

        # read poterntial extended exons
        extended_terminal_exons = pybedtools.BedTool(background_regions_bed_file+"_tmp")

        # read gene coordinates  in memory
        genes = pybedtools.BedTool(genes_bed_file)

        strand = True
        if sequencing_direction == "unstranded":
            strand = False


        # write extended exons that overlap only with 1 gene
        w = open(background_regions_bed_file, 'w')
        for extended_exon in extended_terminal_exons.intersect(genes, s=strand, c=True):

            if int(extended_exon.fields[6]) == 1:

                w.write("\t".join([str(extended_exon.fields[0]),
                                   str(extended_exon.fields[1]),
                                   str(extended_exon.fields[2]),
                                   str(extended_exon.fields[3])+",background",
                                   str(extended_exon.fields[4]),
                                   str(extended_exon.fields[5])+"\n"
                ]))
        w.close()

        # remove unused files
        os.system("rm  %s" % (background_regions_bed_file+"_tmp"))


    def determine_CDS_novel_transcrtips_from_spicing(self, genome):

        """Determine CDS sequence and stop codon for novel transcripts"""

        stop_codons = ['TAA','TAG','TGA']

        # go over all the genes 
        for gene in self.genes:

            # and all the novel transcritps
            for transcript in self.genes[gene].novel_transcripts:

                start_codon_observed = False
                stop_codon_observed  = False
                start_codon_observed = transcript.contains_start_codon_for_novel_transcripts()
                stop_codon_observed  = transcript.contains_stop_codon_for_novel_transcripts()

                # Start and stop codon were observed. This means that the novel 
                # exon falls in the 3'UTR. The existing (annotated) CDS will
                # be used without any changes.
                if start_codon_observed is True and stop_codon_observed is True:
                    
                    transcript.write_CDS = True

                # Only transcsripts that have a start codon and not a stop codon
                # are considered here. These are the ones that might have a novel
                # CDS/stop.
                elif start_codon_observed is True and stop_codon_observed is False:

                    
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
                        # what we do is to change the start of the CDS based one the frame (CDS.start + frame)
                        penultimate_CDS_seq = genome.sequence({'chr': penultimate_exon.chromosome,
                                                               'start': int(penultimate_exon.CDS.start)+int(penultimate_exon.frame),
                                                               'stop': penultimate_exon.CDS.end,
                                                               'strand':penultimate_exon.strand},
                                                               one_based=False)

                        last_exon_seq = genome.sequence({'chr': last_exon.chromosome,
                                                         'start': int(last_exon.start),
                                                         'stop': int(last_exon.end),
                                                         'strand': last_exon.strand},
                                                         one_based=False)

                        # Use the modulo to check if we get 0,1, or 2 (number of bases left)
                        bases_left = len(penultimate_CDS_seq) % 3

                        possible_last_exon_frame = None

                        # determine frame for novel exon
                        if bases_left == 0:

                            possible_last_exon_frame = 0
                            previous_bases = '' # keep the bases (sequence) left from the previous exon

                        elif bases_left == 1:

                            possible_last_exon_frame = 2
                            previous_bases = penultimate_CDS_seq[-bases_left:] # keep the bases (sequence) left from the previous exon

                        elif bases_left == 2:

                            possible_last_exon_frame = 1
                            previous_bases = penultimate_CDS_seq[-bases_left:] # keep the bases (sequence) left from the previous exon

                        # merge the bases left from the previous exon/CDS with the last exon
                        possible_cds = previous_bases + last_exon_seq

                        # Position of stop codon that will be observed
                        first_stop_codon_position = None

                        # Search for a stop codon
                        for i in range(0,len(possible_cds),3):

                            # print(i)
                            # print(possible_cds[i:i+3])

                            # stop codon was detected
                            if possible_cds[i:i+3].upper() in stop_codons:
                                first_stop_codon_position = i
                                break

                        # A stop codon was observed. The CDS for the novel transcript is calculated.
                        if first_stop_codon_position is not None:

                            if bases_left > first_stop_codon_position:
                                
                                transcript.write_CDS =  False
                                sys.stderr.write("-"*80+"\n")
                                sys.stderr.write("[WARNING] Manually check this exon for translation.\n")
                                sys.stderr.write("-"*80+"\n")
                                sys.stderr.write("chromosome: " + str(last_exon.chromosome)+"\n")
                                sys.stderr.write("start: " + str(last_exon.start)+"\n")
                                sys.stderr.write("end: " + str(last_exon.end)+"\n")
                                sys.stderr.write("strand: " + str(last_exon.strand)+"\n")
                                sys.stderr.write("first stop codon position: " + str(first_stop_codon_position)+"\n")
                                sys.stderr.write("bases left: " +str(bases_left)+"\n")
                                sys.stderr.write("possible last exon frame: " + str(possible_last_exon_frame)+"\n")
                                sys.stderr.write("-"*80+"\n\n")
                                continue

                            if last_exon.strand == "+":

                                # CDS for the novel exon
                                CDS = HTSeq.GenomicInterval(last_exon.chromosome, 
                                                            last_exon.start, 
                                                            last_exon.start + first_stop_codon_position - bases_left, 
                                                            last_exon.strand)

                                # stop codon for the novel exon
                                stop_codon = HTSeq.GenomicInterval(last_exon.chromosome, 
                                                                   last_exon.start + first_stop_codon_position - bases_left, 
                                                                   last_exon.start + first_stop_codon_position - bases_left + 2, 
                                                                   last_exon.strand)

                                last_exon.CDS = CDS
                                last_exon.frame = possible_last_exon_frame
                                last_exon.stop_codon = stop_codon

                                transcript.write_CDS = True

                            elif last_exon.strand == "-":

                                # CDS for the novel exon
                                CDS = HTSeq.GenomicInterval(last_exon.chromosome,
                                                            last_exon.end - first_stop_codon_position + bases_left,
                                                            last_exon.end,
                                                            last_exon.strand)

                                # stop codon for the novel exon
                                stop_codon = HTSeq.GenomicInterval(last_exon.chromosome, 
                                                            last_exon.end - first_stop_codon_position + bases_left - 2,
                                                            last_exon.end - first_stop_codon_position + bases_left,
                                                            last_exon.strand)


                                # update information for the last exon
                                last_exon.CDS = CDS
                                last_exon.frame = possible_last_exon_frame
                                last_exon.stop_codon = stop_codon

                                transcript.write_CDS = True

                            # print("**************************************************************************")
                            # print("******************************* STATISTICS *******************************")
                            # print("**************************************************************************")
                            # print("Last Exon      :", last_exon.start, last_exon.end)
                            # print("Last Exon CDS  :", last_exon.CDS.start, last_exon.CDS.end)
                            # print("Last Exon frame:", last_exon.frame)
                            # print("**************************************************************************\n\n")

                        else:

                            transcript.write_CDS = False


                # No stop codon was found, so the CDS should not be written for
                # this transcript
                else:

                    transcript.write_CDS = False

        
    def find_CDS_for_novel_transcripts(self, genome):

        """Determine CDS sequence and stop codon for novel transcripts"""

        stop_codons = ['TAA','TAG','TGA']

        # loop over all genes
        for gene in self.genes:

            # loop over the novel transcripts 
            for transcript in self.genes[gene].novel_transcripts:

                cds_sequence = '' 
                start_codon_observed = False
                stop_codon_observed = False

                # Loop over the exons of the novel transcript
                for exon in transcript.novel_exons: # MAKE SURE THAT THE NOVEL EXONS ARE SORTED PROPERLY (novel exon should be last)

                    if exon.start_codon is not None:
                        start_codon_observed = True

                    if exon.stop_codon is not None:
                        stop_codon_observed = True

                    if exon.CDS is not None:
                        cds_sequence = cds_sequence + genome.sequence({'chr': exon.CDS.chrom, 'start': exon.CDS.start, 'stop': exon.CDS.end, 'strand':exon.CDS.strand}, one_based=False)
                        # print(cds_sequence)

                # Start and stop codon were observed. This means that the novel exon falls in the 3'UTR. The existing (annotated) CDS will be used without any changes.
                if start_codon_observed is True and stop_codon_observed is True:
                    transcript.write_CDS = True

                # Only transcsripts that have a start codon and not a stop codon are considered here. These are the ones that might have a novel CDS/stop.
                elif start_codon_observed is True and stop_codon_observed is False:

                    # Number of bases from previous exon to not disturb the reading frame (0,1,2)
                    number_of_prev_bases = len(cds_sequence) % 3
                    # print(number_of_prev_bases)

                    # Position of stop codon that will be observed
                    first_stop_codon_position = None

                    possible_cds = cds_sequence[-number_of_prev_bases] + genome.sequence({'chr': exon.chromosome, 'start': exon.start, 'stop': exon.end, 'strand':exon.strand}, one_based=False)
                    # print(possible_cds)

                    # Search sequence for stop codon (within the open reading frame)
                    for i in range(0,len(possible_cds),3):

                        # stop codon was detected
                        if possible_cds[i:i+3].upper() in stop_codons:
                            first_stop_codon_position = i
                            break

                    # A stop codon was observed. The CDS for the novel transcript is calculated
                    if first_stop_codon_position is not None:

                        # print("First stop codon", first_stop_codon_position)

                        # This is a case we should skip because no CDS/stop codon can be determined. TEST THIS AGAIN.
                        if first_stop_codon_position == 0:
                            transcript.write_CDS = False
                        else:
                            if exon.strand == "+":
                                # this is the CDS including possible bases from the previous exon
                                # actual_CDS = HTSeq.GenomicInterval(exon.chromosome, exon.start-number_of_prev_bases, exon.start-number_of_prev_bases+first_stop_codon_position, exon.strand)

                                # CDS for the novel exon
                                CDS = HTSeq.GenomicInterval(exon.chromosome, exon.start, exon.start+first_stop_codon_position-1, exon.strand)
                                # stop codon in novel exon
                                stop_codon = HTSeq.GenomicInterval(exon.chromosome, exon.start+first_stop_codon_position-1, exon.start+first_stop_codon_position+2, exon.strand)
                            else:

                                # actual_CDS = HTSeq.GenomicInterval(exon.chromosome, exon.end+number_of_prev_bases-first_stop_codon_position, exon.end+number_of_prev_bases, exon.strand)

                                CDS = HTSeq.GenomicInterval(exon.chromosome, exon.end-first_stop_codon_position+1, exon.end, exon.strand)
                                stop_codon = HTSeq.GenomicInterval(exon.chromosome, exon.end-first_stop_codon_position-2, exon.end-first_stop_codon_position+1, exon.strand)

                            # update information for the last exon
                            exon.CDS = CDS
                            if number_of_prev_bases == 0:
                                exon.frame = 0
                            elif number_of_prev_bases == 1:
                                exon.frame = 2
                            elif number_of_prev_bases == 2:
                                exon.frame = 1
                            exon.stop_codon = stop_codon

                            transcript.write_CDS = True

                        # # # # # 
                        # TESTS #
                        # # # # #

                        # print("************************")
                        # print(">actual_CDS", actual_CDS)
                        # print(genome.sequence({'chr': actual_CDS.chrom, 'start': actual_CDS.start, 'stop':actual_CDS.end, 'strand':actual_CDS.strand}, one_based=False))
                        # print(">CDS", CDS)
                        # print(genome.sequence({'chr': CDS.chrom, 'start': CDS.start, 'stop':CDS.end, 'strand':CDS.strand}, one_based=True))
                        # print(">stop_codon", stop_codon)
                        # print(genome.sequence({'chr': stop_codon.chrom, 'start': stop_codon.start, 'stop':stop_codon.end, 'strand':stop_codon.strand}, one_based=True))
                        # print("************************")

                    # No stop codon was found, so the CDS should not be written for this transcript
                    else:

                        if (len(possible_cds) % 3 == 2) and ((possible_cds[-2] == "TA") or (possible_cds[-2] == "TG")):
                            print("Cool 1")
                            print(gene)
                            transcript.write_CDS = True
                        elif (len(possible_cds) % 3 == 1) and (possible_cds[-1] == "T"):
                            print("COOOOL 2")
                            print(gene)
                            transcript.write_CDS = True
                        else: 
                            transcript.write_CDS = False

                # No CDS or UTR information is available (Do not write CDS)
                else:
                    transcript.write_CDS = False


    def find_CDS_for_novel_exteded_transcripts(self, genome):

        """Determine CDS sequence and stop codon for novel extended transcripts"""

        stop_codons = ['TAA','TAG','TGA']

        # loop over all genes
        for gene in self.genes:

            # loop over the novel transcripts 
            for transcript in self.genes[gene].novel_transcripts_readthrough:

                cds_sequence = '' 
                start_codon_observed = False
                stop_codon_observed = False

                # Loop over the exons of the novel transcript
                for exon in transcript.novel_exons_read_through: # MAKE SURE THAT THE NOVEL EXONS ARE SORTED PROPERLY (novel exon should be last)

                    if exon.start_codon is not None:
                        start_codon_observed = True

                    if exon.stop_codon is not None:
                        stop_codon_observed = True

                    if exon.CDS is not None:
                        cds_sequence = cds_sequence + genome.sequence({'chr': exon.CDS.chrom, 'start': exon.CDS.start, 'stop': exon.CDS.end, 'strand':exon.CDS.strand}, one_based=False)

                # Start and stop codon were observed. This means that the novel exon falls in the 3'UTR. The existing (annotated) CDS will be used without any changes.
                if start_codon_observed is True and stop_codon_observed is True:
                    transcript.write_CDS = True

                # Only transcsripts that have a start codon and not a stop codon are considered here. These are the ones that might have a novel CDS/stop.
                elif start_codon_observed is True and stop_codon_observed is False:

                    # Number of bases from previous exon to not disturb the reading frame (0,1,2)
                    number_of_prev_bases = len(cds_sequence) % 3

                    # Position of stop codon that will be observed
                    first_stop_codon_position = None

                    possible_cds = cds_sequence[-number_of_prev_bases] + genome.sequence({'chr': exon.chromosome, 'start': exon.start, 'stop': exon.end, 'strand':exon.strand}, one_based=False)

                    # Search sequence for stop codon (within the open reading frame)
                    for i in range(0,len(possible_cds),3):

                        # stop codon was detected
                        if possible_cds[i:i+3].upper() in stop_codons:
                            first_stop_codon_position = i
                            break

                    # A stop codon was observed. The CDS for the novel transcript is calculated
                    if first_stop_codon_position is not None:

                        # This is a case we should skip because no CDS/stop codon can be determined. TEST THIS AGAIN.
                        if first_stop_codon_position == 0:
                            transcript.write_CDS = False
                        else:
                            if exon.strand == "+":
                                # this is the CDS including possible bases from the previous exon
                                # actual_CDS = HTSeq.GenomicInterval(exon.chromosome, exon.start-number_of_prev_bases, exon.start-number_of_prev_bases+first_stop_codon_position, exon.strand)

                                # CDS for the novel exon
                                CDS = HTSeq.GenomicInterval(exon.chromosome, exon.start, exon.start+first_stop_codon_position-1, exon.strand)
                                # stop codon in novel exon
                                stop_codon = HTSeq.GenomicInterval(exon.chromosome, exon.start+first_stop_codon_position-1, exon.start+first_stop_codon_position+2, exon.strand)
                            else:

                                # actual_CDS = HTSeq.GenomicInterval(exon.chromosome, exon.end+number_of_prev_bases-first_stop_codon_position, exon.end+number_of_prev_bases, exon.strand)

                                CDS = HTSeq.GenomicInterval(exon.chromosome, exon.end-first_stop_codon_position+1, exon.end, exon.strand)
                                stop_codon = HTSeq.GenomicInterval(exon.chromosome, exon.end-first_stop_codon_position-2, exon.end-first_stop_codon_position+1, exon.strand)

                            # update information for the last exon
                            exon.CDS = CDS
                            exon.frame = number_of_prev_bases
                            exon.stop_codon = stop_codon

                            transcript.write_CDS = True

                        # # # # # 
                        # TESTS #
                        # # # # #

                        # print("************************")
                        # print(">actual_CDS", actual_CDS)
                        # print(genome.sequence({'chr': actual_CDS.chrom, 'start': actual_CDS.start, 'stop':actual_CDS.end, 'strand':actual_CDS.strand}, one_based=False))
                        # print(">CDS", CDS)
                        # print(genome.sequence({'chr': CDS.chrom, 'start': CDS.start, 'stop':CDS.end, 'strand':CDS.strand}, one_based=True))
                        # print(">stop_codon", stop_codon)
                        # print(genome.sequence({'chr': stop_codon.chrom, 'start': stop_codon.start, 'stop':stop_codon.end, 'strand':stop_codon.strand}, one_based=True))
                        # print("************************")

                    # No stop codon was found, so the CDS should not be written for this transcript...
                    else:
                        
                        # ... except there is a novel stop codon creatd with the poly(A) tail used
                        if (len(possible_cds) % 3 == 2) and ((possible_cds[-2] == "TA") or (possible_cds[-2] == "TG")):
                            transcript.write_CDS = True
                        elif (len(possible_cds) % 3 == 1) and (possible_cds[-1] == "T"):
                            transcript.write_CDS = True
                        else: 
                            transcript.write_CDS = False

                # No CDS or UTR information is available (Do not write CDS)
                else:
                    transcript.write_CDS = False


    def write_gtf(self, output_file, with_CDS, accepted_exons_dict="write_all"):

        """ Write out GFT file """

        for gene in self.genes:

            # write gene 
            self.genes[gene].write_gtf(output_file)

            # Find which novel transcripts to use
            novel_transcripts_fitlered = self.genes[gene].decide_novel_transcripts_to_use()

            novel_and_known_merged_sorted = sorted(novel_transcripts_fitlered + self.genes[gene].get_known_transcripts(), key=lambda x: x.start, reverse=False)

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
    
                            exon_id = ":".join([exon.chromosome, str(exon.start), str(exon.end), exon.strand])
    
                            if exon_id in accepted_exons_dict:
                
                                # write transcript
                                transcript.write_transcript_gtf(output_file)
    
                                # write exons
                                transcript.write_exons_gtf(output_file, with_CDS)
    
                                break
                    else:
                        # write transcript
                        transcript.write_transcript_gtf(output_file)
    
                        # write exons
                        transcript.write_exons_gtf(output_file, with_CDS)


    def write_transcript2novel_transcript_mapping_list(self, output_file_handle, accepted_exons_dict="write_all"):

        """ 
            Function that takes as input a list of accepted novel exons 
            and writes out a list of that contains transcripts ids  and
            novel transcipt ids. The transcript ids are all the potential
            "mother" transcripts that generate the novel transcript ids.
        """

        # in case of write all then we should write all novel transcripts
        if accepted_exons_dict == "write_all":

            pass
            # for gene_id in annotation.genes:
            #     for novel_transcript_id in annotation.genes[gene_id].mother_transcripts_of_novel_transcripts:
            #         for transcript_id in annotation.genes[gene_id].mother_transcripts_of_novel_transcripts[novel_transcript_id]:
            #             print("%s \t %s" % (transcript_id, novel_transcript_id))

            # # go over all the genes
            # for gene_id in annotation.genes:
            #     # find all novel transcript ids
            #     for novel_transcript in annotation.genes[gene_id].decide_novel_transcripts_to_use():
            #         # transcript id / novel transcript id
            #         for transcript_id in annotation.genes[gene_id].mother_transcripts_of_novel_transcripts[novel_transcript.transcript_id]:
        
            #             print("%s \t %s" % (transcript_id, novel_transcript.transcript_id))

        
        elif isinstance(accepted_exons_dict, dict):

            # go over the genes
            for gene in self.genes:

                # Find which novel transcripts to use
                novel_transcripts_fitlered = self.genes[gene].decide_novel_transcripts_to_use()

                for novel_transcript in novel_transcripts_fitlered:

                    if "TECtool" in novel_transcript.source:

                        for exon in novel_transcript.novel_exons:

                            exon_id = ":".join([exon.chromosome, str(exon.start), str(exon.end), exon.strand])
    
                            if exon_id in accepted_exons_dict:

                                novel_transcript_id = novel_transcript.transcript_id

                                for transcript_id in self.genes[gene].mother_transcripts_of_novel_transcripts[novel_transcript_id]:

                                    output_file_handle.write("%s \t %s\n" % (transcript_id, novel_transcript_id))

                                break


    def write_transcript2gene_mapping_list(self, transcript2gene_mapping_file_path):

        """Creates a file that contains transcript to gene mappings for 
        all transcripts that exist in the Annotation object."""

        # check whether self.transcripts exists.
        if not self.transcripts:
            sys.stderr.write(("ERROR: 'Annotation.transcripts' is empty! "\
                             +"Thus, writing '%s' is not possible!\n") \
                             % (transcript2gene_mapping_file_path))
            sys.exit(-1)

        # open the file for writing
        mapping_file_handle = open(transcript2gene_mapping_file_path, 'w')
        file_header = "transcript_id\tgene_id\tgene_name\n"
        mapping_file_handle.write(file_header)

        # write out mappings for all transcripts that exist in the Annotation object.
        for transcript_id in self.transcripts:
            
            # create the line to write
            line_to_write = self.transcripts[transcript_id].transcript_id + "\t" \
                          + self.transcripts[transcript_id].gene_id + "\t" \
                          + self.transcripts[transcript_id].gene_name + "\n"
            
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
            fh.write("\t".join([terminal_exon.chromosome, str(terminal_exon.start), str(terminal_exon.end), terminal_exon.gene_id, ".", terminal_exon.strand+"\n"]))
        fh.close()


    def write_all_first_and_intermediate_exons_as_bed(self, outfile):

        """ 
        Writes a file with all first and intermediate exons.
        It can contain duplicate terminal exons.
        """

        all_first_and_intermediate_exons = self.get_all_first_and_intermediate_exons()

        fh = open(outfile, 'w')
        for exon in all_first_and_intermediate_exons:
            fh.write("\t".join([exon.chromosome, str(exon.start), str(exon.end), exon.gene_id, ".", exon.strand+"\n"]))
        fh.close()
