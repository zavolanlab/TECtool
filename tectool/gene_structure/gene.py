import sys
from collections import defaultdict
from itertools import chain
flatten = chain.from_iterable
from interval import Interval, IntervalSet


class Gene(object):

    """This class represents a Gene

        :param gene_id: unique id for a gene (e.g. Ensembl Gene ID).
        :param chromosome: the chromosome on which the gene is located on.
        :param source: the source of the gene annotation 
            (as provided in the GTF/GFF file).
        :param start: the start of the gene
        :param end: the end of the gene
        :param strand: the end of the gene
        :rtype: Gene object

    *Class members*


        *chromosome*
            String. The chromosome at which the gene is located on.
        
        *source*
            String. The name of the program that generated the annotation 
                feature, or the data source (database or project name).

        *feature*
            String. feature type name, e.g. Gene, Transcript, Exon

        *start*
            Int. Start position of the gene.

        *end*
            Int. End position of the gene.

        *score*
            String. Score.

        *strand*
            String. Strand of the gene (defined as + (forward) or - (reverse)).

        *frame*
            String. It can be "."" which is no frame info or 0/1/2.
             '0' indicates that the first base of the feature is 
             the first base of a codon, '1' that the second base is
             the first base of a codon, and so on.

        *gene_id*
            String. Unique id for a gene (e.g. Ensembl ID).            

        *gene_biotype*
            String. The biotype of the gene (e.g. protein_coding).

        *gene_name*
            String. The gene symbol (e.g. FAM188A).

        *exon_list_sorted_by_end_coord*
            dictionary FIXME

        *exon_coordinates_dict*
            Dictionary. 
            key: exon coordinates in the format
            chromosome:start:end:strand. 
            value: list of exons with the doordinates of the key

        *transcripts*
            Dict. Annotated transcsripts.
            key: transcript_id
            value: transcript object

        *novel_transcripts*
            List. Novel transcritps that occur from novel splicing. 
                  Empty in case of no novel transcripts.

        *novel_transcripts_readthrough* 
            List. Novel transcritps that occur from read through
                  Empty in case of no novel transcripts.

        *annotated_terminal_exons*
            List = List of annotated terminal exons.

        *annotated_intermediate_exons*
            List.  List of annotated intermediate exons.

        *background*
            List.  List of background regions.
        
        *annotated_introns_without_polya_site*
            List. List of annotated introns without polya sites.
        
        *annotated_introns_with_polya_site*
            List. List of annotated introns with polya sites.

        *potential_novel_exons*
            List. List of potential novel exons.

        *potential_novel_readthroughs*
            List. List of potential novel readthough exons.
        
        *intron_length*
            Int. Sum of lenhts of the introns.

        *intron_reads*
            Int. Total number of reads that fall in introns.

        *union_exon_length*
            Int. Length of the union exons.

        *total_reads*
            Int. Total number of reads that fall in introns and exons.

        *GeneExpressionPerKBApproximated*
            Float. Gene expression per KB approximated
            (total reads fall in gene loci  / union exon length)*1000.

        *BackgroundPerKB*
            Float. Background per KB.

        *GeneExpressionPerKBInclBackground*
            Float. Gene expression including background.
        
        *GeneExpressionBackgroundPerKB*
            Float. Gene expression without background.

        *GeneExpressionPerKBwithoutBackground*
            Float. Gene expression without background.

        *BackgroundFraction*
            Float. Fraction of background (relative to the total gene expression).

        *overlaps_with_other_gene*
            Boolean. Flag to specify if a gene overlaps with some other. (Default True)

        *mother_transcripts_of_novel_transcripts*
            defaultdict(list). Dictionary that contains as a key the novel transcript 
            ids and as value the list with all transcirpt ids from which the novel 
            transcript can originate from.

    """

    def __init__(self, chromosome, source, feature, start, end, score, strand, frame, gene_id, gene_biotype, gene_name):

        """Constructor for a Gene object."""

        # Basic gtf info [Required]
        self.chromosome = str(chromosome)
        self.source = str(source)
        self.feature = str(feature)
        self.start = int(start)
        self.end = int(end)
        self.score = str(score)
        self.strand = str(strand)
        self.frame = str(frame)

        # Extra gtf info [Required]
        self.gene_id = str(gene_id)
        self.gene_biotype = str(gene_biotype) # in ENSEMBL is gene_biotype, in gencode is gene_type
        self.gene_name = str(gene_name)

        # dictionary that contains as key exon coordinates in the format:
        # chromosome:start:exon:end:strand 
        # and as values a list of exon ids
        self.exon_coordinates_dict = dict()

        # Dictionary of known transcripts. 
        self.transcripts = dict()
        
        # list of novel transcritps (that occur from novel splicing)
        self.novel_transcripts = list()

        # list of novel transcripts (that occur from read through)
        self.novel_transcripts_readthrough = list()

        # list of annotated terminal exons
        self.annotated_terminal_exons = list()

        # list of annotated intermediate exons
        self.annotated_intermediate_exons = list()

        # list of background regions
        self.background = list()

        # list of annotated introns without polya sites
        self.annotated_introns_without_polya_site = list()

        # list of annotated introns with polya sites
        self.annotated_introns_with_polya_site = list()

        # list of potential novel exons
        self.potential_novel_exons = list()

        # list of potential novel readthough exons
        self.potential_novel_readthroughs = list()

        # intron length
        self.intron_length = 0 # IntronLength

        # intronic reads
        self.intron_reads = 0 # IntronicReads

        # union exon length
        self.union_exon_length = 0 # UnionExonLength

        # number of reads that fall in the gene
        self.total_reads = 0 # TotalGeneReads

        # GeneExpressionPerKBApproximated
        self.GeneExpressionPerKBApproximated = 0

        # BackgroundPerKB
        self.BackgroundPerKB = 0

        # GeneExpressionPerKBInclBackground
        self.GeneExpressionPerKBInclBackground = 0

        # GeneExpressionBackgroundPerKB
        self.GeneExpressionBackgroundPerKB = 0

        # GeneExpressionPerKBwithoutBackground
        self.GeneExpressionPerKBwithoutBackground = 0

        # BackgroundFraction
        self.BackgroundFraction = 0

        # Flag to specify if the gene overlaps with some other gene
        # default True
        self.overlaps_with_other_gene = True

        # Dictionary that contains as a key the novel transcript
        # ids and as value the list with all transcirpt ids
        # from which the novel transcript can originate from.
        self.mother_transcripts_of_novel_transcripts = defaultdict(list)


    def __repr__(self):
      
        """Create a representation of the current object."""
      
        repr_string = ""
        for attr in vars(self):
            repr_string += ("%s\tobj.%s = %s\n" 
                         % (type(getattr(self, attr)), 
                            attr, getattr(self, attr)))

        return(repr_string)


    def __str__(self):

        """Create a readable string representation of the current object."""

        str_representation = \
            + 80*"_" + "\n" \
            + 80*"-" + "\n"
        str_representation += \
            ("Gene(object):\t" + self.gene_id + "\n")
        str_representation += \
            + 80*"-" + "\n"
        str_representation += \
            self.__repr__()
        str_representation += \
            + 80*"-" + "\n"

        return(str_representation)


    def __eq__(self, gene):

        """Equality operator."""

        if isinstance(gene, self.__class__):
            return((self.chromosome == gene.chromosome) and 
                   (self.source == gene.source) and            
                   (self.feature == gene.feature) and          
                   (self.start == gene.start) and
                   (self.end == gene.end) and
                   (self.score == gene.score) and              
                   (self.strand == gene.strand) and
                   (self.frame == gene.frame) and
                   (self.gene_id == gene.gene_id) and
                   (self.gene_biotype == gene.gene_biotype) and
                   (self.gene_name == gene.gene_name))

        return NotImplemented


    def __ne__(self, gene):

        """Non-equality operator."""

        if isinstance(gene, self.__class__):
            return(not self.__eq__(gene))
        return NotImplemented


    def extend(self, gene, verbose=False):

        """Method that extends the gene by missing transcripts of an equal gene."""


        # check if we have the same genes here
        if (self != gene):
            sys.stderr.write("ERROR: Gene '" + self.gene_id \
                            +"' cannot be extended by gene '" + gene.gene_id \
                            +"' because the genes differ in annotation.\n")
            sys.exit(-1)
        
        else:

            # check for every transcript in the annotation whether it exists already
            for transcript_id in gene.transcripts:

                if transcript_id in self.transcripts:
                    # if the gene exists already, extend it by the transcripts (if necessary)
                    self.transcripts[transcript_id].extend(gene.transcripts[transcript_id], 
                                                           verbose=verbose)
                
                else:
                    # add the transcript to the gene
                    if verbose:
                        sys.stdout.write(" :: adding transcript '%s'\n" \
                                         % (transcript_id))
                    self.transcripts[transcript_id] = gene.transcripts[transcript_id]


    def determine_earliest_transcript_end(self):

        """Function that parses the transcripts of the genes and returns the earliest end of a transcript"""
        
        earliest_transcript_end = None

        if self.strand is "+":
            for transcript in self.transcripts:
                if (transcript.end < earliest_transcript_end) or (earliest_transcript_end is None):
                    earliest_transcript_end = transcript.end

        elif self.strand is "-":
            
            for transcript in self.transcripts:
                if (transcript.start > earliest_transcript_end) or (earliest_transcript_end is None):
                    earliest_transcript_end = transcript.start

        else:
            sys.stderr.write("Something went wrong with the strand")
            sys.exit(-1)
        
        return(earliest_transcript_end)


    def insert_exon_into_exon_coordinates_dict(self, exon):

        """Insert exons into exon coordinates dictionary"""

        exon_coordinates = ":".join([exon.chromosome, str(exon.start), str(exon.end), exon.strand])

        # exon coordinates already in dictionary
        if str(exon_coordinates) in self.exon_coordinates_dict:
            if exon.exon_id not in self.exon_coordinates_dict[exon_coordinates]:
                self.exon_coordinates_dict[exon_coordinates].append(exon.exon_id)
        else:
            self.exon_coordinates_dict[exon_coordinates] = [exon.exon_id]


    def has_annotated_terminal_exon(self):

        """Gene has annotated terminal exons"""

        if len(self.annotated_terminal_exons) > 0:
            return(True)
        else:
            return(False)



    def has_potential_novel_terminal_exon(self):

        """Gene has potential novel terminal exons"""

        if len(self.potential_novel_exons) > 0:
            return(True)
        else:
            return(False)


    def has_potential_novel_readthrough_exon(self):

        """Gene has potential novel readthrough terminal exons"""

        if len(self.potential_novel_readthroughs) > 0:
            return(True)
        else:
            return(False)


    def get_final_terminal_exons(self):

        """Get last terminal exon of the gene"""

        all_terminal_exons = []

        start, end = self.get_actual_gene_coordinates()

        for transcript in self.transcripts:

            terminal_exon = transcript.get_terminal_exon()

            if transcript.strand == "+":

                if int(terminal_exon.end) == int(end):

                    all_terminal_exons.append(terminal_exon)

            elif transcript.strand == "-":

                if int(terminal_exon.start) == int(start):

                    all_terminal_exons.append(terminal_exon)

        return(all_terminal_exons)


    def get_known_transcripts(self):

        """Return known transcripts of the gene. Objects are stored in a list."""

        known_transcripts = []
        for transcript_id in self.transcripts:
            known_transcripts.append(self.transcripts[transcript_id])

        return(known_transcripts)


    def get_known_transctipt_ids(self):
        
        """Returng known transcript ids of the gene."""

        known_ids = []
        for transcript_id in self.transcripts:
            known_ids.append(transcript_id)

        return(known_ids)


    def order_all_transcripts(self):

        """ Function that sorts old and novel transcripts based on the start, stop and length of the transcsripts"""

        return(sorted(self.transcripts + self.novel_transcripts + self.novel_transcripts_readthrough, key=lambda x: x.start, reverse=False))


    def order_novel_transcripts(self):

        """ Function that sorts novel transcripts based on the start, stop and length of the transcsripts"""

        return(sorted(self.novel_transcripts + self.novel_transcripts_readthrough, key=lambda x: x.start, reverse=False))


    def decide_novel_transcripts_to_use(self):

        """ Function that decides which novel transcripts will be used """

        # dictionary with key the transcript id and values a list of transcripts with the same id
        transcript_id_to_transcript = defaultdict(list)

        for transcript in self.order_novel_transcripts():
            transcript_id_to_transcript.setdefault(transcript.transcript_id, []).append(transcript)

        # list that contains the final transcripts that we will keep
        novel_transcripts_keep = []

        # loop ovel the novel transctipts
        for transcript_id in transcript_id_to_transcript:

            # if we have just one transcript then we use this one
            if len(transcript_id_to_transcript[transcript_id]) == 1:
                novel_transcripts_keep.append(transcript_id_to_transcript[transcript_id][0])

            # If we have more than one transcripts then 
            elif len(transcript_id_to_transcript[transcript_id]) > 1:

                
                transcript_to_keep = None
                # we select the transcrict (and we prefer protein coding cases)
                for transcript in transcript_id_to_transcript[transcript_id]:

                    if transcript_to_keep is None:
                        transcript_to_keep = transcript
                        continue

                    if transcript.write_CDS == True:
                        transcript_to_keep = transcript
                        break

                novel_transcripts_keep.append(transcript_to_keep)

        return(novel_transcripts_keep)


    def get_actual_gene_coordinates(self):

        """ Find gene start and end based on transcript coordinates
            Return (start, end) coordinates
        """

        min_start = None
        max_end = None

        for transcript in self.transcripts:

            if min_start is None:
                min_start = int(self.transcripts[transcript].start)
            if max_end is None:
                max_end = int(self.transcripts[transcript].end)

            if int(self.transcripts[transcript].start) < int(min_start):
                min_start = int(self.transcripts[transcript].start)

            if int(self.transcripts[transcript].end) > int(max_end):
                max_end = int(self.transcripts[transcript].end)

        return(int(min_start), int(max_end))
        

    def get_actual_gene_coordinates_bed(self):

        """ Find gene start and end based on transcript coordinates
            Return line in bed format.
        """

        min_start = None
        max_end = None

        for transcript in self.transcripts:

            if min_start is None:
                min_start = int(self.transcripts[transcript].start)
            if max_end is None:
                max_end = int(self.transcripts[transcript].end)

            if int(self.transcripts[transcript].start) < int(min_start):
                min_start = int(self.transcripts[transcript].start)

            if int(self.transcripts[transcript].end) > int(max_end):
                max_end = int(self.transcripts[transcript].end)

        return("\t".join([self.chromosome, str(min_start), str(max_end), self.gene_id, ".", self.strand]))


    def union_exon_generator(self):

        """ Get union exon. Adapted from: 
            http://stackoverflow.com/questions/24317211/merge-overlapping-numeric-ranges-into-continuous-ranges
        """

        all_exons = []

        LEFT, RIGHT = 1, -1

        offset = 0
        
        for transcript in self.transcripts:

            for exon in self.transcripts[transcript].exon_list_sorted_by_end_coord:

                all_exons.append((int(exon.start), int(exon.end)))

        union_exons = sorted(flatten(((start, LEFT), (stop + offset, RIGHT)) for start, stop in all_exons))

        c = 0
        for value, label in union_exons:
            if c == 0:
                x = value
            c += label
            if c == 0:
                yield x, value - offset


    def get_union_exon_bed(self):

        """Create union exon bed file of the gene"""

        union_exons = self.union_exon_generator()

        union_exons_string = ''

        for start, end in union_exons:

            union_exons_string += "\t".join([self.chromosome, str(start), str(end), self.gene_id, ".", self.strand+"\n"])

        return(union_exons_string)


    def generate_union_introns_per_gene(self):

        """Calculate intronic coordinates by substracting the union exons coordinates from the gene coordinates.
           Based on: http://stackoverflow.com/questions/6462272/subtract-overlaps-between-two-ranges-without-setsi
        """

        # find all exonnic coordinates of the gene
        all_exons = []
        for transcript in self.transcripts:
            for exon in self.transcripts[transcript].exon_list_sorted_by_end_coord:
                all_exons.append((int(exon.start), int(exon.end)))

        # find actual gene start and end based on the annotated transcript coordinates
        gene_start, gene_end = self.get_actual_gene_coordinates()

        # generate intervals
        r1 = IntervalSet([Interval(gene_start, gene_end)])
        r2 = IntervalSet([Interval(start, end) for start, end in all_exons])
        r12 = r1 - r2 # find intronic coordinates

        # create introns list
        introns = []
        for i in r12:
            tmp_coordinates = str(i).strip("\'").strip("(").strip(")").split("..")
            introns.append("\t".join([self.chromosome, str(int(tmp_coordinates[0])), str(int(tmp_coordinates[1])), self.gene_id, '.', self.strand]))

        return(introns)

    def estimate_ExpressionPerKBApproximated(self):

        """calculate expression per KB"""

        if self.union_exon_length > 0:

            self.GeneExpressionPerKBApproximated = ((float(self.total_reads)/float(self.union_exon_length))) * 1000


    def estimate_BackgroundPerKB(self):

        """ calculate background per KB """

        if int(self.intron_length) > 0:

            self.BackgroundPerKB =  (float(self.intron_reads) / float(self.intron_length))*1000


    def estimate_GeneExpressionPerKBInclBackground(self):
        
        """ calculate gene expression including background """

        if self.union_exon_length > 0:
            self.GeneExpressionPerKBInclBackground = ((float(self.total_reads) - float(self.intron_reads))/float(self.union_exon_length)) * 1000

    def estimate_GeneExpressionBackgroundPerKB(self):

        """ calculate gene expression without background """
        
        self.GeneExpressionBackgroundPerKB = (float(self.union_exon_length)/1000)*self.BackgroundPerKB

    def estimate_GeneExpressionPerKBwithoutBackground(self):

        """ calculate the gene expression without background """

        self.GeneExpressionPerKBwithoutBackground = float(self.GeneExpressionPerKBInclBackground) - float(self.GeneExpressionBackgroundPerKB)
        

    def estimate_BackgroundFraction(self):

        """ calculate the fraction of background (relative to the total gene expression) """

        if self.GeneExpressionPerKBInclBackground > 0:
            self.BackgroundFraction = float(self.BackgroundPerKB) / float(self.GeneExpressionPerKBInclBackground)


    def get_annotated_terminal_exons(self):

        terminal_exons_list = []

        for terminal_exon in self.annotated_terminal_exons:

            terminal_exon_list = terminal_exon.get_features()

            if self.union_exon_length > 0:

                terminal_exon_list.append(self.total_reads)
                terminal_exon_list.append(self.union_exon_length)
                terminal_exon_list.append(self.GeneExpressionPerKBApproximated)
                # terminal_exon_list.append(self.intron_reads)
                # terminal_exon_list.append(self.intron_length)
                # terminal_exon_list.append(self.BackgroundPerKB)
                # terminal_exon_list.append(self.GeneExpressionPerKBInclBackground)
                # terminal_exon_list.append(self.GeneExpressionBackgroundPerKB)
                # terminal_exon_list.append(self.GeneExpressionPerKBwithoutBackground)
                # terminal_exon_list.append(self.BackgroundFraction)

                terminal_exons_list.append(terminal_exon_list)

        return(terminal_exons_list)

    def get_annotated_intermediate_exons(self):

        intermediate_exons_list = []

        for intermediate_exon in self.annotated_intermediate_exons:

            intermediate_exon_list = intermediate_exon.get_features()

            if self.union_exon_length > 0:

                intermediate_exon_list.append(self.total_reads)
                intermediate_exon_list.append(self.union_exon_length)
                intermediate_exon_list.append(self.GeneExpressionPerKBApproximated)
                # intermediate_exon_list.append(self.intron_reads)
                # intermediate_exon_list.append(self.intron_length)
                # intermediate_exon_list.append(self.BackgroundPerKB)
                # intermediate_exon_list.append(self.GeneExpressionPerKBInclBackground)
                # intermediate_exon_list.append(self.GeneExpressionBackgroundPerKB)
                # intermediate_exon_list.append(self.GeneExpressionPerKBwithoutBackground)
                # intermediate_exon_list.append(self.BackgroundFraction)

                intermediate_exons_list.append(intermediate_exon_list)

        return(intermediate_exons_list)

    def get_background(self):

        backgrounds_list = []

        for background in self.background:

            background_list = background.get_features()

            if self.union_exon_length > 0:

                background_list.append(self.total_reads)
                background_list.append(self.union_exon_length)
                background_list.append(self.GeneExpressionPerKBApproximated)
                # background_list.append(self.intron_reads)
                # background_list.append(self.intron_length)
                # background_list.append(self.BackgroundPerKB)
                # background_list.append(self.GeneExpressionPerKBInclBackground)
                # background_list.append(self.GeneExpressionBackgroundPerKB)
                # background_list.append(self.GeneExpressionPerKBwithoutBackground)
                # background_list.append(self.BackgroundFraction)

                backgrounds_list.append(background_list)

        return(backgrounds_list)


    def get_potential_novel_exons(self):

        potential_novel_exons_list = []

        for potential_novel_exon in self.potential_novel_exons:

            potential_novel_exon_list = potential_novel_exon.get_features()

            potential_novel_exon_list.append(self.total_reads)
            potential_novel_exon_list.append(self.union_exon_length)
            potential_novel_exon_list.append(self.GeneExpressionPerKBApproximated)
            # potential_novel_exon_list.append(self.intron_reads)
            # potential_novel_exon_list.append(self.intron_length)
            # potential_novel_exon_list.append(self.BackgroundPerKB)
            # potential_novel_exon_list.append(self.GeneExpressionPerKBInclBackground)
            # potential_novel_exon_list.append(self.GeneExpressionBackgroundPerKB)
            # potential_novel_exon_list.append(self.GeneExpressionPerKBwithoutBackground)
            # potential_novel_exon_list.append(self.BackgroundFraction)

            potential_novel_exons_list.append(potential_novel_exon_list)

        return(potential_novel_exons_list)

    def get_potential_novel_readthrough_exons(self):

        """Potential novel readthrough exons"""

        potential_novel_readthough_exons_list = []

        for potential_novel_exon in self.potential_novel_readthroughs:

            potential_novel_exon_list = potential_novel_exon.get_features()

            potential_novel_exon_list.append(self.total_reads)
            potential_novel_exon_list.append(self.union_exon_length)
            potential_novel_exon_list.append(self.GeneExpressionPerKBApproximated)
            # potential_novel_exon_list.append(self.intron_reads)
            # potential_novel_exon_list.append(self.intron_length)
            # potential_novel_exon_list.append(self.BackgroundPerKB)
            # potential_novel_exon_list.append(self.GeneExpressionPerKBInclBackground)
            # potential_novel_exon_list.append(self.GeneExpressionBackgroundPerKB)
            # potential_novel_exon_list.append(self.GeneExpressionPerKBwithoutBackground)
            # potential_novel_exon_list.append(self.BackgroundFraction)

            potential_novel_readthough_exons_list.append(potential_novel_exon_list)

        return(potential_novel_readthough_exons_list)


    def count_all_reads_falling_into_gene_coordinates(self, BAM_file):

        """Method that fetches all reads that entirely fall into the
           genomic coordinates of the gene."""

        return("Implement me!")


    def estimate_gene_expression_background(self, intron_ids_list):

        """Method that estimates the gene expression background from 
           a given set of introns."""
        
        return("Implement me!")


    def estimate_background_corrected_gene_expression(self):

        """Method that estimates the gene expression from union exon length,
           number of reads that fall into the genomic loci,
           number of reads that fall into intronic loci, and
           an estimated gene expression background."""
        
        return("Implement me!")


    def write_gtf(self, output_file):

        """Function that writes in GTF format the gene annotation"""

        output_file.write("\t".join([self.chromosome, self.source, self.feature, str(self.start+1), str(self.end), self.score, self.strand, self.frame, "gene_id \""+self.gene_id+"\"; gene_name \""+self.gene_name+"\"; gene_type \""+self.gene_biotype+"\";\n"]))




