import sys


class Transcript(object):

    """This class represents a Transcript

        :param transcript_id: unique id for a transcript (e.g. Ensembl Transcript ID).
        :param chromosome: the chromosome on which the transcript is located on.
        :param source: the source of the transcript annotation 
            (as provided in the GTF/GFF file).
        :param start: the start of the transcript
        :param end: the end of the transcript
        :param strand: the end of the transcript
        :param gene_id: unique id for the gene to which the transcript belongs to (e.g. Ensembl Gene ID).
        :rtype: Transcript object

    *Class members*

        *chromosome*
            String. The chromosome at which the gene is located on.
        
        *source*
            String. The name of the program that generated the annotation 
                feature, or the data source (database or project name).

        *feature*
            String. feature type name, e.g. Gene, Transcript, Exon

        *start*
            Int. Start position of the transcript.

        *end*
            Int. End position of the transcript.

        *score*
            String. Score.

        *strand*
            String. Strand of the transcript (defined as + (forward) or - (reverse)).

        *frame*
            String. It can be "."" which is no frame info or 0/1/2.
             '0' indicates that the first base of the feature is 
             the first base of a codon, '1' that the second base is
             the first base of a codon, and so on.

        *gene_id*
            String. Unique id for the gene to which the transcript belongs to (e.g. Ensembl Gene ID).

        *transcript_id*
            String. Unique id for a transcript (e.g. Ensembl Transcript ID).

        *gene_name*
            String. The gene symbol (e.g. FAM188A).

        *gene_biotype*
            String. The biotype of the gene (e.g. protein_coding).

        *transcript_name*
            String. The transcript symbol (e.g. CD47-002).

        *transcript_biotype*
            String. The biotype of the transcript (e.g. processed_transctipt).

        *exon_list_sorted_by_end_coord*
            List. List containing known exons of the transcript.

        *novel_exons*
            List. List contatining novel exons (that occur from novel splicing). Empty in case the transcript is not novel.

        *write_CDS*
            Bool. Flag of whether the CDS annotation should be written or not (True/False/None).

    """

    def __init__(self, chromosome, source, feature, start, end, score, strand, frame, gene_id, transcript_id, gene_name, gene_biotype, transcript_name, transcript_biotype):
        
        """Initialise Transcript"""

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
        self.transcript_id = str(transcript_id)
        self.gene_name = str(gene_name)
        self.gene_biotype = str(gene_biotype)
        self.transcript_name = str(transcript_name)
        self.transcript_biotype = str(transcript_biotype)

        #self.transcript_startcodon =
        #self.transcript_stopcodon =

        # List containing known exons of the transcript
        self.exon_list_sorted_by_end_coord = list()

        # List contatining novel exons (that occur from novel splicing). Empty in case the transcript is not novel.
        self.novel_exons = list()

        # Flag of whether the CDS annotation should be written or not (True/False/None)
        self.write_CDS = None


    def __repr__(self):
      
        """Create a string representation of a Transcript object."""
      
        repr_string = ""
        for attr in vars(self):
            repr_string += ("%s\tobj.%s = %s\n" 
                         % (type(getattr(self, attr)), 
                            attr, getattr(self, attr)))

        return(repr_string)


    def __str__(self):

        """Create a readable string representation of a Transcript object."""

        str_representation = \
            + 80*"_" + "\n" \
            + 80*"-" + "\n"
        str_representation += \
            ("Transcript(object):\t" + self.transcript_id + "\n")
        str_representation += \
            + 80*"-" + "\n"
        str_representation += \
            self.__repr__()
        str_representation += \
            + 80*"-" + "\n"

        return(str_representation)


    def __eq__(self, transcript):

        """Equality operator."""

        if isinstance(transcript, self.__class__):
            return((self.chromosome == transcript.chromosome) and 
                   (self.source == transcript.source) and            
                   (self.feature == transcript.feature) and          
                   (self.start == transcript.start) and
                   (self.end == transcript.end) and
                   (self.score == transcript.score) and              
                   (self.strand == transcript.strand) and
                   (self.frame == transcript.frame) and
                   (self.gene_id == transcript.gene_id) and
                   (self.transcript_id == transcript.transcript_id) and
                   (self.gene_name == transcript.gene_name) and
                   (self.gene_biotype == transcript.gene_biotype) and
                   (self.transcript_name == transcript.transcript_name) and
                   (self.transcript_biotype == transcript.transcript_biotype))

        return NotImplemented


    def __ne__(self, transcript):

        """Non-equality operator."""

        if isinstance(transcript, self.__class__):
            return(not self.__eq__(transcript))
        return NotImplemented


    def extend(self, transcript, verbose=False):

        """Method that extends the transcript by missing information from an existing transcript."""

        # check if we have the same genes here
        if (self != transcript):
            sys.stderr.write("ERROR: Transcript '" + self.transcript_id \
                            +"' cannot be extended by transcript '" + transcript.transcript_id \
                            +"' because the transcripts differ in annotation.\n")
            sys.exit(-1)
        
        else:
            
            # check whether we have an equal number of exons
            if len(self.exon_list_sorted_by_end_coord) != len(transcript.exon_list_sorted_by_end_coord):
                    sys.stderr.write("ERROR: Different numbers of exons " \
                                    +" were found for transcript '" \
                                    + self.transcript_id \
                                    +"'.\n")
                    sys.exit(-1)

            # check whether the transcripts are equal according to the exon coordinates
            for idx, exon in enumerate(self.exon_list_sorted_by_end_coord):

                if not exon.has_equal_coordinates(transcript.exon_list_sorted_by_end_coord[idx]):
                    sys.stderr.write("ERROR: Exon '" + exon.exon_id \
                                    +"' has different coordinates than exon '" \
                                    + transcript.exon_list_sorted_by_end_coord[idx].exon_id \
                                    +"' assigned to transcript '" 
                                    + self.transcript_id \
                                    +"'.\n")
                    sys.exit(-1)

            # if our transcript was not annotated to be coding, but the incoming transcript is,
            # let's use the annotation of the incoming transcript.
            if not (self.contains_start_codon and self.contains_stop_codon) and \
              (transcript.contains_start_codon and transcript.contains_stop_codon):

                self.exon_list_sorted_by_end_coord = transcript.exon_list_sorted_by_end_coord


    def insert_into_exon_list_sorted_by_end_coord(self, exon):
        
        """Insert exon in a sorted (by the end position) list, termed 'self.exon_list_sorted_by_end_coord'."""

        index = 0
        if len(self.exon_list_sorted_by_end_coord) == 0:
             self.exon_list_sorted_by_end_coord.insert(0, exon)
        else:
            for idx, info in enumerate(self.exon_list_sorted_by_end_coord):
                if exon.end > info.end:
                    index = idx + 1
            self.exon_list_sorted_by_end_coord.insert(index, exon)


    def get_terminal_exon(self):

        """Get last exon of the transcript"""

        if self.strand == '+':
            return(self.exon_list_sorted_by_end_coord[-1])
        elif self.strand == '-':
            return(self.exon_list_sorted_by_end_coord[0])


    def get_start_exon(self):

        """Get first exon of the transcript"""
        
        if self.strand == '+':
            return(self.exon_list_sorted_by_end_coord[0])
        elif self.strand == '-':
            return(self.exon_list_sorted_by_end_coord[-1])


    def get_intermediate_exons(self):

        """Get intermediate exons of the transcript. All exons but first and last"""

        return(self.exon_list_sorted_by_end_coord[1:-1])


    def get_all_exons_but_last(self):
        
        """Get all exons except the temrinal"""

        if self.strand == "+":
            return(self.exon_list_sorted_by_end_coord[0:-1])
        elif self.strand == "-":
            return(self.exon_list_sorted_by_end_coord[-1:0:-1])


    def contains_start_codon(self):

        """Determine if transcript has start codon"""

        for exon in self.exon_list_sorted_by_end_coord:

            if exon.start_codon is not None:

                return(True)

        return(False)


    def contains_stop_codon(self):

        """Determine if transcript has start codon"""

        for exon in self.exon_list_sorted_by_end_coord:

            if exon.stop_codon is not None:

                return(True)

        return(False)


    def contains_start_codon_for_novel_transcripts(self):

        """Determine if transcript has start codon"""

        for exon in self.novel_exons:

            if exon.start_codon is not None:

                return(True)

        return(False)


    def contains_stop_codon_for_novel_transcripts(self):

        """Determine if transcript has start codon"""

        for exon in self.novel_exons:

            if exon.stop_codon is not None:

                return(True)

        return(False)


    def get_novel_exons_sorted_by_end(self):

        """
        Return list of novel exons sorted by the
        transcript ends
        """

        return(self.sort_exons_by_end(self.novel_exons))


    def sort_exons_by_end(self, list_of_exons):

        """
        Function that gets a list of exons and sorts them by
        the exon end coordinate and returns the sorted list
        """

        return(sorted(list_of_exons, key=lambda x: x.end, reverse=False))


    def get_existing_and_upstream_exons(self, start, end, strand):
        
        """
        Generate the name for the novel transcripts based on 
        the exons found upstream of it
        """

        upstream_exons = []

        if strand == '+':

            for x in self.exon_list_sorted_by_end_coord:

                if int(x.end) <= int(end):

                    upstream_exons.append(x)

        elif strand == '-':

            for x in self.exon_list_sorted_by_end_coord[::-1]:

                if int(start) <= int(x.start):

                    upstream_exons.append(x)

        return(upstream_exons)


    def get_upstream_exons(self, start, end, strand):
        
        """Generate the name for the novel transcripts based on the exons found upstream of it"""

        upstream_exons = []

        if strand == '+':

            for x in self.exon_list_sorted_by_end_coord:

                if int(x.end) < int(end):

                    upstream_exons.append(x)

        elif strand == '-':

            for x in self.exon_list_sorted_by_end_coord[::-1]:

                if int(start) < int(x.start):

                    upstream_exons.append(x)

        return(upstream_exons)


    def write_transcript_gtf(self, output_file):

        """Function that writes in GTF format the gene/transcript/exon annotation (known or/and novel)"""

        extra_field  = "gene_id \""+self.gene_id+"\"; "
        extra_field += "transcript_id \""+self.transcript_id+"\"; "
        extra_field += "gene_name \""+self.gene_name+"\"; "
        extra_field += "gene_biotype \""+self.gene_biotype+"\"; "
        extra_field += "transcript_name \""+self.transcript_name+"\"; "
        extra_field += "transcript_biotype \""+self.transcript_biotype+"\"; \n"

        output_file.write("\t".join([self.chromosome, 
                                    self.source, 
                                    self.feature,
                                    str(self.start+1),
                                    str(self.end),
                                    self.score,
                                    self.strand,
                                    self.frame,
                                    extra_field]))

                                    # "gene_id \""+self.gene_id+"\"; transcript_id \""+self.transcript_id+"\"; gene_name \""+self.gene_name+"\"; gene_source \""+self.gene_source+"\"; gene_biotype \""+self.gene_biotype+"\"; transcript_name \""+self.transcript_name+"\"; transcript_source \""+self.transcript_source+"\";\n"]))


    def write_exons_gtf(self, output_file, with_CDS):
        
        """Write exons of the transcript in gtf format"""

        # Novel exons that have been found in this run
        # (therefore there exist self.novel_exons entries)
        if ("TECtool" in self.source) and (len(self.novel_exons) > 0):

            if self.write_CDS is True:
                for exon in self.novel_exons:
                    if with_CDS is True:
                        exon.write_exon_gtf(output_file, with_CDS=True)
                    else:
                        exon.write_exon_gtf(output_file, with_CDS=False)
            elif self.write_CDS is False:
                for exon in self.novel_exons:
                    exon.write_exon_gtf(output_file, with_CDS=False)

        # Case normal annotation
        else:
            if self.strand == "+":
                for exon in self.exon_list_sorted_by_end_coord:
                    if with_CDS is True:
                        exon.write_exon_gtf(output_file, with_CDS=True)
                    else:
                        exon.write_exon_gtf(output_file, with_CDS=False)
            elif self.strand == "-":
                for exon in self.exon_list_sorted_by_end_coord[::-1]:
                    if with_CDS is True:
                        exon.write_exon_gtf(output_file, with_CDS=True)
                    else:
                        exon.write_exon_gtf(output_file, with_CDS=False)
            else:
                sys.stderr.write("ERROR: '%s' is an unknown strand!\n" % (self.strand))
                sys.exit(-1)