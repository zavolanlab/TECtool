import os

class Exon:

    """
    This class represents an Exon

        :param chromosome: the chromosome on which the exon is located on.
        :param source: the source of the exon annotation
            (as provided in the GTF/GFF file).
        :param start: the start of the exon
        :param end: the end of the exon
        :param strand: the end of the exon
        :param transcript_id: unique id for the transcript to which the
            exon belongs to (e.g. Ensembl Transcript ID).
        :param gene_id: unique id for the gene to which the exon belongs
            to (e.g. Ensembl Gene ID).
        :param exon_id: unique id for the exon (e.g. Ensembl Exon ID).
        :rtype: Exon object

    *Class members*


        *chromosome*
            String. The chromosome at which the exon is located on.

        *source*
            String. The name of the program that generated the annotation
                feature, or the data source (database or project name).

        *feature*
            String. feature type name, e.g. Gene, Transcript, Exon

        *start*
            Int. Start position of the exon.

        *end*
            Int. End position of the exon.

        *score*
            String. Score.

        *strand*
            String. Strand of the exon (defined as + (forward) or - (reverse)).

        *frame*
            String. It can be "."" which is no frame info or 0/1/2.
             '0' indicates that the first base of the feature is
             the first base of a codon, '1' that the second base is
             the first base of a codon, and so on.

        *gene_id*
            String. Unique id for the gene to which the exon belongs to (e.g.
                Ensembl Gene ID).

        *transcript_id*
            String. Unique id for the transcript to which the exon belongs to
                (e.g. Ensembl Transcript ID).

        *exon_number*
            String.

        *gene_name*
            String. The gene symbol (e.g. FAM188A).

        *gene_biotype*
            String. The biotype of the exon (e.g. protein_coding).

        *transcript_name*
            String.

        *transcript_biotype*
            String. The biotype of the exon (e.g. protein_coding).

        *exon_id*
            String. Unique id for the exon (e.g. Ensembl Exon ID).


        *CDS*
            HTSeq Genomic interval. The CDS region. None in case of no CDS.

        *frame*
            Int. Frame of the CDS region. None in case of no CDS.

        *start_codon*
            HTSeq Genomic interval. Start codon region.
                None in case of no start codon.

        *stop_codon*
            HTSeq Genomic interval. Stop codon region.
                None in case of no stop codon.
    """

    def __init__(
        self,
        chromosome,
        source,
        feature,
        start,
        end,
        score,
        strand,
        frame,
        gene_id,
        transcript_id,
        exon_number,
        gene_name,
        gene_biotype,
        transcript_name,
        transcript_biotype,
        exon_id
    ):
        """
        Initialise Transcript
        """

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
        self.exon_number = str(exon_number)
        self.gene_name = str(gene_name)
        self.gene_biotype = str(gene_biotype)
        self.transcript_name = str(transcript_name)
        self.transcript_biotype = str(transcript_biotype)
        self.exon_id = str(exon_id)

        # extra fields important for GTF file
        self.CDS = None
        # self.frame = None
        self.start_codon = None
        self.stop_codon = None

    def __repr__(self):
        """
        Create a string representation of a Exon object.
        """

        repr_string = ""
        for attr in vars(self):
            repr_string += \
                "{}\tobj.{} = {} {}".format(
                    type(getattr(self, attr)),
                    attr,
                    getattr(self, attr),
                    os.linesep
                )

        return(repr_string)

    def __str__(self):
        """
        Create a readable string representation of a Exon object.
        """

        str_representation = \
            + 80 * "_" + os.linesep \
            + 80 * "-" + os.linesep
        str_representation += \
            ("Exon(object):\t" + self.exon_id + "\n")
        str_representation += \
            + 80 * "-" + os.linesep
        str_representation += \
            self.__repr__()
        str_representation += \
            + 80 * "-" + os.linesep

        return(str_representation)

    def has_equal_coordinates(self, exon):
        """
        Method that checks whether the exon has equal
        coordinates to a given exon.
        """

        if isinstance(exon, self.__class__):
            return((self.chromosome == exon.chromosome) and
                   (self.start == exon.start) and
                   (self.end == exon.end) and
                   (self.strand == exon.strand))

        return NotImplemented

    def to_bed(self):
        """
        Return the exon in bed format as string
        """

        return("\t".join([self.chromosome,
                          str(self.start),
                          str(self.end),
                          self.gene_id,
                          "0",
                          self.strand + "\n"]))

    def get_bed_line_with_exonid_geneid_as_list(self):
        """
        Return the exon in bed format as list
        """

        return([self.chromosome,
                str(self.start),
                str(self.end),
                "$".join([self.exon_id, self.gene_id]),
                "0",
                self.strand])

    def get_bed_line_with_geneid_as_list(self):
        """
        Return the exon in bed format as list
        """

        return([self.chromosome,
                str(self.start),
                str(self.end),
                self.gene_id,
                "0",
                self.strand])

    def write_exon_gtf(self, output_file, with_CDS):
        """
        Write exon in gtf format.
        If with_CDS is specified then the start codons,
        CDSs, and stop codons are also written.
        """

        extra_field = "gene_id \"" + self.gene_id + "\"; "
        extra_field += "transcript_id \"" + self.transcript_id + "\"; "
        extra_field += "exon_number \"" + self.exon_number + "\"; "
        extra_field += "gene_name \"" + self.gene_name + "\"; "
        extra_field += "gene_biotype \"" + self.gene_biotype + "\"; "
        extra_field += "transcript_name \"" + self.transcript_name + "\"; "
        extra_field += \
            "transcript_biotype \"" + self.transcript_biotype + "\"; "
        extra_field += "exon_id \"" + self.exon_id + "\"; \n"

        output_file.write("\t".join([
            self.chromosome,
            self.source,
            self.feature,
            str(self.start + 1),
            str(self.end),
            self.score,
            self.strand,
            '.',
            extra_field])
        )

        # Flag to allow or not to write the CDS for the specific transctipt ...
        if with_CDS:

            if self.start_codon is not None:

                output_file.write("\t".join([
                    self.chromosome,
                    self.source,
                    "start_codon",
                    str(self.start_codon.start + 1),
                    str(self.start_codon.end),
                    self.score,
                    self.strand,
                    str(self.frame),
                    extra_field])
                )

            if self.CDS is not None:

                output_file.write("\t".join([
                    self.chromosome,
                    self.source,
                    "CDS",
                    str(self.CDS.start + 1),
                    str(self.CDS.end),
                    self.score,
                    self.strand,
                    str(self.frame),
                    extra_field])
                )

            if self.stop_codon is not None:

                output_file.write("\t".join([
                    self.chromosome,
                    self.source,
                    "stop_codon",
                    str(self.stop_codon.start + 1),
                    str(self.stop_codon.end),
                    self.score,
                    self.strand,
                    str(self.frame),
                    extra_field])
                )

all = ('Exon', )
