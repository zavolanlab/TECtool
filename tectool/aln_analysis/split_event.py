class SplitEvent:

    """ This class represents a split event """

    def __init__(self, chrom, strand, five_prime_ss):
        self.chrom = chrom
        self.strand = strand
        self.five_prime_ss = five_prime_ss
        self.three_prime_ss = None
