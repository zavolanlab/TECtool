# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

try:
  import HTSeq
except(Exception):
  raise("[ERROR] HTSeq was not imported properly. Exiting.")
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
  import os
except(Exception):
  raise("[ERROR] os was not imported properly. Exiting.")
  sys.exit(-1)

try:
    from collections import defaultdict
except(Exception):
    raise("[ERROR] collections dedaultdict was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from itertools import combinations
    from collections import defaultdict
except(Exception):
    raise("[ERROR] collections dedaultdict was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import numpy as np
except(Exception):
    raise("[ERROR] numpy was not imported properly")
    sys.exit(-1)

try:
    import matplotlib.mlab as mlab
except(Exception):
    raise("[ERROR] matplotlib.mlab was not imported properly")
    sys.exit(-1)

try:
    import matplotlib.pyplot as plt
except(Exception):
    raise("[ERROR] matplotlib.pyplot was not imported properly")
    sys.exit(-1)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# PolyASite Class
# -----------------------------------------------------------------------------
class PolyASite():

    """PolyASite Class"""

    def __init__(self, chromosome, start, end, strand, intron_id):

        """Initialise PolyASite"""

        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.intron_id = intron_id

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():

    """ Main function """

    __doc__ = "Create histogram of the distances of the intronic poly(A) sites"

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)


    parser.add_argument("--bed",
                      dest = "bed",
                      help = "BED (bed12) file to read [REQUIRED]",
                      metavar = "FILE")

    parser.add_argument("--out",
                      dest = "out",
                      help = "Output directory [REQUIRED]",
                      metavar="FILE")

    parser.add_argument("-v",
                      "--verbose",
                      action = "store_true",
                      dest = "verbose",
                      default = False,
                      help = "Verbose")

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

    if options.verbose:
        sys.stdout.write("Create histogram of the distances of the intronic poly(A) sites")

    dict_of_introns = defaultdict(list)
    dict_of_observed_polya_sites = dict()

    # parse bed file
    with open(options.bed) as bed:

        for bed_line in bed:

            bed_line_splitted = bed_line.strip().split("\t")

            polya_chr     = bed_line_splitted[0]
            polya_start   = bed_line_splitted[1]
            polya_end     = bed_line_splitted[2]
            polya_strand  = bed_line_splitted[5]

            intron_chr    = bed_line_splitted[6]
            intron_start  = bed_line_splitted[7]
            intron_end    = bed_line_splitted[8]
            intron_strand = bed_line_splitted[11]
            intron_gene   = bed_line_splitted[9]

            # generate intronic id
            intron_id = ":".join([intron_chr, intron_start, intron_end, intron_strand])
            # generate intronic poly(A) site id
            intronic_polya_site_id = ":".join([polya_chr, polya_start, polya_end, polya_strand])

            # if this is an intonic polya site that was not used before
            if intronic_polya_site_id not in dict_of_observed_polya_sites:

                # add it to the dicitionary
                dict_of_observed_polya_sites[intronic_polya_site_id] = intronic_polya_site_id

                # create intronic polya site object
                intronic_polya_site = PolyASite(polya_chr, polya_start, polya_end, polya_strand, intron_id)

                # and append it to the dictionary of introns
                dict_of_introns[intron_id].append(intronic_polya_site)

    # calculate distances
    distances = []
    for intron_id in dict_of_introns:
        for first_intronic_polya_site, second_intronic_polya_site in combinations(dict_of_introns[intron_id], 2):
            distance = -1
            if (first_intronic_polya_site.strand == "+") and (second_intronic_polya_site.strand == "+"):
                distance = abs(first_intronic_polya_site.start - second_intronic_polya_site.start)
                distances.append(distance)
            elif (first_intronic_polya_site.strand == "-") and (second_intronic_polya_site.strand == "-"):
                distance = abs(first_intronic_polya_site.end - second_intronic_polya_site.end)
                distances.append(distance)
#            print(distance)
#        print("--")

    array_distances = np.array(distances)

    # plot
    n, bins, patches = plt.hist(array_distances, 200) #, normed=1, facecolor='green', alpha=0.75)
    #print("Mean", np.mean(array_distances))
    #print("Median", np.median(array_distances))
    #print("Std", np.std(array_distances))
    #print("Min", np.min(array_distances))
    #print("Max", np.max(array_distances))
    plt.savefig(os.path.join(options.out, "histogram_of_intronic_polya_sites.pdf"))


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
