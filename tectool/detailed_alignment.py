class DetailedAlignment:

    """ This class represents a detailed alignment """

    def __init__(self, aln):

        self.aln = aln # this will be dropped in the feature, so that memory does not blow up
        self.number_of_S = 0
        self.split_event_list = list()
        self.regions_set = set()
        self.spans_regions_boarder = False
