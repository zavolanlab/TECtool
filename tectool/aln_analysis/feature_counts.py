class FeatureCounts(object):
#class FeatureCounts:

    """
    :param region: The region of the feature.
    :param annotation: Annotation of the feature. Free text.
    :param gene_id: The gene_id that the region belongs to
    :param splice_in_all: Number of splice reads
                          that splice in the border region.
    :param splice_in_borders: Number of splice reads that splice in the region.
    :param splice_out_all: Number of splice reads that splice
                           out the border the region.
    :param splice_out_borders: Number of splice reads that splice
                           out of the region.
    :param unspliced_feature: Number of unspliced reads that fall
                           within the region.
    :param unspliced_5pSS: Number of unspliced reads that fall in
                           the 5' splice site.
    :param unspliced_3pSS: Number of unspliced reads that fall
                          in the 3' splice site.
    :param profile: Per base profile of the region.
    :rtype: FeatureCounts object
    *Class members*
    *region*
    *annotation
    *gene_id
    *splice_in_all*
    *splice_in_borders*
    *splice_out_all*
    *splice_out_borders*
    *unspliced_feature*
    *unspliced_5pSS*
    *unspliced_3pSS*
    *profile*
    """

    # define class/static variables
    region_colname = 'region'
    annotation_colname = 'annotation'
    gene_id_colname = 'gene_id'
    splice_in_all_colname = 'splice_in_all'
    splice_in_borders_colname = 'splice_in_borders'
    splice_out_all_colname = 'splice_out_all'
    splice_out_borders_colname = 'splice_out_borders'
    unspliced_feature_colname = 'unspliced_feature'
    unspliced_5pSS_colname = 'unspliced_5pSS'
    unspliced_3pSS_colname = 'unspliced_3pSS'
    profile_colname = 'profile'
    total_reads_colname = 'total_reads'
    union_exon_length_colname = 'union_exon_length'
    GeneExpressionPerKBApproximated_colname = 'GeneExpressionPerKBApproximated'

    def __init__(
        self,
        region,
        annotation,
        gene_id,
        splice_in_all,
        splice_in_borders,
        splice_out_all,
        splice_out_borders,
        unspliced_feature,
        unspliced_5pSS,
        unspliced_3pSS,
        profile
    ):
        """Constructor"""

        # initialize variables
        self.region = region
        self.annotation = annotation
        self.gene_id = gene_id
        self.splice_in_all = splice_in_all
        self.splice_in_borders = splice_in_borders
        self.splice_out_all = splice_out_all
        self.splice_out_borders = splice_out_borders
        self.unspliced_feature = unspliced_feature
        self.unspliced_5pSS = unspliced_5pSS
        self.unspliced_3pSS = unspliced_3pSS
        self.profile = profile

        self.total_reads = None
        self.union_exon_length = None
        self.GeneExpressionPerKBApproximated = None

    def to_dict(
        self
    ):
        """
        Return dictionary of the object
        """
        return {
            self.region_colname: self.region,
            self.annotation_colname: self.annotation,
            self.gene_id_colname: self.gene_id,
            self.splice_in_all_colname: self.splice_in_all,
            self.splice_in_borders_colname: self.splice_in_borders,
            self.splice_out_all_colname: self.splice_out_all,
            self.splice_out_borders_colname: self.splice_out_borders,
            self.unspliced_feature_colname: self.unspliced_feature,
            self.unspliced_5pSS_colname: self.unspliced_5pSS,
            self.unspliced_3pSS_colname: self.unspliced_3pSS,
            self.profile_colname: self.profile,
            self.total_reads_colname: self.total_reads,
            self.union_exon_length_colname: self.union_exon_length,
            self.GeneExpressionPerKBApproximated_colname:
                self.GeneExpressionPerKBApproximated
        }


