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
    import os
except(Exception):
    raise("[ERROR] os was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import sys
except(Exception):
    raise("[ERROR] sys was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import copy
except(Exception):
    raise("[ERROR] copy was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import csv
except(Exception):
    raise("[ERROR] csv was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import numpy as np
except(Exception):
    raise("[ERROR] numpy was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import math
except(Exception):
    raise("[ERROR] math was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import itertools
except(Exception):
    raise("[ERROR] itertools was not imported properly")
    sys.exit(-1)

try:
    from sklearn import metrics
    from sklearn.metrics import accuracy_score
    from sklearn.metrics import f1_score
    from sklearn.metrics import confusion_matrix
    from sklearn.metrics import make_scorer
    from sklearn.metrics import roc_curve, auc
    from sklearn.model_selection import GridSearchCV
    from sklearn.preprocessing import label_binarize
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn import linear_model
    from sklearn.cross_validation import train_test_split
    from sklearn import neighbors
    from sklearn.model_selection import StratifiedKFold
except(Exception):
    raise("[ERROR] sklearn was not imported properly")

try:
    from scipy import interp
    from scipy import stats
    from scipy.optimize import curve_fit
except(Exception):
    raise("[ERROR] scipy was not imoprted properly")

try:
    # for feature selection
    # http://rasbt.github.io/mlxtend/user_guide/feature_selection/SequentialFeatureSelector/
    from mlxtend.feature_selection import SequentialFeatureSelector as SFS
except(Exception):
    raise e

try:
    import random
except(Exception):
    raise("[ERROR] random was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import time
except(Exception):
    raise("[ERROR] time was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
except(Exception):
    raise "[ERROR] plt from matplotlib.pyplot was not imported properly. Exiting."
    sys.exit(-1)

try:
    import pandas as pd
except(Exception):
    raise("[ERROR] pandas was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import pybedtools
except(Exception):
    raise("[ERROR] pybedtools was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from ast import literal_eval
except(Exception):
    raise("[ERROR] literal_eval from ast was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import pickle
except(Exception):
    raise("[ERROR] pickle was not imported properly. Exiting.")
    sys.exit(-1)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import our own modules
# -----------------------------------------------------------------------------

from aln_analysis import DetailedAlignment
from aln_analysis import SplitEvent
from aln_analysis import AnalysisUnit
from machine_learning import BayesClassifier

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# CLASSES
# -----------------------------------------------------------------------------

class MachineLearningUnit(object):

    """
        Class for machine learning. A MachineLearningUnit object can collect 
        and filter training data, test data and validation data. Moreover it 
        can train classifiers, select features and perform predictions.

        :rtype: MachineLearningUnit object

    *Class members*

        *terminal_exon_training_data*
            Numpy ndarray. A matrix containing training data.

        *intermediate_exon_training_data*
            Numpy ndarray. A matrix containing training data.

        *background_training_data*
            Numpy ndarray. A matrix containing training data.

        *novel_terminal_exon_candidates_data*
            Numpy ndarray. A matrix containing novel terminal exon candidates.

        *novel_terminal_readthrough_exon_candidates_data*
            Numpy ndarray. A matrix containing novel read through exon candidates.

        *terminal_exon_gene_dict*
            Dictionary. A dictionary that has gene ids as keys.

        *intermediate_exon_gene_dict*
            Dictionary. A dictionary that has gene ids as keys.

        *background_gene_dict*
            Dictionary. A dictionary that has gene ids as keys.

        *terminal_exons_features*
            pandas.DataFrame. Contains one terminal exon data set per line, having one feature per column.

        *intermediate_exons_features*
            pandas.DataFrame. Contains one intermediate exon data set per line, having one feature per column.

        *background_regions_features*
            pandas.DataFrame. Contains one background region data set per line, having one feature per column.

        *novel_terminal_exon_features*
            pandas.DataFrame. Contains one novel terminal exon data set per line, having one feature per column.

        *features*
            List. Contains the names of the calculated features (available in the training dataframes).

        *selected_features*
            List. Contains the names of the selected features (available in the training dataframes).

        *class_col*
            String. Contains the name of the column (available in the sampled dataframes) that holds the information about the class of each data set.

        *terminal_exon_class*
            String. Contains the name of the terminal exon class.

        *intermediate_exon_class*
            String. Contains the name of the intermediate exon class.

        *background_region_class*
            String. Contains the name of the background region class.

        *region_classes*
            List. Contains the region class strings (terminal_exon_class, intermediate_exon_class, background_region_class).

        *training_df*
            pandas.DataFrame. Dataframe that contains all data sets (concatendated terminal_exons_features, intermediate_exons_features, background_regions_features).

        *validation_df*
            pandas.DataFrame. Dataframe that contains all data sets (concatendated terminal_exons_features, intermediate_exons_features, background_regions_features).

        *classifier_dict*
            Dictionary. Dictionary with classifier names as keys and 'sklearn.neighbors.classification.Classifier' objects as values.

        *classifier_funcdict*
            Dictionary. Dictionary with classifier names as keys and function pointers as values. Each function pointer can create a specific type of classifier given X and y data sets.

        *novel_terminal_exons*
            pandas.DataFrame. Candidates that have been classified as novel terminal exons.

        *novel_intermediate_exons*
            pandas.DataFrame. Candidates that have been classified as novel intermediate exons.

        *novel_background_regions*
            pandas.DataFrame. Candidates that have been classified as novel background regions.

        *selected_novel_terminal_exons*
            pandas.DataDrame. Candidates that have been classified as novel terminal exons and were selected based on the probabilities.
    """
    

    def __init__(self):
        
        # further members that are later on crucial 
        # for performing machine learning
        self.terminal_exon_training_data = None
        self.intermediate_exon_training_data = None
        self.background_training_data = None
        self.novel_terminal_exon_candidates_data = None
        self.novel_terminal_readthrough_exon_candidates_data = None

        # dataframes that are used by the machine learning approach
        self.terminal_exons_features = None
        self.intermediate_exons_features = None
        self.background_regions_features = None
        self.novel_terminal_exon_features = None

        # a list of features that are available for each training class
        self.features = None
        
        # list for the features selected
        self.selected_features = None

        # strings that contain
        self.class_col = "class"
        self.terminal_exon_class = "terminal"
        self.intermediate_exon_class = "intermediate"
        self.background_region_class = "background"

        # get a collection in a fixed order of all possible classes
        self.region_classes = \
            [self.terminal_exon_class, 
             self.intermediate_exon_class, 
             self.background_region_class]

        # the dataframe that contains training data of all classes
        self.training_df = None

        # the dataframe that contains validation data of all classes
        self.validation_df = None

        # members for the novel candidates (the intermedidate and 
        # background classified regions are useful for checking the 
        # properties of these regions).
        self.novel_terminal_exons = None
        self.novel_intermediate_exons = None
        self.novel_background_regions = None

        # selected novel terminal exons
        self.selected_novel_terminal_exons = None

        # dictionaries that contain all genes as keys
        # that are used within the according training data
        self.terminal_exon_gene_dict = None
        self.intermediate_exon_gene_dict = None
        self.background_gene_dict = None

        # dictionary for classifiers
        self.classifier_dict = dict()

        # dictionary for classifier functions
        self.classifier_funcdict = \
            {"KNeighbors" : self.create_kneighbors_classifier,
             "multiclass_SVC" : self.create_multiclass_SVC_classifier,
             "Bayes" : self.create_Bayes_classifier}

        # list of training dataframes
        self.list_of_training_dataframes = []

        # list of validation dataframes
        self.list_of_validation_dataframes = []

        # dictionary of terminal probabilities
        self.terminal_probabilities_dict = dict()

        # dictionary of intermediate probabilities
        self.intermediate_probabilities_dict = dict()

        # dictionary of background probabilities
        self.background_probabilities_dict = dict()

        # create a random seed
        random.seed(time.time())


    def write_training_df_to_file(self, 
                                  training_df_file_path,
                                  verbose=True):

        """Method that writes the training data ('training_df') to a file."""

        if verbose: 
            sys.stdout.write("Writing training data ('training_df') to file: "\
                            +"%s\n" % (training_df_file_path))

        self.training_df.to_csv(training_df_file_path, sep="\t")


    def load_training_df_from_file(self, 
                                   training_df_file_path,
                                   verbose=True):

        """Method that reads training data from a file into 'training_df'."""

        if verbose: 
            sys.stdout.write("Reading training data ('training_df') from file: "\
                            +"%s\n" % (training_df_file_path))

        self.training_df = \
            pd.read_csv(training_df_file_path, index_col=0, sep="\t")


    def write_validation_df_to_file(self, 
                                    validation_df_file_path,
                                    verbose=True):

        """Method that writes the validation data ('validation_df') to a file."""

        if verbose: 
            sys.stdout.write("Writing validation data ('validation_df') to file: "\
                            +"%s\n" % (validation_df_file_path))

        self.validation_df.to_csv(validation_df_file_path, sep="\t")


    def load_validation_df_from_file(self, 
                                     validation_df_file_path,
                                     verbose=True):

        """Method that reads validation data from a file into 'validation_df'."""

        if verbose: 
            sys.stdout.write("Reading validation data ('validation_df') from file: "\
                            +"%s\n" % (validation_df_file_path))

        self.validation_df = \
            pd.read_csv(validation_df_file_path, index_col=0, sep="\t")


    def create_terminal_exon_training_set(self, 
                                          terminal_exons_bed_file_path,
                                          sequencing_direction,
                                          max_splice_fuzziness,
                                          output_dir,
                                          genes_to_consider_dict,
                                          bam_file_path,
                                          annotation,
                                          threshold_to_filter):

        """
        Method that creates training data for terminal exons and
        updates the correspoting objects
        
        
        :param terminal_exons_bed_file_path: path to the bed file that contains
          the coordinates of the terminal exons to be considered.

        """

        # ---------------------------------------------------------------------
        # create an AnalysisUnit object for each terminal exon
        # and store it in a dictionary

        # dictionary for the terminal exons
        aunits_terminal_exons_dict = dict()

        with open(terminal_exons_bed_file_path) as annotated_terminal_exons:

            for exon in annotated_terminal_exons:

                exon_sp = exon.strip().split('\t')
                gene_id = str(exon_sp[3].split(",")[0])
                feature_annotation = str(exon_sp[3].split(",")[1])

                if gene_id in genes_to_consider_dict:

                    exon_iv = HTSeq.GenomicInterval(exon_sp[0], 
                                                    int(exon_sp[1]), 
                                                    int(exon_sp[2]), 
                                                    exon_sp[5])
    
                    if exon_iv not in aunits_terminal_exons_dict:
                        aunits_terminal_exons_dict[exon_iv] = \
                            AnalysisUnit(unit_id=exon_iv, 
                                         potential_5pSS_exons=None,
                                         gene_id = gene_id)
                        aunits_terminal_exons_dict[exon_iv].annotation = feature_annotation

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Open the BAM file
        # ---------------------------------------------------------------------
        # now count the things
        bam_file_path = bam_file_path
        bam = HTSeq.BAM_Reader(bam_file_path)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Go over all AnalysisUnit objects for terminal exons, fetch the reads
        # and count
        # ---------------------------------------------------------------------
        sys.stdout.write("Counting annotated terminal exons...\n")
  
        # w = open(terminal_exons_statistics_file_path, 'w')
        # w.write("\t".join(["Region", "Annotation", "GeneId", "SpliceInAll", "SpliceInBorders", 
        #                    "SpliceOutAll", "SpliceOutBorders", "UnsplicedExon", 
        #                    "Unspliced_5pSS", "Unspliced_3pSS", "profile\n"]))
        # w.close()
  
        # go over each unit
        unit_nr = 0
        for unit_id in aunits_terminal_exons_dict.keys():
            unit_nr += 1
  
            # give some feedback about the state of the script
            # (how many units have been analyzed so far?)
            if (unit_nr % 100) == 0:
                sys.stderr.write("Regions processed:\t" + str(unit_nr) + "\n")
  
            # get the AnalysisUnit object
            aunits_terminal_exons_dict[unit_id].analyze_reads_for_annotated_regions(bam=bam,
                                                                                    unit_id = unit_id,
                                                                                    sequencing_direction=sequencing_direction,
                                                                                    splice_fuzziness=max_splice_fuzziness,
                                                                                    count_unique_mapping_reads_only=True,
                                                                                    tmp=output_dir,
                                                                                    threshold_to_filter = threshold_to_filter,
                                                                                    feature_type = "terminal_exon",
                                                                                    annotation = annotation,
                                                                                    verbose=False)
            # free memory
            try:
                del(aunits_terminal_exons_dict[unit_id])
            except(KeyError):
                pass
    

    def create_intermediate_exon_training_set(self, 
                                              intermediate_exons_bed_file_path,
                                              sequencing_direction,
                                              max_splice_fuzziness,
                                              output_dir,
                                              genes_to_consider_dict,
                                              bam_file_path,
                                              annotation,
                                              threshold_to_filter):

        """
        Method that creates training data for intermediate exons and
        updates the correspoting objects

        :param intermediate_exons_bed_file_path: path to the bed file that contains
          the coordinates of the intermediate exons to be considered.
        """

        # ---------------------------------------------------------------------
        # create an AnalysisUnit object for each intermediate exon
        # and store it in a dictionary

        # dictionary for the intermediate exons
        aunits_intermediate_exons_dict = dict()

        with open(intermediate_exons_bed_file_path) as annotated_intermediate_exons:

            for exon in annotated_intermediate_exons:

                exon_sp = exon.strip().split('\t')
                gene_id = str(exon_sp[3].split(",")[0])
                feature_annotation = str(exon_sp[3].split(",")[1])

                if gene_id in genes_to_consider_dict:
    
                    exon_iv = HTSeq.GenomicInterval(exon_sp[0], 
                                                    int(exon_sp[1]),
                                                    int(exon_sp[2]), 
                                                    exon_sp[5])
    
                    if exon_iv not in aunits_intermediate_exons_dict:
                        aunits_intermediate_exons_dict[exon_iv] = \
                            AnalysisUnit(unit_id=exon_iv, 
                                         potential_5pSS_exons=None,
                                         gene_id = gene_id)
                        aunits_intermediate_exons_dict[exon_iv].annotation = feature_annotation

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Open the BAM file
        # -------------------------------------------------------------------------
        # now count the things
        bam_file_path = bam_file_path
        bam = HTSeq.BAM_Reader(bam_file_path)

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Go over all AnalysisUnit objects for intermediate exons, fetch the reads
        # and count
        # -------------------------------------------------------------------------
        sys.stdout.write("Counting annotated intermediate exons...\n")
  
        # # Write the file
        # w = open(intermediate_exons_statistics_file_path, 'w')
        # w.write("\t".join(["Region", "Annotation", "GeneId", "SpliceInAll", "SpliceInBorders", 
        #                    "SpliceOutAll", "SpliceOutBorders", "UnsplicedExon", 
        #                    "Unspliced_5pSS", "Unspliced_3pSS", "profile\n"]))
        # w.close()
  
        # go over each unit
        unit_nr = 0
        for unit_id in aunits_intermediate_exons_dict.keys():
            unit_nr += 1
  
            # give some feedback about the state of the script
            # (how many units have been analyzed so far?)
            if (unit_nr % 100) == 0:
                sys.stderr.write("Regions processed:\t" + str(unit_nr) + "\n")
  
            # get the AnalysisUnit object
            aunits_intermediate_exons_dict[unit_id].analyze_reads_for_annotated_regions(bam = bam,
                                                                                        unit_id = unit_id,
                                                                                        sequencing_direction = sequencing_direction,
                                                                                        splice_fuzziness = max_splice_fuzziness,
                                                                                        count_unique_mapping_reads_only = True,
                                                                                        tmp = output_dir,
                                                                                        threshold_to_filter = threshold_to_filter,
                                                                                        feature_type = "intermediate_exon",
                                                                                        annotation = annotation,
                                                                                        verbose = False)
            # free memory
            try:
                del(aunits_intermediate_exons_dict[unit_id])
            except(KeyError):
                pass


    def create_background_training_set(self, 
                                       background_bed_file_path,
                                       sequencing_direction,
                                       max_splice_fuzziness,
                                       output_dir,
                                       genes_to_consider_dict,
                                       bam_file_path,
                                       annotation,
                                       threshold_to_filter=5):

        """
        Method that creates training data for background regions and
        updates the correspoting objects
        
        :param terminal_exons_bed_file_path: path to the bed file that contains
          the coordinates of the terminal exons to be considered.

        """

        # _________________________________________________________________________
        # ---------------------------------------------------------------------
        # create an AnalysisUnit object for each terminal exon
        # and store it in a dictionary
        # ---------------------------------------------------------------------

        # dictionary for the background regions
        aunits_background_dict = dict()

        with open(background_bed_file_path) as annotated_background:

            for background in annotated_background:

                background_sp = background.strip().split('\t')
                gene_id = str(background_sp[3].split(",")[0])
                feature_annotation = str(background_sp[3].split(",")[1])

                if gene_id in genes_to_consider_dict:

                    # potential 
                    background_iv = HTSeq.GenomicInterval(background_sp[0], 
                                                    int(background_sp[1]), 
                                                    int(background_sp[2]), 
                                                    background_sp[5])

                    if background_iv not in aunits_background_dict:
                        aunits_background_dict[background_iv] = \
                            AnalysisUnit(unit_id=background_iv, 
                                         potential_5pSS_exons=None,
                                         gene_id = gene_id)
                        aunits_background_dict[background_iv].annotation = feature_annotation

                    # # Make sure that this is a terminal exon that was previously selected
                    # for terminal_exon in annotation.genes[gene_id].annotated_terminal_exons:

                    #     terminal_exon_sp = terminal_exon.region.split(":")
                    #     terminal_exon_chr =  terminal_exon_sp[0]
                    #     terminal_exon_start = int(terminal_exon_sp[1])
                    #     terminal_exon_end = int(terminal_exon_sp[2])
                    #     terminal_exon_strand = terminal_exon_sp[3]

                    #     if terminal_exon_strand == "+":

                    #         if  (terminal_exon_chr == background_sp[0]) \
                    #         and (int(terminal_exon_start) == (int(background_sp[1])+1)) \
                    #         and (terminal_exon_strand == background_sp[5]):
        
                    #             if background_iv not in aunits_background_dict:
                    #                 aunits_background_dict[background_iv] = \
                    #                     AnalysisUnit(unit_id=background_iv, 
                    #                                  potential_5pSS_exons=None,
                    #                                  gene_id = gene_id)
                    #                 aunits_background_dict[background_iv].annotation = feature_annotation

                    #     elif terminal_exon_strand == "-":

                    #         if  (terminal_exon_chr == background_sp[0]) \
                    #         and (int(terminal_exon_end) == int(background_sp[2])) \
                    #         and (terminal_exon_strand == background_sp[5]):
        
                    #             if background_iv not in aunits_background_dict:
                    #                 aunits_background_dict[background_iv] = \
                    #                     AnalysisUnit(unit_id=background_iv, 
                    #                                  potential_5pSS_exons=None,
                    #                                  gene_id = gene_id)
                    #                 aunits_background_dict[background_iv].annotation = feature_annotation

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Open the BAM file
        # -------------------------------------------------------------------------
        # now count the things
        bam_file_path = bam_file_path
        bam = HTSeq.BAM_Reader(bam_file_path)

        # _________________________________________________________________________
        # -------------------------------------------------------------------------
        # Go over all AnalysisUnit objects for background, fetch the reads
        # and count
        # -------------------------------------------------------------------------
        sys.stdout.write("Counting background ...\n")

        # go over each unit
        unit_nr = 0
        for unit_id in aunits_background_dict.keys():

            unit_nr += 1
  
            # give some feedback about the state of the script
            # (how many units have been analyzed so far?)
            if (unit_nr % 100) == 0:
                sys.stderr.write("Regions processed:\t" + str(unit_nr) + "\n")
  
            # get the AnalysisUnit object
            aunits_background_dict[unit_id].analyze_reads_for_annotated_regions(bam=bam,
                                                                                unit_id = unit_id,
                                                                                sequencing_direction=sequencing_direction,
                                                                                splice_fuzziness=max_splice_fuzziness,
                                                                                count_unique_mapping_reads_only=True,
                                                                                tmp=output_dir,
                                                                                threshold_to_filter = 5,
                                                                                feature_type = feature_annotation,
                                                                                annotation = annotation,
                                                                                verbose=False)

            # free memory
            try:
                del(aunits_background_dict[unit_id])
            except(KeyError):
                pass

    def estimate_intronic_expression(self, 
                                     intron_regions_bed_file_path,
                                     sequencing_direction,
                                     max_splice_fuzziness,
                                     output_dir,
                                     genes_to_consider_dict,
                                     bam_file_path,
                                     annotation,
                                     threshold_to_filter):

        """
        Method that estimates the intronic expression and
        updates the correspoting objects

        """

        # ---------------------------------------------------------------------
        # create an AnalysisUnit object for each intron
        # and store it in a dictionary

        # dictionary for the introns exons
        aunits_introns_dict = dict()

        with open(intron_regions_bed_file_path) as annotated_introns:

            for intron in annotated_introns:

                intron_sp = intron.strip().split('\t')
                gene_id = str(intron_sp[3].split(",")[0])
                feature_annotation = str(intron_sp[3].split(",")[1])

                if gene_id in genes_to_consider_dict:

                    intron_iv = HTSeq.GenomicInterval(intron_sp[0], 
                                                    int(intron_sp[1]), 
                                                    int(intron_sp[2]), 
                                                    intron_sp[5])
    
                    if intron_iv not in aunits_introns_dict:
                        aunits_introns_dict[intron_iv] = \
                            AnalysisUnit(unit_id=intron_iv, 
                                         potential_5pSS_exons=None,
                                         gene_id = gene_id)
                        aunits_introns_dict[intron_iv].annotation = feature_annotation

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Open the BAM file
        # ---------------------------------------------------------------------
        # now count the things
        bam_file_path = bam_file_path
        bam = HTSeq.BAM_Reader(bam_file_path)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Go over all AnalysisUnit objects for terminal exons, fetch the reads
        # and count
        # ---------------------------------------------------------------------
        sys.stdout.write("Counting annotated introns (that do not contain intronic polya sites)...\n")
  
        # w = open(terminal_exons_statistics_file_path, 'w')
        # w.write("\t".join(["Region", "Annotation", "GeneId", "SpliceInAll", "SpliceInBorders", 
        #                    "SpliceOutAll", "SpliceOutBorders", "UnsplicedExon", 
        #                    "Unspliced_5pSS", "Unspliced_3pSS", "profile\n"]))
        # w.close()
  
        # go over each unit
        unit_nr = 0
        for unit_id in aunits_introns_dict.keys():
            unit_nr += 1
  
            # give some feedback about the state of the script
            # (how many units have been analyzed so far?)
            if (unit_nr % 100) == 0:
                sys.stderr.write("Regions processed:\t" + str(unit_nr) + "\n")
  
            # get the AnalysisUnit object
            aunits_introns_dict[unit_id].analyze_reads_for_annotated_regions(bam=bam,
                                                                            unit_id = unit_id,
                                                                            sequencing_direction=sequencing_direction,
                                                                            splice_fuzziness=max_splice_fuzziness,
                                                                            count_unique_mapping_reads_only=True,
                                                                            tmp=output_dir,
                                                                            threshold_to_filter = threshold_to_filter,
                                                                            feature_type = "intron",
                                                                            annotation = annotation,
                                                                            verbose=False)

            # free memory
            try:
                del(aunits_introns_dict[unit_id])
            except(KeyError):
                pass


    def estimate_gene_expression(self,
                                annotation,
                                genes_to_consider_dict,
                                bam_file_path,
                                sequencing_direction,
                                count_unique_mapping_reads_only = True,
                                verbose = False):

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # create an AnalysisUnit object for each gene region
        # and store it in a dictionary
        # ---------------------------------------------------------------------

        # dictionary for the genes
        aunits_genes_dict = dict()

        for gene_id in annotation.genes:

            # filter genes that do not have expressed terminal exons
            if gene_id in genes_to_consider_dict:
                
                chromosome = annotation.genes[gene_id].chromosome
                start = int(annotation.genes[gene_id].start)
                end = int(annotation.genes[gene_id].end)
                strand = annotation.genes[gene_id].strand

                gene_iv = HTSeq.GenomicInterval(chromosome, start, end, strand)

                if gene_iv not in aunits_genes_dict:

                    aunits_genes_dict[gene_iv] = \
                        AnalysisUnit(unit_id = gene_iv,
                        potential_5pSS_exons = None,
                        gene_id = gene_id)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Open the BAM file
        # ---------------------------------------------------------------------
        bam_file_path = bam_file_path
        bam = HTSeq.BAM_Reader(bam_file_path)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Go over all AnalysisUnit objects for non overlapping genes, fetch
        # the reads and count
        # ---------------------------------------------------------------------

        sys.stdout.write("Estimating gene expression ...\n")

        # go over each unit
        unit_nr = 0
        for unit_id in aunits_genes_dict.keys():

            unit_nr += 1
  
            # give some feedback about the state of the script
            # (how many units have been analyzed so far?)
            if (unit_nr % 100) == 0:
                sys.stderr.write("Regions processed:\t" + str(unit_nr) + "\n")

            gene_reads, gene_id = aunits_genes_dict[unit_id].estimate_gene_expression(bam = bam,
                                                                                      unit_id = unit_id,
                                                                                      sequencing_direction=sequencing_direction,
                                                                                      count_unique_mapping_reads_only=True,
                                                                                      annotation = annotation,
                                                                                      verbose=False)
            annotation.genes[gene_id].total_reads = int(gene_reads)
            annotation.genes[gene_id].estimate_ExpressionPerKBApproximated()            
            annotation.genes[gene_id].estimate_BackgroundPerKB()
            annotation.genes[gene_id].estimate_GeneExpressionPerKBInclBackground()
            annotation.genes[gene_id].estimate_GeneExpressionBackgroundPerKB()
            annotation.genes[gene_id].estimate_GeneExpressionPerKBwithoutBackground()
            annotation.genes[gene_id].estimate_BackgroundFraction()


    def create_training_dataframes(self, annotation, intermediate_output_file, terminal_output_file, background_output_file):

        """ Create training pandas dataframes """

        terminal_exon_training_data_list = []
        intermediate_exon_training_data_list = []
        background_training_data_list = []

        for gene_id in annotation.genes:

            # if (annotation.genes[gene_id].intron_length > 0) and \
            #   (annotation.genes[gene_id].union_exon_length > 0) and \
            #   (annotation.genes[gene_id].GeneExpressionBackgroundPerKB >= 0) and \
            #   (annotation.genes[gene_id].GeneExpressionPerKBInclBackground > 0) and \
            #   (annotation.genes[gene_id].GeneExpressionPerKBwithoutBackground >0) and \
            #   (annotation.genes[gene_id].overlaps_with_other_gene == False):

            if (annotation.genes[gene_id].union_exon_length > 0) and \
               (annotation.genes[gene_id].total_reads >0) and \
               (annotation.genes[gene_id].overlaps_with_other_gene == False):

                # terminal exons
                tmp_terminal_exon_training_data_list = annotation.genes[gene_id].get_annotated_terminal_exons()

                for terminal_exon_stats in tmp_terminal_exon_training_data_list:

                    terminal_exon_training_data_list.append(terminal_exon_stats)

                # intermediate exons
                tmp_intermediate_exon_training_data_list = annotation.genes[gene_id].get_annotated_intermediate_exons()

                for intermediate_exon_stats in tmp_intermediate_exon_training_data_list:

                    intermediate_exon_training_data_list.append(intermediate_exon_stats)

                # background
                tmp_background_training_data_list = annotation.genes[gene_id].get_background()

                for background_stats in tmp_background_training_data_list:

                    background_training_data_list.append(background_stats)

        # List to pandas

        labels = ["Region", "Annotation","GeneId","SpliceInAll", "SpliceInBorders", "SpliceOutAll", "SpliceOutBorders", "UnsplicedExon", "Unspliced_5pSS", "Unspliced_3pSS", "profile", "TotalGeneReads", "UnionExonLength", "GeneExpressionPerKBApproximated"] #, "IntronicReads", "IntronLength", "BackgroundPerKB", "GeneExpressionPerKBInclBackground", "GeneExpressionBackgroundPerKB", "GeneExpressionPerKBwithoutBackground", "BackgroundFraction"]

        self.terminal_exon_training_data = pd.DataFrame.from_records(terminal_exon_training_data_list, columns=labels)
        self.terminal_exon_training_data.set_index("Region", inplace=True)
        self.terminal_exon_training_data.to_csv(terminal_output_file, sep="\t", index=True)

        self.intermediate_exon_training_data = pd.DataFrame.from_records(intermediate_exon_training_data_list, columns=labels)
        self.intermediate_exon_training_data.set_index("Region", inplace=True)
        self.intermediate_exon_training_data.to_csv(intermediate_output_file, sep="\t", index=True)

        self.background_training_data = pd.DataFrame.from_records(background_training_data_list, columns=labels)
        self.background_training_data.set_index("Region", inplace=True)
        self.background_training_data.to_csv(background_output_file, sep="\t", index=True)


    def fit_linear_model_to_profile(self, profile):
        
        """
        Method that fits a linear model to a given profile.
        """

        # create a list with the pins for the profile so that we have always 
        # length 0-1 for each region on the x-axis
        x = np.arange( 0.0, 1.0, 1.0/len(profile) )

        # create the y-axis (=profile)
        y = np.array(profile)

        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        except ValueError:
            slope = np.nan
            intercept = np.nan
            r_value = np.nan

        return([slope, intercept, r_value])

    
    def inverse_cdf_profile(self, profile_probability_norm, quantiles=[0.25,0.5,0.75]):
        
        """Method that creates a profile length normalized inverse cdf profile for given quantiles."""

        # sort the quantiles
        quantiles = sorted(quantiles)
        
        # the output
        quantile_positions = []
        
        # cumulative 
        cum = 0.0
        
        # get the length of the profile
        profile_length = len(profile_probability_norm)

        # go over all nucleotides of the normalized profile
        # and sum up the profile probabilities until we reach the
        # next quantile in quantiles. Then store the nucleotide position
        # and drop the quantil (we reached it already).
        for nt_position in range(profile_length) :
            
            cum += profile_probability_norm[nt_position]
            
            # whenever we reach a quantile, we should store the 
            # (absolute) position and remove the quantile.
            while cum >= quantiles[0]:
                quantile_positions.append(nt_position)
                quantiles.pop(0)
                if len(quantiles) == 0:
                    break

            # when the last quantile was reached, break.
            if len(quantiles) == 0:
                break
        
        # length normalize the inverse cdf profile
        inverse_cdf_profile_norm = [ float(p)/profile_length for p in quantile_positions ]

        return(inverse_cdf_profile_norm)


    def cumulative_fit_funct(self, x, a, b):
        
        """
        Method that returns a second degree polynomial using x and x^2.
        """

        # How to come up with the following equation:
        # 

        ##return ( a * x + (1-a) * pow(x,2) )
        return ( 1.0/(a+b) * (a*x + b*x**2) )


    def fit_polynomial_model_to_profile(self, 
                                        x_axis, 
                                        y_axis, 
                                        diagnostic_plots_dir_path=None,
                                        region_id=None):
        
        """
        Method that fits a polynomial model of second degree to given data
        und creates diagnostic plots to the given directory ('diagnostic_plots_dir_path').
        """
        
        try:
            # fit the parameter
            popt, pcov = curve_fit( self.cumulative_fit_funct, x_axis, y_axis )
            
            # get the estimated y values
            est_y = [ self.cumulative_fit_funct( l, *popt ) for l in x_axis ]
            
            # get the fitted a parameter
            fitted_a_param = popt[0]
            fitted_b_param = popt[1]

            # normalize the parameters
            fitted_a_param_norm = fitted_a_param / (fitted_a_param+fitted_b_param)
            fitted_b_param_norm = fitted_b_param / (fitted_a_param+fitted_b_param)

            # calculate the R squared for the fit
            # https://en.wikipedia.org/wiki/Coefficient_of_determination
            SS_res = 0.0
            for index in range(len(y_axis)):
                SS_res += (y_axis[index] - est_y[index])**2
            
            SS_tot = 0.0
            y_mean = sum(y_axis) / len(y_axis)
            for index in range(len(y_axis)):
                SS_tot += (y_axis[index] - y_mean)**2

            R_squared = 1 - (SS_res / SS_tot)

            # create some plots in order to check what is going on
            if diagnostic_plots_dir_path is not None \
              and region_id is not None \
              and (R_squared > 0.999 or fitted_a_param < 0.0):
                
                file_name = "CDF_"+str(region_id)+".png"
                file_path = os.path.join(diagnostic_plots_dir_path, file_name)

                plt.figure()
                #plt.scatter( x_axis, y_axis , color='blue', alpha=0.5)
                #plt.scatter( x_axis, est_y, color='red', alpha=0.5 )
                plt.plot( x_axis, y_axis, 'ob', alpha=0.5)
                plt.plot( x_axis, est_y, 'or-', alpha=0.5)
                plt.xlim((-0.01,1.0))
                plt.ylim((-0.01,1.0))
                plt.xlabel('Length normalized position')
                plt.ylabel('Cumulative read density (distribution).')
                
                title = "R2=" + str(R_squared) + "; a=" + str(fitted_a_param_norm) + "; b=" + str(fitted_b_param_norm)
                plt.title(title)
                plt.savefig(file_path)
                plt.close('all')

        except ValueError:
            fitted_a_param_norm = np.nan
            fitted_b_param_norm = np.nan
            R_squared = np.nan

        # make a final quality check in order to exclude profiles
        # that make no sense.

        return([fitted_a_param_norm, fitted_b_param_norm, R_squared])


    def calculate_features(self, 
                           row, 
                           profile_col = 'profile', 
                           gene_expr_col = 'GeneExpressionPerKBApproximated', 
                           sep = ',', 
                           results_dir_path = None):

        """
        Method that calculates features for a given training dataframe.
        """

        # create results directory 
        if results_dir_path is not None and not os.path.exists(results_dir_path):
            os.makedirs(results_dir_path)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # PROFILE DEPENDENT FEATURES
        # ---------------------------------------------------------------------
        # split the profile
        profile = [float(x) for x in row[profile_col].split(sep)]
        region_length = len(profile)
        
        if all(nt == 0 for nt in profile):
            # FIXME: we have to get all variables filled in here!

            # we know what the answer should be in such a case...
            ##slope = 0.0
            ##intercept = 0.0
            ##r_value = 0.0
            
            ##feat_list = [slope, intercept, r_value]
            
            # set the feature lists
            ##prof_mean_norm_feat_list = feat_list

            # set ratios to zero for cases with no coverage
            ratio_5p_spliced_in_all = 0.0
            ratio_3p_spliced_out_all = 0.0

            ratio_5p_spliced_in_border = 0.0
            ratio_3p_spliced_out_border = 0.0

            ratio_5p_unspliced = 0.0
            ratio_3p_unspliced = 0.0

            # set the normalized entropy (=entropy efficiency)
            entropy_efficiency = 0.0

        else:

            # normalize by the profile gene expression
            ##prof_expr_norm = [ p*1000/row[gene_expr_col] for p in profile ]

            # fit linear model to normalized by profile gene expression
            ##prof_expr_norm_feat_list = self.fit_linear_model_to_profile(prof_expr_norm)
        
            # normalize the profile by the profile mean
            ##prof_mean_norm = [ p*profile_len/sum(profile) for p in profile ]
            ##prof_mean_norm_feat_list = self.fit_linear_model_to_profile(prof_mean_norm)

            # calculte the normalized entropy (=entropy efficiency)
            prof_probability_norm = [ p/sum(profile) for p in profile ]
            entropy_efficiency = -sum([ p*np.log(p+np.finfo(float).eps) for p in prof_probability_norm ]) / np.log(region_length)

            # -----------------------------------------------------------------
            # Fit polynomial to the inverse CDF
#            inverse_cdf_profile_step_size = 0.05
#            inverse_cdf_quantiles = np.arange(0.0+inverse_cdf_profile_step_size,
#                                              1.0,inverse_cdf_profile_step_size)
#            inverse_cdf_profile_norm = \
#                self.inverse_cdf_profile(profile_probability_norm=prof_probability_norm,
#                                         profile_length=profile_len,
#                                         quantiles=inverse_cdf_quantiles)

            # fit a polynomial to the normalized inverse cdf
#            polynomial_fit_results = \
#                self.fit_polynomial_model_to_profile(x_axis=inverse_cdf_profile_norm,
#                                                     y_axis=inverse_cdf_quantiles,
#                                                     diagnostic_plots_dir_path=results_dir_path,
#                                                     region_id=row.name)

            # -----------------------------------------------------------------
            # FIXME: Get the quartiles at which we reach 
            inverse_cdf_quantiles = np.array([0.05, 0.95])
            # inverse_cdf_profile_step_size = 0.01
            # inverse_cdf_quantiles = np.arange(0.9+inverse_cdf_profile_step_size,
                                              # 1.0,inverse_cdf_profile_step_size)
            inverse_cdf_terminal_profile_norm = \
                self.inverse_cdf_profile(profile_probability_norm=prof_probability_norm,
                                         quantiles=inverse_cdf_quantiles)


            # calculate ratios
            mean_profile_5p = np.mean(profile[:min([10,len(profile)])])+np.finfo(float).eps
            mean_profile_3p = np.mean(profile[max([-10,-len(profile)]):])+np.finfo(float).eps
            
            ratio_5p_spliced_in_all = row['SpliceInAll']/mean_profile_5p
            ratio_3p_spliced_out_all = row['SpliceOutAll']/mean_profile_3p

            ratio_5p_spliced_in_border = row['SpliceInBorders']/mean_profile_5p
            ratio_3p_spliced_out_border = row['SpliceOutBorders']/mean_profile_3p

            ratio_5p_unspliced = row['Unspliced_5pSS']/mean_profile_5p
            ratio_3p_unspliced = row['Unspliced_3pSS']/mean_profile_3p

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # PROFILE-INDEPENDENT FEATURES
        # ---------------------------------------------------------------------
        # Create new feature that adds-up the splice-IN/OUT and crossing-IN/OUT
        # borders
        ReadsIN_borders = row['SpliceInBorders'] + row['Unspliced_5pSS']
        ReadsOUT_borders = row['SpliceOutBorders'] + row['Unspliced_3pSS']
        ReadsOUTvsIN_borders = (ReadsOUT_borders+np.finfo(float).eps)/(ReadsIN_borders+np.finfo(float).eps)
        # all
        ReadsIN_all = row['SpliceInAll'] + row['Unspliced_5pSS']
        ReadsOUT_all = row['SpliceOutAll'] + row['Unspliced_3pSS']
        ReadsOUTvsIN_all = (ReadsOUT_all+np.finfo(float).eps)/(ReadsIN_all+np.finfo(float).eps)

        # Calculate further IN-OUT ratios
        # FIXME: remove the following features
        SpliceOUTvsIN_all = (row['SpliceOutAll']+np.finfo(float).eps)/(row['SpliceInAll']+np.finfo(float).eps)
        SpliceOUTvsIN_borders = (row['SpliceOutBorders']+np.finfo(float).eps)/(row['SpliceInBorders']+np.finfo(float).eps)
        
        ##SpliceOUTbordersVSall = (row['SpliceOutBorders']+np.finfo(float).eps)/(row['SpliceOutAll']+np.finfo(float).eps)
        SpliceINbordersVSall = (row['SpliceInBorders']+np.finfo(float).eps)/(row['SpliceInAll']+np.finfo(float).eps)

        ## 
        GeneExpressionPerKBApproximated = row[gene_expr_col]

        ##
        RegionExpression = ((row['SpliceInAll'] + row['SpliceOutAll'] + row['Unspliced_5pSS'] + row['Unspliced_3pSS'] + row['UnsplicedExon'])/float(region_length))*1000

        ##
        RegionExpressionRatio = RegionExpression/GeneExpressionPerKBApproximated


        
        # Normalization by gene expression
        ##RegionExpressionPerKBInclBackground = np.mean(prof_expr_norm)
        
        # create the results
        # TODO: add maximum splice in per nt feature (should make it possible to 
        #       better distinguish background from terminal exons).
        
        # FIXME: delete the following documentation once we know we do not need this features anymore
        #          'profile_expr_norm_slope' : prof_expr_norm_feat_list[0],
        #          'profile_expr_norm_intercept' : prof_expr_norm_feat_list[1],
        #          'profile_expr_norm_r_value' : prof_expr_norm_feat_list[2], # Delete this feature!
        #          'profile_mean_norm_r_value' : prof_mean_norm_feat_list[2], # Delete this feature!

        # FIXME: 
        #          'profile_mean_norm_slope' : prof_mean_norm_feat_list[0],
        #          'profile_mean_norm_intercept' : prof_mean_norm_feat_list[1],

        

        results = {'SpliceInAll_vs_profile_ratio' : ratio_5p_spliced_in_all,
                   'SpliceOutAll_vs_profile_ratio' : ratio_3p_spliced_out_all,
                   'SpliceInBorders_vs_profile_ratio' : ratio_5p_spliced_in_border,
                   'SpliceOutBorders_vs_profile_ratio' : ratio_3p_spliced_out_border,
                   'CrossingInBorders_vs_profile_ratio' : ratio_5p_unspliced,
                   'CrossingOutBorders_vs_profile_ratio' : ratio_3p_unspliced,
                   'entropy_efficiency' : entropy_efficiency,
                   'region_length' : region_length,
                   'ReadsOUTvsIN_all' : ReadsOUTvsIN_all,
                   'SpliceINbordersVSall' : SpliceINbordersVSall,
                   'RegionExpressionRatio': RegionExpressionRatio}
                   ##'ReadsOUTvsIN_borders' : ReadsOUTvsIN_borders,
                   # 'SpliceOUTvsIN_all' : SpliceOUTvsIN_all,         # FIXME: delete later   
                   # 'SpliceOUTvsIN_borders' : SpliceOUTvsIN_borders, # FIXME: delete later
                   ##'SpliceOUTbordersVSall' : SpliceOUTbordersVSall, 
                   # 'ApproximatedGeneExpressionPerKB' : GeneExpressionPerKBApproximated,
                   # 'RegionExpression': RegionExpression,

                   #'RegionExpressionPerKBInclBackground' : RegionExpressionPerKBInclBackground}, # this feature represents the region length
                   #'PolynomialFitToCDFParameterA' : polynomial_fit_results[0],
                   #'PolynomialFitToCDFParameterB' : polynomial_fit_results[1],
                   #'PolynomialFitToCDFRsquared' : polynomial_fit_results[2]}

        for idx, quantile in enumerate(inverse_cdf_terminal_profile_norm):
            results[("absCDF_quant"+str(inverse_cdf_quantiles[idx]))] = \
                inverse_cdf_terminal_profile_norm[idx]

                    # We do not use them because we use the crossings to filter the positive training set 
                    # and in introns the crossing will not be zero in introns
                    # 'CrossingInBorders_vs_profile_ratio' : ratio_5p_unspliced, 
                    # 'CrossingOutBorders_vs_profile_ratio' : ratio_3p_unspliced,

        # Tell the module which features we have
        self.features = results.keys()
        
        # return the geatures
        return(pd.Series(results))


    def add_features_to_training_dataframes(self, 
                                            output_files_dir, 
                                            nr_to_subsample="all", 
                                            verbose=False):

        """
        Method that calculates features and adds it to the training dataframes.
        """

        # be sure we have converted it to a string
        nr_to_subsample = str(nr_to_subsample)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Check whether we have everything in place needed
        # ---------------------------------------------------------------------
        # Check whether we have the data needed for the next step
        if self.terminal_exon_training_data is None:
            sys.stderr.write("ERROR: no terminal exon training data ('terminal_exon_training_data') available!")
            sys.exit(-1)

        # Check whether we have the data needed for the next step
        if self.intermediate_exon_training_data is None:
            sys.stderr.write("ERROR: no intermediate exon training data ('intermediate_exon_training_data') available!")
            sys.exit(-1)

        # Check whether we have the data needed for the next step
        if self.background_training_data is None:
            sys.stderr.write("ERROR: no background training data ('background_training_data') available!")
            sys.exit(-1)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Get numbers needed for making decisions on the data set size
        # ---------------------------------------------------------------------
        # determine how many data we have available
        nr_terminal_exon_sets = len(self.terminal_exon_training_data)
        nr_intermediate_exon_sets = len(self.intermediate_exon_training_data)
        nr_background_region_sets = len(self.background_training_data)
        nr_available_data_sets = min(nr_terminal_exon_sets, 
                                     nr_intermediate_exon_sets,
                                     nr_background_region_sets)

        # check if we got a number of wished data sets
        if (nr_to_subsample.isdigit()):
            
            # convert to int 
            nr_wanted_data_sets = int(nr_to_subsample)

            # check if we have enough data available
            if nr_available_data_sets >= nr_wanted_data_sets:

                nr_data_sets_to_sample = nr_wanted_data_sets

                if verbose: 
                    sys.stdout.write("Sampling %s data sets from each training class...\n" \
                                    %(nr_data_sets_to_sample))

                # TE
                # random.seed()
                rows = \
                    random.sample(self.terminal_exon_training_data.index, 
                                  (nr_terminal_exon_sets - nr_data_sets_to_sample))
                self.terminal_exon_training_data.drop(rows, inplace = True)

                # IE
                # random.seed()
                rows = \
                    random.sample(self.intermediate_exon_training_data.index, 
                                  (nr_intermediate_exon_sets - nr_data_sets_to_sample))
                self.intermediate_exon_training_data.drop(rows, inplace = True)

                # BG
                # random.seed()
                rows = \
                    random.sample(self.background_training_data.index, 
                                  (nr_background_region_sets - nr_data_sets_to_sample))
                self.background_training_data.drop(rows, inplace = True)

            else:
                nr_data_sets_to_sample = "max_equal_size"

        # use all data
        if str(nr_to_subsample) == "all":
            if verbose:
                sys.stdout.write("Using all data sets from each training class...\n")
        
        # use the maximum possible when choosing equal sized training sets
        elif str(nr_to_subsample) == "max_equal_size":

            if verbose:
                sys.stdout.write("Using maximum possible number (n=%i) of " \
                                +"data sets from each training class so " \
                                +"that all of them have the same size...\n" \
                                %(nr_available_data_sets))

                # TE
                # random.seed()
                rows = \
                    random.sample(self.terminal_exon_training_data.index, 
                                  (nr_terminal_exon_sets - nr_available_data_sets))
                self.terminal_exon_training_data.drop(rows, inplace = True)

                # IE
                # random.seed()
                rows = \
                    random.sample(self.intermediate_exon_training_data.index, 
                                  (nr_intermediate_exon_sets - nr_available_data_sets))
                self.intermediate_exon_training_data.drop(rows, inplace = True)

                # BG
                # random.seed()
                rows = \
                    random.sample(self.background_training_data.index, 
                                  (nr_background_region_sets - nr_available_data_sets))
                self.background_training_data.drop(rows, inplace = True)

        
        # in all other cases we do not understand what the user wants
        elif not nr_to_subsample.isdigit():
            sys.stderr.write(("ERROR: invalid input for 'nr_to_subsample' " \
                             +"parameter in MachineLearningUnit.add_features_to_training_dataframes()."))
            sys.exit(-1)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Calculate features for TERMINAL EXONS
        # ---------------------------------------------------------------------
        # get the features
        if verbose: sys.stdout.write("Calculating features for terminal exon training data...\n")

        # calculate the features
        TE_feat = self.terminal_exon_training_data.merge(self.terminal_exon_training_data.apply(self.calculate_features, axis=1, 
                                                         results_dir_path=os.path.join(output_files_dir, "terminal_exon_training_data")), 
                                                         left_index=True, right_index=True)

        # drop Na values (might occure if it was not possible to calculate one or more features)
        nr_TE_datasets = TE_feat.shape[0]
        TE_feat.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
        
        if verbose: 
            data_fraction_with_features = float(TE_feat.shape[0]) / float(nr_TE_datasets)
            sys.stdout.write(" :: terminal exon training data set fraction for which features could be calculated: %.2f\n" % (data_fraction_with_features))

        # overwrite the old version that lacks features
        self.terminal_exon_training_data = TE_feat.copy()
        
        # clean up...
        del(TE_feat)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Calculate features for INTERMEDIATE EXONS
        # ---------------------------------------------------------------------
        # get the features
        if verbose: sys.stdout.write("Calculating features for intermediate exon training data...\n")

        # calculate the features
        IE_feat = self.intermediate_exon_training_data.merge(self.intermediate_exon_training_data.apply(self.calculate_features, axis=1,
                                                             results_dir_path=os.path.join(output_files_dir, "intermediate_exon_training_data")), 
                                                             left_index=True, right_index=True)
        
        # drop Na values (might occure if it was not possible to calculate one or more features)
        nr_IE_datasets = IE_feat.shape[0]
        IE_feat.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)

        if verbose: 
            data_fraction_with_features = float(IE_feat.shape[0]) / float(nr_IE_datasets)
            sys.stdout.write(" :: intermediate exon training data set fraction for which features could be calculated: %.2f\n" % (data_fraction_with_features))
        
        # overwrite the old version that lacks features
        self.intermediate_exon_training_data = IE_feat.copy()

        # clean up...
        del(IE_feat)

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # Calculate features for BACKGROUND REGIONS
        # ---------------------------------------------------------------------
        # get the features
        if verbose: sys.stdout.write("Calculating features for background regions...\n")

        # calculate the features
        BG_feat = self.background_training_data.merge(self.background_training_data.apply(self.calculate_features, axis=1, 
                                                      results_dir_path=os.path.join(output_files_dir, "background_training_data")), 
                                                      left_index=True, right_index=True)
        
        # drop Na values (might occure if it was not possible to calculate one or more features)
        nr_BG_datasets = BG_feat.shape[0]
        BG_feat.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
        
        if verbose: 
            data_fraction_with_features = float(BG_feat.shape[0]) / float(nr_BG_datasets)
            sys.stdout.write(" :: background region training data set fraction for which features could be calculated: %.2f\n" % (data_fraction_with_features))

        # overwrite the old version that lacks features
        self.background_training_data = BG_feat.copy()

        # clean up...
        del(BG_feat)


    def sample_validation_data_from_training_data(self, fraction, verbose = True):

        """
        Method that samples validation data from the training_df and stores it in validation_df.
        """

        if verbose: 
            sys.stdout.write("Sampling validation data (%.2f%% of training data)...\n" \
                             %(fraction*100))

        # randomly select the speficied fraction from the training data
        # random.seed()
        rows = random.sample(self.training_df.index, int( math.floor( fraction*len(self.training_df.index))))

        #split in training and testing datasets
        self.validation_df = self.training_df.ix[rows]
        self.training_df.drop(rows, inplace = True)
        
        if verbose: 
            sys.stdout.write(" > Final training data set size: %s data entries having %s features.\n"
                             %(str(self.training_df.shape[0]), str(self.training_df.shape[1])))
            sys.stdout.write(" > Final validation data set size: %s data entries having %s features.\n"
                             %(str(self.validation_df.shape[0]), str(self.validation_df.shape[1])))


    def load_training_data(self, 
                           training_data_set_size = "max_equal_size", 
                           validation_data_fraction = 0.2, 
                           output_files_dir = None,
                           run_number = 0,
                           verbose = True):

        """
        Method that samples a specific number of training data from each
        training data set class and then randomly selects the specified fraction
        of validation data.
        """

        training_data_set_size = str(training_data_set_size)

        # check if we have features available
        if (self.features is None) or (len(self.features) < 1):
            sys.stderr.write(("ERROR: no features have been calculated yet.\n"))
            sys.exit(-1)

        # determine how many data we have available
        nr_terminal_exon_sets = len(self.terminal_exon_training_data)
        nr_intermediate_exon_sets = len(self.intermediate_exon_training_data)
        nr_background_region_sets = len(self.background_training_data)
        nr_sets_available = min(nr_terminal_exon_sets, 
                                nr_intermediate_exon_sets,
                                nr_background_region_sets)

        # check if we got a number of wished data sets
        if (training_data_set_size.isdigit()): 
            
            # convert to int 
            nr_wanted_data = int(training_data_set_size)

            # check if we have enough data available
            if nr_sets_available >= nr_wanted_data:
                self.terminal_exons_features = self.terminal_exon_training_data[self.features].sample(n=nr_wanted_data).copy()
                self.intermediate_exons_features = self.intermediate_exon_training_data[self.features].sample(n=nr_wanted_data).copy()
                self.background_regions_features = self.background_training_data[self.features].sample(n=nr_wanted_data).copy()
            else:
                sys.stderr.write(("WARNING: there are not %s data sets " \
                                 +"available for each training class! Thus, " \
                                 +"the maximum possible number (%s) " \
                                 +"will be used.\n") % (str(nr_wanted_data), str(nr_sets_available)))
                training_data_set_size = "max_equal_size"

        # use all data
        if str(training_data_set_size) == "all":
            self.terminal_exons_features = self.terminal_exon_training_data[self.features].copy()
            self.intermediate_exons_features = self.intermediate_exon_training_data[self.features].copy()
            self.background_regions_features = self.background_training_data[self.features].copy()

        # use the maximum possible when choosing equal sized training sets
        elif str(training_data_set_size) == "max_equal_size":
            self.terminal_exons_features = self.terminal_exon_training_data[self.features].sample(n=nr_sets_available).copy()
            self.intermediate_exons_features = self.intermediate_exon_training_data[self.features].sample(n=nr_sets_available).copy()
            self.background_regions_features = self.background_training_data[self.features].sample(n=nr_sets_available).copy()
        
        # in all other cases we do not understand what the user wants
        elif not training_data_set_size.isdigit():
            sys.stderr.write(("ERROR: invalid input for 'training_data_set_size' " \
                             +"parameter in MachineLearningUnit.load_training_data()."))
            sys.exit(-1)

        # add the classes
        self.terminal_exons_features[self.class_col] = self.terminal_exon_class
        self.intermediate_exons_features[self.class_col] = self.intermediate_exon_class
        self.background_regions_features[self.class_col] = self.background_region_class

        # concatenate the training data
        if verbose: sys.stdout.write("Concatenating training data...\n")
        self.training_df = pd.concat([self.terminal_exons_features, 
                                      self.intermediate_exons_features, 
                                      self.background_regions_features])

        # write the training data to a file (in case a output dir was specified).
        if output_files_dir is not None:
            self.write_training_df_to_file(training_df_file_path=os.path.join(output_files_dir, "run_number_" + str(run_number) + "_training_data.tsv"),
                                           verbose=True)

        # select validation data
        if validation_data_fraction > 0.0:
            self.sample_validation_data_from_training_data(fraction = validation_data_fraction, verbose = verbose)

            # write the validation data to a file (in case a output dir was specified).
            if output_files_dir is not None:
                self.write_validation_df_to_file(validation_df_file_path=os.path.join(output_files_dir, "run_number_" + str(run_number) + "_validation_data.tsv"),
                                                 verbose=True)


    def plot_confusion_matrix(self,
                              cm, 
                              file_path,
                              normalize=False,
                              title='Confusion matrix',
                              cmap=plt.cm.Blues):
        """
        Method that plots a confusion matrix.
        Normalization can be applied by setting "normalize=True".
        """

        plt.figure()
        plt.imshow(cm, interpolation='nearest', cmap=cmap)
        plt.title(title)
        plt.colorbar()
        tick_marks = np.arange(len(self.region_classes))
        plt.xticks(tick_marks, self.region_classes, rotation=45)
        plt.yticks(tick_marks, self.region_classes)

        if normalize:
            cm = cm.astype('float') / cm.sum(axis=1)[:,np.newaxis]

        thresh = cm.max() / 2.0
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            plt.text(j, i, str("%0.2f" % cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")

        plt.tight_layout()
        plt.ylabel('True label')
        plt.xlabel('Predicted label')
        plt.savefig(file_path)
        plt.close('all')


    def train_classifier(self, results_dir_path, nr_of_train_vs_test_runs = 25, verbose = False):

        """
        Method to train the classifier. 'nr_of_train_vs_test_runs' runs will be done 
        and the results will be reported in the results_dir_path directory. However, only the
        last classifier will be stored in the MachineLearningUnit object for subsequent use.
        """

        if verbose: sys.stdout.write("Training classifier...\n")

        # Check whether we have data to train the classifier on
        if self.training_df is None:
            sys.stderr.write("ERROR: no training data ('training_df') available!")
            sys.exit(-1)

        # _____________________________________________________________________________
        # -----------------------------------------------------------------------------
        # Training
        # -----------------------------------------------------------------------------
        n_neighbors = 3
        weights = 'uniform'
        ##weights = 'distance'

        # create results directory name
        results_dir = os.path.join(results_dir_path, ('KNeighborsClassifier_%s_%sNodes' % (weights, str(n_neighbors))))
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)

        # create lists that we can use for printing multiple results
        accuracy_scores_list = list()
        precision_scores_list = list()
        recall_scores_list = list()
        f1score_scores_list = list()

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # create multiple runs so that we see how stable our results are
        for i in range(nr_of_train_vs_test_runs):
            
            # _________________________________________________________________
            # -----------------------------------------------------------------
            # Split data
            # -----------------------------------------------------------------
            X_train, X_test, y_train, y_test = \
                train_test_split(self.training_df[self.features], 
                                 self.training_df[self.class_col], 
                                 test_size = 0.2, 
                                 random_state = random.randint(0,1000))

            # _________________________________________________________________
            # -----------------------------------------------------------------
            # Model training
            # -----------------------------------------------------------------
            # get the classifier
            self.classifier = neighbors.KNeighborsClassifier(n_neighbors, 
                                                             weights=weights)
                
            # fit the classifier
            self.classifier.fit(X_train, y_train)

            # _________________________________________________________________
            # -----------------------------------------------------------------
            # Model validation
            # -----------------------------------------------------------------
            # TODO: 
            # Suggestion from Andrea
            # Use predict_proba() -> returns a list of the propabilities for each class
            # -----------------------------------------------------------------
            # perform predictions on the test set
            y_pred = self.classifier.predict(X_test)
            y_true = y_test.tolist()
            
            # -----------------------------------------------------------------
            # calculate the accuracy
            # -----------------------------------------------------------------
            accuracy = accuracy_score(y_true = y_true, y_pred = y_pred, 
                                      normalize = True)
            accuracy_scores_list.append(accuracy)

            # -----------------------------------------------------------------
            # create confusion matrixes
            # -----------------------------------------------------------------
            cnf_matrix = confusion_matrix(y_true = y_true, y_pred = y_pred, 
                                          labels = self.region_classes)

            # Plot non-normalized confusion matrix
            # -----------------------------------------------------------------
            # TODO:
            # Suggestion from Andrea
            # use .roc_score
            # use .roc_auc_score
            # -----------------------------------------------------------------
            cm_file_name = ("normalized_confusion_matrix_RUN_%s.png" % str(i))
            self.plot_confusion_matrix(cnf_matrix, 
                                       file_path=os.path.join(results_dir, cm_file_name),
                                       normalize=True,
                                       title='Confusion matrix')
            
            # -----------------------------------------------------------------
            # create precission, recall and F1-scores
            # -----------------------------------------------------------------
            precision_scores_list.append(metrics.precision_score(y_true = y_true, y_pred = y_pred, average='macro'))
            recall_scores_list.append(metrics.recall_score(y_true = y_true, y_pred = y_pred, average='micro'))
            f1score_scores_list.append(metrics.f1_score(y_true = y_true, y_pred = y_pred, average='weighted'))

        # _____________________________________________________________________
        # ---------------------------------------------------------------------
        # print accuracy
        plt.hist(accuracy_scores_list, bins = np.arange(0.0, 1.005, 0.005))
        plt.savefig(os.path.join(results_dir, "accuracy.png"))
        plt.close('all')

        # print precission
        plt.hist(precision_scores_list, bins = np.arange(0.0, 1.005, 0.005))
        plt.savefig(os.path.join(results_dir, "precission.png"))
        plt.close('all')

        # print recall
        plt.hist(recall_scores_list, bins = np.arange(0.0, 1.005, 0.005))
        plt.savefig(os.path.join(results_dir, "recall.png"))
        plt.close('all')

        # print F1-score
        plt.hist(f1score_scores_list, bins = np.arange(0.0, 1.005, 0.005))
        plt.savefig(os.path.join(results_dir, "f1.png"))
        plt.close('all')


    # perfomance function on specific subset of features, number_of_randomization times
    def training(self, 
                 classifier,
                 features_for_training, 
                 number_of_randomization = 10):

        # create lists that we can use for printing multiple results
        training_scores_list = list()
        
        # create multiple runs so that we see how stable our results are
        for i in range(number_of_randomization):

            X_train, X_test, y_train, y_test = \
                train_test_split(self.training_df[features_for_training], 
                                 self.training_df[self.class_col], 
                                 test_size = 0.2, 
                                 random_state = random.randint(0,65000))
            
            # use the classifier with the best performing parameters
            # perform predictions on the test set
            y_pred = \
                self.classifier_funcdict[classifier](X_train = X_train,
                                                     y_train = y_train).predict(X_test)
            y_true = y_test.tolist()

            # training score from my_score_func
            #training_scores_list.append(my_score_func(y_true, y_pred, exponent=exponent))
            training_scores_list.append(metrics.f1_score(y_true, y_pred, average="macro"))
        
        return([ np.mean(training_scores_list), np.std(training_scores_list, ddof=1) ])


    def greedy_feature_selection(self, classifier, manually_selected_features = list(), number_of_randomization = 10, verbose = True):

        """
        Method that performs Greedy Feature Selection.
        """

        if verbose: 
            sys.stdout.write(("Performing Greedy Feature Selection " \
                             +"(using %i independent runs per " \
                             +"feature)...\n") \
                             % (number_of_randomization))

        # Calculate the t-value that we consider to be significant
        significant_t_value = stats.t.ppf(q = 0.9, df = number_of_randomization, loc = 0, scale = 1)

        # start from the full list of features and select always 
        # the one that performs best
        features_to_test = copy.copy(self.features)

        # initialize the list for selected features
        self.selected_features = copy.copy(manually_selected_features)

        # remove features that are already selected
        for feature in self.selected_features:
            features_to_test.remove(feature)


        # create lists for the feature scores and the stdevs
        feature_scores = list()
        feature_stdevs = list()

        # do the greedy
        while len(features_to_test) > 0:
            
            # initialize variables for investigating features
            max_score = 0.0
            max_stdev = 0.0
            best_feature = ""

            # iterate over all the features and find out which one
            # performs best
            for feat in features_to_test:
                
                # add the current feature to the already selected features
                features = self.selected_features + [feat]
                
                # train
                score, score_stdev = \
                  self.training(classifier = classifier,
                                features_for_training = features, 
                                number_of_randomization = number_of_randomization)
                
                if verbose:
                    sys.stdout.write((" :: :: Score: %.2f%% by using features: %s\n") \
                                     %(score, features))

                if max_score < score:
                    max_score = score
                    max_stdev = score_stdev
                    best_feature = feat
            
            # in case we found already the best performing feature (of all)
            # we want to know now, whether the additional features contribute
            # significantly to the performance.
            if len(feature_scores) > 0:
                                
                # Add features only in case they contribute sugnificantly to the predictions.
                t_value = ( max_score - feature_scores[-1] ) / math.sqrt( ( max_stdev**2 + feature_stdevs[-1]**2 ) / number_of_randomization )

                # Print the t-value of the feature
                if verbose:
                    sys.stdout.write((" :: :: t-value: %.2f by using features: %s\n") \
                                     %(t_value, best_feature))
                
                # in case the feature does not contribute significantly to 
                # the performance we stop searching for additional features.
                if significant_t_value > t_value:
                    break

            # store the max_score and the max_stdev so that we can later on find out
            # whether a novel feature significantly improves the perormance
            feature_scores.append(max_score)
            feature_stdevs.append(max_stdev)
            
            # add the best feature to the selected features
            self.selected_features.append(best_feature)
            
            # remove the best feature from the features we still want to test
            features_to_test.remove(best_feature)

            if verbose:
                sys.stdout.write((" :: Selected Features Score: %.2f%% by using features: %s\n") \
                                 %(max_score, self.selected_features))
                sys.stdout.flush()

        # give some final feedback
        if verbose:
            sys.stdout.write((" :: Finally Selected Features Score: %.2f%% by using features: %s\n") \
                             %(feature_scores[-1], self.selected_features))
            sys.stdout.flush()


    def create_kneighbors_classifier(self, X_train, y_train):


        """ 
        Method that finds the optimal parameters for a KNeighbors Classifier
        making use of the selected features and returns the best estimator.
        """

        #scorer = make_scorer(my_score_func, exponent = 2.0)
        scorer = make_scorer(metrics.f1_score, average = "macro")

        parameters = {'n_neighbors': range(2,15),
                      'algorithm': ['auto'], #['ball_tree', 'kd_tree', 'brute'],
                      'weights': ['uniform', 'distance'],
                      'p': range(1,3) }

        # create the KNeighbors Classifier
        kn = neighbors.KNeighborsClassifier()

        # get classifier for all parameter combinations
        clf = GridSearchCV(kn, parameters, scoring = scorer)
        clf.fit(X_train, y_train)

        # return the best performing estimator
        return(clf.best_estimator_)


    def create_multiclass_SVC_classifier(self, X_train, y_train):


        """ 
        Method that finds the optimal parameters for a create_multi-class SVC
        Classifier making use of the selected features and returns the best estimator.
        """


        # return the classifier
        return("Implement me!")


    def create_Bayes_classifier(self, X_train, y_train):

        """ 
        Method that finds the optimal parameters for a Bayes
        Classifier making use of the selected features and returns the best estimator.
        """

        clf = BayesClassifier()
        clf.fit(X_train, y_train)

        # return the classifier
        return(clf)


    def train_classifier_on_selected_features(self, classifier, results_dir_path, nr_of_train_vs_test_runs = 25, verbose = False):

        """
        Method to trains the classifier on greedy selected features. 'nr_of_train_vs_test_runs' runs will be done 
        and the results will be reported in the results_dir_path directory. However, only the
        last classifier will be stored in the MachineLearningUnit object for subsequent use.
        """

        if verbose: sys.stdout.write("Training classifier on greedy selected features: %s\n" % (self.selected_features))

        # -----------------------------------------------------------------
        # Check that everything is in place in order to get started
        # -----------------------------------------------------------------
        # Check whether we have data to train the classifier on
        if self.selected_features is None:
            sys.stderr.write("ERROR: no greedy selected features ('greedy_selected_features') available!")
            sys.exit(-1)

        # Check whether we have data to train the classifier on
        if self.training_df is None:
            sys.stderr.write("ERROR: no training data ('training_df') available!")
            sys.exit(-1)

        # Check whether we have data to validate the classifier on
        if self.validation_df is None:
            sys.stderr.write("ERROR: no validation data ('validation_df') available!")
            sys.exit(-1)

        # -----------------------------------------------------------------
        # Get the best performing classifier.
        # -----------------------------------------------------------------        
        self.classifier_dict[classifier] = \
            self.classifier_funcdict[classifier](X_train = self.training_df[self.selected_features],
                                                   y_train = self.training_df[self.class_col])

        # -----------------------------------------------------------------
        # Create predictions for the validation data
        # -----------------------------------------------------------------        
        y_pred = self.classifier_dict[classifier].predict(self.validation_df[self.selected_features])
        
        # -----------------------------------------------------------------
        # Get the actual classes
        # -----------------------------------------------------------------        
        y_true = self.validation_df[self.class_col].tolist()
        
        # -----------------------------------------------------------------
        # create and plot confusion matrixes
        # -----------------------------------------------------------------
        cnf_matrix = confusion_matrix(y_true = y_true, y_pred = y_pred, 
                                      labels = self.region_classes)

        cm_file_name = "normalized_confusion_matrix.png"
        cm_file_path = os.path.join(results_dir_path, cm_file_name)
        self.plot_confusion_matrix(cnf_matrix, 
                                   file_path = cm_file_path,
                                   normalize = True,
                                   title = 'Confusion matrix')
        
        if verbose: 
            sys.stdout.write("Writing confusion matrix: %s\n" % (cm_file_path))

        # -----------------------------------------------------------------
        # Calculate validation scores
        # -----------------------------------------------------------------        
        accuracy = accuracy_score(y_true = y_true, y_pred = y_pred, normalize = True)
        f1_result = metrics.f1_score(y_true = y_true, y_pred = y_pred, average="macro")

        if verbose: 
            sys.stdout.write(" :: Accuracy: %.2f\n" % (accuracy))
            sys.stdout.write(" :: F1 score: %.2f\n" % (f1_result))


    def create_terminal_exon_candidates_dataframe(self, 
                                                  annotation, 
                                                  novel_terminal_output_file,
                                                  verbose=False):
        """
        Create terminal exon candidates dataframe
        """

        novel_terminal_exon_candidates_data_list = []

        for gene_id in annotation.genes:

            # if (annotation.genes[gene_id].intron_length > 0) and \
            #    (annotation.genes[gene_id].union_exon_length > 0) and \
            #    (annotation.genes[gene_id].GeneExpressionBackgroundPerKB >= 0) and \
            #    (annotation.genes[gene_id].GeneExpressionPerKBInclBackground > 0) and \
            #    (annotation.genes[gene_id].GeneExpressionPerKBwithoutBackground >0) and \
            #    (annotation.genes[gene_id].has_potential_novel_terminal_exon()):

            if (annotation.genes[gene_id].union_exon_length > 0) and \
               (annotation.genes[gene_id].total_reads > 0) and \
               (annotation.genes[gene_id].has_potential_novel_terminal_exon()):

                # terminal exons
                tmp_novel_terminal_exon_candidates_data_list = annotation.genes[gene_id].get_potential_novel_exons()

                for novel_terminal_exon_stats in tmp_novel_terminal_exon_candidates_data_list:

                    novel_terminal_exon_candidates_data_list.append(novel_terminal_exon_stats)

        labels = ["Region", "Annotation","GeneId","SpliceInAll", "SpliceInBorders", "SpliceOutAll", "SpliceOutBorders", "UnsplicedExon", "Unspliced_5pSS", "Unspliced_3pSS", "profile", "TotalGeneReads", "UnionExonLength", "GeneExpressionPerKBApproximated"] #, "IntronicReads", "IntronLength", "BackgroundPerKB", "GeneExpressionPerKBInclBackground", "GeneExpressionBackgroundPerKB", "GeneExpressionPerKBwithoutBackground", "BackgroundFraction"]

        self.novel_terminal_exon_candidates_data = pd.DataFrame.from_records(novel_terminal_exon_candidates_data_list, columns=labels)
        self.novel_terminal_exon_candidates_data.set_index(["Region","GeneId"], inplace=True)
        self.novel_terminal_exon_candidates_data.to_csv(novel_terminal_output_file, sep='\t', index=True)


    def filter_terminal_exon_candidates_that_overlap_with_annotated_exons(self,
                                                                          annotation,
                                                                          novel_terminal_output_file,
                                                                          sequencing_direction,
                                                                          exons_per_gene_bed_file,
                                                                          verbose):

        """ 
        Filter potential novel terminal exons that overlap with annotated exons
        """

        w = open(novel_terminal_output_file+"_tmp.bed", 'w')

        # reset index for novel terminal exons
        self.novel_terminal_exon_candidates_data.reset_index(inplace=True)

        for novel_exon in self.novel_terminal_exon_candidates_data[["Region","GeneId"]].iterrows():
            
            # get region
            region=novel_exon[1]["Region"]
            # split
            chrom,start,end,strand=region.split(":")
            # get gene id
            gene=novel_exon[1]["GeneId"]

            # write the bed file
            w.write("\t".join([chrom, str(int(start)-1), end, gene, ".", strand+"\n"]))

        w.close()

        # Read the novel terminal exons bed file as a bedtools object
        novel_terminal_exon_candidates_bed = pybedtools.BedTool(novel_terminal_output_file+"_tmp.bed")

        # Read the exons bed file as a bedtools object
        exons_per_gene_bed = pybedtools.BedTool(exons_per_gene_bed_file)

        strand_option=True
        if sequencing_direction == "unstranded":
            strand_option=False

        # Intetsect (bedtools intersect -wa -v -s -a <> -b <>)
        novel_terminal_exon_candidates = novel_terminal_exon_candidates_bed.intersect(exons_per_gene_bed, s=strand_option, v=True)

        # Find accepted regions
        accepted_novel_regions = []
        for novel_exon in novel_terminal_exon_candidates:
            chrom = novel_exon[0]
            start = str(int(novel_exon[1])+1)
            end = novel_exon[2]
            strand = novel_exon[5]

            region = ":".join([chrom, start, end, strand])
            accepted_novel_regions.append(region)

        accepted_novel_regions = list(set(accepted_novel_regions))

        # novel terminal exons 
        self.novel_terminal_exon_candidates_data = self.novel_terminal_exon_candidates_data[self.novel_terminal_exon_candidates_data['Region'].isin(accepted_novel_regions)]
        # set index
        self.novel_terminal_exon_candidates_data.set_index(["Region","GeneId"], inplace=True)
        # write out files
        self.novel_terminal_exon_candidates_data.to_csv(novel_terminal_output_file, sep='\t', index=True)


    def create_terminal_exon_readthough_candidates_dataframe(self, 
                                                             annotation, 
                                                             novel_readthough_output_file,
                                                             verbose = False):

        """
        Create terminal exon readthrough candidates dataframes
        """
        novel_terminal_exon_readthrough_candidates_data_list = []

        for gene_id in annotation.genes:

            if (annotation.genes[gene_id].union_exon_length > 0) and \
               (annotation.genes[gene_id].total_reads > 0) and \
               (annotation.genes[gene_id].has_potential_novel_readthrough_exon()):

                # readthrough exons extended
                tmp_novel_terminal_exon_readthrough_candidates_data_list = annotation.genes[gene_id].get_potential_novel_readthrough_exons()

                for novel_terminal_exon_stats in tmp_novel_terminal_exon_readthrough_candidates_data_list:

                    novel_terminal_exon_readthrough_candidates_data_list.append(novel_terminal_exon_stats)

        labels = ["Region", "Annotation","GeneId","SpliceInAll", "SpliceInBorders", "SpliceOutAll", "SpliceOutBorders", "UnsplicedExon", "Unspliced_5pSS", "Unspliced_3pSS", "profile", "TotalGeneReads", "UnionExonLength", "GeneExpressionPerKBApproximated"] #, "IntronicReads", "IntronLength", "BackgroundPerKB", "GeneExpressionPerKBInclBackground", "GeneExpressionBackgroundPerKB", "GeneExpressionPerKBwithoutBackground", "BackgroundFraction"]

        self.novel_terminal_readthrough_exon_candidates_data = pd.DataFrame.from_records(novel_terminal_exon_readthrough_candidates_data_list, columns=labels)
        self.novel_terminal_readthrough_exon_candidates_data.set_index("Region", inplace=True)
        self.novel_terminal_readthrough_exon_candidates_data.to_csv(novel_readthough_output_file, sep='\t', index=True)


    def min_profile_coverage_fraction_reached(self, 
                                              row, 
                                              profile_col = 'profile', 
                                              min_profile_coverage_fraction = 0.0, 
                                              sep = ','):

        """
        Method that checks whether profiles have low coverage.
        """

        # split the profile
        profile = [float(x) for x in row[profile_col].split(sep)]
        region_length = len(profile)

        # determine the fraction of covered nucleotides
        nr_zero_bases = profile.count(0)
        profile_coverage_fraction = 1 - (float(nr_zero_bases) / float(region_length))
        
        # check first if we want 
        if (profile_coverage_fraction < min_profile_coverage_fraction):
            return(False)
        else:
            return(True)


    def add_features_to_terminal_exon_candidates_dataframe(self, output_files_dir, verbose=True):

        """
        Method that calculates features and adds it to the novel terminal exon dataframe.
        """

        if verbose: sys.stdout.write("Calculating features for the novel terminal exon candidate regions...\n")

        # check whether we have data to train the classifier on
        if self.novel_terminal_exon_candidates_data is None:
            sys.stderr.write("ERROR: no novel terminal exon candidates data ('novel_terminal_exon_candidates_data') available!")
            sys.exit(-1)

        # filter out candidates with low profile coverage
        min_profile_coverage_fraction = 0.80
        sufficiently_covered_candidates_idx = \
          self.novel_terminal_exon_candidates_data.apply(self.min_profile_coverage_fraction_reached, axis=1, 
                                                         min_profile_coverage_fraction = min_profile_coverage_fraction)

        # only get the candidates that are sufficiently covered
        NTE_sufficiently_covered = \
          self.novel_terminal_exon_candidates_data.loc[sufficiently_covered_candidates_idx,:].copy()

        # give some feedback to the user
        if verbose: 
            
            data_fraction_with_coverage = \
              float(NTE_sufficiently_covered.shape[0]) / float(self.novel_terminal_exon_candidates_data.shape[0])
            
            sys.stdout.write((" :: terminal exon candidate data set " \
                             +"fraction for which sufficient coverage " \
                             +"(>=%.2f) is available and therefore will " \
                             +"be considered: %.2f (=%i candidates)\n") \
                             % (min_profile_coverage_fraction*100,
                                data_fraction_with_coverage, 
                                NTE_sufficiently_covered.shape[0]))

        # calculate the features
        NTE_feat = NTE_sufficiently_covered.merge(NTE_sufficiently_covered.apply(self.calculate_features, axis=1, 
                                                  results_dir_path=os.path.join(output_files_dir, "terminal_exon_candidates_data")), 
                                                  left_index=True, right_index=True)
        
        # drop Na values (might occure if it was not possible to calculate one or more features)
        nr_NTE_datasets = NTE_feat.shape[0]
        NTE_feat.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
        
        if verbose: 
            data_fraction_with_features = float(NTE_feat.shape[0]) / float(nr_NTE_datasets)
            sys.stdout.write(" :: fraction of considered terminal exon candidates for which features could be calculated: %.2f\n" % (data_fraction_with_features))

        # overwrite the old version that lacks features
        self.novel_terminal_exon_candidates_data = NTE_feat.copy()
        
        # clean up...
        del(NTE_feat)


    def load_terminal_exon_candidates(self, verbose=True):

        """
        Creates 'self.novel_terminal_exon_features', which contains features and a class column
        for novel terminal exons.
        """

        # get only the features
        self.novel_terminal_exon_features = self.novel_terminal_exon_candidates_data[self.selected_features].copy()

        # add the class
        self.novel_terminal_exon_features[self.class_col] = self.terminal_exon_class


    def fill_probabilities_dictionary(self, classifier, results_dir, verbose=True):

        """
        Use the classifier to classify given regions based on their features. 
        """

        if verbose: sys.stdout.write("Calculating probabilities for terminal exon candidates...\n")

        # Check whether we have data to train the classifier on
        if classifier not in self.classifier_dict:
            sys.stderr.write("ERROR: no classifier available!")
            sys.exit(-1)

        # read in the file with the novel terminal exon candidates
        X = self.novel_terminal_exon_features[self.selected_features]

        for index, row in X.iterrows():

            tmp_probabilities = self.classifier_dict[classifier].predict_proba(row)

            sys.stderr.write("row\n")
            sys.stderr.write(str(row))
            sys.stderr.write("\n")

            if index in self.terminal_probabilities_dict:
                self.terminal_probabilities_dict[index].append(tmp_probabilities['terminal_probability'])
            else:
                self.terminal_probabilities_dict[index] = [tmp_probabilities['terminal_probability']]

            if index in self.intermediate_probabilities_dict:
                self.intermediate_probabilities_dict[index].append(tmp_probabilities['intermediate_probability'])
            else:
                self.intermediate_probabilities_dict[index] = [tmp_probabilities['intermediate_probability']]
            
            if index in self.background_probabilities_dict:
                self.background_probabilities_dict[index].append(tmp_probabilities['background_probability'])
            else:
                self.background_probabilities_dict[index] = [tmp_probabilities['background_probability']]

        pickle.dump( self.terminal_probabilities_dict, open( os.path.join(results_dir, 'DEBUG_terminal_probabilities_dict.p'), "wb" ) )

        



    def classify_terminal_exon_candidates(self, classifier, results_dir, verbose=True):

        mean_terminal_probabilities_dict = dict()
        mean_intermediate_probabilities_dict = dict()
        mean_background_probabilities_dict = dict()
        classification_dict = dict()

        for index in self.terminal_probabilities_dict:

            terminal_probability = np.mean(self.terminal_probabilities_dict[index])
            mean_terminal_probabilities_dict[index] = terminal_probability

            intermediate_probability = np.mean(self.intermediate_probabilities_dict[index])
            mean_intermediate_probabilities_dict[index] = intermediate_probability

            background_probability = np.mean(self.background_probabilities_dict[index])
            mean_background_probabilities_dict[index] = background_probability

            if ((terminal_probability >= intermediate_probability) and \
                (terminal_probability >= background_probability)):
                classification_dict[index] = 'terminal'
            elif ((intermediate_probability >= terminal_probability) and \
                (intermediate_probability >= background_probability)):
                classification_dict[index] = 'intermediate'
            else:
                classification_dict[index] = 'background'

        terminal_probability_df = pd.DataFrame.from_dict(mean_terminal_probabilities_dict, orient='index')
        terminal_probability_df.reset_index(inplace=True)
        terminal_probability_df.columns = ["ind","terminal_probability"]
        # terminal_probability_df["actual_indexes"] = [ literal_eval(a) for a in terminal_probability_df["ind"].tolist() ]
        terminal_probability_df[['Region', 'GeneId']] = terminal_probability_df["ind"].apply(pd.Series)
        terminal_probability_df = terminal_probability_df[["Region", "GeneId", "terminal_probability"]]
        terminal_probability_df.set_index(["Region","GeneId"], inplace=True)

        intermediate_probability_df = pd.DataFrame.from_dict(mean_intermediate_probabilities_dict, orient='index')
        intermediate_probability_df.reset_index(inplace=True)
        intermediate_probability_df.columns=["ind","intermediate_probability"]
        # intermediate_probability_df["actual_indexes"] = [ literal_eval(a) for a in intermediate_probability_df["ind"].tolist() ]
        intermediate_probability_df[['Region', 'GeneId']] = intermediate_probability_df["ind"].apply(pd.Series)
        intermediate_probability_df = intermediate_probability_df[["Region", "GeneId", "intermediate_probability"]]
        intermediate_probability_df.set_index(["Region","GeneId"], inplace=True)

        background_probability_df = pd.DataFrame.from_dict(mean_background_probabilities_dict, orient='index')
        background_probability_df.reset_index(inplace=True)
        background_probability_df.columns=["ind","background_probability"]
        # background_probability_df["actual_indexes"] = [ literal_eval(a) for a in background_probability_df["ind"].tolist() ]
        background_probability_df[['Region', 'GeneId']] = background_probability_df["ind"].apply(pd.Series)
        background_probability_df = background_probability_df[["Region", "GeneId", "background_probability"]]
        background_probability_df.set_index(["Region","GeneId"], inplace=True)

        classification_df = pd.DataFrame.from_dict(classification_dict, orient='index')
        classification_df.reset_index(inplace=True)
        classification_df.columns=["ind","classification"]
        # classification_df["actual_indexes"] = [ literal_eval(a) for a in classification_df["ind"].tolist() ]
        classification_df[['Region', 'GeneId']] = classification_df["ind"].apply(pd.Series)
        classification_df = classification_df[["Region", "GeneId", "classification"]]
        classification_df.set_index(["Region","GeneId"], inplace=True)

        # read in the file with the novel terminal exon candidates
        X = self.novel_terminal_exon_features[self.selected_features]

        X_with_probabilities = reduce(lambda left,right: pd.merge(left,right, left_index=True, right_index=True), [X, terminal_probability_df, intermediate_probability_df, background_probability_df, classification_df])

        X_tmp = X_with_probabilities.copy()
        X_tmp.reset_index(inplace=True)
        X_tmp["chromosome"], X_tmp["start"], X_tmp["end"], X_tmp["strand"] = X_tmp["Region"].str.split(':',3).str

        self.selected_novel_terminal_exons = pd.DataFrame()

        # group by strand
        for strand_group in X_tmp[X_tmp["classification"]=="terminal"].groupby(["strand"]):

            if strand_group[0] == "+":
                # group by chromosome, start and gene id and concatenate to the final dataframe
                for geneid_chromosome_start in strand_group[1].groupby(["chromosome", "start", "GeneId"]):
                    current_terminal = geneid_chromosome_start[1].loc[[geneid_chromosome_start[1]["terminal_probability"].argmax()]]
                    self.selected_novel_terminal_exons = pd.concat([self.selected_novel_terminal_exons, current_terminal])

            if strand_group[0] == "-":
                # group by chromosome, end and gene id and concatenate to the final dataframe
                for geneid_chromosome_end in strand_group[1].groupby(["chromosome", "end", "GeneId"]):
                    current_terminal = geneid_chromosome_end[1].loc[[geneid_chromosome_end[1]["terminal_probability"].argmax()]]
                    self.selected_novel_terminal_exons = pd.concat([self.selected_novel_terminal_exons, current_terminal])

        # Write out the final terminal exons that we will use
        if not self.selected_novel_terminal_exons.empty:
            self.selected_novel_terminal_exons.set_index(["Region", "GeneId"], inplace=True)
            self.selected_novel_terminal_exons.to_csv(os.path.join(results_dir, 'classified_as_terminal_with_probabilities.tsv'), sep='\t', index=True)
        else:
            sys.stderr.write("[WARNING] No novel terminal exons were detected\n")
            self.selected_novel_terminal_exons = pd.DataFrame(columns=["Region","GeneId"])
            self.selected_novel_terminal_exons.to_csv(os.path.join(results_dir, 'classified_as_terminal_with_probabilities.tsv'), sep='\t', index=False)


    def classify_terminal_exon_candidates_original(self, classifier, results_dir, verbose=True):

        """
        Use the classifier to classify given regions based on their features. 
        """

        if verbose: sys.stdout.write("Classifying terminal exon candidates...\n")

        # Check whether we have data to train the classifier on
        if classifier not in self.classifier_dict:
            sys.stderr.write("ERROR: no classifier available!")
            sys.exit(-1)

        # read in the file with the novel terminal exon candidates
        X = self.novel_terminal_exon_features[self.selected_features]

        # classify the candidates
        # FIXME: here we have to use the predict_proba method and then
        #        filter out the class with the highest probability
        y_pred = self.classifier_dict[classifier].predict(X)

        # # determine the accuracy (even though this does not make too much sense here)
        # y_true = self.novel_terminal_exon_features[self.class_col]
        # accuracy = accuracy_score(y_true = y_true, y_pred = y_pred, normalize = True)
        
        # # create a confusion matrix
        # cnf_matrix = confusion_matrix(y_true = y_true, y_pred = y_pred, labels = self.region_classes)

        # # Plot the normalized confusion matrix
        # cm_file_name = "normalized_confusion_matrix_TERMINAL_EXON_CANDIDATES.png"
        # self.plot_confusion_matrix(cnf_matrix, 
        #                            file_path=os.path.join(results_dir, cm_file_name),
        #                            normalize=True,
        #                            title='Confusion matrix')

        # # write out files with the candidates that have been classified as novel terminal exons
        # self.novel_terminal_exons = X.loc[(y_pred == self.terminal_exon_class),:]
        # # FIXME: remove writing to a file
        # self.novel_terminal_exons.to_csv(os.path.join(results_dir,'classified_as_novel_terminal_exon.tsv'), sep='\t', index=True)

        # # write out files with the candidates that have been classified as intermediate exons
        # self.novel_intermediate_exons = X.loc[(y_pred == self.intermediate_exon_class),:]
        # # FIXME: remove writing to a file
        # self.novel_intermediate_exons.to_csv(os.path.join(results_dir,'classified_as_intermediate_exon.tsv'), sep='\t', index=True)

        # # write out files with the candidates that have been classified as background
        # self.novel_background_regions = X.loc[(y_pred == self.background_region_class),:]
        # # FIXME: remove writing to a file
        # self.novel_background_regions.to_csv(os.path.join(results_dir,'classified_as_background_region.tsv'), sep='\t', index=True)

        # ---------------------------------------------------------------------        
        # write out dataframe with probabilities
        X_with_probabilies = X.merge(X.apply(self.classifier_dict[classifier].predict_proba, axis=1), left_index=True, right_index=True)
        X_with_probabilies.to_csv(os.path.join(results_dir, 'classified_with_probabilities.tsv'), sep='\t', index=True)

        # Select terminal exons with the highest probability
        # For each 5pSS we select the exon with the highest probability to be terminal
        # This means that for each 5pSS we select only one poly(A) site
        X_tmp = X_with_probabilies.copy()
        X_tmp.reset_index(inplace=True)
        X_tmp["chromosome"], X_tmp["start"], X_tmp["end"], X_tmp["strand"] = X_tmp["Region"].str.split(':',3).str

        self.selected_novel_terminal_exons = pd.DataFrame()

        # group by strand
        for strand_group in X_tmp[X_tmp["classification"]=="terminal"].groupby(["strand"]):

            if strand_group[0] == "+":
                # group by chromosome, start and gene id and concatenate to the final dataframe
                for geneid_chromosome_start in strand_group[1].groupby(["chromosome", "start", "GeneId"]):
                    current_terminal = geneid_chromosome_start[1].loc[[geneid_chromosome_start[1]["terminal_probability"].argmax()]]
                    self.selected_novel_terminal_exons = pd.concat([self.selected_novel_terminal_exons, current_terminal])

            if strand_group[0] == "-":
                # group by chromosome, end and gene id and concatenate to the final dataframe
                for geneid_chromosome_end in strand_group[1].groupby(["chromosome", "end", "GeneId"]):
                    current_terminal = geneid_chromosome_end[1].loc[[geneid_chromosome_end[1]["terminal_probability"].argmax()]]
                    self.selected_novel_terminal_exons = pd.concat([self.selected_novel_terminal_exons, current_terminal])

        # Write out the final terminal exons that we will use
        if not self.selected_novel_terminal_exons.empty:
            self.selected_novel_terminal_exons.set_index(["Region", "GeneId"], inplace=True)
            self.selected_novel_terminal_exons.to_csv(os.path.join(results_dir, 'classified_as_terminal_with_probabilities.tsv'), sep='\t', index=True)
        else:
            sys.stderr.write("[WARNING] No novel terminal exons were detected\n")
            self.selected_novel_terminal_exons = pd.DataFrame(columns=["Region","GeneId"])
            self.selected_novel_terminal_exons.to_csv(os.path.join(results_dir, 'classified_as_terminal_with_probabilities.tsv'), sep='\t', index=False)



    def filter_training_data(raw_data_file_path, rownames_col="Region", profile_col="profile", min_feature_reads=5):
        
        """Filters raw data."""
        
        # read in the file
        df_raw = pd.io.parsers.read_csv(raw_data_file_path, sep="\t", index_col=rownames_col, header=0)

        # select the terminal exons we are interested in
        keep_rows_idx = (df_raw.ix[:,df_raw.columns != profile_col] >= min_feature_reads).any(axis=1)
        df_filtered = df_raw.loc[keep_rows_idx,:]

        # how many are we left with?
        print("Number of data sets before filtering:\t" + str(df_raw.shape[0]))
        print("Number of data sets after filtering:\t" + str(df_filtered.shape[0]))
        
        # return the filtered data
        return(df_filtered)

