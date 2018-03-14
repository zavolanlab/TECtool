#!/usr/bin/env Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Jun 8, 2017
### Updates: March, 14 2018
### Author: Foivos Gypas
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Requirements: Gviz, rtracklayer, biomaRt, optparse, GenomicFeatures
### R version used:
#==================#
### Description: Plots novel terminal exons as identified by TECtool
### Output: A directory with novel plots
#==================#
#    HEADER END    #
#==================#

#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD LIBRARIES <---#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("Gviz"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("GenomicFeatures"))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
	make_option(c("--gtf"), action="store", type="character", default="", help="REQUIRED: GTF file with annotated transcripts", metavar="file"),
	make_option(c("--polyasites"), action="store", type="character", default="", help="REQUIRED: BED file with polya sites", metavar="file"),
	make_option(c("--bam"), action="store", type="character", default="", help="REQUIRED: Alignment file in BAM format", metavar="file"),
	make_option(c("--tectool_exons"), action="store", type="character", default="", help="REQUIRED: TECtool exons file in tsv format. Output of TECtool (classified_as_terminal_with_probabilities.tsv).", metavar="file"),
	make_option(c("--output_dir"), action="store", type="character", default="", help="REQUIRED: Output directory", metavar="directory"),
	make_option(c("--help"), action="store_true", default=FALSE, help="Show this information and die"),
	make_option(c("--verbose"), action="store_true", default=FALSE, help="Be Verbose")
)
#	make_option(c("-i", "--genome_version", action="store", type="character", default="hg38", help="REQUIRED: Genome Version. Example hg38, hg19, mm10"))
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf [FILE] --polyasites [FILE] --tectool_exons [FILE] --bam [FILE] --output_dir [DIRECTORY]", option_list = option_list, add_help_option=FALSE, description="\n Plots novel terminal exons as identified by TECtool")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if ( opt$gtf== "" || opt$polyasites=="" || opt$bam=="" || opt$tectool_exons==""  || opt$output_dir=="") {
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
	stop(print_help(opt_parser))
}

#==========================#
#    PRE-REQUISITES END    #
#==========================#

#================#
#   MAIN START   #
#================#

#---> START MESSAGE <---#
if ( opt$verbose  ) cat("Starting script'", "'...\n\n", sep="")

#---> Configuration for non ucsc chromosomes (chomosome that start with chr) <---
options(ucscChromosomeNames=FALSE)

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose  ) cat("Reading annotation file '", basename(opt$gtf), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object
# In BioC 3.1, all the makeTranscriptDbFrom*() functions were renamed makeTxDbFrom*(). The old names are still working but are now deprecated.
#txdb <- makeTranscriptDbFromGFF(opt$gtf)
txdb <- makeTxDbFromGFF(opt$gtf)
geneTrack <- GeneRegionTrack(txdb)

#---> IMPORT Polya sites <---#
# Print status message
if ( opt$verbose  ) cat("Reading bed file '", basename(opt$polyasites), "'...\n", sep="")
polyAtrack <- AnnotationTrack(opt$polyasites, shape="box")

#---> IMPORT Polya sites <---#
if ( opt$verbose   ) cat("Reading alignment file '", basename(opt$bam), "'...\n", sep="")
bam <- AlignmentsTrack(opt$bam, isPaired=FALSE) #, genome=opt$genome_version)

#---> Create output directory <---#
if ( opt$verbose    ) cat("Creating output directory'", basename(opt$output_dir), "'...\n", sep="")
dir.create(opt$output_dir, showWarnings = FALSE)

#---> IMPORT novel exons file <---#
# Print status message
if ( opt$verbose  ) cat("Reading tsv file of novel exons'", basename(opt$tectool_exons), "'...\n", sep="")
tectool_exons <- read.table(opt$tectool_exons, stringsAsFactors=FALSE, header=TRUE)

if ("terminal_probability" %in% colnames(tectool_exons)){

	# Use rtracklayer::import method to import GTF file to GRanges object 
	gr <- import(con=opt$gtf, format="gtf")
	# Subset EXONS (discards all other categories, e.g. CDS, start_codon etc.)
	gr <- gr[values(gr)[["type"]] == "exon"]
	# keep genes that we are interested in
	gr_subset <- gr[gr$gene_id %in% tectool_exons[["gene_id"]]]
	# create a granges object from exons
	tectool_exons_granges <- makeGRangesFromDataFrame(tectool_exons)

	# find upstream and downstream regions
	for(i in seq_along(tectool_exons_granges)) {
		
		# http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/nearest-methods.html

		# find previous exon
		precede_exon <- gr_subset[follow(tectool_exons_granges[i], gr_subset, select=c("arbitrary"), ignore.strand=FALSE)]
		precede_start <- start(precede_exon)
		precede_end <- end(precede_exon)
		# find next exon
		follow_exon <- gr_subset[precede(tectool_exons_granges[i], gr_subset, select=c("arbitrary"), ignore.strand=FALSE)]
		follow_start <- start(follow_exon)
		follow_end <- end(follow_exon)
		min_start <- min(precede_start, precede_end, follow_start, follow_end)
		max_end <- max(precede_start, precede_end, follow_start, follow_end)

		tectool_exons[i, "min_start"] <- min_start
		tectool_exons[i, "min_end"] <- max_end
	}
}

plot_novel_exon <- function(x) {

	chromosome <- x["chromosome"]
	# length <- abs(strtoi(x["end"]) - strtoi(x["start"]))
	# norm <- round(0.1*length)
	# start <- strtoi(x["start"]) - norm
	# end <- strtoi(x["end"]) + norm
	start <- strtoi(x["start"])
	end <- strtoi(x["end"])

	min_start <- strtoi(x["min_start"])
	min_end <- strtoi(x["min_end"])

	ht <- HighlightTrack(
		trackList = list(geneTrack, polyAtrack, bam),
		start = c(start),
		end=c(end),
		chromosome = chromosome
	)

	plotTracks(c(ht),
			   from=min_start,
			   to=min_end,
			   chromosome=chromosome,
			   type=c("coverage","sashimi"),
			   sashimiNumbers=FALSE,
			   main=paste(x["chromosome"],
			   	          x["min_start"],
			   	          x["min_end"],
			   	          x["terminal_probability"],
			   	          sep=" ")
	)
}

if ("terminal_probability" %in% colnames(tectool_exons)) {
	# sort by terminal probability
	tectool_exons <- tectool_exons[order(-tectool_exons["terminal_probability"]),]
	if ( opt$verbose  ) cat("Generating plot for novel terminal exons")
	
	if ("chromosome" %in% colnames(tectool_exons)) {

		pdf(file.path(opt$output_dir,"plots.pdf"))
		for (i in 1:nrow(tectool_exons)){
				plot_novel_exon(tectool_exons[i,])
		}
		dev.off()
	}
} else{
	x <- data.frame()
	write.table(x, file=file.path(opt$output_dir,"plots.pdf"), col.names=FALSE)
}

#================#
#   MAIN END     #
#================#