#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("stringr")
library("dplyr")
library("readr")
library("VariantAnnotation")
library("parallel")
library("stringr")
library("data.table")

cols_tsv <- 
cols(
  .default = col_character(),
  case_id = col_character(),
  normal_id = col_character(),
  tumour_id = col_character(),
  chromosome = col_character(),
  start = col_double(),
  stop = col_double(),
  gene = col_character(),
  gene_id = col_character(),
  type = col_character(),
  filter = col_character(),
  ref = col_character(),
  alt = col_character(),
  gt = col_character(),
  pl = col_character(),
  mut_pr = col_double(),
  tr = col_double(),
  ta = col_double(),
  nr = col_double(),
  na = col_double(),
  dbsnp = col_character(),
  thousand_genomes = col_character(),
  cosmic = col_character(),
  caller = col_character(),
  trinucleotide_ref = col_character(),
  trinucleotide_alt = col_character(),
  amino_acid_change = col_character(),
  functional_class = col_character(),
  gene_coding = col_character(),
  project = col_character(),
  substitution = col_character(),
  impact = col_character(),
  dna_change = col_character(),
  biotype = col_character()
)

cols_pair_strelka <- 
cols(
  case_id = col_character(),
  normal_id = col_character(),
  tumour_id = col_character(),
  chromosome = col_double(),
  start = col_double(),
  stop = col_double(),
  gene = col_character(),
  gene_id = col_character(),
  type = col_character(),
  ref = col_character(),
  alt = col_character(),
  tr = col_double(),
  ta = col_double(),
  nr = col_double(),
  na = col_double(),
  dbsnp = col_character(),
  thousand_genomes = col_logical(),
  cosmic = col_logical(),
  caller = col_character(),
  amino_acid_change = col_character(),
  functional_class = col_character(),
  gene_coding = col_character(),
  project = col_character(),
  dna_change = col_character(),
  biotype = col_character()
)


#args <- c("--random", "--targetfile=/scratch/shahlab_tmp/danlai/APARICIO-590/PAIR_STRE/strelka_pipeline/OUTPUT/RUN/SA1259NC_strelka/outputs/results/TASK_12_STRELKA_PARSE_INDL_strelka_indels_parsed.tsv", "--outputfile=testfile.csv")

#args <- c("--random", "--targetfile=/scratch/shahlab_tmp/danlai/APARICIO-590/PAIR_STRE/strelka_pipeline/OUTPUT/RUN/SA1259NC_strelka/outputs/results/TASK_12_STRELKA_PARSE_INDL_strelka_indels_parsed.tsv", "--outputfile=testfile.csv")


target_file <- args[str_detect(args,"--targetfile")] %>% str_replace("--targetfile=","")
#!str_detect(paste(target_file), paste(getwd()))
output_file <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")
target_file_type <-gsub( "vcf.gz","vcf",basename(target_file)) %>% str_split("\\.") %>% unlist() %>% last()


#target_file_2 <- "/scratch/shahlab_tmp/danlai/APARICIO-590/SING_MUTA/mutationseq_pipeline/OUTPUT/RUN/SA1228N_museq/outputs/TASK_8_PARSE_museq_parsed_premapp.tsv"

#import referenc sequences
gencode_reference_col_types = cols(
  X1 = col_character(),
  X2 = col_character(),
  X3 = col_character(),
  X4 = col_double(),
  X5 = col_double(),
  X6 = col_character(),
  X7 = col_character(),
  X8 = col_character(),
  X9 = col_character()
)

gene_from_add_info <- function(x){
		output <- NA
		add_info_type <- unlist(strsplit(x, ";")) %>% trimws()
		matches <- str_detect(add_info_type, "gene_name")
		if(length(matches == 1)){
			output <- add_info_type[str_detect(add_info_type, "gene_name")] %>% str_replace("gene_name ","")
		}
		output 
	} 

load_parsed_reference_data <- TRUE
if(load_parsed_reference_data == FALSE){
	#gencode_reference <- read_tsv("/shahlab/archive/misc/sbeatty/reference/gencode.v19.annotation.gtf_withproteinids", skip=5, col_names=FALSE, col_types=gencode_reference_col_types)
	#names(gencode_reference) <- c("chr", "source", "feature_type","start","stop", "score", "strand", "phase", "add_info")
	
	clinvar <- read_tsv("/shahlab/archive/misc/sbeatty/reference/variant_summary.txt")
	clinvar_GRCh37 <- clinvar[c(clinvar[,"Assembly"] == "GRCh37"),] ; rm(clinvar)
	
	
	
	gencode_reference <- read_tsv("/shahlab/archive/misc/sbeatty/reference/gencode.v19.annotation.gtf_withproteinids", skip=5, col_names=FALSE,quote="XXX", col_types=gencode_reference_col_types)
	names(gencode_reference) <- c("Chromosome", "source", "feature_type","start","stop", "score", "strand", "phase", "add_info")
	gencode_reference$Chromosome <- gsub("chr","", gencode_reference$Chromosome)
	gencode_reference$add_info <- gsub(regex('\\"'),"",gencode_reference$add_info)
	gencode_reference$genecode_gene_name <-mcmapply(gencode_reference$add_info, FUN=gene_from_add_info, USE.NAMES = FALSE, mc.cores=20)



clinvar <- read_tsv("/shahlab/archive/misc/sbeatty/reference/variant_summary.txt")
	clinvar_GRCh37 <- clinvar[c(clinvar[,"Assembly"] == "GRCh37"),] ; rm(clinvar)
	#eliminate annotations for missing a chromosome number 
	clinvar_GRCh37 <- clinvar_GRCh37[!is.na(clinvar_GRCh37$Chromosome),]
	#eliminate rows stop position less than 0.
	clinvar_GRCh37 <- clinvar_GRCh37[!c(clinvar_GRCh37$Stop < 0),]
	#reduce clinvar to only GRCh37 reference
	clin_var_ranges <-  makeGRangesFromDataFrame(data.frame(clinvar_GRCh37), 
	       seqnames.field ="Chromosome" , start.field="Start", end.field="Stop", 
	    keep.extra.columns=TRUE)
	gencode_ranges <-  makeGRangesFromDataFrame(data.frame(gencode_reference), 
	       seqnames.field ="Chromosome" , start.field="start", end.field="stop", 
	    keep.extra.columns=TRUE)
	gencode_ranges$add_info <- gsub(regex('\\"'),"",gencode_ranges$add_info)
	gencode_ranges$genecode_gene_name <-mcmapply(gencode_ranges$add_info, FUN=gene_from_add_info, USE.NAMES = FALSE, mc.cores=20)

	save(gencode_ranges, file="//scratch/shahlab_tmp/sbeatty/ind231/reference_data/gencode_ranges.Rdata")
	save(clin_var_ranges,file="/scratch/shahlab_tmp/sbeatty/ind231/reference_data/clin_var_ranges.Rdata")
	save(clinvar_GRCh37, file="/scratch/shahlab_tmp/sbeatty/ind231/reference_data/clinvar_GRCh37.Rdata")
	save(gencode_reference,file="/scratch/shahlab_tmp/sbeatty/ind231/reference_data/gencode_reference.Rdata")

} else {
	load("//scratch/shahlab_tmp/sbeatty/ind231/reference_data/gencode_ranges.Rdata")
	load("//scratch/shahlab_tmp/sbeatty/ind231/reference_data/clin_var_ranges.Rdata")
	load("//scratch/shahlab_tmp/sbeatty/ind231/reference_data/clinvar_GRCh37.Rdata")
	load("//scratch/shahlab_tmp/sbeatty/ind231/reference_data/gencode_reference.Rdata")
}
get_gene_ids <- function(row, match_list, first_only=TRUE){
    output <- NA
    matched_results <- data.frame(gencode_ranges[match_list@to[which(match_list@from == row)],"add_info"])[,"add_info"]
    if(first_only == TRUE & length(matched_results) > 0){
        matched_results <- matched_results[1]
    }
    if(length(matched_results) > 0){
        output <- paste(mapply(matched_results, USE.NAMES=FALSE, FUN=function(x) {get_gencode(string=x, key="gene_id")}), collapse="; ")
    }
    output
}

genome_range_to_dataframe <- function(genome_range){
	df <- genome_range
	df_rowRanges <- rowRanges(df)
	range_start <- ranges(df)@start
	range_end <- ranges(df)@start
	range_width <- ranges(df)@width
	df_seqnames <-  as.character(seqnames(df))
	df_elementMetadata <- elementMetadata(rowRanges(df))
	strand_info <- as.character(df_rowRanges@strand)
	df_paramRangeID <- as.character(df_elementMetadata[,"paramRangeID"])
	df_REF <- as.character(df_elementMetadata[,"REF"])
	df_ALT <- rownames(df_elementMetadata) %>% str_sub(start=-1) %>% as.character
	df_QUAL <- as.character(df_elementMetadata[,"QUAL"])
	df_FILTER <- as.character(df_elementMetadata[,"FILTER"])
	df_out <- data.frame(chr=df_seqnames, pos=range_start, start = range_start,
		end=range_end, width=range_width, strand=strand_info,
		paramRangeID=df_paramRangeID, REF=df_REF, ALT=df_ALT, QUAL=df_QUAL, Filter=df_FILTER)
	df_out
}

if(target_file_type == "vcf"){
	df <- readVcf(paste(target_file), genome="hg19")
	ranges <- df
	annotation_df <- genome_range_to_dataframe(df)
	annotation_df[,"genecode_gene_name"] <- NA
	names(annotation_df) <- gsub("stop", "end", names(annotation_df))


} else if (target_file_type == "tsv") {
	if(str_detect(target_file,"mutationseq_pipeline")){
		df <- read_tsv(paste(target_file),  col_types = cols_tsv)
	} else if(str_detect(target_file,"PAIR_STRE") & str_detect(target_file,"indels")) {
		df <- read_tsv(paste(target_file),  col_types = cols_pair_strelka)
	}
	else {
		df <- read_tsv(paste(target_file))
	}

	df <- df[!is.na(df$chromosome),]
	ranges <- makeGRangesFromDataFrame(df, 
       seqnames.field ="chromosome" , start.field="start", end.field="stop", keep.extra.columns=TRUE)
	annotation_df <- df
	annotation_df <- data.frame(annotation_df)
	annotation_df[,"genecode_gene_name"] <- NA
	names(annotation_df) <- gsub("chromosome", "chr", names(annotation_df))
	names(annotation_df) <- gsub("stop", "end", names(annotation_df))

} else {
	print("error stuff is wack yo")
}


gencode_matches <- suppressWarnings(findOverlaps(ranges, gencode_ranges, type="any"))
annotation_df[,"genecode_gene_name"][gencode_matches@from] <- gencode_ranges$genecode_gene_name[gencode_matches@to]
annotation_df[,"unique_id"] <- paste0(annotation_df$chr,"@",annotation_df$end)
clin_var_matches <- findOverlaps(ranges, clin_var_ranges)

#subset clinvar to only those with matches
match_list_short <- clin_var_matches[match(unique( clin_var_matches@to), clin_var_matches@to),]
clin_var_relevant <- clin_var_ranges[match_list_short@to]
clin_var_matches_subset <- findOverlaps(ranges, clin_var_relevant)
clin_var_relevant$unique_id <- NA
clin_var_relevant$unique_id[clin_var_matches_subset@to] <- annotation_df$unique_id[clin_var_matches_subset@from]

joined_df <- left_join(data.frame(annotation_df), data.frame(clin_var_relevant), by="unique_id", suffix=c("museq","clinvar"), KEEP=TRUE)

fwrite(joined_df, file=paste(output_file))

