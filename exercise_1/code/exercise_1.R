print_point <- function(txt) {
    open <- "\n\n\n>>>>>>>>>>>>>>>>>>>>>>>>>"
    close <- "<<<<<<<<<<<<<<<<<<<<<<<<\n\n\n"
    cat(open, " ", txt, " ", close)
}
print_task <- function(txt) {
    open <- "\n\n##########"
    close <- "##########\n\n"
    cat(open, " ", txt, " ", close)
}

tryCatch(
    { # nolint
        setwd("exercise_1/code")
    },
    error = function(cond) {

    }
)

############## Point 3

print_point("Point 3")

# install VariantAnnotation if not installed
print_task("install VariantAnnotation if not installed")

if (!requireNamespace("BiocManager", quietly = TRUE, verbose = FALSE)) {
    install.packages("BiocManager")
}

# BiocManager::install("VariantAnnotation")
# BiocManager::install("snpStats")
# BiocManager::install("GenomicRanges")
# BiocManager::install("biomaRt")

# Load the ’VariantAnnotation’ package which will also be used later on in the script
print_task("Load the ’VariantAnnotation’ package which will also be used later on in the script")

library(VariantAnnotation, quietly = TRUE)

# Get a list of the functions within the package
print_task("Get a list of the functions within the package")

available_funs <- lsf.str("package:VariantAnnotation")
print(available_funs)

# Check out the documentation for the function "readVcf"
print_task("Check out the documentation for the function \"readVcf\"")

?readVcf
help("readVcf")

# Get your current working directory
print_task("Get your current working directory")

cwd <- getwd()
print(cwd)

# Assign the sum of 2,3 and 4 to variable x
print_task("Assign the sum of 2,3 and 4 to variable x")

x <- 2 + 3 + 4
print(x)
x <- sum(c(2, 3, 4))
print(x)
x <- sum(2:4)
print(x)

# Make a character vector of the gene names PAX6, ZIC2, OCT4 and SOX2 and a second numeric countvector of the same length containing randomly sampled numbers between 1 and 10 (set the seed 42)
print_task("Make a character vector of the gene names PAX6, ZIC2, OCT4 and SOX2 and a second numeric countvector of the same length containing randomly sampled numbers between 1 and 10 (set the seed 42)")

genes <- c("PAX6", "ZIC2", "OCT4", "SOX2")
set.seed(24)
counts <- sample(1:10, length(genes))
print(genes)
print(counts)

# Subset the gene - vector using [] notation, and get the 2nd and 4th element
print_task("Subset the gene - vector using [] notation, and get the 2nd and 4th element")

subset <- genes[c(2, 4)]
print(subset)

# Generate a dataframe out of the two generated vectors
print_task("Generate a dataframe out of the two generated vectors")

genes_df <- data.frame(genes, counts)
print(genes_df)

# Select the genes of the generated dataframe for which the corresponding count value is greater than 5
print_task("Select the genes of the generated dataframe for which the corresponding count value is greater than 5")

selection <- genes_df[genes_df$counts > 5, ]
print(selection)

# Make a boxplot of the distribution of the generated count values in the dataframe
print_task("Make a boxplot of the distribution of the generated count values in the dataframe")

# boxplot(genes_df$counts)

# Write a function that takes a gene name and a dataframe as input with one column named
# "genes" and searches whether this gene name occurs in the column. It returns TRUE in case the
# genename occurs in the dataset and FALSE in case it doesn’t. Test the function with the above generated dataset.
print_task("Write a function that takes a gene name and a dataframe as input with one column named \"genes\" and searches whether this gene name occurs in the column. It returns TRUE in case the genename occurs in thedataset and FALSE in case it doesn’t. Test the function with the above generated dataset.")

is_gene_present <- function(gene_name, genes_df) {
    return(gene_name %in% genes_df$genes)
}
print(
    is_gene_present("BURUNDU", genes_df)
)
print(
    is_gene_present("OCT4", genes_df)
)

############## Point 4

print_point("Point 4")

vcf <- readVcf("data/e-geuv-1_filtered.vcf")

# Get the number of samples and variants.
print_task("Get the number of samples and variants.")

tmp <- dim(vcf)
n_samples <- tmp[1]
n_variants <- tmp[2]
print(n_samples)
print(n_variants)

# Get the first 5 SNPs from the first 3 sample
print_task("Get the first 5 SNPs from the first 3 sample")

gt <- geno(vcf)$GT[1:3, 1:5]
print(gt)

# Get the reference and alternative alleles for these first 5 SNPs.
print_task("Get the reference and alternative alleles for these first 5 SNPs.")

r_alleles <- ref(vcf)[1:3, 1:5]
print(r_alleles)

alt_alleles <- alt(vcf)[1:3, 1:5]
print(alt_alleles)

# Get the genotypes of samples HG00351, HG00353 and HG00355 for the variant rs17042098
print_task("Get the genotypes of samples HG00351, HG00353 and HG00355 for the variant rs17042098")

tmp <- geno(vcf)$GT["rs17042098", c("HG00351", "HG00353", "HG00355")]
print(tmp)

# Get the frequencies of the genotypes for SNP rs17042098
print_task("Get the frequencies of the genotypes for SNP rs17042098")

freqs <- table(geno(vcf)$GT["rs17042098", ])
print(freqs)

# Convert the genotypes 0/0, 0/1, 1/1 for SNP rs17042098 to 0, 1, 2
print_task("Convert the genotypes 0/0, 0/1, 1/1 for SNP rs17042098 to 0, 1, 2")

print(head(as(genotypeToSnpMatrix(vcf["rs17042098", ])$genotypes, "numeric")))

# Get all SNPs in the region of chromsome 4, 2 000 000 - 3 000 000 bp.

grange_search <- GRanges(seqnames = "4", ranges = IRanges(start = 2000000, width = 3000000))
grange <- rowRanges(vcf)
overlaps <- findOverlaps(grange, grange_search)
grange[overlaps@from]

############## Point 6

print_point("Point 6")

# Use this to look up the position (chromosome_name,start_position,end_position)
# and HGNC gene symbols forthe genes with Ensembl gene ids (ensembl_gene_id)
# ENSG00000196620, ENSG00000109787, ENSG00000241163and ENSG00000000938
print_task("Use this to look up the position (chromosome_name,start_position,end_position) and HGNC gene symbols forthe genes with Ensembl gene ids (ensembl_gene_id) ENSG00000196620, ENSG00000109787, ENSG00000241163and ENSG00000000938")

library(biomaRt, quietly = TRUE)

# get the basic ensembl annotations based on GRCh37
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
# define genes of interest with ensembl ids
genes <- c("ENSG00000196620", "ENSG00000109787", "ENSG00000241163", "ENSG00000000938")
# filter for desired gene ids
annotations <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes, mart = ensembl)
print("DDD")
print(annotations)