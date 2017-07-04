#'
#'
#' Quantify interactions between CNEs and target/bystander genes within GRBs
#' 
#'

# source("https://bioconductor.org/biocLite.R")
# biocLite(c("biomaRt", "rtracklayer", "InteractionSet", "TxDb.Hsapiens.UCSC.hg19.knownGene"))

require(biomaRt)        # to retrieve human paralogs from Ensembl
require(rtracklayer)    # to parse .bed files
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(chromloop)      # devtools::install_github("ibn-salem/chromloop")
require(RColorBrewer)
require(tidyverse)      # for tidy data

#-------------------------------------------------------------------------------
# set some parameters:
#-------------------------------------------------------------------------------

GENES_FILE <- "data/GRB_Interactions/targets_bystanders_hg38.tsv"
GRB_FILE <- "data/GRB_Interactions/grbsHg38.bed.hg19.bed"
CNE_FILE <- "data/GRB_Interactions/human_mouse_cnes_hg38.bed.hg19.bed"

CaptureHiC_Files <- c(
  "data/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt",
  "data/Mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt"
)

COL_GENE_GROUP <- brewer.pal(8, "Set1")[c(1, 5)]

outPrefix <- "results/GRB_interactions"
# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# parse Capture Hi-C data
#-------------------------------------------------------------------------------
seqInfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)

chicList <- lapply(CaptureHiC_Files,
       chromloop::parseCaptureHiC, seqinfo = seqInfo)

captureHiC <- do.call("c", chicList)

#-------------------------------------------------------------------------------
# parse GRBs and CNEs
#-------------------------------------------------------------------------------
grbGR <- import.bed(GRB_FILE, seqinfo = seqInfo)
cneGR <- import.bed(CNE_FILE, seqinfo = seqInfo)

#-------------------------------------------------------------------------------
# get target and bystander gene in hg19 as GenomicRanges object
#-------------------------------------------------------------------------------
genesDF <- read_tsv(GENES_FILE, col_types = cols(
  ENSEMBLGeneID = col_character(),
  ucscCoordinates = col_character(),
  geneName = col_character(),
  originalSupport = col_double(),
  normalisedSupport = col_double(),
  prediction = col_character()
))


#-------------------------------------------------------------------------------
# get genes for ENSG
#-------------------------------------------------------------------------------
ensemblGRCh37 <- useMart(host = "grch37.ensembl.org", 
                         biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl", 
                         verbose = FALSE)

geneAttributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", 
                   "start_position", "end_position", "strand", "status", 
                   "gene_biotype")

genes = as_tibble(getBM(
  attributes = geneAttributes, 
  mart = ensemblGRCh37, 
  filters = list("ensembl_gene_id" = genesDF$ENSEMBLGeneID, 
                 "chromosome_name" = c(1:22, "X", "Y"))
  ))

# add target gene prediction to downloaded ensambl data frame
genes <- genes %>%
  left_join(genesDF, by = c("ensembl_gene_id" = "ENSEMBLGeneID"))

# make GRanges object for tss of all known prot coding genes
tssCoord <- ifelse(genes$strand == 1, genes$start_position, genes$end_position)
tssGR <- GRanges(paste0("chr", genes$chromosome_name),
        IRanges(tssCoord, tssCoord),
        strand = ifelse(genes$strand == 1, '+', '-'), 
        names = genes$ensembl_gene_id, 
        select(genes, -strand),
        seqinfo = seqInfo
        )

#-------------------------------------------------------------------------------
# map CNEs to GRBS and GRBs to genes
#-------------------------------------------------------------------------------

#' Goal is a dataset with following columns
#' gene (Idx)
#' GRB (IDx)
#' CNE (IDx)

geneToGRB <- as_tibble(as.data.frame(findOverlaps(tssGR, grbGR))) %>% 
  rename(gene_ID = queryHits, grb_ID = subjectHits)

GrbToCne <- as_tibble(as.data.frame(findOverlaps(grbGR, cneGR))) %>% 
  rename(grb_ID = queryHits, cne_ID = subjectHits)

gene_cne <- genes %>% 
  dplyr::select(-chromosome_name, -start_position, -end_position, -strand, 
                -status, -gene_biotype, -ucscCoordinates) %>% 
  mutate(gene_ID = 1:nrow(genes)) %>% 
  left_join(geneToGRB, by = "gene_ID") %>% 
  left_join(GrbToCne, by = "grb_ID") %>% 
  filter(!is.na(gene_ID) & !is.na(cne_ID)) %>% 
  mutate(gene_cne_ID = 1:n())

#-------------------------------------------------------------------------------
# interactions between genes and CNEs
#-------------------------------------------------------------------------------

# build GI of Gene-CNE pairs as query
gene_cne_GI <- GInteractions(tssGR[gene_cne$gene_ID], cneGR[gene_cne$cne_ID]) 

captureHiC_DF <- as_tibble(as.data.frame(mcols(captureHiC))) %>% 
  mutate(CHiC_ID = 1:length(captureHiC))

# get mapping of Gene-CNE pairs to supporting (overlapping) interactions
gene_cne_interactions <- as_tibble(data.frame(
  findOverlaps(gene_cne_GI, captureHiC)
  )) %>% 
  rename(gene_cne_ID = queryHits, CHiC_ID = subjectHits) %>% 
  left_join(captureHiC_DF, by = "CHiC_ID") %>% 
  group_by(gene_cne_ID) %>% 
  summarise(
    raw_count_mean = mean(raw_count),
    log_observed_expected_mean = mean(log_observed_expected),
    interactions = n()
  )


#-------------------------------------------------------------------------------
# annotate gene_cne data.frame with interactions
#-------------------------------------------------------------------------------
gene_cne <- gene_cne %>%
  left_join(gene_cne_interactions, by = "gene_cne_ID")

# save result table
write_tsv(gene_cne, paste0(outPrefix, ".gene_cne.tsv"))
save(gene_cne, file = paste0(outPrefix, ".gene_cne.Rdata"))


#===============================================================================
# Analyse interactions between target and bystander genes
#===============================================================================

#-------------------------------------------------------------------------------
# Number of CNEs per gene for target vs. bystander
#-------------------------------------------------------------------------------
countDF <- gene_cne %>% 
  dplyr::group_by(ensembl_gene_id, prediction) %>% 
  dplyr::summarize(n_cne = n())

p <- ggplot(countDF, aes(x = n_cne, y = ..count.., fill = prediction, color = prediction)) + 
  geom_density(alpha = .2) + 
  scale_fill_manual(values = COL_GENE_GROUP) +
  scale_color_manual(values = COL_GENE_GROUP) +
  labs(x = "Number of CNEs per gene")

#-------------------------------------------------------------------------------
# Number of CNEs per gene for target vs. bystander
#-------------------------------------------------------------------------------
