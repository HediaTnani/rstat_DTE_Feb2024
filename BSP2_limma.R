# Rstat club 02/02/2024
# Load the libraries
library(SummarizedExperiment)
library(limma)
library(purrr)
library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(edgeR)
library(tximport)

# Load the data
load("/dcs04/lieber/lcolladotor/chessBrain_LIBD4085/processed_data/rda/rse_tx.bsp2_dlpfc_gencode25.txi.n500.Rdata", verbose = TRUE)

# Filtering the data
set.seed(123)
cutoff <- jaffelab::expression_cutoff(assays(rse_tx)$tpm,max_cut = 1)
# the suggested expression cutoff is 0.29
# Subset the rse_tx object by the expression cutoff
rse_tx<-rse_tx[rowMeans(assays(rse_tx)$tpm)>=cutoff,]

# Check the assays
assays(rse_tx)
# List of length 3
# names(3): counts tpm length

# Note: Transcripts Per Million (TPM) is a normalization method for RNA-seq data that facilitates comparisons both within and between samples. It adjusts for gene length and sequencing depth, making it a helpful measure for relative gene expression

# Remove the "|" from the rownames
rownames(rse_tx) <- rowData(rse_tx)$transcript_id

## Subset the rse_tx object by DLPFC region
rse_tx <- rse_tx[, rse_tx$Region == "DLPFC"]

# Creates a model matrix for regression analysis using the predictors Dx, Age, Sex and Race from the column data of the rse_tx object.
rse_tx$Dx <- as.factor(rse_tx$Dx)
mod <- model.matrix(~ Dx + Age + Sex + Race,
                    data = colData(rse_tx)
)

# Fit the linear model
## Extract the transcript expression values and put them in the
## log2(TPM + 0.5) scale
txExprs = (assay(rse_tx,"tpm")+0.5) %>% log2()

## Run the standard linear model for differential expression
fitTx <- lmFit(txExprs, mod)
eBTx <- eBayes(fitTx)
## Extract the differential expression results
sigTx_logtpm <- topTable(eBTx,
                        coef = "DxSCZD", # test for differential expression with respect to the "DxSchizo" coefficient.
                        p.value = 1, # all genes will be returned, regardless of their p-values
                        number = nrow(rse_tx) # all genes to be returned
)

#######################################################
# Trend=T model
arrayw  <- arrayWeights(txExprs, mod)  # calculate weights
fitTx_trend   <- lmFit(txExprs, mod, weights=arrayw) # fit the model
# Bayes shrinkage
eBTx_trend <- eBayes(fitTx_trend, trend = TRUE)
## Extract the differential expression results
sigTx_trend <- topTable(eBTx_trend,
                  coef = "DxSCZD", # test for differential expression with respect to the "DxSchizo" coefficient.
                  p.value = 1, # all genes will be returned, regardless of their p-values
                  number = nrow(rse_tx) # all genes to be returned
)

#######################################################
# Volcano plot
# Create the two volcano plots
p_logtpm <- EnhancedVolcano(sigTx_logtpm,
                            lab = rownames(sigTx_logtpm),
                            x = 'logFC',
                            y = 'adj.P.Val',
                            pCutoff = 0.05,
                            selectLab = FALSE,
                            FCcutoff = 0.3,
                            xlim = c(-1.2, 1.2),
                            ylim = c(0, 4.3),
                            title = 'logtpm model')

p_trend <- EnhancedVolcano(sigTx_trend,
                                  lab = rownames(sigTx_trend),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  pCutoff = 0.05,
                                  FCcutoff = 0.3,
                                  selectLab = FALSE,
                                  xlim = c(-1.2, 1.2),
                                  ylim = c(0, 4.3),
                                  title = 'Trend model')

# Combine the plots
combined_plot <- p_logtpm + p_trend

# Print the combined plot
combined_plot

############ top1500 ############
# Function to get the top n significant transcripts
topn_sig <- function(data, significance_level = 0.05, top_n) {
  data %>%
    tibble::rownames_to_column("Transcripts") %>%
    dplyr::select(Transcripts, logFC, t, P.Value, adj.P.Val) %>%
    dplyr::filter(adj.P.Val <= significance_level) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice(1:top_n)
}

# Get the top 1500 significant transcripts
sigTx_top1500 <- topn_sig(sigTx_logtpm, significance_level = 0.05, top_n = 1500)
sigTx_trend_top1500 <- topn_sig(sigTx_trend, significance_level = 0.05, top_n = 1500)

# Prepare the data for the Venn diagram
x <- list(
  sigTx_top1500 = sigTx_top1500$Transcripts, 
  sigTx_trend_top1500 = sigTx_trend_top1500$Transcripts
)
# Create the Venn diagram
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

##########################
# Limma-voom
# Function to run the limma-voom
limma_voom <- function(rse) {
  scaledCounts <- tximport::makeCountsFromAbundance(assays(rse)$counts, assays(rse)$tpm,assays(rse)$length, countsFromAbundance = "lengthScaledTPM")
  y <- DGEList(scaledCounts)
  ## model design here:
  design <- model.matrix(~Dx + Age + Sex + Race, data=colData(rse_tx))
  y <- y[filterByExpr(y, design), ]
  y <- calcNormFactors(y)
  v <- voom(y, design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  topTable(fit, coef = "DxSCZD", p.value = 1, number = nrow(rse_tx))
} 

# Run the limma-voom
# Remove the "|" from the rownames
rownames(rse_tx) <- rowData(rse_tx)$transcript_id
sigTx_limma_voom <- limma_voom(rse_tx)

# Volcano plot
p_voom <- EnhancedVolcano(sigTx_limma_voom,
                                  lab = rownames(sigTx_limma_voom),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  pCutoff = 0.05,
                                  FCcutoff = 0.3,
                                  selectLab = FALSE,
                                  xlim = c(-1.2, 1.2),
                                  ylim = c(0, 4),
                                  title = 'Voom model')

# Combine the plots
combined_plot <- p_logtpm + p_voom

# Get the top1500 significant transcripts
sigTx_top1500_voom <- topn_sig(sigTx_limma_voom, significance_level = 0.05, top_n = 1500)

# Prepare the data for the Venn diagram
x <- list(
  limma_logtpm_top1500 = sigTx_top1500$Transcripts, 
  limma_trend_top1500 = sigTx_trend_top1500$Transcripts,
  limma_voom_top1500 = sigTx_top1500_voom$Transcripts
)
# Create the Venn diagram
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)
