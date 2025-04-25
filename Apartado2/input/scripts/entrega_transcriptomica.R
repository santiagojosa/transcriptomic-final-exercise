library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(tidyverse)


# Carga de datos ----------------------------------------------------------

rawcounts <- read_tsv("input/rawcounts.tsv", col_names = TRUE)
colnames(rawcounts) <- c("gene_id", colnames(rawcounts)[-1])
rawcounts <- column_to_rownames(rawcounts, "gene_id")

metadata <-  read_tsv("input/metadata.tsv", col_names = TRUE)
colnames(metadata) <- c("sample_id", colnames(metadata)[-1])
metadata$agent <- as.factor(metadata$agent)
metadata$time <- as.factor(metadata$time)
metadata$patient <- as.factor(metadata$patient)

# Creacion de nueva variable que agrupe time y agent
metadata <- metadata |> mutate(group = paste(agent, time, sep = "_"))
metadata$group <- as.factor(metadata$group)
# creacion de objeto SummarizedExperiment
summarized_experiment <- SummarizedExperiment(assays = list(counts = as.matrix(rawcounts)),
                                              colData = metadata)


# Exploracion de batch effect ---------------------------------------------

dds_exploratorio <- DESeqDataSet(summarized_experiment, design = ~group)  # aún sin modelo
keep <- rowSums(counts(dds_exploratorio)) >= 10 ## Seleccionar genes con más de 10 cuentas en todos los samples
dds_exploratorio <- dds_exploratorio[keep, ]

vsd_exploratorio <- vst(dds_exploratorio, blind = TRUE)

plotPCA(vsd_exploratorio, intgroup = "sample_id") +
  labs(color = "Sample ID") 
plotPCA(vsd_exploratorio, intgroup = "patient") +
  labs(color = "patient")
plotPCA(vsd_exploratorio, intgroup = "agent") +
  labs(color = "agent")
plotPCA(vsd_exploratorio, intgroup = "time") +
  labs(color = "time")
plotPCA(vsd_exploratorio, intgroup = "group") +
  labs(color = "group")
plotPCA(vsd_exploratorio, intgroup = c("patient", "agent")) + 
  labs(color = "patient:agent")

## Calculamos las distancias a partir de las cuentas normalizadas y variance-stabilized (vst)
sampleDists <- dist(t(assay(vsd_exploratorio))) # Calcula la distancia de unas muestras a otras

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_exploratorio$patient, vsd_exploratorio$time, vsd_exploratorio$agent, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# Creacion de objeto dds con modelo ---------------------------------------

dds <- DESeqDataSet(summarized_experiment, design = ~patient + group)

keep <- rowSums(counts(dds)) >= 10 ## Seleccionar genes con más de 10 cuentas en todos los samples
dds <- dds[keep, ]
vsd <- vst(dds, blind = TRUE)

dds <- DESeq(dds, test = "Wald")

plotDispEsts(dds)

plotMA(dds)

# DPNvsControl - 24h ----------------------------------------------

results_OHTvsControl <- results(object = dds,
                                contrast = c("group", "OHT_24h", "Control_24h"), # elegir las dos que nos interesa comparar. El control va el segundo.
                                pAdjustMethod = "BH",
                                tidy = TRUE)

genesID <-mygene::queryMany(results_OHTvsControl$row,
                            scopes="ensembl.gene",
                            fields="symbol",
                            species="human")
genesID <- genesID[!duplicated(genesID$query),]
results_OHTvsControl$row <- ifelse(is.na(genesID$symbol),
                                   genesID$query,
                                   genesID$symbol)

EnhancedVolcano(results_OHTvsControl,
                lab = results_OHTvsControl$row,
                x = "log2FoldChange",
                y = "padj",
                title = "DEG OHT vs Control",
                FCcutoff = 1,
                pCutoff = 0.05,
                subtitle = NULL,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                labSize = 3.0,
                pointSize = 2,
                legendPosition = "none",
                widthConnectors = 1)

# Selecciono genes significativos para heatmap
OHTvsControl_padj05_FC2 <-results_OHTvsControl |> filter(padj < 0.05, abs(log2FoldChange) > 1)
write_tsv(OHTvsControl_padj05_FC2, "out/tabs/DEG_OHTvsControl.tsv")

## Heatmap de los genes TOP DGE por p-valor ajustado
heatmap_results_OHTvsControl <- results(object = dds,
                                contrast = c("group", "OHT_24h", "Control_24h"), # elegir las dos que nos interesa comparar. El control va el segundo.
                                pAdjustMethod = "BH",
                                tidy = TRUE)
heatmap_results_OHTvsControl_padj05_FC2 <-heatmap_results_OHTvsControl |> filter(padj < 0.05, abs(log2FoldChange) > 1)
heatmap_genes_selected_OHTvsControl_padj05_FC2 <- heatmap_results_OHTvsControl_padj05_FC2 |> pull(row)

vsd |> colnames() |> sort()
mat <- assay(vsd)[heatmap_genes_selected_OHTvsControl_padj05_FC2, ] 
pheatmap(mat)


# DPNvsControl - 24h ------------------------------------------------------

results_DPNvsControl <- results(object = dds,
                                contrast = c("group", "DPN_24h", "Control_24h"), # elegir las dos que nos interesa comparar. El control va el segundo.
                                pAdjustMethod = "BH",
                                tidy = TRUE)

genesID <-mygene::queryMany(results_DPNvsControl$row,
                            scopes="ensembl.gene",
                            fields="symbol",
                            species="human")
genesID <- genesID[!duplicated(genesID$query),]
results_DPNvsControl$row <- ifelse(is.na(genesID$symbol),
                                   genesID$query,
                                   genesID$symbol)

EnhancedVolcano(results_DPNvsControl,
                lab = results_DPNvsControl$row,
                x = "log2FoldChange",
                y = "padj",
                title = "DEG DPN vs Control",
                FCcutoff = 1,
                pCutoff = 0.05,
                subtitle = NULL,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                labSize = 3.0,
                pointSize = 2,
                legendPosition = "none",
                widthConnectors = 1)

# Selecciono genes significativos para heatmap
DPNvsControl_padj05_FC2 <-results_DPNvsControl |> filter(padj < 0.05, abs(log2FoldChange) > 1)
write_tsv(DPNvsControl_padj05_FC2, "out/tabs/DEG_DPNvsControl.tsv")

## Heatmap de los genes TOP DGE por p-valor ajustado
heatmap_results_DPNvsControl <- results(object = dds,
                                        contrast = c("group", "DPN_24h", "Control_24h"), # elegir las dos que nos interesa comparar. El control va el segundo.
                                        pAdjustMethod = "BH",
                                        tidy = TRUE)
heatmap_results_DPNvsControl_padj05_FC2 <-heatmap_results_DPNvsControl |> filter(padj < 0.05, abs(log2FoldChange) > 1)
heatmap_genes_selected_DPNvsControl_padj05_FC2 <- heatmap_results_DPNvsControl_padj05_FC2 |> pull(row)

vsd |> colnames() |> sort()
mat <- assay(vsd)[heatmap_genes_selected_DPNvsControl_padj05_FC2, ] 
pheatmap(mat)


# Obtener tabla para GSEA
GSEAresults_DPNvsControl <- results(object = dds,
                                contrast = c("group", "DPN_24h", "Control_24h"), # elegir las dos que nos interesa comparar. El control va el segundo.
                                pAdjustMethod = "BH",
                                tidy = TRUE)

datos_GSEA <- GSEAresults_DPNvsControl |> select(row, log2FoldChange) |> arrange(desc(log2FoldChange))
write_tsv(datos_GSEA, "out/tabs/DEA.rnk", col_names = FALSE)
