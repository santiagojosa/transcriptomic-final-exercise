# intialize ----------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(readxl)
library(VennDiagram)
library(magick)


# Este metodo usa la función DESeqDataSetFromMatrix a partir de una matriz de 
# rawcounts (reads por gene)

# Loading Data ----------------------------------------------------------------
# Primero se extrae un subset de la matriz de rawcounts con la info de metadatos
# Cargamos archivos con los rawcounts de RNAseq por gene para todos los PDX del laboratorio y metadatos de todos esos PDX
rawcounts_file <- read_tsv(file = "data/RNAseq_PDX_pase2_rawcounts_18nov24.tsv")
metadata_file <- read_tsv(file = "data/PDX_PASE2_metadata_corrected_18nov24.tsv")

# Descargo las tablas de surfaceome y las guardo en res/
surfaceome_all <- read_xlsx("res/table_S3_surfaceome.xlsx",
                            sheet = "SurfaceomeMasterTable",
                            skip = 1,
                            col_names = TRUE,
                            guess_max = 10000) # guess_max es mas alto que por defecto porque me detectaba mal las columnas
surfaceome_surface <- read_xlsx("res/table_S3_surfaceome.xlsx",
                                sheet = "in silico surfaceome only",
                                skip = 1,
                                col_names = TRUE,
                                guess_max = 10000)
surfaceome_ID <- read_table("res/surfaceome_ids.txt",
                            col_names = "ID")
# Compruebo que los IDs de las tablas son los mismos. Usaré la tabla más completa (archivo "table_S3_surfaceome", hoja "SurfaceomeMasterTable")
all(sort(surfaceome_surface$`UniProt name`) == sort(surfaceome_ID$ID))

# De la tabla de surfaceoma, selecciono las columnas de interes para el analisis posterior
surface_genes_info <- 
  surfaceome_surface |> 
  select(`UniProt gene`,
         GeneID,
         `Surfaceome Label`,
         `Membranome Almen main-class`,
         `Membranome Almen sub-class`, 
         `UniProt subcellular`,
         `COMPARTMENTS benchmark pos`,
         `COMPARTMENTS benchmark neg`)
# Tengo que limpiar tabla
# # # # # # # #

# SCLC_vs_ADC ----------------------------------------------------------------
### Comparación SCLC vs ADC all genes
# Selecciono los metadatos de mis muestras
metadata <- metadata_file |>
  filter(ID != "TP191") |> # quito el TP191 porque no está bien clasificado
  filter(Type == "SCLC" | Subtype == "ADC") |> # selecciono los SCLC y los ADC
  select(ID, Type, Subtype) |> 
  arrange(Type, Subtype) |> 
  column_to_rownames("ID")
metadata$Type <- factor(metadata$Type) # Convierto Type en un factor
metadata$Subtype <- factor(metadata$Subtype) # Convierto Subtype en un factor

# Indico las muestras a estudiar
# muestras <- c("CDX1",
#               "CDX5",
#               "CDX9",
#               "TP135",
#               "TP171",
#               "TP176",
#               "TP185",
#               "TP200",
#               "RTP202_10",
#               "TP103",
#               "TP118",
#               "TP126",
#               "TP134",
#               "TP143",
#               "TP148",
#               "TP166",
#               "TP181",
#               "TP184")

# Subset de las rawcounts de las muestras de interes
rawcounts <- rawcounts_file |> 
  column_to_rownames("Gene") |> 
  select(rownames(metadata))

# Revisa si hay NA en la tabla
sum(colSums(is.na(rawcounts)))
colSums(is.na(rawcounts))

# Revisar por qué hay NAs, no debería.
# Si no son muchos NAs, sustituyo esos valores por cero. Si son muchos, revisar a fondo.
# rawcounts[is.na(rawcounts[, "RTP202_10"]), "RTP202_10"] <- 0
# colSums(is.na(rawcounts))

# Compruebo que el orden de rawcounts y de metadata es el mismo
# SI NO ESTAN EN EL MISMO ORDEN Y SE LLAMAN IGUAL, DESEQ2 NO FUNCIONA Y DA ERROR
all(names(rawcounts) == metadata$PDX_pase2)
### Aquí debería haber un warning customizado que avise de lo que ocurre y permita al usuario volver a empezar
### aunque haya ordenado las filas desde el principio

# Construyo el DESeqDataSet (es un Dataset con toda la info que vamos a usar luego)
dds <- DESeqDataSetFromMatrix(countData = rawcounts,
                              colData = metadata,
                              design = ~ Type)
# Comprobacion de que esta bien
dds # Resumen del objeto dds
head(assay(dds)) # Datos de los counts
colData(dds) # Metadatos
design(dds) # Variable cualitativa del diseño
validObject(dds) # Comprueba que cumple los requisitos para seguir
summary(assay(dds)) # Resumen de cada muestra

vsd <- vst(dds, blind=FALSE) # Normalizacion
PCA_type <- plotPCA(vsd, intgroup="Type") # PCA de Type
# Guardo el PCA
path_PCA_type = file.path("out", "img", "21_SCLC_ADC_PCA_type.pdf")
pdf(file = path_PCA_type)
PCA_type
dev.off()

PCA_subtype <- plotPCA(vsd, intgroup="Subtype") # PCA de Subtype
# Guardo el PCA
path_PCA_subtype = file.path("out", "img", "21_SCLC_ADC_PCA_subtype.pdf")
pdf(file = path_PCA_subtype)
PCA_subtype
dev.off()

levels(dds$Type) # Comprueba que los niveles de la variable cualitativa son correctos
levels(dds$Subtype)

# Filtrado de conteos bajos: pre-filtering to keep only rows that have a count 
# of at least 10 for a minimal number of samples
# 10 por defecto, como valor minimo de conteos, razonable para bulk
# minimo numero de muestras en los que tiene que haber el minimo conteo es 1 
# (con que haya una de las muestras con conteos, vale)
smallestGroupSize <- 1
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
dds <- dds[keep,]

# Aqui ya hago el análisis DESeq para obtener DEG.
# Lo guarda en la misma variable, en el apartado results()
dds <- DESeq(dds)
# Extraigo los resultados, y me aseguro que el contraste es el que me interesa. Limpio la tabla
# log2FC = log2(SCLC/NSCLC)
SCLC_ADC <- results(dds, contrast = c("Type", "SCLC", "ADC"), tidy = T) |>
  as_tibble() |> dplyr::rename(gene=row) |> na.omit()
# SCLC_ADC <- left_join(SCLC_ADC, surface_genes_info, by = join_by(gene == `UniProt gene`))
SCLC_ADC

# Guardo la tabla para entregar en informe
DEG_all_file_path = file.path("out", "21_SCLC_ADC_DEG_allgenes.csv")
write_csv(SCLC_ADC, file = DEG_all_file_path)

# Creo un volcano plot con todos
# DESeq2 me da los log2FC. Esto me indica que por cada unidad de log2FC que aumenta, la expresión aumenta en el doble
ev <- EnhancedVolcano(
  SCLC_ADC,
  lab = SCLC_ADC$gene,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'SCLC vs ADC',
  pCutoff = 10e-2,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 6.0,
  subtitle = "",
  legendPosition = "none",
  selectLab = c("SEZ6", "DLL3", "H7B3", "HER2"),
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = "green"
)
ev
# Guardo el volcano plot
path_volcano_all = file.path("out", "img", "21_SCLC_ADC_volcano_plot_allgenes.pdf")
pdf(file=path_volcano_all)
ev
dev.off()

# Guardo tabla con DEG a padj < 0.05 y FC > 2
SCLC_ADC_pval05_FC2 <- SCLC_ADC |> 
  filter(padj < 0.05, abs(log2FoldChange) > 1)
SCLC_ADC_pval05_FC2

DEG_pval05_FC2_file_path = file.path("out", "21_SCLC_ADC_DEG_allgenes_pval05_FC2.csv")
write_csv(SCLC_ADC_pval05_FC2, file = DEG_pval05_FC2_file_path)

# Selecciono los genes con ese padj y FC
genes_selected_SCLC_ADC_all <- SCLC_ADC_pval05_FC2 |> pull(gene)

### CREO HEATMAP

# normalizo los raw counts
vst_counts <- dds |> vst(blind = FALSE)

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-1.5, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 1.5, length.out=floor(palette_length/2)))

vst_counts |> colnames() |> sort()


# selecciono los counts normalizados que correspondan a los DEG (genes_selected) y hago heatmap
hm_out <- assay(vst_counts)[genes_selected_SCLC_ADC_all, ] |> pheatmap(scale = "row",
                                                     cluster_cols = T,
                                                     border_color = NA,
                                                     color = my_color,
                                                     breaks = my_breaks,
                                                     annotation_col = metadata,
                                                     fontsize_row = 12,
                                                     fontsize_col = 12,
                                                     gaps_col = c(9),
                                                     show_rownames = F,
                                                     treeheight_row = 0
                                                     )

# Guardo el heatmap
path_heatmap = file.path("out", "img", "21_SCLC_ADC_heatmap_allgenes.pdf")
pdf(file=path_heatmap)
hm_out
dev.off()

### FILTRAR POR GENES DE SUPERFICIE
### Comparación SCLC vs ADC surface genes
SCLC_ADC_surface <- SCLC_ADC |> filter(gene %in% surface_genes_info$`UniProt gene`)

nrow(SCLC_ADC_surface)

# Guardo la tabla para entregar en informe
SCLC_ADC_surface_all_file_path = file.path("out", "21_SCLC_ADC_DEG_surfacegenes.csv")
write_csv(SCLC_ADC_surface, file = SCLC_ADC_surface_all_file_path)

# Creo un volcano plot con los surface
ev_surface <- EnhancedVolcano(
  SCLC_ADC_surface,
  lab = SCLC_ADC_surface$gene,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'SCLC vs ADC surface genes',
  pCutoff = 10e-2,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 6.0,
  subtitle = "",
  legendPosition = "none",
  selectLab = c("SEZ6", "DLL3", "H7B3", "HER2"),
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = "green"
)
ev_surface

# Guardo el volcano plot
path_volcano_surface_all = file.path("out", "img", "21_SCLC_ADC_volcano_plot_surfacegenes.pdf")
pdf(file=path_volcano_surface_all)
ev_surface
dev.off()

# Guardo tabla con DEG a padj < 0.05 y FC > 2
SCLC_ADC_surface_pval05_FC2 <- SCLC_ADC_surface |> 
  filter(padj < 0.05, abs(log2FoldChange) > 1)
SCLC_ADC_surface_pval05_FC2

SCLC_ADC_surface_pval05_FC2_file_path = file.path("out", "21_SCLC_ADC_DEG_surfacegenes_pval05_FC2.csv")
write_csv(SCLC_ADC_surface_pval05_FC2, file = SCLC_ADC_surface_pval05_FC2_file_path)

# Selecciono los genes de superficie con ese padj y FC
genes_SCLC_ADC_surface_selected <- SCLC_ADC_surface_pval05_FC2 |> pull(gene)

### CREO HEATMAP

# normalizo los raw counts
vst_counts <- dds |> vst(blind = FALSE)

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-1.5, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 1.5, length.out=floor(palette_length/2)))

vst_counts |> colnames() |> sort()


# selecciono los counts normalizados que correspondan a los DEG (genes_selected) y hago heatmap
hm_surface_out <- assay(vst_counts)[genes_SCLC_ADC_surface_selected, ] |>
  pheatmap(
    scale = "row",
    cluster_cols = T,
    border_color = NA,
    color = my_color,
    breaks = my_breaks,
    annotation_col = metadata,
    fontsize_row = 12,
    fontsize_col = 12,
    show_rownames = F,
    treeheight_row = 0
  )

# Guardo el heatmap
path_heatmap_surface = file.path("out", "img", "21_SCLC_ADC_heatmap_surfacegenes.pdf")
pdf(file=path_heatmap_surface)
hm_surface_out
dev.off()


# SCLC_vs_SCC ----------------------------------------------------------------
### Comparación SCLC vs SCC all genes
# Selecciono los metadatos de mis muestras
metadata <- metadata_file |>
  filter(ID != "TP191") |> # quito el TP191 porque no está bien clasificado
  filter(Type == "SCLC" | Subtype == "SCC") |> # selecciono los SCLC y los ADC
  select(ID, Type, Subtype) |> 
  arrange(Type, Subtype) |> 
  column_to_rownames("ID")
metadata$Type <- factor(metadata$Type) # Convierto Type en un factor
metadata$Subtype <- factor(metadata$Subtype) # Convierto Subtype en un factor

# Indico las muestras a estudiar
# muestras <- c("CDX1",
#               "CDX5",
#               "CDX9",
#               "TP135",
#               "TP171",
#               "TP176",
#               "TP185",
#               "TP200",
#               "RTP202_10",
#               "TP103",
#               "TP118",
#               "TP126",
#               "TP134",
#               "TP143",
#               "TP148",
#               "TP166",
#               "TP181",
#               "TP184")

# Subset de las rawcounts de las muestras de interes
rawcounts <- rawcounts_file |> 
  column_to_rownames("Gene") |> 
  select(rownames(metadata))

# Revisa si hay NA en la tabla
sum(colSums(is.na(rawcounts)))
colSums(is.na(rawcounts))

# Revisar por qué hay NAs, no debería.
# Si no son muchos NAs, sustituyo esos valores por cero. Si son muchos, revisar a fondo.
# rawcounts[is.na(rawcounts[, "RTP202_10"]), "RTP202_10"] <- 0
# colSums(is.na(rawcounts))

# Compruebo que el orden de rawcounts y de metadata es el mismo
# SI NO ESTAN EN EL MISMO ORDEN Y SE LLAMAN IGUAL, DESEQ2 NO FUNCIONA Y DA ERROR
all(names(rawcounts) == metadata$PDX_pase2)
### Aquí debería haber un warning customizado que avise de lo que ocurre y permita al usuario volver a empezar
### aunque haya ordenado las filas desde el principio

# Construyo el DESeqDataSet (es un Dataset con toda la info que vamos a usar luego)
dds <- DESeqDataSetFromMatrix(countData = rawcounts,
                              colData = metadata,
                              design = ~ Type)
# Comprobacion de que esta bien
dds # Resumen del objeto dds
head(assay(dds)) # Datos de los counts
colData(dds) # Metadatos
design(dds) # Variable cualitativa del diseño
validObject(dds) # Comprueba que cumple los requisitos para seguir
summary(assay(dds)) # Resumen de cada muestra

vsd <- vst(dds, blind=FALSE) # Normalizacion
PCA_type <- plotPCA(vsd, intgroup="Type") # PCA de Type
# Guardo el PCA
path_PCA_type = file.path("out", "img", "21_SCLC_SCC_PCA_type.pdf")
pdf(file = path_PCA_type)
PCA_type
dev.off()

PCA_subtype <- plotPCA(vsd, intgroup="Subtype") # PCA de Subtype
# Guardo el PCA
path_PCA_subtype = file.path("out", "img", "21_SCLC_SCC_PCA_subtype.pdf")
pdf(file = path_PCA_subtype)
PCA_subtype
dev.off()

levels(dds$Type) # Comprueba que los niveles de la variable cualitativa son correctos
levels(dds$Subtype)

# Filtrado de conteos bajos: pre-filtering to keep only rows that have a count 
# of at least 10 for a minimal number of samples
# 10 por defecto, como valor minimo de conteos, razonable para bulk
# minimo numero de muestras en los que tiene que haber el minimo conteo es 1 
# (con que haya una de las muestras con conteos, vale)
smallestGroupSize <- 1
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
dds <- dds[keep,]

# Aqui ya hago el análisis DESeq para obtener DEG.
# Lo guarda en la misma variable, en el apartado results()
dds <- DESeq(dds)
# Extraigo los resultados, y me aseguro que el contraste es el que me interesa. Limpio la tabla
# log2FC = log2(SCLC/NSCLC)
SCLC_SCC <- results(dds, contrast = c("Type", "SCLC", "SCC"), tidy = T) |>
  as_tibble() |> dplyr::rename(gene=row) |> na.omit()
# SCLC_SCC <- left_join(SCLC_SCC, surface_genes_info, by = join_by(gene == `UniProt gene`))
SCLC_SCC

# Guardo la tabla para entregar en informe
DEG_all_file_path = file.path("out", "21_SCLC_SCC_DEG_allgenes.csv")
write_csv(SCLC_SCC, file = DEG_all_file_path)

# Creo un volcano plot con todos
# DESeq2 me da los log2FC. Esto me indica que por cada unidad de log2FC que aumenta, la expresión aumenta en el doble
ev <- EnhancedVolcano(
  SCLC_SCC,
  lab = SCLC_SCC$gene,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'SCLC vs SCC',
  pCutoff = 10e-2,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 6.0,
  subtitle = "",
  legendPosition = "none",
  selectLab = c("SEZ6", "DLL3", "H7B3", "HER2"),
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = "green"
)
ev
# Guardo el volcano plot
path_volcano_all = file.path("out", "img", "21_SCLC_SCC_volcano_plot_allgenes.pdf")
pdf(file=path_volcano_all)
ev
dev.off()

# Guardo tabla con DEG a padj < 0.05 y FC > 2
SCLC_SCC_pval05_FC2 <- SCLC_SCC |> 
  filter(padj < 0.05, abs(log2FoldChange) > 1)
SCLC_SCC_pval05_FC2

DEG_pval05_FC2_file_path = file.path("out", "21_SCLC_SCC_DEG_allgenes_pval05_FC2.csv")
write_csv(SCLC_SCC_pval05_FC2, file = DEG_pval05_FC2_file_path)

# Selecciono los genes con ese padj y FC
genes_selected_SCLC_SCC_all <- SCLC_SCC_pval05_FC2 |> pull(gene)

### CREO HEATMAP

# normalizo los raw counts
vst_counts <- dds |> vst(blind = FALSE)

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-1.5, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 1.5, length.out=floor(palette_length/2)))

vst_counts |> colnames() |> sort()


# selecciono los counts normalizados que correspondan a los DEG (genes_selected) y hago heatmap
hm_out <- assay(vst_counts)[genes_selected_SCLC_SCC_all, ] |> pheatmap(scale = "row",
                                                                       cluster_cols = T,
                                                                       border_color = NA,
                                                                       color = my_color,
                                                                       breaks = my_breaks,
                                                                       annotation_col = metadata,
                                                                       fontsize_row = 12,
                                                                       fontsize_col = 12,
                                                                       gaps_col = c(9),
                                                                       show_rownames = F,
                                                                       treeheight_row = 0
)

# Guardo el heatmap
path_heatmap = file.path("out", "img", "21_SCLC_SCC_heatmap_allgenes.pdf")
pdf(file=path_heatmap)
hm_out
dev.off()

### FILTRAR POR GENES DE SUPERFICIE
### Comparación SCLC vs ADC surface genes
SCLC_SCC_surface <- SCLC_SCC |> filter(gene %in% surface_genes_info$`UniProt gene`)

nrow(SCLC_SCC_surface)

# Guardo la tabla para entregar en informe
SCLC_SCC_surface_all_file_path = file.path("out", "21_SCLC_SCC_DEG_surfacegenes.csv")
write_csv(SCLC_SCC_surface, file = SCLC_SCC_surface_all_file_path)

# Creo un volcano plot con los surface
ev_surface <- EnhancedVolcano(
  SCLC_SCC_surface,
  lab = SCLC_SCC_surface$gene,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'SCLC vs SCC surface genes',
  pCutoff = 10e-2,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 6.0,
  subtitle = "",
  legendPosition = "none",
  selectLab = c("SEZ6", "DLL3", "H7B3", "HER2"),
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = "green"
)
ev_surface

# Guardo el volcano plot
path_volcano_surface_all = file.path("out", "img", "21_SCLC_SCC_volcano_plot_surfacegenes.pdf")
pdf(file=path_volcano_surface_all)
ev_surface
dev.off()

# Guardo tabla con DEG a padj < 0.05 y FC > 2
SCLC_SCC_surface_pval05_FC2 <- SCLC_SCC_surface |> 
  filter(padj < 0.05, abs(log2FoldChange) > 1)
SCLC_SCC_surface_pval05_FC2

SCLC_SCC_surface_pval05_FC2_file_path = file.path("out", "21_SCLC_SCC_DEG_surfacegenes_pval05_FC2.csv")
write_csv(SCLC_SCC_surface_pval05_FC2, file = SCLC_SCC_surface_pval05_FC2_file_path)

# Selecciono los genes de superficie con ese padj y FC
genes_SCLC_SCC_surface_selected <- SCLC_SCC_surface_pval05_FC2 |> pull(gene)

### CREO HEATMAP

# normalizo los raw counts
vst_counts <- dds |> vst(blind = FALSE)

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-1.5, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 1.5, length.out=floor(palette_length/2)))

vst_counts |> colnames() |> sort()


# selecciono los counts normalizados que correspondan a los DEG (genes_selected) y hago heatmap
hm_surface_out <- assay(vst_counts)[genes_SCLC_SCC_surface_selected, ] |>
  pheatmap(
    scale = "row",
    cluster_cols = T,
    border_color = NA,
    color = my_color,
    breaks = my_breaks,
    annotation_col = metadata,
    fontsize_row = 12,
    fontsize_col = 12,
    show_rownames = F,
    treeheight_row = 0
  )

# Guardo el heatmap
path_heatmap_surface = file.path("out", "img", "21_SCLC_SCC_heatmap_surfacegenes.pdf")
pdf(file=path_heatmap_surface)
hm_surface_out
dev.off()


## Intersección entre las dos comparaciones -----------------------------------
genes_SCLC_both <- intersect(SCLC_ADC_pval05_FC2$gene, SCLC_SCC_pval05_FC2$gene)
length(genes_SCLC_both)
genes_SCLC_surface_both <- intersect(SCLC_ADC_surface_pval05_FC2$gene, SCLC_SCC_surface_pval05_FC2$gene)
length(genes_SCLC_surface_both)

genes_SCLC_both_path = file.path("out", "21_genes_SCLC_both.txt")
write.table(genes_SCLC_both, file = genes_SCLC_both_path, row.names = FALSE, col.names = FALSE)

genes_SCLC_surface_both_path = file.path("out", "21_genes_SCLC_surface_both.txt")
write.table(genes_SCLC_surface_both, file = genes_SCLC_surface_both_path, row.names = FALSE, col.names = FALSE)

# Saco la lista de genes up y down. También los up y down en ambos
# (hay genes que están up en uno y down en otro. Esos no estan)
up_genes_surface_SCLC_ADC <- 
  SCLC_ADC_surface_pval05_FC2 |> filter(log2FoldChange > 0)

down_genes_surface_SCLC_ADC <- 
  SCLC_ADC_surface_pval05_FC2 |> filter(log2FoldChange < 0)

up_genes_surface_SCLC_SCC <- 
  SCLC_SCC_surface_pval05_FC2 |> filter(log2FoldChange > 0)

down_genes_surface_SCLC_SCC <- 
  SCLC_SCC_surface_pval05_FC2 |> filter(log2FoldChange < 0)

up_genes_surface_SCLC_both <- 
  intersect(up_genes_surface_SCLC_ADC$gene,
            up_genes_surface_SCLC_SCC$gene)
length(up_genes_surface_SCLC_both)

down_genes_surface_SCLC_both <- 
  intersect(down_genes_surface_SCLC_ADC$gene,
            down_genes_surface_SCLC_SCC$gene)
length(down_genes_surface_SCLC_both)

venn.diagram(x = list(up_genes_surface_SCLC_ADC$gene, up_genes_surface_SCLC_SCC$gene), 
             category.names = c("SCLCvsADC", "SCLCvsSCC"), 
             filename = "out/img/21_venn_ADC_SCC_up.tiff", 
             disable.logging = TRUE,
             main = "Overexpressed genes",
             lwd = 2,                      # Grosor del borde
             lty = 'solid',                # Tipo de línea
             fill = c("lightblue", "pink"),# Colores de fondo
             cex = 2,                      # Tamaño de número en el diagrama
             fontface = "bold",            # Negrita
             fontfamily = "sans",
             cat.cex = 2,                  # Tamaño del texto de categorías
             cat.fontface = "bold",
             cat.default.pos = "outer",   # Posición del texto
             cat.pos = c(-20, 20),         # Ángulos del texto de categorías
             cat.dist = 0.05               # Distancia del texto respecto al círculo)
)

venn.diagram(x = list(down_genes_surface_SCLC_ADC$gene, down_genes_surface_SCLC_SCC$gene), 
             category.names = c("SCLCvsADC", "SCLCvsSCC"), 
             filename = "out/img/21_venn_ADC_SCC_down.tiff", 
             disable.logging = TRUE,
             main = "Underexpresed genes",
             lwd = 2,                      # Grosor del borde
             lty = 'solid',                # Tipo de línea
             fill = c("lightblue", "pink"),# Colores de fondo
             cex = 2,                      # Tamaño de número en el diagrama
             fontface = "bold",            # Negrita
             fontfamily = "sans",
             cat.cex = 2,                  # Tamaño del texto de categorías
             cat.fontface = "bold",
             cat.default.pos = "outer",   # Posición del texto
             cat.pos = c(-20, 20),         # Ángulos del texto de categorías
             cat.dist = 0.05               # Distancia del texto respecto al círculo)
)

# Creo un tib con los nombres de los genes de superficie comunes,
# y un score aleatorio y personal para puntuar cada gen en función de cada análisis

# SCLC_surface_both <- tibble(gene = genes_SCLC_surface_both) |> 
#   left_join(SCLC_ADC_surface_pval05_FC2 |> select(gene, log2FoldChange)) |> 
#   dplyr::rename(log2FoldChange_ADC = log2FoldChange) |> 
#   mutate(score_ADC = (log2FoldChange_ADC-1)*100/max(log2FoldChange_ADC)) |> 
#   left_join(SCLC_SCC_surface_pval05_FC2 |> select(gene, log2FoldChange)) |> 
#   dplyr::rename(log2FoldChange_SCC = log2FoldChange) |>
#   mutate(score_SCC = (log2FoldChange_SCC-1)*100/max(log2FoldChange_SCC))
#   

# Protein from paper ------------------------------------------------------

protein_dif_expr <- read_excel("res/1-s2.0-S0092867423013351-mmc3.xlsx", sheet = "Table S3A", col_names = TRUE)

paper_protein_DEG_all_path = file.path("out", "21_paper_DEG_allproteins.csv")
write_csv(protein_dif_expr, file = paper_protein_DEG_all_path)

protein_dif_expr_pval05_FC2 <- protein_dif_expr |> filter(FDR < 0.05, abs(Log2.median.foldchange) > 1)

paper_protein_DEG_pval05_FC2_all_path = file.path("out", "21_paper_DEG_allproteins_pval05_FC2.csv")
write_csv(protein_dif_expr_pval05_FC2, file = paper_protein_DEG_pval05_FC2_all_path)

ev_protein <- EnhancedVolcano(
  protein_dif_expr,
  lab = protein_dif_expr$Gene.symbol,
  x = "Log2.median.foldchange",
  y = "FDR",
  pCutoff = 10e-2,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 6.0,
  subtitle = "",
  legendPosition = "none",
  selectLab = c("SEZ6", "DLL3", "H7B3", "HER2"),
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = "green"
)
ev_protein

path_volcano_paper_protein = file.path("out", "img", "21_paper_volcano_plot_allproteins.pdf")
pdf(file=path_volcano_paper_protein)
ev_protein
dev.off()

protein_surface <- protein_dif_expr |> filter(Gene.symbol %in% surface_genes_info$`UniProt gene`)

paper_protein_surface_path = file.path("out", "21_paper_DEG_surfaceproteins.csv")
write_csv(protein_surface, file = paper_protein_surface_path)

ev_protein_surface <- EnhancedVolcano(
  protein_surface,
  lab = protein_surface$Gene.symbol,
  x = "Log2.median.foldchange",
  y = "FDR",
  pCutoff = 10e-2,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 6.0,
  subtitle = "",
  legendPosition = "none",
  selectLab = c("SEZ6", "DLL3", "H7B3", "HER2"),
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = "green"
)

ev_protein_surface

path_volcano_paper_protein_surface = file.path("out", "img", "21_paper_volcano_plot_surfaceproteins.pdf")
pdf(file=path_volcano_paper_protein_surface)
ev_protein_surface
dev.off()

protein_surface_pval05_FC2 <- protein_surface |> 
  filter(FDR < 0.05, abs(Log2.median.foldchange) > 1)
protein_surface_pval05_FC2

paper_protein_DEG_pval05_FC2_surface_path = file.path("out", "21_paper_DEG_surfaceproteins_pval05_FC2.csv")
write_csv(protein_surface_pval05_FC2, file = paper_protein_DEG_pval05_FC2_surface_path)

genes_protein_surface_selected <- protein_surface_pval05_FC2 |> pull(Gene.symbol)

up_protein_surface <- 
  protein_surface_pval05_FC2 |> filter(Log2.median.foldchange > 0)

down_protein_surface <- 
  protein_surface_pval05_FC2 |> filter(Log2.median.foldchange < 0)

# Comparación con datos inhouse
up_surface_inhouse_paper <- 
  intersect(up_genes_surface_SCLC_both, up_protein_surface$Gene.symbol)
length(up_surface_inhouse_paper)
up_surface_inhouse_paper

up_surface_inhouse_paper_path = file.path("out", "21_upgenesproteins_surface_inhouse_paper.txt")
write.table(up_surface_inhouse_paper, file = up_surface_inhouse_paper_path, row.names = FALSE, col.names = FALSE)

venn.diagram(x = list(up_genes_surface_SCLC_ADC$gene, up_genes_surface_SCLC_SCC$gene, up_protein_surface$Gene.symbol), 
             category.names = c("SCLCvsADC", "SCLCvsSCC", "SCLCvsHealthy"), 
             filename = "out/img/21_venn_healthy_ADC_SCC_up.tiff", 
             disable.logging = TRUE,
             main = "Overexpressed genes",
             lwd = 2,                      # Grosor del borde
             lty = 'solid',                # Tipo de línea
             fill = c("lightblue", "pink", "#823"),# Colores de fondo
             cex = 2,                      # Tamaño de número en el diagrama
             fontface = "bold",            # Negrita
             fontfamily = "sans",
             cat.cex = 2,                  # Tamaño del texto de categorías
             cat.fontface = "bold",
             cat.default.pos = "outer",   # Posición del texto
             cat.pos = c(-20, 20, 0),         # Ángulos del texto de categorías
             cat.dist = 0.05               # Distancia del texto respecto al círculo)
)

down_surface_inhouse_paper <- 
  intersect(down_genes_surface_SCLC_both, down_protein_surface$Gene.symbol)
length(down_surface_inhouse_paper)
down_surface_inhouse_paper

down_surface_inhouse_paper_path = file.path("out", "21_downgenesproteins_surface_inhouse_paper.txt")
write.table(down_surface_inhouse_paper, file = down_surface_inhouse_paper_path, row.names = FALSE, col.names = FALSE)

venn.diagram(x = list(down_genes_surface_SCLC_ADC$gene, down_genes_surface_SCLC_SCC$gene, down_protein_surface$Gene.symbol), 
             category.names = c("SCLCvsADC", "SCLCvsSCC", "SCLCvsHealthy"), 
             filename = "out/img/21_venn_healthy_ADC_SCC_down.tiff", 
             disable.logging = TRUE,
             main = "Underexpressed genes",
             lwd = 2,                      # Grosor del borde
             lty = 'solid',                # Tipo de línea
             fill = c("lightblue", "pink", "#823"),# Colores de fondo
             cex = 2,                      # Tamaño de número en el diagrama
             fontface = "bold",            # Negrita
             fontfamily = "sans",
             cat.cex = 2,                  # Tamaño del texto de categorías
             cat.fontface = "bold",
             cat.default.pos = "outer",   # Posición del texto
             cat.pos = c(-20, 20, -200),         # Ángulos del texto de categorías
             cat.dist = 0.05               # Distancia del texto respecto al círculo)
)
