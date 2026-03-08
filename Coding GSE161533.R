##PENGOLAHAN DATA

#A: Persiapan Lingkungan Kerja

#Pemanggilan library
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(hgu133a.db)



#B: Pengambilan Data Dari GEO

gset <- getGEO("GSE161533", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)



#C: Pre-Processing Data Ekspresi

#Transformasi data agar data bersifat vektor dan Log2. Hal ini untuk:
#1. Menstabilkan varians
#2. Mendekati asumsi model linear
#3. Memudahkan interpretasi log fold change

ex <- exprs(gset)

qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}



#D: Definisi Kelompok Sampel

group_info <- pData(gset)[["characteristics_ch1"]]

groups <- make.names(group_info)

gset$group <- factor(groups)

nama_grup <- levels(gset$group)
print(nama_grup)



#E: Definisi Matrix (Kerangka Statistik)

design <- model.matrix(~0 + gset$group)

colnames(design) <- levels(gset$group)

grup_normal <- nama_grup[1]
grup_paratumor <- nama_grup [2]
grup_tumor <- nama_grup[3]

contrast_formula1 <- paste(grup_paratumor, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula1))

contrast_formula2 <- paste(grup_tumor, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula2))

#lmFit():
#Membangun model linear untuk setiap gen
fit <- lmFit(ex, design)

#makeContrasts(): mendefinisikan perbandingan grup normal dengan paratumor
contrast_matrix1 <- makeContrasts(contrasts = contrast_formula1, levels = design)

#contrasts.fit(): menerapkan kontras ke model grup normal dengan paratumor
fit2a <- contrasts.fit(fit, contrast_matrix1)

#makeContrasts(): mendefinisikan perbandingan grup normal dengan tumor
contrast_matrix2 <- makeContrasts(contrasts = contrast_formula2, levels = design)

#contrasts.fit(): menerapkan kontras ke model grup normal dengan tumor
fit2b <- contrasts.fit(fit, contrast_matrix2)

#eBayes():
#Empirical Bayes untuk menstabilkan estimasi varians
fit2a <- eBayes(fit2a)
fit2b <- eBayes(fit2b)

#topTable():
#Mengambil hasil akhir DEG
#adjust = "fdr" -> koreksi multiple testing
#p.value = 0.01  -> gen sangat signifikan

#grup normal dengan paratumor
topTableResults_A <- topTable(
  fit2a,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.1
)

#grup normal dengan tumor
topTableResults_B <- topTable(
  fit2b,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.1
)

head(topTableResults_A)
head(topTableResults_B)



#G:  Anotasi Nama Gen

#Pengambilan ID probe dari hasil DEG
prob_ids_A <- rownames(topTableResults_A)
prob_ids_B <- rownames(topTableResults_B)

#Mapping probe antara gene symbol dengan gene name
#Grup normal dengan paratumor
gene_annotation_A <- AnnotationDbi::select(
  hgu133a.db,
  keys = prob_ids_A,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

#Grup normal dengan tumor
gene_annotation_B <- AnnotationDbi::select(
  hgu133a.db,
  keys = prob_ids_B,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

#Penggabungan dengan hasil limma
topTableResults_A$PROBEID <- rownames(topTableResults_A)
topTableResults_B$PROBEID <- rownames(topTableResults_B)

topTableResults_A <- merge(
  topTableResults_A,
  gene_annotation_A,
  by = "PROBEID",
  all.x = TRUE
)

topTableResults_B <- merge(
  topTableResults_B,
  gene_annotation_B,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil anotasi
head(topTableResults_A[, c("PROBEID", "SYMBOL", "GENENAME")])
head(topTableResults_B[, c("PROBEID", "SYMBOL", "GENENAME")])



##VISUALISASI DATA

#H.1 Boxplot Distribusi Nilai Ekspresi

group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)


#H.2 Distribusi Nilai Ekspresi (Density Plot)

expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x= Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() + 
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )


#H.3 UMAP

umap_input <- t(ex)

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )


#I.1 Volcano Plot

#Karena gen hasil DEG antara Normal vs Paratumor hanya satu,
#maka hanya dibuat untuk Normal vs Tumor
#Definisi volcano plot 
volcano_data_B <- data.frame(
  logFC = topTableResults_B$logFC,
  adj.P.Val = topTableResults_B$adj.P.Val,
  Gene = topTableResults_B$SYMBOL
)

#Klasifikasi status gen
volcano_data_B$status <- "NO"
volcano_data_B$status[volcano_data_B$logFC > 1 & volcano_data_B$adj.P.Val < 0.01] <- "UP"
volcano_data_B$status[volcano_data_B$logFC < -1 & volcano_data_B$adj.P.Val < 0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data_B, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Normal dengan Tumor")


#I.2 Visualisasi Heatmap

#Karena gen hasilDEG untuk Normal vs Paratumor hanya satu,
#maka hanya akan dibuat heatmap untuk Normal vs Tumor

#Pemilihan 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults_B <- topTableResults_B[
  order(topTableResults_B$adj.P.Val),
]

top50_B <- head(topTableResults_B, 50)

#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50_B$PROBEID, ]

#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label_B <- ifelse(
  is.na(top50_B$SYMBOL) | top50_B$SYMBOL == "",
  top50_B$PROBEID,      # jika SYMBOL kosong → probe ID
  top50_B$SYMBOL        # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label_B

#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi heatmap 
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)



##PENUTUP

#J: Menyimpan Hasil
write.csv(topTableResults_A, "Hasil_GSE161533_DEG_NormalvsParatumor.csv")
write.csv(topTableResults_B, "Hasil_GSE161533_DEG_NormalvsTumor_seluruh.csv")
write.csv(top50_B, "Hasil_GSE161533_DEG_NormalvsTumor_top50.csv")

message("Analisis selesai. File hasil telah disimpan.")


#Hasil GO Top50_B http://biit.cs.ut.ee/gplink/l/avkP_vq2jS9