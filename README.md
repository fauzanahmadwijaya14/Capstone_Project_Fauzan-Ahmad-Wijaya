## Capstone_Project_Fauzan-Ahmad-Wijaya

**Pendahuluan**

Percobaan ini dilakukan untuk menganalisis perbedaan ekspresi gen pada data GSE 161533 milik Qiu et al. (2020). GSE ini berisi data perbandingan Esophageal Squamous Cell Carcinoma (ESCC) dengan kondisi normal serta paratumor. ESCC merupakan salah satu jenis kanker paling agresif di dunia dan sebagian besar kasus berasal dari Asia Timur, Iran Utara, Asia Tengah, dan Cina (Chen et al., 2025; Weidenbaum & Gibson, 2022). ESCC merupakan jenis paling dominan dari kanker esofagus (>90%) dengan tingkat deteksi awal rendah, keluasan genomik, resistensi pengobatan, penyebaran lokal yang cepat, dan tingkat harapan hidup 5 tahun hanya sekitar 15-25% (Chen et al., 2025; Wang et al., 2025). ESCC sendiri terjadi di bagian atas hingga tengah dari esofagus, sedangkan Esophageal Adenocarcinoma (EAC) terjadi di bagian distal serta persimpangan gastroesofagus (Weidenbaum & Gibson, 2022). Pada proyek ini dilakukan analisis DEG dengan peranti R untuk mengetahui gen-gen yang terespresikan serta relasinya dengan jalur metabolisme khususnya pada kasus ESCC.

**Metode**

Data set yang dipilih merupakan [GSE161533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161533). Pada analisis ini, data diakses langsung menggunakan R versi 4.5.2. Pemilihan R sebagai peranti dalam analisis karena R memfasilitasi analisis serta visualisasi data yang lebih luas dan fleksibel. Sebelumnya, R Studio yang digunakan sudah terinstal BiocManager dengan sejumlah tools tambahan untuk mendukung analisis seperti GEOquery, limma, dplyr, ggplot2, pheatmap, gplots, dan RcolorBrewer.


```sh
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

```

