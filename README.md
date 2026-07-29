# rnaseq-hfd-vs-lfd-urothelium
This repository contains my mini project completed as part of the RNA-seq Bioinformatics Bootcamp organized by INBIO Indonesia. The project reproduces and extends the transcriptomic analysis of mouse urothelium exposed to a high-fat diet using publicly available GEO datasets.

## Comparative RNA-seq Analysis of Mouse Urothelium Following High-Fat Diet

Tentang: Pengaruh Obesitas terhadap Infeksi Saluran Kemih (ISK/UTI) melalui treatment 
High Fat Diet (HFD) vs Low Fat Diet (LFD)
menggunakan dataset publik GEO Kode GSE294660

### 01 DATA PREPARATION 

#### Instalasi Package
- **Install BiocManager (Pengelola Packages Bioconductor)**
  ```r
  install.packages("BiocManager")
  ```
- **Install semua packages untuk RNA-seq**
  ```r
  BiocManager::install(c("DESeq2", "ggplot2", "pheatmap", "GEOquery", "dplyr", "tidyr"))
  ```
- **Load library setiap kali mau pakai**
  ```r
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(GEOquery)
  library(dplyr)
  library(tidyr)
  ```

#### Data Preparation
#####  a. Download Data
- **Load packages**
  ```r
  library(GEOquery)
  ```
- **download series_matrix GSE**
  ```r
  gse <- getGEO("GSE294660", GSEMatrix = TRUE)
  ```
- **download semua file supplementary dari GSE**
  ```r
  getGEOSuppFiles("GSE294660")
  ```

##### b. Data Cleaning
- **Ambil data ekspresi dan metadata dari series matrix**
  ```r
  gse_data<-gse[[1]]
  expression_data<-exprs(gse_data)
  metadata<-pData(gse_data)
  ```

- **Cek struktur datanya**
  ```r
  dim(expression_data)
  ```
- **Melihat supplementary file yang sudah didownload** 
  ```r
  list.files("GSE294660/")
  ```
- **Baca isi file, dan simpan sebagai objek "counts". Sesuaikan nama file dengan output list.files tadi**
  ```r
  library(readr)
  counts <- read.csv(
    "GSE294660_dio_12WOD_urothelium_baseline_counts.csv",
    row.names = 1,
    check.names = FALSE
  )
  dim(counts)
  head(counts)
  ```

#### Membuat metadata
- **sederhanakan nama sampel8**
  ```r
  colnames(counts) <- c(
    "F_LFD_1", "F_LFD_2", "F_LFD_3", "F_LFD_4",
    "F_HFD_1", "F_HFD_2", "F_HFD_3", "F_HFD_4",
    "M_LFD_1", "M_LFD_2", "M_LFD_3", "M_LFD_4",
    "M_HFD_1", "M_HFD_2", "M_HFD_3", "M_HFD_4"
  )
  ```
-  **cek nama sampel**
  ```r
  colnames(counts)
  ```

- **membuat metadata**
  ```r
  metadata <- data.frame(
    sample = colnames(counts),
    sex = factor(c(rep("Female", 8),rep("Male", 8)
    )),
    condition = factor(c(rep("LFD", 4),rep("HFD", 4),rep("LFD", 4),rep("HFD", 4)
    ))
  )
  ```
- **Jadikan nama sampel sebagai row names**
  ```r
  rownames(metadata) <- metadata$sample
  ```
- **Tampilkan metadata**
  ```r
  metadata
  ```
- **Verifikasi dengan**
  ```r
  colnames(counts) == rownames(metadata)
  semua harus TRUE
  ```
##### Distribusi Data
- **Load library untuk visualisasi**
  ```r
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  ```
- **ubah data ke format panjang (long format) untuk ggplot**
  ```r
  counts_long <- counts %>%
    as.data.frame() %>%
    mutate(gene = rownames(.)) %>%
    pivot_longer(
      cols = -gene,
      names_to = "sample",
      values_to = "count"
    ) %>%
    mutate(type = "Data Mentah")
  ```
- **Boxplot (dalam skala log2 agar lebih jelas)**
  ```r
  ggplot(counts_long, aes(x = sample, y = log2(count + 1), fill = sample)) +
    geom_boxplot() +
    theme_minimal() +
    labs(
      title = "Distribusi Ekspresi Gen per Sampel (Data Mentah)",
      x = "Sampel",
      y = "Log2(Count + 1)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ```
##### a. Filter gen dengan reads rendah (opsional)
- **Filter: gen dengan total reads >= 10 di semua sampel**
  ```r
  keep <- rowSums(counts) >= 10
  counts_filtered <- counts[keep, ]
  ```
- **Bandingkan jumlah gen sebelum dan sesudah filter**
  ```r
  dim(counts)
  dim(counts_filtered)
  ```
#### Korelasi antar sampel
- **install dan load pheatmap (jika belum)**
  ```r
  if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")
  library(pheatmap)
  ```
- **Hitung korelasi antar sampel (pakai data yang sudah difilter)**
  ```r
  cor_matrix <- cor(counts_filtered)
  ```
- **Tampilkan sebagai heatmap**
  ```r
  pheatmap(cor_matrix,
           main = "Korelasi Antar Sampel (Data Mentah)",
           display_numbers = TRUE,
           number_format = "%.2f",
           color = colorRampPalette(c("blue", "white", "red"))(50))
  ```
#### PCA (Principal Component Analysis)
- **PCA dengan data yang sudah di-log (karena data counts sangat skewed)**
  ```r
  log_data <- log2(counts_filtered + 1)
  pca_results <- prcomp(t(log_data), scale. = TRUE)
  ```
- **Buat data frame untuk plotting**
  ```r
  pca_df <- data.frame(
    PC1 = pca_results$x[,1],
    PC2 = pca_results$x[,2],
    condition = metadata$condition
  )
  ```
- **Hitung persentase varians**
  ```r
  var_explained <- summary(pca_results)$importance[2, 1:2]*100
  ```
- **Plot PCA**
  ```r
  ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 5) +
    theme_minimal() +
    labs(title = "PCA Plot (Data Mentah - Log2)",
         x = paste0("PC1: ", round(var_explained[1], 1), "% variance"),
         y = paste0("PC2: ", round(var_explained[2], 1), "% variance")) +
    theme(legend.position = "bottom")
  ```

### 02 NORMALIZATION
#### Normalisasi
- **Load library**
  ```r
  library(DESeq2)
  ```

- **Buat DESeq2 object (dari count matrix dan metadata)**
  ```r
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = metadata,
    design = ~ sex + condition
  )
  ```
- **mulai normalisasi**
  ```r
  dds <- DESeq(dds)
  ```
- **ambil size factors (faktor normalisasi)**
  ```r
  sizeFactors(dds)
  ```
- **ambil data yang sudah dinormalisasi. Data ini bisa digunakan untuk visualisasi (heatmap, PCA)**
  ```r
  rld <- rlog(dds, blind = FALSE)
  vsd <- vst(dds, blind = FALSE)
  ```
#### Persiapan visualisasi
- **ubah data normalisasi (rlog) ke format panjang**
  ```r
  rld_data <- assay(rld)
  rld_long <- rld_data %>%
    as.data.frame() %>%
    mutate(gene = rownames(.)) %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
    mutate(type = "Data Normalisasi (rlog)")
  ```
- **gabungkan kedua data**
  ```r
  combine_data <- bind_rows(
    counts_long %>% mutate(value = log2(count + 1)), 
    rld_long %>% mutate(value = expression)
  )
  ```
#### Visualisasi dengan Boxplot
  ```r
  library(ggplot2)
  ggplot(combine_data, aes(x = sample, y = value, fill = type)) +
    geom_boxplot() +
    facet_wrap(~ type, scales = "free_y", ncol = 2) + 
    theme_minimal() +
    labs(title = "Perbandingan Distribusi Ekspresi Gen", 
         subtitle = "Data Mentah vs Data Normalisasi (rlog)",
         x = "Sample",
         y = "Log2 Ekspresi") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "none")
  ```
### 03 Differential Gene Expression (DEG) Analysis
#### DESeq2
- **membuat object dds**
  ```r
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = metadata,
    design = ~ sex + condition
  )
  ```
- **lihat level condition dan jalankan DESEq2**
  ```r
  dds$condition <- relevel(dds$condition, ref = "LFD")
  dds <- DESeq(dds)
  ```
- **ambil hasil perbandingan**
  ```r
  res <- results(dds,
                 contrast = c("condition", "HFD", "LFD"))
  res <- res[order(res$padj), ]
  head(res)
  ```
#### Filter Gen
- **ubah data jadi format data frame**
  ```r
  res_df <-as.data.frame(res)
  ```
- **tambahkan kolom status signifikansi**
  ```r
  res_df$significant <- ifelse(
    !is.na(res_df$padj)&
      res_df$padj <0.05&
      abs(res_df$log2FoldChange) > 1,
    "Signifikan",
    "Tidak Signifikan"
  )
  ```
- **tambahkan kolom arah perubahan**
  ```r
  res_df$direction <- ifelse(
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1,
    "Naik",
    ifelse(
      !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1,
      "Turun",
      "Tidak Signifikan"
    )
  )
  ```
- **filter gen signifikan**
  ```r
  sig_genes <- res_df[res_df$significant == "Signifikan", ]
  nrow(sig_genes)
  ```
- **lihat gen naik dan turun**
  ```r
  table(sig_genes$direction)
  ```
#### Export Hasil
- **export hasil**
  ```r
  write.csv(res_df, "DEG_results_GSE294660_12WOD.csv")
  ```
- **lihat top 10 gen naik, diurutkan berdasarkan Log2FoldChange** 
  ```r
  top_up <- res_df[res_df$direction == "Naik",]
  top_up <- top_up[order(top_up$log2FoldChange,decreasing = TRUE),]
  head(top_up[,c("log2FoldChange","pvalue","padj")],10)
  ```
- **lihat top 10 gen turun**
  ```r
  top_down <- res_df[res_df$direction == "Turun",]
  top_down <- top_down[order(top_down$log2FoldChange),]
  head(top_down[,c("log2FoldChange", "pvalue", "padj")], 10)
  ```
### 04 VISUALIZATION
#### Volcano plot
  ```r
  library(dplyr)
  library(ggplot2)
  ```
- **ambil 10 gen dengan padj terkecil (paling signifikan)**
  ```r
  top_gene <- res_df %>%
    filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >1) %>%
    arrange(padj) %>%
    head(10)
  ```
- **tambahkan kolom tabel (hanya untuk gen top)**
  ```r
  top_genes <- res_df[order(res_df$padj), ][1:10, ]
  res_df$label <- ifelse(
    rownames(res_df) %in% rownames(top_genes),
    rownames(res_df),
    ""
  )
  ```
#### Buat volcano plot dengan label
  ```r
  volcano_plot_labeled <- ggplot(
    res_df,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = direction
    )
  ) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c("Naik" = "red","Turun" = "blue","Tidak Signifikan" = "gray"),
      name = "Regulasi"
    ) +
    geom_vline(xintercept = c(-1, 1),linetype = "dashed",color = "black",alpha = 0.5
    ) +
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",color = "black",alpha = 0.5
    ) +
    geom_text(aes(label = label),vjust = -0.5,size = 3,check_overlap = TRUE
    ) +
    labs(
      title = "Volcano Plot: HFD vs LFD (Top 10 DEG)",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  print(volcano_plot_labeled)
  ggsave(
    "volcano_plot_labeled.png",volcano_plot_labeled,width = 10,height = 7,dpi = 300
  )
  ```
#### Heatmap
  ```r
  library(pheatmap)
  library(DESeq2)
  ```
- **ambil top 30 gen signifikan (berdasarkan padj)**
  ```r
  res_sorted <- res[order(res$padj), ]
  top_genes <- rownames(res_sorted)[!is.na(res_sorted$padj)][1:30]
  ```
- **ambil data ekspresi untuk gen gen tersebut dari rlog**
  ```r
  heatmap_data <- assay(rld)[top_genes, ]
  ```
- **buat annotation untuk sampel (warna berdasarkan kondisi)**
  ```r
  annotation_col <- data.frame(
    Sex = metadata$sex,
    Condition = metadata$condition
  )
  rownames(annotation_col) <- colnames(heatmap_data)
  ```
- **buat heatmap**
  ```r
  heatmap_plot <- pheatmap(
    heatmap_data,
    scale = "row",
    main = "Heatmap Top 30 DEG",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 8,
    annotation_col = annotation_col,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    border_color = NA
    )
    ```

Ini untuk 8 weeks
Tentang: Pengaruh Obesitas terhadap Infeksi Saluran Kemih (ISK/UTI) melalui treatment High Fat Diet (HFD) vs Low Fat Diet (LFD) (8 weeks)
menggunakan dataset publik GEO Kode GSE294660

### 01 DATA PREPARATION 
#### Instalasi Package
- **Install BiocManager (Pengelola Packages Bioconductor)**
  ```r
  install.packages("BiocManager")
  ```
- **Install semua packages untuk RNA-seq**
  ```r
  BiocManager::install(c("DESeq2", "ggplot2", "pheatmap", "GEOquery", "dplyr", "tidyr"))
  ```
- **Load library setiap kali mau pakai**
  ```r
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(GEOquery)
  library(dplyr)
  library(tidyr)
  ```
#### Data Preparation
##### a. Download Data
- **Load packages**
  ```r
  library(GEOquery)
  ```
- **download series_matrix GSE**
  ```r
  gse <- getGEO("GSE294660", GSEMatrix = TRUE)
  ```
- download semua file supplementary dari GSE
  ```r
  getGEOSuppFiles("GSE294660")
  ```
##### b. Data Cleaning
- **Ambil data ekspresi dan metadata dari series matrix**
  ```r
  gse_data<-gse[[1]]
  expression_data<-exprs(gse_data)
  metadata<-pData(gse_data)
  ```
- **Cek struktur datanya**
  ```r
  dim(expression_data)
  ```
- **Melihat supplementary file yang sudah didownload**
  ```r
  list.files("GSE294660/")
  ```
- **Baca isi file, dan simpan sebagai objek "counts". Sesuaikan nama file dengan output list.files tadi**
  ```r
  library(readr)
  counts <- read.csv(
    "GSE294660_dio_8WOD_urothelium_baseline_counts.csv",
    row.names = 1,
    check.names = FALSE
  )
  dim(counts)
  head(counts)
  ```
#### Membuat Metadata, Urutkan kolom menjadi:
##### a. Female LFD -> Female HFD -> Male LFD -> Male HFD
  ```r
  counts <- counts[, c(
      1, 2, 3,      # Female LFD
      4, 8, 9,      # Female HFD
      6, 7, 12,     # Male LFD
      5, 10, 11     # Male HFD
    )]
  ```
##### b. Beri nama baru
  ```r
  colnames(counts) <- c(
    "F_LFD_1", "F_LFD_2", "F_LFD_3",
    "F_HFD_1", "F_HFD_2", "F_HFD_3",
    "M_LFD_1", "M_LFD_2", "M_LFD_3",
    "M_HFD_1", "M_HFD_2", "M_HFD_3"
  )
  ```
-  **cek nama sampel**
  ```r
  colnames(counts)
  ```
- **membuat metadata**
  ```r
  metadata <- data.frame(
    sample = colnames(counts),
    sex = factor(c(rep("Female", 6), rep("Male", 6))),
    condition = factor(c(
      rep("LFD", 3),
      rep("HFD", 3),
      rep("LFD", 3),
      rep("HFD", 3)
    ))
  )
  ```
- **Jadikan nama sampel sebagai row names**
  ```r
  rownames(metadata) <- metadata$sample
  ```
- **Tampilkan metadata**
  ```r
  metadata
  ```
- **Verifikasi dengan**
  ```r
  colnames(counts) == rownames(metadata)
  semua harus TRUE
  ```
#### Distribusi Data
- **Load library untuk visualisasi**
  ```r
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  ```
- **ubah data ke format panjang (long format) untuk ggplot**
  ```r
  counts_long <- counts %>%
    as.data.frame() %>%
    mutate(gene = rownames(.)) %>%
    pivot_longer(
      cols = -gene,
      names_to = "sample",
      values_to = "count"
    ) %>%
    mutate(type = "Data Mentah")
  ```
- **Boxplot (dalam skala log2 agar lebih jelas)**
  ```r
  ggplot(counts_long, aes(x = sample, y = log2(count + 1), fill = sample)) +
    geom_boxplot() +
    theme_minimal() +
    labs(
      title = "Distribusi Ekspresi Gen per Sampel (Data Mentah)",
      x = "Sampel",
      y = "Log2(Count + 1)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ```
##### a. Filter gen dengan reads rendah (opsional)
- **Filter: gen dengan total reads >= 10 di semua sampel**
  ```r
  keep <- rowSums(counts) >= 10
  counts_filtered <- counts[keep, ]
  ```
- **Bandingkan jumlah gen sebelum dan sesudah filter**
  ```r
  dim(counts)
  dim(counts_filtered)
  ```
#### Korelasi antar sampel
- **install dan load pheatmap (jika belum)**
  ```r
  if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")
  library(pheatmap)
  ```
- **Hitung korelasi antar sampel (pakai data yang sudah difilter)**
  ```r
  cor_matrix <- cor(counts_filtered)
  ```
- **Tampilkan sebagai heatmap**
  ```r
  pheatmap(cor_matrix,
           main = "Korelasi Antar Sampel (Data Mentah)",
           display_numbers = TRUE,
           number_format = "%.2f",
           color = colorRampPalette(c("blue", "white", "red"))(50))
  ```
#### PCA (Principal Component Analysis)
- **PCA dengan data yang sudah di-log (karena data counts sangat skewed)**
  ```r
  log_data <- log2(counts_filtered + 1)
  pca_results <- prcomp(t(log_data), scale. = TRUE)
  ```
- **Buat data frame untuk plotting**
  ```r
  pca_df <- data.frame(
    PC1 = pca_results$x[,1],
    PC2 = pca_results$x[,2],
    condition = metadata$condition
  )
  ```
- **Hitung persentase varians**
  ```r
  var_explained <- summary(pca_results)$importance[2, 1:2]*100
  ```
- **Plot PCA**
  ```
  ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 5) +
    theme_minimal() +
    labs(title = "PCA Plot (Data Mentah - Log2)",
         x = paste0("PC1: ", round(var_explained[1], 1), "% variance"),
         y = paste0("PC2: ", round(var_explained[2], 1), "% variance")) +
    theme(legend.position = "bottom")
  ```
### 02 NORMALIZATION 
#### Normalisasi
- **Load library**
  ```r
  library(DESeq2)
  ```
- **Buat DESeq2 object (dari count matrix dan metadata)**
  ```r
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = metadata,
    design = ~ sex + condition
  )
  ```
- **mulai normalisasi**
  ```r
  dds <- DESeq(dds)
  ```
- **ambil size factors (faktor normalisasi)**
  ```r
  sizeFactors(dds)
  ```
- **ambil data yang sudah dinormalisasi. Data ini bisa digunakan untuk visualisasi (heatmap, PCA)**
  ```r
  rld <- rlog(dds, blind = FALSE)
  vsd <- vst(dds, blind = FALSE)
  ```
#### Persiapan visualisasi
- **ubah data normalisasi (rlog) ke format panjang**
  ```r
  rld_data <- assay(rld)
  rld_long <- rld_data %>%
    as.data.frame() %>%
    mutate(gene = rownames(.)) %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
    mutate(type = "Data Normalisasi (rlog)")
  ```
- **gabungkan kedua data**
  ```r
  combine_data <- bind_rows(
    counts_long %>% mutate(value = log2(count + 1)), 
    rld_long %>% mutate(value = expression)
  )
  ```
#### Visualisasi dengan Boxplot
  ```r
  library(ggplot2)
  ggplot(combine_data, aes(x = sample, y = value, fill = type)) +
    geom_boxplot() +
    facet_wrap(~ type, scales = "free_y", ncol = 2) + 
    theme_minimal() +
    labs(title = "Perbandingan Distribusi Ekspresi Gen", 
         subtitle = "Data Mentah vs Data Normalisasi (rlog)",
         x = "Sample",
         y = "Log2 Ekspresi") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "none")
  ```
### 03 Differential Gene Expression (DEG) Analysis
#### DESeq2
- **membuat object dds**
  ```r
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = metadata,
    design = ~ sex + condition
  )
  ```
- **lihat level condition dan jalankan DESEq2**
  ```r
  dds$condition <- relevel(dds$condition, ref = "LFD")
  dds <- DESeq(dds)
  ```
- **ambil hasil perbandingan**
  ```r
  res <- results(dds,
                 contrast = c("condition", "HFD", "LFD"))
  res <- res[order(res$padj), ]
  head(res)
  ```
#### Filter Gen
- **ubah data jadi format data frame**
  ```r
  res_df <-as.data.frame(res)
  ```
- **tambahkan kolom status signifikansi**
  ```r
  res_df$significant <- ifelse(
    !is.na(res_df$padj)&
      res_df$padj <0.05&
      abs(res_df$log2FoldChange) > 1,
    "Signifikan",
    "Tidak Signifikan"
  )
  ```
- **tambahkan kolom arah perubahan**
  ```r
  res_df$direction <- ifelse(
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1,
    "Naik",
    ifelse(
      !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1,
      "Turun",
      "Tidak Signifikan"
    )
  )
  ```
- **filter gen signifikan**
  ```r
  sig_genes <- res_df[res_df$significant == "Signifikan", ]
  nrow(sig_genes)
  ```
- **lihat gen naik dan turun**
  ```r
  table(sig_genes$direction)
  ```
#### Export Hasil
- **export hasil**
  ```r
  write.csv(res_df, "DEG_results_GSE294660_8WOD.csv")
  ```
- **lihat top 10 gen naik, diurutkan berdasarkan Log2FoldChange**
  ```r
  top_up <- res_df[res_df$direction == "Naik",]
  top_up <- top_up[order(top_up$log2FoldChange,decreasing = TRUE),]
  head(top_up[,c("log2FoldChange","pvalue","padj")],10)
  ```
- **lihat top 10 gen turun**
  ```r
  top_down <- res_df[res_df$direction == "Turun",]
  top_down <- top_down[order(top_down$log2FoldChange),]
  head(top_down[,c("log2FoldChange", "pvalue", "padj")], 10)
  ```
### 04 VISUALIZATION
#### Volcano plot
  ```r
  library(dplyr)
  library(ggplot2)
  ```
- **ambil 10 gen dengan padj terkecil (paling signifikan)**
  top_gene <- res_df %>%
    filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >1) %>%
    arrange(padj) %>%
    head(10)

- **tambahkan kolom tabel (hanya untuk gen top)**
```r
  top_genes <- res_df[order(res_df$padj), ][1:10, ]
  res_df$label <- ifelse(
    rownames(res_df) %in% rownames(top_genes),
    rownames(res_df),
    ""
  )
```
#### Buat volcano plot dengan label
  ```r
  volcano_plot_labeled <- ggplot(
    res_df,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = direction
    )
  ) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c("Naik" = "red","Turun" = "blue","Tidak Signifikan" = "gray"),
      name = "Regulasi"
    ) +
    geom_vline(xintercept = c(-1, 1),linetype = "dashed",color = "black",alpha = 0.5
    ) +
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",color = "black",alpha = 0.5
    ) +
    geom_text(aes(label = label),vjust = -0.5,size = 3,check_overlap = TRUE
    ) +
    labs(title = "Volcano Plot: HFD vs LFD (Top 10 DEG)",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  print(volcano_plot_labeled)
  ggsave(
    "volcano_plot_labeled 8 weeks.png",volcano_plot_labeled, width = 10,height = 7,dpi = 300
  )
  ```
##### Heatmap
  ```r
  library(pheatmap)
  library(DESeq2)
  ```
- **ambil top 30 gen signifikan (berdasarkan padj)**
  ```r
  res_sorted <- res[order(res$padj), ]
  top_genes <- rownames(res_sorted)[!is.na(res_sorted$padj)][1:30]
  ```
- **ambil data ekspresi untuk gen gen tersebut dari rlog**
  ```r
  heatmap_data <- assay(rld)[top_genes, ]
  ```
- **buat annotation untuk sampel (warna berdasarkan kondisi)**
  ```r
  annotation_col <- data.frame(
    Sex = metadata$sex,
    Condition = metadata$condition
  )
  rownames(annotation_col) <- colnames(heatmap_data)
  ```
- **buat heatmap**
  ```r
  heatmap_plot <- pheatmap(
    heatmap_data,
    scale = "row",
    main = "Heatmap Top 30 DEG",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 8,
    annotation_col = annotation_col,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    border_color = NA
  )
  ```

  

