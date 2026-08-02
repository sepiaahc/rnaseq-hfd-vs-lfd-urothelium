Tentang: Pengaruh Obesitas terhadap Infeksi Saluran Kemih (ISK/UTI) melalui treatment 
High Fat Diet (HFD) vs Low Fat Diet (LFD)
menggunakan dataset publik GEO Kode GSE294660

#### 01 DATA PREPARATION 

# Instalasi Package
- Install BiocManager (Pengelola Packages Bioconductor)
install.packages("BiocManager")
- Install semua packages untuk RNA-seq
BiocManager::install(c("DESeq2", "ggplot2", "pheatmap", "GEOquery", "dplyr", "tidyr"))
- Load library setiap kali mau pakai
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(GEOquery)
library(dplyr)
library(tidyr)

# Data Preparation
## Download Data
- Load packages
library(GEOquery)
- download series_matrix GSE

gse <- getGEO("GSE294660", GSEMatrix = TRUE)
- download semua file supplementary dari GSE
getGEOSuppFiles("GSE294660")

## Data Cleaning
- Ambil data ekspresi dan metadata dari series matrix
gse_data<-gse[[1]]
expression_data<-exprs(gse_data)
metadata<-pData(gse_data)

- Cek struktur datanya
dim(expression_data)
- Melihat supplementary file yang sudah didownload 
list.files("GSE294660/")

- Baca isi file, dan simpan sebagai objek "counts". Sesuaikan nama file dengan output list.files tadi
library(readr)
counts <- read.csv(
  "GSE294660_dio_12WOD_urothelium_baseline_counts.csv",
  row.names = 1,
  check.names = FALSE
)

# Ambil hanya sampel female
counts_female <- counts[, grepl("^female", colnames(counts))]

dim(counts_female)

# Membuat metadata
- sederhanakan nama sampel
colnames(counts_female) <- c(
  "F_LFD", "F_LFD", "F_LFD", "F_LFD",
  "F_HFD", "F_HFD", "F_HFD", "F_HFD"
)

-  cek nama sampel
colnames(counts_female)

- membuat metadata
metadata <- data.frame(
  sample = colnames(counts_female),
  sex = factor(rep("Female", 8)),
  condition = factor(c(rep("LFD", 4), rep("HFD", 4)))
)

- Jadikan nama sampel sebagai row names
rownames(metadata) <- metadata$sample

- Tampilkan metadata
metadata

- Verifikasi dengan
colnames(counts_female) == rownames(metadata)
semua harus TRUE

# Distribusi Data
- Load library untuk visualisasi
library(ggplot2)
library(tidyr)
library(dplyr)

- ubah data ke format panjang (long format) untuk ggplot
counts_long <- counts_female %>%
  as.data.frame() %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "count"
  ) %>%
  mutate(type = "Data Mentah")

- Boxplot (dalam skala log2 agar lebih jelas)
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
## Filter gen dengan reads rendah (opsional)
- Filter: gen dengan total reads >= 10 di semua sampel
keep <- rowSums(counts_female) >= 10
counts_female_filtered <- counts_female[keep, ]
- Bandingkan jumlah gen sebelum dan sesudah filter

dim(counts_female)
dim(counts_female_filtered)

# Korelasi antar sampel
- install dan load pheatmap (jika belum)
if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

- Hitung korelasi antar sampel (pakai data yang sudah difilter)
cor_matrix <- cor(counts_female_filtered)

- Tampilkan sebagai heatmap
pheatmap(cor_matrix,
         main = "Korelasi Antar Sampel (Data Mentah)",
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("blue", "white", "red"))(50))

# PCA (Principal Component Analysis)
- PCA dengan data yang sudah di-log (karena data counts sangat skewed)
log_data <- log2(counts_female_filtered + 1)
pca_results <- prcomp(t(log_data), scale. = TRUE)

- Buat data frame untuk plotting
pca_df <- data.frame(
  PC1 = pca_results$x[,1],
  PC2 = pca_results$x[,2],
  condition = metadata$condition
)

- Hitung persentase varians
var_explained <- summary(pca_results)$importance[2, 1:2]*100

- Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 5) +
  theme_minimal() +
  labs(title = "PCA Plot (Data Mentah - Log2)",
       x = paste0("PC1: ", round(var_explained[1], 1), "% variance"),
       y = paste0("PC2: ", round(var_explained[2], 1), "% variance")) +
  theme(legend.position = "bottom")

#### 02 NORMALIZATION (baru sampai sini)
# Normalisasi
- Load library
library(DESeq2)

- Buat DESeq2 object (dari count matrix dan metadata)
dds <- DESeqDataSetFromMatrix(
  countData = counts_female_filtered,
  colData = metadata,
  design = ~ condition
)

- mulai normalisasi
dds <- DESeq(dds)

- ambil size factors (faktor normalisasi)
sizeFactors(dds)

- ambil data yang sudah dinormalisasi. Data ini bisa digunakan untuk visualisasi (heatmap, PCA)
rld <- rlog(dds, blind = FALSE)
vsd <- vst(dds, blind = FALSE)

# Persiapan visualisasi
- ubah data normalisasi (rlog) ke format panjang
rld_data <- assay(rld)
rld_long <- rld_data %>%
  as.data.frame() %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  mutate(type = "Data Normalisasi (rlog)")

- gabungkan kedua data
combine_data <- bind_rows(
  counts_long %>% mutate(value = log2(count + 1)), 
  rld_long %>% mutate(value = expression)
)

# Visualisasi dengan Boxplot
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

#### 03 Differential Gene Expression (DEG) Analysis
# DESeq2
- membuat object dds
dds <- DESeqDataSetFromMatrix(
  countData = counts_female_filtered,
  colData = metadata,
  design = ~ condition
)

- lihat level condition dan jalankan DESEq2
dds$condition <- relevel(dds$condition, ref = "LFD")
dds <- DESeq(dds)

- ambil hasil perbandingan
res <- results(dds,
               contrast = c("condition", "HFD", "LFD"))
res <- res[order(res$padj), ]
head(res)

# Filter Gen
- ubah data jadi format data frame
res_df <-as.data.frame(res)
res_df$gene <- rownames(res_df)

- tambahkan kolom status signifikansi
res_df$significant <- ifelse(
  !is.na(res_df$padj) &
    res_df$padj < 0.05,
  "Signifikan",
  "Tidak Signifikan"
)

- tambahkan kolom arah perubahan
res_df$direction <- ifelse(
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange > 0,
  "Naik",
  ifelse(
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange < 0,
    "Turun",
    "Tidak Signifikan"
  )
)
- Hitung jumlah DEG
up <- sum(res_df$direction=="Naik")

down <- sum(res_df$direction=="Turun")

- filter gen signifikan
sig_genes <- res_df[res_df$significant == "Signifikan", ]
nrow(sig_genes)

- lihat gen naik dan turun
table(sig_genes$direction)

# Export Hasil
- export hasil
write.csv(res_df, "DEG_results_GSE294660_12WOD_female.csv")

- lihat top 10 gen naik, diurutkan berdasarkan Log2FoldChange 

top_up <- res_df[res_df$direction == "Naik",]

top_up <- top_up[order(top_up$log2FoldChange,decreasing = TRUE),]

head(top_up[,c("log2FoldChange","pvalue","padj")],10)


- lihat top 10 gen turun
top_down <- res_df[res_df$direction == "Turun",]

top_down <- top_down[order(top_down$log2FoldChange),]

head(top_down[,c("log2FoldChange", "pvalue", "padj")], 10)

#### 04 VISUALIZATION
# Volcano plot
library(dplyr)
library(ggplot2)

- ambil 10 gen dengan padj terkecil (paling signifikan)
top_gene <- res_df %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >1) %>%
  arrange(padj) %>%
  head(10)

- tambahkan kolom tabel (hanya untuk gen top)
top_genes <- res_df[order(res_df$padj), ][1:10, ]
res_df$label <- ifelse(
  rownames(res_df) %in% rownames(top_genes),
  rownames(res_df),
  ""
)

# Buat volcano plot dengan label
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
    values = c(
      "Naik" = "red",
      "Turun" = "blue",
      "Tidak Signifikan" = "gray"
    ),
    name = "Regulasi"
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.5,
    size = 3,
    check_overlap = TRUE
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
  "volcano_plot_labeled_female 12 weeks.png",
  volcano_plot_labeled,
  width = 10,
  height = 7,
  dpi = 300
)

# Heatmap
library(pheatmap)
library(DESeq2)

- ambil top 30 gen signifikan (berdasarkan padj)
res_sorted <- res[order(res$padj), ]
top_genes <- rownames(res_sorted)[!is.na(res_sorted$padj)][1:30]

- ambil data ekspresi untuk gen gen tersebut dari rlog
heatmap_data <- assay(rld)[top_genes, ]

- buat annotation untuk sampel (warna berdasarkan kondisi)
annotation_col <- data.frame(
  Condition = metadata$condition
)

rownames(annotation_col) <- colnames(heatmap_data)

- buat heatmap
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

