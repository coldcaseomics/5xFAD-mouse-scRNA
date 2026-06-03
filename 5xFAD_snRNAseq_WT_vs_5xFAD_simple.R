# =============================================================================
# 5xFAD snRNA-seq 분석 파이프라인
# 샘플: WT_1 (정상) vs WT_5XFAD_1 (알츠하이머 모델)
# =============================================================================


# ── 0. 패키지 로드 ────────────────────────────────────────────────────────────

library(Seurat)
library(ggplot2)
library(ggrepel)


# ── 1. 경로 설정 ──────────────────────────────────────────────────────────────

setwd("/Users/ijeonglyul/Desktop/5xFAD_project")

data_dir  <- "./data"

dir.create("./plots", showWarnings = FALSE)
dir.create("./rds",   showWarnings = FALSE)


# ── 2. WT 샘플 로드 ───────────────────────────────────────────────────────────

# 임시 폴더 만들기
tmp_wt <- file.path(tempdir(), "WT_1")
dir.create(tmp_wt, showWarnings = FALSE)

# 파일 복사 (이름을 10X 표준 형식으로 변경)
file.copy(from = "./data/GSM4173504_WT_1_barcodes.tsv.gz",
          to   = file.path(tmp_wt, "barcodes.tsv.gz"), overwrite = TRUE)
file.copy(from = "./data/GSM4173504_WT_1_features.tsv.gz",
          to   = file.path(tmp_wt, "features.tsv.gz"), overwrite = TRUE)
file.copy(from = "./data/GSM4173504_WT_1_matrix.mtx.gz",
          to   = file.path(tmp_wt, "matrix.mtx.gz"), overwrite = TRUE)

# count matrix 읽기
counts_wt <- Read10X(data.dir = tmp_wt)

# Seurat object 생성
obj_wt <- CreateSeuratObject(counts = counts_wt, project = "WT_1",
                              min.cells = 3, min.features = 200)

# 메타데이터 추가
obj_wt$sample_id <- "WT_1"
obj_wt$genotype  <- "WT"

# 미토콘드리아 % 계산 (마우스: 소문자 mt-)
obj_wt[["percent.mt"]] <- PercentageFeatureSet(obj_wt, pattern = "^mt-")

# QC 필터링
obj_wt <- subset(obj_wt,
                 subset = nFeature_RNA > 200 &
                          nFeature_RNA < 6000 &
                          nCount_RNA   > 500  &
                          percent.mt   < 5)

# 확인
obj_wt


# ── 3. 5xFAD 샘플 로드 ────────────────────────────────────────────────────────

# 임시 폴더 만들기
tmp_5xfad <- file.path(tempdir(), "WT_5XFAD_1")
dir.create(tmp_5xfad, showWarnings = FALSE)

# 파일 복사
file.copy(from = "./data/GSM4173510_WT_5XFAD_1_barcodes.tsv.gz",
          to   = file.path(tmp_5xfad, "barcodes.tsv.gz"), overwrite = TRUE)
file.copy(from = "./data/GSM4173510_WT_5XFAD_1_features.tsv.gz",
          to   = file.path(tmp_5xfad, "features.tsv.gz"), overwrite = TRUE)
file.copy(from = "./data/GSM4173510_WT_5XFAD_1_matrix.mtx.gz",
          to   = file.path(tmp_5xfad, "matrix.mtx.gz"), overwrite = TRUE)

# count matrix 읽기
counts_5xfad <- Read10X(data.dir = tmp_5xfad)

# Seurat object 생성
obj_5xfad <- CreateSeuratObject(counts = counts_5xfad, project = "WT_5XFAD_1",
                                 min.cells = 3, min.features = 200)

# 메타데이터 추가
obj_5xfad$sample_id <- "WT_5XFAD_1"
obj_5xfad$genotype  <- "WT_5XFAD"

# 미토콘드리아 % 계산
obj_5xfad[["percent.mt"]] <- PercentageFeatureSet(obj_5xfad, pattern = "^mt-")

# QC 필터링
obj_5xfad <- subset(obj_5xfad,
                    subset = nFeature_RNA > 200 &
                             nFeature_RNA < 6000 &
                             nCount_RNA   > 500  &
                             percent.mt   < 5)

# 확인
obj_5xfad


# ── 4. QC 시각화 ──────────────────────────────────────────────────────────────

# WT QC violin plot
p_qc_wt <- VlnPlot(obj_wt,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3, pt.size = 0)
ggsave("./plots/QC_WT.jpg", plot = p_qc_wt, width = 12, height = 5)

# 5xFAD QC violin plot
p_qc_5xfad <- VlnPlot(obj_5xfad,
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       ncol = 3, pt.size = 0)
ggsave("./plots/QC_5xFAD.jpg", plot = p_qc_5xfad, width = 12, height = 5)


# ── 5. 두 샘플 합치기 ─────────────────────────────────────────────────────────

seurat_merged <- merge(
  obj_wt,
  y            = obj_5xfad,
  add.cell.ids = c("WT", "5XFAD")
)

seurat_merged


# ── 6. 정규화 (SCTransform) ───────────────────────────────────────────────────

options(future.globals.maxSize = Inf)

seurat_merged <- SCTransform(
  seurat_merged,
  vst.flavor      = "v2",
  vars.to.regress = "percent.mt",
  verbose         = FALSE
)


# ── 7. PCA ────────────────────────────────────────────────────────────────────

seurat_merged <- RunPCA(seurat_merged, npcs = 50, verbose = FALSE)

# ElbowPlot 확인 후 PC 수 결정
p_elbow <- ElbowPlot(seurat_merged, ndims = 50)
p_elbow
ggsave("./plots/elbow_plot.jpg", plot = p_elbow, width = 8, height = 5)

# ElbowPlot 보고 평평해지는 지점 선택 → 여기서는 30 사용


# ── 8. UMAP ───────────────────────────────────────────────────────────────────

seurat_merged <- RunUMAP(seurat_merged, dims = 1:30, verbose = FALSE)


# ── 9. 클러스터링 ─────────────────────────────────────────────────────────────

seurat_merged <- FindNeighbors(seurat_merged, dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)

# 클러스터 수 확인
table(seurat_merged$seurat_clusters)

# UMAP 시각화
p_umap <- DimPlot(seurat_merged, label = TRUE, label.size = 4)
p_umap
ggsave("./plots/umap_clusters.jpg", plot = p_umap, width = 10, height = 8)


# ── 10. 세포 타입 마커 확인 ───────────────────────────────────────────────────

seurat_merged <- PrepSCTFindMarkers(seurat_merged)
all_markers <- FindAllMarkers(
  seurat_merged,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25
)

# 파일저장
all_markers$p_val <- format(all_markers$p_val, scientific = TRUE)
all_markers$p_val_adj <- format(all_markers$p_val_adj, scientific = TRUE)

write.csv(all_markers, "./rds/all_cluster_markers.csv", row.names = FALSE)

library(dplyr)
top5 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(top5, n = 105)

# 각 클러스터별로 유전자 찾기 (어떤 세포인지 판단하기 위함!)
all_markers[all_markers$cluster == "0", ] %>%
head(20)


cluster_annotation <- c(
  "0"  = "Oligodendrocyte",
  "1"  = "Excitatory_Neuron",
  "2"  = "Excitatory_Neuron",
  "3"  = "Excitatory_Neuron",
  "4"  = "Excitatory_Neuron",
  "5"  = "Inhibitory_Neuron",
  "6"  = "Inhibitory_Neuron",
  "7"  = "Microglia",
  "8"  = "Astrocyte",
  "9"  = "Inhibitory_Neuron",
  "10" = "Unknown",
  "11" = "Inhibitory_Neuron",
  "12" = "Inhibitory_Neuron",
  "13" = "Inhibitory_Neuron",
  "14" = "OPC",
  "15" = "Excitatory_Neuron",
  "16" = "Vascular",
  "17" = "Excitatory_Neuron",
  "18" = "Excitatory_Neuron",
  "19" = "Oligodendrocyte",
  "20" = "OPC"
)

seurat_merged <- RenameIdents(seurat_merged, cluster_annotation)
seurat_merged$cell_type <- Idents(seurat_merged)
table(seurat_merged$cell_type)

p_celltype <- DimPlot(seurat_merged,
                      group.by   = "cell_type",
                      label      = TRUE,
                      label.size = 4,
                      repel      = TRUE)
p_celltype

ggsave("./plots/umap_celltype.jpg", plot = p_celltype, width = 10, height = 8)


p_split <- DimPlot(seurat_merged,
                   group.by  = "cell_type",
                   split.by  = "genotype",
                   label     = TRUE,
                   label.size = 3,
                   repel     = TRUE)

p_split
ggsave("./plots/umap_wt_vs_5xfad.jpg", plot = p_split, width = 14, height = 6)

prop_table <- prop.table(
  table(seurat_merged$genotype, seurat_merged$cell_type),
  margin = 1) * 100

round(prop_table, 2)


prop_data <- as.data.frame(prop_table)
colnames(prop_data) <- c("genotype", "cell_type", "proportion")

p_prop <- ggplot(prop_data, aes(x = cell_type, y = proportion, fill = genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell type proportion: WT vs 5xFAD",
       x = "Cell type", y = "Proportion (%)")

p_prop
ggsave("./plots/celltype_proportion.jpg", plot = p_prop, width = 8, height = 5)


#DEG 분석!
# microglia 부터
microglia <- subset(seurat_merged, subset = cell_type == "Microglia")

DefaultAssay(microglia) <- "RNA"
microglia <- NormalizeData(microglia)
microglia <- JoinLayers(microglia)

Idents(microglia) <- "genotype"

deg_microglia <- FindMarkers(
  microglia,
  ident.1  = "WT_5XFAD",
  ident.2  = "WT",
  test.use = "wilcox",
  min.pct  = 0.1
)

deg_microglia$gene <- rownames(deg_microglia)
head(deg_microglia, 20)



p_volcano_mg <- ggplot(deg_microglia,
                       aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5, size = 1, color = "gray") +
  geom_point(data = subset(deg_microglia, p_val_adj < 0.05 & avg_log2FC >  0.5),
             color = "#E24B4A", size = 1.5) +
  geom_point(data = subset(deg_microglia, p_val_adj < 0.05 & avg_log2FC < -0.5),
             color = "#185FA5", size = 1.5) +
  geom_text_repel(
    data = subset(deg_microglia, p_val_adj < 1e-10 & abs(avg_log2FC) > 1),
    aes(label = gene), size = 3, max.overlaps = 15) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "Microglia DEG: 5xFAD vs WT",
       x = "avg log2FC", y = "-log10(adj.p)") +
  theme_classic()

p_volcano_mg
ggsave("./plots/volcano_microglia.jpg", plot = p_volcano_mg, width = 8, height = 6)


# astrocyte

astrocyte <- subset(seurat_merged, subset = cell_type == "Astrocyte")

DefaultAssay(astrocyte) <- "RNA"
astrocyte <- NormalizeData(astrocyte)
astrocyte <- JoinLayers(astrocyte)

Idents(astrocyte) <- "genotype"

deg_astrocyte <- FindMarkers(
  astrocyte,
  ident.1  = "WT_5XFAD",
  ident.2  = "WT",
  test.use = "wilcox",
  min.pct  = 0.1
)

deg_astrocyte$gene <- rownames(deg_astrocyte)
head(deg_astrocyte, 20)

p_volcano_ast <- ggplot(deg_astrocyte,
                        aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5, size = 1, color = "gray") +
  geom_point(data = subset(deg_astrocyte, p_val_adj < 0.05 & avg_log2FC >  0.5),
             color = "#E24B4A", size = 1.5) +
  geom_point(data = subset(deg_astrocyte, p_val_adj < 0.05 & avg_log2FC < -0.5),
             color = "#185FA5", size = 1.5) +
  geom_text_repel(
    data = subset(deg_astrocyte, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5),
    aes(label = gene), size = 3, max.overlaps = 15) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  xlim(-2.5, 3) +
  labs(title = "Astrocyte DEG: 5xFAD vs WT",
       x = "avg log2FC", y = "-log10(adj.p)") +
  theme_classic()

p_volcano_ast
ggsave("./plots/volcano_astrocyte.jpg", plot = p_volcano_ast, width = 8, height = 6)


# 최종 저장
write.csv(deg_microglia,  "./rds/DEG_microglia_5xFAD_vs_WT.csv",  row.names = FALSE)
write.csv(deg_astrocyte,  "./rds/DEG_astrocyte_5xFAD_vs_WT.csv",  row.names = FALSE)
saveRDS(seurat_merged,    "./rds/GSE140510_2samples_final.rds")

message("=== 분석 완료 ===")


