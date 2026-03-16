
library(survminer)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(gridExtra)
library(viridis)
library(survival)
library(ggpubr)
library(broom)
library(forestplot)
library(stringr)
library(GSEABase)
library(patchwork)

sci_if_needed <- function(x) {
  sapply(x, function(v) {
    
    if (abs(v) > 100 | abs(v) < 0.01) {
      formatC(v, format = "e", digits = 2)
    } else {
      sprintf("%.2f", v)
    }
  })
}

threshold_continous = function(nums,min=0,max=300){
  nums[nums<min] = min
  nums[nums>max] = max
  return(nums)
}
options(future.globals.maxSize = 100 * 1000 * 1024^2)
SaveDir = "/mnt/Data1/groupjuan/Ming/PanBrainCancer/Results/step5_Analysis_scRNAseq_Tcell/Output"
DataDir = "/mnt/Data1/groupjuan/Ming/PanBrainCancer/Results/step4_AnnotationFullmeta/Output"
gmtDir = "/mnt/Data1/groupjuan/Ming/PanBrainCancer/Results/step5_Analysis_scRNAseq_Tcell/Data"
dir(gmtDir)

sample_meta_scRNAseq = read.csv2(file = file.path(gmtDir,"scRNAseq_cohort_sample_meta.csv"),sep = ",")
cell_meta_scRNAseq = readRDS(file = file.path(DataDir,"full_meta.rds"))
rownames(sample_meta_scRNAseq) = sample_meta_scRNAseq$sample_id
# table(cell_meta_scRNAseq$dataset)
#This analysis is for CD8 Tcell infiltration paper. Only GBM sample is used
table(sample_meta_scRNAseq$condition_top)
sample_meta_scRNAseq = sample_meta_scRNAseq[sample_meta_scRNAseq$condition_top == "Glioblastoma",]
cell_meta_scRNAseq = cell_meta_scRNAseq[cell_meta_scRNAseq$orig.ident %in% names(table(sample_meta_scRNAseq$sample_id)),]

#check annotation dot plot####
TcellDir = "/mnt/Data1/groupjuan/Ming/PanBrainCancer/Results/step5_Analysis_scRNAseq_Tcell/Data"
obj_immune = readRDS(file = file.path(TcellDir,"immune_integrated_harmony_sketch.rds"))
obj_immune = obj_immune[,(obj_immune$dataset %in% c("GSE165080","GSE174554"))]
gc()
obj_nonimmune = readRDS(file = file.path(TcellDir,"nonimmune_integrated_harmony_sketch.rds"))
obj_nonimmune = obj_nonimmune[,(obj_nonimmune$dataset %in% c("GSE165080","GSE174554"))]
gc()
obj_immune<- JoinLayers(obj_immune)
obj_nonimmune<- JoinLayers(obj_nonimmune)
obj_immune@active.assay = "RNA"
obj_nonimmune@active.assay = "RNA"
common_genes = rownames(obj_immune)[rownames(obj_immune) %in% rownames(obj_nonimmune)]
obj_immune = obj_immune[common_genes,]
obj_nonimmune = obj_nonimmune[common_genes,]
gc()

all(rownames(obj_immune) == rownames(obj_nonimmune))
obj = merge(obj_immune, y = obj_nonimmune)


obj@meta.data[obj$united_annotation %in% c("CD4_Naive_Tcells","CD4_Treg"),"general_annotation"]= "CD4Tcell"
obj@meta.data[obj$united_annotation %in% c("CD8_Effector_Tcells","CD8_Effector_Tcells_Exhausted",
                                           "CD8_Effector_Tcells_Exhausted_Proliferating_IDH_TP53",
                                           "CD8_Naive_Tcells"),"general_annotation"]= "CD8Tcell"
obj@meta.data[obj$united_annotation %in% c("DN_NKT","DN_MAIT"),"general_annotation"]= "DN_T"

obj$general_annotation2 = obj$general_annotation
obj@meta.data[obj$united_annotation %in% c("Astrocytes"),"general_annotation2"]= "Astrocytes"
obj@meta.data[obj$united_annotation %in% c("Oligodendrocytes"),"general_annotation2"]= "Oligodendrocytes"
obj@meta.data[obj$united_annotation %in% c("OPC"),"general_annotation2"]= "OPC"
obj@meta.data[obj$united_annotation %in% c("Endothelials"),"general_annotation2"]= "Endothelials"
obj@meta.data[obj$united_annotation %in% c("Fibroblasts"),"general_annotation2"]= "Fibroblasts"

obj$general_annotation2 = factor(obj$general_annotation2,levels = c("CD4Tcell","CD8Tcell","DN_T","NKcells","ILCs","Bcells","Myeloids",
                                                                    "Endothelials","Fibroblasts","Tumorcells","Astrocytes","OPC","Oligodendrocytes","Neurons"))

marker_sets = getGmt(file.path(gmtDir,"Presentation_markers_simple.gmt"),sep = ",")
all_genes = rownames(obj)
p_list = list()
width_list = c()
Dotplot_with_label = function(object, features, geneset_name, 
                              geom_segment_x = 0.5, geom_segment_xend = 9, 
                              annotate_x = 4, show_y_axis = F,group.by="immune_based_clusters"){
  
  colour_paletts = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026")
  
  p <- DotPlot(object = object,features = features, group.by=group.by, scale.by = "size", scale=T) +
    RotatedAxis() +
    scale_color_gradientn(values = seq(0,1,0.1),colours = colour_paletts)+ 
    xlab(label = "") + ylab(label = "") + 
    scale_size_continuous(range = c(0.5, 5),breaks = c(10, 20, 30, 40)) + 
    geom_segment(x=geom_segment_x, xend=geom_segment_xend,
                 y=-1.2,yend=-1.2,linewidth = 0.45) + 
    annotate(geom="text",x = annotate_x, y = -Inf,label=geneset_name, hjust = 1, vjust = 8, angle = 45, size=6.3) + 
    coord_cartesian(clip="off")
  # labs(title = "Dotplot") + 
  if(show_y_axis){
    p = p + theme( plot.margin = margin(r = 0,b = 30) ,
                   plot.title = element_text(size = 15),
                   axis.text = element_text(size = 14) )
  }else{
    p = p + theme(
      axis.text.y = element_blank(),  #
      axis.ticks.y = element_blank(), # 
      axis.title.y = element_blank(),  # 
      axis.line.y = element_blank(),  # 
      axis.text = element_text(size = 14),
      plot.margin = margin(l = 0,b = 30))
  }
  p = p + theme( legend.position = "bottom")
  return(p)
}

for (i in 1:length(marker_sets)) {
  print(marker_sets[[i]]@setName)
  features = marker_sets[[i]]@geneIds[marker_sets[[i]]@geneIds != ""]
  features = features[features %in% all_genes]
  width_list = c(width_list, length(features))
  if(i == 1){
    p_list[[i]] = Dotplot_with_label(object = obj, features = features, geneset_name = marker_sets[[i]]@setName,
                                     geom_segment_x = 0.5, geom_segment_xend = length(features) , annotate_x = length(features) / 6, 
                                     show_y_axis = T,group.by="general_annotation2")
  }else{
    p_list[[i]] = Dotplot_with_label(object = obj, features = features, geneset_name = marker_sets[[i]]@setName,
                                     geom_segment_x = 0.5, geom_segment_xend = length(features) , annotate_x = length(features) / 6, 
                                     show_y_axis = F,group.by="general_annotation2")
  }
}
{
  p <- (p_list[[1]] + p_list[[2]] + p_list[[3]] + p_list[[4]] + 
          p_list[[5]] + p_list[[6]] + p_list[[7]] + p_list[[8]] + 
          p_list[[9]] + p_list[[10]] + p_list[[11]] + p_list[[12]] + 
          p_list[[13]] + p_list[[14]] + p_list[[15]]) + 
    plot_layout(guides = "collect", widths = width_list[1:38]) &
    theme( legend.position = "bottom", 
           legend.margin = margin(t = 26, r = 0, b = 0, l = 0) ) &
    guides(
      color = guide_colorbar(title.position = "top", title.hjust = 0.5,title = "Avg.exp."),
      size = guide_legend(title.position = "top", title.hjust = 0.5,
                          label.position = "bottom", 
                          direction = "horizontal",
                          nrow = 1,
                          keywidth = unit(0.5, "cm"),
                          keyheight = unit(0.8, "cm"),
                          override.aes = list(color = "black")) ) 
  pdf(file = file.path(SaveDir,paste0("dotplot_allcelltypes_","general_annotation2","_Presentation_simple.pdf")),width = 30,height = 6)
  p
  dev.off()
}

# retrive Healthy control data ####
sample_meta_scRNAseq = read.csv2(file = file.path(gmtDir,"scRNAseq_cohort_sample_meta.csv"),sep = ",")
cell_meta_scRNAseq = readRDS(file = file.path(DataDir,"full_meta.rds"))
rownames(sample_meta_scRNAseq) = sample_meta_scRNAseq$sample_id
table(sample_meta_scRNAseq$condition_top)
sample_meta_scRNAseq = sample_meta_scRNAseq[sample_meta_scRNAseq$condition_top == "Healthy",]
cell_meta_scRNAseq = cell_meta_scRNAseq[cell_meta_scRNAseq$orig.ident %in% names(table(sample_meta_scRNAseq$sample_id)),]

df = cell_meta_scRNAseq
df[df$united_annotation %in% c("CD4_Naive_Tcells","CD4_Treg"),"general_annotation"] = "CD4Tcell"
df[df$united_annotation %in% c("CD8_Effector_Tcells","CD8_Effector_Tcells_Exhausted",
                               "CD8_Effector_Tcells_Exhausted_Proliferating_IDH_TP53",
                               "CD8_Naive_Tcells"),"general_annotation"] = "CD8Tcell"
df[df$united_annotation %in% c("DN_NKT","DN_MAIT"),"general_annotation"] = "DN_T"
table(df$general_annotation)
df_count = df %>% 
  group_by(orig.ident,general_annotation) %>%
  summarise(n = n())
df_percent = df_count  %>% 
  group_by(orig.ident) %>%
  summarise(all_count = sum(n))

df_count = as.data.frame(df_count)
df_count[((df_count$orig.ident == df_percent[i,"orig.ident"][[1]])&(df_count$general_annotation=="CD8Tcell")),"n"]
df_count = df_count[df_count$general_annotation =="CD8Tcell",]
rownames(df_count)=df_count$orig.ident
df_percent = as.data.frame(df_percent)
rownames(df_percent)=df_percent$orig.ident
df_percent[rownames(df_count),"CD8Tcell_count"] = df_count$n
df_percent$CD8Tcell_percent = df_percent$CD8Tcell_count / df_percent$all_count


sample_meta_scRNAseq[rownames(df_percent),"CD8Tcell_percent"] = df_percent$CD8Tcell_percent
sample_meta_scRNAseq$CD8Tcell_cat = NA
saveRDS(sample_meta_scRNAseq,file.path(SaveDir,"sample_meta_scRNAseq_Healthy.rds"))
saveRDS(cell_meta_scRNAseq,file.path(SaveDir,"cell_meta_scRNAseq_Healthy.rds"))


# category CD8high/low infiltration####
df = cell_meta_scRNAseq
df[df$united_annotation %in% c("CD4_Naive_Tcells","CD4_Treg"),"general_annotation"] = "CD4Tcell"
df[df$united_annotation %in% c("CD8_Effector_Tcells","CD8_Effector_Tcells_Exhausted","CD8_Effector_Tcells_Exhausted_Proliferating_IDH_TP53"),"general_annotation"] = "CD8Tcell"
df[df$united_annotation %in% c("DN_NKT"),"general_annotation"] = "DN_NKT"
table(df$general_annotation)
colnames(cell_meta_scRNAseq)
df_count = df %>% 
  group_by(orig.ident,general_annotation) %>%
  summarise(n = n())
df_percent = df_count  %>% 
  group_by(orig.ident) %>%
  summarise(all_count = sum(n))

df_count = as.data.frame(df_count)
df_count[((df_count$orig.ident == df_percent[i,"orig.ident"][[1]])&(df_count$general_annotation=="CD8Tcell")),"n"]
df_count = df_count[df_count$general_annotation =="CD8Tcell",]
rownames(df_count)=df_count$orig.ident
df_percent = as.data.frame(df_percent)
rownames(df_percent)=df_percent$orig.ident
df_percent[rownames(df_count),"CD8Tcell_count"] = df_count$n
df_percent$CD8Tcell_percent = df_percent$CD8Tcell_count / df_percent$all_count


sample_meta_scRNAseq[rownames(df_percent),"CD8Tcell_percent"] = df_percent$CD8Tcell_percent
pdf(file = file.path(SaveDir,"hist_CD8Tcell_percent.pdf"),width = 5,height = 5)
print(hist(sample_meta_scRNAseq$CD8Tcell_percent,breaks = 30))
dev.off()

sample_meta_scRNAseq$CD8Tcell_percent_100 = sample_meta_scRNAseq$CD8Tcell_percent * 100
p = ggplot(sample_meta_scRNAseq, aes(x = CD8Tcell_percent_100)) +
  # 1. Increase bins for more detail; set border color for clarity
  geom_histogram(bins = 100, fill = "steelblue", color = "white", alpha = 0.8) +
  
  # 2. Add the vertical line (using linewidth for modern ggplot2)
  geom_vline(aes(xintercept = 0.41061), color = "firebrick", linetype = "dashed", linewidth = 1) +
  
  # 3. Add a label for the line
  annotate("text", x = 0.41061 + 0.2, y = 55, label = paste("Median:", round(0.41061, 2),"%"), 
           color = "firebrick", fontface = "bold", hjust = 0) +
  coord_cartesian(ylim = c(0, 15))+
  # 4. Apply a clean theme and labels
  theme_minimal() +
  labs(
    title = "Distribution of Percent of CD8 TIL",
    x = "Percent of CD8 TIL",
    y = "Frequency"
  )
pdf(file = file.path(SaveDir,"hist_CD8Tcell_percent_100.pdf"),width = 5,height = 4)
print(p)
dev.off()

sample_meta_scRNAseq$CD8Tcell_cat = "Low"
sample_meta_scRNAseq[sample_meta_scRNAseq$CD8Tcell_percent>= 0.0041061 ,"CD8Tcell_cat"] = "High"

saveRDS(sample_meta_scRNAseq,file.path(SaveDir,"sample_meta_scRNAseq_GBM.rds"))
saveRDS(cell_meta_scRNAseq,file.path(SaveDir,"cell_meta_scRNAseq_GBM.rds"))

sample_meta_scRNAseq$PFS_time
tmp = sample_meta_scRNAseq[,c("patients_id","PRM","CD8Tcell_percent","OS_time","PFS_time")]

# survival analysis ####

sample_meta = readRDS(file.path(SaveDir,"sample_meta_scRNAseq_GBM.rds"))
cell_meta = readRDS(file.path(SaveDir,"cell_meta_scRNAseq_GBM.rds"))

sample_meta$CD8Tcell_percent = sample_meta$CD8Tcell_percent * 100 # so that the interpretation of regression coef is change per 1 %
#CD8Tcell_percent
model = coxph(Surv(OS_time, OS_status) ~ CD8Tcell_percent + Gender + Age, data=sample_meta[sample_meta$PRM == "Primary",])
{
  df <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      Variable = recode(term,
                        "CD8Tcell_percent" = "CD8 percentage",
                        "Age" = "Age (years)",
                        "GenderMale" = "Male (vs Female)"
      ),
      HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      p_value = sprintf("%.3f", p.value)
    ) %>%
    select(Variable, HR_CI, p_value, estimate, conf.low, conf.high)
  
  df <- df %>%
    mutate(
      HR_CI = sprintf(
        "%s (%s–%s)",
        sci_if_needed(estimate),
        sci_if_needed(conf.low),
        sci_if_needed(conf.high)
      )
    )
  table_text <- rbind(
    c("Variable", "HR (95% CI)", "p"),
    as.matrix(df[, c("Variable", "HR_CI", "p_value")])
  )
  
  p = forestplot(
    labeltext = table_text,
    mean  = c(NA, df$estimate),
    lower = c(NA, df$conf.low),
    upper = c(NA, df$conf.high),
    
    zero = 1,
    graph.pos = 2,
    xlog = TRUE,
    
    boxsize = 0.20,
    lwd.ci = 1.8,
    ci.vertices = TRUE,
    ci.vertices.height = 0.06,
    
    col = fpColors(
      box  = "#1A4E8A",     
      lines = "#1A4E8A",
      zero = "#333333"    
    ),
    
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.20, fontfamily = "Helvetica"),
      ticks = gpar(cex = 1.10, fontfamily = "Helvetica"),
      xlab  = gpar(cex = 1.25, fontfamily = "Helvetica")
    ),
    
    lineheight = unit(6.5, "mm"),
    colgap = unit(7, "mm"),
    
    xticks = c(0.2, 0.5, 1, 2, 5, 20),
    
    xlab = "Hazard Ratio (log scale)",
    title ="GBM OS Primary cohort"
  )
  pdf(file = paste0(SaveDir,"/Forestplot_Primarycohort_OS_CD8percent_fullmodel.pdf"),width = 8,height = 6)
  print(p)
  dev.off()
  dev.off()
}

library(survminer)
library(ggsurvfit)
library(survival)
data = sample_meta[sample_meta$PRM == "Primary",]
fit <- survfit(Surv(OS_time, OS_status) ~ CD8Tcell_cat, data = data)
p = ggsurvplot(
    fit, 
    data = data, 
    size = 1,                 # change line size
    palette = c("#3ec1d3", "#ff9a00"),# custom color palettes
    # conf.int = TRUE,          # Add confidence interval
    pval = TRUE,              # Add p-value
    risk.table = TRUE,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c("Low", "High"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw(),      # Change ggplot2 theme
    title = "Overall Survival in Primary GBM",
    subtitle = "group by CD8 tensity"
  )
pdf(file = paste0(SaveDir,"/SurvCurve_OS_CD8_PrimaryGBM.pdf"),width = 4,height = 5,onefile = F)
print(p)
dev.off()

# check CD8Tcell phenotype ####
## prepare meta dfs and seurat obj ####
table(cell_meta_GBM$united_annotation)
sample_meta_GBM = readRDS(file.path(SaveDir,"sample_meta_scRNAseq_GBM.rds"))
cell_meta_GBM = readRDS(file.path(SaveDir,"cell_meta_scRNAseq_GBM.rds"))
sample_meta_Healthy = readRDS(file.path(SaveDir,"sample_meta_scRNAseq_Healthy.rds"))
cell_meta_Healthy = readRDS(file.path(SaveDir,"cell_meta_scRNAseq_Healthy.rds"))
TcellDir = "/mnt/Data1/groupjuan/Ming/PanBrainCancer/Results/step5_Analysis_scRNAseq_Tcell/Data"
obj = readRDS(file = file.path(TcellDir,"obj_Tcells.rds"))

sample_meta = rbind(sample_meta_GBM,sample_meta_Healthy)
cell_meta = rbind(cell_meta_GBM,cell_meta_Healthy)
df = cell_meta
df[df$united_annotation %in% c("CD4_Naive_Tcells","CD4_Treg"),"general_annotation"] = "CD4Tcell"
df[df$united_annotation %in% c("CD8_Effector_Tcells","CD8_Effector_Tcells_Exhausted",
                               "CD8_Effector_Tcells_Exhausted_Proliferating_IDH_TP53",
                               "CD8_Naive_Tcells"),"general_annotation"] = "CD8Tcell"
df[df$united_annotation %in% c("DN_NKT","DN_MAIT"),"general_annotation"] = "DN_T"
cell_meta = df
cell_meta = cell_meta[(cell_meta$general_annotation%in%c("CD8Tcell","CD4Tcell")),]
obj = obj[,cell_meta$fullname]
gc()
saveRDS(obj, file = file.path(SaveDir,"obj_CD84Tcell.rds"))
saveRDS(cell_meta, file = file.path(SaveDir,"cell_meta.rds"))
saveRDS(sample_meta, file = file.path(SaveDir,"sample_meta.rds"))


## process seurat obj for T cells####
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$dataset)
obj <- SetIdent(obj, value = "dataset")#some samples only has few cells. Cannot integrate at sample level. integrate at dataset level instead.

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
gc()
# visualize before intedgration
pdf(file = file.path(SaveDir,"umap.unintegrated_dataset.pdf"),width = 12,height = 10)
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("dataset"),label = T,repel = T)
dev.off()
length(table(obj$orig.ident))#[1] 223
pdf(file = file.path(SaveDir,"umap.unintegrated_orig.ident.pdf"),width = 16,height = 10)
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("orig.ident"),label = T,repel = T)
dev.off()

#integration
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
# obj <- FindClusters(obj, resolution = 0.5, cluster.name = "harmony_clusters")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
obj<- JoinLayers(obj)
obj$general_annotation = cell_meta[colnames(obj),"general_annotation"]

obj$condition = obj$dataset
obj@meta.data[obj$dataset=="GSE165080","condition"]="HealthyPBMC"
obj@meta.data[obj$dataset=="GSE174554","condition"]="GBMTissue"
obj$condition = paste0(obj$condition,"_",obj$general_annotation)

saveRDS(obj,file = file.path(SaveDir,"obj_CD84Tcell.rds"))


## visualize processed seurat obj ####

obj = readRDS(file = file.path(SaveDir,"obj_CD84Tcell.rds"))
obj@active.assay = "sketch"
pdf(file = file.path(SaveDir,"umap.integrated.harmony_dataset.pdf"),width = 12,height = 10)
DimPlot(obj, reduction = "umap.harmony", group.by = c("dataset"),label = T,repel = T)
dev.off()
pdf(file = file.path(SaveDir,"umap.integrated.harmony_orig.ident.pdf"),width = 16,height = 10)
DimPlot(obj, reduction = "umap.harmony", group.by = c("orig.ident"),label = T,repel = T)
dev.off()
pdf(file = file.path(SaveDir,"umap.integrated.harmony_general_annotation.pdf"),width = 12,height = 10)
DimPlot(obj, reduction = "umap.harmony", group.by = c("general_annotation"),label = T,repel = T)
dev.off()
pdf(file = file.path(SaveDir,"umap.integrated.harmony_condition.pdf"),width = 12,height = 10)
DimPlot(obj, reduction = "umap.harmony", group.by = c("condition"),label = T,repel = T)
dev.off()

obj@active.assay = "sketch"
pdf(file = file.path(SaveDir,"FeaturePlot_lineagemarker.pdf"),width = 15,height = 15)
FeaturePlot(obj, features = c("PTPRC","CD3D","CD3E","CD3G","CD247","CD8A","CD8B","CD4"),pt.size = 0.1)
dev.off()
pdf(file = file.path(SaveDir,"FeaturePlot_functionmarker.pdf"),width = 20,height = 15)
FeaturePlot(obj, features = c("GZMB","GZMH","GZMK","PRF1","PDCD1","HAVCR2","CTLA4","LAG3","CD274","TIGIT"),pt.size = 0.1)
dev.off()

obj@active.assay = "RNA"
pdf(file = file.path(SaveDir,"VlnPlot_functionmarker.pdf"),width = 15,height = 15)
VlnPlot(obj, features = c("GZMB","GZMH","GZMK","PRF1","PDCD1","HAVCR2","CTLA4","LAG3","CD274","TIGIT"),group.by = c("condition"))
dev.off()
pdf(file = file.path(SaveDir,"RidgePlot_functionmarker.pdf"),width = 15,height = 15)
RidgePlot(obj, features = c("GZMB","GZMH","GZMK","PRF1","PDCD1","HAVCR2","CTLA4","LAG3","CD274","TIGIT"),group.by = c("condition"), ncol = 2)
dev.off()


## calculate positive cell percentage ####
obj = readRDS(file = file.path(SaveDir,"obj_CD84Tcell.rds"))
cell_meta = readRDS(file = file.path(SaveDir,"cell_meta.rds"))

obj@active.assay = "RNA"
pheno_df = FetchData(object = obj,c("condition","GZMB","GZMH","GZMK","PRF1","PDCD1","HAVCR2","CTLA4","LAG3","CD274","TIGIT"))
all(rownames(cell_meta) %in% rownames(pheno_df))
pheno_df = pheno_df[rownames(cell_meta),]
cell_meta = cbind(cell_meta,pheno_df)
saveRDS(cell_meta, file = file.path(SaveDir,"cell_meta.rds"))

table(cell_meta$dataset)
Samples_Healthy = names(table(cell_meta[cell_meta$dataset=="GSE165080","orig.ident"]))
Samples_GBM = names(table(cell_meta[cell_meta$dataset=="GSE174554","orig.ident"]))
Celltype = names(table(cell_meta$condition))
# [1] "GBMTissue_CD4Tcell"   "GBMTissue_CD8Tcell"   "HealthyPBMC_CD4Tcell" "HealthyPBMC_CD8Tcell"
feature = c("GZMB","GZMH","GZMK","PRF1","PDCD1","HAVCR2","CTLA4","LAG3","CD274","TIGIT")
df = as.data.frame(matrix(nrow = 2* (length(Samples_Healthy) + length(Samples_GBM)),ncol = 13))
colnames(df) = c("Sample","condition","count",feature)

for (s in Samples_Healthy) {
  for (con in c("HealthyPBMC_CD4Tcell","HealthyPBMC_CD8Tcell")) {
    df[i,"Sample"] = s
    df[i,"condition"] = con
    tmp = cell_meta[((cell_meta$orig.ident==s) &(cell_meta$condition==con)),]
    df[i,"count"] = nrow(tmp)
    if(nrow(tmp) == 0){
      df[i,feature] = 0
    }else{
      for (f in feature) {
        df[i,f] = sum(tmp[[f]]>0)
      }
    }
    i = i + 1
  }
}
for (s in Samples_GBM) {
  for (con in c("GBMTissue_CD4Tcell","GBMTissue_CD8Tcell")) {
    df[i,"Sample"] = s
    df[i,"condition"] = con
    tmp = cell_meta[((cell_meta$orig.ident==s) &(cell_meta$condition==con)),]
    df[i,"count"] = nrow(tmp)
    if(nrow(tmp) == 0){
      df[i,feature] = 0
    }else{
      for (f in feature) {
        df[i,f] = sum(tmp[[f]]>0)
      }
    }
    i = i + 1
  }
}
feature_count_df = df
feature_percent_df = df
for (f in feature) {
  feature_percent_df[,f] = feature_percent_df[,f] / feature_percent_df[,"count"]
}



library(ggplot2)
library(ggpubr)

my_colors <- c("GBMTissue_CD4Tcell"= "#c0ebd7",   "GBMTissue_CD8Tcell"= "#177cb0",  
               "HealthyPBMC_CD4Tcell"= "#ffc773", "HealthyPBMC_CD8Tcell"= "#f36838")
my_comparisons <- list( c("GBMTissue_CD8Tcell", "HealthyPBMC_CD8Tcell") )

for (f in feature) {
  p = ggplot(feature_percent_df, aes(x = condition, y = .data[[f]], fill = condition)) +
    
    # 2. Add Jittered Points
    geom_point(position = position_jitterdodge(jitter.width = 0.25, seed = 42),
               shape = 21, size = 3, color = "black", stroke = 0.5, alpha = 0.8) +
    # 1. Add Mean and Error Bars (SEM)
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.5) +
    stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2, color = "black") +
    
    scale_fill_manual(values = my_colors) +
    
    # 3. Fixed P-value logic
    stat_compare_means(
      comparisons = my_comparisons,
      method = "t.test",
      tip.length = 0.02,
      aes(label = after_stat(paste0("p = ", p.format))) # This ensures 'p =' is rendered
    ) +
    
    labs(x = NULL, y = paste0(f, " + cells (%)"), 
         caption = "p value of T.test comparison; % in CD4+/CD8+ cells") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.position = "none",
      plot.caption = element_text(size = 8, face = "italic", color = "grey30", hjust = 0)
    )
  pdf(file = file.path(SaveDir,paste0("grouped_dotplot_",f,".pdf")),width = 3,height = 4)
  print(p)
  dev.off()
}

table(cell_meta$general_annotation)
nrow(cell_meta)


#screen for good predictive markers associated with CD8 ####
## prepare data ####
obj = readRDS(file = file.path(SaveDir,"obj_CD84Tcell.rds"))
cell_meta = readRDS(file = file.path(SaveDir,"cell_meta.rds"))
sample_meta = readRDS(file.path(SaveDir,"sample_meta_scRNAseq_GBM.rds"))
obj@active.assay = "RNA"

table(cell_meta$condition)
cell_meta = cell_meta[cell_meta$condition == "GBMTissue_CD8Tcell",]#only focus on CD8T in GBM
obj = obj[,rownames(cell_meta)]
genes = rownames(obj)
Samples = names(table(cell_meta$orig.ident))

exp_tb = FetchData(object = obj,c("orig.ident",genes))
exp_tb_mean = exp_tb %>%
  group_by(orig.ident) %>%
  summarise(across(all_of(genes),~mean(.x,na.rm = T)))
remove(obj)
remove(exp_tb)
gc()
saveRDS(exp_tb_mean, file = file.path(SaveDir,"exp_tb_mean.rds"))

## clean data ####
cell_meta = readRDS(file = file.path(SaveDir,"cell_meta.rds"))
sample_meta = readRDS(file.path(SaveDir,"sample_meta_scRNAseq_GBM.rds"))
exp_tb_mean = readRDS(file = file.path(SaveDir,"exp_tb_mean.rds"))

exp_tb_mean = as.data.frame(exp_tb_mean)
rownames(exp_tb_mean) = exp_tb_mean$orig.ident
exp_tb_mean$orig.ident = NULL
colnames(sample_meta)
sample_meta = sample_meta[rownames(exp_tb_mean),c("PRM","Age","Gender","OS_status","OS_time","PFS_status","PFS_time")]
sample_meta = sample_meta[!(is.na(sample_meta$OS_status) & is.na(sample_meta$PFS_status)),]
exp_tb_mean = exp_tb_mean[rownames(sample_meta),]
exp_tb_mean = exp_tb_mean[,(colSums(exp_tb_mean !=0)>0)]
genes = colnames(exp_tb_mean)

res_tb = as.data.frame(matrix(nrow = length(genes),ncol = 3))
colnames(res_tb) = c("gene","beta","p")
rownames(res_tb) = genes
res_tb$gene = genes

##Primary, OS ####
tmp_sample_meta = sample_meta
tmp_sample_meta = tmp_sample_meta[(tmp_sample_meta$PRM=="Primary"),c("Age","Gender","OS_status","OS_time")]
tmp_sample_meta = tmp_sample_meta[(!is.na(tmp_sample_meta$OS_status)),]
tmp_exp_tb_mean = exp_tb_mean
tmp_exp_tb_mean = tmp_exp_tb_mean[rownames(tmp_sample_meta),]
tmp_res_tb = res_tb

for (i in 1:length(genes)) {
  if(i%%100 == 0){
    print(i)
  }
  g = genes[i]
  tmp_df = cbind(tmp_sample_meta,tmp_exp_tb_mean[,g])
  colnames(tmp_df)[5] = "exp"
  model = coxph(Surv(OS_time, OS_status) ~ exp + Gender + Age, data=tmp_df)
  res = summary(model)
  tmp_res_tb[g,"beta"] = res[["coefficients"]]["exp","coef"]
  tmp_res_tb[g,"p"] = res[["coefficients"]]["exp","Pr(>|z|)"]
}
write.table(tmp_res_tb,file = file.path(SaveDir,"Cox_OS_PrimaryGBM.csv"),sep = ",")
saveRDS(tmp_res_tb,file = file.path(SaveDir,"Cox_OS_PrimaryGBM.rds"))
res_tb_OS_primary = tmp_res_tb

#Primary, PFS ####
tmp_sample_meta = sample_meta
tmp_sample_meta = tmp_sample_meta[(tmp_sample_meta$PRM=="Primary"),c("Age","Gender","PFS_status","PFS_time")]
tmp_sample_meta = tmp_sample_meta[(!is.na(tmp_sample_meta$PFS_status)),]
tmp_exp_tb_mean = exp_tb_mean
tmp_exp_tb_mean = tmp_exp_tb_mean[rownames(tmp_sample_meta),]
tmp_res_tb = res_tb

for (i in 1:length(genes)) {
  if(i%%100 == 0){
    print(i)
  }
  g = genes[i]
  tmp_df = cbind(tmp_sample_meta,tmp_exp_tb_mean[,g])
  colnames(tmp_df)[5] = "exp"
  model = coxph(Surv(PFS_time, PFS_status) ~ exp + Gender + Age, data=tmp_df)
  res = summary(model)
  tmp_res_tb[g,"beta"] = res[["coefficients"]]["exp","coef"]
  tmp_res_tb[g,"p"] = res[["coefficients"]]["exp","Pr(>|z|)"]
}
write.table(tmp_res_tb,file = file.path(SaveDir,"Cox_PFS_PrimaryGBM.csv"),sep = ",")
saveRDS(tmp_res_tb,file = file.path(SaveDir,"Cox_PFS_PrimaryGBM.rds"))
res_tb_PFS_primary = tmp_res_tb

#Recurrent, OS ####
tmp_sample_meta = sample_meta
tmp_sample_meta = tmp_sample_meta[(tmp_sample_meta$PRM=="Recurrent"),c("Age","Gender","OS_status","OS_time")]
tmp_sample_meta = tmp_sample_meta[(!is.na(tmp_sample_meta$OS_status)),]
tmp_exp_tb_mean = exp_tb_mean
tmp_exp_tb_mean = tmp_exp_tb_mean[rownames(tmp_sample_meta),]
tmp_res_tb = res_tb

for (i in 1:length(genes)) {
  if(i%%100 == 0){
    print(i)
  }
  g = genes[i]
  tmp_df = cbind(tmp_sample_meta,tmp_exp_tb_mean[,g])
  colnames(tmp_df)[5] = "exp"
  model = coxph(Surv(OS_time, OS_status) ~ exp + Gender + Age, data=tmp_df)
  res = summary(model)
  tmp_res_tb[g,"beta"] = res[["coefficients"]]["exp","coef"]
  tmp_res_tb[g,"p"] = res[["coefficients"]]["exp","Pr(>|z|)"]
}
write.table(tmp_res_tb,file = file.path(SaveDir,"Cox_OS_RecurrentGBM.csv"),sep = ",")
saveRDS(tmp_res_tb,file = file.path(SaveDir,"Cox_OS_RecurrentGBM.rds"))
res_tb_OS_Recurrent = tmp_res_tb
#Recurrent, PFS ####
tmp_sample_meta = sample_meta
tmp_sample_meta = tmp_sample_meta[(tmp_sample_meta$PRM=="Recurrent"),c("Age","Gender","PFS_status","PFS_time")]
tmp_sample_meta = tmp_sample_meta[(!is.na(tmp_sample_meta$PFS_status)),]
tmp_exp_tb_mean = exp_tb_mean
tmp_exp_tb_mean = tmp_exp_tb_mean[rownames(tmp_sample_meta),]
tmp_res_tb = res_tb

for (i in 1:length(genes)) {
  if(i%%100 == 0){
    print(i)
  }
  g = genes[i]
  tmp_df = cbind(tmp_sample_meta,tmp_exp_tb_mean[,g])
  colnames(tmp_df)[5] = "exp"
  model = coxph(Surv(PFS_time, PFS_status) ~ exp + Gender + Age, data=tmp_df)
  res = summary(model)
  tmp_res_tb[g,"beta"] = res[["coefficients"]]["exp","coef"]
  tmp_res_tb[g,"p"] = res[["coefficients"]]["exp","Pr(>|z|)"]
}
write.table(tmp_res_tb,file = file.path(SaveDir,"Cox_PFS_RecurrentGBM.csv"),sep = ",")
saveRDS(tmp_res_tb,file = file.path(SaveDir,"Cox_PFS_RecurrentGBM.rds"))
res_tb_PFS_Recurrent = tmp_res_tb

## COLLECTION ####

res_tb_PFS_Recurrent = readRDS(file = file.path(SaveDir,"Cox_PFS_RecurrentGBM.rds"))
res_tb_OS_Recurrent = readRDS(file = file.path(SaveDir,"Cox_OS_RecurrentGBM.rds"))
res_tb_PFS_primary = readRDS(file = file.path(SaveDir,"Cox_PFS_PrimaryGBM.rds"))
res_tb_OS_primary = readRDS(file = file.path(SaveDir,"Cox_OS_PrimaryGBM.rds"))

p_thre = 0.05
g_os_p = res_tb_OS_primary[(!is.na(res_tb_OS_primary$p))&(res_tb_OS_primary$p<p_thre),"gene"]
g_pfs_p = res_tb_PFS_primary[(!is.na(res_tb_PFS_primary$p))&(res_tb_PFS_primary$p<p_thre),"gene"]
g_os_r = res_tb_OS_Recurrent[(!is.na(res_tb_OS_Recurrent$p))&(res_tb_OS_Recurrent$p<p_thre),"gene"]
g_pfs_r = res_tb_PFS_Recurrent[(!is.na(res_tb_PFS_Recurrent$p))&(res_tb_PFS_Recurrent$p<p_thre),"gene"]


g_pfs_r[g_pfs_r %in% g_os_p[g_os_p %in% g_pfs_p]]
# [1] "TMEM128" "IFITM3" 
g_os_r[g_os_r %in% g_os_p[g_os_p %in% g_pfs_p]]
# [1] "TMEM128"
g_os_r[g_os_r %in% g_os_p[g_os_p %in% g_pfs_r]]
# [1] "MFN2"    "GSTM4"   "PHTF1"   "CD8B"    "MAP4"    "RHOA"    "TMCC1"   "TMEM128" "RPL9"    "TNFAIP3" "AKAP12"  "ERMARD" 
# [13] "CCDC120" "GPR107"  "TUBB4B"  "FTH1"    "DUSP16"  "RPL41"   "DEPDC4"  "FNDC3A"  "BAZ1A"   "MNAT1"   "UBE3A"   "CNOT1"  
# [25] "TRPV1"   "CHD3"    "SLFN12L" "USP36"   "TPGS2"   "ATP9B"   "DTD1"   
g_os_r[g_os_r %in% g_pfs_p[g_pfs_p %in% g_pfs_r]]
# [1] "STXBP3"     "NBPF11"     "CTSS"       "HDGF"       "CD244"      "RAB29"      "RRP15"      "B3GNT2"     "METTL21A"  
# [10] "TIGIT"      "H1FX"       "GPR171"     "MECOM"      "TMEM128"    "AIMP1"      "IRF1"       "SMAD5"      "PKD2L2"    
# [19] "RING1"      "EEF1A1"     "TBP"        "CLDN12"     "ARPC1B"     "ALKBH4"     "PIM2"       "TSPYL2"     "FAM120C"   
# [28] "HTATSF1"    "OTUD6B-AS1" "RPL30"      "SET"        "TMEM179B"   "TAGLN"      "MCMBP"      "NACA"       "PIP4K2C"   
# [37] "OAS2"       "EXOSC8"     "ZSCAN2"     "YPEL3"      "EIF1"       "PLTP"       "CEBPB"      "ZBP1"       "SH2D3A"    
# [46] "COX6B1"     "CCDC97"    
g_os_p[g_os_p%in% g_os_r[g_os_r %in% g_pfs_p[g_pfs_p %in% g_pfs_r]]]
# [1] "TMEM128"


g_sum = table(c(g_os_p,g_os_r,g_pfs_p,g_pfs_r))
g_sum = g_sum[order(g_sum,decreasing = T)]

#PIM2 has better specificity
exp_tb_mean[,"SMAD5"]

library(VennDiagram)
x <- list(Set_A = g_os_p, Set_B = g_os_r, Set_C = g_pfs_p, Set_D = g_pfs_r)
p = venn.diagram(
  x = x,
  category.names = c("OS Primary GBM", "OS Recurrent GBM", "PFS Primary GBM", "PFS Recurrent GBM"),
  # filename = file.path(SaveDir,"venn_plot_screen_gene_survival.pdf"),
  filename = NULL,
  output = TRUE,
  # Optional styling
  fill = c("skyblue", "pink1", "mediumorchid", "orange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)
pdf(file = file.path(SaveDir,"venn_plot_screen_gene_survival.pdf"),width = 5,height = 5)
print(p)
dev.off()

plot_forest = function(model,g,title,SaveDir){
  df <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      Variable = recode(term,
                        "exp" = g,
                        "Age" = "Age (years)",
                        "GenderMale" = "Male (vs Female)"
      ),
      HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      p_value = sprintf("%.3f", p.value)
    ) %>%
    select(Variable, HR_CI, p_value, estimate, conf.low, conf.high)
  
  df_text <- df %>%
    mutate(
      HR_CI = sprintf(
        "%s (%s–%s)",
        sci_if_needed(estimate),
        sci_if_needed(conf.low),
        sci_if_needed(conf.high)
      )
    )
  table_text <- rbind(
    c("Variable", "HR (95% CI)", "p"),
    as.matrix(df_text[, c("Variable", "HR_CI", "p_value")])
  )
  
  cap_upper <- 10      
  cap_lower <- 0.001    
  
  df <- df %>%
    mutate(
      estimate_plot  = estimate,
      conf.low_plot  = conf.low,
      conf.high_plot = conf.high,

      # handel Inf and big number
      estimate_plot  = ifelse( estimate_plot > cap_upper, cap_upper, estimate_plot),
      conf.low_plot  = ifelse( conf.low_plot > cap_upper, cap_upper, conf.low_plot),
      conf.high_plot = ifelse( conf.high_plot > cap_upper, cap_upper, conf.high_plot),
      
      # handel -Inf
      estimate_plot  = ifelse(is.infinite(estimate_plot) & estimate_plot < 0, cap_lower, estimate_plot),
      conf.low_plot  = ifelse(is.infinite(conf.low_plot) & conf.low_plot < 0, cap_lower, conf.low_plot),
      conf.high_plot = ifelse(is.infinite(conf.high_plot) & conf.high_plot < 0, cap_lower, conf.high_plot),
      
      # handel 0 
      estimate_plot  = ifelse(estimate_plot == 0, cap_lower, estimate_plot),
      conf.low_plot  = ifelse(conf.low_plot == 0, cap_lower, conf.low_plot),
      conf.high_plot = ifelse(conf.high_plot == 0, cap_lower, conf.high_plot)
    )
  p = forestplot(
    labeltext = table_text,
    mean  = c(NA, df$estimate_plot),
    lower = c(NA, df$conf.low_plot),
    upper = c(NA, df$conf.high_plot),
    
    zero = 1,
    graph.pos = 2,
    xlog = TRUE,
    
    boxsize = 0.20,
    lwd.ci = 1.8,
    ci.vertices = TRUE,
    ci.vertices.height = 0.06,
    
    col = fpColors(
      box  = "#1A4E8A",     
      lines = "#1A4E8A",
      zero = "#333333"    
    ),
    
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.20, fontfamily = "Helvetica"),
      ticks = gpar(cex = 1.10, fontfamily = "Helvetica"),
      xlab  = gpar(cex = 1.25, fontfamily = "Helvetica")
    ),
    
    lineheight = unit(6.5, "mm"),
    colgap = unit(7, "mm"),
    
    xticks = c(0.2, 0.5, 1, 2, 5, 10),
    
    xlab = "Hazard Ratio (log scale)",
    title = title
  )
  pdf(file = paste0(SaveDir,"/",g,"_Forestplot_",title,".pdf"),width = 8,height = 6)
  print(p)
  dev.off()
  dev.off()
  
  return("")
}
plot_forest_group = function(model,g,title,SaveDir){
  df <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      Variable = recode(term,
                        "groupHigh" = paste0(g," High (vs Low)"),
                        "Age" = "Age (years)",
                        "GenderMale" = "Male (vs Female)"
      ),
      HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      p_value = sprintf("%.3f", p.value)
    ) %>%
    select(Variable, HR_CI, p_value, estimate, conf.low, conf.high)
  
  df_text <- df %>%
    mutate(
      HR_CI = sprintf(
        "%s (%s–%s)",
        sci_if_needed(estimate),
        sci_if_needed(conf.low),
        sci_if_needed(conf.high)
      )
    )
  table_text <- rbind(
    c("Variable", "HR (95% CI)", "p"),
    as.matrix(df_text[, c("Variable", "HR_CI", "p_value")])
  )
  
  cap_upper <- 10      
  cap_lower <- 0.001    
  
  df <- df %>%
    mutate(
      estimate_plot  = estimate,
      conf.low_plot  = conf.low,
      conf.high_plot = conf.high,
      
      # handel Inf and big number
      estimate_plot  = ifelse( estimate_plot > cap_upper, cap_upper, estimate_plot),
      conf.low_plot  = ifelse( conf.low_plot > cap_upper, cap_upper, conf.low_plot),
      conf.high_plot = ifelse( conf.high_plot > cap_upper, cap_upper, conf.high_plot),
      
      # handel -Inf
      estimate_plot  = ifelse(is.infinite(estimate_plot) & estimate_plot < 0, cap_lower, estimate_plot),
      conf.low_plot  = ifelse(is.infinite(conf.low_plot) & conf.low_plot < 0, cap_lower, conf.low_plot),
      conf.high_plot = ifelse(is.infinite(conf.high_plot) & conf.high_plot < 0, cap_lower, conf.high_plot),
      
      # handel 0 
      estimate_plot  = ifelse(estimate_plot == 0, cap_lower, estimate_plot),
      conf.low_plot  = ifelse(conf.low_plot == 0, cap_lower, conf.low_plot),
      conf.high_plot = ifelse(conf.high_plot == 0, cap_lower, conf.high_plot)
    )
  p = forestplot(
    labeltext = table_text,
    mean  = c(NA, df$estimate_plot),
    lower = c(NA, df$conf.low_plot),
    upper = c(NA, df$conf.high_plot),
    
    zero = 1,
    graph.pos = 2,
    xlog = TRUE,
    
    boxsize = 0.20,
    lwd.ci = 1.8,
    ci.vertices = TRUE,
    ci.vertices.height = 0.06,
    
    col = fpColors(
      box  = "#1A4E8A",     
      lines = "#1A4E8A",
      zero = "#333333"    
    ),
    
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.20, fontfamily = "Helvetica"),
      ticks = gpar(cex = 1.10, fontfamily = "Helvetica"),
      xlab  = gpar(cex = 1.25, fontfamily = "Helvetica")
    ),
    
    lineheight = unit(6.5, "mm"),
    colgap = unit(7, "mm"),
    
    xticks = c(0.2, 0.5, 1, 2, 5, 10),
    
    xlab = "Hazard Ratio (log scale)",
    title = title
  )
  pdf(file = paste0(SaveDir,"/",g,"_Forestplot_",title,".pdf"),width = 8,height = 6)
  print(p)
  dev.off()
  dev.off()
  
  return("")
}
# dev.off()
dir.create(file.path(SaveDir,"SurvivalCurve_RiskGene"))
g_sum = g_sum[g_sum>2]
# i = 3
sink(file.path(SaveDir,"SurvivalCurve_RiskGene","log.txt"), append=TRUE, split=TRUE)
for (i in c(1:length(g_sum))) {
  g = names(g_sum)[i]
  print(paste0("running for ",g," ..."))
  tmp_exp_tb_mean = exp_tb_mean
  tmp_exp_tb_mean = tmp_exp_tb_mean[rownames(sample_meta),]
  tmp_df = cbind(sample_meta,tmp_exp_tb_mean[,g])
  colnames(tmp_df)[ncol(tmp_df)] = "exp"
  tmp_df$group = "Low"
  tmp_df[tmp_df$exp>median(tmp_df$exp),"group"] = "High"
  tmp_df$group = factor(tmp_df$group,levels = c("Low","High"))
  
  #forest exp
  print(paste0("Cox of OS for ",g, " in Primary GBM: "))
  model = coxph(Surv(OS_time, OS_status) ~ exp + Gender + Age, data=tmp_df[tmp_df$PRM == "Primary",])
  res = summary(model)
  print(res)
  plot_forest(model=model,g = g,
              title = "Cox of OS in Primary GBM",
              SaveDir = file.path(SaveDir,"SurvivalCurve_RiskGene"))
    
  print(paste0("Cox of PFS for ",g, " in Primary GBM: "))
  model = coxph(Surv(PFS_time, PFS_status) ~ exp + Gender + Age, data=tmp_df[tmp_df$PRM == "Primary",])
  res = summary(model)
  print(res)
  plot_forest(model=model,g = g,
              title = "Cox of PFS in Primary GBM",
              SaveDir = file.path(SaveDir,"SurvivalCurve_RiskGene"))
  
  print(paste0("Cox of OS for ",g, " in Recurrent GBM: "))
  model = coxph(Surv(OS_time, OS_status) ~ exp + Gender + Age, data=tmp_df[tmp_df$PRM == "Recurrent",])
  res = summary(model)
  print(res)
  plot_forest(model=model,g = g,
              title = "Cox of OS in Recurrent GBM",
              SaveDir = file.path(SaveDir,"SurvivalCurve_RiskGene"))
  
  print(paste0("Cox of PFS for ",g, " in Recurrent GBM: "))
  model = coxph(Surv(PFS_time, PFS_status) ~ exp + Gender + Age, data=tmp_df[tmp_df$PRM == "Recurrent",])
  res = summary(model)
  print(res)
  plot_forest(model=model,g = g,
              title = "Cox of PFS in Recurrent GBM",
              SaveDir = file.path(SaveDir,"SurvivalCurve_RiskGene"))
  #forest group
  print(paste0("Cox of OS for ",g, " in Primary GBM: "))
  model = coxph(Surv(OS_time, OS_status) ~ group + Gender + Age, data=tmp_df[tmp_df$PRM == "Primary",])
  res = summary(model)
  print(res)
  plot_forest_group(model=model,g = g,
              title = "Cox of OS in Primary GBM_bygroup",
              SaveDir = file.path(SaveDir,"SurvivalCurve_RiskGene"))
  
  print(paste0("Cox of PFS for ",g, " in Primary GBM: "))
  model = coxph(Surv(PFS_time, PFS_status) ~ group + Gender + Age, data=tmp_df[tmp_df$PRM == "Primary",])
  res = summary(model)
  print(res)
  plot_forest_group(model=model,g = g,
              title = "Cox of PFS in Primary GBM_bygroup",
              SaveDir = file.path(SaveDir,"SurvivalCurve_RiskGene"))
  
  print(paste0("Cox of OS for ",g, " in Recurrent GBM: "))
  model = coxph(Surv(OS_time, OS_status) ~ group + Gender + Age, data=tmp_df[tmp_df$PRM == "Recurrent",])
  res = summary(model)
  print(res)
  plot_forest_group(model=model,g = g,
              title = "Cox of OS in Recurrent GBM_bygroup",
              SaveDir = file.path(SaveDir,"SurvivalCurve_RiskGene"))
  
  print(paste0("Cox of PFS for ",g, " in Recurrent GBM: "))
  model = coxph(Surv(PFS_time, PFS_status) ~ group + Gender + Age, data=tmp_df[tmp_df$PRM == "Recurrent",])
  res = summary(model)
  print(res)
  plot_forest_group(model=model,g = g,
              title = "Cox of PFS in Recurrent GBM_bygroup",
              SaveDir = file.path(SaveDir,"SurvivalCurve_RiskGene"))
  
  print(paste0("Plotting survival curve for grouped ",g, "(based on Median) in GBM: "))
  
  data = tmp_df[tmp_df$PRM == "Primary",]
  fit <- survfit(Surv(OS_time, OS_status) ~ group, data = data)
  p_os_primary = ggsurvplot(
    fit, 
    data = data, 
    size = 1,                 # change line size
    palette = c("#3ec1d3", "#ff9a00"),# custom color palettes
    # conf.int = TRUE,          # Add confidence interval
    pval = TRUE,              # Add p-value
    risk.table = TRUE,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c("Low", "High"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw(),      # Change ggplot2 theme
    title = "Overall Survival in Primary GBM",
    subtitle = paste0("group by ",g," tensity")
  )
  
  data = tmp_df[tmp_df$PRM == "Primary",]
  fit <- survfit(Surv(PFS_time, PFS_status) ~ group, data = data)
  p_pfs_primary = ggsurvplot(
    fit, 
    data = data, 
    size = 1,                 # change line size
    palette = c("#3ec1d3", "#ff9a00"),# custom color palettes
    # conf.int = TRUE,          # Add confidence interval
    pval = TRUE,              # Add p-value
    risk.table = TRUE,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c("Low", "High"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw(),      # Change ggplot2 theme
    title = "Progression Free Survival in Primary GBM",
    subtitle = paste0("group by ",g," tensity")
  )
  
  data = tmp_df[tmp_df$PRM == "Recurrent",]
  fit <- survfit(Surv(OS_time, OS_status) ~ group, data = data)
  p_os_recurrent = ggsurvplot(
    fit, 
    data = data, 
    size = 1,                 # change line size
    palette = c("#3ec1d3", "#ff9a00"),# custom color palettes
    # conf.int = TRUE,          # Add confidence interval
    pval = TRUE,              # Add p-value
    risk.table = TRUE,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c("Low", "High"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw(),      # Change ggplot2 theme
    title = "Overall Survival in Recurrent GBM",
    subtitle = paste0("group by ",g," tensity")
  )
  
  data = tmp_df[tmp_df$PRM == "Recurrent",]
  fit <- survfit(Surv(PFS_time, PFS_status) ~ group, data = data)
  p_pfs_recurrent = ggsurvplot(
    fit, 
    data = data, 
    size = 1,                 # change line size
    palette = c("#3ec1d3", "#ff9a00"),# custom color palettes
    # conf.int = TRUE,          # Add confidence interval
    pval = TRUE,              # Add p-value
    risk.table = TRUE,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c("Low", "High"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw(),      # Change ggplot2 theme
    title = "Progression Free Survival in Recurrent GBM",
    subtitle = paste0("group by ",g," tensity")
  )
  
  list_of_plots <- list(p_os_primary, p_pfs_primary, p_os_recurrent, p_pfs_recurrent)
  pdf(file = file.path(SaveDir, "SurvivalCurve_RiskGene", paste0(g,".pdf")), width = 10, height = 10,onefile = F)
  print(arrange_ggsurvplots(list_of_plots, print = F, ncol = 2, nrow = 2))
  dev.off()
}
sink()

#check specificity
TcellDir = "/mnt/Data1/groupjuan/Ming/PanBrainCancer/Results/step5_Analysis_scRNAseq_Tcell/Data"
obj_immune = readRDS(file = file.path(TcellDir,"immune_integrated_harmony_sketch.rds"))
obj_immune = obj_immune[,(obj_immune$dataset == "GSE174554")]
gc()
table(obj_immune$dataset)
obj_nonimmune = readRDS(file = file.path(TcellDir,"nonimmune_integrated_harmony_sketch.rds"))
obj_nonimmune = obj_nonimmune[,(obj_nonimmune$dataset == "GSE174554")]
gc()

all(names(g_sum) %in% rownames(obj_nonimmune))
# [1] TRUE
all(names(g_sum) %in% rownames(obj_immune))
# [1] TRUE
obj_immune = obj_immune[names(g_sum),]
obj_nonimmune = obj_nonimmune[names(g_sum),]
gc()
obj_immune@active.assay = "RNA"
df_immune1 = FetchData(object = obj_immune,c("orig.ident","united_annotation","general_annotation"))
df_immune2 = t(GetAssayData(obj_immune,assay = "RNA",layer ="counts" ))
all(rownames(df_immune1) ==rownames(df_immune2))
df_immune = cbind(df_immune1,df_immune2)

obj_nonimmune@active.assay = "RNA"
df_nonimmune1 = FetchData(object = obj_nonimmune,c("orig.ident","united_annotation","general_annotation"))
df_nonimmune2 = t(GetAssayData(obj_nonimmune,assay = "RNA",layer ="counts" ))
all(rownames(df_nonimmune1) ==rownames(df_nonimmune2))
df_nonimmune = cbind(df_nonimmune1,df_nonimmune2)
all(colnames(df_immune) == colnames(df_nonimmune))

df_all = rbind(df_immune,df_nonimmune)
table(df_all$united_annotation)


df_all[df_all$united_annotation %in% c("CD4_Naive_Tcells","CD4_Treg"),"general_annotation"] = "CD4Tcell"
df_all[df_all$united_annotation %in% c("CD8_Effector_Tcells","CD8_Effector_Tcells_Exhausted",
                               "CD8_Effector_Tcells_Exhausted_Proliferating_IDH_TP53",
                               "CD8_Naive_Tcells"),"general_annotation"] = "CD8Tcell"
df_all[df_all$united_annotation %in% c("DN_NKT","DN_MAIT"),"general_annotation"] = "DN_T"

df_all = df_all[!(df_all$general_annotation %in% c("Bcells","DN_T","ILCs")),]
saveRDS(df_all,file.path(SaveDir,"df_GSE174554_count_sig_genes.rds"))


dir.create(file.path(SaveDir,"Specificity_RiskGene"))
obj <- CreateSeuratObject(counts = t(df_all[,names(g_sum)]), project = "GSE174554", min.cells = 0, min.features = 0)
all(colnames(obj) == rownames(df_all))
# [1] TRUE
obj$general_annotation = df_all$general_annotation 


pdf(file = file.path(SaveDir,"Specificity_RiskGene","GBM_cohort_piechart.pdf"),width = 6,height = 5,onefile=FALSE)
slices <- round(100* table(df_all$general_annotation) / nrow(df_all),digits = 2)
slices <- slices[c("Tumorcells","Glial_cells","Myeloids","Neurons","Stromalcells","CD8Tcell","NKcells","CD4Tcell")]
lbls <- paste0(names(slices),": ",slices,"%")
pie(slices, labels = lbls, main=paste0("Cell proporton in GBM Cohort")) #label shows the label names
dev.off()
df_all$orig.ident


res = df_all %>%
  group_by(general_annotation) %>%
  summarise(across(all_of(names(g_sum)),~mean(.x,na.rm = T)))
res = as.data.frame(res)
rownames(res) = res$general_annotation
res$general_annotation = NULL
res = t(res)

library(pheatmap)
colnames(res)
p = pheatmap(res[,c("CD8Tcell","Glial_cells","Myeloids","Neurons","Stromalcells","Tumorcells")], 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         main = "Specificity of Risk Genes")
pdf(file = file.path(SaveDir,"Specificity_RiskGene","heatmap_specificity_withoutCD4TandNK.pdf"),width = 5,height = 15)
print(p)
dev.off()
dev.off()

p = pheatmap(res, 
             scale = "row", 
             cluster_rows = TRUE, 
             cluster_cols = TRUE,
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D2",
             main = "Specificity of Risk Genes")
pdf(file = file.path(SaveDir,"Specificity_RiskGene","heatmap_specificity_complete.pdf"),width = 5,height = 15)
print(p)
dev.off()
dev.off()
# 
# names(g_sum)[4]
# pdf(file = file.path(SaveDir,"Specificity_RiskGene","test.pdf"),width = 15,height = 15)
# VlnPlot(obj, features = names(g_sum)[4] ,assay = "RNA",layer = "counts",group.by = c("general_annotation"))
# dev.off()
# colnames(obj)
library(scico)
g = names(g_sum)[61]
for (g in names(g_sum)) {
  data <- data.frame(category = colnames(res), 
                     value = res[g,])
  
  p = ggplot(data, aes(x = category, y = value, fill = category)) +
    geom_col() +
    labs(title = g, 
         x = "Cell type", 
         y = "Mean expression") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")+
    scale_fill_scico_d(palette = "batlow")
  pdf(file = file.path(SaveDir,"Specificity_RiskGene",paste0(g,".pdf")),width = 4,height = 4)
  print(p)
  dev.off()
  
  tmp_res = df_all %>% 
    group_by(orig.ident,general_annotation) %>%
    summarise(across(all_of(g),~mean(.x,na.rm = T)))
  tmp_res = as.data.frame(tmp_res)
  colnames(tmp_res) = c("sample","celltype","exp")
  
  p = ggplot(tmp_res, aes(x = celltype, y = exp, color = celltype)) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                 geom = "errorbar", width = 0.2, color = "black") +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
    labs(title = g, 
         x = "", 
         y = "Mean expression",
         caption = "Each dot represents one sample") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none",
          plot.caption = element_text(size = 8, face = "italic", hjust = 0))+
    scale_color_scico_d(palette = "batlow")
  pdf(file = file.path(SaveDir,"Specificity_RiskGene",paste0(g,"_dotplot.pdf")),width = 4,height = 4)
  print(p)
  dev.off()
}

# DEG and Enrichment analysis between SLFN12L+/- CD8 T cells####
dir.create(file.path(SaveDir,"DEG_Enrich_SLFN12L"))
obj = readRDS(file = file.path(SaveDir,"obj_CD84Tcell.rds"))
cell_meta = readRDS(file = file.path(SaveDir,"cell_meta.rds"))
sample_meta = readRDS(file.path(SaveDir,"sample_meta_scRNAseq_GBM.rds"))
obj@active.assay = "RNA"

table(cell_meta$condition)
cell_meta = cell_meta[cell_meta$condition == "GBMTissue_CD8Tcell",]#only focus on CD8T in GBM
obj = obj[,rownames(cell_meta)]
genes = rownames(obj)
Samples = names(table(cell_meta$orig.ident))

exp_tb = FetchData(object = obj,c("orig.ident","SLFN12L"))
rownames(exp_tb) = colnames(obj)
obj$SLFN12L_group = "Low"
obj@meta.data[(rownames(exp_tb)[exp_tb$SLFN12L > median(exp_tb$SLFN12L)]),"SLFN12L_group"] = "High"
DEG_res <- FindMarkers(obj, group.by = "SLFN12L_group",ident.1 = "High", ident.2 = "Low", verbose = FALSE)
write.table(DEG_res,file = file.path(SaveDir,"DEG_Enrich_SLFN12L","DEG_res_group_by_SLFN12L.csv"),sep = ",")
saveRDS(DEG_res,file = file.path(SaveDir,"DEG_Enrich_SLFN12L","DEG_res_group_by_SLFN12L.rds"))
DEG_res[DEG_res$p_val_adj<=0.05,]
# p_val avg_log2FC pct.1 pct.2    p_val_adj
# SLFN12L     0.000000e+00  3.6581672 1.000 0.257 0.000000e+00
# XAF1        6.623978e-12  0.8243187 0.274 0.151 1.023670e-07
# STAT1       1.518454e-10  0.6805728 0.292 0.175 2.346618e-06
# KMT2A       1.049107e-09  0.5450130 0.359 0.237 1.621290e-05
# THEMIS      1.417775e-09  0.3745322 0.703 0.595 2.191029e-05
# SP100       3.010457e-09  0.3910125 0.633 0.505 4.652360e-05
# RAP1GDS1    1.421119e-08  0.3236496 0.478 0.342 2.196198e-04
# PLCXD2      7.441615e-08  1.3966985 0.093 0.037 1.150027e-03
# UTRN        8.349703e-08  0.2426156 0.700 0.585 1.290363e-03
# MIR4435-2HG 1.746849e-07  0.5144115 0.366 0.258 2.699580e-03
# TOX         3.300992e-07  0.3530737 0.635 0.525 5.101354e-03
# TRPM7       4.324222e-07  0.5205651 0.308 0.209 6.682653e-03
# OAS3        4.409112e-07  1.5436530 0.056 0.016 6.813842e-03
# DOCK10      4.935368e-07  0.2338568 0.675 0.543 7.627118e-03
# DDX60       7.353612e-07  0.7160266 0.162 0.090 1.136427e-02
# MX1         7.792696e-07  1.0048207 0.142 0.077 1.204283e-02
# CD96        8.383080e-07  0.2783613 0.689 0.557 1.295521e-02
# CD53        8.976464e-07  0.3390267 0.444 0.330 1.387223e-02
# PYHIN1      1.168857e-06  0.3359081 0.416 0.311 1.806352e-02
# NPIPB5      1.230721e-06  0.3677261 0.213 0.131 1.901956e-02
# SLFN5       1.310526e-06  0.3369177 0.321 0.221 2.025287e-02
# SAMD9L      1.320609e-06  0.8184930 0.212 0.137 2.040869e-02
# DNAJB14     2.588385e-06  0.3893047 0.283 0.191 4.000090e-02
#all upregulated in SLFN12L High group

# library(ggplot2)
# library(ggrepel)
# library(dplyr)
# 
# # 1. Prepare the data
# plot_data <- DEG_res %>%
#   mutate(
#     gene = rownames(.),
#     # Create a significance column for coloring
#     significance = case_when(
#       avg_log2FC >= 1 & p_val_adj <= 0.05 ~ "Up-regulated",
#       avg_log2FC <= -1 & p_val_adj <= 0.05 ~ "Down-regulated",
#       TRUE ~ "Not significant"
#     )
#   )
# 
# # 2. Define the Top 10 genes to label (by significance)
# top_genes <- plot_data %>%
#   filter(p_val_adj < 0.05) %>%
#   slice_min(order_by = p_val_adj, n = 20)
# 
# # 3. Generate the plot
# ggplot(plot_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
#   # Use smaller points with transparency to handle density
#   geom_point(alpha = 0.6, size = 1.2) +
#   
#   # Custom Colors (Journal Standards)
#   scale_color_manual(values = c("Up-regulated" = "#E64B35FF",    # Red (NPG Style)
#                                 "Down-regulated" = "#4DBBD5FF",  # Blue (NPG Style)
#                                 "Not significant" = "grey80")) +
#   
#   # Add Threshold Lines
#   geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", size = 0.3) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.3) +
#   
#   # Smart Labels
#   geom_text_repel(data = top_genes, aes(label = gene),
#                   size = 3.5, color = "black", fontface = "italic",
#                   box.padding = 0.1, max.overlaps = Inf) +
#   
#   # Professional Styling
#   theme_classic() +
#   labs(title = "Differential Expression: Cluster A vs B",
#        subtitle = paste0("Adjusted p-value < 0.05 | |log2FC| > 1"),
#        x = expression(log[2]~Fold~Change),
#        y = expression(-log[10]~Adjusted~P-value)) +
#   theme(
#     text = element_text(family = "Arial", size = 12),
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5),
#     legend.position = "right",
#     legend.title = element_blank(),
#     axis.line = element_line(colour = 'black', size = 0.5)
#   )


upregulated_genes = rownames(DEG_res[DEG_res$p_val_adj<=0.05,])
library(ggplot2)
library(cowplot)
library(ggrepel)
# theme_set(theme_cowplot())

aggregate_obj <- AggregateExpression(obj, group.by = c("SLFN12L_group"), return.seurat = TRUE)
aggregate_obj <- NormalizeData(aggregate_obj)

# Extract the expression data into a data frame
plot_data <- as.data.frame((LayerData(aggregate_obj, assay = "RNA")))

# Create a column to identify which genes to highlight
plot_data$gene <- rownames(plot_data)
plot_data$is_highlight <- ifelse(plot_data$gene %in% upregulated_genes, "Selected", "Unselected")

# Ensure 'Selected' points are plotted LAST so they appear on top
plot_data <- plot_data[order(plot_data$is_highlight, decreasing = FALSE), ]

p <- ggplot(plot_data, aes(x = Low, y = High, color = is_highlight)) +
  geom_point(aes(size = is_highlight, alpha = is_highlight)) +
  # Custom Colors, Sizes, and Transparency
  scale_color_manual(values = c("Selected" = "red", "Unselected" = "grey80")) +
  scale_size_manual(values = c("Selected" = 1.5, "Unselected" = 0.5)) +
  scale_alpha_manual(values = c("Selected" = 1, "Unselected" = 0.5)) +
  # Add Labels
  geom_text_repel(
    data = subset(plot_data, is_highlight == "Selected"),
    aes(label = gene),
    size = 3,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.5
  ) +
  # Limits, Titles, and Theme
  xlim(0, 4) + 
  ylim(0, 4) +
  labs(
    title = "Gene Expression Comparison",
    # subtitle = "SLFN12L Low vs High",
    x = "SLFN12L Low CD8 T Cells",
    y = "SLFN12L High CD8 T Cells",
    caption = "Red dot: Differentially Expressed Genes"
  ) +
  theme_classic() +
  theme(
    legend.position = "none", # Hides the 'Selected/Unselected' legend
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(size = 8, face = "italic", color = "grey30", hjust = 0)
  )

pdf(file = file.path(SaveDir,"DEG_Enrich_SLFN12L","Gene_Expression_Comparison.pdf"),width = 5,height = 5)
print(p)
dev.off()

#Enrichment analysis
#note: linux has issue install clusterProfiler. Use laptop instead.

# Correlation between SLFN12L and other genes in CD8 T cells ####
dir.create(file.path(SaveDir,"Correlation_SLFN12L"))
obj = readRDS(file = file.path(SaveDir,"obj_CD84Tcell.rds"))
cell_meta = readRDS(file = file.path(SaveDir,"cell_meta.rds"))
sample_meta = readRDS(file.path(SaveDir,"sample_meta_scRNAseq_GBM.rds"))
obj@active.assay = "RNA"

table(cell_meta$condition)
cell_meta = cell_meta[cell_meta$condition == "GBMTissue_CD8Tcell",]#only focus on CD8T in GBM
obj = obj[,rownames(cell_meta)]
genes = rownames(obj)
Samples = names(table(cell_meta$orig.ident))

exp_tb = FetchData(object = obj,c("orig.ident",genes))
exp_tb_mean = exp_tb %>%
  group_by(orig.ident) %>%
  summarise(across(all_of(genes),~mean(.x,na.rm = T)))

cor_matrix_data <- exp_tb_mean %>%
  tibble::column_to_rownames("orig.ident") %>%
  as.matrix()
target_vector <- cor_matrix_data[, "SLFN12L"]

stats_list <- apply(cor_matrix_data, 2, function(x) {
  # We use tryCatch in case a gene has zero variance (which breaks cor.test)
  tryCatch({
    test <- cor.test(x, target_vector, method = "spearman")
    return(c(correlation = test$estimate[[1]], p_value = test$p.value))
  }, error = function(e) return(c(correlation = NA, p_value = NA)))
})
#Reformat into a clean table
results_df <- as.data.frame(t(stats_list)) %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is.na(p_value)) %>%
  # Filter out the target gene comparing to itself
  filter(gene != "SLFN12L")

#Apply the Multiple Testing Correction (FDR)
results_df$padj <- p.adjust(results_df$p_value, method = "BH")
results_df <- results_df %>% arrange(desc(correlation))

write.table(results_df,file = file.path(SaveDir,"Correlation_SLFN12L","Correlation_SLFN12L.csv"),sep = ",")
saveRDS(results_df,file = file.path(SaveDir,"Correlation_SLFN12L","Correlation_SLFN12L.rds"))
results_df[results_df$padj<=0.05,]

#scatterplot
sig_genes_df <- results_df %>%
  filter(padj <= 0.05) 
genes_to_plot <- sig_genes_df$gene 
#following genes are correlated under spearman and pearson.
# [1] "SIK3"       "AKAP13"     "PRKCB"      "STAT4"      "TSPAN5"     "PRKCQ"      "FTX"        "KANSL1"    
# [9] "CASP8"      "RBMS1"      "DENND1B"    "STAG1"      "NIPBL"      "FAM13B"     "RPRD2"      "VOPP1"     
# [17] "GPATCH8"    "PPP6R3"     "PIBF1"      "TRAPPC8"    "PHRF1"      "ADK"        "EVL"        "ALG13"     
# [25] "PIGL"       "FCHSD2"     "ST6GAL1"    "RAD51B"     "MCPH1"      "CYFIP2"     "TRIM56"     "AP3B1"     
# [33] "MICU2"      "FAM214A"    "WASF2"      "PITPNB"     "TXLNG"      "CSGALNACT2" "TBL1X"      "ABI2"      
# [41] "CENPC"      "SLC39A9"    "MIER1"      "MAP3K1"     "STX18"      "IKZF3"  
for (g in genes_to_plot) {
  
  # Get padj for the current gene to put in the subtitle
  current_padj <- sig_genes_df$padj[sig_genes_df$gene == g]
  
  # Create the plot
  p <- ggplot(exp_tb_mean, aes_string(x = "SLFN12L", y = g)) +
    geom_point(aes(), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
    # This adds the raw R and p-value automatically
    stat_cor(method = "spearman", label.x.npc = "center") + 
    theme_minimal() +
    labs(
      title = paste("Correlation:", "SLFN12L", "vs", g),
      subtitle = paste0("Adjusted p-value (BH): ", formatC(current_padj, format = "e", digits = 2)),
      x = "SLFN12L (Mean Expression)",
      y = paste(g, "(Mean Expression)"),
      caption = "Each dot represents CD8+ TILs in one sample"
    ) +
    theme(legend.position = "bottom",
          plot.caption = element_text(size = 8, face = "italic", color = "grey30", hjust = 0))
  
  # Print the plot to the PDF device
  pdf(file.path(SaveDir,"Correlation_SLFN12L",paste0("Corplot_SLFN12L_",g,".pdf")), width = 5, height = 5)
  print(p)
  dev.off()
}

#corplot with selected genes
select_genes = c("GZMB","GZMH","GZMK","PRF1","PDCD1","HAVCR2","CTLA4","LAG3","CD274","TIGIT","MKI67","STMN1","TOP2A")
for (g in select_genes) {
  
  # Get padj for the current gene to put in the subtitle
  current_padj <- sig_genes_df$padj[sig_genes_df$gene == g]
  
  # Create the plot
  p <- ggplot(exp_tb_mean, aes_string(x = "SLFN12L", y = g)) +
    geom_point(aes(), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
    # This adds the raw R and p-value automatically
    stat_cor(method = "spearman", label.x.npc = "center") + 
    theme_minimal() +
    labs(
      title = paste("Correlation:", "SLFN12L", "vs", g),
      subtitle = paste0("Adjusted p-value (BH): ", formatC(current_padj, format = "e", digits = 2)),
      x = "SLFN12L (Mean Expression)",
      y = paste(g, "(Mean Expression)"),
      caption = "Each dot represents CD8+ TILs in one sample"
    ) +
    theme(legend.position = "bottom",
          plot.caption = element_text(size = 8, face = "italic", color = "grey30", hjust = 0))
  
  # Print the plot to the PDF device
  pdf(file.path(SaveDir,"Correlation_SLFN12L",paste0("Corplot_SLFN12L_",g,".pdf")), width = 5, height = 5)
  print(p)
  dev.off()
}
#Enrichment analysis
#note: linux has issue install clusterProfiler. Use laptop instead.

# Correlation and DEG between GPR171/EEF1A1/ZBP1/IRF1 and other genes in CD8 T cells ####
target_gene = "ZBP1"#EEF1A1, GPR171, ZBP1, IRF1


obj = readRDS(file = file.path(SaveDir,"obj_CD84Tcell.rds"))
cell_meta = readRDS(file = file.path(SaveDir,"cell_meta.rds"))
sample_meta = readRDS(file.path(SaveDir,"sample_meta_scRNAseq_GBM.rds"))
obj@active.assay = "RNA"

table(cell_meta$condition)
cell_meta = cell_meta[cell_meta$condition == "GBMTissue_CD8Tcell",]#only focus on CD8T in GBM
obj = obj[,rownames(cell_meta)]
genes = rownames(obj)
Samples = names(table(cell_meta$orig.ident))

#Correlation
dir.create(file.path(SaveDir,paste0("Correlation_",target_gene)))
{
exp_tb = FetchData(object = obj,c("orig.ident",genes))
exp_tb_mean = exp_tb %>%
  group_by(orig.ident) %>%
  summarise(across(all_of(genes),~mean(.x,na.rm = T)))
exp_tb_mean[[target_gene]]

cor_matrix_data <- exp_tb_mean %>%
  tibble::column_to_rownames("orig.ident") %>%
  as.matrix()
target_vector <- cor_matrix_data[, target_gene]

stats_list <- apply(cor_matrix_data, 2, function(x) {
  # We use tryCatch in case a gene has zero variance (which breaks cor.test)
  tryCatch({
    test <- cor.test(x, target_vector, method = "spearman")
    return(c(correlation = test$estimate[[1]], p_value = test$p.value))
  }, error = function(e) return(c(correlation = NA, p_value = NA)))
})
#Reformat into a clean table
results_df <- as.data.frame(t(stats_list)) %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is.na(p_value)) %>%
  # Filter out the target gene comparing to itself
  filter(gene != target_gene)

#Apply the Multiple Testing Correction (FDR)
results_df$padj <- p.adjust(results_df$p_value, method = "BH")
results_df <- results_df %>% arrange(desc(correlation))

write.table(results_df,file = file.path(SaveDir,paste0("Correlation_",target_gene),paste0("Correlation_",target_gene,".csv")),sep = ",")
saveRDS(results_df,file = file.path(SaveDir,paste0("Correlation_",target_gene),paste0("Correlation_",target_gene,".rds")))
results_df[results_df$padj<=0.05,]

#scatterplot
sig_genes_df <- results_df %>%
  filter(padj <= 0.05) 
genes_to_plot <- sig_genes_df$gene 
dir.create(file.path(SaveDir,paste0("Correlation_",target_gene),"correlated_genes"))
for (g in genes_to_plot[1:min(50,length(genes_to_plot))]) {
  
  # Get padj for the current gene to put in the subtitle
  current_padj <- sig_genes_df$padj[sig_genes_df$gene == g]
  
  # Create the plot
  p <- ggplot(exp_tb_mean, aes_string(x = target_gene, y = g)) +
    geom_point(aes(), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
    # This adds the raw R and p-value automatically
    stat_cor(method = "spearman", label.x.npc = "center") + 
    theme_minimal() +
    labs(
      title = paste("Correlation:", target_gene, "vs", g),
      subtitle = paste0("Adjusted p-value (BH): ", formatC(current_padj, format = "e", digits = 2)),
      x = paste(target_gene, "(Mean Expression)"),
      y = paste(g, "(Mean Expression)"),
      caption = "Each dot represents CD8+ TILs in one sample"
    ) +
    theme(legend.position = "bottom",
          plot.caption = element_text(size = 8, face = "italic", color = "grey30", hjust = 0))
  
  # Print the plot to the PDF device
  pdf(file.path(SaveDir,paste0("Correlation_",target_gene),"correlated_genes",paste0("Corplot_",target_gene,"_",g,".pdf")), width = 5, height = 5)
  print(p)
  dev.off()
}

#corplot with selected genes
select_genes = c("GZMB","GZMH","GZMK","PRF1","PDCD1","HAVCR2","CTLA4","LAG3","CD274","TIGIT","MKI67","STMN1","TOP2A")
dir.create(file.path(SaveDir,paste0("Correlation_",target_gene),"select_genes"))
for (g in select_genes) {
  
  # Get padj for the current gene to put in the subtitle
  current_padj <- results_df$padj[results_df$gene == g]
  
  # Create the plot
  p <- ggplot(exp_tb_mean, aes_string(x = target_gene, y = g)) +
    geom_point(aes(), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
    # This adds the raw R and p-value automatically
    stat_cor(method = "spearman", label.x.npc = "center") + 
    theme_minimal() +
    labs(
      title = paste("Correlation:", target_gene, "vs", g),
      subtitle = paste0("Adjusted p-value (BH): ", formatC(current_padj, format = "e", digits = 2)),
      x = paste(target_gene, "(Mean Expression)"),
      y = paste(g, "(Mean Expression)"),
      caption = "Each dot represents CD8+ TILs in one sample"
    ) +
    theme(legend.position = "bottom",
          plot.caption = element_text(size = 8, face = "italic", color = "grey30", hjust = 0))
  
  # Print the plot to the PDF device
  pdf(file.path(SaveDir,paste0("Correlation_",target_gene),"select_genes",paste0("Corplot_",target_gene,"_",g,".pdf")), width = 5, height = 5)
  print(p)
  dev.off()
}
#Enrichment analysis
#note: linux has issue install clusterProfiler. Use laptop instead.
}
#DEG 
dir.create(file.path(SaveDir,paste0("DEG_",target_gene)))
{
exp_tb = FetchData(object = obj,c("orig.ident",target_gene))
rownames(exp_tb) = colnames(obj)
obj@meta.data[[paste0(target_gene,"_group")]] = "Low"

obj@meta.data[(rownames(exp_tb)[exp_tb[[target_gene]] > median(exp_tb[[target_gene]])]),paste0(target_gene,"_group")] = "High"
DEG_res <- FindMarkers(obj, group.by = paste0(target_gene,"_group"),ident.1 = "High", ident.2 = "Low", verbose = FALSE)
write.table(DEG_res,file = file.path(SaveDir,paste0("DEG_",target_gene),paste0("DEG_res_group_by_",target_gene,".csv")),sep = ",")
saveRDS(DEG_res,file = file.path(SaveDir,paste0("DEG_",target_gene),paste0("DEG_res_group_by_",target_gene,".rds")))
}
{
upregulated_genes = rownames(DEG_res[DEG_res$p_val_adj<=0.05,])
upregulated_genes = upregulated_genes[1:min(30,length(upregulated_genes))]
library(ggplot2)
library(cowplot)
library(ggrepel)
# theme_set(theme_cowplot())

aggregate_obj <- AggregateExpression(obj, group.by = c(paste0(target_gene,"_group")), return.seurat = TRUE)
aggregate_obj <- NormalizeData(aggregate_obj)

# Extract the expression data into a data frame
plot_data <- as.data.frame((LayerData(aggregate_obj, assay = "RNA")))

# Create a column to identify which genes to highlight
plot_data$gene <- rownames(plot_data)
plot_data$is_highlight <- ifelse(plot_data$gene %in% upregulated_genes, "Selected", "Unselected")

# Ensure 'Selected' points are plotted LAST so they appear on top
plot_data <- plot_data[order(plot_data$is_highlight, decreasing = FALSE), ]

p <- ggplot(plot_data, aes(x = Low, y = High, color = is_highlight)) +
  geom_point(aes(size = is_highlight, alpha = is_highlight)) +
  # Custom Colors, Sizes, and Transparency
  scale_color_manual(values = c("Selected" = "red", "Unselected" = "grey80")) +
  scale_size_manual(values = c("Selected" = 1.5, "Unselected" = 0.5)) +
  scale_alpha_manual(values = c("Selected" = 1, "Unselected" = 0.5)) +
  # Add Labels
  geom_text_repel(
    data = subset(plot_data, is_highlight == "Selected"),
    aes(label = gene),
    size = 3,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.5
  ) +
  # Limits, Titles, and Theme
  xlim(0, 4) + 
  ylim(0, 4) +
  labs(
    title = "Gene Expression Comparison",
    x = paste0(target_gene," Low CD8 T Cells"),
    y = paste0(target_gene," High CD8 T Cells"),
    caption = "Red dot: Differentially Expressed Genes"
  ) +
  theme_classic() +
  theme(
    legend.position = "none", # Hides the 'Selected/Unselected' legend
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(size = 8, face = "italic", color = "grey30", hjust = 0)
  )

pdf(file = file.path(SaveDir,paste0("DEG_",target_gene),"Gene_Expression_Comparison.pdf"),width = 5,height = 5)
print(p)
dev.off()

#Enrichment analysis
#note: linux has issue install clusterProfiler. Use laptop instead.
}

dev.off()

