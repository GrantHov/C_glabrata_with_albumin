library(DESeq2)
library("RColorBrewer")
colors <-  brewer.pal(8, "Set1")
source("https://bioconductor.org/biocLite.R")
library(reshape2)
library(ggplot2)
library(ggrepel)
library(tibble)   ### for inserting coulmn in dataframe
library(EDASeq)
library(VennDiagram)
library(openxlsx)
library(RUVSeq)
library(patchwork)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
#BiocManager::install("org.Hs.eg.db")
##### setting wd
local<-"~/users/tg/current/hhovhannisyan/jena_project/"
remote<-"/nfs/users/tg/hhovhannisyan/jena_project"

setwd(local)



################### GLABRATA DATA ####################


yeast<-NULL
human<-NULL
yeast_names<-NULL
human_names<-NULL
count_file_list_round1<-read.table("albumin/mapping/round_1/count_files_list_round1", header = F)  ### count_files_list_round1 was reordered manualy
count_file_list_round1<-as.vector(count_file_list_round1$V1)

count_file_list_round2<-read.table("albumin/mapping/round_2/count_files_list_round2", header = F)  ### count_files_list_round2 was reordered manualy
count_file_list_round2<-as.vector(count_file_list_round2$V1)

for (f in c(count_file_list_round1, count_file_list_round2))
{
  file<-read.table(paste0("albumin/mapping/",f), sep = "\t", row.names = 1)
  file<-file[-c(1,2,3,4),]
  yeast_col<-as.data.frame(subset(file, !grepl("ENS", row.names(file))))
  yeast_col<-yeast_col[,3,drop=FALSE]
  yeast<-cbind(yeast,yeast_col[,1])
  id_yeast<-strsplit(f,split = "/")[[1]][2]
  id_yeast<-strsplit(id_yeast,split = "_")[[1]][1]
  yeast_names<-c(yeast_names,id_yeast)
  human_data<-as.data.frame(subset(file, grepl("ENS", row.names(file))))
  if (dim(human_data)[1]>0) {
    human_col<-human_data[,3,drop=FALSE]
    human<-cbind(human,human_col[,1])
    id_human<-strsplit(f,split = "/")[[1]][2]
    id_human<-strsplit(id_human,split = "_")[[1]][1]
    human_names<-c(human_names,id_human)
  }
  
}
colnames(yeast)<-yeast_names
colnames(human)<-human_names

rownames(yeast)<-rownames(yeast_col)
rownames(human)<-rownames(human_col)

################### TEST FOR BATCH EFFECT IN CGLAB FOR ROUND1 and ROUND2  ###################################


colData_yeast<-data.frame(time=factor(c(rep("0",3), rep("3",3), rep("24", 3), rep("24c",3), rep("0",3), rep("3",3),rep("24",2),rep("24c",3)),levels=c("0","3","24", "24c")), 
                          batch=c(rep("round1",12), rep("round2",11)), names=yeast_names)


pre_dds_cglab <- DESeqDataSetFromMatrix(yeast, colData=colData_yeast, ~time) ### change yeast to human to see the data for cglab

dds_cglab <- DESeq(pre_dds_cglab)


all_cglab_data_norm<-vst(dds_cglab)
all_cglab_data_norm<-assay(all_cglab_data_norm)

pca <- prcomp(t(all_cglab_data_norm))

percentVar <- pca$sdev^2/sum(pca$sdev^2)

df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],colData_yeast) 

df_plotting$time<-factor(df_plotting$time, levels = c("0","3","24","24c"))
ggplot(df_plotting, aes(PC1,PC2,label=df_plotting$names))+geom_point(aes(colour=time, shape=batch), size=4)+
  geom_text_repel(label=df_plotting$names, size = 3, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + 
  coord_fixed()+ theme_bw()+guides(colour = guide_legend(order = 1),size = guide_legend(order = 2))

# NO BATCH EFFECT OBSERVED


###############  CGLAB  #########################################################


dds_cglab_collapsed <- collapseReplicates(dds_cglab, dds_cglab$names)


cglab_collapsed_counts<-counts(dds_cglab_collapsed, normalized=F)
cglab_collapsed_counts<-cglab_collapsed_counts[,c(1,5:12,2:4)]
#write.table(cglab_collapsed_counts, file="albumin/mapping/results_and_tables/count_data_cglab.txt", sep = "\t", quote = FALSE)

cglab_collapsed_coldata<-data.frame(time=factor(c(rep("0",3), rep("3",3), rep("24", 3), rep("24_control",3)), levels=c("0","3","24", "24_control")))


pre_dds_cglab <- DESeqDataSetFromMatrix(cglab_collapsed_counts, colData=cglab_collapsed_coldata, ~time)

dds_cglab <- DESeq(pre_dds_cglab)

all_cglab_data_norm<-vst(dds_cglab)
all_cglab_data_norm<-assay(all_cglab_data_norm)

pca <- prcomp(t(all_cglab_data_norm))

percentVar <- pca$sdev^2/sum(pca$sdev^2)

df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],cglab_collapsed_coldata) 

df_plotting$time<-factor(df_plotting$time, levels = c("0","3","24","24_control"))
ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=time), size=4)+
  geom_text_repel(label=rownames(df_plotting), size = 3, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + 
  coord_fixed()+ theme_bw()+guides(colour = guide_legend(order = 1),size = guide_legend(order = 2))+
  ggsave("~/Dropbox/Project_albumin/albumin_results_and_tables/fungal_data_PCA.png", width = 6, height = 4, units = "in", dpi = 300)





res_cglab_0_3 <- results(dds_cglab, cooksCutoff=FALSE, contrast = c("time","3","0"))
write.table(res_cglab_0_3, file="albumin/mapping/results_and_tables/res_cglab_0_vs_3_albumin.txt", sep = "\t", quote = FALSE)

res_cglab_0_24 <- results(dds_cglab, cooksCutoff=FALSE, contrast = c("time","24","0"))
write.table(res_cglab_0_24, file="albumin/mapping/results_and_tables/res_cglab_0_vs_24_albumin.txt", sep = "\t", quote = FALSE)

res_cglab_0_24c <- results(dds_cglab, cooksCutoff=FALSE, contrast = c("time","24_control","0"))
write.table(res_cglab_0_24c, file="albumin/mapping/results_and_tables/res_cglab_0_vs_24c_albumin.txt", sep = "\t", quote = FALSE)


res_cglab_24c_24 <- results(dds_cglab, cooksCutoff=FALSE, contrast = c("time","24","24_control"))
write.table(res_cglab_24c_24, file="albumin/mapping/results_and_tables/res_cglab_24c_vs_24.txt", sep = "\t", quote = FALSE)


#######################################  TEST FOR BATCH EFFECT IN HUMAN FOR ROUND1 and ROUND2   #########################

colData_human<-data.frame(time=factor(c(rep("0",3), rep("3",3), rep("24", 3), rep("24c",3), rep("0",3), rep("3",3),rep("24",2),rep("24c",3)),levels=c("0","3","24", "24c")), 
                          batch=c(rep("round1",12), rep("round2",11)), names=human_names)


pre_dds_human <- DESeqDataSetFromMatrix(human, colData=colData_human, ~time) 

dds_human <- DESeq(pre_dds_human)


all_human_data_norm<-vst(dds_human)
all_data_data_norm<-assay(all_human_data_norm)

pca <- prcomp(t(all_data_data_norm))

percentVar <- pca$sdev^2/sum(pca$sdev^2)

df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],colData_human) 

df_plotting$time<-factor(df_plotting$time, levels = c("0","3","24","24c"))
ggplot(df_plotting, aes(PC1,PC2,label=df_plotting$names))+geom_point(aes(colour=time, shape=batch), size=4)+
  geom_text_repel(label=df_plotting$names, size = 3, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + 
  coord_fixed()+ theme_bw()+guides(colour = guide_legend(order = 1),size = guide_legend(order = 2))

#############   HUMAN  ################################################


dds_human_collapsed <- collapseReplicates(dds_human, dds_human$names)


human_collapsed_counts<-counts(dds_human_collapsed, normalized=F)
human_collapsed_counts<-human_collapsed_counts[,c(1,5:12,2:4)]
#write.table(human_collapsed_counts, file="albumin/mapping/results_and_tables/count_data_human.txt", sep = "\t", quote = FALSE)

human_collapsed_coldata<-data.frame(time=factor(c(rep("0",3), rep("3",3), rep("24", 3), rep("24_control",3)), levels=c("0","3","24", "24_control")))


pre_dds_human <- DESeqDataSetFromMatrix(human_collapsed_counts, colData=human_collapsed_coldata, ~time)

dds_human <- DESeq(pre_dds_human)


all_human_data_norm<-vst(dds_human)
all_data_data_norm<-assay(all_human_data_norm)

pca <- prcomp(t(all_data_data_norm))

percentVar <- pca$sdev^2/sum(pca$sdev^2)

df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],human_collapsed_coldata) 

df_plotting$time<-factor(df_plotting$time, levels = c("0","3","24","24_control"))
ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=time), size=4)+
  geom_text_repel(label=rownames(df_plotting), size = 3, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + 
  coord_fixed()+ theme_bw()+guides(colour = guide_legend(order = 1),size = guide_legend(order = 2))+
  ggsave("~/Dropbox/Project_albumin/albumin_results_and_tables/human_data_PCA.png", width = 6, height = 4, units = "in", dpi = 300)




res_human_0_3 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","3","0"))
write.table(res_human_0_3, file="albumin/mapping/results_and_tables/res_human_0_vs_3_albumin.txt", sep = "\t", quote = FALSE)

res_human_0_24 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","24","0"))
write.table(res_human_0_24, file="albumin/mapping/results_and_tables/res_human_0_vs_24_albumin.txt", sep = "\t", quote = FALSE)

res_human_0_24c <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","24_control","0"))
write.table(res_human_0_24c, file="albumin/mapping/results_and_tables/res_human_0_vs_24c_albumin.txt", sep = "\t", quote = FALSE)


res_human_24c_24 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","24","24_control"))
write.table(res_human_24c_24, file="albumin/mapping/results_and_tables/res_human_24c_vs_24.txt", sep = "\t", quote = FALSE)


as.data.frame(res_human_0_3[c("ENSG00000198886","ENSG00000198888","ENSG00000198763",
                              "ENSG00000198840","ENSG00000212907","ENSG00000198786","ENSG00000198695","ENSG00000198727",
                              "ENSG00000198804","ENSG00000198712","ENSG00000198938","ENSG00000198899","ENSG00000228253"),])



#####   CREATING EXCELL DOCUMENT   ########


wb <- createWorkbook("cglab_albumin_data")

header<-c("gene", "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
fungal_data_files<-c("res_cglab_0_3","res_cglab_0_24","res_cglab_0_24c","res_cglab_24c_24")

human_data_files<-c("res_human_0_3","res_human_0_24","res_human_0_24c","res_human_24c_24")

for (file in fungal_data_files){
  assign("fung",as.data.frame(eval(parse(text=paste0(file)))))
  fung <- cbind("genes"=rownames(fung), data.frame(fung, row.names=NULL))
  addWorksheet(wb, file)
  writeData(wb, file, fung,keepNA = T)
}

for (file in human_data_files){
  assign("hum",as.data.frame(eval(parse(text=paste0(file)))))
  hum<-hum[complete.cases(hum),]
  hum <- cbind("genes"=rownames(hum), data.frame(hum, row.names=NULL))
  addWorksheet(wb, file)
  writeData(wb, file, hum, keepNA = T)
}

saveWorkbook(wb, "albumin/mapping/results_and_tables/cglab_albumin_data.xlsx", overwrite = TRUE)


### GO TERMS CGLAB



### Finding MF,BP and CC terms
all_go_tems<-as.data.frame(GOTERM)

all_go_tems<-all_go_tems[!(all_go_tems$go_id %in% c("GO:0003674","GO:0008150","GO:0005575")),]


MF_terms<-all_go_tems$go_id[all_go_tems$Ontology=="MF"]
BP_terms<-all_go_tems$go_id[all_go_tems$Ontology=="BP"]
CC_terms<-all_go_tems$go_id[all_go_tems$Ontology=="CC"]


### Association of GOID with description
goterms <- Term(GOTERM)
a<-as.data.frame(goterms)
go_names<-cbind(row.names(a),a)


### CALB GOIDS
cglab_go<-read.table("~/Dropbox (Personal)/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab_go.txt")


MF_universe<-cglab_go[cglab_go$V1 %in% MF_terms,]
BP_universe<-cglab_go[cglab_go$V1 %in% BP_terms,]
CC_universe<-cglab_go[cglab_go$V1 %in% CC_terms,]



fungal_goterms<-function(file,name,ont){
  ego<-enricher(file, pvalueCutoff = 0.05, 
                pAdjustMethod = "BH", universe = eval(parse(text=paste0("as.character(",ont,"_universe$V2)"))),minGSSize = 2, 
                maxGSSize = 10000, TERM2GENE = eval(parse(text=paste0(ont,"_universe"))),TERM2NAME = go_names)
  
  if (!is.null(ego)){
    p<-dotplot(ego,showCategory = 30)
    p+ggplot2::ggsave(sprintf("~/Dropbox/Project_albumin/albumin_results_and_tables/GO_terms/fungi/%s_%s.png",name,ont),
                      units="in", width=10, heigh=7, dpi=600)
    
    write.table(ego,sprintf("~/Dropbox/Project_albumin/albumin_results_and_tables/GO_terms/fungi/%s_%s.txt",name,ont), sep="\t")
  }
  
}


for (comp in fungal_data_files){
  for (ontology in c("BP","MF","CC")){
    data<- eval(parse(text = c(comp)))
    
    assign(paste0(comp,"_up"),row.names(data[!is.na(data$padj) & data$padj<0.01 & data$log2FoldChange>1.5,]))
    assign(paste0(comp,"_down"),row.names(data[!is.na(data$padj) & data$padj<0.01 & data$log2FoldChange< -1.5,]))
    
    fungal_goterms(eval(parse(text = paste0(comp,"_up"))),paste0(comp,"_up"),ontology)
    fungal_goterms(eval(parse(text = paste0(comp,"_down"))),paste0(comp,"_down"),ontology)
  }
}



############ HUMAN GO TERMS  #####################

human_goterms<-function(file,name,ontology){
  ego<-enrichGO(file, 
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                pvalueCutoff = 0.05, 
                pAdjustMethod = "BH",
                ont = ontology)
  if (!is.null(ego)){
    p<-dotplot(ego,showCategory = 30)
    p+ggplot2::ggsave(sprintf("~/Dropbox/Project_albumin/albumin_results_and_tables/GO_terms/human/%s_%s.png",name,ontology),
                      units="in", width=10, heigh=7, dpi=600)
    write.table(ego,sprintf("~/Dropbox/Project_albumin/albumin_results_and_tables/GO_terms/human/%s_%s.txt",name,ontology), sep="\t")
    
  }
}


for (comp in human_data_files){
  for (ontology in c("BP","MF","CC")){
    data <- eval(parse(text = c(comp)))
    
    assign(paste0(comp,"_up"),row.names(data[!is.na(data$padj) & data$padj<0.01 & data$log2FoldChange>1.5,]))
    assign(paste0(comp,"_down"),row.names(data[!is.na(data$padj) & data$padj<0.01 & data$log2FoldChange< -1.5,]))
    
    human_goterms(eval(parse(text = paste0(comp,"_up"))),paste0(comp,"_up"),ontology)
    human_goterms(eval(parse(text = paste0(comp,"_down"))),paste0(comp,"_down"),ontology)
  }
}


#### GO comparisons with bigRNASeq data


#### Human


## bigRNASeq                     
bigRNASeq_0_3 <-read.table("~/Dropbox (Personal)/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab/res_human_cglab_0_vs_3.txt")
bigRNASeq_0_24 <-read.table("~/Dropbox (Personal)/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab/res_human_cglab_0_vs_24.txt")

bigRNASeq_up_0_3 <- bigRNASeq_0_3[!is.na(bigRNASeq_0_3$padj) & bigRNASeq_0_3$log2FoldChange>1.5 & bigRNASeq_0_3$padj<0.01,]
bigRNASeq_up_0_24 <- bigRNASeq_0_24[!is.na(bigRNASeq_0_24$padj) & bigRNASeq_0_24$log2FoldChange>1.5 & bigRNASeq_0_24$padj<0.01,]

bigRNASeq_down_0_3 <- bigRNASeq_0_3[!is.na(bigRNASeq_0_3$padj) & bigRNASeq_0_3$log2FoldChang < -1.5 & bigRNASeq_0_3$padj<0.01,]
bigRNASeq_down_0_24 <- bigRNASeq_0_24[!is.na(bigRNASeq_0_24$padj) & bigRNASeq_0_24$log2FoldChange < -1.5 & bigRNASeq_0_24$padj<0.01,]

## current  
current_up_0_3 <- res_human_0_3[!is.na(res_human_0_3$padj) & res_human_0_3$log2FoldChange>1.5 & res_human_0_3$padj<0.01,]
current_up_0_24 <- res_human_0_24[!is.na(res_human_0_24$padj) & res_human_0_24$log2FoldChange>1.5 & res_human_0_24$padj<0.01,]

current_down_0_3 <- res_human_0_3[!is.na(res_human_0_3$padj) & res_human_0_3$log2FoldChange < -1.5 & res_human_0_3$padj<0.01,]
current_down_0_24 <- res_human_0_24[!is.na(res_human_0_24$padj) & res_human_0_24$log2FoldChange < -1.5 & res_human_0_24$padj<0.01,]

grouped_up <- rbind.data.frame(cbind.data.frame("gene"=rownames(bigRNASeq_up_0_3),"l2fc"=bigRNASeq_up_0_3$log2FoldChange,"spp"="cglab","time"="3hpi"),
                               cbind.data.frame("gene"=rownames(bigRNASeq_up_0_24),"l2fc"=bigRNASeq_up_0_24$log2FoldChange,"spp"="cglab","time"="24hpi"),
                               cbind.data.frame("gene"=rownames(current_up_0_3),"l2fc"=current_up_0_3$log2FoldChange,"spp"="cglab_albumin","time"="3hpi"),
                               cbind.data.frame("gene"=rownames(current_up_0_24),"l2fc"=current_up_0_24$log2FoldChange,"spp"="cglab_albumin","time"="24hpi"))

grouped_down <- rbind.data.frame(cbind.data.frame("gene"=rownames(bigRNASeq_down_0_3),"l2fc"=bigRNASeq_down_0_3$log2FoldChange,"spp"="cglab","time"="3hpi"),
                                 cbind.data.frame("gene"=rownames(bigRNASeq_down_0_24),"l2fc"=bigRNASeq_down_0_24$log2FoldChange,"spp"="cglab","time"="24hpi"),
                                 cbind.data.frame("gene"=rownames(current_down_0_3),"l2fc"=current_down_0_3$log2FoldChange,"spp"="cglab_albumin","time"="3hpi"),
                                 cbind.data.frame("gene"=rownames(current_down_0_24),"l2fc"=current_down_0_24$log2FoldChange,"spp"="cglab_albumin","time"="24hpi"))



formula_res_up <- compareCluster(gene~spp+time,OrgDb = org.Hs.eg.db, data=grouped_up, ont = "BP",keyType="ENSEMBL", fun="enrichGO")
clusterProfiler::simplify(formula_res_up,
                          cutoff = 0.7,
                          by = "p.adjust",
                          select_fun = min,
                          measure = "Wang",
                          semData = NULL)


formula_res_up_to_plot<-formula_res_up
df<-as.data.frame(formula_res_up_to_plot)
fake_1 <- df[1,]
fake_2 <- df[1,]

fake_1$Cluster <- "cglab_albumin.3hpi"
fake_1$Description <- ""
fake_1$GeneRatio <- 0
fake_2$Cluster <- "cglab.24hpi"
fake_2$Description <- ""
fake_2$time <- "24hpi"
fake_2$GeneRatio <- 0
df <- rbind.data.frame(df,fake_1,fake_2)

#df$GeneRatio<-1
df$time<-factor(df$time, levels = c("3hpi","24hpi"))
df$Cluster<-factor(df$Cluster, levels = c("cglab_albumin.3hpi","cglab.3hpi","cglab_albumin.24hpi","cglab.24hpi"))
formula_res_up_to_plot@compareClusterResult<-df
pup<-dotplot(formula_res_up_to_plot,showCategory=20)+ ggplot2::facet_grid(~time, scales = "free")


formula_res_down <- compareCluster(gene~spp+time,OrgDb = org.Hs.eg.db, data=grouped_down, ont = "BP",keyType="ENSEMBL", fun="enrichGO")

formula_res_down_to_plot<-formula_res_down
df<-as.data.frame(formula_res_down_to_plot)
#df$GeneRatio<-1
df$time<-factor(df$time, levels = c("3hpi","24hpi"))
df$Cluster<-factor(df$Cluster, levels = c("cglab_albumin.3hpi","cglab.3hpi","cglab_albumin.24hpi","cglab.24hpi"))
formula_res_down_to_plot@compareClusterResult<-df
pdown<-dotplot(formula_res_down_to_plot,showCategory=20)+ ggplot2::facet_grid(~time, scales = "free")


png("~/Dropbox/Project_albumin/albumin_results_and_tables/GO_human_up.png", width = 12,height = 14,units = "in", res = 300)
pup
dev.off()

png("~/Dropbox/Project_albumin/albumin_results_and_tables/GO_human_down.png", width = 12,height = 14,units = "in", res = 300)
pdown
dev.off()


## fungus

current_0_3_up <-read.table("~/Dropbox (Personal)/Project_albumin/albumin_results_and_tables/GO_terms/fungi/res_cglab_0_3_up_BP.txt")
current_0_24_up <-read.table("~/Dropbox (Personal)/Project_albumin/albumin_results_and_tables/GO_terms/fungi/res_cglab_0_24_up_BP.txt")

current_0_3_down <-read.table("~/Dropbox (Personal)/Project_albumin/albumin_results_and_tables/GO_terms/fungi/res_cglab_0_3_down_BP.txt")
current_0_24_down <-read.table("~/Dropbox (Personal)/Project_albumin/albumin_results_and_tables/GO_terms/fungi/res_cglab_0_24_down_BP.txt")



big_0_3_up <- read.table("~/Dropbox (Personal)/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab/fungi/res_cglab_0_3_up_BP.txt")
big_0_24_up <- read.table("~/Dropbox (Personal)/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab/fungi/res_cglab_0_24_up_BP.txt")

big_0_3_down <- read.table("~/Dropbox (Personal)/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab/fungi/res_cglab_0_3_down_BP.txt")
big_0_24_down <- read.table("~/Dropbox (Personal)/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab/fungi/res_cglab_0_24_down_BP.txt")


up_to_plot <- rbind(cbind.data.frame("spp"="cglab_albumin","time"="3hpi",current_0_3_up),
                    cbind.data.frame("spp"="cglab_albumin","time"="24hpi",current_0_24_up),
                    cbind.data.frame("spp"="cglab","time"="3hpi",big_0_3_up),
                    cbind.data.frame("spp"="cglab","time"="24hpi",big_0_24_up))
#up_to_plot$GeneRatio <- 1
up_to_plot<-cbind.data.frame("Cluster"=paste0(as.character(up_to_plot$spp),".",as.character(up_to_plot$time)),up_to_plot)

down_to_plot <- rbind(cbind.data.frame("spp"="cglab_albumin","time"="3hpi",current_0_3_down),
                      cbind.data.frame("spp"="cglab_albumin","time"="24hpi",current_0_24_down),
                      cbind.data.frame("spp"="cglab","time"="3hpi",big_0_3_down),
                      cbind.data.frame("spp"="cglab","time"="24hpi",big_0_24_down))
#down_to_plot$GeneRatio <- 1
down_to_plot<-cbind.data.frame("Cluster"=paste0(as.character(down_to_plot$spp),".",as.character(down_to_plot$time)),down_to_plot)



### hacking clusterProfile
formula_res_up_to_plot_fungal<-formula_res_up

formula_res_up_to_plot_fungal@compareClusterResult<-up_to_plot

df<-as.data.frame(formula_res_up_to_plot_fungal)
#df$GeneRatio<-1
df$time<-factor(df$time, levels = c("3hpi","24hpi"))
formula_res_up_to_plot_fungal@compareClusterResult<-df
pup_fung<-dotplot(formula_res_up_to_plot_fungal,showCategory=20)+ ggplot2::facet_grid(~time,scales = "free")

### only albumin

only_alb<-formula_res_up_to_plot_fungal@compareClusterResult[formula_res_up_to_plot_fungal@compareClusterResult$spp=="cglab_albumin",]
only_alb <- only_alb[-c(13,26,28,29,35),]

formula_res_up_to_plot_fungal@compareClusterResult <- only_alb
enrich
bp2 <- clusterProfiler::simplify(formula_res_up_to_plot_fungal,   cutoff = 0.7,
                                 by = "p.adjust",
                                 select_fun = min,
                                 measure = "Wang",
                                 semData = NULL)


class(formula_res_up_to_plot_fungal)
pup_fung_alb<-dotplot(formula_res_up_to_plot_fungal,showCategory=20)+ ggplot2::facet_grid(~time,scales = "free")
pup_fung_alb+ggsave("~/Dropbox (Personal)/Project_albumin/albumin_results_and_tables/GO_fungal_up_only_albumin_no_redund.pdf",device = "pdf",dpi = 600,width = 12,height = 10)


formula_res_down_to_plot_fungal<-formula_res_down
formula_res_down_to_plot_fungal@compareClusterResult<-down_to_plot
df<-as.data.frame(formula_res_down_to_plot_fungal)
#df$GeneRatio<-1
df$time<-factor(df$time, levels = c("3hpi","24hpi"))
formula_res_down_to_plot_fungal@compareClusterResult<-df
pdown_fung<-dotplot(formula_res_down_to_plot_fungal,showCategory=20)+ ggplot2::facet_grid(~time,scales = "free")


png("~/Dropbox/Project_albumin/albumin_results_and_tables/GO_fungal_up.png", width = 12,height = 14,units = "in", res = 300)
pup_fung
dev.off()

png("~/Dropbox/Project_albumin/albumin_results_and_tables/GO_fungal_down.png", width = 20,height = 12,units = "in", res = 300)
pdown_fung
dev.off()



##### human 24 BP up

to_plot<-res_human_0_24[!is.na(res_human_0_24$padj) & res_human_0_24$padj<0.01 & res_human_0_24$log2FoldChange>1.5,]

ego<-enrichGO(rownames(to_plot), 
              OrgDb = org.Hs.eg.db,
              keyType = "ENSEMBL",
              pvalueCutoff = 0.05, 
              pAdjustMethod = "BH",
              ont = "BP")


a<- as.data.frame(ego@result)
ego@result  <- ego@result[-c(1,2,3,7,9,10,11),]  


p<-dotplot(ego,showCategory = 30)
p+ggplot2::ggsave("~/Dropbox (Personal)/Project_albumin/albumin_results_and_tables/human_24UP_BP.pdf",units="in", width=10, heigh=7, dpi=600)

#### PCAs

### human PCA all samples


human_with_e_and_D_calb<-read.table("~/Dropbox/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/calb_human_data.txt",check.names = FALSE)
human_with_e_and_D_cglab<-read.table("~/Dropbox/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab_human_data.txt",check.names = FALSE)
human_with_e_and_D_cpar<-read.table("~/Dropbox/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cpar_human_data.txt",check.names = FALSE)
human_with_e_and_D_ctrop<-read.table("~/Dropbox/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/ctrop_human_data.txt",check.names = FALSE)

all_human_data<-cbind.data.frame(human_with_e_and_D_calb,human_with_e_and_D_cglab,human_with_e_and_D_cpar,human_with_e_and_D_ctrop)

#### collapsing technical replicates
all_human_data<-add_column(all_human_data, "3+reseq" = all_human_data$`3`+all_human_data$`3reseq`, .after = "2")
all_human_data$`3reseq`<-NULL
all_human_data$`3`<-NULL

all_human_data<-add_column(all_human_data, "24+reseq" = all_human_data$`24`+all_human_data$`24reseq`, .after = "D6")
all_human_data$`24reseq`<-NULL
all_human_data$`24`<-NULL

all_human_data<-add_column(all_human_data, "36+reseq" = all_human_data$`36`+all_human_data$`36reseq`, .after = "35")
all_human_data$`36reseq`<-NULL
all_human_data$`36`<-NULL

all_human_data<-add_column(all_human_data, "37+reseq" = all_human_data$`37`+all_human_data$`37reseq`, .after = "D9")
all_human_data$`37reseq`<-NULL
all_human_data$`37`<-NULL

all_human_data<-add_column(all_human_data, "52+2reseq" = all_human_data$`52`+all_human_data$`52reseq`+all_human_data$`52reseq2`, .after = "51")
all_human_data$`52reseq`<-NULL
all_human_data$`52reseq2`<-NULL
all_human_data$`52`<-NULL

all_human_data<-add_column(all_human_data, "53+reseq" = all_human_data$`53`+all_human_data$`53reseq`, .after = "52+2reseq")
all_human_data$`53reseq`<-NULL
all_human_data$`53`<-NULL

all_human_data <- cbind.data.frame(all_human_data,human_collapsed_counts)


### all data, no batch correction

all_human_data_to_plot<-all_human_data[,-c(3:5,9:14,27:32,39:44,51:56)]
species<-c(rep("control",2),rep("calb",10), rep("control",2),rep("cglab",6),rep("cpar",6),rep("ctrop",6), rep("cglab_albumin",9),rep("control_albumin",3))
time_point<-c(c("0_control","0_control","3","3","3","24","24","24","24","24_ece1","24_ece1","24_ece1","24_control","24_control"),     ## calb
              c("3","3","3","24","24","24"),   ### cglab
              c("3","3","3","24","24","24"),   ### cpar
              c("3","3","3","24","24","24"),   ### ctrop
              c("0_control","0_control","0_control","3","3","3","24","24","24","24_control","24_control","24_control"))  ### cglab_albumin

colData<-data.frame(species=species, time_point=time_point)
pre_dds <- DESeqDataSetFromMatrix(all_human_data_to_plot, colData=colData,~time_point)
dds_human <- DESeq(pre_dds)
all_data_norm<-vst(dds_human)
vst_data<-assay(all_data_norm)
pca <- prcomp(t(vst_data))

percentVar <- pca$sdev^2/sum(pca$sdev^2)

df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2]*-1,colData) ### do PC2 = pca$x[, 2]*-1 to invert values

myColors <- brewer.pal(8,"Dark2")
names(myColors) <- c("0_control","1.5","3","24","24_control","24_ece1")
colScale <- scale_colour_manual(name = "Time point (h)",values = myColors)


df_plotting$time_point <- factor(df_plotting$time_point, levels = c("0_control","3","24","24_ece1","24_control"))
df_plotting$species <- factor(df_plotting$species, levels = c("calb","cglab","cpar","ctrop","control","cglab_albumin","control_albumin"))

ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=time_point, shape=species), size=6)+
  geom_text_repel(label=as.character(rownames(df_plotting)), size = 4, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
  coord_fixed()+ theme_bw()+colScale+
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2))+
  labs(shape="Species", colour="Time point")+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size=16),
        legend.text = element_text(size=14),legend.title = element_text(size=14))+
  scale_shape_manual(values = c(4,8,15,16,17,18,21))+
  ggsave("~/Dropbox/Project_albumin/albumin_results_and_tables/human_pca_all_data_no_batch_correction.png", units="in", width=10, heigh=7, dpi=600)




##### cglab data only, no batch corrections

all_human_data_cglab<-all_human_data[,-c(3:21,27:29,36:59)]

all_human_data_cglab<- all_human_data_cglab[,-c(8:10)]

species<-c(rep("no albumin",10) ,rep("with albumin",12))
time_point<-c(c("0_control","0_control","24_control","24_control"),     ## calb
              c("3","3","3","24","24","24"),   ### cglab
              c("0_control","0_control","0_control","3","3","3","24","24","24","24_control","24_control","24_control"))  ### cglab_albumin

colData<-data.frame(species=species, time_point=time_point)
pre_dds <- DESeqDataSetFromMatrix(all_human_data_cglab, colData=colData,~time_point)
dds_human <- DESeq(pre_dds)
all_data_norm<-vst(dds_human)
vst_data<-assay(all_data_norm)
pca <- prcomp(t(vst_data))

percentVar <- pca$sdev^2/sum(pca$sdev^2)

df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],colData) ### do PC2 = pca$x[, 2]*-1 to invert values

myColors <- brewer.pal(8,"Dark2")
names(myColors) <- c("0_control","1.5","3","24","24_control","24_ece1")
colScale <- scale_colour_manual(name = "Time point (h)",values = myColors)


df_plotting$time_point <- factor(df_plotting$time_point, levels = c("0_control","3","24","24_ece1","24_control"))


ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=time_point, shape=species), size=6)+
  geom_text_repel(label=as.character(rownames(df_plotting)), size = 4, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
  coord_fixed()+ theme_bw()+colScale+
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2))+
  labs(shape="Species", colour="Time point")+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size=16),
        legend.text = element_text(size=14),legend.title = element_text(size=14))+
  ggsave("~/Dropbox/Project_albumin/albumin_results_and_tables/human_pca_cglab_data_no_batch_correction.png", units="in", width=10, heigh=7, dpi=600)

########batch corrections

##### NON-DE GENES DURING INFECTIONS (ideally we also need non-DE due to albumin addition)
for (specie in c("calb","cglab","cpar","ctrop")){
  data<-read.table(paste0("~/Dropbox/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/",specie,"/res_lrt_human_",specie,".txt"))
  assign(paste0("human_",specie,"_rlt_no_DE"),data[!is.na(data$padj) & data$padj>0.05 & data$baseMean>10,]) 
}

overlap_for_batch<-as.data.frame(Reduce(intersect,list(rownames(human_calb_rlt_no_DE),
                                                       rownames(human_cglab_rlt_no_DE),
                                                       rownames(human_cpar_rlt_no_DE),
                                                       rownames(human_ctrop_rlt_no_DE))))
colnames(overlap_for_batch)<-"genes"


deseq2_human <- function(cond) {
  myColors <- brewer.pal(8,"Dark2")
  names(myColors) <- c("0_control","1.5","3","24","24_control","24_ece1")
  colScale <- scale_colour_manual(name = "Time point (h)",values = myColors)
  
  
  
  if (cond == c("with dead and ece1")){
    ## uncomment below to exclude 1.5 h samples
    all_human_data<-all_human_data[,-c(3:5,9:14,27:32,39:44,51:56)]
    species<-c(rep("control",2),rep("calb",10), rep("control",2),rep("cglab",6),rep("cpar",6),rep("ctrop",6), rep("cglab_albumin",9),rep("control_albumin",3))
    time_point<-c(c("0_control","0_control","3","3","3","24","24","24","24","24_ece1","24_ece1","24_ece1","24_control","24_control"),     ## calb
                  c("3","3","3","24","24","24"),   ### cglab
                  c("3","3","3","24","24","24"),   ### cpar
                  c("3","3","3","24","24","24"),   ### ctrop
                  c("0_control","0_control","0_control","3","3","3","24","24","24","24_control","24_control","24_control"))  ### cglab_albumin
    
    batch<-c(rep("1",9),rep("2",3),rep("1",2),
             rep("1",6),
             rep("1",6),
             rep("1",6),
             rep("3",12))
    
    
    colData<-data.frame(batch=batch)
    
    set_to_remove_batch <- newSeqExpressionSet(as.matrix(all_human_data),phenoData = as.data.frame(colData, row.names=colnames(all_human_data)))
    
    removed_batch_effect <<- RUVg(set_to_remove_batch, as.character(overlap_for_batch$genes), k=1)
    pData(removed_batch_effect)
    
    pre_dds_no_batch <- DESeqDataSetFromMatrix(countData = normCounts(removed_batch_effect),colData = pData(removed_batch_effect), design = ~W_1+batch)
    
    dds_no_batch<-DESeq(pre_dds_no_batch)
    
    
    all_data_nobatch_norm<-vst(dds_no_batch,blind = T)
    vst_data_no_batch<-assay(all_data_nobatch_norm)
    pca <- prcomp(t(vst_data_no_batch))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    
    colData_mod<-data.frame(species=species, time_point=time_point, batch_effect=batch)
    df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2]*-1,colData_mod) ### do PC2 = pca$x[, 2]*-1 to invert values
    
    
    
    df_plotting$time_point <- factor(df_plotting$time_point, levels = c("0_control","3","24","24_ece1","24_control"))
    df_plotting$species <- factor(df_plotting$species, levels = c("calb","cglab","cpar","ctrop","control","cglab_albumin","control_albumin"))
    
    
    ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=time_point, shape=species), size=6)+
      geom_text_repel(label=as.character(rownames(df_plotting)), size = 4, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
      xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
      coord_fixed()+ theme_bw()+colScale+
      guides(colour = guide_legend(order = 1),
             size = guide_legend(order = 2))+
      labs(shape="Species", colour="Time point")+
      theme(axis.text = element_text(size = 16),axis.title = element_text(size=16),
            legend.text = element_text(size=14),legend.title = element_text(size=14))+
      scale_shape_manual(values = c(4,8,15,16,17,18,21))
    ggsave("~/Dropbox/Project_albumin/albumin_results_and_tables/human_pca_all_data_k_1.png", units="in", width=10, heigh=7, dpi=600)
    
  }
  
  else{
    
    ## uncomment below to exclude 1.5 h samples
    all_human_data<-all_human_data[,-c(3:21,27:29,36:59)]
    all_human_data<- all_human_data[,-c(8:10)]
    
    species<-c(rep("no albumin",10) ,rep("with albumin",12))
    time_point<-c(c("0_control","0_control","24_control","24_control"),     ## calb
                  c("3","3","3","24","24","24"),   ### cglab
                  c("0_control","0_control","0_control","3","3","3","24","24","24","24_control","24_control","24_control"))  ### cglab_albumin
    
    batch<-c(rep("1",10), rep("2",12))
    
    colData<-data.frame(batch=batch)
    
    set_to_remove_batch <- newSeqExpressionSet(as.matrix(all_human_data),phenoData = as.data.frame(colData, row.names=colnames(all_human_data)))
    
    removed_batch_effect <<- RUVg(set_to_remove_batch, as.character(overlap_for_batch$genes), k=1)
    pData(removed_batch_effect)
    
    pre_dds_no_batch <- DESeqDataSetFromMatrix(countData = normCounts(removed_batch_effect),colData = pData(removed_batch_effect), design = ~W_1+batch)
    
    dds_no_batch<-DESeq(pre_dds_no_batch)
    
    
    all_data_nobatch_norm<-vst(dds_no_batch,blind = T)
    vst_data_no_batch<-assay(all_data_nobatch_norm)
    pca <- prcomp(t(vst_data_no_batch))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    
    colData_mod<<-data.frame(species=species, time_point=time_point)
    df_plotting<<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2]*-1,colData_mod) ### do PC2 = pca$x[, 2]*-1 to invert values
    
    
    
    df_plotting$time_point <- factor(df_plotting$time_point, levels = c("0_control","3","24","24_control"))
    
    ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=time_point, shape=species), size=6)+
      geom_text_repel(label=as.character(rownames(df_plotting)), size = 4, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
      xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
      coord_fixed()+ theme_bw()+colScale+
      guides(colour = guide_legend(order = 1),
             size = guide_legend(order = 2))+
      labs(shape="Species", colour="Time point")+
      theme(axis.text = element_text(size = 16),axis.title = element_text(size=16),
            legend.text = element_text(size=14),legend.title = element_text(size=14))+
      ggsave("~/Dropbox/Project_albumin/albumin_results_and_tables/human_pca_cglab_data_k_1.png", units="in", width=10, heigh=7, dpi=600)
    
  }
  
}
### This generates PCA with RUVg k=1. To make the same for k=2 and k=3, change the value of k in :
### RUVg(set_to_remove_batch, as.character(overlap_for_batch$genes), k=1)
### and change "~ W_1 + batch" to "~ W_1 + W_2 + batch" for k=2, and "~ W_1+W_2+W_3+ batch" for k=3

### this generates PCA for all samples
deseq2_human("with dead and ece1") 

### this generates PCA for cglab samples only
deseq2_human("only cglab")


### fungal data PCAs


yeast_with_e_and_D<-read.table("~/Dropbox/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab/CGLAB_counts.txt",check.names = FALSE)

all_cglab <- cbind.data.frame(yeast_with_e_and_D,cglab_collapsed_counts)

all_cglab <- all_cglab[,-c(6,7,8)]

cglab_coldata<-data.frame(time=factor(c(rep("0_control",2), rep("3",3),rep("24", 3), rep("24_control",2),
                                        rep("0_control",3), rep("3",3), rep("24", 3), rep("24_control",3)), levels=c("0_control","3","24","24_control")))

species<-c(rep("no albumin",10) ,rep("with albumin",12))



pre_dds_cglab <- DESeqDataSetFromMatrix(all_cglab, colData=cglab_coldata, ~time)
dds_cglab <- DESeq(pre_dds_cglab)
all_data_norm<-vst(dds_cglab, blind = TRUE)
vst_data<-assay(all_data_norm)
pca <- prcomp(t(vst_data))
percentVar <- pca$sdev^2/sum(pca$sdev^2)


colData_cglab_to_plot <-cbind.data.frame(cglab_coldata,species)
df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],colData_cglab_to_plot) ### do PC2 = pca$x[, 2]*-1 to invert values

myColors <- brewer.pal(8,"Dark2")
names(myColors) <- c("0_control","3","24","24_control")
colScale <- scale_colour_manual(name = "Time point (h)",values = myColors)

ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=time, shape=species), size=6)+
  geom_text_repel(label=as.character(rownames(df_plotting)), size = 4, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
  coord_fixed()+ theme_bw()+colScale+
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2))+
  labs(shape="Species", colour="Time point")+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size=16),
        legend.text = element_text(size=14),legend.title = element_text(size=14))+
  ggsave("~/Dropbox/Project_albumin/albumin_results_and_tables/cglab_pca_no_batch_correction.png", units="in", width=10, heigh=7, dpi=300)

### batch corrected pca

res_lrt_cglab <- read.table("~/Dropbox/Project_Jena/Results_shared/Plots_and_Suppl/codes_for_github_tested_FINAL_to_save_for_website/cglab/res_lrt_cglab.txt")

non_DE_cglab<-rownames(res_lrt_cglab[!is.na(res_lrt_cglab$padj) & res_lrt_cglab$padj>0.5 & res_lrt_cglab$baseMean>10,])


batch<-c(rep("1",10) ,rep("2",12))

colData<-data.frame(batch=batch)

set_to_remove_batch <- newSeqExpressionSet(as.matrix(all_cglab),phenoData = as.data.frame(colData, row.names=colnames(all_cglab)))

removed_batch_effect <- RUVg(set_to_remove_batch, as.character(non_DE_cglab), k=1)
pData(removed_batch_effect)

pre_dds_no_batch <- DESeqDataSetFromMatrix(countData = normCounts(removed_batch_effect),colData = pData(removed_batch_effect), design = ~W_1+batch)

dds_no_batch<-DESeq(pre_dds_no_batch)


all_data_nobatch_norm<-vst(dds_no_batch,blind = T)
vst_data_no_batch<-assay(all_data_nobatch_norm)
pca <- prcomp(t(vst_data_no_batch))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

#colData_mod<-data.frame(species=species, time_point=time, batch_effect=batch)
df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2]*-1,colData_cglab_to_plot) ### do PC2 = pca$x[, 2]*-1 to invert values



df_plotting$time_point<-factor(df_plotting$time, levels = c("0_control","3","24","24_control"))

ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=time_point, shape=species), size=6)+
  geom_text_repel(label=as.character(rownames(df_plotting)), size = 4, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
  coord_fixed()+ theme_bw()+colScale+
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2))+
  labs(shape="Species", colour="Time point")+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size=16),
        legend.text = element_text(size=14),legend.title = element_text(size=14))+
  ggsave("~/Dropbox/Project_albumin/albumin_results_and_tables/cglab_pca_all_data_k_1.png", units="in", width=10, heigh=7, dpi=600)

