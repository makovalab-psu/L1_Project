#Clustering with additional features: MSPI and TaqI
require(IWTomics)
require(dendextend)

setwd("~/Google_Drive/L1_Project/Analysis/Additional_analysis/Enzyme_site_analysis/clustering/")
source('plot_only_means.r') # function plot_only_means to plot only means (no boxplots)

# Load .RData with MSPI and TaqI included as additional features
load('L1_complete_enzyme_clustering.RData')
# Select only autosomes
index=lapply(regionsFeatures@regions,function(region) seqnames(region)!='chrX')
regionsFeatures@metadata$region_datasets$size=unlist(lapply(index,sum))
regionsFeatures@regions=GRangesList(mapply(function(region,ind) region[ind,],regionsFeatures@regions,index,SIMPLIFY=FALSE))
regionsFeatures@features=lapply(regionsFeatures@features,function(feature) mapply(function(feat,ind) feat[,which(ind)],feature,index,SIMPLIFY=FALSE))
regionsFeatures@length_features=lapply(regionsFeatures@length_features,function(feature) mapply(function(feat,ind) feat[which(ind)],feature,index,SIMPLIFY=FALSE))
validObject(regionsFeatures)
# number of windows in each dataset
lengthRegions(regionsFeatures)
save(regionsFeatures,file='L1_autosomes_enzyme_clustering.RData')

################################################ 
#########Plot MspI and TaqI separately##########
################################################ 
# load data
load('L1_autosomes_enzyme_clustering.RData')
# select only feature we want to analyze: testis expression and sperm-hypomethylation
idFeatures_select=c("H2AFZ_signal","H3K27ac_signal","H4K20me1_signal","H3K36me3_signal","H3K4me1_signal",
                    "H3K4me2_signal","H3K4me3_signal","H3K79me2_signal","H3K9ac_signal","H3K9me3_signal",
                    "H3K27me3_signal","CTCF_signal","DNase_DHS_signal","RNA_PolII",
                    "Quadruplex","A_Phased","Direct_Repeats","Inverted_Repeats","Mirror_Repeats",
                    "Z_DNA","Most_Conserved","Exons","Introns","GC_Content",
                    #(remove AT_content from the analysis)"AT_Content",
                    "Mononucleotide","Morethan1_nuc","DNA_Transposons","SINE_Alu","SINE_MIR",
                    "LTR_Elements","L1_Targets","LINE_L2&L3","CpG_Islands",
                    "5hMc","Sperm_hypometh","Rep_origin","Recombination_Hot",
                    #"Exon_Expression",
                    "Gene_Expression",
                    #"Tracript_Expression", "Testis_Expression",
                    "chh_meth","chg_meth","cpg_meth","MSPI_Count","TaqI_Count")

#idFeatures_select=c("Sperm_hypometh","Testis_Expression","Gene_Expression")
regionsFeatures=regionsFeatures[,idFeatures_select]

# clustering using spearman correlation
clustering_spearman <- function(data,file_name,h,labels){
  correlation <- cor(data,method="spearman");
  abs.corr=as.dist(1-abs(correlation))
  fit <- hclust(abs.corr)
  pdf(paste0("clustering_spearman",file_name,".pdf"),width=12,height=12)
  par(mar=c(4,1,1,14)+0.1)
  library("dendextend")
  dend=as.dendrogram(fit,hang=0.05)
  dend_cut=cutree(dend,h=h)[order.dendrogram(dend)]
  dend_col=rep(1,length(dend_cut))
  cluster=which(table(dend_cut)>1)
  for(k in seq_along(cluster)){
    dend_col[dend_cut==cluster[k]]=k+2
  }
  labels_colors(dend)=dend_col
  labels(dend)=labels[labels(dend)]
  plot(dend,xlab= "1-|Spearman's correlation|", ylab="",main=NULL,cex=0.7,cex.lab=1.5,horiz=TRUE,xlim=c(1,0))
  lines(c(h,h),c(-2,nrow(correlation)+1),col='red',lty='dashed',lwd=2)
  dev.off()
}
# Plot all regions
features_pairs=plot(regionsFeatures,type='pairs',N_regions=lengthRegions(regionsFeatures),plot=FALSE)
clustering_spearman(features_pairs$features_plot,"_autosomes_all_separate",h=0.2,labels=nameFeatures(regionsFeatures))

# # de novo
# features_pairs=plot(regionsFeatures,type='pairs',id_regions_subset='L1denovo',N_regions=lengthRegions(regionsFeatures),plot=FALSE)
# clustering_spearman(features_pairs$features_plot,"_autosomes_denovo_separate",h=0.2,labels=nameFeatures(regionsFeatures))
# # polymorphic
# features_pairs=plot(regionsFeatures,type='pairs',id_regions_subset='L1Pol',N_regions=lengthRegions(regionsFeatures),plot=FALSE)
# clustering_spearman(features_pairs$features_plot,"_autosomes_pol_separate",h=0.2,labels=nameFeatures(regionsFeatures))
# # human specific
# features_pairs=plot(regionsFeatures,type='pairs',id_regions_subset='L1HS',N_regions=lengthRegions(regionsFeatures),plot=FALSE)
# clustering_spearman(features_pairs$features_plot,"_autosomes_hs_separate",h=0.2,labels=nameFeatures(regionsFeatures))
# # control
# features_pairs=plot(regionsFeatures,type='pairs',id_regions_subset='Control',N_regions=lengthRegions(regionsFeatures),plot=FALSE)
# clustering_spearman(features_pairs$features_plot,"_autosomes_controls_separate",h=0.2,labels=nameFeatures(regionsFeatures))
# 
