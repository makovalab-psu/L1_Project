require(IWTomics)
require(dendextend)

setwd("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/L1_polyA_filter/IWT_polyT_subset/")
#setwd("/Users/Bruce/Desktop/L1/IWTomics analysis_1000control_1000denovo")

# files with datasets and features
datasets=read.table("datasets_polyT.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasets=read.table("features_datasets_polyT.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# exclude low-resoution features (analyzed in IWT_L1_low_resolution.r)
features_datasets=features_datasets[1:3,]

# Save data
#regionsFeatures=IWTomicsData(datasets$file,features_datasets[,datasets$id],'center',
#                             datasets$id,datasets$name,features_datasets$id,features_datasets$name,
#                             path='files')
#save(regionsFeatures,file='L1_complete_filtered.RData')

# select only autosomes
load('L1_complete_filtered.RData')
index=lapply(regionsFeatures@regions,function(region) seqnames(region)!='chrX')
regionsFeatures@metadata$region_datasets$size=unlist(lapply(index,sum))
regionsFeatures@regions=GRangesList(mapply(function(region,ind) region[ind,],regionsFeatures@regions,index,SIMPLIFY=FALSE))
regionsFeatures@features=lapply(regionsFeatures@features,function(feature) mapply(function(feat,ind) feat[,which(ind)],feature,index,SIMPLIFY=FALSE))
regionsFeatures@length_features=lapply(regionsFeatures@length_features,function(feature) mapply(function(feat,ind) feat[which(ind)],feature,index,SIMPLIFY=FALSE))
validObject(regionsFeatures)
# number of windows in each dataset
lengthRegions(regionsFeatures)
#save(regionsFeatures,file='L1_autosomes.RData')


# number of 0 in the different features
# require that at least 10% of data are not zeros... 
# if this is not the case, smooth them with kernel smoothing and increasing bandwidth
zero_count=lapply(regionsFeatures@features,
                  function(feature){
                    counts=Reduce(rbind,lapply(feature,function(feat) c(sum(feat==0),sum(feat!=0),length(as.vector(feat)),length(unique(as.vector(feat))))))
                    colnames(counts)=c('0','>0','tot','distinct')
                    rownames(counts)=names(feature)
                    return(counts)
                  })
zero_count
zero_count_tot=Reduce(rbind,lapply(regionsFeatures@features,
                                   function(feature){
                                     feat=unlist(feature)
                                     count=c(sum(feat==0),sum(feat!=0),length(feat),length(unique(feat)))
                                     names(count)=c('0','>0','tot','distinct')
                                     return(count)
                                   }))
rownames(zero_count_tot)=names(zero_count)
zero_count_tot

many_zeros=names(which((zero_count_tot[,2]/zero_count_tot[,3]*100)<10))
many_zeros
#No actions needed

###Take into consideration stran info###
# de novo
L1_flanking_100=read.table('files/L1denovo_100kb_no_gaps_no_blacklisted_method2_15polyT.bed',header=FALSE,sep='\t')
names(L1_flanking_100)=c('chr','start','end','strand','barcode','site','annot','overlap_same','overlap_other')
L1_flanking_100=makeGRangesFromDataFrame(L1_flanking_100,starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
L1_flanking_100=L1_flanking_100[match(regionsFeatures@regions$L1denovo,L1_flanking_100,ignore.strand=TRUE)]
L1denovo_new=L1_flanking_100

# polymorphic
L1_flanking_100=read.table('files/L1Pol_100kb_no_gaps_no_blacklisted_method2.interval',header=FALSE,sep='\t')
L1_flanking_100=cbind(L1_flanking_100[,1:4],NA,NA,L1_flanking_100[,5:7])
names(L1_flanking_100)=c('chr','start','end','strand','barcode','site','annot','overlap_same','overlap_other')
L1_flanking_100=makeGRangesFromDataFrame(L1_flanking_100,starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
L1_flanking_100=L1_flanking_100[match(regionsFeatures@regions$L1Pol,L1_flanking_100,ignore.strand=TRUE)]
L1Pol_new=L1_flanking_100

# human specific
L1_flanking_100=read.table('files/L1HS_100kb_no_gaps_no_blacklisted_method2_withL1seq.interval',header=FALSE,sep='\t')
L1_flanking_100=cbind(L1_flanking_100[,1:4],NA,NA,L1_flanking_100[,5:7])
names(L1_flanking_100)=c('chr','start','end','strand','barcode','site','annot','overlap_same','overlap_other')
L1_flanking_100=makeGRangesFromDataFrame(L1_flanking_100,starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
L1_flanking_100=L1_flanking_100[match(regionsFeatures@regions$L1HS,L1_flanking_100,ignore.strand=TRUE)]
L1HS_new=L1_flanking_100

Control_new=regionsFeatures@regions$Control
mcols(Control_new)=data.frame(barcode=NA,site=NA,annot=NA,overlap_same=0,overlap_other=0)


# reverse the 100-kb region measurements for minus strand
regionsFeatures@regions=GRangesList(L1denovoFlasch=L1denovo_new_flasch,L1denovoSultana=L1denovo_new_sultana, L1denovo=L1denovo_new,L1Pol=L1Pol_new,L1HS=L1HS_new,Control=Control_new)
regionsFeatures@features=lapply(regionsFeatures@features,
                                function(feature){
                                  for(region in idRegions(regionsFeatures)){
                                    minus=which(strand(regions(regionsFeatures)[[region]])=='-')
                                    feature[[region]][,minus]=feature[[region]][seq(nrow(feature[[region]]),1),minus]
                                  }
                                  return(feature)
                                })
validObject(regionsFeatures)
save(regionsFeatures,file='L1_autosomes_zero_fixed_with_strand_filtered.RData')



load('L1_autosomes_zero_fixed_with_strand_filtered.RData')
# plot
# pdf('curves_autosomes_with_strand.pdf',width=10,height=8)
# plot(regionsFeatures,type='curves',id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
#      col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
# dev.off()

##Visulazation with boxplot
pdf('boxplot_autosomes_with_strand_filtered.pdf',width=10,height=7)
plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
     col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
dev.off()

# pdf('boxplot_autosomes_with_strand_denovo.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1denovo'),
#      col=c('black','red'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_pol.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1Pol'),
#      col=c('black','blue'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_hs.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1HS'),
#      col=c('black','green'),xlab='kb',ask=FALSE)
# dev.off()

# smooth all the features
already_smoothed=c("")
band=2
regionsFeatures_smoothed=smooth(regionsFeatures,id_features_subset=setdiff(idFeatures(regionsFeatures),already_smoothed),type='kernel',bandwidth=band)
save(regionsFeatures_smoothed,file='L1_autosomes_zero_fixed_with_strand_smoothed_filtered.RData')


load('L1_autosomes_zero_fixed_with_strand_smoothed_filtered.RData')
# plot
# pdf('curves_autosomes_with_strand_smoothed.pdf',width=10,height=8)
# plot(regionsFeatures_smoothed,type='curves',id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
#      col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
# dev.off()

pdf('boxplot_autosomes_with_strand_smoothed_filtered.pdf',width=10,height=7)
plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
     col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
dev.off()

# pdf('boxplot_autosomes_with_strand_smoothed_denovo.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1denovo'),
#      col=c('black','red'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_smoothed_pol.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1Pol'),
#      col=c('black','blue'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_smoothed_hs.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1HS'),
#      col=c('black','green'),xlab='kb',ask=FALSE)
# dev.off()


#################### IWT test ####################
#### all regions, zero fixed, with strand ####
##############################################
load('L1_autosomes_zero_fixed_with_strand_smoothed_filtered.RData')

result_mean=IWTomicsTest(regionsFeatures_smoothed,
                         id_region1=c("L1denovo","L1Pol","L1HS","L1denovo","L1denovo","L1Pol"),
                         id_region2=c("Control","Control","Control","L1Pol","L1HS","L1HS"),
                         statistics='mean',B=10000)
save(result_mean,file='L1_autosomes_results_smoothed_mean.RData')
load('L1_autosomes_results_smoothed_mean.RData')
pdf('IWT_autosomes_smoothed_mean.pdf',width=7,height=10)
plotTest(result_mean,col=c('red','blue','green','black'),
         scale_threshold=unlist(lapply(result_mean@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
dev.off()
plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),".pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',scale_threshold=10,
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),"_scale10.pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),".pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',scale_threshold=10,
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),"_scale10.pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)







