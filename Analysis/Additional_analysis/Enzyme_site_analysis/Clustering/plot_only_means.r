plot_only_means <- function(x,type='boxplot',method='pearson',include_zero=TRUE,
                              N_regions=pmin(lengthRegions(x),ifelse(type=='curves',10,ifelse(type=='pairs',1000,+Inf))),
                              probs=c(0.25,0.5,0.75),average=TRUE,size=TRUE,
                              id_regions_subset=idRegions(x),id_features_subset=idFeatures(x),
                              log_scale=FALSE,log_shift=0,col=1+seq_along(id_regions_subset),
                              plot=TRUE,ask=TRUE,xlab='Windows',ylim=NULL,...){
  if(sum(!(id_regions_subset %in% idRegions(x))))
    stop('invalid id_regions_subset. The region datasets provided are not listed in x.')
  if(sum(!(id_features_subset %in% idFeatures(x))))
    stop('invalid id_features_subset. The features provided are not listed in x.')
  if(!(type %in% c('curves','boxplot','pairs','pairsSmooth')))
    stop('invalid plot type \'',type,'\'. Available types are \'curves\', \'boxplot\', \'pairs\' and \'pairsSmooth\'.')
  if(!(method %in% c('pearson','kendall','spearman')))
    stop('invalid method type \'',method,'\'. Available methods are \'pearson\', \'kendall\' and \'spearman\'.')
  
  x=x[id_regions_subset,id_features_subset]
  features_plot=features(x)
  if(log_scale){
    features_plot=lapply(features_plot,lapply,'+',log_shift)
    features_plot=lapply(features_plot,lapply,log)
    if(sum(unlist(lapply(features_plot,lapply,is.infinite))))
      stop('logarithm of 0.')
  }
  if(type %in% c('pairs','pairsSmooth')){
    if(length(unique(resolution(x)))!=1)
      stop('type\'',type,'\' but selected features with different resolution. Smooth data first to have the same resolution.')
    N_regions=rep(N_regions,length.out=length(id_regions_subset))
    N_regions=pmin(lengthRegions(x),N_regions)
    index_plot=mapply(sample,lengthRegions(x),N_regions,SIMPLIFY=FALSE)
    shuffle=sample(sum(N_regions))
    features_plot=lapply(features_plot,function(feature_plot) Reduce(cbind,mapply(function(feature,index_plot) feature[,index_plot],feature_plot,index_plot,SIMPLIFY=FALSE))[,shuffle])
    col=rep(col,length.out=length(N_regions))
    col_plot=rep(rep(col,N_regions)[shuffle],each=nrow(features_plot[[1]]))
    features_plot=do.call(cbind,lapply(features_plot,as.vector))
    features_cor=cor(features_plot[!is.na(rowSums(features_plot)),],method=method)
    z=list(features_plot=features_plot,features_cor=features_cor,type=type)
    if(plot){
      panel.cor <- function(x,y,digits=2,method_cor=method,prefix="",cex.cor, ...){
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        index=(!is.na(x))&(!is.na(y))
        r <- cor(x[index], y[index],method=method_cor)
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste0(prefix, txt)
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor*(0.4+0.8*abs(r)))
      }
      devAskNewPage(ask)
      if(type=='pairs'){
        upper.panel=points
      }else{
        upper.panel=function(...) {par(new=TRUE);smoothScatter(...,nrpoints=0,colramp=colorRampPalette(c("white",col_plot[1])),add=TRUE)}
      }
      pairs(features_plot,main='Features correlation',labels=nameFeatures(x),col=col_plot,pch=3,lower.panel=panel.cor,upper.panel=upper.panel,...)
      invisible(z)
    }else{
      z
    }
  }else{
    if(alignment(x)=='left')
      x_plot=lapply(features_plot,function(feature_plot) 0.5:(nrow(feature_plot[[1]])-0.5))
    if(alignment(x)=='right')
      x_plot=lapply(features_plot,function(feature_plot) -(nrow(feature_plot[[1]])-0.5):(-0.5))
    if(alignment(x)=='center')
      x_plot=lapply(features_plot,function(feature_plot) seq_len(nrow(feature_plot[[1]]))-nrow(feature_plot[[1]])/2-0.5)
    if(alignment(x)=='scale'){
      length_features_plot=mapply(function(feature_plot,length_feature) mapply(function(feature,length) pmax(apply(feature,2,function(feature) rev(which(!is.na(feature)))[1]),length),
                                                                               feature_plot,length_feature,SIMPLIFY=FALSE),
                                  features_plot,x@length_features,SIMPLIFY=FALSE)
      if(sum(unlist(lapply(length_features_plot,function(length_feature) length(unique(unlist(length_feature)))!=1)))){
        if(type=='boxplot')
          stop('type \'boxplot\' is incompatible with \'scale\' alignment and regions of different length. Smooth data first.')
        if(average){
          warning('average=TRUE is incompatible with \'scale\' alignment and regions of different length. Setting average=FALSE.')
          average=FALSE
        }
        if(size){
          warning('size=TRUE is incompatible with \'scale\' alignment and regions of different length. Setting size=FALSE.')
          size=FALSE
        }
        x_plot=lapply(length_features_plot,
                      function(length_feature){
                        length_max=max(unlist(length_feature))
                        x_plot=lapply(length_feature,function(length) do.call(cbind,lapply(length,function(length) c(seq(0,1,length.out=length),rep(NA,length_max-length)))))
                        return(x_plot)
                      })
      }else{
        x_plot=lapply(features_plot,function(feature_plot) seq(0,1,length.out=nrow(feature_plot[[1]])))
      }
    }
    if(length(col)!=length(id_regions_subset)){
      warning('number of colors in \'col\' different from the number of region datasets considered.')
      col=rep(col,length.out=length(id_regions_subset))
    }
    if(average)
      features_average=lapply(features_plot,function(feature_plot) Reduce(cbind,lapply(feature_plot,rowMeans,na.rm=TRUE)))
    if(size)
      features_position_size=lapply(features_plot,function(feature_plot) do.call(cbind,lapply(rev(feature_plot),function(feature) rowSums(!is.na(feature)))))
    if(type=='curves'){
      N_regions=rep(N_regions,length.out=length(id_regions_subset))
      N_regions=pmin(lengthRegions(x)[id_regions_subset],N_regions)
      index_plot=mapply(sample,lengthRegions(x)[id_regions_subset],N_regions,SIMPLIFY=FALSE)
      shuffle=sample(sum(N_regions))
      features_plot=lapply(features_plot,function(feature_plot) Reduce(cbind,mapply(function(feature,index_plot) as.matrix(feature[,index_plot]),feature_plot,index_plot,SIMPLIFY=FALSE))[,shuffle])
      if(is.list(x_plot[[1]]))
        x_plot=lapply(x_plot,function(x_plot) Reduce(cbind,mapply(function(x,index_plot) x[,index_plot],x_plot,index_plot,SIMPLIFY=FALSE))[,shuffle])
      col_plot=rep(col,N_regions)[shuffle]
    }
    if(type=='boxplot'){
      #features_plot=lapply(features_plot,function(feature_plot) lapply(feature_plot,function(feature) t(apply(feature,1,quantile,na.rm=TRUE,probs=probs))))
      col_plot=lapply(col,function(col) rgb(colorRamp(c('white',col))((1:4)/4)[-1,],alpha=c(80,230),maxColorValue=255))
      names(col_plot)=id_regions_subset
    }
    
    if(plot){
      devAskNewPage(ask)
      if(size){
        layout(matrix(1:3,nrow=3),heights=c(5,1,1))
        mar.left=7
      }else{
        mar.left=4
      }
      par(oma=c(0,0,0,8))
      if(is.null(ylim)){
        if(average){
          #ylim=mapply(function(feature,average) pmin(range(c(unlist(feature),unlist(average)),na.rm=TRUE),c(0,+Inf)),features_plot,features_average,SIMPLIFY=FALSE)
          if(include_zero){
            ylim=lapply(features_average,function(average) pmin(range(unlist(average),na.rm=TRUE),c(0,+Inf)))
          }else{
            ylim=lapply(features_average,function(average) pmin(range(unlist(average),na.rm=TRUE),c(+Inf,+Inf)))
          }
        }else{
          ylim=lapply(features_plot,function(feature) pmin(range(unlist(feature),na.rm=TRUE),c(0,+Inf)))
        }
      }else{
        ylim=lapply(features_plot,function(feature) ylim)
      }
      for(id_feature in id_features_subset){
        par(mar=c(5,mar.left,4,1.5))
        if(type=='curves'){
          matplot(x_plot[[id_feature]],features_plot[[id_feature]],type='l',col=col_plot,ylim=ylim[[id_feature]],
                  main=nameFeatures(x)[id_feature],xlab=xlab,
                  ylab=paste0(ifelse(log_scale,'log ',''),nameFeatures(x)[id_feature]),...)
        }
        if(type=='boxplot'){
          plot(1,type="n",xlim=range(x_plot[[id_feature]]),ylim=ylim[[id_feature]],
               main=nameFeatures(x)[id_feature],xlab=xlab,
               ylab=paste0(ifelse(log_scale,'log ',''),nameFeatures(x)[id_feature]),...)
          for(id_region in id_regions_subset){
            #polygon(c(x_plot[[id_feature]],rev(x_plot[[id_feature]])),c(features_plot[[id_feature]][[id_region]][,1],rev(features_plot[[id_feature]][[id_region]][,length(probs)])),
            #        col=col_plot[[id_region]][1],border=FALSE)
            #for(i in seq_along(probs))
            #  lines(x_plot[[id_feature]],features_plot[[id_feature]][[id_region]][,i],col=col_plot[[id_region]][2],lty=2,lwd=1.5)
          }
        }
        if(average)
          matplot(x_plot[[id_feature]],features_average[[id_feature]],type='l',col=col,lty=1,lwd=2,add=TRUE)
        args=as.list(match.call())
        if(is.null(args$cex)){
          cex=ifelse(is.null(args$cex.lab),1,args$cex.lab)
        }else{
          cex=args$cex
        }
        legend(par('usr')[2],mean(par('usr')[3:4]),legend=nameRegions(x)[id_regions_subset],xpd=NA,bty='n',lty=1,lwd=2,col=col,yjust=0.5,cex=cex)
        if(size){
          par(mar=c(1,mar.left,2,1.5))
          image(x_plot[[id_feature]],seq_along(id_regions_subset),features_position_size[[id_feature]],
                col=cm.colors(101),xlim=par('usr')[1:2],ylim=range(seq_along(id_regions_subset))+c(-0.5,0.5),axes=FALSE,xlab='',ylab='',...)
          axis(side=3,...)
          axis(side=2,at=seq_along(id_regions_subset),labels=rev(nameRegions(x)[id_regions_subset]),tick=FALSE,las=1,line=-1,...)
          x.rect=range(x_plot[[id_feature]])+c(-1,1)*diff(x_plot[[id_feature]])[1]/2
          rect(x.rect[1],seq_along(id_regions_subset)-0.5,x.rect[2],seq_along(id_regions_subset)+0.5,border='black')
          par(mar=c(5,mar.left,0.5,1.5))
          image(-50:50,1,as.matrix(seq.int(101)),xlim=c(-100,100),axes=FALSE,xlab='Sample size',ylab='',col=cm.colors(101),...)
          labels=seq(min(features_position_size[[id_feature]]),max(features_position_size[[id_feature]]),length.out=5)
          if(length(unique(labels))==1)
            labels=-2:2+labels
          axis(side=1,at=c(-50,-25,0,25,50),labels=labels,...)
          rect(-50.5,par('usr')[3],50.5,par('usr')[4],border='black')
        }
      }
      if(!average)
        features_average=NULL
      if(!size)
        features_position_size=NULL
      z=list(x_plot=x_plot,features_plot=features_plot,features_average=features_average,features_position_size=features_position_size,type=type,col=col,col_plot=col_plot)
      invisible(z)
    }else{
      if(!average)
        features_average=NULL
      if(!size)
        features_position_size=NULL
      z=list(x_plot=x_plot,features_plot=features_plot,features_average=features_average,features_position_size=features_position_size,type=type,col=col,col_plot=col_plot)
      z
    }
  }
}
