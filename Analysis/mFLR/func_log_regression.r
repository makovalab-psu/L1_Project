source('funLASSO_model_modified_for_L1.r')

require(IWTomics)
require(fda)
require(fda.usc)



organize_data_for_feature_selection <- function(IWTomics_object_low,IWTomics_object,index_dataset1,index_dataset2,
                                                index_low,index_scalars,index_functionals){
  # Create design matrix: intercept, low-resolution predictors, scalar predictors and functional predictors
  Z=NULL # predictors
  for(i in index_low){
    l1<-IWTomics_object_low@features[[i]][[index_dataset1]]
    ct<-IWTomics_object_low@features[[i]][[index_dataset2]]
    
    # standardize
    l1_ct_std<-scale(c(l1,ct))
    
    # add feature to design matrix
    Z=cbind(Z,l1_ct_std)
  }
  for(i in index_scalars){
    l1<-IWTomics_object@features[[i]][[index_dataset1]]
    ct<-IWTomics_object@features[[i]][[index_dataset2]]
    
    # compute mean in each region
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    
    # standardize
    l1_ct_mean_std<-scale(c(l1_mean,ct_mean))
    
    # add feature to design matrix
    Z=cbind(Z,l1_ct_mean_std)
  }
  colnames(Z)=c(idFeatures(IWTomics_object_low)[index_low],idFeatures(IWTomics_object)[index_scalars])
  for(i in index_functionals){
    l1<-IWTomics_object@features[[i]][[index_dataset1]]
    ct<-IWTomics_object@features[[i]][[index_dataset2]]
    x<-cbind(l1,ct)
    
    # standardize
    x_std=matrix(scale(as.vector(x)),nrow=nrow(x),ncol=ncol(x))
    
    # define the same basis for both x and b
    basis_common=create.bspline.basis(norder=3,breaks=seq(from=0.5,to=100.5,length.out=6))
    
    # create cross product matrix
    m_Jphi=inprod(basis_common,basis_common)
    
    # create cross product matrix for second derivatives (only needed in Gertheiss norm)
    m_Jphi2=inprod(basis_common,basis_common,Lfdobj1=int2Lfd(2),Lfdobj2=int2Lfd(2))
    
    # create functional object for x
    x.fd=Data2fd(argvals=1:100,y=x_std,basisobj=basis_common)
    
    # add feature to design matrix
    Z_new=t(x.fd$coefs)%*%m_Jphi
    colnames(Z_new)=paste0(idFeatures(IWTomics_object)[i],c('',seq_len(basis_common$nbasis)[-1]))
    Z=cbind(Z,Z_new)
  }
  Z=cbind(1, Z) #add intercept
  
  # create response
  y=rep(c(1,0),c(ncol(l1),ncol(ct)))
  
  # create groups of predictors, corresponding to different features
  groups=c(seq_len(length(index_low)+length(index_scalars)),
           rep(seq_along(index_functionals)+length(index_low)+length(index_scalars),each=basis_common$nbasis))
  
  return(list(Z=Z,y=y,groups=groups,m_Jphi=m_Jphi,m_Jphi2=m_Jphi2))
}

feature_selection_and_plot <- function(Z,y,groups,m_Jphi,m_Jphi2,type,lambda_all,lambda2_all=0,name=''){
  m_K=2 # number of classes
  prob=matrix(0, nr=nrow(Z), nc=m_K)  # class probability
  yhat=y
  
  # Penalty P1 (in our case P1 and P2 are the same penalty since we only have 2 classes)
  if(type %in% 1:2){
    bic_all=rep(NA,length(lambda_all))
    err_all=rep(NA,length(lambda_all))
    minbic=Inf
    for(i in seq_along(lambda_all)) {
      lambda=lambda_all[i]
      message('lambda ',signif(lambda,3))
      est = gLQA_mlogistic2_mod(Z, y, groups, lambda, type, lam2=0, a=1, b0=TRUE, m_Jphi=m_Jphi, m_Jphi2=m_Jphi2) # Multiclass functional logistic + group lasso (between variable only) estimated by LQA
      bic = BIC_multilog2_L1_mod(Z, y, groups, lambda, est$beta, type, lam2=0, a=1, b0=TRUE, m_Jphi=m_Jphi, m_Jphi2=m_Jphi2) # BIC for group SCAD
      bic_all[i]=bic
      err = 0 # error
      for(ii in 1:nrow(Z)){	
        for(l in 1:(m_K-1)){
          prob[ii,l]<-exp(Z[ii,]%*%est$beta[,l])/(1+sum(exp(Z[ii,]%*%est$beta)))
        }
        prob[ii, m_K] = 1 - sum(prob[ii, 1:(m_K-1)])
        switch (which.max(prob[ii, ]),
                "1" = (yhat[ii] = 1),
                "2" = (yhat[ii] = 0),
                (yhat[ii] = 1)
        )
        # classification error
        if (y[ii]!=yhat[ii]) {
          err = err+1
        }
      }
      err = err / nrow(Z)
      err_all[i]=err
      if (bic < minbic) {
        minbic = bic
        resbeta = est$beta # coefficients
        reserase = est$erase #a ctive set
        reslambda = lambda # regularization parameter
        reslambda2 = 0 # regularization parameter
        reserr = err # classification error
      }
      #message('bic ',round(bic,0))
      #message('err ',round(err,2)*100)
    } 
    rownames(reserase)=colnames(Z)
    
    plot(lambda_all,bic_all,type='l',xlab='Lambda',ylab='BIC',main=paste('BIC',name))
    abline(v=reslambda,col='red')
    plot(lambda_all,err_all*100,type='l',xlab='Lambda',ylab='Classification error',main=paste('Error',name))
    abline(v=reslambda,col='red')
  }
  if(type %in% 3){
    bic_all=matrix(NA,nrow=length(lambda_all),ncol=length(lambda2_all))
    err_all=matrix(NA,nrow=length(lambda_all),ncol=length(lambda2_all))
    minbic=Inf
    for(i in seq_along(lambda_all)) {
      lambda=lambda_all[i]
      message('lambda ',signif(lambda,3))
      for(j in seq_along(lambda2_all)){
        lambda2=lambda2_all[j]
        message('lambda2 ',signif(lambda2,3))
        est = gLQA_mlogistic2_mod(Z, y, groups, lambda, type, lam2=lambda2, a=1, b0=TRUE, m_Jphi=m_Jphi, m_Jphi2=m_Jphi2) # Multiclass functional logistic + group lasso (between variable only) estimated by LQA
        bic = BIC_multilog2_L1_mod(Z, y, groups, lambda, est$beta, type, lam2=lambda2, a=1, b0=TRUE, m_Jphi=m_Jphi, m_Jphi2=m_Jphi2) # BIC for group SCAD
        bic_all[i,j]=bic
        err = 0 # error
        for(ii in 1:nrow(Z)){	
          for(l in 1:(m_K-1)){
            prob[ii,l]<-exp(Z[ii,]%*%est$beta[,l])/(1+sum(exp(Z[ii,]%*%est$beta)))
          }
          prob[ii, m_K] = 1 - sum(prob[ii, 1:(m_K-1)])
          switch (which.max(prob[ii, ]),
                  "1" = (yhat[ii] = 1),
                  "2" = (yhat[ii] = 0),
                  (yhat[ii] = 1)
          )
          # classification error
          if (y[ii]!=yhat[ii]) {
            err = err+1
          }
        }
        err = err / nrow(Z)
        err_all[i,j]=err
        if (bic < minbic) {
          minbic = bic
          resbeta = est$beta # coefficients
          reserase = est$erase # active set
          reslambda = lambda # regularization parameter
          reslambda2 = lambda2 # regularization parameter
          reserr = err # classification error
        }
        #message('bic ',round(bic,0))
        #message('err ',round(err,2)*100)
      }
    } 
    rownames(reserase)=colnames(Z)
    
    contour(lambda_all,lambda2_all,resbic_all,xlab='Lambda',ylab='Lambda2',main=paste('BIC',name))
    points(reslambda,reslambda2,col='red')
    contour(lambda_all,lambda2_all,err_all*100,xlab='Lambda',ylab='Lambda2',main=paste('Error',name))
    points(reslambda,reslambda2,col='red')
  }
  
  return(list(minbic=minbic,resbeta=beta,reserase=reserase,reslambda=reslambda,reslambda2=reslambda2,reserr=reserr,
              bic_all=bic_all,err_all=err_all))
}

fit_reduced_model <- function(IWTomics_object_low,IWTomics_object,index_dataset1,index_dataset2,
                              index_low,index_scalars,index_functionals,index_selected_low,index_selected){
  # Fit model with the selected variables
  formula='y~'
  df=data.frame(y=rep(c(1,0),IWTomics_object@metadata$region_datasets$size[c(index_dataset1,index_dataset2)]))
  names_df=names(df)
  # Scalar predictors
  if(length(intersect(index_low,index_selected_low))>0){
    for(i in intersect(index_low,index_selected_low)){
      l1<-IWTomics_object_low@features[[i]][[index_dataset1]]
      ct<-IWTomics_object_low@features[[i]][[index_dataset2]]
      
      # add feature to data.frame
      df=cbind(df,c(l1,ct))
    }
    names(df)=c(names_df,idFeatures(IWTomics_object_low)[intersect(index_low,index_selected_low)])
    formula=paste0(formula,paste(idFeatures(IWTomics_object_low)[intersect(index_low,index_selected_low)],collapse='+'))
  }
  names_df=names(df)
  if(length(intersect(index_scalars,index_selected))>0){
    if(length(intersect(index_low,index_selected_low))>0)
      formula=paste0(formula,'+')
    for(i in intersect(index_scalars,index_selected)){
      l1<-IWTomics_object@features[[i]][[index_dataset1]]
      ct<-IWTomics_object@features[[i]][[index_dataset2]]
      
      # compute mean in each region
      l1_mean<-colMeans(l1)
      ct_mean<-colMeans(ct)
      
      # add feature to data.frame
      df=cbind(df,c(l1_mean,ct_mean))
    }
    names(df)=c(names_df,idFeatures(IWTomics_object)[intersect(index_scalars,index_selected)])
    formula=paste0(formula,paste(idFeatures(IWTomics_object)[intersect(index_scalars,index_selected)],collapse='+'))
  }
  data_regr=list(df=df)
  
  # Functional predictors
  basis=NULL
  if(length(intersect(index_functionals,index_selected))>0){
    if(length(intersect(index_scalars,index_selected))>0)
      formula=paste0(formula,'+')
    basis=vector('list')
    for(i in intersect(index_functionals,index_selected)){
      l1<-IWTomics_object@features[[i]][[index_dataset1]]
      ct<-IWTomics_object@features[[i]][[index_dataset2]]
      x<-cbind(l1,ct)
      
      # define the same basis for both x and b
      # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
      basis_common=create.bspline.basis(norder=3,breaks=seq(from=0.5,to=100.5,length.out=6))
      
      # add feature to data_regr
      data_regr[[paste(idFeatures(IWTomics_object)[i])]]=fdata(t(x),argvals=1:100,names=list(main=idFeatures(IWTomics_object)[i]))
      basis[[paste(idFeatures(IWTomics_object)[i])]]=basis_common
    }
    formula=paste0(formula,paste(idFeatures(IWTomics_object)[intersect(index_functionals,index_selected)],collapse='+'))
  }
  # Fit model
  flr_output=fregre.glm(formula(formula),data=data_regr,family=binomial(link="logit"),basis.x=basis,basis.b=basis)
  
  # deviance explained
  null_deviance=flr_output$null.deviance
  residual_deviance=flr_output$deviance
  dev.exp=(1-residual_deviance/null_deviance)
  dev.exp
  
  # compute RCDE (relative contribution to the deviance explained)
  pred_names=strsplit(as.character(flr_output$formula.ini)[3],'[+]')[[1]]
  RCDE_vector=NULL
  if(length(pred_names)>1){
    for(i in seq_along(pred_names)){
      reduced_formula=paste0('y~',paste0(pred_names[-i],collapse='+'))
      flr_output_reduced=fregre.glm(formula(reduced_formula),data=data_regr,family=binomial(link="logit"),basis.x=basis,basis.b=basis)
      deviance_reduced=flr_output_reduced$deviance
      RCDE_vector[i]=(((null_deviance-residual_deviance)-(null_deviance-deviance_reduced))/(null_deviance-residual_deviance))
    }
  }else{
    reduced_formula='y~1'
    flr_output_reduced=fregre.glm(formula(reduced_formula),data=data_regr,family=binomial(link="logit"),basis.x=basis,basis.b=basis)
    deviance_reduced=flr_output_reduced$deviance
    RCDE_vector=(((null_deviance-residual_deviance)-(null_deviance-deviance_reduced))/(null_deviance-residual_deviance))
  }
  names(RCDE_vector)=pred_names
  RCDE_vector=as.matrix(RCDE_vector)
  
  return(list(data_regr=data_regr,flr_output=flr_output,dev.exp=dev.exp,RCDE=RCDE_vector))
}

fit_two_predictors <- function(IWTomics_object_low,IWTomics_object,index_dataset1,index_dataset2,
                               index_low,index_scalars,index_functionals,index_selected_low,index_selected,index_predictor1){
  
  pred_name1=idFeatures(IWTomics_object)[index_predictor1]
  # Fit model with the selected variables
  formula='y~'
  df=data.frame(y=rep(c(1,0),IWTomics_object@metadata$region_datasets$size[c(index_dataset1,index_dataset2)]))
  names_df=names(df)
  # Scalar predictors
  if(length(intersect(index_low,index_selected_low))>0){
    for(i in intersect(index_low,index_selected_low)){
      l1<-IWTomics_object_low@features[[i]][[index_dataset1]]
      ct<-IWTomics_object_low@features[[i]][[index_dataset2]]
      
      # add feature to data.frame
      df=cbind(df,c(l1,ct))
    }
    names(df)=c(names_df,idFeatures(IWTomics_object_low)[intersect(index_low,index_selected_low)])
    formula=paste0(formula,paste(idFeatures(IWTomics_object_low)[intersect(index_low,index_selected_low)],collapse='+'))
  }
  names_df=names(df)
  if(length(intersect(index_scalars,index_selected))>0){
    if(length(intersect(index_low,index_selected_low))>0)
      formula=paste0(formula,'+')
    for(i in intersect(index_scalars,index_selected)){
      l1<-IWTomics_object@features[[i]][[index_dataset1]]
      ct<-IWTomics_object@features[[i]][[index_dataset2]]
      
      # compute mean in each region
      l1_mean<-colMeans(l1)
      ct_mean<-colMeans(ct)
      
      # add feature to data.frame
      df=cbind(df,c(l1_mean,ct_mean))
    }
    names(df)=c(names_df,idFeatures(IWTomics_object)[intersect(index_scalars,index_selected)])
    formula=paste0(formula,paste(idFeatures(IWTomics_object)[intersect(setdiff(index_scalars,index_predictor1),index_selected)],collapse='+'))
  }
  data_regr=list(df=df)
  # Functional predictors
  basis=NULL
  if(length(intersect(index_functionals,index_selected))>0){
    if(length(intersect(index_scalars,index_selected))>0)
      formula=paste0(formula,'+')
    basis=vector('list')
    for(i in intersect(index_functionals,index_selected)){
      l1<-IWTomics_object@features[[i]][[index_dataset1]]
      ct<-IWTomics_object@features[[i]][[index_dataset2]]
      x<-cbind(l1,ct)
      
      # define the same basis for both x and b
      # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
      basis_common=create.bspline.basis(norder=3,breaks=seq(from=0.5,to=100.5,length.out=6))
      
      # add feature to data_regr
      data_regr[[paste(idFeatures(IWTomics_object)[i])]]=fdata(t(x),argvals=1:100,names=list(main=idFeatures(IWTomics_object)[i]))
      basis[[paste(idFeatures(IWTomics_object)[i])]]=basis_common
    }
    formula=paste0(formula,paste(idFeatures(IWTomics_object)[intersect(setdiff(index_functionals,index_predictor1),index_selected)],collapse='+'))
  }
  
  # Fit model
  flr_output=fregre.glm(formula(formula),data=data_regr,family=binomial(link="logit"),basis.x=basis,basis.b=basis)
  
  # compute R-squared for predictor1 alone, and predictor1 plus another predictor
  pred_names=strsplit(as.character(flr_output$formula.ini)[3],'[+]')[[1]]
  
  flr_output_pred1=fregre.glm(formula(paste0('y~',pred_name1)),data=data_regr,family=binomial(link="logit"),basis.x=basis,basis.b=basis)
  null_deviance=flr_output_pred1$null.deviance
  residual_deviance=flr_output_pred1$deviance
  R2_vector=(1-residual_deviance/null_deviance)
  if(length(pred_names)>1){
    for(i in seq_along(pred_names)){
      reduced_formula=paste0('y~',pred_name1,'+',pred_names[i])
      flr_output_reduced=fregre.glm(formula(reduced_formula),data=data_regr,family=binomial(link="logit"),basis.x=basis,basis.b=basis)
      null_deviance=flr_output_reduced$null.deviance
      residual_deviance=flr_output_reduced$deviance
      R2_vector[i+1]=(1-residual_deviance/null_deviance)
      
      reduced_formula=paste0('y~',pred_names[i])
      flr_output_reduced=fregre.glm(formula(reduced_formula),data=data_regr,family=binomial(link="logit"),basis.x=basis,basis.b=basis)
      null_deviance=flr_output_reduced$null.deviance
      residual_deviance=flr_output_reduced$deviance
      R2_vector[i+1]=R2_vector[i+1]-(1-residual_deviance/null_deviance)
    }
  }
  names(R2_vector)=c(pred_name1,pred_names)
  R2_vector=as.matrix(R2_vector)
  
  
  return(list(data_regr=data_regr,flr_output_pred1=flr_output_pred1,R2=R2_vector))
}




# Select functional predictors by the localization table created based on the boxplots & IWT results
localization_all<-read.table("localization_table.txt",header=TRUE)

func_1<-as.vector(localization_all$test1)
func_2<-as.vector(localization_all$test2)
func_3<-as.vector(localization_all$test3)
func_4<-as.vector(localization_all$test4)
func_5<-as.vector(localization_all$test5)
func_6<-as.vector(localization_all$test6)

# Select low resolution predictors by the localization table created based on the test results
localization_low<-read.table("localization_table_LowFeatures.txt",header=TRUE)

low_1<-as.vector(localization_low$test1)
low_2<-as.vector(localization_low$test2)
low_3<-as.vector(localization_low$test3)
low_4<-as.vector(localization_low$test4)
low_5<-as.vector(localization_low$test5)
low_6<-as.vector(localization_low$test6)






setwd('./random1')
load(paste0("L1_transformed_random_",r,'.RData'))
load(paste0("L1_lowFeatures_transformed_random",r,'.RData'))

## Change the feature names to the same as IWT plot
feature_names<-as.vector(result_mean@metadata$feature_datasets$name)

feature_num<- as.vector(c(1:45))
feature<-as.data.frame(matrix(ncol=1))
for (i in 1:45){
  feature<-rbind(feature,paste("feature_", i, sep = ""))
}
feature<-feature[-1,]
feature_table<-data.frame(feature_num,feature,feature_names)

# change "LINE_L2&L3" to "LINE_L2_L3"
row.names(result_mean@metadata$feature_datasets)[33]="LINE_L2_L3"
names(result_mean@features)[33]="LINE_L2_L3"







####### Multiple logistic regression (both scalar and functional predictors)



### Comparison 1 (denovo vs control)
####################################
comp_scal<-which(func_1=="n")
comp_func<-which(func_1=="y")
comp_low<-which(low_1=="n")
index1=1
index2=4

# Create design matrix, response and groups of predictors
data=organize_data_for_feature_selection(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func)


## Variable selection
type=1 # 1:Euclidean 2:Aguilera-Morillo 3:Gertheiss (this need the second parameter lambda2)

# type 1 and 2
lambda_all=seq(0.02,0.06,0.001)
lambda2_all=0

pdf(paste0('comp1_feature_selection_type',type,'.pdf'),width=5,height=5)
selection_results=feature_selection_and_plot(data$Z,data$y,data$groups,data$m_Jphi,data$m_Jphi2,
                                             type,lambda_all,lambda2_all,name='De novo vs Control')
dev.off()

# Selected variables
comp_select_low=which(idFeatures(regionsFeatures) %in% rownames(selection_results$reserase)[selection_results$reserase==1])
comp_select=which(idFeatures(result_mean) %in% rownames(selection_results$reserase)[selection_results$reserase==1])

# Fit model with the selected variables, compute deviance explained and RCDE
flr_results=fit_reduced_model(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select)
summary(flr_results$flr_output)
round(flr_results$dev.exp*100,2)
round(flr_results$RCDE*100,2)



# Compute R2 added by L1 targets to model with one other predictor
comp1=32 # L1 targets
L1targets_effect=fit_two_predictors(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select,comp1)
round(L1targets_effect$R2*100,2)

# Compute R2 added by GC content to model with one other predictor
comp1=24 # GC content
L1targets_effect=fit_two_predictors(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select,comp1)
round(L1targets_effect$R2*100,2)

save(data,lambda_all,lambda2_all,selection_results,comp_select_low,comp_select,flr_results,file=paste0('comp1_final_model_type',type,'.RData'))

#---------------------------





### Comparison 2 (pol vs control)
###################################
comp_scal<-which(func_2=="n")
comp_func<-which(func_2=="y")
comp_low<-which(low_2=="n")
index1=2
index2=4

# Create design matrix, response and groups of predictors
data=organize_data_for_feature_selection(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func)


## Variable selection
type=1 # 1:Euclidean 2:Aguilera-Morillo 3:Gertheiss (this need the second parameter lambda2)

# type 1 and 2
lambda_all=seq(0.02,0.06,0.001)
lambda2_all=0

pdf(paste0('comp2_feature_selection_type',type,'.pdf'),width=5,height=5)
selection_results=feature_selection_and_plot(data$Z,data$y,data$groups,data$m_Jphi,data$m_Jphi2,
                                             type,lambda_all,lambda2_all,name='Polymorphic vs Control')
dev.off()

# Selected variables
comp_select_low=which(idFeatures(regionsFeatures) %in% rownames(selection_results$reserase)[selection_results$reserase==1])
comp_select=which(idFeatures(result_mean) %in% rownames(selection_results$reserase)[selection_results$reserase==1])

# Fit model with the selected variables, compute deviance explained and RCDE
flr_results=fit_reduced_model(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select)
summary(flr_results$flr_output)
round(flr_results$dev.exp*100,2)
round(flr_results$RCDE*100,2)



# Compute R2 added by L1 targets to model with one other predictor
comp1=32 # L1 targets
L1targets_effect=fit_two_predictors(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select,comp1)
round(L1targets_effect$R2*100,2)

# Compute R2 added by GC content to model with one other predictor
comp1=24 # GC content
L1targets_effect=fit_two_predictors(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select,comp1)
round(L1targets_effect$R2*100,2)

save(data,lambda_all,lambda2_all,selection_results,comp_select_low,comp_select,flr_results,file=paste0('comp2_final_model_type',type,'.RData'))
#---------------------------





### Comparison 3 (hs vs control)
###################################
comp_scal<-which(func_3=="n")
comp_func<-which(func_3=="y")
comp_low<-which(low_3=="n")
index1=3
index2=4

# Create design matrix, response and groups of predictors
data=organize_data_for_feature_selection(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func)


## Variable selection
type=1 # 1:Euclidean 2:Aguilera-Morillo 3:Gertheiss (this need the second parameter lambda2)

# type 1 and 2
lambda_all=seq(0.06,0.10,0.001)
lambda2_all=0

pdf(paste0('comp3_feature_selection_type',type,'.pdf'),width=5,height=5)
selection_results=feature_selection_and_plot(data$Z,data$y,data$groups,data$m_Jphi,data$m_Jphi2,
                                             type,lambda_all,lambda2_all,name='Human specific vs Control')
dev.off()

# Selected variables
comp_select_low=which(idFeatures(regionsFeatures) %in% rownames(selection_results$reserase)[selection_results$reserase==1])
comp_select=which(idFeatures(result_mean) %in% rownames(selection_results$reserase)[selection_results$reserase==1])

# Fit model with the selected variables, compute deviance explained and RCDE
flr_results=fit_reduced_model(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select)
summary(flr_results$flr_output)
round(flr_results$dev.exp*100,2)
round(flr_results$RCDE*100,2)



# Compute R2 added by L1 targets to model with one other predictor
comp1=32 # L1 targets
L1targets_effect=fit_two_predictors(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select,comp1)
round(L1targets_effect$R2*100,2)

# Compute R2 added by GC content to model with one other predictor
comp1=24 # GC content
L1targets_effect=fit_two_predictors(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select,comp1)
round(L1targets_effect$R2*100,2)

save(data,lambda_all,lambda2_all,selection_results,comp_select_low,comp_select,flr_results,file=paste0('comp3_final_model_type',type,'.RData'))
#---------------------------





### Comparison 4 (pol vs denovo)
###################################
comp_scal<-which(func_4=="n")
comp_func<-which(func_4=="y")
comp_low<-which(low_4=="n")
index1=2
index2=1

# Create design matrix, response and groups of predictors
data=organize_data_for_feature_selection(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func)


## Variable selection
type=1 # 1:Euclidean 2:Aguilera-Morillo 3:Gertheiss (this need the second parameter lambda2)

# type 1 and 2
lambda_all=seq(0.14,0.18,0.001)
lambda2_all=0

pdf(paste0('comp4_feature_selection_type',type,'.pdf'),width=5,height=5)
selection_results=feature_selection_and_plot(data$Z,data$y,data$groups,data$m_Jphi,data$m_Jphi2,
                                             type,lambda_all,lambda2_all,name='Polymorphic vs De novo')
dev.off()

# Selected variables
comp_select_low=which(idFeatures(regionsFeatures) %in% rownames(selection_results$reserase)[selection_results$reserase==1])
comp_select=which(idFeatures(result_mean) %in% rownames(selection_results$reserase)[selection_results$reserase==1])

# Fit model with the selected variables, compute deviance explained and RCDE
flr_results=fit_reduced_model(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select)
summary(flr_results$flr_output)
round(flr_results$dev.exp*100,2)
round(flr_results$RCDE*100,2)

save(data,lambda_all,lambda2_all,selection_results,comp_select_low,comp_select,flr_results,file=paste0('comp4_final_model_type',type,'.RData'))
#---------------------------





### Comparison 5 (hs vs denovo)
###################################
comp_scal<-which(func_5=="n")
comp_func<-which(func_5=="y")
comp_low<-which(low_5=="n")
index1=3
index2=1

# Create design matrix, response and groups of predictors
data=organize_data_for_feature_selection(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func)


## Variable selection
type=1 # 1:Euclidean 2:Aguilera-Morillo 3:Gertheiss (this need the second parameter lambda2)

# type 1 and 2
lambda_all=seq(0.13,0.15,0.001)
lambda2_all=0

pdf(paste0('comp5_feature_selection_type',type,'.pdf'),width=5,height=5)
selection_results=feature_selection_and_plot(data$Z,data$y,data$groups,data$m_Jphi,data$m_Jphi2,
                                             type,lambda_all,lambda2_all,name='Human specific vs De novo')
dev.off()

# Selected variables
comp_select_low=which(idFeatures(regionsFeatures) %in% rownames(selection_results$reserase)[selection_results$reserase==1])
comp_select=which(idFeatures(result_mean) %in% rownames(selection_results$reserase)[selection_results$reserase==1])

# Fit model with the selected variables, compute deviance explained and RCDE
flr_results=fit_reduced_model(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select)
summary(flr_results$flr_output)
round(flr_results$dev.exp*100,2)
round(flr_results$RCDE*100,2)

save(data,lambda_all,lambda2_all,selection_results,comp_select_low,comp_select,flr_results,file=paste0('comp5_final_model_type',type,'.RData'))
#---------------------------





### Comparison 6 (hs vs pol)
###################################
comp_scal<-which(func_6=="n")
comp_func<-which(func_6=="y")
comp_low<-which(low_6=="n")
index1=3
index2=2

# Create design matrix, response and groups of predictors
data=organize_data_for_feature_selection(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func)


## Variable selection
type=1 # 1:Euclidean 2:Aguilera-Morillo 3:Gertheiss (this need the second parameter lambda2)

# type 1 and 2
lambda_all=seq(0.14,0.18,0.001)
lambda2_all=0

pdf(paste0('comp6_feature_selection_type',type,'.pdf'),width=5,height=5)
selection_results=feature_selection_and_plot(data$Z,data$y,data$groups,data$m_Jphi,data$m_Jphi2,
                                             type,lambda_all,lambda2_all,name='Human specific vs Polymorphic')
dev.off()

# Selected variables
comp_select_low=which(idFeatures(regionsFeatures) %in% rownames(selection_results$reserase)[selection_results$reserase==1])
comp_select=which(idFeatures(result_mean) %in% rownames(selection_results$reserase)[selection_results$reserase==1])

# Fit model with the selected variables, compute deviance explained and RCDE
flr_results=fit_reduced_model(regionsFeatures,result_mean,index1,index2,comp_low,comp_scal,comp_func,comp_select_low,comp_select)
summary(flr_results$flr_output)
round(flr_results$dev.exp*100,2)
round(flr_results$RCDE*100,2)

save(data,lambda_all,lambda2_all,selection_results,comp_select_low,comp_select,flr_results,file=paste0('comp6_final_model_type',type,'.RData'))
#---------------------------

