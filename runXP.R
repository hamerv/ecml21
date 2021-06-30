library(glmnet)
library(doParallel)
library(datamicroarray)
library(ggalt)
library(viridis)
library(e1071)
library(LiblineaR)
library(CORElearn)
library(grplasso)
library(randomForest)
library(deepboost)
library(stringr)
library(class)

if(!exists("knnClass"))
  knnClass = knn

library(infotheo)
library(praznik)
library(mRMRe)


source("phi_lin.R")
source("stab_vs_acc.R")


logistic_LASSO = function(x,y,nfeat=20,frac=0.2,lambda=0.1,lambda_i=5000,train = sample(1:nrow(x),replace=T,size=1*nrow(x)),test = setdiff(1:nrow(x),train), lambda_finals=c(0.1),alpha_finals=c(1), imp_method="linear"){
  
  feature_initial_indices = 1:ncol(x)
  norm_train = x[train,]
  norm_test = x[test,]
  
  no_norm_train = norm_train
  norm_train = scale(x[train,])
  norm_test = scale(x[test,],attr(norm_train,"scaled:center"),attr(norm_train,"scaled:scale"))
  
  pre_removed_features = which(is.na(norm_train[1,]))
  if(length(pre_removed_features)>0)
    feature_initial_indices = feature_initial_indices[-pre_removed_features]
  
  norm_train = norm_train[,feature_initial_indices]
  no_norm_train = no_norm_train[,feature_initial_indices]
  
  g = variances(no_norm_train,y[train])
  selected_features = order(g,decreasing=T)[1:min(lambda_i,length(feature_initial_indices))]
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]
  
  
  coeffs = glmnet(norm_train,y[train],family="binomial",alpha=alpha_finals[1],lambda=c(lambda))
  selected_features = which(as.vector(abs(coeffs$beta)) != 0)
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]
  
  model_beta = coeffs$beta[selected_features]
  beta = rep(0,length(model_beta))
  a0 = rep(0,1)
  for(i in 1:length(beta))
    beta[i] = model_beta[i]
  a0[1] = -coeffs$a0[1]
  eva = evaluate_model(beta,a0,norm_test[,feature_initial_indices],y[test],str="on OOB sample")
  
  list(list(alpha_final=alpha_finals[1],lambda_final=lambda_finals[1],feat=feature_initial_indices,acc=eva$acc,bcr=eva$bcr,beta=beta))
}

#for simulations only
logistic_gLASSO = function(x,y,nfeat=20,frac=0.2,lambda=0.1,lambda_i=20,train = sample(1:nrow(x),replace=T,size=1*nrow(x)),test = setdiff(1:nrow(x),train), lambda_finals=c(0.1),alpha_finals=c(0),imp_method="linear"){
  
  print("gLASSO")
  y[which(y==-1)] = 0
  feature_initial_indices = 1:ncol(x)
  norm_train = x[train,]
  norm_test = x[test,]
  
  no_norm_train = norm_train
  norm_train = scale(x[train,])
  norm_test = scale(x[test,],attr(norm_train,"scaled:center"),attr(norm_train,"scaled:scale"))
  
  pre_removed_features = which(is.na(norm_train[1,]))
  if(length(pre_removed_features)>0)
    feature_initial_indices = feature_initial_indices[-pre_removed_features]
  
  norm_train = norm_train[,feature_initial_indices]
  no_norm_train = no_norm_train[,feature_initial_indices]
  
  t = y[train]
  fit <<- grplasso(norm_train, t, clusters, c(lambda), model = LogReg(),weights=rep(1,length(t)),offset = rep(0,length(t)),standardize = F)
  selected_features = which(as.vector(abs(fit$coefficients)) != 0)
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]
  
  model_beta = fit$coefficients[selected_features]
  beta = rep(0,length(model_beta))
  a0 = rep(0,1)
  for(i in 1:length(beta))
    beta[i] = model_beta[i]
  a0[1] = 0
  y[which(y==0)] = -1
  eva = evaluate_model(beta,a0,norm_test[,feature_initial_indices],y[test],str="on OOB sample")
  
  list(list(alpha_final=alpha_finals[1],lambda_final=lambda_finals[1],feat=feature_initial_indices,acc=eva$acc,bcr=eva$bcr,beta=beta))
}

RF = function(x,y,nfeat=20,frac=0.2,lambda=0.1,lambda_i=20,train = sample(1:nrow(x),replace=T,size=1*nrow(x)),test = setdiff(1:nrow(x),train), lambda_finals=c(0.1),alpha_finals=c(0), imp_method="linear"){
  print("Random forest")
  
  feature_initial_indices = 1:ncol(x)
  norm_train = x[train,]
  norm_test = x[test,]
  
  no_norm_train = norm_train
  norm_train = scale(x[train,])
  norm_test = scale(x[test,],attr(norm_train,"scaled:center"),attr(norm_train,"scaled:scale"))
  
  pre_removed_features = which(is.na(norm_train[1,]))
  if(length(pre_removed_features)>0)
    feature_initial_indices = feature_initial_indices[-pre_removed_features]
  
  norm_train = norm_train[,feature_initial_indices]
  no_norm_train = no_norm_train[,feature_initial_indices]
  
  g = variances(no_norm_train,y[train])
  selected_features = order(g,decreasing=T)[1:min(lambda_i,length(feature_initial_indices))]
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]#feature_initial_indices]
  
  
  forest = randomForest(norm_train,as.factor(y[train]),importance = T, ntree=lambda)
  #to ensure that in case of ties, features are selected randomly
  varImp = forest$importance[,"MeanDecreaseAccuracy"] + runif(ncol(norm_train),-10^-12,10^-12)
  selected_features = order(varImp,decreasing = T)[1:nfeat]
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]
  
  
  return_list = list()
  #myCluster <- makeCluster(8, # number of cores to use
                      #     type = "FORK")
  #registerDoParallel(myCluster)
  return_list = foreach(lambda_final=lambda_finals)%do%{
    forest = randomForest(norm_train,as.factor(y[train]),importance = T,xtest=norm_test[,feature_initial_indices],ytest=as.factor(y[test]),keep.forest=T,ntree=lambda_final)
    predicted = forest$test$predicted
    
    acc = mean(predicted == as.factor(y[test]))
    
    predicted = as.numeric(as.vector(predict(forest,scale(x)[,feature_initial_indices])))
    if(imp_method == "linear")
      beta = get_importance_pred_change_v2(forest,scale(x)[,feature_initial_indices],y,1:nrow(x),predicted)
    else
      beta = length(selected_features)*abs(forest$importance[,"MeanDecreaseAccuracy"])/sum(abs(forest$importance[,"MeanDecreaseAccuracy"]))
    beta = sapply(beta, function(i) { max(i,10^-10) })
    
    list(alpha_final=alpha_finals[1],lambda_final=lambda_final,feat=feature_initial_indices,acc=acc,bcr=0,beta=beta)
  }
  return_list 
}

Relief = function(x,y,nfeat=20,frac=0.2,lambda=0.1,lambda_i=20,train = sample(1:nrow(x),replace=T,size=1*nrow(x)),test = setdiff(1:nrow(x),train), lambda_finals=c(0.1),alpha_finals=c(0),imp_method="linear"){
  print("Relief")
  feature_initial_indices = 1:ncol(x)
  norm_train = x[train,]
  norm_test = x[test,]
  
  no_norm_train = norm_train
  norm_train = scale(x[train,])
  norm_test = scale(x[test,],attr(norm_train,"scaled:center"),attr(norm_train,"scaled:scale"))
  
  pre_removed_features = which(is.na(norm_train[1,]))
  if(length(pre_removed_features)>0)
    feature_initial_indices = feature_initial_indices[-pre_removed_features]
  
  norm_train = norm_train[,feature_initial_indices]
  no_norm_train = no_norm_train[,feature_initial_indices]
  
  g = variances(no_norm_train,y[train])
  selected_features = order(g,decreasing=T)[1:min(lambda_i,length(feature_initial_indices))]
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]
  
  splits = str_split(lambda, pattern="-")
  meth = splits[[1]][1]
  k = as.numeric(splits[[1]][2])
  
  rel = relief(norm_train,y[train],meth,k=k)
  selected_features = order(rel,decreasing = T)[1:nfeat]
  beta = rel[selected_features]
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]
  
  nn = knnClass(norm_train,norm_test[,feature_initial_indices], y[train],k=k)
  eva = evaluate_model(xtest=norm_test[,feature_initial_indices],ytest=y[test],test_indices=test,predictionsTest=nn)
  
  if(imp_method == "linear"){
    customPredict = function(model, test_row){
        return(knnClass(norm_train,test_row, y[train],k=k))
    }
    predicted = knnClass(norm_train,scale(x)[,feature_initial_indices], y[train],k=k)
    beta = get_importance_pred_change_v2(NULL,scale(x)[,feature_initial_indices],y,1:nrow(x),predicted,customPredict)
  }
  list(list(alpha_final=alpha_finals[1],lambda_final=lambda_finals[1],feat=feature_initial_indices,acc=eva$acc,bcr=eva$bcr,beta=beta))
}


logistic_RFE_multiples = function(x,y,nfeat=50,frac=0.2,lambda=0.1,lambda_i=50000,train = sample(1:nrow(x),replace=T,size=1*nrow(x)),test = setdiff(1:nrow(x),train), lambda_finals=c(0.1),alpha_finals=c(0),imp_method="linear-logistic"){
  
  feature_initial_indices = 1:ncol(x)
  norm_train = x[train,]
  norm_test = x[test,]
  y_train = y[train]
  y_test = y[test]
  
  no_norm_train = norm_train
  norm_train = scale(norm_train)
  norm_test = scale(norm_test,attr(norm_train,"scaled:center"),attr(norm_train,"scaled:scale"))
  
  pre_removed_features = which(is.na(norm_train[1,]))
  if(length(pre_removed_features)>0)
    feature_initial_indices = feature_initial_indices[-pre_removed_features]
  
  norm_train = norm_train[,feature_initial_indices]
  no_norm_train = no_norm_train[,feature_initial_indices]
  
  lambda_i = min(lambda_i,length(feature_initial_indices))
  
  g = variances(no_norm_train,y_train)
  selected_features = order(g,decreasing=T)[1:min(lambda_i,length(feature_initial_indices))]
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]#feature_initial_indices]
  
  while(ncol(norm_train)>nfeat){
    if(imp_method == "linear-logistic"){
      coeffs = glmnet(norm_train,y_train,family="binomial",alpha=0,lambda=c(lambda),thresh=10^-5)
    }else{
      coeffs = svm(norm_train,y[train],cost=1/lambda,kernel="linear")
      coeffs$beta = (t(coeffs$coefs) %*% coeffs$SV)
    }
    selected_features = order(as.vector(abs(coeffs$beta)),decreasing = T)[1:max(ncol(norm_train)*(1-frac),nfeat)]
    norm_train = norm_train[,selected_features]
    feature_initial_indices = feature_initial_indices[selected_features]
  }
  
  res = list()
  for(lambda_final in lambda_finals){
    
      if(imp_method=="linear-logistic"){
        model = glmnet(norm_train,y_train,family="binomial",alpha=0,lambda=c(lambda_final))
      }else{
        model = svm(norm_train,y[train],cost=1/lambda_final,kernel="linear")
        model$beta = (t(model$coefs) %*% model$SV)
        model$a0 = c(model$rho)
      }
      
      a0 = rep(0,1)
      beta = rep(0,length(model$beta))
      for(i in 1:length(beta))
        beta[i] = model$beta[i]
      a0[1] = -model$a0[1]
      
      eva = evaluate_model(beta,a0,norm_test[,feature_initial_indices],y_test,str="on OOB sample",test)
      res[[length(res)+1]] = list(alpha_final=0,lambda_final=lambda_final,feat=feature_initial_indices,acc=eva$acc,bcr=eva$bcr,beta=beta)
  }
  res
}

filterThenLogReg = function(x,y,nfeat=50,frac=0.2,lambda=0.1,lambda_i=50000,train = sample(1:nrow(x),replace=T,size=1*nrow(x)),test = setdiff(1:nrow(x),train), lambda_finals=c(0.1),alpha_finals=c(0),imp_method="linear"){
  
  feature_initial_indices = 1:ncol(x)
  norm_train = x[train,]
  norm_test = x[test,]
  y_train = y[train]
  y_test = y[test]
  
  no_norm_train = norm_train
  norm_train = scale(norm_train)
  norm_test = scale(norm_test,attr(norm_train,"scaled:center"),attr(norm_train,"scaled:scale"))
  
  pre_removed_features = which(is.na(norm_train[1,]))
  if(length(pre_removed_features)>0)
    feature_initial_indices = feature_initial_indices[-pre_removed_features]
  
  norm_train = norm_train[,feature_initial_indices]
  no_norm_train = no_norm_train[,feature_initial_indices]
  
  fil = get(METHODS[lambda])(x=norm_train,y=y[train],nfeat=nfeat)
  selected_features = order(fil,decreasing = T)[1:nfeat]
  
  feature_initial_indices = feature_initial_indices[selected_features]
  norm_train = norm_train[,selected_features]
  
  res = list()
  
  for(lambda_final in lambda_finals){
    model = glmnet(norm_train,y_train,family="binomial",alpha=0,lambda=lambda_final)
    beta = rep(0,length(model$beta))
    a0 = rep(0,1)
    for(i in 1:length(beta))
      beta[i] = model$beta[i]
    a0[1] = -model$a0[1]
    eva = evaluate_model(beta,a0,norm_test[,feature_initial_indices],y_test,str="on OOB sample",test)
    imps = rep(1,length(feature_initial_indices))
    res[[length(res)+1]] = list(alpha_final=0,lambda_final=lambda_final,feat=feature_initial_indices,acc=eva$acc,bcr=eva$bcr,beta=beta)
}
  res
}

runEnsemble_RFE_multiple = function(x,y,Ms=1:10,nfeat=50,frac=0.2,lambda=0.2,alpha_finals=c(1),lambda_finals=c(0.1),lambda_i=5000,method=logistic_RFE_multiples,imp_method="linear",stabs=TRUE){
  
  M = length(Ms)
  d = ncol(x)
  subsets = matrix(data=0,nrow=M,ncol=d)
  coeffs = matrix(data=0,nrow=M,ncol=d)
  accuracies = rep(0,M)
  bcrs = rep(0,M)
  minM = Ms[1]
  
  accuracies = list()
  bcrs = list()
  phi_msis = list()
  coeffs = list()
  for(lambda_final in lambda_finals){
    for(alpha_final in alpha_finals){
      accuracies[[length(accuracies)+1]] = 0
      bcrs[[length(bcrs)+1]] = 0
      coeffs[[length(coeffs)+1]] = matrix(data=0,nrow=M,ncol=d)
    }
  }
  NCLUSTER = 1
  myCluster <- makeCluster(NCLUSTER, type = "FORK")
  registerDoParallel(myCluster)
  clusters = rep(0,length(Ms))
  for(i in 1:(length(clusters))){
    clusters[i] = ((i-1)%%NCLUSTER) + 1
  }
  
  resAll = foreach(ind=1:min(NCLUSTER,length(Ms)))%do%{
    
    subMs = Ms[which(clusters==ind)]
    subsets2 = matrix(data=0,nrow=length(subMs),ncol=d)
    accuracies2 = list()
    bcrs2 = list()
    coeffs2 = list()
    
    for(lambda_final in lambda_finals){
      for(alpha_final in alpha_finals){
        accuracies2[[length(accuracies2)+1]] = 0
        bcrs2[[length(bcrs2)+1]] = 0
        coeffs2[[length(coeffs2)+1]] = matrix(data=0,nrow=length(subMs),ncol=d)
      }
    }
    
    for(mind in 1:length(subMs)){
      m = subMs[mind]
      set.seed(m)
      res = method(x,y,nfeat=nfeat,frac=frac,lambda=lambda,lambda_i=lambda_i,lambda_finals = lambda_finals,alpha_finals=alpha_finals,imp_method=imp_method)
      for(i in 1:length(res)){
        accuracies2[[i]] = accuracies2[[i]] + res[[i]][["acc"]]/length(subMs)
        bcrs2[[i]] = bcrs2[[i]] + res[[i]][["bcr"]]/length(subMs)
        coeffs2[[i]][mind,res[[1]]$feat] = res[[i]][["beta"]]
      }
      subsets2[mind,res[[1]]$feat] = T
    }
    list(subsets=subsets2,coeffs=coeffs2,accuracies=accuracies2,bcrs=bcrs2)
  }
  stopCluster(myCluster)
  
  Mcounter = 1
  for(i in 1:length(resAll)){
    Mplus = nrow(resAll[[i]]$subsets)
    subsets[Mcounter:(Mcounter+Mplus-1),] = resAll[[i]]$subsets
    for(j in 1:length(resAll[[i]]$coeffs)){
      coeffs[[j]][Mcounter:(Mcounter+Mplus-1),] = resAll[[i]]$coeffs[[j]]
      accuracies[[j]][Mcounter:(Mcounter+Mplus-1)] = resAll[[i]]$accuracies[[j]]
      bcrs[[j]][Mcounter:(Mcounter+Mplus-1)] = resAll[[i]]$bcrs[[j]]
    }
    Mcounter = Mcounter + Mplus
  }
  
  simMatrix = abs(cor(x,method="spearman"))
  simMatrix[which(is.na(simMatrix))] = 0
  
  for(i in 1:length(accuracies)){
      if(stabs){
        phi_msis[[i]] = phi_lin(coeffs[[i]],simMatrix)
      }else{
        phi_msis[[i]] = 0
      }
  }
  
  phi_S = 0
  POGR = 0
  phi = 0
  if(stabs){
    phi_S = phi_C(subsets,simMatrix)
    POGR = POGr(subsets,simMatrix)
    phi = stability(subsets,d)
  }
  
  list(acc=accuracies,bcr=bcrs,stab=phi,phi_msis=phi_msis, phi_S=phi_S, bcrs=bcrs,subsets=subsets, coeffs=coeffs, POGR=POGR)
  
}


compute_frame_multipleks = function(x,y,Ms,lambdas=0.2,alpha_finals=c(0),lambda_finals=c(0.1),frac=0.2,nfeats=c(20),method=logistic_RFE_multiples,imp_method="linear",dfil=5000){
  i = 1
  frame <- foreach(nfeat=nfeats, .combine=rbind)%do%{
    foreach(lambda=lambdas, .combine=rbind)%do%{
      print(paste0("------", as.character(i), "/",as.character(length(lambdas))))
      i = i + 1
      compute_frame_multiple(x,y,Ms,lambda,alpha_finals,lambda_finals,frac,nfeat,method,imp_method,dfil=dfil)$frame
    }
  }
  frame
}

compute_frame_multiple = function(x,y,Ms,lambda=0.2,alpha_finals=c(0),lambda_finals=c(0.1),frac=0.2,nfeat=20,method=logistic_RFE_multiples,imp_method="linear",dfil=5000,stabs=TRUE){
  
  M = length(Ms)
  res = runEnsemble_RFE_multiple(x,y,Ms=Ms,lambda=lambda,frac=frac,nfeat=nfeat,alpha_finals=alpha_finals,lambda_finals=lambda_finals,lambda_i=dfil,method=method,imp_method=imp_method,stabs=stabs)
  
  frame = data.frame(matrix(nrow=0,ncol=17))
  i = 1
  nfeat = sum(res$coeffs[[i]]!=0)/length(Ms)
  for(lambda_final in lambda_finals){
    for(alpha_final in alpha_finals){
      frame = rbind(frame,c(mean(res[["acc"]][[i]]),mean(res[["bcr"]][[i]]),res$stab,res$phi_msis[[i]],res$phi_S,res$POGR,lambda,lambda_final,alpha_final,dfil,nfeat))
      i = i+1
    }
  }
  colnames(frame) <- c("accuracy","bcr","phi","phi_msi","phi_S","POGR","lambda","lambda_f","alpha_f","dfil","k")
  list(frame=frame,sub=res$subsets,coeffs=res$coeffs)
}

#CODE TO GENERATE THE RESULTS OF SECTION 6.2

logrfe_frame <<- list()
logrfes = function(datas,M){
  lambdas = c(1)
  lambda_finals = c(0.01)
  for(i in 1:length(datas)){
    print(paste0(datas[[i]]$name," " ,as.character(i), "/", as.character(length(datas))))
    nfeats=c(min(20,round(sqrt(ncol(datas[[i]]$x)))))
    logrfe_frame[[datas[[i]]$name]] <<- compute_frame_multipleks(datas[[i]]$x,datas[[i]]$y,1:M,lambdas,method=logistic_RFE_multiples,nfeats=nfeats,lambda_finals = lambda_finals, imp_method = "linear-logistic")
  }
  logrfe_frame
}

svmrfe_frame <<- list()
svmrfes = function(datas,M){
  lambdas = c(1000)
  lambda_finals = c(300)
  for(i in 1:length(datas)){
    print(paste0(datas[[i]]$name," " ,as.character(i), "/", as.character(length(datas))))
    nfeats=c(min(20,round(sqrt(ncol(datas[[i]]$x)))))
    svmrfe_frame[[datas[[i]]$name]] <<- compute_frame_multipleks(datas[[i]]$x,datas[[i]]$y,1:M,lambdas,method=logistic_RFE_multiples,nfeats=nfeats,lambda_finals = lambda_finals, imp_method = "linear-svm")
  }
  svmrfe_frame
}

lasso_frame <<- list()
lassos = function(datas,M, alphas=c(1,0.8)){
  
  for(i in 1:length(datas)){
    print(datas[[i]]$name)
    
    for(alpha in alphas){
      print(paste0("alpha: ", toString(alpha)))
      lambda_min = 10^-8
      lambda_max = 1
      
      nfeat = 0
      target_nfeat = min(20, round(sqrt(ncol(datas[[i]]$x))))
      lambda = 0
      while(abs(nfeat - target_nfeat) > 0.1){
        lambda = (lambda_min+lambda_max)/2
        print(lambda)
        tryCatch({
          lasso_fr = compute_frame_multiple(datas[[i]]$x,datas[[i]]$y,Ms=1:M,lambda=lambda,method=logistic_LASSO,alpha_finals=c(alpha),stabs=FALSE)
          nfeat = lasso_fr$frame[,"k"]
        },
        error = function(e) {  }
        )
        print(nfeat)
        if(abs(nfeat - target_nfeat) > 0.1){#} && (lambda_max-lambda_min) > 10^-7){
          if(nfeat > target_nfeat) lambda_min = lambda
          else lambda_max = lambda
        }
      }
      f = compute_frame_multiple(datas[[i]]$x,datas[[i]]$y,Ms=1:M,lambda=lambda,method=logistic_LASSO,alpha_finals=c(alpha))
      if(alpha == alphas[1]){
          lasso_frame[[datas[[i]]$name]] <<- f$frame
      }else{
          lasso_frame[[datas[[i]]$name]] <<- rbind(lasso_frame[[datas[[i]]$name]], f$frame)
        }
      }
  }
  lasso_frame
}

rf_frame  <<- list()

forests = function(datas, M){
  lambdas = c(1000)
  lambda_finals = c(1000)
  for(i in 1:length(datas)){
    print(datas[[i]]$name)
    #for reproducibility, problem occured for some runs in range 30:100
    Ms = 1:M
    if(M>30) Ms = c(1:30,101:(100+M-30))
    rf_frame[[datas[[i]]$name]] <<- compute_frame_multipleks(datas[[i]]$x,datas[[i]]$y,1:M,lambdas,method=RF,nfeat=min(20,round(sqrt(ncol(datas[[i]]$x)))),lambda_finals = lambda_finals,imp_method = "linear")
    save(rf_frame,file=paste0("resultsRF-", datas[[i]]$name, ".RData"))
  }
  save(rf_frame,file=paste0("resultsRF-", datas[[i]]$name, ".RData"))
  rf_frame
}

relief_frame  <<- list()

reliefs = function(datas, M){
  lambda_finals = c(1)
  lambdas = c("ReliefFequalK-5")
  for(i in 1:length(datas)){
    corrs = datas[[i]]$x
    print(datas[[i]]$name)
    relief_frame[[datas[[i]]$name]] <<- compute_frame_multipleks(datas[[i]]$x,datas[[i]]$y,1:M,lambdas,method=Relief,nfeat=min(20,round(sqrt(ncol(datas[[i]]$x)))),imp_method = "linear")
  }
  relief_frame
}


filter_frame  <<- list()

METHODS <<- c("t_test","MIM","mrmr")

filters = function(datas, M, method="t_test"){
  lambdas = c(which(METHODS==method))
  lambda_finals = c(0.01)
  for(i in 1:length(datas)){
    print(datas[[i]]$name)
    filter_frame[[datas[[i]]$name]] <<- compute_frame_multipleks(datas[[i]]$x,datas[[i]]$y,1:M,lambdas,lambda_finals=lambda_finals,method=filterThenLogReg,nfeat=min(20,round(sqrt(ncol(datas[[i]]$x)))),imp_method = "linear")
  }
  filter_frame
}

#for paper results
#takes more than 10 hours
paper = function(){
  logrfe = logrfes(datasetsECML,100)
  svmrfe = svmrfes(datasetsECML,100)
  rf = forests(datasetsECML,100)
  rel = reliefs(datasetsECML,100)
  mr=filters(datasetsECML,100,"mrmr")
  mi=filters(datasetsECML,100,"MIM")
  tt=filters(datasetsECML,100,"t_test")
  
  plot_multiple_mci(list(logrfe,svmrfe,rf,rel,mr,mi,tt),c("phi"),c("logrfe","svmrfe","rf","relief","mrmr","mim","ttest"), list(alon), "paper-phi.pdf")
  plot_multiple_mci(list(logrfe,svmrfe,rf,rel,mr,mi,tt),c("phi_msi"),c("logrfe","svmrfe","rf","relief","mrmr","mim","ttest"), list(alon), "paper-phi_msi.pdf")
}


#for demo
#takes +- 2 minutes
demo = function(){
  alon = datasetsECML[[1]]
  logrfe=logrfes(list(alon),5)
  svmrfe=svmrfes(list(alon),5)
  rf=forests(list(alon),5)
  rel=reliefs(list(alon),5)
  mr=filters(list(alon),5,"mrmr")
  mi=filters(list(alon),5,"MIM")
  tt=filters(list(alon),5,"t_test")

  plot_multiple_mci(list(logrfe,svmrfe,rf,rel,mr,mi,tt),c("phi"),c("logrfe","svmrfe","rf","relief","mrmr","mim","ttest"), list(alon), "test-phi.pdf")
  plot_multiple_mci(list(logrfe,svmrfe,rf,rel,mr,mi,tt),c("phi_msi"),c("logrfe","svmrfe","rf","relief","mrmr","mim","ttest"), list(alon), "test-phi_msi.pdf")
}
