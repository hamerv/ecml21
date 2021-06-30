library(glmnet)
library(doParallel)
#library(datamicroarray)
library(ggalt)
library(viridis)
library(latex2exp)

POGr = function(coeffs,C,thresh=0.5){
  
  score = 0
  M = nrow(coeffs)
  C = C>=thresh
  
  rands = list()
  
  for(i in 1:M){
    fi = which(coeffs[i,]!=0)
    for(j in 1:M){
      fj = which(coeffs[j,]!=0)
      if(i!=j){
        r = 0
        for(f in fi){
          if(sum(C[fj,f])>0){
            r = r + 1
          }
        }
        score = score + (r)/(length(fi)) 
      }
    }
  }
  return(score/(M*(M-1)))
}

phi_C = function(coeffs,C){
  
  d = ncol(coeffs)
  M = nrow(coeffs)
  k = sum(coeffs!=0)/M
  
  coeffs = abs(coeffs) != 0
  S = matrix(data=0,nrow=d,ncol=d)
  pfs = apply(coeffs,2,mean)
  pfs_non_zero = which(pfs!=0)
  ks = apply(coeffs,1,sum)
  
  kC = 0
  for(f in pfs_non_zero){
    for(f2 in pfs_non_zero){
      pff2 = mean((coeffs[,f]*coeffs[,f2])!=0)
      S[f,f2] = (M/(M-1)) * (pff2 - pfs[f]*pfs[f2]  )
      kC = kC + pff2*C[f,f2]
    }
  }

  num = sum(diag(C%*%S))
  diago = (1-(k/d))*k/d
  offdiago = (mean(ks^2) - k)/(d^2-d) - (k^2)/(d^2)
  S0 = matrix(data=offdiago,nrow=d,ncol=d)
  for(i in 1:d){
    S0[i,i] = diago
  }
  denom = sum(diag(C%*%S0))
  phi = 1- (num/denom)
  return(phi)
}

get_importance_pred_change_v2 = function(model,data,labels,test,predictions=NULL, customPredict=NULL){
  
  if(is.null(predictions))
    predictions = as.numeric(as.vector(predict(model,data)))
  
  if(!is.null(customPredict))
    predict = customPredict
  
  print("V2")
  importances = rep(1,ncol(data))
  #change_table = seq(0,100,by=0.01)
  tol = 0.001
  
  changes = matrix(data=0,nrow=length(test),ncol=ncol(data))
  
  for(tt in 1:length(test)){
    test_index = test[tt]
    for(f in 1:ncol(data)){
      changes_min = 0
      changes_max = 100
      
      ch = 1
      while((changes_max-changes_min) > tol){
        
        switch = F
        test_ex = data[test_index,]
        #print(test_ex)
        test_ex[f] = test_ex[f] + (changes_max+changes_min)/2
        new_pred = as.numeric(as.vector(predict(model,matrix(test_ex,nrow=1))))
        if(new_pred != predictions[test_index])
          switch = T
        
        test_ex[f] = test_ex[f] - 2*(changes_max+changes_min)/2
        new_pred = as.numeric(as.vector(predict(model,matrix(test_ex,nrow=1))))
        if(new_pred != predictions[test_index])
          switch = T
        
        
        if(switch){
          changes_max = (changes_max + changes_min)/2
        }else{
          changes_min = (changes_max + changes_min)/2
        }
      }
      changes[tt,f] = 1/abs((changes_max+changes_min)/2)
    }
    changes[tt,] = changes[tt,]*ncol(data)/sum(changes[tt,])
    #changes = changes/length(test)
    #print(importances)
  }
  for(f in 1:ncol(data)){
    importances[f] = mean(changes[,f])
  }
  #importances = ncol(data)*importances/sum(importances)
  return(importances)
}

evaluate_model = function(w,b,xtest,ytest,str="test set:", test_indices=1:length(ytest),predictionsTest=NULL){
  if(is.null(predictionsTest)){
    predictionsTest = c()
    for(i in 1:nrow(xtest)){
      predictionsTest[i] = (2*((xtest[i,]%*%w - b)>0)-1)
    }
  }
  accuracy = sum(predictionsTest==ytest)/length(ytest)
  errors = which(predictionsTest!=ytest)
  pos_labels = which(predictionsTest == 1)
  
  TN = sum(predictionsTest==-1 & ytest==-1)
  TP = sum(predictionsTest==1 & ytest==1)
  FN = sum(predictionsTest==-1 & ytest==1)
  FP = sum(predictionsTest==1 & ytest==-1)
  bcr = 0.5*(TP/(TP+FN) + TN/(FP+TN))
  
  list(acc=accuracy,errors=test_indices[errors],bcr=bcr,pos_labels=test_indices[pos_labels])
}


stability = function(subsets,d){
  subsets = subsets != 0
  phi = 0
  M = nrow(subsets)
  k = sum(subsets)/(d*M)
  for(i in 1:ncol(subsets)){
    p_i = mean(subsets[,i])
    phi = phi + (M/(M-1))*p_i*(1-p_i)
  }
  phi = phi/d
  phi = phi/(k*(1-k))
  return(1-phi)
}

phi_pears = function(coeffs){
  sims = 0
  M = nrow(coeffs)
  coeffs = abs(coeffs)
  for(i in 1:M){
    for(j in 1:M){
      if(i!=j){
          sims = sims + cor(coeffs[i,],coeffs[j,]) 
      }
    }
  }
  sims/(M*(M-1))
} 


get_pareto_points = function(frame,ini="def",measure="stability"){
  init = which(frame[,"accuracy"] == max(frame[,"accuracy"]))[[1]]
  p = frame[1,]
  pareto = data.frame(matrix(nrow=0,ncol=ncol(frame)))
  colnames(pareto) = colnames(frame)
  pareto[1,] = frame[init,]
  frame = frame[-init,]
  frame = as.data.frame(frame)
  while(is.data.frame(frame) && nrow(frame)>0){
    candidate = which(frame[,"accuracy"]==max(frame[,"accuracy"]))
    if(frame[candidate,measure]>pareto[nrow(pareto),measure]){
      pareto[nrow(pareto)+1,] = frame[candidate,]
    }
    frame = frame[-candidate,]
  }
  pareto
}


t_test = function(x,y,nfeat){
  pos = which(y==1)
  neg = which(y==-1)
  ts = rep(0,ncol(x))
  nas = c()
  for(i in 1:length(ts)){
    ts[i] = abs(mean(x[pos,i]) - mean(x[neg,i]))/(sqrt(var(x[pos,i])/length(pos) + var(x[neg,i])/length(neg)))
    if(is.na(ts[i])){
      ts[i] = 10^(-10)
      nas[length(nas)+1] = i 
    }
    if(is.infinite(ts[i])){
      nas[length(nas)+1] = i
    }  
  }
  if(length(nas)==0){
    ts/mean(ts)
  }else{
    ts/mean(ts[-nas])
  }
}

relief = function(x,y,method="ReliefFbestK",iter=200,k=3){
  print(method)
  merged = as.data.frame(cbind(x,y))
  colnames(merged)[length(colnames(merged))] = "label"
  abs(attrEval(label ~ ., merged, estimator=method, ReliefIterations=iter, kNearestEqual=k))
}


variances = function(x,y){
  apply(X=x,MARGIN=2,FUN=var)
}

MIM = function(x,y,nfeat){
  dis = discretize(x)
  sapply(X=1:ncol(x),FUN=function(i){
    mutinformation(dis[,i],y)
  })
}

mrmr = function(x,y,nfeat){
  scores = rep(0,ncol(x))
  data <- mRMR.data(data = data.frame(target=y, x))
  r = mRMR.classic(data=data, target_indices=1, feature_count=nfeat)
  scores[as.vector(r@filters$"1")-1] = 1
  scores
}

draw_map = function(m,groups=NULL,name="map.pdf", file=T){
  library(RColorBrewer)
  k_max = 0
  kbar = sum(m!=0)/nrow(m)
  for(i in 1:nrow(m)){
    k_i = sum(m[i,]!=0)
    k_max = max(kbar,k_max)
    m[i,] = 100*abs(m[i,])/sum(abs(m[i,which(m[i,]!=0)]))
  }
  d = ncol(m)
  occurences = rep(0,d)
  importances = rep(0,d)
  variances = rep(0,d)
  for(f in 1:d){
    occurences[f] = sum(m[,f]!=0)
    importances[f] = mean(abs(m[which(m[,f]!=0),f]))
    variances[f] = sqrt(var(m[which(m[,f]!=0),f]))/mean(m[which(m[,f]!=0),f])
  }
  
  
  vals = occurences-0.00000001*variances
  f_sorted = order(vals,decreasing = T)#[c(2,1,3:length(occurences))]
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if(file)
    pdf(name, width = 30/1.3 / 2.54, height = 20/1.3 / 2.54)#35
  plot(c(0, 100), c(0, nrow(m)), type = "n", xlab = "Cumulative importance (%)", ylab = "i",
       main = "Feature map",cex.lab=2, cex.axis=2, cex.main=2)
  m = m*1
  sum_group = 0
  for(l in 1:nrow(m)){
    k_l = sum(m[l,])
    #print(k_l)
    previous_sum = 0#0.7*5/9 * 100
    group_index = 1
    for(i in 1:d){
      f = f_sorted[i]
      if(!is.null(groups) && i>1 && groups[f_sorted[i]]==groups[f_sorted[i-1]]){
        #print(i)
        i = group_index
      }else if(i>1){
        group_index = group_index + 1
      }
      i = group_index + 0
      if(i <= 5){
        sum_group = sum_group + m[l,f]
      }
      rect(previous_sum, (l-1), previous_sum+m[l,f], l, col = col_vector[i%%length(col_vector)+0])
      previous_sum = previous_sum + m[l,f]
    }
  }
  if(file)
    dev.off()
}


plot_multiple_mci = function(results, measures=c("phi"), names=dnames, datasets=datasetsECML, file="file.pdf"){
  returns = list()
  pdf(file, width = 35/1.3 / 2.54, height = 20/1.3 / 2.54, #50/1.3
  )#colormodel = "grey")
  for(dataset in datasets){
    frame = as.data.frame(matrix(data=0,nrow=0,ncol=ncol(results[[1]][[dataset$name]]+3)))
    
    for(measure in measures){
      paretos = list()
      frameAll = as.data.frame(matrix(data=0,nrow=0,ncol=ncol(results[[1]][[dataset$name]])+1))
      colnames(frameAll) = c(colnames(results[[1]][[dataset$name]]),"groups")
      for(i in 1:length(results)){
        if(!is.null(results[[i]][[dataset$name]])){
          #results[[i]][[dataset$name]] = results[[i]][[dataset$name]][complete.cases(results[[i]][[dataset$name]]),]
          results[[i]][[dataset$name]][,"accuracy"] = as.numeric(as.character(results[[i]][[dataset$name]][,"accuracy"]))
          results[[i]][[dataset$name]][,measure] = as.numeric(as.character(results[[i]][[dataset$name]][,measure]))
          print(results[[i]][[dataset$name]])
          paretos[[i]] = get_pareto_points(results[[i]][[dataset$name]], measure=measure)
          paretos[[i]] = cbind(paretos[[i]], rep(paste(names[[i]],dataset$name,measure,sep="-"),nrow(paretos[[i]])))
          paretos[[i]] = cbind(rep(names[[i]],nrow(paretos[[i]])),paretos[[i]])
          frameAll = rbind(frameAll,paretos[[i]])
        }
      }
      frameAll = cbind(frameAll,frameAll[,measure],rep(measure,nrow(frameAll)))#,rep(paste(dataset$name,measure,sep="-"),nrow(frameAll)))
      colnames(frameAll)[length(colnames(frameAll))-1] = "measure"
      colnames(frameAll)[length(colnames(frameAll))] = "type"
      colnames(frameAll)[length(colnames(frameAll))-2] = "groups"
      frameAll[,"type"] = as.character(frameAll[,"type"])
      frameAll[,"groups"] = as.character(frameAll[,"groups"])
      
      frameAll = cbind(frameAll,rep("no",nrow(frameAll)))
      colnames(frameAll)[length(colnames(frameAll))] = "pareto"
      frameAll[,"pareto"] = as.character(frameAll[,"pareto"])
      for(i in 1:nrow(frameAll)){
        w = which(frameAll[,"measure"]>frameAll[i,"measure"] & frameAll[,"accuracy"]>frameAll[i,"accuracy"])
        if(length(w)==0) frameAll[i,"pareto"] = "yes"
      }
      
      frame = rbind(frame,frameAll)
    }
    
    frame[,"type"] = as.character(frame[,"type"])
    frame[which(frame[,"type"]=="lin_stab"),"type"] = "phi_pears"
    frame[which(frame[,"type"]=="int_stab"),"type"] = "phi_cic"
    frame[which(frame[,"type"]=="stability"),"type"] = "phi"
    colnames(frame)[1] = "method"
    
    p <- ggplot(data=frame,mapping=aes(x=measure,y=accuracy, group = method, shape=type, col=method)) + theme_grey(base_size = 30,base_rect_size = 3)
    p <- p + geom_point(size=10,shape=17) + scale_shape(solid = T) #+ geom_line(mapping=aes(x=measure,y=accuracy, group = groups),linetype="dotted",size=0.75) #+ xlim(min(frameAll[,measure])-0.05,max(frameAll[,measure])+0.05) + ylim(min(frameAll[,"accuracy"])-0.05,max(frameAll[,"accuracy"])+0.05)
    p <- p + xlab("Stability") + ylab("Accuracy") + labs(title=dataset$name)
    
    p <- p + geom_line(data=frame[which(frame[,'pareto']=="yes"),], mapping=aes(x=measure,y=accuracy, group=type, color=type), linetype="dotted", color="black", size=2)
    p <- p + geom_text(aes(label=method),hjust=-0.2, vjust=-0.2, size=10)
    p <- p + ylim(min(frame[,"accuracy"])-0.005,max(frame[,"accuracy"])+0.005) + xlim(min(frame[,"measure"])-0.02,max(frame[,"measure"])+0.08)
    
    plot(p)
    
    returns[[dataset$name]] = frame
  }
  dev.off()
  return(returns)
}