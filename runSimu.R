library(MASS)
source("runXP.R")

build_simu_dataset = function(MU1, INTRA_GROUP_COR, SEED, REL_SIZE){
  set.seed(SEED)
  mu1 = MU1
  mu2 = 0.05
  N=100
  
  NGROUPS = 1000
  NRELEVANT_GROUPS = 5
  GROUP_SIZES = c(rep(REL_SIZE,NRELEVANT_GROUPS),rep(1,NGROUPS-NRELEVANT_GROUPS))
  d = sum(GROUP_SIZES)
  
  
  RELEVANT_GROUPS = 1:NRELEVANT_GROUPS
  drel = sum(GROUP_SIZES[RELEVANT_GROUPS])
  
  #INTRA_GROUP_COR = 0.95
  INTER_GROUP_COR = 0
  
  get_group_index = function(f){
    
    if(f <= GROUP_SIZES[1]){
      return(1)
    }
    current_sum = GROUP_SIZES[1]
    for(i in 2:(NGROUPS-1)){
      current_sum = current_sum + GROUP_SIZES[i]
      if(f <= current_sum){
        return(i)
      }
    }
    
    return(NGROUPS)
  }
  
  clusters = rep(0,d)
  for(i in 1:d){
    clusters[i] = get_group_index(i)
  }
  
  covariance = matrix(data=0,nrow=d,ncol=d)
  
  for(i in 1:nrow(covariance)){
    for(j in 1:ncol(covariance)){
      if(i == j){
        covariance[i,j] = 1
      }else if(clusters[i] == clusters[j]){
        covariance[i,j] = INTRA_GROUP_COR
      }else{
        covariance[i,j] = INTER_GROUP_COR
      }
    }
  }
  
  dataPos = mvrnorm(n=N,mu=c(rep(mu1,drel),rep(mu2,d-drel)),Sigma=covariance)
  dataNeg = mvrnorm(n=N,mu=c(rep(-mu1,drel),rep(-mu2,d-drel)),Sigma=covariance)
  data = rbind(dataPos,dataNeg)
  labels = c(rep(1,N),rep(-1,N))
  list("labels"=labels,"data"=data,"clusters"=clusters)
}

runSimu = function(index){
  
  INTRA_GROUP_COR = 0.8
  MU1 = 0.5
  if(index == 1){
    MU1 = 0.35
  }
  if(index == 3){
    INTRA_GROUP_COR = 0.95
  }
  
  ms = seq(1,15,by=1)
  
  phis = rep(0,length(ms))
  phis_msi = rep(0,length(ms))
  phis_pears = rep(0,length(ms))
  phis_S = rep(0,length(ms))
  phis_POGR = rep(0,length(ms))
  
  ks = rep(0,length(ms))
  ks_out = rep(0,length(ms))
  ds = list()
  M=30
  ind = 1
  
  for(REL_SIZE in ms){
    print("SIZE:")
    print(REL_SIZE)
    NSEED = 10
    for(SEED in 1:NSEED){
      r = build_simu_dataset(MU1, INTRA_GROUP_COR, SEED, REL_SIZE)
      
      clusters <<- r$clusters
      if(index == 3)
        d = compute_frame_multiple(r$data,r$labels,Ms=1:M,lambda=0.1,alpha_finals = c(1),lambda_finals=c(1),method = logistic_LASSO)
      else
        d = compute_frame_multiple(r$data,r$labels,Ms=1:M,lambda=10,alpha_finals = c(0),lambda_finals=c(1),method = logistic_gLASSO)
      
      phis_msi[ind] =  phis_msi[ind] + d$frame[1,"phi_msi"]/NSEED#d$frame[1,"iw_stab"]
      phis_S[ind] = phis_S[ind] + d$frame[1,"phi_S"]/NSEED
      phis[ind] = phis[ind] + d$frame[1,"phi"]/NSEED
      phis_pears[ind] = phis_pears[ind] + phi_pears(d$coeffs[[1]])/NSEED 
      phis_POGR[ind] = phis_POGR[ind] + d$frame[1,"POGR"]/NSEED
    }
    
    plot_simu(ms, ind, phis, phis_pears, phis_S, phis_msi, phis_POGR, index)
    ind=ind+1
  }
}

plot_simu = function(ms, ind, phis, phis_pears, phis_S, phis_msi, phis_POGR, index){
  library(ggplot2)
  xs = ms[1:ind]
  frame = matrix(data=0,nrow=5*ind,ncol=3)
  frame[,1] = c(xs,xs,xs,xs,xs)
  frame[,2] = c(phis[1:ind],phis_pears[1:ind],phis_S[1:ind], phis_msi[1:ind], phis_POGR[1:ind])
  frame[,3] = c(rep("phi",ind),rep("phi_pears",ind),rep("phi_S",ind),rep("phi_msi",ind), rep("POGR",ind))
  
  frame = as.data.frame(frame)
  
  colnames(frame) = c("q","Stability","Type")
  frame[,1] = as.numeric(as.character(frame[,1]))
  frame[,2] = as.numeric(as.character(frame[,2]))
  
  p <- ggplot(data=frame,mapping=aes(x=q,y=Stability,group=Type,col=Type)) + theme_grey(base_size = 28,base_rect_size = 3)
  p <- p + geom_point(size=4) + geom_line()
  plot(p)
  
  pdf(paste0("simu",toString(index),".pdf"), width = 48.25/1.3 / 2.54, height = 20/1.3 / 2.54)#width = 35/1.3 / 2.54, height = 20/1.3 / 2.54)#width = 48.25/1.3 / 2.54, height = 20/1.3 / 2.54)
  #colormodel = "grey")
  plot(p)
  dev.off()
}
