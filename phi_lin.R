library(lpSolve)

sim_lin = function(feats1,feats2,corrs,imps1,imps2, m1, m2){
  
  cors = abs(corrs[feats1,feats2])
  
  f.obj = as.vector(t(cors))
  
  f.constraints = matrix(data=0,nrow=length(feats1)+length(feats2),ncol=length(feats2)*length(feats1))
  
  for(i in 1:length(feats1)){
    start = (i-1)*length(feats2)+1
    end = start + length(feats2)-1
    f.constraints[i,start:end] = 1
  }
  
  for(i in 1:length(feats2)){
    indexes = seq(i,length(feats2)*length(feats1),by=length(feats2))
    f.constraints[length(feats1)+i,indexes] = 1
  }
  
  f.rhs = c(imps1[feats1],imps2[feats2])
  
  f.dir = rep("<=", nrow(f.constraints))
  
  r = lp ("max", f.obj, f.constraints, f.dir, f.rhs)
  solution = matrix(r$solution,nrow=length(feats1),ncol=length(feats2),byrow = T)
  
  return(r$objval)
}

normalize = function(coeffs){
  d = ncol(coeffs)
  M = nrow(coeffs)
  kbar = sum(coeffs!=0)/M
  for(m in 1:M){
    sum_coef = sum(abs(coeffs[m,which(coeffs[m,]!=0)]))
    coeffs[m,] = kbar*abs(coeffs[m,])/sum_coef
  }
  coeffs
}

phi_lin = function(coeffs, corrs){
  
  M = nrow(coeffs)
  kbar = sum(coeffs!=0)/M
  coeffs = normalize(coeffs)
  
  sims = matrix(data=0,nrow=M,ncol=M)
  for(m1 in 1:M){
    for(m2 in 1:M){
      if(m1 < m2){
        sims[m1,m2] = sim_lin(which(coeffs[m1,]!=0), which(coeffs[m2,]!=0), corrs, coeffs[m1,], coeffs[m2,], m1, m2)
      }
    }
  }
  return(2*sum(sims)/(kbar*M*(M-1)))
}