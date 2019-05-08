###### Functions ######
computeBNDistance <- function(ob1,ob2,d){
  
  
  mat1 <- matrix(as.numeric(ob1),nrow = dim(ob1)[1], ncol = dim(ob1)[2])
  mat2 <- matrix(as.numeric(ob2),nrow = dim(ob2)[1], ncol = dim(ob2)[2])
  #print(typeof(mat1))
  #print(typeof(mat2))
  bndist = bottleneck(mat1,mat2,dimension =d)
  return(bndist)
}


calculateLandscape <- function(PD,d,kk,nNodes){
  Diag1 <- matrix(as.numeric(PD),nrow = dim(PD)[1], ncol = dim(PD)[2])
  land <- landscape(Diag1,dimension = d, KK = kk, tseq = 1:nNodes)
  return(land)
}

#makeDiagFromPIs <- function(allPIs,graphNumber,nNodes){
makeDiagFromPIs <- function(pis){  
  # this is a cell now a list need to extract matrices from it
  graphNumber <- 1
  nNodes <- 70
  
 
  pi1 <- pis[graphNumber,1]
  pi2 <- pis[graphNumber,2]
  pi3 <- pis[graphNumber,3]
  
  pi1[pi1>nNodes] <- nNodes+1
  pi2[pi2>nNodes] <- nNodes+1
  pi3[pi3>nNodes] <- nNodes+1
  
  
  # Need to make into PD that TDA functions expect
  if(dim(pi1)[1]>0){
  pi1_1 <- cbind(matrix(1,dim(pi1)[1]),pi1)
  } else {
    pi1_1 <- c()
  }
  
  if(dim(pi2)[1]>0){
  pi2_1 <- cbind(matrix(2,dim(pi2)[1]),pi2)
  } else {
    pi2_1 <- c()
  }
  
  if(dim(pi3)[1]>0){
  pi3_1 <- cbind(matrix(3,dim(pi3)[1]),pi3)
  } else {
    pi3_1 <- c()
  }
  
  Diag1 <- rbind(pi1_1,pi2_1,pi3_1)
  
  return(Diag1)
  
}






makeDiagFromPIs_2 <- function(PIs){
  
  # this is a cell now a list need to extract matrices from it
  

  
  pi1 <- PIs[[1]][[1]]
  pi2 <- PIs[[1]][[2]]
  pi3 <- PIs[[1]][[3]]
  
  # Need to make into PD that TDA functions expect
  if(dim(pi1)[1]>0){
    pi1_1 <- cbind(matrix(1,dim(pi1)[1]),pi1)
  } else {
    pi1_1 <- c()
  }
  
  if(dim(pi2)[1]>0){
    pi2_1 <- cbind(matrix(2,dim(pi2)[1]),pi2)
  } else {
    pi2_1 <- c()
  }
  
  if(dim(pi3)[1]>0){
    pi3_1 <- cbind(matrix(3,dim(pi3)[1]),pi3)
  } else {
    pi3_1 <- c()
  }
  
  Diag1 <- rbind(pi1_1,pi2_1,pi3_1)
  
  
  if(is.null(Diag1)){
    Diag1 <- matrix(0,1,3)
  }
  
  
  return(Diag1)
  
}
