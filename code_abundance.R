#*******************************************************************************
#********************************************************************************
#
## R scripts "DiP" for computing Chiu et al.'s (2014)'s and Routledge's (1979) alpha, beta, gamma diversities and associated 
# dissimilarity measures (Sorensen-type 1-C and Jaccard-type 1-U). See the following paper for details:
# Chao et al. (2019) Comparing two classes of alpha diversities and their corresponding beta and (dis)similarity measures: 
# with an application to the Formosan sika deer (Cervus nippon taiouanus) reintroduction program 
#
# This code includes three parts:
# (1) Taxonomic version
# (2) Phylogenetic version
# (3) Plot dissimilarity function
# (4) Example
#
# Note: The packages "phytools", "ade4", "dplyr", "reshape2" and "ggplot2",must be installed and loaded before running the scripts
####################################################################################
#
# (1). Taxonomic alpha, beta, gamma diversities and dissimilarity measures
#
####################################################################################
library(phytools)
library(ade4)
library(dplyr)
library(reshape2)
library(ggplot2)

#' Taxo_dissimilarity(data, q) computes taxonomic alpha, beta, gamma diversities and dissimilarity measures.
#' @param data a SxN dataframe of species raw abundance, where S is the number of species in the pooled assemblage and N is the number
#'   of assemblages.
#' @param q a non-negative vector specifying the diversity order
#' @return a list consisting of three objects:
#' (1) $Routledge: a table of alpha, beta, gamma diversities and dissimilarity based on Routledge's framework.
#' (2) $Chiu: a table of alpha, beta, gamma diversities and dissimilarity based on Chiu et al.'s framework.
#' (3) $Size_weight: size weight of each assemblage, which is used in Routledge's decomposition. 
Taxo_dissimilarity <- function(data, q, datatype){
  wi <- colSums(data)/sum(data)
  if(length(wi) != ncol(data)) stop("please check your weight, whose length should be the same as the number of assemblages.")
  if(sum(wi)!=1) wi/sum(wi) else wi
  if(colSums(data)[1] != 1) data <- apply(data,MARGIN = 2,function(x) x/sum(x))
  N <- ncol(data)
  weight_data <- sapply(1:ncol(data), function(i) wi[i]*data[,i] )
  gamma <- sapply(q , function(x) Diversity_profile(rowSums(weight_data), x) )
  alpha_C <- sapply(q , function(x) Diversity_profile(c(weight_data), x)/N )
  beta_C <- gamma/alpha_C
  alpha_R <- sapply(q , function(x) alpha_r(data, x , wi) )
  beta_R <- gamma/alpha_R
  beta_max_R <- sapply(q , function(x) Diversity_profile(weight_data, x)/alpha_r(data, x , wi) )
  beta_max_C <- N
  modified_CqN <- sapply(q, function(x) {
    if(x == 1){
      log(beta_R[which(q==1)])/log(beta_max_R[which(q==1)])
    }else{
      (beta_R[which(q==x)]^(1-x)-1)/(beta_max_R[which(q==x)]^(1-x)-1)
    }
  })
  
  modified_UqN <- sapply(q, function(x) {
    if(x == 1){
      log(beta_R[which(q==1)])/log(beta_max_R[which(q==1)])
    }else{
      (beta_R[which(q==x)]^(x-1)-1)/(beta_max_R[which(q==x)]^(x-1)-1)
    }
  })
  modified_VqN <- sapply(q, function(x)  (beta_R[which(q==x)]-1)/(beta_max_R[which(q==x)]-1) )
  
  modified_SqN <- sapply(q, function(x)  (beta_R[which(q==x)]^(-1)-1)/(beta_max_R[which(q==x)]^(-1)-1) )
  
  
  CqN_C <- sapply(q , function(x) {
    if(x==1) log(beta_C[which(q==1)])/log(ncol(data)) else (beta_C[which(q==x)]^(1-x)-1)/(N^(1-x)-1)
  })
  UqN_C <- sapply(q , function(x) {
    if(x==1) log(beta_C[which(q==1)])/log(ncol(data)) else (beta_C[which(q==x)]^(x-1)-1)/(N^(x-1)-1)
  })
  
  VqN_C <- sapply(q , function(x) (beta_C[which(q==x)]-1)/(N-1) )
  SqN_C <- sapply(q , function(x) (beta_C[which(q==x)]^(-1)-1)/(N^(-1)-1) )
  
  output <- list(Routledge = data.frame(q, gamma,alpha_R,beta_R,beta_max_R,modified_CqN,
                                        modified_UqN,modified_VqN,modified_SqN), 
                 Chiu = data.frame(q,gamma,alpha_C,beta_C,beta_max_C,CqN_C,UqN_C,VqN_C,SqN_C),
                 Size_weight = c(wi))
  colnames(output[[1]]) <- c("Order","gamma","alpha","beta","beta max","1-CqN","1-UqN","1-VqN","1-SqN")
  colnames(output[[2]]) <- c("Order","gamma","alpha","beta","beta max","1-CqN","1-UqN","1-VqN","1-SqN")
  names(wi) = colnames(data)
  output[[1]] <- round(output[[1]],4)
  output[[2]] <- round(output[[2]],4)
  return(output)
}

####################################################################################
#
# (2).  Phylogenetic alpha, beta, gamma diversities and dissimilarity measures
#
####################################################################################
#' Phylo_dissimilarity(data, tree, q) computes the phylogenetic alpha, beta, gamma decomposition and dissimilarity measures.
#' @param data a SxN dataframe of species raw abundance, where S is the number of species in the pooled assembalges and N is
#'  the number of assemblage.
#' @param tree the phylogenetic tree spanned by S species of the pooled assemblage
#' @param q a non-negative vector specifying the diversity order.
#' @return a list consisting of two objects:
#' (1) $Routledge: a table of phylogenetic alpha, beta, gamma diversities and dissimilarity measures based on Routledge's framework.
#' (2) $Chiu: a table of phylogenetic alpha, beta, gamma diversities and dissimilarity measures based on Chiu et al's framework.
#' (3) $Size_weight: size weight of each assemblage, which is used in Routledge's decomposition. 
Phylo_dissimilarity = function(data, tree, q){
  outputR = Beta.Equal(data, "abundance", tree, q)[[1]]
  outputC = Beta.Size(data, "abundance", tree, q)
  weis <- colSums(data)/sum(data)
  output = list( Routledge = outputR, Chiu = outputC, Size_weight = weis)
  return( output )
}
Beta.Equal = function(data, datatype, tree, q.Order){
  tmp = TranMul(data, datatype, tree, method = "relative") 
  G <- function(L,a,q, T_bar){
    if(q==0){
      G = sum(L[a>0])
    }else if(q==1){
      G <- -sum( L[a>0]*a[a>0]/T_bar*log(a[a>0]/T_bar))  
    }else{
      G <- (sum(L*(a/T_bar)^q))
    }
    return(G)
  }
  if(datatype == "abundance"){
    N <- ncol(data)
    pij <- apply(data, 2 , function(x) x/sum(x))
    n = sum(data)
    weights <- colSums(data)/n
  }
  if(datatype == "incidence_raw"){
    N <- length(data)
    size_each <- sapply(data,sum)
    n = sum(size_each)
    weights <- size_each/n
    weights <- rep(1/N,N)
  }
  pool.p <- 0
  aik <- matrix(NA,nrow(tmp$Gamma), N)
  position1 <- which(q.Order==1)
  position_else = c(1:length(q.Order))[-c(position1)]
  q.Orderelse <- q.Order[-c(position1)]
  for(k in 1:N){
    zik <- tmp$Alpha[[k]]$branch_abun
    #zik = zzz[,k]
    # aik[ ,k] <- zik/tail(tmp$Alpha[[k]]$branch_abun,1)
    aik[ ,k] <- zik/max(tmp$Alpha[[k]]$branch_abun)
    etc = weights[k]*aik[ ,k]
    pool.p <- pool.p+etc
  }
  LL <- tmp$Gamma$branch_length
  Tbar=sum(LL%*%pool.p)
  #LL = rep(tmp$Gamma$branch_length,2)
  G1 <- sapply(q.Order,function(q) G(L = LL, a = pool.p, q = q, T_bar = Tbar))
  G1[position1] <- exp(G1[position1])
  G1[position_else] <- G1[position_else]^(1/(1-q.Order[position_else]))
  A1 <- sapply(q.Order,function(q) sum(apply(aik, 2, G, L = LL, q = q, T_bar = Tbar)*weights))
  A1[position1] <- exp(A1[position1])
  A1[position_else] <- A1[position_else]^(1/(1-q.Order[position_else]))
  b <- G1/A1
  Beta_max_Rout <- function(pij,wt,q,N,li,Tbar){
    p <- pij
    sub <- function(q){
      if(q != 1){
        a <- sum(sapply(1:N,function(i) {
          pi <- p[ ,i][p[ ,i]>0]
          lis <- li[p[ ,i]>0]
          wt[i]*sum(lis*pi^q)
        }))^(1/(1-q))
        b <- sum(sapply(1:N,function(i) {
          pi <- p[ ,i][p[ ,i]>0]
          lis <- li[p[ ,i]>0]
          wt[i]^q*sum(lis*pi^q)
        }))^(1/(1-q))
        b/a
      } else {
        exp(-li%*%pij%*%(wt*log(wt))/Tbar)
        #exp(-sum(wt*log(wt)))
      }
    }
    sapply(q, sub)
  }
  
  ans = rep(0, length(q))
  ans2 = rep(0, length(q))
  Beta_max_q <- sapply(q.Order, function(q) Beta_max_Rout(aik,weights,q,N,li = LL,Tbar = Tbar))
  if ( length(position1)!=0 ){
    ans[position1] <- log(b[position1])/log(Beta_max_q[position1])
    ans2[position1] <- log(b[position1])/log(Beta_max_q[position1])
  }
  if (length(q.Orderelse)!=0 ){
    temp <- sapply(position_else, function(i){
      ((b[i])^(1-q.Order[i])-1)/(Beta_max_q[i]^(1-q.Order[i])-1)
    })
    ans[position_else] = temp
    temp <- sapply(position_else, function(i){
      ((b[i])^(q.Order[i]-1)-1)/(Beta_max_q[i]^(q.Order[i]-1)-1)
    })
    ans2[position_else] = temp
  }
  est <- data.frame( b, ans, ans2, G1, A1)
  colnames(est) <- c("Beta", "1-CqN", "1-UqN", "Gamma", "Alpha")
  output = data.frame(q.Order, est)
  colnames(output) = c("Order", "Beta", "1-CqN", "1-UqN", "Gamma", "Alpha")
  output <- list(Mle=output)
  return(output)
}
Beta.Size = function(data, datatype, tree, q){
  tmp = TranMul(data, datatype, tree, "absolute")
  if(datatype == "abundance"){
    N = ncol(data)
    n = sum(data)
    G1 = PhD.q.mle(tmp$Gamma, n, "abundance", q)
    A1 = PhD.q.mle(tmp$Alpha, n, "abundance", q)/N
  }
  if(datatype == "incidence_raw"){
    N = length(data)
    n = sum(unlist(lapply(data, ncol)))
    G1 = PhD.q.mle(tmp$Gamma, n, "incidence_raw", q)
    A1 = PhD.q.mle(tmp$Alpha, n, "incidence_raw", q)/N
  }
  b <- G1/A1
  if( sum(b<1) > 0) b[b<1] = 1 ;   if( sum(b>N) > 0) b[b>N] = N
  trans = TranSim(b, N, q=q)
  est = cbind(b, A1, G1, trans)
  output = data.frame(q, est)
  colnames(output) = c("Order", "Beta","Alpha","Gamma",
                       "1-CqN","1-UqN","1-VqN","1-SqN")
  
  return(output)
}
convToNewick <- function(tree){
  tree<-reorder.phylo(tree,"cladewise")
  n<-length(tree$tip)
  string<-vector(); string[1]<-"("; j<-2
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,2]<=n){
      string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
      if(!is.null(tree$edge.length)){
        string[j]<-paste(c(":",round(tree$edge.length[i],10)), collapse="")
        j<-j+1
      }
      v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
      while(length(v)>0&&k==v[length(v)]){
        string[j]<-")"; j<-j+1
        w<-which(tree$edge[,2]==tree$edge[k,1])
        if(!is.null(tree$edge.length)){
          string[j]<-paste(c(":",round(tree$edge.length[w],10)), collapse="")
          j<-j+1
        }
        v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w
      }
      string[j]<-","; j<-j+1
    } else if(tree$edge[i,2]>=n){
      string[j]<-"("; j<-j+1
    }
  }
  if(is.null(tree$edge.length)) string<-c(string[1:(length(string)-1)], ";")
  else string<-c(string[1:(length(string)-2)],";")
  string<-paste(string,collapse="")
  return(string)
}
PhD.q.mle = function(tmp, n, datatype, q){ # proposed
  PD_obs <- sum(tmp[tmp[,1]>0,2])
  t_bar <- sum(tmp[, 1]*tmp[, 2]/n)
  position0 <- which(q==0)
  position1 <- which(q==1)
  position_else = c(1:length(q))[-c(position0,position1)]
  ans = rep(0, length(q))
  
  if ( length(position0)!=0 ){
    ans[position0] <- PD_obs
  } 
  if ( length(position1)!=0 ){
    tmp <- tmp[tmp[,1]>0,]
    ans[position1] <- exp(-sum(tmp[,2]*(tmp[,1]/t_bar/n)*log((tmp[,1]/t_bar/n))))
  }
  if (length(position_else)!=0 ){
    temp <- sapply(1:length(position_else), function(i){
      qq <- q[position_else[i]]
      ( sum(tmp[,2]*(tmp[,1]/t_bar/n)^(qq)) ) ^ (1/(1-qq))
    })
    ans[position_else] = temp
  }
  return( ans )
}
TranMul = function(data, datatype, tree, method){
  if(datatype == "abundance"){
    rtree = newick2phylog(convToNewick(tree))
    data = data[names(rtree$leaves), ]
    rdata = apply(data, 1, sum)
    GammaTmp = choose_data(rdata, "abundance", rtree)
    GammaTbar = sum(GammaTmp[, 1]*GammaTmp[, 2]) / sum(rdata)
    AlphaTmp = list()
    for(i in 1:ncol(data)){
      adata = aadata = data[, i]
      names(aadata) = rownames(data)
      names(adata) = tree$tip.label
      # tip = tree$tip.label[-match(names(adata[adata>0]), tree$tip.label)]
      # subtree = drop.tip(tree, tip)
      if(length(tree$tip.label)<=2){
        abun <- c(adata[names(adata)%in%tree$tip.label], root = sum(adata))
        tmp = data.frame('branch_abun'=abun, "branch_length" = c(tree$edge.length,0))
      } else{
        treeA = newick2phylog(convToNewick(tree))
        # tmp = choose_data(aadata[aadata>0], "abundance", treeA)
        tmp = choose_data(aadata, "abundance", treeA)#yhc20190219
      }
      # AlphaTbar = sum (tmp[, 1]*tmp[, 2] / sum(adata))
      # if(AlphaTbar < GammaTbar) tmp[nrow(tmp),2] = GammaTbar - AlphaTbar
      AlphaTmp[[i]] = tmp
    }
  }
  if(datatype == "incidence_raw"){
    rtree = newick2phylog(convToNewick(tree))
    data = lapply(data, function(x) x[names(rtree$leaves), ])
    rdata = do.call(cbind, data)
    GammaTmp = choose_data(rdata, "incidence_raw", rtree)
    GammaTbar = sum(GammaTmp[, 1]*GammaTmp[, 2]) / ncol(rdata)
    AlphaTmp = list()
    for(i in 1:length(data)){
      adata = aadata = data[[i]]
      adata = rowSums(adata)
      #tip = tree$tip.label[-match(names(adata[adata>0]), tree$tip.label)]
      #subtree = drop.tip(tree, tip)
      #if(length(subtree$tip.label)<=2){
      if(length(tree$tip.label)<=2){
        #abun <- c(adata[names(adata)%in%subtree$tip.label], root = ncol(aadata))
        abun <- c(adata, root = ncol(aadata))
        #tmp = data.frame('branch_abun'=abun, "branch_length" = c(subtree$edge.length,0))
        tmp = data.frame('branch_abun'=abun, "branch_length" = c(tree$edge.length,0))
      } else{
        #treeA = newick2phylog(convToNewick(subtree))
        #tmp = choose_data(aadata[adata>0, ], "incidence_raw", treeA)
        tmp = choose_data(aadata, "incidence_raw", rtree)
      }
      # AlphaTbar = sum (tmp[, 1]*tmp[, 2] / ncol(aadata))
      # if(AlphaTbar < GammaTbar) tmp[nrow(tmp),2] = GammaTbar - AlphaTbar
      AlphaTmp[[i]] = tmp
    }
  }
  if(method == "absolute"){
    AlphaTmp = do.call(rbind, AlphaTmp)
  }
  
  output = list(Alpha=AlphaTmp, Gamma=GammaTmp)
  return(output)
}#modify abundance case 
TranSim = function(b, N, q=NULL, m=NULL, Q=NULL){
  if(is.null(q) == F){
    position1 <- which(q == 1)
    position_else = c(1:length(q))[-c(position1)]
    ans.C = rep(0, length(q))
    ans.U = rep(0, length(q))
    if ( length(position1)!= 0 ){
      ans.C[position1] = 1 - log(b[position1])/log(N)
      ans.U[position1] = 1 - log(b[position1])/log(N)
    }
    if (length(position_else)!= 0 ){
      temp = sapply(1:length(position_else), function(i){
        ((b[position_else[i]])^(1-q[position_else[i]]) - N^(1-q[position_else[i]]))/(1- N^(1-q[position_else[i]]))
      })
      ans.C[position_else] = temp
      temp = sapply(1:length(position_else), function(i){
        ((b[position_else[i]])^(q[position_else[i]]-1) - N^(q[position_else[i]]-1))/(1- N^(q[position_else[i]]-1))
      })
      ans.U[position_else] = temp
    }
    ans.C = 1 - ans.C
    ans.U = 1 - ans.U
  }
  if(is.null(m) == F){
    if( Q == 1 ) ans.C = ans.U = log(b)/log(N)
    if( Q != 1 ){
      ans.C =  1 - (b^(1-Q) - N^(1-Q))/(1- N^(1-Q))
      ans.U =  1 - (b^(Q-1) - N^(Q-1))/(1- N^(Q-1))
    }
  }
  ans.V = 1 - (N - b) / (N - 1)
  ans.S = 1 - (1/b - 1/N) / (1 - 1/N)
  cbind(ans.C, ans.U, ans.V, ans.S)
}
choose_data = function(data, datatype, tree){
  if(datatype == 'abundance'){
    tmp <- data[names(tree$leaves)]  
    for(i in 1:length(tree$parts)){
      tmp[1+length(tmp)] <- sum(tmp[tree$parts[[i]]])
      names(tmp)[length(tmp)] <- names(tree$parts)[i]
    }
    tmp <- data.frame('branch_abun' = tmp,"branch_length" = c(tree$leaves,tree$nodes))
  }
  if(datatype == 'incidence_raw'){
    data <- data[names(tree$leaves),]
    t <- apply(apply(data, 2, function(x){
      tmp <- x ; names(tmp) <- names(tree$leaves)
      for(i in 1:length(tree$parts)){
        tmp[1+length(tmp)] <- sum(tmp[tree$parts[[i]]])
        names(tmp)[length(tmp)] <- names(tree$parts)[i]
      }
      tmp[tmp>0] <- 1
      return(tmp)
    }), 1, sum )
    tmp <- data.frame('branch_abun'=t, "branch_length" = c(tree$leaves,tree$nodes)) 
  }
  return(tmp)
}
Diversity_profile <- function(x,q){
  pi <- if (sum(x) != 1) x/sum(x) else x
  pi <- pi[pi>0]
  Sub <- function(q){
    if (q == 1) exp(-sum(pi*log(pi))) else exp(1/(1-q)*log(sum(pi^q)))
  }
  sapply(q, Sub)
}
alpha_r <- function(data , q, wi){
  if(sum(wi)!=1) wi/sum(wi) else wi
  if(q==1) exp(sum(wi*log(apply(data, MARGIN = 2, function(x) Diversity_profile(x,1)))))
  else sum(wi*(apply(data, MARGIN = 2, function(x) Diversity_profile(x,q)))^(1-q))^(1/(1-q))
}


####################################################################################
#
# (3). Plot of dissimilarity measures 
#
####################################################################################
#' dis_plot(outcome) plots the output of Taxo_disimilarity or Phylo_disimilarity
#' @param outcome the output of Taxo_disimilarity or Phylo_disimilarity
#' @return a list consisting of two ggplot objects:
#' (1) $Sorsen: plot for Routledge and Chiu et al.'s 1-C q-profile
#' (2) $Jaccard: plot for Routledge and Chiu et al.'s 1-U q-profile
dis_plot <- function(outcome){
  two <- data.frame(q = outcome$Routledge$Order, 
                    Absolute = outcome$Chiu[,"1-CqN"], Relative = outcome$Routledge[,"1-CqN"]) %>% melt(.,id=c("q"))
  dis = "Sorensen-type dissimilarity"
  
  Sorsen <- ggplot(data = two)+geom_line(aes(x = q, y = value, col = variable, lty = variable), size = 2) + theme_bw() + 
    scale_color_manual(values=c(2,4)) + 
    theme(legend.position="bottom", legend.key.width = unit(5,"cm"),
          legend.title = element_blank(), plot.title = element_text(size = 16),
          axis.title=element_text(size=12), axis.text = element_text(size = 12))+
    xlab("Order q") + ylab(dis) 
  
  two <- data.frame(q = outcome$Routledge$Order, 
                    Absolute = outcome$Chiu[,"1-UqN"], Relative = outcome$Routledge[,"1-UqN"]) %>% melt(.,id=c("q"))
  dis = "Jaccard-type dissimilarity"
  
  Jaccard <- ggplot(data = two)+geom_line(aes(x = q, y = value, col = variable, lty = variable), size = 2) + theme_bw() + 
    scale_color_manual(values=c(2,4)) + 
    theme(legend.position="bottom", legend.key.width = unit(5,"cm"),
          legend.title = element_blank(), plot.title = element_text(size = 16),
          axis.title=element_text(size=12), axis.text = element_text(size = 12))+
    xlab("Order q") + ylab(dis) 
  return(list(Sorsen = Sorsen, Jaccard = Jaccard))
}



####################################################################################
#
# (4). Example (Kenting 2008 and 2013 censuses data) 
#
####################################################################################
data <- read.csv("Kenting_census.csv",row.names = 1,header = T)
tree <- read.newick(file = "phylogenetic_tree.txt")
tax <- Taxo_dissimilarity(data = data, q = seq(0,2,0.2), datatype = "abundance")
phyloge <- Phylo_dissimilarity(data = data, tree = tree, q = seq(0,2,0.2))
dis_plot(tax)
dis_plot(phyloge)
