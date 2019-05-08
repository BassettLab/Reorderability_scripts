#### Running stats on local reordering graphs



calculate_transitivity <- function(badj){

G <- graph_from_adjacency_matrix(badj, mode = "undirected")  
tr <- transitivity(G, type = "undirected")

return(tr) 
  
}





calculate_communities_bu <- function(badj){
  
  G <- graph_from_adjacency_matrix(badj, mode = "undirected")
  comms <- cluster_louvain(G)
  comm_membership <- membership(comms)
  return(comm_membership)
  
  
}


calculate_communities_wu <- function(adj){
  
  G <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
  comms <- cluster_louvain(G)
  comm_membership <- membership(comms)
  Q <- modularity(comms)
  print(paste("Modularity = ",Q,"for the weighted graph", sep = ""))
  
  return(list(comm_membership,Q))
  
}



calculate_connected_components <- function(badj){
  
  G <- graph_from_adjacency_matrix(badj, mode = "undirected")
  comps <- components(G)

  g_comps <- comps$membership

  
  return(g_comps)
  
}



calculate_top_overlap <- function(badj){
  
  nNodes <- ncol(badj)
  tolap_mat <- matrix(1,nNodes,nNodes)


  for(i in 1:nNodes){
    ki <- sum(badj[i, ])
   


    for(j in i:nNodes){
      if(j>i){

        kj <- sum(badj[j, ])
        
        ci <- badj[i,-c(i,j)]
        cj <- badj[j,-c(i,j)]
        
        l_ij <- sum(ci*cj)  # R only does element-wise multiplication
  
        num <- l_ij + badj[i,j]
        den <- max(ki,kj)
   
        tolap_mat[i,j] <- num/den
        tolap_mat[j,i] <- num/den
     
      }
      
      
    }
    
  }
  
  
  return(tolap_mat)
  
}


lm_eqn <- function(df){
  
  # found at https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
  
  m <- lm(yvar ~ xvar, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


plotGGplot <- function(Xdata,Ydata,filename){
  
  df <- data.frame(xvar = Xdata,yvar = Ydata)
  
  plt1 <- ggplot(df, aes(x = xvar, y = yvar)) + 
    geom_smooth(method = lm) +
    theme_bw()
  #pdf(paste(filename,"1.pdf",sep = ""), width = 5, height = 5)
  #print(plt1)
  #dev.off()
  ggplot2::ggsave(paste(filename,"1.pdf",sep = ""), device = "pdf", width = 4, height = 4)
  print("saved 1")
  
  plt2 <- ggplot(df, aes(x = xvar, y = yvar)) + 
    stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) + 
    scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) + coord_fixed() +
    geom_smooth(method = lm) +
    theme_bw()

  
  pdf(paste(filename,"2.pdf",sep = ""))
  print(plt2)
  dev.off()
  print("saved 2")
  
  
  
}