library('ggplot2')

createUMAPsPerTopic <- function (cellData, pathToPlotsDir){
  #' Generates UMAP plots for each topic showing the scores of each cell

  pdf(file.path(pathToPlotsDir, 'Umaps.pdf'))
  rangeTopics = ncol(cellData[ , grepl( "Topic" , names( cellData ) )]) 
  
  plot <- ggplot(data = cellData, 
                 aes(x = UMAP1, 
                     y = UMAP2, 
                     color = cellType))+
    theme_bw()+
    geom_point()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x = "UMAP1",
         y = "UMAP2",
         subtitle = 'UMAP plot colored by cell type')
  print (plot)
  for (i in 1:rangeTopics ){
  colours = as.numeric(unlist(cellData[paste(c('Topic', i), collapse = "")]))
  plot <- ggplot(data = cellData, 
                 aes(x = UMAP1, 
                     y = UMAP2,
                     color = colours))+
          theme_bw()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
          geom_point()+
          scale_color_gradient(name = paste(c("Topic", i, 'score'),  collapse = " "),
                              low = "grey", high = "red")+
          labs(x = "UMAP1",
               y = "UMAP2",
           subtitle = paste(c("UMAP plot Topic", i),  collapse = " "))
  
  print (plot)
  }
  dev.off()
}

createtSNEsPerTopic <- function (cellData, pathToPlotsDir){
  #' Generates tSNE plots for each topic showing the scores of each cell
  
  pdf(file.path(pathToPlotsDir, 'tSNEs.pdf'))
  rangeTopics = ncol(cellData[ , grepl( "Topic" , names( cellData ) )]) 
  
  plot <- ggplot(data = cellData, 
                 aes(x = tSNE1, 
                     y = tSNE2, 
                     color = cellType))+
    theme_bw()+
    geom_point()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x = "tSNE1",
         y = "tSNE2",
         subtitle = 'tSNE plot colored by cell type')
  print (plot)
  for (i in 1:rangeTopics ){
    colours = as.numeric(unlist(cellData[paste(c('Topic', i), collapse = "")]))
    plot <- ggplot(data = cellData, 
                   aes(x = tSNE1, 
                       y = tSNE2,
                       color = colours))+
      theme_bw()+
      geom_point()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      scale_color_gradient(name = paste(c("Topic", i, 'score'),  collapse = " "),
                           low = "grey", high = "red")+
      labs(x = "tSNE1",
           y = "tSNE2",
           subtitle = paste(c("tSNE plot Topic", i),  collapse = " "))
    
    print (plot)
  }
  dev.off()
}

plotCNNhistory <- function (cnnHistory, pathToPlotsDir){
  
  pdf(file = file.path(pathToPlotsDir, 'lossMetric.pdf'),
      width = 7,
      height =10)
  plot(cnnHistory )
  
  dev.off()
}
