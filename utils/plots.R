library(ggplot2)
library(ggtext)
library(RColorBrewer)

createUMAPsPerTopic <- function (cellData, pathToPlotsDir){
  #' Generates UMAP plots for each topic showing the scores of each cell

  colorVector = c('#D73027', '#F46D43', '#C51B7D', '#FDAE61', '#FEE090','#FFFFBF',
                  '#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#DE77AE','#F1B6DA', 
                  '#D8DAEB', '#B2ABD2', '#8073AC', '#542788', '#E08214','#FDB863',
                  '#8C510A', '#BF812D')
  
  pdf(file.path(pathToPlotsDir, 'Umaps.pdf'))
  rangeTopics = ncol(cellData[ , grepl( "Topic" , names( cellData ) )]) 
  
  plot <- ggplot(data = cellData, 
                 aes(x = UMAP1, 
                     y = UMAP2, 
                     color = cellType))+
    theme_bw()+
    geom_point()+
    scale_colour_manual(values = colorVector)+
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
  
  colorVector = c('#D73027', '#F46D43', '#C51B7D', '#FDAE61', '#FEE090', '#FFFFBF',
                  '#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#DE77AE','#F1B6DA', 
                  '#D8DAEB', '#B2ABD2', '#8073AC', '#542788', '#E08214','#FDB863',
                  '#8C510A', '#BF812D')

  rangeTopics = ncol(cellData[ , grepl( "Topic" , names( cellData ) )]) 
  
  pdf(file.path(pathToPlotsDir, 'tSNEs.pdf'))
  plot <- ggplot(data = cellData, 
                 aes(x = tSNE1, 
                     y = tSNE2, 
                     color = cellType))+
    theme_bw()+
    geom_point()+
    scale_colour_manual(values = colorVector)+
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

plotPredictedvsTrue <- function (yPred, yTrue, pathToPlotsDir){
  #' Generates scatter plots for predicted vs True Values for each variable.
  
  predictedValues = colnames(yPred)

  pdf(file.path(pathToPlotsDir, 'Predictions.pdf'))
  
  for (i in  predictedValues){
    corr <- cor.test(unlist(yPred[i]), unlist(yTest[i]), method=c("pearson"))
    p <- ggplot() +
      geom_point(data = data.frame(x = unlist(yPred[i]), y = unlist(yTest[i])), aes(x = x, y = y))+
      ggtitle(i) +  xlab("Predicted score") + ylab("True score")+
      annotate( geom='richtext', label = paste ("correlation", round(corr$estimate, 2)), x= -Inf, y= Inf,  hjust = -0.5, vjust = 1.25)  
    print (p)
  }
  dev.off()

}