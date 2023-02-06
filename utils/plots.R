library('ggplot2')

createUMAPsPerTopic <- function (cellData, pathToPlotsDir){
  #' Generates UMAP plots for each topic showing the scores of each cell
  #' 

pdf(file.path(pathToPlotsDir, 'Umaps.pdf'))
for (i in 1:9){

  ggplot(data = cellData, 
         aes(x = UMAP1, 
             y = UMAP2, 
             color = paste(c('Topic', i), collapse = "")))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") +
  scale_color_gradient()

}
dev.off()
}