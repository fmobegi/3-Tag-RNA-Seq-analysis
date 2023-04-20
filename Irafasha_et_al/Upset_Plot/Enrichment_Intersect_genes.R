#6hai
Six_hai <- read_csv("Enrichment_Common_6hai_PPP_Rev_1.csv")
Six_hai_plot <- ggplot(data = Six_hai, 
                             aes(x = Six_hai$`Fold Enrichment`, 
                                 y = fct_reorder(Six_hai$Pathway, 
                                                 Six_hai$`Fold Enrichment`,.desc = FALSE))) + 
  geom_point(aes(color=`Enrichment FDR`,size= nGenes)) + 
  geom_segment(aes(x=Six_hai$`Fold Enrichment`, 
                   xend=0, y=Six_hai$Pathway, yend=Six_hai$Pathway, 
                   color= Six_hai$`Enrichment FDR`))  + 
  labs(x= "Fold Enrichment", y="",title = "6_hai")
Six_hai_plot=plot(Six_hai_plot)
Six_hai_plot

#Turquise_module
Turquise_module <- read_csv("Enrichment_Common_12hai_PPP_Rev_1.csv")
Turquise_module_plot <- ggplot(data = Turquise_module, 
                               aes(x = Turquise_module$`Fold Enrichment`, 
                                   y = fct_reorder(Turquise_module$Pathway, 
                                                   Turquise_module$`Fold Enrichment`,.desc = FALSE))) + 
  geom_point(aes(color=`Enrichment FDR`,size= nGenes)) + 
  geom_segment(aes(x=Turquise_module$`Fold Enrichment`, 
                   xend=0, y=Turquise_module$Pathway, yend=Turquise_module$Pathway, 
                   color= Turquise_module$`Enrichment FDR`)) +  
  labs(x= "Fold Enrichment", y="",title = "ME_Turquise")
Plot_Turquise_module=plot(Turquise_module_plot)
Plot_Turquise_module

#plot grid
hubs<- plot_grid(Plot_green_module, Plot_yellow_module, Plot_Turquise_module, nrow = 3, ncol = 1)
ggsave(plot = hubs, filename = "Enrichment_Plots_Grid.pdf", height = 30, width = 25)
dev.off
