
### INDVAL ####


require(labdsv)

indic.sp <- indval(sp.mat[, MySpecies], gr)

# Table of the significant indicator species
pAdj <- p.adjust(indic.sp$pval, "holm")
maxiv <- indic.sp$maxcls #[pAdj <= 0.05]
iv <- indic.sp$indcls #[pAdj <= 0.05]
pv <- indic.sp$pval #[pAdj <= 0.05]
fr <- apply(sp.mat[,MySpecies] > 0, 2, sum)#[pAdj <= 0.05]
fidg <- data.frame(group=maxiv, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(-fidg$group, -fidg$indval),]
fidg$group <- factor(fidg$group, 1:length(gr.name), gr.name)


# Heatmap

indic <- indic.sp$indval %>% 
  tibble::rownames_to_column(var = "Species")%>% 
  arrange(Species) %>% 
  melt(id.vars = "Species") 

indic$Species <- factor(x = indic$Species, levels = rev(row.names(fidg)))
indic$variable <- factor(x = indic$variable, levels = rev(unique(indic$variable)))

ggplot(data = indic, aes(x = variable, y = Species)) +
  geom_tile(aes(fill = value),colour="white") +
  #low= "#CED5F0" , high = "#072175"; low="#E5E5E5", high = "steelblue" #104E8B
  scale_fill_gradient(low="#E5E5E5", high = "red3", limits=c(0.05,0.9), na.value ="white") +
  theme_classic(base_size=9)+
  theme(axis.text.x= element_text(size = 8.5)) +
  labs(fill='IndVal', x= "Community states") 
