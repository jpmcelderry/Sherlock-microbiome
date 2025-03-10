# alpha diversity with mutations
diversity_vs_genomics<-
  RNAseq_diversity%>%
  left_join(rna_annots_full)%>%
  filter(Type=="Tumor")%>%
  select(-Type)%>%
  left_join(sherlock_data_full,by=c("Subject"))%>%
  filter(Type %in% c("Mutation_Driver","Signature_Cosmic")|
           str_detect(Type,"SCNA")|
           Gene %in% c("Kataegis","HLA_LOH","WGD_Status")
  )%>%
  group_by(Gene,Type)%>%
  filter(nlevels(as.factor(Alteration))>1,sum(Alteration=="Yes")>10)%>%
  mutate(Alteration=if_else(Alteration=="Yes",1,0))%>%
  drop_na(`500`)%>%
  reframe(glm(Alteration~Site_REF+`500`,family = "binomial")%>%tidy())%>%
  filter(str_detect(term,"`500`"))%>%
  ungroup()%>%
  mutate(fdr=p.adjust(p.value,"fdr"))%>%
  arrange(p.value)
diversity_vs_genomics_Plot<-
  diversity_vs_genomics%>%
  mutate(Type=str_replace_all(Type,"_"," "),experiment="RNA-seq Shannon diversity")%>%
  ggplot(aes(x=estimate,y=-log10(p.value),
             label=if_else(p.value<0.05,Gene,""),
             color=Type))+
  geom_point()+
  geom_texthline(yintercept=-log10(0.05),color="red",
                 label="p=0.05",lty=2,family="Roboto Condensed",fontface="italic")+
  geom_texthline(yintercept=-log10(0.000532),color="blue",
                 label="FDR=0.05",lty=2,family="Roboto Condensed",fontface="italic")+
  geom_text_repel(min.segment.length = 0,show.legend = F,nudge_y = 0.1)+
  scale_color_manual(values = SBScolor)+
  labs(x="Binomial regression coefficient")+
  ggh4x::facet_wrap2(~experiment,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=RNA_color)))
ggsave(diversity_vs_genomics_Plot,filename="plots/current_plots/genomics_v_diversity_nonVolcano.png",device = png,height = 4,width=5)

diversity_vs_genomics%>%
  write_tsv("results_tables/Adiversity_v_genomics_RNA.tsv")

richness_v_genomics<-
  rna_richness%>%
  filter(RNAseq_SampleID %in% names(rna.adj.vars))%>%
  left_join(rna_annots)%>%
  filter(Type=="Tumor")%>%
  select(-Type)%>%
  left_join(sherlock_data_full,by=c("Subject"))%>%
  filter(Type %in% c("Mutation_Driver",
                     "Signature_Cosmic")|str_detect(Type,"SCNA")|Gene %in% c("Kataegis","HLA_LOH","WGD_Status"))%>%
  group_by(Gene,Type)%>%
  filter(nlevels(as.factor(Alteration))>1,sum(Alteration=="Yes")>10)%>%
  mutate(Alteration=if_else(Alteration=="Yes",1,0))%>%
  reframe(glm(Alteration~Site_REF+N500,family = "binomial")%>%tidy())%>%
  filter(str_detect(term,"N500"))%>%
  ungroup()%>%
  mutate(fdr=p.adjust(p.value,"fdr"))%>%
  arrange(p.value)
richness_v_genomicsPlot<-
  richness_v_genomics%>%
  mutate(Type=str_replace_all(Type,"_"," "),experiment="RNA-seq genus richness")%>%
  ggplot(aes(x=estimate,y=-log10(p.value),
             label=if_else(p.value<0.05,Gene,""),
             color=Type))+
  geom_point()+
  geom_texthline(yintercept=-log10(0.05),color="red",
                 label="p=0.05",lty=2,family="Roboto Condensed",fontface="italic",hjust=.99)+
  geom_texthline(yintercept=-log10(0.000532),color="blue",
                 label="FDR=0.05",lty=2,family="Roboto Condensed",fontface="italic",hjust=.99)+
  geom_text_repel(min.segment.length = 0,show.legend = F,nudge_y = 0.1)+
  scale_color_manual(values = SBScolor)+
  labs(x="Binomial regression coefficient")+
  ggh4x::facet_wrap2(~experiment,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=RNA_color)))


###### Associate various bacteria with mutational signatures
bacteria_x_mut_signatures<-
  rna_combatd_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>50)>=.1)%>%
  t()%>%
  clr(.+.1)%>%
  as.data.frame()%>%
  rownames_to_column("RNAseq_SampleID")%>%
  pivot_longer(where(is.numeric),names_to = "tax_id",values_to = "clr")%>%
  mutate(zscore=(clr-mean(clr)/sd(clr)),
         .by=tax_id)%>%
  left_join(rna_annots_full[,c("RNAseq_SampleID","Site_REF","Subject","Type","AGE_AT_DIAGNOSIS")])%>%
  filter(Type=="Tumor")%>%
  left_join(sherlock_data_full,by=c("Subject"))%>%
  filter(Type.y=="Signature_Cosmic")%>%
  drop_na(Alteration)%>%
  mutate(Alteration=if_else(Alteration=="Yes",1,0))%>%
  group_by(tax_id,Gene)%>%
  filter(sum(Alteration)>=20)%>%
  summarise(glm(Alteration~Site_REF+AGE_AT_DIAGNOSIS+clr,family = "binomial")%>%tidy(),)%>%
  ungroup()%>%
  filter(term=="clr")%>%
  mutate(fdr=p.adjust(p.value))%>%
  left_join(kraken_taxonomy)
bacteria_x_mut_signatures

# bacteria_x_mut_signatures%>%
#   ggplot(aes(x=name,y=Gene,fill=estimate))+
#   geom_tile(alpha=0.5)+
#   scale_fill_gradient2(low="blue",high="red",mid ="white")+
#   theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45,size=7))

bacteria_x_mut_signatures_plot<-
  bacteria_x_mut_signatures%>%
  mutate(estimate=case_when(estimate< -10 ~ -5,
    estimate>5 ~ 5,
    TRUE ~ estimate
  ))%>%
  mutate(experiment="RNA-seq")%>%
  ggplot(aes(x=Gene,y=name,fill=estimate))+
  geom_tile(color="grey")+
  theme(axis.text.y = element_text(size=8))+
  theme(axis.text.x = element_text(size=8))+
  scale_fill_gradient2(low="blue",high="red",mid="white",
                       breaks=c(-5.0,-2.5,0,2.3),
                       labels=c("<=-5","-2.5","0.0","2.3"))+
  labs(x=NULL,y=NULL,fill="Binomial Regression coefficient")+
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=9),
        legend.position = "top")
bacteria_x_mut_signatures_plot

#omics_correlation_legend<-get_legend(diversity_vs_genomics_Plot)
plot_grid(
  plot_grid(richness_v_genomicsPlot+theme(legend.position = "bottom"),
            diversity_vs_genomics_Plot+guides(color="none"),
            # omics_correlation_legend,
            labels = c("a","b",""),rel_heights = c(1.2,1),ncol=1),
  bacteria_x_mut_signatures_plot,labels=c("","c"),
  ncol=2,rel_widths = c(0.45,0.55))

ggsave(filename="plots/current_plots/omics_v_diversity_nonVolcanos.png",
       device = png,height = 9,width=11,bg="white")
ggsave(filename="plots/current_plots/omics_v_diversity_nonVolcanos.pdf",
       device = cairo_pdf,height = 9,width=11,bg="white")

