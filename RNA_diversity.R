paired_subjects_rna<-
  rna_annots%>%
  count(Subject,Type)%>%
  pivot_wider(names_from = "Type",values_from = "n")%>%
  drop_na()%>%
  pull(Subject)

map_Adiversity<-
  function(adiversity,annotations,testvar){
  testFormula<-as.formula(paste0("`500`~Site_REF+",testvar))
  dataset<-adiversity%>%
    left_join(annotations)%>%
    drop_na(`500`,{{testvar}})%>%
    filter(Type=="Tumor")%>%
    mutate( testvar = fct_lump_min({{testvar}},min = 20,other_level = "Other"))
  glm(data = dataset,formula=testFormula,)%>%
    tidy(control=T)%>%
    mutate(term=str_replace(term,testvar,paste0(testvar," ")))
}
RNA_aDiversity_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","ANCESTRY_DERIVED","PASSIVE_SMOKE",
                 "SEX_DERIVED","AGE_AT_DIAGNOSIS","ASTHMA","ANY_PREVIOUS_LUNG_DISEASE"),map_Adiversity,
          adiversity=RNAseq_diversity,
          annotations=(rna_annots_full%>%
                         filter(RNAseq_SampleID %in% colnames(rna_combatd_decontamd))%>%
                         replace(.=="99999",NA)%>%
                         mutate(ANCESTRY_DERIVED=fct_lump_min(ANCESTRY_DERIVED,40,other_level = "Other"))%>%
                         mutate_at("ANCESTRY_DERIVED",relevel,ref="EUR")))%>%
  filter(!str_detect(term,"Site_REF"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))

map_richness<-
  function(richness,annotations,testvar){
    testFormula<-as.formula(paste0("N500~Site_REF+",testvar))
    dataset<-richness%>%
      left_join(annotations)%>%
      drop_na(N500,{{testvar}})%>%
      filter(Type=="Tumor")%>%
      mutate( testvar = fct_lump_min({{testvar}},min = 20,other_level = "Other"))
    glm(data = dataset,formula=testFormula,)%>%
      tidy()%>%
      mutate(term=str_replace(term,testvar,paste0(testvar," ")))
  }
RNA_richness_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","ANCESTRY_DERIVED","PASSIVE_SMOKE",
                 "SEX_DERIVED","AGE_AT_DIAGNOSIS","ASTHMA","ANY_PREVIOUS_LUNG_DISEASE"),map_richness,
          richness=rna_richness,
          annotations=(rna_annots_full%>%
                         filter(RNAseq_SampleID %in% colnames(rna_combatd_decontamd))%>%
                         replace(.=="99999",NA)%>%
                         mutate(ANCESTRY_DERIVED=fct_lump_min(ANCESTRY_DERIVED,40,other_level = "Other"))%>%
                         mutate_at("ANCESTRY_DERIVED",relevel,ref="EUR")))%>%
  filter(!str_detect(term,"Site_REF"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))

RNA_aDiversity_clinical%>%
  write_tsv("results_tables/Adiversity_v_clinicalVars_RNA.txt")

RNA_richness_diversity_results<-
  rbind(RNA_richness_clinical%>%mutate(var="RNA-seq genus richness"),
        RNA_aDiversity_clinical%>%mutate(var="RNA-seq Shannon diversity"))%>%
  mutate_at("term",str_replace,pattern="SMOKE",replacement="SMOKING")%>%
  mutate(term=str_to_title(str_replace(term,"_simple","")%>%str_replace_all("_"," ")))%>%
  mutate(term=str_to_title(str_replace(term,"Derived ","")))%>%
  mutate_at("term",str_replace,"ii$","II")%>%
  mutate_at("term",str_replace,"i$","I")%>%
  mutate_at("term",str_replace,"v$","V")%>%
  mutate(fdr=p.adjust(p.value),term=str_replace(term,"Eas$","EAS"),term=str_replace(term,"Amr","AMR"))
richness_diversity_Forest<-
  RNA_richness_diversity_results%>%
  forestplot(name = term,
             estimate = estimate,
             se=std.error)+
  labs(x="Correlation coefficient")+
  geom_vline(xintercept = 0,lty=3,alpha=0.5)+
  geom_text(mapping = aes(x=estimate,y=term,label=paste0("p=",round(p.value,digits = 3))),
            nudge_y = 0.25,family="Roboto Condensed",size=2.5)+
  theme_bw(base_family = "Roboto Condensed")+
  theme(text = element_text(size = 10))+
  ggh4x::facet_wrap2(~var,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=RNA_color)))
richness_diversity_Forest$layers[[2]]<-NULL
richness_diversity_Forest

rna_richness_country<-
  rna_richness%>%
  left_join(rna_annots_full)%>%
  filter(N500!=N200)%>%
  mutate(Country=case_when(Site_REF=="Mayo"~"Minnesota",
                           Site_REF=="Roswell Park"~"New York",
                           Site_REF=="INCAN"~"Mexico City",
                           Site_REF=="Yale"~"Connecticut",
                           .default=Site_REF
  ))%>%
  ggplot(aes(x=reorder(Country,FUN = median,X = N200),y=N200))+
  geom_violin()+
  geom_quasirandom(alpha=0.075,fill=NA)+
  geom_boxplot(width=0.2,outlier.shape = NA)+
  labs(x=NULL,y="Genus richness")+
  stat_compare_means(label.x.npc = 0.4)+
  theme(axis.text.x = element_text(vjust = 1,hjust=1,angle = 45))

hospital_v_aDiversity<-
  RNAseq_diversity%>%
  left_join(rna_annots_full)%>%
  drop_na(`500`)%>%
  mutate(Country=case_when(Site_REF=="Mayo"~"Minnesota",
                           Site_REF=="Roswell Park"~"New York",
                           Site_REF=="INCAN"~"Mexico City",
                           Site_REF=="Yale"~"Connecticut",
                           .default=Site_REF
  ))%>%
  ggplot(aes(x=reorder(Country,FUN = median,X = `500`),y=`500`))+
  geom_violin()+
  geom_quasirandom(alpha=0.075,fill=NA)+
  geom_boxplot(width=0.2,outlier.shape = NA)+
  labs(x=NULL,y="Shannon diversity")+
  stat_compare_means(label.x.npc = 0.4)+
  theme(axis.text.x = element_text(vjust = 1,hjust=1,angle = 45))
  
# plot_grid(rna_richness_country,hospital_v_aDiversity,
#           rna_richness_histology,histology_v_aDiversity,
#           labels = c("a","b","c","d"))
# ggsave(filename = "plots/current_plots/suppFigure_adiverAbundance.pdf",device=cairo_pdf,width=8,height=8)
# ggsave(filename = "plots/current_plots/suppFigure_adiverAbundance.png",device=png,width=8,height=8)

# # SCATTER PLOT CODE
Danaher_immuneScore<-
  read_delim("extra_data/deconvolution/Sherlock1+2_edge_normalized_logCPM_Danaher_median.txt")%>%
  pivot_longer(-1,names_to = "RNAseq_SampleID",values_to = "Score")%>%
  rename(`Immune Cell`="Genelist")
mean_cell_props2<-
  Danaher_immuneScore%>%
  left_join(rna_annots[,c("RNAseq_SampleID","Type")])%>%
  filter(RNAseq_SampleID %in% colnames(rna_combatd_decontamd))%>%
  summarise(mean=mean(Score),.by=c(`Immune Cell`,Type))

###### Richness v immune component
richness_v_immuneCells<-
  rna_richness%>%
  left_join(rna_annots_full)%>%
  left_join(Danaher_immuneScore,by=c("RNAseq_SampleID"))%>%
  filter(nlevels(as.factor(Type))>1,
         .by=`Subject ID`)%>%
  rename(Cell="Immune Cell")%>%
  drop_na(Cell)%>%
  summarise(glm(Score~Site_REF+N500)%>%tidy(),
            .by=c(Cell,Type))%>%
  filter(term=="N500")

richness_v_immuneCells_plot<-
  richness_v_immuneCells%>%
  mutate(fdr=p.adjust(p.value),experiment="RNA-seq genus richness")%>%
  left_join(mean_cell_props2,by=c("Type","Cell"="Immune Cell"))%>%
  ggplot(aes(x=estimate,y=-log10(p.value),
             color=factor(Type,levels=c("Tumor","Normal")),
             label=if_else(fdr<0.05,Cell,"")
  ))+
  geom_point(aes(size=mean),alpha=0.3)+
  geom_texthline(yintercept = -log10(0.05),lty=2,label="p=0.05",
                 color="red",fontface="italic",hjust=.99,size=3)+
  geom_texthline(yintercept = -log10(0.0023),lty=2,label="FDR=0.05",
                 color="blue",fontface="italic",hjust=.99,size=3)+
  geom_text_repel(max.overlaps = 7,show.legend = F,min.segment.length = 0)+
  geom_vline(xintercept = 0,lty=3)+
  labs(y="-log10(p.value)",x="Correlation coefficient",color="Tissue")+
  guides(label="none")+
  ggh4x::facet_wrap2(~experiment,scales = "free_x",
                     strip = strip_themed(background_x = element_rect(fill=RNA_color)))
ggsave(richness_v_immuneCells_plot,filename = "plots/current_plots/richness_v_cells_Danaher.pdf",device=cairo_pdf,width=6,height=5)

diversity_v_immuneCells<-
  RNAseq_diversity%>%
  left_join(rna_annots_full)%>%
  left_join(Danaher_immuneScore,by=c("RNAseq_SampleID"))%>%
  filter(nlevels(as.factor(Type))>1,.by=`Subject ID`)%>%
  rename(Cell="Immune Cell")%>%
  drop_na(Cell)%>%
  summarise(glm(Score~Site_REF+`500`)%>%tidy(),
            .by=c(Cell,Type))%>%
  filter(term=="`500`")

diversity_v_immuneCells_plot<-
  diversity_v_immuneCells%>%
  mutate(fdr=p.adjust(p.value),experiment="RNA-seq Shannon diversity")%>%
  left_join(mean_cell_props2,by=c("Type","Cell"="Immune Cell"))%>%
  ggplot(aes(x=estimate,y=-log10(p.value),
             color=factor(Type,levels=c("Tumor","Normal")),
             label=if_else(p.value<0.05,Cell,"")
  ))+
  geom_point(aes(size=mean),alpha=0.3)+
  geom_texthline(yintercept = -log10(0.00175),lty=2,label="FDR=0.05",
                 color="blue",fontface="italic",hjust=.99,size=3)+
  geom_texthline(yintercept = -log10(0.05),lty=2,label="p=0.05",
                 color="red",fontface="italic",hjust=.99,size=3)+
  geom_text_repel(max.overlaps = 7,show.legend = F,min.segment.length = 0)+
  geom_vline(xintercept = 0,lty=3)+
  labs(y="-log10(p.value)",x="Correlation coefficient",color="Tissue")+
  guides(label="none")+
  ggh4x::facet_wrap2(~experiment,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=RNA_color)))
ggsave(diversity_v_immuneCells_plot,filename = "plots/current_plots/diversity_v_cells_Danaher.pdf",device=cairo_pdf,width=6,height=5)

###############################################################
########################### BETA-DIVERSITY

big_model_rnaAnnos<-
  rename(.data = rna_annots_full,Sample_Source="Type")%>%
  column_to_rownames("RNAseq_SampleID")%>%
  .[labels(rna_combatd_avgdist),]%>%
  mutate(BMI=cut(BMI, breaks=c(-Inf,18.5,24.9,30,Inf),labels=c("under","normal","overweight","obese"))%>%
           as.character(),
         AGE_AT_DIAGNOSIS=cut(AGE_AT_DIAGNOSIS,breaks=5)%>%as.character())%>%
  replace(.=="99999",NA)%>%
  drop_na(Site_REF,AGE_AT_DIAGNOSIS,STAGE_simple,HISTOLOGY_COMPOSITE,VITAL_STATUS,ANCESTRY_DERIVED,SEX_DERIVED)

big_model_rna2<-
  adonis2(formula=as.matrix(rna_combatd_avgdist)[rownames(big_model_rnaAnnos),
                                                 rownames(big_model_rnaAnnos)]~Site_REF+Sample_Source+AGE_AT_DIAGNOSIS+STAGE_simple+HISTOLOGY_COMPOSITE+VITAL_STATUS+ANCESTRY_DERIVED+SEX_DERIVED,
          data = big_model_rnaAnnos,by = "margin",permutations = 999,parallel = 4)

# betadisper_rna_annos<-
#   big_model_rnaAnnos%>%
#   rownames_to_column("RNAseq_SampleID")%>%
#   filter(RNAseq_SampleID %in% labels(rna_combatd_avgdist))%>%
#   pull(var = VITAL_STATUS,name = RNAseq_SampleID)
# betadisper(as.dist(as.matrix(rna_combatd_avgdist)[names(betadisper_rna_annos),
#                                             names(betadisper_rna_annos)]),
#            group = betadisper_rna_annos)%>%anova()
# 
# rna_combatd_cmdscale<-cmdscale(rna_combatd_avgdist,list. = T,eig=T,k=10)
# rna_combatd_cmdscale$eig<-(rna_combatd_cmdscale$eig/sum(rna_combatd_cmdscale$eig))*100
# loadings<-
#   rna_combatd_decontamd%>%
#   # select(tax_id,any_of(names(rna.adj.vars)))%>%
#   filter(tax_id %in% all_bact_genera)%>%
#   column_to_rownames("tax_id")%>%
#   filter(rowMeans(.>10)>0.05)%>%
#   .[,labels(rna_combatd_avgdist)]%>%
#   apply(MARGIN = 2,FUN = function(x){x/sum(x)})%>%
#   apply(MARGIN = 1,FUN = function(x){x/sum(x)})%>%
#   t()%>%
#   `%*%`(rna_combatd_cmdscale$points)%>%
#   as.data.frame()%>%
#   `colnames<-`(paste0("PC",seq(1:10)))%>%
#   rownames_to_column("tax_id")%>%
#   left_join(kraken_taxonomy)%>%
#   mutate(magnitude=sqrt(PC1^2+PC2^2))%>%
#   arrange(-magnitude)
# rna_combatd_cmdscale$points%>%
#   as.data.frame()%>%
#   `colnames<-`(paste0("PC",seq(1:10)))%>%
#   rownames_to_column("RNAseq_SampleID")%>%
#   left_join(rna_annots_full%>%
#               left_join(read_delim("extra_data/RNA_extraction_batch.txt")%>%
#                           select(`Sherlock PID`,`Tissue Attribute`,`Created By`,`Source Material`)%>%
#                           unique(),
#                         by=c("Sherlock_PID"="Sherlock PID","Type"="Tissue Attribute")))%>%
#   ggplot(aes(x=PC1,y=PC2,label=NULL))+
#   geom_point(aes(fill=`Created By`),color="gray80",shape=21)+ 
#   stat_ellipse(aes(group=`Created By`),geom="polygon", # code for ellipse outlines
#                level=0.80,fill=NA,color="black",size=1.02)+
#   stat_ellipse(aes(color=`Created By`),geom="polygon", # code for colored ellipses
#                level=0.80,fill=NA)+
#   labs(x=paste0("PC1 (",round(rna_combatd_cmdscale$eig[1],1),"%)"),
#        y=paste0("PC2 (",round(rna_combatd_cmdscale$eig[2],1),"%)"))+
#   geom_point(data=loadings[1:15,],aes(x=PC1,y=PC2),inherit.aes = F,shape=17)+ # code for genera dots
#   geom_text_repel(data=loadings[1:15,],aes(x=PC1,y=PC2,label=name), # code for genera labels
#                   inherit.aes = F,bg.colour="white",
#                   max.overlaps = 18,force=4,min.segment.length = 0)+
#   scale_color_manual(values = SBScolor)+
#   scale_fill_manual(values = SBScolor)
# 
# 

  
