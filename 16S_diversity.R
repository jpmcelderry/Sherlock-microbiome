########## ALPHA DIVERSITY
s16_paired_filter<-
  s16_metadata%>%
  filter(`Tumor-NormalStatus` %in% c("Tumor","Normal"))%>%
  count(AdditionalAttributes,`Tumor-NormalStatus`)%>%
  pivot_wider(names_from = "Tumor-NormalStatus",values_from = "n")%>%
  drop_na()%>%
  pull(AdditionalAttributes)
map_16s_Adiversity<-
  function(adiversity,annotations,testvar){
    testFormula<-as.formula(paste0("`Shannon diversity`~Study+",testvar))
    dataset<-adiversity%>%
      left_join(annotations)%>%
      drop_na(`Shannon diversity`,{{testvar}})%>%
      filter(`Tumor-NormalStatus`=="Tumor")%>%
      mutate( testvar = fct_lump_min({{testvar}},min = 20,other_level = "Other"))
    glm(data = dataset,formula=testFormula,)%>%
      tidy(control=T)%>%
      mutate(term=str_replace(term,testvar,paste0(testvar," ")))
  }
s16_aDiversity_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","METASTASIS","ANCESTRY_DERIVED","PASSIVE_SMOKE",
                 "SEX_DERIVED","AGE_AT_DIAGNOSIS","ANY_PREVIOUS_LUNG_DISEASE"),map_16s_Adiversity,
          adiversity=s16_diversity%>%
            pivot_longer(`100`:`750`,names_to = "depth",values_to = "Shannon diversity")%>%
            filter(depth==250)%>%
            select(SampleID,`Shannon diversity`),
          annotations=(s16_clinical_annos%>%
                         filter(SampleID %in% colnames(s16_scrubbed))%>%
                         replace(.=="99999",NA)%>%
                         mutate(ANCESTRY_DERIVED=fct_lump_min(ANCESTRY_DERIVED,40,other_level = "Other"))%>%
                         mutate_at("ANCESTRY_DERIVED",relevel,ref="EUR")))%>%
  filter(!str_detect(term,"Site_REF"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))

# RICHNESS
map_16s_richness<-
  function(adiversity,annotations,testvar){
    testFormula<-as.formula(paste0("N250~Study+",testvar))
    dataset<-adiversity%>%
      left_join(annotations)%>%
      drop_na(N250,{{testvar}})%>%
      filter(`Tumor-NormalStatus`=="Tumor")%>%
      mutate( testvar = fct_lump_min({{testvar}},min = 20,other_level = "Other"))
    glm(data = dataset,formula=testFormula,)%>%
      tidy(control=T)%>%
      mutate(term=str_replace(term,testvar,paste0(testvar," ")))
  }
s16_richness_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","METASTASIS","ANCESTRY_DERIVED","PASSIVE_SMOKE",
                 "SEX_DERIVED","AGE_AT_DIAGNOSIS","ANY_PREVIOUS_LUNG_DISEASE"),map_16s_richness,
          adiversity=s16_richness%>%select(SampleID,N250),
          annotations=(s16_clinical_annos%>%
                         filter(SampleID %in% colnames(s16_scrubbed))%>%
                         replace(.=="99999",NA)%>%
                         mutate(ANCESTRY_DERIVED=fct_lump_min(ANCESTRY_DERIVED,40,other_level = "Other"))%>%
                         mutate_at("ANCESTRY_DERIVED",relevel,ref="EUR")))%>%
  filter(!str_detect(term,"Site_REF"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))

s16_richness_diversity_results<-
  rbind(s16_richness_clinical%>%mutate(var="16S genus richness"),
      s16_aDiversity_clinical%>%mutate(var="16S Shannon diversity"))%>%
  filter(!str_detect(term,"Study"))%>%
  mutate_at("term",str_replace,pattern="SMOKE",replacement="SMOKING")%>%
  mutate(term=str_to_title(str_replace(term,"_simple","")%>%str_replace_all("_"," ")))%>%
  mutate(term=str_to_title(str_replace(term,"Derived ","")))%>%
  mutate_at("term",str_replace,"ii$","II")%>%
  mutate_at("term",str_replace,"i$","I")%>%
  mutate_at("term",str_replace,"v$","V")%>%
  mutate_at("term",str_replace,"Diff Final ","")%>%
  mutate(fdr=p.adjust(p.value),term=str_replace(term,"Eas$","EAS"),term=str_replace(term,"Amr","AMR"))

s16_richness_diversity_forest<-
  s16_richness_diversity_results%>%
  forestplot(name = term,
             estimate = estimate,
             se=std.error)+
  labs(x="Correlation coefficient")+
  geom_vline(xintercept = 0,lty=3,alpha=0.5)+
  geom_text(mapping = aes(x=estimate,y=term,label=paste0("p=",round(p.value,digits = 3))),
            nudge_y = 0.25,family="Roboto Condensed",size=2.5)+
  theme_bw(base_family = "Roboto Condensed")+
  theme(text = element_text(size = 10))+
  ggh4x::facet_wrap2(~var,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=s16_color)))
s16_richness_diversity_forest$layers[[2]]<-NULL
s16_richness_diversity_forest

########################################################################
##### VS IMMUNE CELLS
s16_immune_cells<-
  s16_richness%>%
  filter(N100!=N250)%>%
  left_join(s16_metadata)%>%
  left_join(rna_annots[,c("Sherlock_EAGLE_PID","Type","RNAseq_SampleID","Site_REF")],by=c("AdditionalAttributes"="Sherlock_EAGLE_PID","Tumor-NormalStatus"="Type"))%>%
  drop_na(RNAseq_SampleID)%>%
  filter(nlevels(as.factor(RNAseq_SampleID))==1,.by = SampleID)%>%
  left_join(Danaher_immuneScore,by=c("RNAseq_SampleID"))%>%
  filter(nlevels(as.factor(`Tumor-NormalStatus`))>1,.by = AdditionalAttributes)%>%
  rename(Cell="Immune Cell")%>%
  drop_na(Cell)%>%
  summarise(glm(Score~Study+N250)%>%tidy(),
            .by=c(`Tumor-NormalStatus`,Cell))%>%
  filter(term=="N250")

s16_diversity_immuneCells<-
  s16_diversity%>%
  left_join(s16_metadata)%>%
  left_join(rna_annots[,c("Sherlock_EAGLE_PID","Type","RNAseq_SampleID","Site_REF")],by=c("AdditionalAttributes"="Sherlock_EAGLE_PID","Tumor-NormalStatus"="Type"))%>%
  drop_na(RNAseq_SampleID)%>%
  filter(nlevels(as.factor(RNAseq_SampleID))==1,.by = SampleID)%>%
  left_join(Danaher_immuneScore,by=c("RNAseq_SampleID"))%>%
  filter(nlevels(as.factor(`Tumor-NormalStatus`))>1,.by = AdditionalAttributes)%>%
  rename(Cell="Immune Cell")%>%
  drop_na(Cell)%>%
  summarise(glm(Score~Site_REF+`250`)%>%tidy(),.by=c(`Tumor-NormalStatus`,Cell))%>%
  filter(term=="`250`")

s16_immune_richness_plot<-
  s16_immune_cells%>%
  mutate(fdr=p.adjust(p.value),experiment="16S genus richness")%>%
  left_join(mean_cell_props2,by=c("Tumor-NormalStatus"="Type","Cell"="Immune Cell"))%>%
  ggplot(aes(x=estimate,y=-log10(p.value),
             color=factor(`Tumor-NormalStatus`,levels=c("Tumor","Normal")),
             label=if_else(p.value<0.05,Cell,"")
  ))+
  geom_point(aes(size=mean),alpha=0.3)+
  geom_texthline(yintercept = -log10(0.05),lty=2,label="p=0.05",
                 color="red",fontface="italic",hjust=.99,size=3)+
  # geom_texthline(yintercept = -log10(.00175),lty=2,label="FDR=0.05",
  #                color="blue",fontface="italic",hjust=.99,size=3)+
  geom_text_repel(max.overlaps = 7,show.legend = F,min.segment.length = 0)+
  geom_vline(xintercept = 0,lty=3)+
  labs(y="-log10(p.value)",x="Correlation coefficient",color="Tissue")+
  guides(label="none")+
  ggh4x::facet_wrap2(~experiment,scales = "free",
                     strip = strip_themed(background_x = element_rect(fill=s16_color)))

s16_immune_diversity_plot<-
  s16_diversity_immuneCells%>%
  mutate(fdr=p.adjust(p.value),experiment="16S Shannon diversity")%>%
  left_join(mean_cell_props2,by=c("Tumor-NormalStatus"="Type","Cell"="Immune Cell"))%>%
  ggplot(aes(x=estimate,y=-log10(p.value),
             color=factor(`Tumor-NormalStatus`,levels=c("Tumor","Normal")),
             label=if_else(p.value<0.05,Cell,"")
  ))+
  geom_point(aes(size=mean),alpha=0.3)+
  geom_texthline(yintercept = -log10(0.05),lty=2,label="p=0.05",
                 color="red",fontface="italic",hjust=.99,size=3)+
  # geom_texthline(yintercept = -log10(.0035),lty=2,label="FDR=0.05",
  #                color="blue",fontface="italic",hjust=.99,size=3)+
  geom_text_repel(max.overlaps = 7,show.legend = F,min.segment.length = 0)+
  geom_vline(xintercept = 0,lty=3)+
  labs(y="-log10(p.value)",x="Correlation coefficient",color="Tissue")+
  guides(label="none")+
  ggh4x::facet_wrap2(~experiment,scales = "free",
                     strip = strip_themed(background_x = element_rect(fill=s16_color)))

########################################################################
################## BETA DIVERSITY
paste0(c("s16_avgdist_mat~Study","Sample_Source","STAGE_simple","HISTOLOGY_COMPOSITE",
         "ANY_PREVIOUS_LUNG_DISEASE","VITAL_STATUS",
         "RECURRENCE","PASSIVE_SMOKE","METASTASIS","ANCESTRY_DERIVED","SEX_DERIVED"),collapse = "+")
big_model_16sAnnos<-
  rename(.data = s16_clinical_annos,Sample_Source="Tumor-NormalStatus")%>%
  column_to_rownames("SampleID")%>%
  .[labels(s16_avgdist),]%>%
  mutate(BMI=cut(BMI, breaks=c(-Inf,18.5,24.9,30,Inf),labels=c("under","normal","overweight","obese"))%>%
           as.character(),
         AGE_AT_DIAGNOSIS=cut(AGE_AT_DIAGNOSIS,breaks=5)%>%as.character())%>%
  replace(.=="99999",NA)%>%
  drop_na(Study,AGE_AT_DIAGNOSIS,STAGE_simple,HISTOLOGY_COMPOSITE,VITAL_STATUS,ANCESTRY_DERIVED,SEX_DERIVED,METASTASIS)

big_model_16s<-
  adonis2(formula=as.matrix(s16_avgdist)[rownames(big_model_16sAnnos),
                                         rownames(big_model_16sAnnos)]~Study+Sample_Source+AGE_AT_DIAGNOSIS+STAGE_simple+HISTOLOGY_COMPOSITE+METASTASIS+VITAL_STATUS+ANCESTRY_DERIVED+SEX_DERIVED,
          data = big_model_16sAnnos,by = "margin",permutations = 999,parallel = 4)


