###### is blood different from tissue?
blood_v_tissue<-
  wgs_full_annotations%>%
  filter(Barcode %in% colnames(wgs_decontamd))%>%
  mutate(Tissue_Type=if_else(Sample_Source=="Blood","Blood","Lung/Tumor"))%>%
  filter(nlevels(as.factor(Tissue_Type))>1,.by=SUBJECT_ID)%>%
  column_to_rownames("Barcode")
blood_v_tissue_avgdist<-
  wgs_decontamd%>%
  column_to_rownames("tax_id")%>%
  select(any_of(rownames(blood_v_tissue)))%>%
  filter(rowMeans(.>0)>0.01)%>%
  t()%>%
  avgdist(100,iterations = 50)
blood_v_tissue<-
  blood_v_tissue[labels(blood_v_tissue_avgdist),]
blood_v_tissue_results<-
  adonis2(formula=as.matrix(blood_v_tissue_avgdist)[rownames(blood_v_tissue),
                                                    rownames(blood_v_tissue)]~site_study+Tissue_Type,
        data = blood_v_tissue,
        by = "margin",
        permutations = 999,
        parallel = 4)

wgs_cmdscale<-
  cmdscale(blood_v_tissue_avgdist,
           list. = T,eig=T,k=10)
wgs_cmdscale$eig<-
  (wgs_cmdscale$eig/sum(wgs_cmdscale$eig))*100

# Plot
TissueVBlood_betadiversity_plot<-
  wgs_cmdscale$points%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("Barcode")%>%
  left_join(wgs_full_annotations)%>%
  mutate(Tissue_Type=if_else(Sample_Source=="Blood","Blood","Lung/Tumor"))%>%
  filter(nlevels(as.factor(Tissue_Type))>1,.by=SUBJECT_ID)%>%
  # plotting
  ggplot(aes(x=PC1,y=PC2,group=Sample_Source,label=NULL))+
  geom_point(aes(fill=Sample_Source),color="gray80",shape=21)+
  stat_ellipse(geom="polygon",level=0.80,fill=NA,color="black",size=1.02)+
  stat_ellipse(aes(color=Sample_Source),geom="polygon",level=0.80,fill=NA)+
  labs(x=paste0("PC1 (",round(wgs_cmdscale$eig[1],1),"%)"),
       y=paste0("PC2 (",round(wgs_cmdscale$eig[2],1),"%)"))+
  labs(fill="Type",color="Type")+
  facet_wrap2(~"WGS",strip = strip_themed(background_x = element_rect(fill=WGS_color)))
TissueVBlood_betadiversity_plot

ggsave("plots/current_plots/blood_betadiversity.pdf",device=cairo_pdf,width=7,height=6)
ggsave("plots/current_plots/blood_betadiversity.png",device=png,width=7,height=6)

########## ALPHA DIVERISITY
# Isolate blood samples from diversity/richness tables
blood_diversity<-
  WGS_diversity%>%
  filter(Barcode %in% (annotations_withClinical%>%
                         filter(Sample_Source=="Blood")%>%
                         pull(Barcode)))
blood_richness<-
  WGS_richness%>%
  filter(Barcode %in% (annotations_withClinical%>%
                         filter(Sample_Source=="Blood")%>%
                         pull(Barcode)))
blood_ng232_diversity<-
  NG232_diversity%>%
  filter(Barcode %in% (annotations_withClinical%>%
                         filter(Sample_Source=="Blood")%>%
                         pull(Barcode)))
blood_ng232_richness<-
  NG232_richness%>%
  filter(Barcode %in% (annotations_withClinical%>%
                         filter(Sample_Source=="Blood")%>%
                         pull(Barcode)))

# diversity-clinical associations
map_Adiversity_blood<-
  function(adiversity,annos,testvar){
    testFormula<-as.formula(paste0("`100`~site_study+",testvar))
    dataset<-adiversity%>%
      left_join(annos)%>%
      drop_na(`100`,{{testvar}})%>%
      mutate( testvar = fct_lump_min({{testvar}},min = 20,other_level = "Other"))
    glm(data = dataset,formula=testFormula) %>% 
      tidy() %>% 
      mutate(term=str_replace(term,testvar,paste0(testvar," ")))
  }
blood_aDiversity<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple",
                 "RECURRENCE","VITAL_STATUS"),
          map_Adiversity_blood,
          adiversity=blood_diversity,
          annos=(wgs_full_annotations%>%
                         replace(.=="99999",NA)))%>%
  filter(!str_detect(term,"site_study"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))%>%
  mutate(term=toupper(str_replace(term,"_simple","")%>%str_replace("_"," ")))
blood_ng232_aDiversity<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","VITAL_STATUS"),
          map_Adiversity_blood,
          adiversity=blood_ng232_diversity,
          annos=(wgs_full_annotations%>%
                   replace(.=="99999",NA)))%>%
  filter(!str_detect(term,"site_study"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))%>%
  mutate(term=toupper(str_replace(term,"_simple","")%>%str_replace("_"," ")))
blood_meta_aDiversity<-
  rbind(blood_aDiversity%>%mutate(data="NG232"),
      blood_ng232_aDiversity%>%mutate(data="WGS"))%>%
  summarise(estimate=sum(estimate/std.error^2)/sum(1/std.error^2),
            std.error=sqrt(1/sum(1/std.error^2)),
            p.value= 1 - pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)%>%
  mutate(meta_fdr=p.adjust(p.value))

# richness-clinical association
map_richness<-
  function(richness,annos,testvar){
    testFormula<-as.formula(paste0("N100~site_study+",testvar))
    dataset<-richness%>%
      left_join(annos)%>%
      drop_na(N100,{{testvar}})%>%
      mutate( testvar = fct_lump_min({{testvar}},min = 20,other_level = "Other"))
    glm(data = dataset,formula=testFormula)%>%
      tidy()%>%
      mutate(term=str_replace(term,testvar,paste0(testvar," ")))
}
blood_richness_clinical<-
  map_dfr(.x = c("STAGE_simple",
                 "HISTOLOGY_simple",
                 "RECURRENCE",
                 "VITAL_STATUS"),
          .f = map_richness,
          richness=blood_richness,
          annos=(wgs_full_annotations%>%replace(.=="99999",NA)))%>%
  filter(!str_detect(term,"site_study"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))%>%
  mutate(term=toupper(str_replace(term,"_simple","")%>%str_replace("_"," ")))

blood_ng232_richness_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","VITAL_STATUS"),
          .f = map_richness,
          richness=blood_ng232_richness,
          annos=(wgs_full_annotations%>%replace(.=="99999",NA)))%>%
  filter(!str_detect(term,"site_study"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))%>%
  mutate(term=toupper(str_replace(term,"_simple","")%>%str_replace("_"," ")))

blood_meta_richness<-
  rbind(blood_ng232_richness_clinical%>%mutate(data="NG232"),
      blood_richness_clinical%>%mutate(data="WGS"))%>%
  summarise(estimate=sum(estimate/std.error^2)/sum(1/std.error^2),
            std.error=sqrt(1/sum(1/std.error^2)),
            p.value= 1 - pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)%>%
  mutate(meta_fdr=p.adjust(p.value))

# create diversity/richness forest plots
bloodAdiver_forest<-
  # Join tables and clean up labels
  rbind(blood_meta_aDiversity%>%mutate(metric="Blood Shannon diversity"),
          blood_meta_richness%>%mutate(metric="Blood genus richness"))%>%
  mutate_at("term",str_replace,pattern="SMOKE",replacement="SMOKING")%>%
  mutate_at("term",str_replace,pattern="GRADE DIFF_FINAL",replacement="Grade")%>%
  mutate(term=str_to_title(str_replace(term,"_simple","")%>%str_replace_all("_"," ")))%>%
  mutate(term=str_to_title(str_replace(term,"Derived ","")))%>%
  mutate_at("term",str_replace,"ii$","II")%>%
  mutate_at("term",str_replace,"i$","I")%>%
  mutate_at("term",str_replace,"v$","V")%>%
  filter(!str_detect(term,"99999"))%>%
  # plot
  forestplot(name = term,
             estimate = estimate,
             se=std.error)+
  labs(x="Correlation coefficient")+
  geom_vline(xintercept = 0,lty=3,alpha=0.5)+
  geom_text(mapping = aes(x=estimate,y=term,label=paste0("p=",round(p.value,digits = 3))),
            nudge_y = 0.25,family="Roboto Condensed",size=2.5)+
  theme_bw(base_family = "Roboto Condensed")+
  theme(text = element_text(size = 10))+
  ggh4x::facet_wrap2(~metric,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=WGS_color)))
bloodAdiver_forest$layers[[2]]<-NULL
bloodAdiver_forest

########################################################################
### BETA DIVERSITY
blood_avgdist<-
  wgs_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  select(-tax_id)%>%
  select_if(colnames(.) %in% (annotations_withClinical%>%
                                pull(Barcode)))%>%
  filter(rowMeans(.>0)>0.01)%>%
  t()%>%
  avgdist(100,iterations = 50)
ng232_blood_avgdist<-
  NG232_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  select(-tax_id)%>%
  select_if(colnames(.) %in% (annotations_withClinical%>%
                                pull(Barcode)))%>%
  filter(rowMeans(.>0)>0.01)%>%
  t()%>%
  avgdist(100,iterations = 50)
########## Wrangle annotations
big_model_bloodAnnos<-
  wgs_full_annotations%>%
  column_to_rownames("Barcode")%>%
  .[labels(blood_avgdist),]%>%
  filter(Sample_Source=="Blood")%>%
  mutate(BMI=cut(BMI, breaks=c(-Inf,18.5,24.9,30,Inf),
                 labels=c("under","normal","overweight","obese"))%>%
           as.character(),
         AGE_AT_DIAGNOSIS=cut(AGE_AT_DIAGNOSIS,breaks=5)%>%as.character())%>%
  drop_na(VITAL_STATUS,STAGE_simple,HISTOLOGY_simple)%>%
  replace(is.na(.),"99999")
big_model_bloodAnnos_ng232<-
  wgs_full_annotations%>%
  column_to_rownames("Barcode")%>%
  .[labels(ng232_blood_avgdist),]%>%
  filter(Sample_Source=="Blood")%>%
  mutate(BMI=cut(BMI, breaks=c(-Inf,18.5,24.9,30,Inf),
                 labels=c("under","normal","overweight","obese"))%>%
           as.character(),
         AGE_AT_DIAGNOSIS=cut(AGE_AT_DIAGNOSIS,breaks=5)%>%as.character())%>%
  drop_na(VITAL_STATUS,STAGE_simple,HISTOLOGY_simple)%>%
  replace(is.na(.),"99999")

########## Run adonis model
big_model_blood<-
  adonis2(formula=as.matrix(blood_avgdist)[rownames(big_model_bloodAnnos),
                                           rownames(big_model_bloodAnnos)]~site_study+STAGE_simple+HISTOLOGY_simple+RECURRENCE+VITAL_STATUS,
          data = big_model_bloodAnnos,by = "margin",permutations = 999,parallel = 4)

big_model_ng232_blood<-
  adonis2(formula=as.matrix(ng232_blood_avgdist)[rownames(big_model_bloodAnnos_ng232),
                                           rownames(big_model_bloodAnnos_ng232)]~site_study+STAGE_simple+HISTOLOGY_simple+RECURRENCE+VITAL_STATUS,
          data = big_model_bloodAnnos_ng232,by = "margin",permutations = 999,parallel = 4)

blood_beta_metaP<-
  rbind(tidy(big_model_blood)%>%mutate(data="WGS"),
                        tidy(big_model_ng232_blood)%>%mutate(data="NG232"))%>%
  summarise(meta.p.value= 1 - pchisq(-2*sum(log(p.value)),df=2*n()),
            R2=mean(R2),
            .by=term)%>%
  mutate(meta.fdr=p.adjust(meta.p.value))

# Plot adonis results
betadiversity_blood<-
  blood_beta_metaP%>%
  filter(!term %in% c("Total","Residual"))%>%
  mutate(term=str_to_title(str_replace(term,"_simple","")%>%str_replace("_"," ")))%>%
  ggplot(aes(x=term,y=R2))+
  geom_bar(stat="identity",width=0.6,position = position_dodge())+
  geom_text(aes(label = case_when(meta.p.value<0.0001 ~ "****",
                                  meta.p.value<0.001 ~ "***",
                                  meta.p.value<0.01 ~ "**",
                                  meta.p.value<0.05 ~ "*",
                                  .default = ""), y=(R2+.01*R2)),
            position = position_dodge(width = .6), size = 12/.pt)+
  theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1),text = element_text(size = 9))+
  scale_fill_brewer(type="qual")+
  labs(y="Proportion of variance explained",x=NULL)+
  facet_wrap2(~"Blood beta diversity",strip = strip_themed(background_x = element_rect(fill=WGS_color)))
betadiversity_blood

plot_grid(bloodAdiver_forest,
          betadiversity_blood,
          labels=c("a","b"),nrow=1,rel_widths = c(1.5,1))

ggsave("plots/current_plots/blood_supplement.pdf",device=cairo_pdf,height=5.5,width=8)
ggsave("plots/current_plots/blood_supplement.png",height=5.5,width=8)