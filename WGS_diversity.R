########################################################################
################## BETA DIVERSITY
big_model_wgsAnnos<-
  wgs_full_annotations%>%
  column_to_rownames("Barcode")%>%
  .[labels(wgs_avgdist),]%>%
  filter(Sample_Source=="Tumor")%>%
  mutate(BMI=cut(BMI, breaks=c(-Inf,18.5,24.9,30,Inf),
                 labels=c("under","normal","overweight","obese"))%>%as.character(),
         AGE_AT_DIAGNOSIS=cut(AGE_AT_DIAGNOSIS,breaks=5)%>%as.character())%>%
  replace(.=="99999",NA)%>%
  drop_na(site_study,AGE_AT_DIAGNOSIS,STAGE_simple,HISTOLOGY_COMPOSITE,METASTASIS,VITAL_STATUS,ANCESTRY_DERIVED,SEX_DERIVED,METASTASIS)
big_model_ng232Annos<-
  wgs_full_annotations%>%
  column_to_rownames("Barcode")%>%
  .[labels(NG232_avgdist),]%>%
  filter(Sample_Source=="Tumor")%>%
  mutate(BMI=cut(BMI, breaks=c(-Inf,18.5,24.9,30,Inf),labels=c("under","normal","overweight","obese"))%>%
           as.character(),
         AGE_AT_DIAGNOSIS=cut(AGE_AT_DIAGNOSIS,breaks=5)%>%as.character())%>%
  replace(.=="99999",NA)%>%
  drop_na(site_study,AGE_AT_DIAGNOSIS,STAGE_simple,HISTOLOGY_COMPOSITE,VITAL_STATUS,ANCESTRY_DERIVED,SEX_DERIVED,METASTASIS)

big_model_wgs<-
  adonis2(formula=as.matrix(wgs_avgdist)[rownames(big_model_wgsAnnos),
                                         rownames(big_model_wgsAnnos)]~site_study+AGE_AT_DIAGNOSIS+STAGE_simple+HISTOLOGY_COMPOSITE+METASTASIS+VITAL_STATUS+ANCESTRY_DERIVED+SEX_DERIVED,
                           data = big_model_wgsAnnos,by = "margin",permutations = 999,parallel = 4)
big_model_ng232<-
  adonis2(formula=as.matrix(NG232_avgdist)[rownames(big_model_ng232Annos),
                                         rownames(big_model_ng232Annos)]~site_study+AGE_AT_DIAGNOSIS+STAGE_simple+HISTOLOGY_COMPOSITE+METASTASIS+VITAL_STATUS+ANCESTRY_DERIVED+SEX_DERIVED,
          data = big_model_ng232Annos,by = "margin",permutations = 999,parallel = 4)

wgs_betadiver_meta_results<-
  rbind(tidy(big_model_wgs)%>%mutate(dataset="This study"),
      tidy(big_model_ng232)%>%mutate(dataset="Zhang et al 2021"))%>%
  summarise(meta.p.value= 1 - pchisq(-2*sum(log(p.value)),df=2*n()),
            R2=mean(R2),
            .by=term)%>%
  mutate(meta.fdr=p.adjust(meta.p.value))

################################ A-DIVERSITY ################################
map_wgs_Adiversity<-
  function(adiversity,annotations,testvar){
    testFormula<-
      as.formula(paste0("`100`~site_study+",testvar))
    dataset<-
      adiversity%>%
      left_join(annotations)%>%
      drop_na(`100`,{{testvar}})%>%
      filter(Sample_Source=="Tumor")%>%
      mutate( testvar = fct_lump_min({{testvar}},min = 20,other_level = "Other"))
    glm(data = dataset,formula=testFormula,)%>%
      tidy(control=T)%>%
      mutate(term=str_replace(term,testvar,paste0(testvar," ")))
  }
wgs_aDiversity_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","METASTASIS","ANCESTRY_DERIVED","PASSIVE_SMOKE",
                 "SEX_DERIVED","AGE_AT_DIAGNOSIS","ANY_PREVIOUS_LUNG_DISEASE"),map_wgs_Adiversity,
          adiversity=WGS_diversity,
          annotations=(wgs_full_annotations%>%
                         filter(Barcode %in% colnames(wgs_decontamd))%>%
                         replace(.=="99999",NA)%>%
                         mutate(ANCESTRY_DERIVED=fct_lump_min(ANCESTRY_DERIVED,40,other_level = "Other"))%>%
                         mutate_at("ANCESTRY_DERIVED",relevel,ref="EUR")))%>%
  filter(!str_detect(term,"site_study"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))
NG232_aDiversity_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","METASTASIS","ANCESTRY_DERIVED","PASSIVE_SMOKE",
                 "SEX_DERIVED","AGE_AT_DIAGNOSIS"),map_wgs_Adiversity,
          adiversity=NG232_diversity,
          annotations=(wgs_full_annotations%>%
                         filter(Barcode %in% colnames(NG232_decontamd))%>%
                         replace(.=="99999",NA)%>%
                         mutate(ANCESTRY_DERIVED=fct_lump_min(ANCESTRY_DERIVED,40,other_level = "Other"))%>%
                         mutate_at("ANCESTRY_DERIVED",relevel,ref="EUR")))%>%
  filter(!str_detect(term,"site_study"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))

################################ RICHNESS ################################
map_wgs_richness<-
  function(adiversity,annotations,testvar){
    testFormula<-as.formula(paste0("N100~site_study+",testvar))
    dataset<-adiversity%>%
      left_join(annotations)%>%
      drop_na(N100,{{testvar}})%>%
      filter(Sample_Source=="Tumor")%>%
      mutate( testvar = fct_lump_min({{testvar}},min = 20,other_level = "Other"))
    glm(data = dataset,formula=testFormula,)%>%
      tidy(control=T)%>%
      mutate(term=str_replace(term,testvar,paste0(testvar," ")))
  }
wgs_richness_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","METASTASIS","ANCESTRY_DERIVED","PASSIVE_SMOKE",
                 "SEX_DERIVED","AGE_AT_DIAGNOSIS","ANY_PREVIOUS_LUNG_DISEASE"),map_wgs_richness,
          adiversity=WGS_richness%>%select(Barcode,N100),
          annotations=(wgs_full_annotations%>%
                         filter(Barcode %in% colnames(wgs_decontamd))%>%
                         replace(.=="99999",NA)%>%
                         mutate(ANCESTRY_DERIVED=fct_lump_min(ANCESTRY_DERIVED,40,other_level = "Other"))%>%
                         mutate_at("ANCESTRY_DERIVED",relevel,ref="EUR")))%>%
  filter(!str_detect(term,"site_study"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))
NG232_richness_clinical<-
  map_dfr(.x = c("STAGE_simple","HISTOLOGY_simple","METASTASIS","ANCESTRY_DERIVED","PASSIVE_SMOKE",
                 "SEX_DERIVED","AGE_AT_DIAGNOSIS"),map_wgs_richness,
          adiversity=NG232_richness%>%select(Barcode,N100),
          annotations=(wgs_full_annotations%>%
                         filter(Barcode %in% colnames(NG232_decontamd))%>%
                         replace(.=="99999",NA)%>%
                         mutate(ANCESTRY_DERIVED=fct_lump_min(ANCESTRY_DERIVED,40,other_level = "Other"))%>%
                         mutate_at("ANCESTRY_DERIVED",relevel,ref="EUR")))%>%
  filter(!str_detect(term,"site_study"),!str_detect(term,"Intercept"))%>%
  mutate(fdr=p.adjust(p.value,method="fdr"))

meta_wgs_aDiversity<-
  wgs_aDiversity_clinical%>%
  mutate(dataset="WGS")%>%
  rbind(NG232_aDiversity_clinical%>%mutate(dataset="NG232"))%>%
  summarise(estimate=sum(estimate/std.error^2)/sum(1/std.error^2),
            std.error=sqrt(1/sum(1/std.error^2)),
            p.value= 1 - pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)%>%
  mutate(fdr=p.adjust(p.value))
meta_wgs_richness<-
  wgs_richness_clinical%>%
  mutate(dataset="WGS")%>%
  rbind(NG232_richness_clinical%>%mutate(dataset="NG232"))%>%
  summarise(estimate=sum(estimate/std.error^2)/sum(1/std.error^2),
            std.error=sqrt(1/sum(1/std.error^2)),
            p.value= 1 - pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)%>%
  mutate(fdr=p.adjust(p.value))

wgs_richness_diversity_results<-
  rbind(meta_wgs_richness%>%mutate(var="WGS genus richness"),
        meta_wgs_aDiversity%>%mutate(var="WGS Shannon diversity"))%>%
  rename()%>%
  filter(!str_detect(term,"site_study"))%>%
  mutate_at("term",str_replace,pattern="SMOKE",replacement="SMOKING")%>%
  mutate(term=str_to_title(str_replace(term,"_simple","")%>%str_replace_all("_"," ")))%>%
  mutate(term=str_to_title(str_replace(term,"Derived ","")))%>%
  mutate_at("term",str_replace,"ii$","II")%>%
  mutate_at("term",str_replace,"i$","I")%>%
  mutate_at("term",str_replace,"v$","V")%>%
  mutate_at("term",str_replace,"Diff Final ","")%>%
  mutate(fdr=p.adjust(p.value),term=str_replace(term,"Eas$","EAS"),term=str_replace(term,"Amr","AMR"))

wgs_richness_diversity_forest<-
  wgs_richness_diversity_results%>%
  forestplot(name = term,
             estimate = estimate,
             se=std.error)+
  labs(x="Correlation coefficient")+
  geom_vline(xintercept = 0,lty=3,alpha=0.5)+
  geom_text(mapping = aes(x=estimate,y=term,label=paste0("p=",round(p.value,digits = 3))),
            nudge_y = 0.25,family="Roboto Condensed",size=2.5)+
  theme_bw(base_family = "Roboto Condensed")+
  theme(text = element_text(size = 10))+
  ggh4x::facet_wrap2(~var,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=WGS_color)))
wgs_richness_diversity_forest$layers[[2]]<-NULL
wgs_richness_diversity_forest