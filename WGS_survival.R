################################# WGS ABUNDANCE SURVIVAL #################################
wgs_bacteria_survival<-
  function(abundance_matrix,Source){
    samples_list<-
      wgs_full_annotations%>%
      filter(Sample_Source==Source)%>%
      pull(Barcode)
    
    wgs_clr<-
      abundance_matrix%>%
      select(any_of(c("tax_id",samples_list)))%>%
      # filter(tax_id %in% all_bact_genera)%>%
      left_join(unified_taxonomy[,c("tax_id","name")])%>%
      column_to_rownames("name")%>%
      select(-tax_id)%>%
      select_if(colSums(.)>=100)%>%
      filter(rowMeans(.>=10)>0.1)%>%
      t()%>%clr(.+0.1)%>%t()%>%
      as.data.frame()
    
    # SURVIVAL WITH GENUS ABUNDANCE
    setup_wgs_surv_relabund<-
      wgs_clr%>%
      t()%>%
      as.data.frame()%>%
      rownames_to_column("Barcode")%>%
      left_join(wgs_full_annotations)%>%
      drop_na(VITAL_STATUS,STAGE_simple,SURVIVAL_TIME_WEEKS_DERIVED)%>%
      mutate(stage_binned=if_else(STAGE_simple %in% c("II","III","IV"),"late",STAGE_simple),
             age_over_65=AGE_AT_DIAGNOSIS>65,
             age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
      filter(VITAL_STATUS!="99999",STAGE_simple!="99999")
    
    survival_v_genera_wgs<-
      intersect(rownames(wgs_clr),unified_taxonomy$name)
    
    genera_wgs_survival<-
      map_dfr(survival_v_genera_wgs,
              function(x){
                formula<-as.formula(paste0("survival~strata(site_study)+strata(age_over_65)+strata(stage_binned)+age_by10+HISTOLOGY_simple+`",x,"`"));
                survival<-Surv(setup_wgs_surv_relabund$survival_weeks_10Y, setup_wgs_surv_relabund$vital_status_10Y);
                coxph(formula,data = setup_wgs_surv_relabund)%>%tidy()
              })%>%
      mutate(fdr=p.adjust(p.value,method="BH"))%>%
      filter(!str_detect(term,"_"),!term=="AlterationYes")%>%
      arrange(p.value)%>%
      mutate(log2_estimate=log2(exp(estimate)))
    
    return(genera_wgs_survival)
  }

################################# RUN THE ABOVE FOR WGS and NG232 #################################
genera_wgs_survival_tumor<-wgs_bacteria_survival(wgs_decontamd,"Tumor")
genera_wgs_survival_Blood<-wgs_bacteria_survival(wgs_decontamd,"Blood")
genera_NG232_survival_tumor<-wgs_bacteria_survival(NG232_decontamd,"Tumor")
genera_NG232_survival_Blood<-wgs_bacteria_survival(NG232_decontamd,"Blood")

meta_genera_results_tumor<-
  genera_NG232_survival_tumor%>%
  mutate(dataset="NG232")%>%
  rbind(genera_wgs_survival_tumor%>%mutate(dataset="WGS"))%>%
  summarise(meta_est=sum(estimate/std.error^2)/sum(1/std.error^2),
            meta_err=sqrt(1/sum(1/std.error^2)),
            meta_p= 1-pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)%>%
  mutate(meta_fdr=p.adjust(meta_p))

meta_genera_results_blood<-
  genera_NG232_survival_Blood%>%
  mutate(dataset="NG232")%>%
  rbind(genera_wgs_survival_Blood%>%mutate(dataset="WGS"))%>%
  summarise(meta_est=sum(estimate/std.error^2)/sum(1/std.error^2),
            meta_err=sqrt(1/sum(1/std.error^2)),
            meta_p= 1-pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)%>%
  mutate(meta_fdr=p.adjust(meta_p))

meta_wgs_survival_plot<-
  meta_genera_results_tumor%>%mutate(`Sample Type`="Tumor")%>%
  rbind(meta_genera_results_blood%>%mutate(`Sample Type`="Blood"))%>%
  arrange(-meta_est)%>%
  mutate(experiment="WGS")%>%
  forestplot(name = term,
             estimate = meta_est,
             se=meta_err,
             colour=`Sample Type`)+
  theme_bw()+
  facet_wrap(~experiment)+
  labs(x="ln(Cox hazard ratio)")+
  theme(axis.text.y=element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 10),
        strip.background = element_rect(fill=WGS_color))+
  scale_color_manual(values = c(Tumor="red",Normal="blue",Blood="goldenrod"),drop=F)+
  coord_cartesian(xlim=c(-7,7))
meta_wgs_survival_plot

################################# WGS SURVIVAL AND RICHNESS #################################
wgs_richness_survival<-
  function(richness_dat,Source){
    WGS_richness_survival<-
      richness_dat%>%
      left_join(wgs_full_annotations)%>%
      #mutate(ng232=Barcode %in% natgen_232_barcodes)%>%
      filter(Sample_Source==Source)%>%
      filter(VITAL_STATUS!="99999",
             !is.na(SURVIVAL_TIME_WEEKS_DERIVED),
             STAGE_simple!="99999")%>%
      mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
             age_over_65=AGE_AT_DIAGNOSIS>65,
             age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
      mutate(rarefied_100=ntile(N100,n=2),
             .by = c(site_study,STAGE_DERIVED))
    
    richness_cox<-
      coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(site_study)+strata(stage_binned)+strata(age_over_65)+HISTOLOGY_simple+age_by10+N100,
          data = WGS_richness_survival)
    
    richness_cox_p_wgs<-
      richness_cox%>%
      tidy()%>%
      filter(term=="N100")%>%
      pull(p.value)%>%
      round(digits = 3)
    
    survivalRichness_wgs<-
      survfit(formula=Surv(WGS_richness_survival$survival_weeks_10Y,
                           WGS_richness_survival$vital_status_10Y)~rarefied_100,
              data=WGS_richness_survival)%>%
      ggsurvplot(#conf.int = TRUE,          # Add confidence interval
                 #pval = paste0("Cox p=",richness_cox_p_wgs),              # Add p-value
                 legend.labs =c("Low richness","High richness"),palette = c("gold","blue"),
                 ggtheme = theme_bw(base_family = "Roboto Condensed"),
                 xlab="Time (weeks)")%>%
      .$plot+
      theme(legend.position = "bottom")
    
    return(list(richness_cox,survivalRichness_wgs))
  }
wgs_tumor_richness<-wgs_richness_survival(WGS_richness,"Tumor")
wgs_blood_richness<-wgs_richness_survival(WGS_richness,"Blood")
ng232_tumor_richness<-wgs_richness_survival(NG232_richness,"Tumor")
ng232_blood_richness<-wgs_richness_survival(NG232_richness,"Blood")

tumor_meta_richness_results<-
  rbind(tidy(wgs_tumor_richness[[1]])%>%mutate(data="WGS"),
        tidy(ng232_tumor_richness[[1]])%>%mutate(data="NG232"))%>%
  summarise(meta_est=sum(estimate/std.error^2)/sum(1/std.error^2),
            meta_err=sqrt(1/sum(1/std.error^2)),
            meta_p= 1-pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)
blood_meta_richness_results<-
  rbind(tidy(wgs_blood_richness[[1]])%>%mutate(data="WGS"),
        tidy(ng232_blood_richness[[1]])%>%mutate(data="NG232"))%>%
  summarise(meta_est=sum(estimate/std.error^2)/sum(1/std.error^2),
            meta_err=sqrt(1/sum(1/std.error^2)),
            meta_p= 1-pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)

survivalRichness_wgs<-
  wgs_tumor_richness[[2]]+
  facet_wrap2(~"WGS tumor richness",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(tumor_meta_richness_results$meta_p[nrow(tumor_meta_richness_results)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(tumor_meta_richness_results$meta_est[nrow(tumor_meta_richness_results)]),digits = 2)),
           x = 50,y=0.25)


Blood_survivalRichness_wgs<-
  wgs_blood_richness[[2]]+
  facet_wrap2(~"WGS blood richness",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(blood_meta_richness_results$meta_p[nrow(blood_meta_richness_results)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(blood_meta_richness_results$meta_est[nrow(blood_meta_richness_results)]),digits = 2)),
           x = 50,y=0.25)

################################# REPEAT WITH DIVERSITY #################################
wgs_diversity_survival<-
  function(diversity_data,Source){
    WGS_diversity_survival<-
      diversity_data%>%
      pivot_longer(`100`:`750`,names_to = "depth",values_to = "Shannon diversity")%>%
      left_join(wgs_full_annotations)%>%
      #mutate(ng232=Barcode %in% natgen_232_barcodes)%>%
      filter(Sample_Source==Source,depth==100)%>%
      filter(VITAL_STATUS!="99999",
             !is.na(SURVIVAL_TIME_WEEKS_DERIVED),
             STAGE_simple!="99999")%>%
      mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
             age_over_65=AGE_AT_DIAGNOSIS>65,
             age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
      mutate(rarefied_100=ntile(`Shannon diversity`,n=2),
             .by = c(site_study,STAGE_DERIVED))
    
    diversity_cox<-
      coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(site_study)+strata(stage_binned)+strata(age_over_65)+HISTOLOGY_simple+age_by10+rarefied_100,
          data = WGS_diversity_survival)
    
    diversity_cox_p_wgs<-
      diversity_cox%>%
      tidy()%>%
      filter(term=="rarefied_100")%>%
      pull(p.value)%>%
      round(digits = 3)
    
    survivalDiversity_wgs<-
      survfit(formula=Surv(WGS_diversity_survival$survival_weeks_10Y,
                           WGS_diversity_survival$vital_status_10Y)~rarefied_100,
              data=WGS_diversity_survival)%>%
      ggsurvplot(#conf.int = TRUE,          # Add confidence interval
                 #pval = paste0("Cox p=",diversity_cox_p_wgs),              # Add p-value
                 legend.labs =c("Low diversity","High diversity"),palette = c("gold","blue"),
                 ggtheme = theme_bw(base_family = "Roboto Condensed"),
                 xlab="Time (weeks)")%>%
      .$plot+
      #facet_wrap2(~"WGS tumor diversity",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
      theme(legend.position = "bottom")
    
    return(list(diversity_cox,survivalDiversity_wgs))
  }

wgs_tumor_diversity<-wgs_diversity_survival(WGS_diversity,"Tumor")
wgs_blood_diversity<-wgs_diversity_survival(WGS_diversity,"Blood")
ng232_tumor_diversity<-wgs_diversity_survival(NG232_diversity,"Tumor")
ng232_blood_diversity<-wgs_diversity_survival(NG232_diversity,"Blood")

tumor_meta_diversity_results<-
  rbind(tidy(wgs_tumor_diversity[[1]])%>%mutate(data="WGS"),
        tidy(ng232_tumor_diversity[[1]])%>%mutate(data="NG232"))%>%
  summarise(meta_est=sum(estimate/std.error^2)/sum(1/std.error^2),
            meta_err=sqrt(1/sum(1/std.error^2)),
            meta_p= 1-pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)
blood_meta_diversity_results<-
  rbind(tidy(wgs_blood_diversity[[1]])%>%mutate(data="WGS"),
        tidy(ng232_blood_diversity[[1]])%>%mutate(data="NG232"))%>%
  summarise(meta_est=sum(estimate/std.error^2)/sum(1/std.error^2),
            meta_err=sqrt(1/sum(1/std.error^2)),
            meta_p= 1-pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)

survivalDiversity_wgs<-
  wgs_tumor_diversity[[2]]+
  facet_wrap2(~"WGS tumor diversity",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(tumor_meta_diversity_results$meta_p[nrow(tumor_meta_diversity_results)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(tumor_meta_diversity_results$meta_est[nrow(tumor_meta_diversity_results)]),digits = 2)),
           x = 50,y=0.25)
           
  
Blood_survivalDiversity_wgs<-
  wgs_blood_diversity[[2]]+
  facet_wrap2(~"WGS blood diversity",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(blood_meta_diversity_results$meta_p[nrow(blood_meta_diversity_results)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(blood_meta_diversity_results$meta_est[nrow(blood_meta_diversity_results)]),digits = 2)),
           x = 50,y=0.25)
