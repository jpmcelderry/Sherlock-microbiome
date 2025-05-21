############################# RELATE SURVIVAL TO TOTAL ABUNDANCE, ALPHA DIVERSITY
ADiversity_survival<-
  RNAseq_diversity%>%
  left_join(rna_annots_full)%>%
  drop_na(STAGE_simple,SURVIVAL_TIME_WEEKS_DERIVED)%>%
  filter(Type=="Tumor",STAGE_simple!="99999",!VITAL_STATUS=="99999",HISTOLOGY_simple=="ADENOCARCINOMA")%>%
  mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  mutate(rarefied_500=ntile(`500`,n=2),rarefied_750=ntile(`750`,n=2),
         rarefied_1000=ntile(`1000`,n=2),rarefied_1500=ntile(`1500`,n=2),.by = c(Site_REF,stage_binned))%>%
  rename(r500="500",r750="750")%>%
  drop_na(rarefied_500)

diversity_cox_p<-
  coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(Site_REF)+strata(stage_binned)+strata(age_over_65)+age_by10+r500,
        data = ADiversity_survival)%>%
  tidy()%>%
  filter(term=="r500")

survivalDiversity<-
  survfit(formula=Surv(survival_weeks_10Y,vital_status_10Y)~rarefied_500,
          data=ADiversity_survival)%>%
  ggsurvplot(#pval = paste0("Cox p=",diversity_cox_p),              # Add p-value
             legend.labs =c("Low diversity","High diversity"),palette = c("gold","blue"),
             ggtheme = theme_bw(base_family = "Roboto Condensed"),
             xlab="Time (weeks)")%>%
  .$plot+
  facet_wrap2(~"RNA-seq LUAD diversity",strip = strip_themed(background_x = element_rect(fill=RNA_color)))+
  theme(legend.position = "bottom")+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(diversity_cox_p$p.value[nrow(diversity_cox_p)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(diversity_cox_p$estimate[nrow(diversity_cox_p)]),digits = 2)),
           x = 50,y=0.25)
  

richness_survival<-
  rna_richness%>%
  left_join(rna_annots_full)%>%
  filter(Type=="Tumor")%>%
  filter(VITAL_STATUS!="99999",
         !is.na(SURVIVAL_TIME_WEEKS_DERIVED),
         STAGE_simple!="99999",
         HISTOLOGY_simple=="ADENOCARCINOMA")%>%
  mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  mutate(rarefied_500=ntile(N500,n=2),
         rarefied_750=ntile(N750,n=2),
         rarefied_1000=ntile(N1000,n=2),
         rarefied_1500=ntile(N1500,n=2),
         .by = c(Site_REF,STAGE_DERIVED))

richness_cox_p<-
  coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(Site_REF)+strata(stage_binned)+strata(age_over_65)+age_by10+N500,
        data = richness_survival)%>%
  tidy()%>%
  filter(term=="N500")

survivalRichness<-
  survfit(formula=Surv(richness_survival$survival_weeks_10Y,
                       richness_survival$vital_status_10Y)~rarefied_500,
          data=richness_survival)%>%
  ggsurvplot(#pval = paste0("Cox p=",round(richness_cox_p$p.value,digits = 3)),##"; HR=",round(1/exp(richness_cox_p$estimate),digits = 2)),
             legend.labs =c("Low richness","High richness"),palette = c("gold","blue"),
             ggtheme = theme_bw(base_family = "Roboto Condensed"),
             xlab="Time (weeks)")%>%
  .$plot+
  facet_wrap2(~"RNA-seq LUAD richness",strip = strip_themed(background_x = element_rect(fill=RNA_color)))+
  theme(legend.position = "bottom")+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(richness_cox_p$p.value[nrow(richness_cox_p)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(richness_cox_p$estimate[nrow(richness_cox_p)]),digits = 2)),
           x = 50,y=0.25)

########################################################################### 16S
s16_richness_survival<-
  s16_richness%>%
  left_join(s16_clinical_annos)%>%
  filter(`Tumor-NormalStatus`=="Tumor")%>%
  filter(VITAL_STATUS!="99999",
         !is.na(SURVIVAL_TIME_WEEKS_DERIVED),
         STAGE_simple!="99999",HISTOLOGY_simple=="ADENOCARCINOMA")%>%
  mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  mutate(rarefied_250=ntile(N250,n=2),
         .by = c(Study,STAGE_DERIVED))

richness_cox_p_16s<-
  coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(Study)+strata(stage_binned)+strata(age_over_65)+age_by10+N250,
        data = s16_richness_survival)%>%
  tidy()%>%
  filter(term=="N250")

survivalRichness_16S<-
  survfit(formula=Surv(s16_richness_survival$survival_weeks_10Y,
                       s16_richness_survival$vital_status_10Y)~rarefied_250,
          data=s16_richness_survival)%>%
  ggsurvplot(#pval = paste0("Cox p=",richness_cox_p_16s),              # Add p-value
    legend.labs =c("Low richness","High richness"),palette = c("gold","blue"),
    ggtheme = theme_bw(base_family = "Roboto Condensed"),
    xlab="Time (weeks)")%>%
  .$plot+
  facet_wrap2(~"16S LUAD richness",strip = strip_themed(background_x = element_rect(fill=s16_color)))+
  theme(legend.position = "bottom")+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(richness_cox_p_16s$p.value[nrow(richness_cox_p_16s)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(richness_cox_p_16s$estimate[nrow(richness_cox_p_16s)]),digits = 2)),
           x = 50,y=0.25)

#### repeat with diversity
s16_diversity_survival<-
  s16_diversity%>%
  left_join(s16_clinical_annos)%>%
  pivot_longer(`100`:`750`,names_to = "depth",values_to = "Shannon diversity")%>%
  filter(`Tumor-NormalStatus`=="Tumor",depth==250)%>%
  filter(VITAL_STATUS!="99999",
         !is.na(SURVIVAL_TIME_WEEKS_DERIVED),
         STAGE_simple!="99999",HISTOLOGY_simple=="ADENOCARCINOMA")%>%
  mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  mutate(rarefied_250=ntile(`Shannon diversity`,n=2),
         .by = c(Study,STAGE_DERIVED))

diversity_cox_p_16s<-
  coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(Study)+strata(stage_binned)+strata(age_over_65)+age_by10+rarefied_250,
        data = s16_diversity_survival)%>%
  tidy()%>%
  filter(term=="rarefied_250")

survivalDiversity_16s<-
  survfit(formula=Surv(s16_richness_survival$survival_weeks_10Y,
                       s16_richness_survival$vital_status_10Y)~rarefied_250,
          data=s16_richness_survival)%>%
  ggsurvplot(#pval = paste0("Cox p=",diversity_cox_p_16s),              # Add p-value
    legend.labs =c("Low diversity","High diversity"),palette = c("gold","blue"),
    ggtheme = theme_bw(base_family = "Roboto Condensed"),
    xlab="Time (weeks)")%>%
  .$plot+
  facet_wrap2(~"16S LUAD diversity",strip = strip_themed(background_x = element_rect(fill=s16_color)))+
  theme(legend.position = "bottom")+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(diversity_cox_p_16s$p.value[nrow(diversity_cox_p_16s)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(diversity_cox_p_16s$estimate[nrow(diversity_cox_p_16s)]),digits = 2)),
           x = 50,y=0.25)

############################################################################ WGS
wgs_richness_survival<-
  function(richness_dat,Source){
    
    WGS_richness_survival<-
      richness_dat%>%
      left_join(wgs_full_annotations)%>%
      #mutate(ng232=Barcode %in% natgen_232_barcodes)%>%
      filter(Sample_Source==Source)%>%
      filter(VITAL_STATUS!="99999",
             !is.na(SURVIVAL_TIME_WEEKS_DERIVED),
             STAGE_simple!="99999",
             HISTOLOGY_simple=="ADENOCARCINOMA")%>%
      mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
             age_over_65=AGE_AT_DIAGNOSIS>65,
             age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
      mutate(rarefied_100=ntile(N100,n=2),
             .by = c(site_study,STAGE_DERIVED))
    
    richness_cox<-
      coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(site_study)+strata(stage_binned)+strata(age_over_65)+age_by10+N100,
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
            .by=term)%>%
  filter(term=="N100")
blood_meta_richness_results<-
  rbind(tidy(wgs_blood_richness[[1]])%>%mutate(data="WGS"),
        tidy(ng232_blood_richness[[1]])%>%mutate(data="NG232"))%>%
  summarise(meta_est=sum(estimate/std.error^2)/sum(1/std.error^2),
            meta_err=sqrt(1/sum(1/std.error^2)),
            meta_p= 1-pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)%>%
  filter(term=="N100")

survivalRichness_wgs<-
  wgs_tumor_richness[[2]]+
  facet_wrap2(~"WGS LUAD richness",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(tumor_meta_richness_results$meta_p[nrow(tumor_meta_richness_results)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(tumor_meta_richness_results$meta_est[nrow(tumor_meta_richness_results)]),digits = 2)),
           x = 50,y=0.25)


Blood_survivalRichness_wgs<-
  wgs_blood_richness[[2]]+
  facet_wrap2(~"WGS LUAD blood richness",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
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
             STAGE_simple!="99999",
             HISTOLOGY_simple=="ADENOCARCINOMA")%>%
      mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
             age_over_65=AGE_AT_DIAGNOSIS>65,
             age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
      mutate(rarefied_100=ntile(`Shannon diversity`,n=2),
             .by = c(site_study,STAGE_DERIVED))
    
    diversity_cox<-
      coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(site_study)+strata(stage_binned)+strata(age_over_65)+age_by10+rarefied_100,
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
            .by=term)%>%
  filter(term=="rarefied_100")
blood_meta_diversity_results<-
  rbind(tidy(wgs_blood_diversity[[1]])%>%mutate(data="WGS"),
        tidy(ng232_blood_diversity[[1]])%>%mutate(data="NG232"))%>%
  summarise(meta_est=sum(estimate/std.error^2)/sum(1/std.error^2),
            meta_err=sqrt(1/sum(1/std.error^2)),
            meta_p= 1-pchisq(-2*sum(log(p.value)),df=2*n()),
            .by=term)%>%
  filter(term=="rarefied_100")

survivalDiversity_wgs<-
  wgs_tumor_diversity[[2]]+
  facet_wrap2(~"WGS LUAD diversity",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(tumor_meta_diversity_results$meta_p[nrow(tumor_meta_diversity_results)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(tumor_meta_diversity_results$meta_est[nrow(tumor_meta_diversity_results)]),digits = 2)),
           x = 50,y=0.25)


Blood_survivalDiversity_wgs<-
  wgs_blood_diversity[[2]]+
  facet_wrap2(~"WGS LUAD blood diversity",strip = strip_themed(background_x = element_rect(fill=WGS_color)))+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(blood_meta_diversity_results$meta_p[nrow(blood_meta_diversity_results)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(blood_meta_diversity_results$meta_est[nrow(blood_meta_diversity_results)]),digits = 2)),
           x = 50,y=0.25)

plot_grid(survivalRichness+theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalDiversity+theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalRichness_16S+theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalDiversity_16s+theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalRichness_wgs+theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalDiversity_wgs+theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          Blood_survivalRichness_wgs+theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          Blood_survivalDiversity_wgs+theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          labels = c("a","b","c","d","e","f","g","h"),ncol=2)
ggsave(filename = "plots/current_plots/survival_adenosOnly.pdf",device = cairo_pdf,height=9,width=8,bg="white")
ggsave(filename = "plots/current_plots/survival_adenosOnly.png",device = png,height=9,width=8,bg="white")
