####### 16S survival
s16_clr<-
  s16_scrubbed%>%
  filter(tax_id %in% all_bact_genera)%>%
  left_join(unified_taxonomy[,c("tax_id","name")])%>%
  column_to_rownames("name")%>%
  select(-tax_id)%>%
  select_if(colSums(.)>=250)%>%
  filter(rowMeans(.>=10)>0.1)%>%
  t()%>%clr(.+0.1)%>%t()%>%
  as.data.frame()

# SURVIVAL WITH GENUS ABUNDANCE 2-TILES
setup_16S_surv_relabund_normal<-
  s16_clr%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("SampleID")%>%
  left_join(s16_clinical_annos)%>%
  filter(`Tumor-NormalStatus`=="Normal")%>%
  drop_na(VITAL_STATUS,STAGE_simple,SURVIVAL_TIME_WEEKS_DERIVED)%>%
  mutate(stage_binned=if_else(STAGE_simple %in% c("II","III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  filter(VITAL_STATUS!="99999",STAGE_simple!="99999")

setup_16S_surv_relabund_tumor<-
  s16_clr%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("SampleID")%>%
  left_join(s16_clinical_annos)%>%
  filter(`Tumor-NormalStatus`=="Tumor")%>%
  drop_na(VITAL_STATUS,STAGE_simple,SURVIVAL_TIME_WEEKS_DERIVED)%>%
  mutate(stage_binned=if_else(STAGE_simple %in% c("II","III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  filter(VITAL_STATUS!="99999",STAGE_simple!="99999")

survival_v_genera_16S<-
  intersect(colnames(setup_16S_surv_relabund_tumor),unified_taxonomy$name)

genera_16S_survival_tumor<-
  map_dfr(survival_v_genera_16S,
          function(x){
            formula<-as.formula(paste0("survival~strata(Study)+strata(age_over_65)+strata(stage_binned)+HISTOLOGY_simple+age_by10+`",x,"`"));
            survival<-Surv(setup_16S_surv_relabund_tumor$survival_weeks_10Y, setup_16S_surv_relabund_tumor$vital_status_10Y);
            coxph(formula,data = setup_16S_surv_relabund_tumor)%>%tidy()
          })%>%
  mutate(fdr=p.adjust(p.value,method="BH"))%>%
  filter(!str_detect(term,"_"),!term=="AlterationYes")%>%
  arrange(p.value)%>%
  mutate(estimate2=log2(exp(estimate)))

genera_16S_survival_normal<-
  map_dfr(survival_v_genera_16S,
          function(x){
            formula<-as.formula(paste0("survival~strata(Study)+strata(age_over_65)+strata(stage_binned)+HISTOLOGY_simple+age_by10+`",x,"`"));
            survival<-Surv(setup_16S_surv_relabund_normal$survival_weeks_10Y, setup_16S_surv_relabund_normal$vital_status_10Y);
            coxph(formula,data = setup_16S_surv_relabund_normal)%>%tidy()
          })%>%
  mutate(fdr=p.adjust(p.value,method="BH"))%>%
  filter(!str_detect(term,"_"),!term=="AlterationYes")%>%
  arrange(p.value)%>%
  mutate(estimate2=log2(exp(estimate)))

genera_16S_survival_plot<-
  genera_16S_survival_tumor%>%mutate(`Sample Type`="Tumor")%>%
  rbind(genera_16S_survival_normal%>%mutate(`Sample Type`="Normal"))%>%
  # left_join(survival_prevalences,by=c("term"="name"))%>%
  arrange(-estimate)%>%
  mutate(experiment="16S")%>%
  forestplot(name = term,
             estimate = estimate,
             se=std.error,
             # pvalue = p.value,
             colour=`Sample Type`,)+
  theme_bw()+
  facet_wrap(~experiment)+
  labs(x="ln(Cox hazard ratio)")+
  theme(axis.text.y=element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 10),
        strip.background = element_rect(fill=s16_color))+
  scale_color_manual(values = c(Tumor="red",Normal="blue"))
genera_16S_survival_plot+
  scale_fill_manual(name = "Unadjusted p-value",
                    labels = c("<=0.05", ">0.05"),
                    values = c(1, 22))+guides(fill=guide_legend(order=2))

################################# RICHNESS ############################################
s16_richness_survival<-
  s16_richness%>%
  left_join(s16_clinical_annos)%>%
  filter(`Tumor-NormalStatus`=="Tumor")%>%
  filter(VITAL_STATUS!="99999",
         !is.na(SURVIVAL_TIME_WEEKS_DERIVED),
         STAGE_simple!="99999")%>%
  mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  mutate(rarefied_250=ntile(N250,n=2),
         .by = c(Study,STAGE_DERIVED))

richness_cox_p_16s<-
  coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(Study)+strata(stage_binned)+strata(age_over_65)+HISTOLOGY_simple+age_by10+N250,
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
  facet_wrap2(~"16S tumor richness",strip = strip_themed(background_x = element_rect(fill=s16_color)))+
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
         STAGE_simple!="99999")%>%
  mutate(stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  mutate(rarefied_250=ntile(`Shannon diversity`,n=2),
         .by = c(Study,STAGE_DERIVED))

diversity_cox_p_16s<-
  coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(Study)+strata(stage_binned)+strata(age_over_65)+HISTOLOGY_simple+age_by10+rarefied_250,
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
  facet_wrap2(~"16S tumor diversity",strip = strip_themed(background_x = element_rect(fill=s16_color)))+
  theme(legend.position = "bottom")+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(diversity_cox_p_16s$p.value[nrow(diversity_cox_p_16s)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(diversity_cox_p_16s$estimate[nrow(diversity_cox_p_16s)]),digits = 2)),
           x = 50,y=0.25)
