rna_combatd_clr<-
  rna_combatd_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  left_join(kraken_taxonomy[,c("tax_id","name")])%>%
  column_to_rownames("name")%>%
  select(-tax_id)%>%
  select_if(colSums(.)>=500)%>%
  filter(rowMeans(.>=50)>0.1)%>%
  t()%>%clr(.+0.1)%>%t()%>%
  as.data.frame()

# SURVIVAL WITH GENUS ABUNDANCE 2-TILES
setup_rna_surv_relabund_normal<-
  rna_combatd_clr%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("RNAseq_SampleID")%>%
  left_join(rna_annots_full)%>%
  filter(Type=="Normal")%>%
  drop_na(VITAL_STATUS,STAGE_simple,SURVIVAL_TIME_WEEKS_DERIVED)%>%
  mutate(VITAL_STATUS=(VITAL_STATUS%%2)+1)%>%
  mutate(vital_status_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,1,VITAL_STATUS),
         survival_weeks_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,520.0,SURVIVAL_TIME_WEEKS_DERIVED),
         stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  filter(VITAL_STATUS!="99999",STAGE_simple!="99999")

setup_rna_surv_relabund_tumor<-
  rna_combatd_clr%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("RNAseq_SampleID")%>%
  left_join(rna_annots_full)%>%
  filter(Type=="Tumor")%>%
  drop_na(VITAL_STATUS,STAGE_simple,SURVIVAL_TIME_WEEKS_DERIVED)%>%
  mutate(VITAL_STATUS=(VITAL_STATUS%%2)+1)%>%
  mutate(vital_status_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,1,VITAL_STATUS),
         survival_weeks_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,520.0,SURVIVAL_TIME_WEEKS_DERIVED),
         stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  filter(VITAL_STATUS!="99999",STAGE_simple!="99999")

survival_v_genera<-
  intersect(colnames(setup_rna_surv_relabund_tumor),kraken_taxonomy$name)

genera_survival_tumor<-
  map_dfr(survival_v_genera,
          function(x){
            formula<-as.formula(paste0("survival~strata(Site_REF)+strata(age_over_65)+strata(stage_binned)+HISTOLOGY_simple+age_by10+`",x,"`"));
            survival<-Surv(setup_rna_surv_relabund_tumor$survival_weeks_10Y, setup_rna_surv_relabund_tumor$vital_status_10Y);
            coxph(formula,data = setup_rna_surv_relabund_tumor)%>%tidy()
          })%>%
  mutate(fdr=p.adjust(p.value,method="BH"))%>%
  filter(!str_detect(term,"_"),!term=="AlterationYes")%>%
  arrange(p.value)

genera_survival_normal<-
  map_dfr(survival_v_genera,
          function(x){
            formula<-as.formula(paste0("survival~strata(Site_REF)+strata(age_over_65)+strata(stage_binned)+HISTOLOGY_simple+age_by10+`",x,"`"));
            survival<-Surv(setup_rna_surv_relabund_normal$survival_weeks_10Y, setup_rna_surv_relabund_normal$vital_status_10Y);
            coxph(formula,data = setup_rna_surv_relabund_normal)%>%tidy()
          })%>%
  mutate(fdr=p.adjust(p.value,method="BH"))%>%
  filter(!str_detect(term,"_"),!term=="AlterationYes")%>%
  arrange(p.value)

genera_survival_plot<-
  genera_survival_tumor%>%mutate(`Sample Type`="Tumor")%>%
  rbind(genera_survival_normal%>%mutate(`Sample Type`="Normal"))%>%
  arrange(-estimate)%>%
  mutate(experiment="RNA-seq")%>%
  forestplot(name = term,
             estimate = estimate,
             se=std.error,
             # pvalue = p.value,
             colour=`Sample Type`)+
  theme_bw()+
  facet_wrap2(~experiment,strip=strip_themed(background_x = element_rect(fill = RNA_color)))+
  labs(x="ln(Cox hazard ratio)")+
  theme(axis.text.y=element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 10))+
  scale_color_manual(values = c(Tumor="red",Normal="blue"))
genera_survival_plot

### Kaplan-Meier curves
setup_genus_KMeier<-
  rna_combatd_clr%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("RNAseq_SampleID")%>%
  left_join(rna_annots_full[,c("RNAseq_SampleID","Site_REF")])%>%
  pivot_longer(where(is.numeric),values_to = "clr",names_to = "name")%>%
  mutate(clr=ntile(clr,n=2),.by = Site_REF)%>%
  left_join(rna_annots_full)%>%
  drop_na(VITAL_STATUS,STAGE_simple)%>%
  mutate(VITAL_STATUS=(VITAL_STATUS%%2)+1)%>%
  mutate(vital_status_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,1,VITAL_STATUS),
         survival_weeks_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,520.0,SURVIVAL_TIME_WEEKS_DERIVED),
         stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65)%>%
  left_join(kraken_taxonomy[,c("tax_id","name")])%>%
  select(-tax_id)%>%
  pivot_wider(names_from = "name",values_from = "clr")%>%
  filter(VITAL_STATUS!="99999",STAGE_simple!="99999")

# KM_1<-
#   survfit(Surv(survival_weeks_10Y,vital_status_10Y)~Microbacterium,
#         data=setup_genus_KMeier,)%>%
#   ggsurvplot(conf.int = TRUE,          # Add confidence interval
#              pval = FALSE,
#              legend.labs =c("Low","High"),palette = c("blue","red"),
#              legend.title="Microbacterium\nrelative abundance",
#              ggtheme = theme_bw(base_family = "Roboto Condensed"),
#              xlab="Time (weeks)")
# KM_1

########################### RELATE SURVIVAL TO DIVERSITY ############################################
ADiversity_survival<-
  RNAseq_diversity%>%
  left_join(rna_annots_full)%>%
  drop_na(STAGE_simple,SURVIVAL_TIME_WEEKS_DERIVED)%>%
  filter(Type=="Tumor",STAGE_simple!="99999",!VITAL_STATUS=="99999")%>%
  mutate(VITAL_STATUS=(VITAL_STATUS%%2)+1)%>%
  mutate(vital_status_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,1,VITAL_STATUS),
         survival_weeks_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,520.0,SURVIVAL_TIME_WEEKS_DERIVED),
         stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  mutate(rarefied_500=ntile(`500`,n=2),rarefied_750=ntile(`750`,n=2),
         rarefied_1000=ntile(`1000`,n=2),rarefied_1500=ntile(`1500`,n=2),.by = c(Site_REF,stage_binned))%>%
  rename(r500="500",r750="750")%>%
  drop_na(rarefied_500)

diversity_cox_p<-
  coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(Site_REF)+strata(stage_binned)+strata(age_over_65)+HISTOLOGY_simple+age_by10+r500,
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
  facet_wrap2(~"RNA-seq tumor diversity",strip = strip_themed(background_x = element_rect(fill=RNA_color)))+
  theme(legend.position = "bottom")+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(diversity_cox_p$p.value[nrow(diversity_cox_p)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(diversity_cox_p$estimate[nrow(diversity_cox_p)]),digits = 2)),
           x = 50,y=0.25)

############################################ RICHNESS #################################

richness_survival<-
  rna_richness%>%
  left_join(rna_annots_full)%>%
  filter(Type=="Tumor")%>%
  filter(VITAL_STATUS!="99999",
         !is.na(SURVIVAL_TIME_WEEKS_DERIVED),
         STAGE_simple!="99999")%>%
  mutate(VITAL_STATUS=(VITAL_STATUS%%2)+1)%>%
  mutate(vital_status_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,1,VITAL_STATUS),
         survival_weeks_10Y=if_else(SURVIVAL_TIME_WEEKS_DERIVED>520,520.0,SURVIVAL_TIME_WEEKS_DERIVED),
         stage_binned=if_else(STAGE_simple %in% c("II,III","IV"),"late",STAGE_simple),
         age_over_65=AGE_AT_DIAGNOSIS>65,
         age_by10=cut_interval(AGE_AT_DIAGNOSIS,length = 10))%>%
  mutate(rarefied_500=ntile(N500,n=2),
         rarefied_750=ntile(N750,n=2),
         rarefied_1000=ntile(N1000,n=2),
         rarefied_1500=ntile(N1500,n=2),
         .by = c(Site_REF,STAGE_DERIVED))

richness_cox_p<-
  coxph(Surv(survival_weeks_10Y,vital_status_10Y)~strata(Site_REF)+strata(stage_binned)+strata(age_over_65)+HISTOLOGY_simple+age_by10+N500,
        data = richness_survival)%>%
  tidy()%>%
  filter(term=="N500")

survivalRichness<-
  survfit(formula=Surv(richness_survival$survival_weeks_10Y,
                       richness_survival$vital_status_10Y)~rarefied_500,
          data=richness_survival)%>%
  ggsurvplot(#pval = paste0("Cox p=",round(richness_cox_p$p.value,digits = 3)),
             legend.labs =c("Low richness","High richness"),palette = c("gold","blue"),
             ggtheme = theme_bw(base_family = "Roboto Condensed"),
             xlab="Time (weeks)",)%>%
  .$plot+
  facet_wrap2(~"RNA-seq tumor richness",strip = strip_themed(background_x = element_rect(fill=RNA_color)))+
  theme(legend.position = "bottom")+
  annotate(geom = "text", 
           label=paste0("p-value=",
                        round(richness_cox_p$p.value[nrow(richness_cox_p)],digits = 2),
                        "\n",
                        "Cox HR=",
                        round(exp(richness_cox_p$estimate[nrow(richness_cox_p)]),digits = 2)),
           x = 50,y=0.25)



