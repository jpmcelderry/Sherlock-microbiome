# library(SKIT)
library(snm)
library(vegan)
library(phyloseq)
library(ggrepel)
library(gridExtra)
library(tidyverse)
library(broom)
library(ggpubr)
library(ggsci)
library(hrbrthemes)
source("ancom_v2.1.R")
source("~/Downloads/CCLasso-master/R/cclasso.R")
source("~/Downloads/CCLasso-master/R/SparCC.R")
library(doParallel)
library(ConQuR)
library(ggh4x)
library(survminer)
library(survival)
library(gtsummary)
library(fast.adonis)
library(shadowtext)
# library(plyr,include.only = c(".","ddply"))
library(rempsyc) #### for pretty publication tables
library(sva)
library(ggbeeswarm)
library(ggforestplot)
library(cowplot)
library(patchwork)
library(geomtextpath)
library(ggupset)
library(ggsankey)
theme_set(theme_bw(base_family = "Roboto Condensed"))

transpose_frame<-
  function(x){
    return(as.data.frame(t(x)))
  }

graph_phyla_composition_barplot<-
  function(pivoted_relabund_dataset,ordering_dataset=pivoted_relabund_dataset,level="phylum",
           ordering_taxa="Actinobacteria",otus=5,formula="~Sample_Source+Smoking",
           abundance_cutoff=0,
           legend_size=8){
    actinobacter_order<-
      ordering_dataset%>%
      filter(type==level,name==ordering_taxa)%>%
      arrange(Relabund)%>%
      pull(Barcode)
    pivoted_relabund_dataset%>%
      filter(type==level)%>%
      filter(Relabund>abundance_cutoff)%>%
      ggplot(aes(y=Relabund,
                 x=factor(Barcode,levels=actinobacter_order),
                 fill=fct_lump(name,n=otus)))+
      geom_bar(position = "fill",stat="identity")+
      ggsci::scale_fill_igv(name=level)+
      facet_grid2(as.formula(formula),scales = "free_x",independent = "x")+
      scale_x_discrete(name=NULL,labels=NULL)+
      theme_bw(base_family="Roboto Condensed")+
      theme(legend.text = element_text(size=legend_size),
            panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
  }

# s16_phylum_bargraph<-
#   function(pivoted_dataset,level="phylum",ordering_taxa="Actinobacteria",otus=5){
#     bacter_order_16s<-
#       pivoted_dataset%>%
#       filter(type==level,name==ordering_taxa)%>%
#       arrange(Relabund)%>%`$`(index)
#     pivoted_dataset%>%
#       replace_na(list(`Tumor-NormalStatus`="negative ctrl"))%>%
#       filter(type==level,`Tumor-NormalStatus`!="Peritumoral",Relabund>0)%>%
#       mutate(name=fct_lump(name,n=otus))%>%
#       ggplot(mapping = aes(y=Relabund,
#                            x=factor(index,levels=bacter_order_16s),
#                            fill=factor(name,sort(name)%>%unique())))+
#       geom_bar(position = "fill",stat="identity")+
#       scale_x_discrete(labels=NULL,name="Samples")+
#       ggsci::scale_fill_igv(name=level)+
#       theme(legend.text = element_text(size=rel(0.5)))+
#       facet_wrap(~`Tumor-NormalStatus`,scales = "free_x")
#   }

# annotate_WGS<-
#   function(data){
#     taxalevel <- c("superkingdom","phylum","class","order","family","genus","species")
#     return(data %>% 
#              left_join(annotations,by=c("file"="Barcode")) %>% 
#              mutate(type=factor(type,levels = taxalevel),
#                     Sample_Source =factor(Sample_Source,levels = c('Blood','Lung','Tumor')))%>%
#              rename(Abundance="value",Barcode="file")
#     )
#   }

# WGS alpha diversity
# plot_WGS_alpha_diversity<-
#   function(dataset,x=Smoking,y=Sample_Source,
#            level="genus",compare_means=NULL,seqs=200,iters=100){
#     print("rarefying")
#     rarefied_diversity<-function(dataset,n){
#       rrarefy(dataset%>%
#                 filter(type==level)%>%
#                 select(Barcode,Reads,taxonomy)%>%
#                 pivot_wider(values_from = "Reads",names_from = "taxonomy")%>%
#                 replace(is.na(.),0)%>%
#                 column_to_rownames("Barcode")%>%
#                 ceiling()%>%
#                 filter(rowSums(.)>n),
#               sample=n)%>%
#         as.data.frame()%>%
#         summarise(Barcode=rownames(.),
#                   shannon=diversity(.,index = "shannon"),
#                   inv_simpson=diversity(.,index="invsimpson"))
#     }
#     all_samplings<-map_dfr(1:iters,~rarefied_diversity(dataset,seqs),.id="iters")
#     print("summarising rarefied results")
#     rare<-all_samplings%>%
#       group_by(Barcode)%>%
#       summarise(shannon=median(shannon),inv_simpson=median(inv_simpson))%>%
#       left_join(dataset%>%select(Barcode,{{x}},{{y}})%>%unique())%>%
#       pivot_longer(c("shannon","inv_simpson"),values_to = "diversity",names_to = "metric")
#     # print("summarising non-rarefied data")
#     # nonrare<-dataset%>%
#     #   filter(type=="genus",Barcode %in% (dataset%>%filter(type=="superkingdom",Reads>seqs))$Barcode)%>%
#     #   group_by(Barcode,{{y}},{{x}})%>%
#     #   summarise(shannon=diversity(Reads,index="shannon"),inv_simpson=diversity(Reads,index="invsimpson"))%>%
#     #   pivot_longer(c("shannon","inv_simpson"),values_to = "diversity",names_to = "metric")
#     print("graphing")
#     # rbind(rare%>%mutate(rarefied=paste0("Rarefied (",seqs," reads)")),
#     #       nonrare%>%mutate(rarefied="Non-rarefied"))%>%
#     rare%>%
#       mutate(rarefied=paste0("Rarefied (",seqs," reads)"))%>%
#       ggplot(mapping=aes(x=paste({{x}},{{y}},sep = "_"),y=diversity,color={{y}}))+
#       stat_summary(fun = mean,geom = "crossbar",width=1,col=1,show.legend = FALSE)+
#       geom_jitter(aes(alpha=0.1),size=0.5)+
#       facet_grid(metric~rarefied,scales = "free_y")+
#       theme_bw(base_family="Roboto Condensed")+
#       theme(legend.position = 'none',axis.text.x = element_text(angle=45, vjust=1,hjust=1),axis.title.x = element_blank())+
#       stat_compare_means(comparisons=compare_means)
#   }

# updated alpha diversity function, only does rarefied but at several rarefaction levels
# plot_alphaDiver_WGS_2<-
#   function(dataset,x=Smoking,y=Sample_Source,level="genus",compare_means=NULL,
#            iters=100,rarelevels=c(100,200,500,1000),violin=F,metric="shannon"){
#     data_as_matrix<-dataset%>%
#       filter(type==level)%>%
#       select(Barcode,Reads,taxonomy)%>%
#       pivot_wider(values_from = "Reads",names_from = "taxonomy")%>%
#       replace(is.na(.),0)%>%
#       column_to_rownames("Barcode")%>%
#       ceiling()
#     print("rarefying")
#     
#     results<-
#       map_dfr(rarelevels,
#               function(x) map_dfr(1:iters,~(data_as_matrix%>%filter(rowSums(.)>x)%>%
#                                               rrarefy(sample = x)%>%
#                                               diversity(index=metric)))%>%
#                 apply(2,median,na.rm=TRUE),.id = "sample_depth")%>%
#       mutate(sample_depth=factor(paste0("Min. ",rarelevels," reads"),levels=paste0("Min. ",rarelevels," reads")))%>%
#       pivot_longer(where(is.numeric),names_to = "Barcode",values_to = "diversity")%>%
#       drop_na()%>%
#       left_join(dataset%>%select(Barcode,Subject,{{x}},{{y}})%>%unique())%>%
#       left_join(annotations_withClinical)%>%
#       group_by(Subject,sample_depth)%>%
#       filter(n()>1)
#     if(violin){
#     plot<-
#       results%>%
#       group_by({{x}},{{y}},sample_depth)%>%
#       mutate(n=paste0(n()," paired samples"))%>%
#       ggplot(mapping=aes(x=paste({{x}},{{y}},sep = " "),y=diversity))+
#       geom_violin(aes(fill=paste({{x}},{{y}},sep = " ")),alpha=0.7,show.legend = F)+
#       geom_boxplot(width=0.1)+
#       facet_grid(~sample_depth,scales = "free")+
#       theme_bw(base_family="Roboto Condensed")+
#       theme(legend.position = 'bottom',legend.title = element_blank(),
#             axis.text.x = element_blank(),
#             axis.title.x = element_blank())+
#       stat_compare_means(comparisons=compare_means)+
#       ylab(paste0(metric," diversity"))
#     }
#     else{
#       plot<-
#         results%>%
#         group_by({{x}},{{y}},sample_depth)%>%
#         mutate(n=paste0(n()," samples"))%>%
#         ggplot(mapping=aes(x=paste({{x}},{{y}},sep = " "),y=diversity))+
#         geom_jitter(aes(fill=paste({{x}},{{y}},sep = " ")),size=2,pch=21,col="gray95",show.legend = F)+
#         geom_boxplot(fill=NA,outlier.shape = NA)+
#         facet_wrap(~sample_depth,scales = "free")+
#         theme_bw(base_family="Roboto Condensed")+
#         theme(legend.position = 'bottom',legend.title = element_blank(),
#               axis.title.x = element_blank())+
#         stat_compare_means(comparisons=compare_means)+
#         ylab(paste0(metric," diversity"))
#     }
#     return(list(results=results,plot=plot))
#   }


# run_adonis2<-
#   function(column,dataset,
#            variables=genetic_alterations_wide,
#            adjust_variables="Study"){
#     formula<-as.formula(paste0("dataset~",adjust_variables,"+`",column,"`"))
#     sherlock_annos<-
#       variables%>%
#       filter(Tumor_Barcode %in% labels(dataset))%>%
#       left_join(annotations,by=c("Tumor_Barcode"="Barcode"))%>%
#       left_join(clinical_data%>%replace_na(list(Histology="NA")))%>%
#       column_to_rownames("Tumor_Barcode")
#     # only run adonis if more than 2 observed values
#     if(sherlock_annos[column]%>%unique()%>%nrow()>1){
#       print(column)
#       adonis2(formula,data = sherlock_annos,by="terms",na.rm=TRUE)%>%
#         tidy()%>%
#         filter(!is.na(p.value))
#     }
#   }

# FUNCTION FOR VISUALIZING BETA DIVERSITY
calc_pcoa<-
  function(distances,annos,var,join_column="Barcode",var1=1,var2=2,dims=5){
    annotation<-
      annos%>%
      rename(Barcode=join_column)%>%
      filter(Barcode %in% labels(distances))
    distances<-
      as.matrix(distances)[annotation$Barcode, annotation$Barcode]%>%
      as.dist()
    return(cmdscale(distances,k=dims,eig = TRUE,list. = TRUE))
  }
graph_betadiver<-
  function(distances,annos,var,join_column="Barcode",var1=1,var2=2,dims=5){
    pcoa_out<-
      calc_pcoa(distances,annos,var,join_column,var1,var2,dims)  
    percent_explained<-
      round((pcoa_out$eig/sum(pcoa_out$eig))*100,digits = 2)
    pcoa_plot<-pcoa_out$points[,c(var1,var2)]%>%  
      as.data.frame()%>%
      rownames_to_column(join_column)%>%
      left_join(annos)%>%
      ggplot(mapping=aes(x=V1,y=V2,color={{var}}))+
      geom_point()+
      stat_ellipse(geom="polygon",aes(fill={{var}}),
                   alpha=0.1,level=0.75,show.legend = FALSE)+
      geom_point(.%>%group_by({{var}})%>%summarise(V1=mean(V1,na.rm=TRUE),V2=mean(V2,na.rm=TRUE)),
                 mapping=aes(fill={{var}}),size=5,color=1,pch=21)+
      xlab(paste0("V",var1," (",percent_explained[var1],"%)"))+
      ylab(paste0("V",var2," (",percent_explained[var2],"%)"))+
      theme_bw(base_family="Roboto Condensed")
    return(list(data=pcoa_out,plot=pcoa_plot))
  }

# PREPARE INPUT FOR RUNNING ANCOM
# prepare_ancom<-
#   function(data, lib_cut, taxonomy_in=taxadata_all, metadata=annotations, level="genus"){
#     ancom_data<-
#       data%>%
#       filter(tax_id %in% (taxonomy_in%>%filter(type==level))$tax_id)%>%
#       select(taxonomy,Barcode,Reads)%>%
#       pivot_wider(names_from = "Barcode",values_from = "Reads")%>%
#       column_to_rownames("taxonomy")
#     meta<-metadata%>%
#       filter(Barcode %in% colnames(ancom_data))
#     return(feature_table_pre_process(ancom_data,
#                                      meta_data = metadata,
#                                      sample_var = "Barcode",
#                                      lib_cut=lib_cut))
#   }

## construct a taxonomy table
build_taxTable<-function(base_taxonomy,levels=c("superkingdom","phylum","class","order","family","genus","species")){
  # working_taxonomy = the current level being added to the taxa table
  # next_table = children of the working level
  
  # initiate the table
  taxa_table<-setNames(data.frame(matrix(nrow=0,ncol=length(levels))),levels)
  
  for(level_indx in 1:(length(levels))){
    #isolate current level and next level
    level<-levels[level_indx]
    next_level<-levels[level_indx+1]
    
    # at the top of the taxonomy, add node to the taxonomy with no parent
    if(level_indx==1){
      working_taxonomy<-
        base_taxonomy%>%
        filter(type==level)%>%
        mutate(parent=NA)
      next_table<-data.frame()
    }
    # copy current level from previous "next level"
    else{
      working_taxonomy<-next_table
      next_table<-data.frame()
    }
    # Loop to add nodes to the taxa table
    for(otu_indx in 1:nrow(working_taxonomy)){
      # extract an otu from the working table
      otu<-working_taxonomy[otu_indx,]
      
      # copy info from otu parent
      if(!is.na(otu$parent)){
        taxa_table[otu$tax_id,]<-taxa_table[otu$parent,]
      }
      # add current otu information at appropriate taxonomic level
      taxa_table[otu$tax_id,level]<-otu$tax_id
      
      # if level is not a leaf, find current otu's children and store in table next_level
      if(!is.na(next_level)){
        next_table<-base_taxonomy%>%
          filter(type==next_level)%>%
          filter(str_detect(pattern = otu$taxonomy,string = taxonomy))%>%
          mutate(parent=otu$tax_id)%>%
          bind_rows(next_table)
      }
    }
  }
  
  return(taxa_table)
}

## recursively propagates decontaminations
decontaminate_up<-function(data,contaminants,decontam_level,taxa_table=tax_table){
  # get levels from taxonomy table, reverse order
  levels<-rev(colnames(taxa_table))
  
  # move taxa ids to rownames
  if("tax_id" %in% colnames(data)){
    # narrow the taxonomy to only the necessary taxa
    taxa_table<-taxa_table[data$tax_id,]
    data<-data%>%
      column_to_rownames("tax_id")
  }
  else{
    taxa_table<-taxa_table[rownames(data),]
  }
  
  # initialize the contaminant frame from the decontam level, all non-contaminants set to 0
  contam_frame<-data%>%
    filter(rownames(.) %in% taxa_table[,decontam_level])
  contam_frame[!(rownames(contam_frame)%in%contaminants),]<-0
  
  # initialize the decontam frame from the decontam level, all contaminants set to 0
  decontam_frame<-data%>%
    filter(rownames(.) %in% taxa_table[,decontam_level])
  decontam_frame[(rownames(decontam_frame)%in%contaminants),]<-0
  
  # iterate over the levels
  for(level in levels[which(levels==decontam_level):length(levels)]){
    
    print(level)
    # decontamination at the decontam level is already done, skip
    if(level==decontam_level){
      next
    }
    
    # get the otus that need decontaming from the overlap of data and taxa at "level"
    otus<-data%>%
      filter(rownames(.) %in% taxa_table[,level])%>%
      rownames()
    # get the sublevel to use for decontamination via the next taxonomy header
    sub_level<-
      colnames(taxa_table)[(which(colnames(taxa_table)==level)+1)]
    #print(paste0(level,",",sub_level))
    
    for(otu in otus){
      # get children (@ sub-level index) from taxa table
      sub_otus<-
        taxa_table[(taxa_table[,level]==otu),sub_level]%>%
        .[!is.na(.)]%>%
        unique()
      # sum the child reads in the decontam frame
      decontam_frac<-decontam_frame[sub_otus,]%>%
        replace(is.na(.),0)%>%
        colSums(na.rm = T)
      # sum the child reads in the contam frame
      contam_frac<-contam_frame[sub_otus,]%>%
        replace(is.na(.),0)%>%
        colSums(na.rm = T)
      # calculate sample-wise percent reads in the contam frame
      contam_read_proportion<-
        contam_frac/(contam_frac+decontam_frac)
      # replace any nans with 0, happens when there's no child reads
      contam_read_proportion[is.nan(contam_read_proportion)]<-0
      # multiply appropriate contam/decontam read percents with the original read count and add to final frame
      decontam_frame<-
        (data[otu,]*(1-contam_read_proportion))%>%
        bind_rows(decontam_frame)
      contam_frame<-
        (data[otu,]*(contam_read_proportion))%>%
        bind_rows(contam_frame)    
    }
  }
  return(list(decontam_frame,contam_frame))
}

calculate_avg_aDiversity<-function(x,depth,index,iters=100){
  map_dfr(1:iters,~(x%>%
                      rrarefy(sample=depth)%>%
                      diversity(index = index)
  ))%>%
    apply(2,median,na.rm=TRUE)%>%
    return()
}