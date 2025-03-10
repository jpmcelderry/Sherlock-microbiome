tax_table_16S<-
  taxonomy_16S%>%
  filter(tax_id %in% s16_kraken$tax_id)%>%
  mutate(tax_id=as.character(tax_id))%>%
  build_taxTable(levels=c("superkingdom","phylum","class","order","family","genus"))

### this method is the same as the other decontamination, 
### but handles fractional decontamination rather than binary "yes/no" contaminant removal
### Uses the raw and decontaminated dataframes to calculate contaminant percent at the decontam level
decontaminate_up_scrub<-function(data,orig_data,decontam_level,taxa_table=tax_table){
  # get levels from taxonomy table, reverse order
  levels<-rev(colnames(taxa_table))
  # narrow the taxa table to only the necessary taxa
  taxa_table<-taxa_table[data$tax_id,]
  
  # set up the data to numbers only
  if("tax_id" %in% colnames(data)){
    # narrow the taxa table to only the necessary taxa
    taxa_table<-taxa_table[data$tax_id,]
    data<-data%>%
      column_to_rownames("tax_id")
  }
  else{
    taxa_table<-taxa_table[rownames(data),]
  }
  if("tax_id" %in% colnames(orig_data)){
    orig_data<-orig_data%>%
      column_to_rownames("tax_id")
  }
  
  
  # initialize the decontam frame
  decontam_frame<-data%>%
    filter(rownames(.) %in% taxa_table[,decontam_level])
  
  # subtract to find 
  contam_frame<-orig_data%>%
    .[rownames(decontam_frame),colnames(decontam_frame)]
  contam_frame<-(contam_frame-decontam_frame)
  
  # iterate over the levels
  for(level in levels[which(levels==decontam_level):length(levels)]){
    
    print(level)
    # decontamination at the lowest level is already done, skip
    if(level==decontam_level){
      next
    }
    
    # get the otus that need de-contaming from the overlap of data taxa and taxa at "level"
    otus<-data%>%
      filter(rownames(.) %in% taxa_table[,level])%>%
      rownames()
    # get the sublevel to use for decontamination by the next header in the table
    sub_level<-
      colnames(taxa_table)[(which(colnames(taxa_table)==level)+1)]
    #print(paste0(level,",",sub_level))
    
    for(otu in otus){
      # get children (@ sub-level index) from taxa table where level=otu
      sub_otus<-
        taxa_table[taxa_table[,level]==otu,sub_level]%>%
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
      # calculate sample-wise percent total reads in the contam frame
      contam_read_proportion<-
        contam_frac/(contam_frac+decontam_frac)
      # replace any nans with 0, happens when no species-level reads for the genus
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

# filter to only genus-level taxa
s16_preScrub<-
  s16_kraken%>%
  filter(tax_id %in% (taxonomy_16S%>%filter(type=="genus")%>%pull(tax_id)))%>%
  column_to_rownames("tax_id")%>%
  select_if(colSums(.)>0)%>%
  t()

# reformat 16s metadata for use with scrub
scrub16s_metadata<-
  s16_metadata%>%
  select(SampleID,`Sample Type`,PCR_Plate,Well)%>%
  mutate(contam_source=(!`Sample Type`=="Study"),
         Well=str_replace(Well,"_",""))%>%
  column_to_rownames("SampleID")%>%
  select(contam_source,`Sample Type`,Well,PCR_Plate)%>%
  .[rownames(s16_preScrub),]

# pull list of plates to iterate over
scrub_plates<-
  scrub16s_metadata%>%
  filter(!contam_source)%>%
  pull(PCR_Plate)%>%
  unique()

scrubbed_data<-list()
scrub_results_full<-list()
for (plate in scrub_plates) {
  working_meta<-
    scrub16s_metadata%>%
    filter(PCR_Plate==plate)%>%
    select(-PCR_Plate)
  
  scrub_results_full[[plate]]<-SCRuB::SCRuB(s16_preScrub[rownames(working_meta),],
                                       working_meta)
  scrubbed_data[[plate]]<-scrub_results_full[[plate]]$decontaminated_samples
}
# merge scrubbed data by plate into one data frame
scrubbed_data<-do.call("rbind", scrubbed_data)
scrubbed_data<-
  scrubbed_data%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("tax_id")

# combine the scrubbed genus data back with raw family+ data
s16_postScrub<-
  s16_kraken%>%
  select(colnames(scrubbed_data))%>%
  filter(!tax_id %in% (taxonomy_16S%>%filter(type %in% c("genus","species"))%>%pull(tax_id)))%>%
  rbind(.,scrubbed_data)

# propagate the decontamination upward
s16_scrubbed_part1<-decontaminate_up_scrub(data = s16_postScrub,
                       orig_data = s16_kraken,
                       decontam_level = "genus",taxa_table = tax_table_16S)%>%.[[1]]
# remove known contaminators
s16_scrubbed<-
  decontaminate_up(data = s16_scrubbed_part1,
                               contaminants = c("1912216","561","1386",salter_list_nonHuman$tax_id),
                               decontam_level = "genus",
                   taxa_table = tax_table_16S)%>%.[[1]]%>%
  rownames_to_column("tax_id")
# remove excluded sherlock samples (determined here based on missing study site info)
bad_16s_samples<-
  s16_clinical_annos%>%
  filter(is.na(STUDY_SITE)|str_detect(HISTOLOGY_COMPOSITE,"EXCLUDE")|AdditionalAttributes %in% sherlock_excluded$`Sherlock PID`)%>%
  pull(SampleID)
s16_scrubbed<-
  s16_scrubbed%>%
  select(-any_of(bad_16s_samples))%>%
  filter(tax_id %in% kraken_microbeTaxonomy$tax_id)

### 16S batch effect beta-diversity plots
s16_batchBefore<-
  s16_kraken%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>0)>0.05)%>%
  select(any_of(colnames(s16_scrubbed)))%>%
  t()%>%
  avgdist(sample = 200,iterations=50)
s16_batchBefore_cmdscale<-cmdscale(s16_batchBefore,list. = T,eig=T,k=10)
s16_batchBefore_cmdscale$eig<-(s16_batchBefore_cmdscale$eig/sum(s16_batchBefore_cmdscale$eig))*100
s16_batchBefore_plot<-
  s16_batchBefore_cmdscale$points%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("SampleID")%>%
  left_join(s16_clinical_annos,by=c("SampleID"))%>%
  ggplot(aes(x=PC1,y=PC2,group=PCR_Plate,label=NA))+
  stat_ellipse(aes(color=PCR_Plate),show.legend = F)+
  geom_point(aes(fill=PCR_Plate),shape=21,color="gray80")+
  labs(x=paste0("PC1 (",round(s16_batchAfter_cmdscale$eig[1],1),"%)"),
       y=paste0("PC2 (",round(s16_batchAfter_cmdscale$eig[2],1),"%)"),
       fill="PCR Plate")+
  scale_fill_manual(values=rev(SBScolor))+
  scale_color_manual(values=rev(SBScolor))+
  theme(legend.text = element_text(size=8))+
  guides(fill=guide_legend(ncol=2))
s16_batchBefore_plot

s16_batchAfter<-
  s16_scrubbed%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>0)>0.05)%>%
  select(any_of(colnames(s16_scrubbed)))%>%
  t()%>%
  avgdist(sample = 100,iterations=50)
s16_batchAfter_cmdscale<-cmdscale(s16_batchAfter,list. = T,eig=T,k=10)
s16_batchAfter_cmdscale$eig<-(s16_batchAfter_cmdscale$eig/sum(s16_batchAfter_cmdscale$eig))*100
s16_batchAfter_plot<-
  s16_batchAfter_cmdscale$points%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("SampleID")%>%
  left_join(s16_clinical_annos,by=c("SampleID"))%>%
  ggplot(aes(x=PC1,y=PC2,group=PCR_Plate,label=NA))+
  stat_ellipse(aes(color=PCR_Plate),show.legend = F)+
  geom_point(aes(fill=PCR_Plate),shape=21,color="gray80")+
  labs(x=paste0("PC1 (",round(s16_batchAfter_cmdscale$eig[1],1),"%)"),
       y=paste0("PC2 (",round(s16_batchAfter_cmdscale$eig[2],1),"%)"),
       fill="PCR Plate")+
  scale_fill_manual(values=rev(SBScolor))+
  scale_color_manual(values=rev(SBScolor))+
  theme(legend.text = element_text(size=8))+
  guides(fill=guide_legend(ncol=2))
s16_batchAfter_plot
