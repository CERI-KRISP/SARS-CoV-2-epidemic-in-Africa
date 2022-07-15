#######################################################################################################
###########        Analysis of first 75,000 SARS-CoV-2 Sequences from Africa     ######################
###########     Author: Houriiyah Tegally, CERI/KRISP, Stellenbosch Uni & UKZN   ######################
#######################################################################################################

#Please also note that the genomic metadata shared in this training
#belong to the individual data generators in each African country and
#should only be used while complying to the GISAID terms and conditions

library(ggplot2)
library("readxl")
library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library('gridExtra')
library('data.table')
library('scales')

library(sf)
library(raster)
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
library(ggnewscale)
library(scatterpie)

#####################################
############# FIG 2 #################
#####################################

### Loading Data and Data Curation

global_metadata<-read.csv('metadata_2022-03-30_23-13.tsv',sep = "\t")
### have to select only humans and only original passage
africa_metadata<-subset(global_metadata, region=='Africa')
africa_metadata<-subset(africa_metadata, country!='Tanzania')

df_africa<-africa_metadata
df_africa$date<-as.Date(df_africa$date)

df_africa$days<-as.Date(cut(df_africa$date,breaks = "day",start.on.monday = FALSE))
df_africa$date2<-as.Date(cut(df_africa$date,breaks = "2 weeks",start.on.monday = FALSE))
df_africa$date3<-as.Date(cut(df_africa$date,breaks = "1 month",start.on.monday = FALSE))
df_africa$date<-as.Date(cut(df_africa$date,breaks = "1 week",start.on.monday = FALSE))


df_africa$date_submitted<-as.Date(df_africa$date_submitted)

df_africa$days_submitted<-as.Date(cut(df_africa$date_submitted,breaks = "day",start.on.monday = FALSE))


df_africa$submission_lag=df_africa$days_submitted-df_africa$days



df_africa['country'][df_africa['country'] == "Republic of the Congo"] <- "Congo"
df_africa['country'][df_africa['country'] == "Democratic Republic of the Congo"] <- "DRC"




owid_data<-read_excel('owid-covid-data_1April2022.xlsx')
owid_data$date<-as.Date(owid_data$date)
owid_data$date2<-as.Date(cut(owid_data$date,breaks = "1 week",start.on.monday = FALSE))


owid_data['location'][owid_data['location'] == "Democratic Republic of Congo"] <- "DRC"

owid_data_africa<-subset(owid_data, location=='Africa')



owid_data$location[owid_data$location == "Cabo Verde"] <- "Cape Verde"
owid_data$location[owid_data$location == "Côte d'Ivoire"] <- "Cote d'Ivoire"
owid_data$location[owid_data$location == "Union of the Comoros"] <- "Comoros"


df_africa$country[df_africa$country == "Cabo Verde"] <- "Cape Verde"
df_africa$country[df_africa$country == "Côte d'Ivoire"] <- "Cote d'Ivoire"
df_africa$country[df_africa$country == "Union of the Comoros"] <- "Comoros"

### PANEL D

country_list<-unique(df_africa$country)
df_sequencing_stats <- matrix(ncol=5, nrow=0)
colnames(df_sequencing_stats) <- c("country", "epiweeks_genomes","epiweeks_cases", "total_genomes","total_cases")
for(i in 1:length(unique(df_africa$country))){
  temp_df<-subset(df_africa,country==country_list[i])
  temp_df_2<-subset(owid_data,location==country_list[i])
  country = country_list[i]
  epiweeks_genomes = length(unique(temp_df$date))
  epiweeks_cases = length(unique(temp_df_2$date2))
  total_genomes<-length(temp_df$strain)
  total_cases<-max(subset(temp_df_2, !is.na(total_cases))$total_cases)
  #z = 3
  df_sequencing_stats = rbind(df_sequencing_stats, data.frame(country, epiweeks_genomes,epiweeks_cases,total_genomes,total_cases))
  
}

df_sequencing_stats

df_sequencing_stats$percentage_cases_sequenced<-(df_sequencing_stats$total_genomes/df_sequencing_stats$total_cases)*100
df_sequencing_stats$proportion_weeks_sequenced<-df_sequencing_stats$epiweeks_genomes/df_sequencing_stats$epiweeks_cases

df_sequencing_stats


df_sequencing_stats = df_sequencing_stats %>% 
  left_join(RCC_countries, by = c("country" = "country"))



p_sequencing_regularity<-ggplot(data=df_sequencing_stats,aes(x=as.numeric(proportion_weeks_sequenced),y=as.numeric(percentage_cases_sequenced)))+
  theme_minimal()+
  geom_point(aes(color=log10(as.numeric(total_cases)), size=log10(total_genomes)), alpha=0.8)+
  scale_color_distiller(palette='YlGnBu', direction=1,trans = 'log10', name='Total Cases', 
                        guide = guide_colourbar(frame.colour = "white", 
                                                ticks.colour = "white", # you can also remove the ticks with NA
                                                barwidth=4, barheight=0.8,
                                                direction = "horizontal",
                                                title.position = "top"),
                        
                        breaks=c(4,6), labels=c('10,000','1,000,000'))+
  geom_text_repel(label=df_sequencing_stats$country, size=4)+
  scale_size_continuous(name='Total Genomes', breaks=c(2,3,4), labels=c('100','1000','10,000'))+
  ylab('Percentage of cases sequenced overall (log scale)')+
  xlab('Proportion of epiweeks with sequencing data')+
  #labs(size='Total Genomes (Log)')+
  theme(legend.position = c(0.9,0.25))+
  #theme(legend.direction = 'horizontal')+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(axis.text=element_text(size=12))+
  theme(axis.title = element_text(size=14))+
  scale_y_continuous(trans = 'log10', breaks=c(0.01,0.1,1,7))+
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
  annotate("segment", x = 0.9, xend = 1, y = 12, yend = 12,
           colour = "grey40", size = 1, arrow = arrow(length = unit(.2,"cm"))) +
  annotate("segment", x = 0.1, xend = 0, y = 12, yend = 12,
           colour = "grey40", size = 1, arrow = arrow(length = unit(.2,"cm")))+
  annotate("text",label='Frequent Sequencing',x=0.8,y=14, size=5)+
  annotate("text",label='Sparse Sequencing',x=0.2,y=14, size=5)

p_sequencing_regularity




##### Additional data curation

df_sequencing_stats['country'][df_sequencing_stats['country'] == "Congo"] <- "Republic of Congo"
df_sequencing_stats['country'][df_sequencing_stats['country'] == "DRC"] <- "Democratic Republic of the Congo"
df_sequencing_stats['country'][df_sequencing_stats['country'] == "Cote d'Ivoire"] <- "Ivory Coast"

df_africa['country'][df_africa['country'] == "Congo"] <- "Republic of Congo"
df_africa['country'][df_africa['country'] == "Republic of the Congo"] <- "Republic of Congo"

df_africa['country'][df_africa['country'] == "DRC"] <- "Democratic Republic of the Congo"
df_africa['country'][df_africa['country'] == "Cote d'Ivoire"] <- "Ivory Coast"
df_africa['country'][df_africa['country'] == "Eswatini"] <- "Swaziland"


### Checking if we have info for all submitting labs
submitting_labs_info<-read_excel('submitting_lab_info_1APRIL.xlsx')
submitting_labs_current<-data.frame(unique(df_africa$submitting_lab))
names(submitting_labs_current)<-c('submitting_lab')
submitting_labs_info_prelim_15Feb<-right_join(submitting_labs_info,submitting_labs_current)
unique(subset(submitting_labs_info_prelim_15Feb,is.na(SequencingCountry)))$submitting_lab ### Some not able to get info
#write.csv(submitting_labs_info_prelim_15Feb,'submitting_labs_info_prelim_1April_v1.csv')

submission_locations<-read_excel('submitting_lab_info_1APRIL.xlsx')
submission_locations['SequencingCountry'][submission_locations['SequencingCountry'] == "DRC"] <- "Democratic Republic of the Congo"
submission_locations['SequencingCountry'][submission_locations['SequencingCountry'] == "Congo"] <- "Republic of Congo"


#df_africa <- merge(df_africa,submission_locations,by=c("submitting_lab"="submitting_lab"),all.x = TRUE)
df_africa <- unique(merge(df_africa,submission_locations,by=c("submitting_lab"="submitting_lab"),all.x = TRUE))
names(df_africa)[names(df_africa) == 'duplicate tracking (will give a uniform name for submitting lab'] <- 'SequencingLab'


sequencing_strategy = function(x,y,z)
{  
  
  if ( identical(x,y)) {
    return (as.character('Local'))
  } else if ( identical(z,"Africa")) {
    return (as.character('Africa'))
  } else {
    return (as.character('Outside Africa'))
  }
  
}

df_africa$SequencingStrategy<-mapply(sequencing_strategy, df_africa$country, df_africa$SequencingCountry, df_africa$SequencingContinent)


##### PANEL E


stat_box_data <- function(y, upper_limit = max(iris$Sepal.Length) * 1.15) {
  return( 
    data.frame(
      y = 270,
      label = paste('n =', length(y), '\n',
                    'median =', round(median(y), 1), '\n'
      )
    )
  )
}

my_comparisons = list(c("Africa","Local"), c("Local","Outside Africa"), c('Africa','Outside Africa'))

p_seq_location<-ggplot(data=subset(df_africa, !is.na(SequencingStrategy) & submission_lag<366), aes(x=reorder(SequencingStrategy,as.numeric(submission_lag),mean),y=as.numeric(submission_lag), fill=SequencingStrategy))+
  theme_classic()+
  geom_boxplot()+
  #stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "black" )+
  scale_fill_brewer(palette = 'Spectral', name='Sequencing Location')+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  ylab('Days from Collection to Submission')+
  xlab('Sequencing Strategy')+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text",
    size=3.5,
    #hjust = 0.5,
     vjust = 0.2
  )+
  theme(axis.text=element_text(size=12))+
  theme(axis.title = element_text(size=14))+
  theme(legend.position = 'none')+
  coord_flip()
p_seq_location

#### Panel A inset graph

df_count_seq_lab_all<-subset(df_africa,!is.na(submission_lag)) %>%
  group_by(country) %>%
  mutate(seq_lab_count=length(unique(SequencingLab)),submission_lag_mean=mean(submission_lag))

df_count_seq_lab_local<-subset(subset(df_africa,!is.na(submission_lag)),SequencingStrategy=="Local") %>%
  group_by(country,SequencingStrategy) %>%
  mutate(seq_lab_count=length(unique(SequencingLab)),submission_lag_mean=mean(submission_lag))

mydf<-subset(df_africa, SequencingStrategy=='Local')
mydf = mydf[order(mydf$date3, mydf$country), ]
mydf$country_id = as.integer(factor(mydf$country, levels = unique(mydf$country)))
mydf$cum_n_country = cummax(mydf$country_id)

df_count_seq_country_local<-mydf %>%
  group_by(date3) %>%
  summarize(cum_n_country = max(cum_n_country)) 


df_sequencing_regularity2<-left_join(df_sequencing_regularity,df_count_seq_lab_all)
df_sequencing_regularity2<-df_sequencing_regularity2[!duplicated(df_sequencing_regularity2), ]
df_sequencing_regularity2


p_local_increase<- ggplot(data=df_count_seq_country_local,aes(x = date3, y = cum_n_country)) + 
  theme_minimal()+
  #geom_point(color='grey40',size=2)+
  geom_line(size=0.5,color='grey40')+
  ylab('')+
  xlab('')+
  #ggtitle("Number of Countries with Local Sequencing")+
  scale_y_continuous(breaks=c(10,20,30),position='right')+
  scale_x_date(date_labels = "%b\n%y",date_breaks = "9 month")+
  theme(plot.title = element_text(hjust=0.5))

p_local_increase


#### Panel F

cortest2df <- function(cortest,x,y) {
  #cortest=cor.test(as.numeric(first_VOC_country_table_gamma$Gamma_arrival), first_VOC_country_table_gamma$total, method = "pearson", use = "complete.obs")
  #x="Gamma"
  #y="Sampling"
  temp_df <- matrix(ncol=6, nrow=0)
  
  a=cortest$estimate
  b=cortest$conf.int[1]
  c=cortest$conf.int[2]
  
  temp_df<-data.frame(x, y,a,b,c)
  colnames(temp_df) <- c("test","method","correlation", "upper","lower")
  rownames(temp_df) <- NULL
  
  return(temp_df)
}	

#number of seq labs vs proportion of weeks sequenced
cor_labs_weeks_sequenced_p<-cortest2df(cor.test(method='pearson',df_sequencing_regularity2$seq_lab_count,as.numeric(df_sequencing_regularity2$proportion_weeks_sequenced), use = "complete.obs"),"No of Labs vs Weeks Sequenced","pearson")
#cor_labs_weeks_sequenced_s<-cortest2df(cor.test(method='spearman',df_sequencing_regularity2$seq_lab_count,as.numeric(df_sequencing_regularity2$proportion_weeks_sequenced), use = "complete.obs"),"No of Labs vs Weeks Sequenced","spearman")

cor_labs_total_genomes_p<-cortest2df(cor.test(method='pearson',df_sequencing_regularity2$seq_lab_count,as.numeric(df_sequencing_regularity2$total_genomes), use = "complete.obs"),"No of Labs vs Total Genomes","pearson")
cor.test(method='spearman',df_sequencing_regularity2$seq_lab_count,as.numeric(df_sequencing_regularity2$total_genomes))

cor_labs_submission_lag_p<-cortest2df(cor.test(method='pearson',subset(df_sequencing_regularity2,!is.na(as.numeric(submission_lag_mean)))$seq_lab_count,as.numeric(subset(df_sequencing_regularity2,!is.na(as.numeric(submission_lag_mean)))$submission_lag_mean), use = "complete.obs"),"No of Labs vs Sequencing\nTurn-around time","pearson")
cor.test(method='spearman',subset(df_sequencing_regularity2,!is.na(as.numeric(submission_lag_mean)))$seq_lab_count,as.numeric(subset(df_sequencing_regularity2,!is.na(as.numeric(submission_lag_mean)))$submission_lag_mean))



all<-rbind(cor_labs_weeks_sequenced_p,cor_labs_total_genomes_p,cor_labs_submission_lag_p)
#alpha_travel_correlations<-rbind(corr_alpha_sampling_uk_p,corr_alpha_inferred_uk_p)


num_labs_corr<-ggplot(all)+
  theme_bw()+
  geom_vline(xintercept = 0,linetype='dashed',size=1,color='grey')+
  
  geom_point(aes(x=abs(correlation),y=test),size=4, shape=21,color='grey40',fill='white', stroke=1)+
  #geom_segment(aes(x=abs(lower),y=test,xend=abs(upper),yend=test), color='darkseagreen4', size=2)+
  
  ylab("")+
  xlab("Strength of Association")+
  xlim(0,1)+
  ggtitle("No. of Sequencing Labs vs Sequencing Output")+
  scale_y_discrete(position = "right")+
  theme(plot.title = element_text(hjust=0.5))
num_labs_corr


##### PANEL C

africa_cdc_map<-read.csv('AfricaCDC_shapefile/Africa_fortified_Name_0.tsv',sep = "\t")
africa_cdc_map['id'][africa_cdc_map['id'] == "Côte d'Ivoire"] <- "Ivory Coast"

world3 <- map_data("world")

world4<- world3 %>% 
  left_join(df_sequencing_stats, by = c("region" = "country"))


world5<- africa_cdc_map %>% 
  left_join(df_sequencing_stats, by = c("id" = "country"))

seq_location<-as.data.frame.matrix(table(df_africa %>% dplyr::select("country","SequencingStrategy")))
seq_location <- cbind(rownames(seq_location), seq_location)
rownames(seq_location) <- NULL
colnames(seq_location) <- c("country","Africa","Local",'Outside Africa')
seq_location<- seq_location %>% 
  left_join(world4, by = c("country" = "region"))



scatterpie_data <- seq_location %>%
  group_by(country,total_genomes) %>%
  summarise(long = mean(as.numeric(long)),
            lat = mean(as.numeric(lat)))


final_data <- left_join(scatterpie_data, seq_location, by = "country") 

final_data <- unique( final_data[ , 1:7 ] )

###Manually changing the number breakdowns for Egypt, Algeria, Tunisia, Gabon and Ghana to reflect ACEGID sequencing
###(Not correctly captured in the submitting lab column on GISAID as the originating labs submitted the data themselves)
final_data[final_data$country=='Egypt', "Africa"] <- final_data[final_data$country=='Egypt', "Africa"]+15
final_data[final_data$country=='Egypt', "Local"] <- final_data[final_data$country=='Egypt', "Local"]-15

final_data[final_data$country=='Algeria', "Africa"] <- final_data[final_data$country=='Algeria', "Africa"]+50
final_data[final_data$country=='Algeria', "Local"] <- final_data[final_data$country=='Algeria', "Local"]-50

final_data[final_data$country=='Tunisia', "Africa"] <- final_data[final_data$country=='Tunisia', "Africa"]+50
final_data[final_data$country=='Tunisia', "Local"] <- final_data[final_data$country=='Tunisia', "Local"]-50


final_data[final_data$country=='Gabon', "Africa"] <- final_data[final_data$country=='Gabon', "Africa"]+77
final_data[final_data$country=='Gabon', "Local"] <- final_data[final_data$country=='Gabon', "Local"]-77


final_data[final_data$country=='Ghana', "Africa"] <- final_data[final_data$country=='Ghana', "Africa"]+50
final_data[final_data$country=='Ghana', "Local"] <- final_data[final_data$country=='Ghana', "Local"]-50


###Manually changing the coordinates of the pie-chart for South Africa to not overlap with Lesotho's
final_data[final_data$country=='South Africa', "long.x"] <- 22
final_data[final_data$country=='South Africa', "lat.x"] <- -31

map_mix_strategy<-ggplot(world5, aes(long, lat)) +
  theme_void()+
  geom_map(map=world5, aes(map_id=id, fill=total_genomes), color="grey40",size=0.2) +
  scale_fill_distiller(palette = "PuBuGn", direction = 1,na.value = "white",trans = "log",
                       breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000), name='Genomes') + 
  
  new_scale("fill") +
  
  geom_scatterpie(data = subset(final_data,!is.na(long.x)), aes(x=long.x, y=lat.x, group=country, r=2), 
                                   cols=c("Africa", "Local", "Outside Africa"),legend_name ='Sequencing',size=0.1)+
  scale_fill_brewer(palette = 'Spectral')+
  
  coord_fixed()+
  scale_x_continuous(limits = c(-26, 60))+
  scale_y_continuous(limits = c(-40, 38))
#map_mix_strategy

#ggsave("Africa3.pdf", width = 50, height = 40, units = "cm", limitsize = FALSE)


##### PANEL A

df_seq_locations2<-read_excel('local_sequencing_data.xlsx')



map_local_sequencing<-ggplot() +
  theme_void()+
  geom_map(data=world5,map=world5, aes(long, lat,map_id=id),fill='white', color="grey40",size=0.2)+
  geom_map(data=subset(world5,id=='South Africa'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Lesotho'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='white',size=0.2)+
  
   geom_map(data=subset(world5,id=='Nigeria'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Ghana'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Senegal'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Uganda'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Republic of Congo'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Democratic Republic of the Congo'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Zimbabwe'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  
  geom_map(data=subset(world5,id=='Kenya'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Gambia'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Sudan'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Algeria'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Guinea'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Botswana'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Algeria'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Burkina Faso'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Cameroon'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Central African Republic'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Egypt'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Ethiopia'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Reunion'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Gabon'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Reunion'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Libya'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Madagascar'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Malawi'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Mali'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Mauritius'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Morocco'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Rwanda'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Sierra Leone'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Tunisia'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Zambia'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Ivory Coast'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Mozambique'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Seychelles'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Niger'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Sao Tome and Principe'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Cape Verde'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  geom_map(data=subset(world5,id=='Mauritania'),map=world5, aes(long, lat,map_id=id), color="black", fill='grey60')+
  
  geom_point(data = subset(df_seq_locations2,!is.na(as.numeric(long_submitting_lab))),
             aes(x = as.numeric(long_submitting_lab), y = as.numeric(lat_submitting_lab)),
             size=2, color='black', shape=21, fill='red1')+
  #geom_label_repel(data = distinct(subset(df_seq_locations2,!is.na(as.numeric(long_submitting_country))),country,.keep_all = TRUE),
   #                aes(x = as.numeric(long_submitting_lab), y = as.numeric(lat_submitting_lab), label=stringr::str_wrap(country, 20)),
   #                size=2, color='black',label.padding = unit(0.1, "lines"),
   #                label.r = unit(0.15, "lines"),
    #               label.size = 0.1, alpha=0.7,nudge_y = -3)+
  scale_x_continuous(limits = c(-26, 60))+
  scale_y_continuous(limits = c(-40, 38))+
  coord_fixed()

#map_local_sequencing



##### PANEL B


df_seq_locations3<-read_excel('Africa_sequencing_data.xlsx')

map_africa_sequencing_network<-ggplot() +
  theme_void()+
  geom_map(data=world5,map=world5, aes(long, lat,map_id=id), color="grey40", fill='white',size=0.2)+
  geom_map(data=world5,map=world5, aes(long, lat,map_id=id), color="grey40", fill='white',size=0.2)+
  
  geom_map(data=subset(world5,id=='Ethiopia'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Malawi'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Angola'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Zimbabwe'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Mozambique'),map=world5, aes(long, lat,map_id=id), color="grey40",fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Uganda'),map=world5, aes(long, lat,map_id=id), color="grey40",fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Namibia'),map=world5, aes(long, lat,map_id=id), color="grey40",fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Cameroon'),map=world5, aes(long, lat,map_id=id), color="grey40",fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Botswana'),map=world5, aes(long, lat,map_id=id), color="grey40",fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Mauritius'),map=world5, aes(long, lat,map_id=id), color="grey40",fill='blue4', alpha=0.5,size=0.2)+
  
  geom_map(data=subset(world5,id=='Sudan'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Swaziland'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Lesotho'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Swaziland'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  
  geom_map(data=subset(world5,id=='Burkina Faso'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Cameroon'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Rwanda'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Somalia'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Central African Republic'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Tunisia'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Algeria'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Egypt'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Gabon'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Ghana'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkgreen', alpha=0.5,size=0.2)+
  
  geom_map(data=subset(world5,id=='Togo'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='gold2', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Sierra Leone'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='gold2', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Benin'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='gold2', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Nigeria'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='gold2', alpha=0.5,size=0.2)+
  
  geom_map(data=subset(world5,id=='Tunisia'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='burlywood4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Cape Verde'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='burlywood4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Madagascar'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='burlywood4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Guinea'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='burlywood4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Central African Republic'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='burlywood4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Cameroon'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='burlywood4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Niger'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='burlywood4', alpha=0.5,size=0.2)+
  
  geom_map(data=subset(world5,id=='South Sudan'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='deeppink4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Burundi'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='deeppink4', alpha=0.5,size=0.2)+
  
  geom_map(data=subset(world5,id=='Republic of Congo'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkorange3', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Chad'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkorange3', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Central African Republic'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkorange3', alpha=0.5,size=0.2)+
  #geom_map(data=subset(world5,id=='Democratic Republic of the Congo'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkorange3', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Cameroon'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='darkorange3', alpha=0.5,size=0.2)+
  
  geom_map(data=subset(world5,id=='Seychelles'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='purple4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Comoros'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='purple4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Ethiopia'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='purple4', alpha=0.5,size=0.2)+
  geom_map(data=subset(world5,id=='Democratic Republic of the Congo'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='purple4', alpha=0.5,size=0.2)+
  
  geom_map(data=subset(world5,id=='Guinea-Bissau'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='black', alpha=0.5,size=0.2)+

  geom_map(data=subset(world5,id=='South Africa'),map=world5, aes(long, lat,map_id=id), color="black", fill='blue4')+
  geom_map(data=subset(world5,id=='Nigeria'),map=world5, aes(long, lat,map_id=id), color="black", fill='darkgreen')+
  geom_map(data=subset(world5,id=='Ghana'),map=world5, aes(long, lat,map_id=id), color="black", fill='gold2')+
  geom_map(data=subset(world5,id=='Senegal'),map=world5, aes(long, lat,map_id=id), color="black", fill='burlywood4')+
  geom_map(data=subset(world5,id=='Uganda'),map=world5, aes(long, lat,map_id=id), color="black", fill='deeppink4')+
  geom_map(data=subset(world5,id=='Democratic Republic of the Congo'),map=world5, aes(long, lat,map_id=id), color="black", fill='darkorange3')+
  geom_map(data=subset(world5,id=='Kenya'),map=world5, aes(long, lat,map_id=id), color="black", fill='purple4')+
  geom_map(data=subset(world5,id=='Gambia'),map=world5, aes(long, lat,map_id=id), color="black", fill='black')+
  geom_map(data=subset(world5,id=='Lesotho'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='white',size=0.2)+
  geom_map(data=subset(world5,id=='Lesotho'),map=world5, aes(long, lat,map_id=id), color="grey40", fill='blue4', alpha=0.5,size=0.2)+
  
  
  geom_curve(data = subset(df_seq_locations3,!is.na(long_sequencing_country)),
             aes(color=duplicate_submitting_lab2,x = long_originating_country, y = lat_originating_country, xend = long_sequencing_lab, yend = lat_sequencing_lab),
             size=0.7)+
  geom_point(data = subset(df_seq_locations3,!is.na(long_sequencing_country)),
             aes(x = long_sequencing_lab, y = lat_sequencing_lab),
             size=2, color='black', shape=21, fill='red1')+
  #geom_label(data = subset(subset(df_seq_locations3,!is.na(long_sequencing_country)),duplicate_submitting_lab2!='MRCG'),
   #          aes(x = long_sequencing_lab, y = lat_sequencing_lab, label=stringr::str_wrap(duplicate_submitting_lab2)),
   #          size=2, color='black',label.padding = unit(0.1, "lines"),
   #          label.r = unit(0.15, "lines"),
   #          label.size = 0.1, alpha=0.2,nudge_y = -3)+
  
  #geom_label(data = subset(subset(df_seq_locations3,!is.na(long_sequencing_country)),duplicate_submitting_lab2=='MRCG'),
  #           aes(x = long_sequencing_lab, y = lat_sequencing_lab, label=stringr::str_wrap(duplicate_submitting_lab2)),
   #          size=2, color='black',label.padding = unit(0.1, "lines"),
   #          label.r = unit(0.15, "lines"),
    #         label.size = 0.1, alpha=0.2,nudge_y = -4)+
  
  scale_color_manual(values=c('springgreen4','dodgerblue1', 'plum4','darkorange2','wheat4','purple3','deeppink3','black','blue1','gold3'))+
  theme(legend.position = 'none')+
  scale_x_continuous(limits = c(-26, 60))+
  scale_y_continuous(limits = c(-40, 38))+
  coord_fixed()
  
#map_africa_sequencing_network

temp<-plot_grid(map_local_sequencing,map_africa_sequencing_network,map_mix_strategy,ncol=3,rel_widths = c(0.3,0.3,0.4))


ggsave(file="fig2_maps.pdf",plot=temp, width = 25, height = 10, units = "cm", limitsize = FALSE)



