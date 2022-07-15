#######################################################################################################
###########        Analysis of first 100,000 SARS-CoV-2 Sequences from Africa     ######################
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
############# FIG 1 #################
#####################################

global_metadata<-read.csv('metadata_2022-03-30_23-13.tsv',sep = "\t")

### have to select only humans and only original passage
global_metadata<-subset(global_metadata, host=='Human')
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


#count_prop_variant<-data.frame(prop.table(table(df_africa$Nextstrain_clade)))

country_South<-c('South Africa', 'Botswana','Zambia')
RCC_South<-c('South','South','South')

country_East<-c('Kenya','Uganda', 'Ethiopia')
RCC_East<-c('East','East','East')


country_West<-c('Nigeria','Ghana', 'Senegal')
RCC_West<-c('West','West','West')


country_Central<-c('DRC','Congo', 'Cameroon')
RCC_Central<-c('Central','Central','Central')

country_North<-c('Egypt','Morocco','Tunisia')
RCC_North<-c('North','North','North')



country<-c(country_South,country_East,country_West,country_Central,country_North)
RCC<-c(RCC_South,RCC_East,RCC_West,RCC_Central,RCC_North)

RCC_countries <- data.frame(country, RCC)

df_africa = df_africa %>% 
  left_join(RCC_countries, by = c("country" = "country"))

df_africa_second = df_africa_second %>% 
  left_join(RCC_countries, by = c("country" = "country"))

owid_data<-read_excel('owid-covid-data_1April2022.xlsx')
owid_data$date<-as.Date(owid_data$date)

#owid_data[owid_data == "Republic of Congo"] <- "Congo"
owid_data['location'][owid_data['location'] == "Democratic Republic of Congo"] <- "DRC"


owid_data <- owid_data %>% 
  dplyr::mutate(new_cases_7 = zoo::rollmean(new_cases, k = 7, fill = NA)) %>% 
  dplyr::ungroup()

owid_data_africa<-subset(owid_data, location=='Africa')





owid_data = owid_data %>% 
  left_join(RCC_countries, by = c("location" = "country"))

owid_data$date2<-as.Date(cut(owid_data$date,breaks = "1 week",start.on.monday = FALSE))

owid_data$location[owid_data$location == "Cabo Verde"] <- "Cape Verde"
owid_data$location[owid_data$location == "Côte d'Ivoire"] <- "Cote d'Ivoire"
owid_data$location[owid_data$location == "Union of the Comoros"] <- "Comoros"


epi_curve<-ggplot(data=owid_data_africa, aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='darkolivegreen3', alpha=0.3) +
  
  geom_line(color='darkolivegreen3', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  ggtitle("Africa Total")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.1,40))
epi_curve

panelD<-ggplot(data=df_africa)  + theme_minimal_hgrid()+
   geom_count(data=subset(df_africa, Nextstrain_clade %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(df_africa, Nextstrain_clade %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(df_africa, Nextstrain_clade %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(df_africa, Nextstrain_clade %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  
  geom_count(data=subset(df_africa, Nextstrain_clade %like% 'Eta'),aes(fill='Eta',x=date3, y='Eta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(df_africa, !Nextstrain_clade %like% 'Beta' & !Nextstrain_clade %like% 'Alpha' & !Nextstrain_clade %like% 'Delta' & !Nextstrain_clade %like% 'Omicron' & !Nextstrain_clade %like% 'Eta'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  scale_size(range=c(0,15), limits = c(5,10000), breaks=c(50,500,5000))+
  ylab('')+ xlab('month')+
  labs(size='Genomes Sampled Monthly')+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4","#2B83BA",'hotpink2','grey'), name='')+
  theme(legend.position="bottom") +
  theme(legend.justification = 0.5) +
  theme(legend.direction="horizontal") +
  theme(legend.box="vertical") +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=10))+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  guides(fill = guide_legend(override.aes = list(size=4),nrow=1, order=2), size=guide_legend(order=1))+
  xlab('Sampling Dates')

panelD


Africa_fig<-plot_grid(epi_curve+theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.x = element_blank()),
                      panelD,ncol=1,align='v', rel_heights = c(0.35,0.65))+
  theme(plot.background = element_rect(fill = 'grey95'))
Africa_fig

epi_curve_SA<-ggplot(data=subset(owid_data, RCC=='South'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='grey70', alpha=0.3) +
  
  geom_line(color='grey70', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))

epi_curve_SA

genomes_SA<-ggplot(data=subset(df_africa, RCC=='South'))  + theme_minimal_hgrid()+
  geom_count(data=subset(subset(df_africa, RCC=='South'), Nextstrain_clade %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='South'), Nextstrain_clade %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='South'), Nextstrain_clade %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='South'), Nextstrain_clade %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
 geom_count(data=subset(subset(df_africa, RCC=='South'), pango_lineage=='C.1.2'),aes(fill='C.1.2',x=date3, y='C.1.2'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  
  geom_count(data=subset(subset(df_africa, RCC=='South'), !Nextstrain_clade %like% 'Beta' & !Nextstrain_clade %like% 'Alpha' & !Nextstrain_clade %like% 'Delta' & !Nextstrain_clade %like% 'Omicron' & pango_lineage!='C.1.2'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_size(range=c(0,15), limits = c(0,5000), breaks=c(5,50,500,5000))+
  
  theme(legend.position="bottom") +
  theme(legend.justification = 0.5) +
  
  theme(legend.direction="horizontal") +
  theme(legend.box="vertical") +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=8)) +
  theme(axis.title.x = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.title.y = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF",'indianred3',"#ABDDA4",'hotpink2','grey'), name='Lineages', guide='none')+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())

genomes_SA

SA_fig<-plot_grid(epi_curve_SA+theme(axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank()),
                  genomes_SA,ncol=1,align='v', rel_heights = c(0.32,0.65))
SA_fig


epi_curve_East<-ggplot(data=subset(owid_data, RCC=='East'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='indianred3', alpha=0.3) +
  
  geom_line(color='indianred3', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))
epi_curve_East

genomes_East<-ggplot(data=subset(df_africa, RCC=='East'))  + theme_minimal_hgrid()+
 geom_count(data=subset(subset(df_africa, RCC=='East'), Nextstrain_clade %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='East'), Nextstrain_clade %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='East'), Nextstrain_clade %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='East'), Nextstrain_clade %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
    geom_count(data=subset(subset(df_africa, RCC=='East'), pango_lineage=='A.23.1'),aes(fill='A.23.1',x=date3, y='A.23.1'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  
  geom_count(data=subset(subset(df_africa, RCC=='East'), !Nextstrain_clade %like% 'Beta' & !Nextstrain_clade %like% 'Alpha' & !Nextstrain_clade %like% 'Delta' & !Nextstrain_clade %like% 'Omicron' & pango_lineage!='A.23.1'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_fill_manual(values=c('indianred3',"#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name='Lineages', guide='none')+
  scale_size(range=c(0,15), limits = c(0,1000), breaks=c(5,50,500))+
  
  theme(legend.position="bottom") +
  theme(legend.justification = 0.5) +
  
  theme(legend.direction="horizontal") +
  theme(legend.box="vertical") +
  #theme(legend.box.just = 0.7) +
  
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.title.y = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())


genomes_East

East_fig<-plot_grid(epi_curve_East+theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.title.x = element_blank()),
                    genomes_East,ncol=1,align='v', rel_heights = c(0.35,0.65))
East_fig



epi_curve_West<-ggplot(data=subset(owid_data, RCC=='West'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='grey20', alpha=0.3) +
  
  geom_line(color='grey20', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))
epi_curve_West

genomes_West<-ggplot(data=subset(df_africa, RCC=='West'))  + theme_minimal_hgrid()+
 geom_count(data=subset(subset(df_africa, RCC=='West'), Nextstrain_clade %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='West'), Nextstrain_clade %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='West'), Nextstrain_clade %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='West'), Nextstrain_clade %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  
  geom_count(data=subset(subset(df_africa, RCC=='West'), Nextstrain_clade %like% 'Eta'),aes(fill='Eta',x=date3, y='Eta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='West'), !Nextstrain_clade %like% '20H (Beta, V2)' & !Nextstrain_clade %like% '20I (Alpha, V1)' & !Nextstrain_clade %like% 'Delta' & !Nextstrain_clade %like% 'Omicron' & !Nextstrain_clade %like% 'Eta'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4","#2B83BA",'hotpink2','grey'), name='Lineages', guide='none')+
  scale_size(range=c(0,15), limits = c(0,500), breaks=c(5,50,500))+
  
  theme(legend.position="bottom") +
  theme(legend.justification = 0.5) +
  
  theme(legend.direction="horizontal") +
  theme(legend.box="vertical") +
  #theme(legend.box.just = 0.7) +
  
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.title.y = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())


genomes_West

West_fig<-plot_grid(epi_curve_West+theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.title.x = element_blank()),
                    genomes_West,ncol=1,align='v', rel_heights = c(0.35,0.65))
West_fig



epi_curve_Central<-ggplot(data=subset(owid_data, RCC=='Central'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='green4', alpha=0.3) +
  
  geom_line(color='green4', size=0.8)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))
epi_curve_Central

genomes_Central<-ggplot(data=subset(df_africa, RCC=='Central'))  + theme_minimal_hgrid()+
 geom_count(data=subset(subset(df_africa, RCC=='Central'), Nextstrain_clade %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='Central'), Nextstrain_clade %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='Central'), Nextstrain_clade %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='Central'), Nextstrain_clade %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='Central'), pango_lineage=='B.1.620'),aes(fill='B.1.620',x=date3, y='B.1.620'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  
  geom_count(data=subset(subset(df_africa, RCC=='Central'), !Nextstrain_clade %like% 'Beta' & !Nextstrain_clade %like% 'Alpha' & !Nextstrain_clade %like% 'Delta' & !Nextstrain_clade %like% 'Delta'  & pango_lineage!='B.1.620'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
 
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_fill_manual(values=c("#FDAE61",'indianred3',"#FFFFBF","#ABDDA4",'hotpink2','grey'), name='Lineages', guide='none')+
  scale_size(range=c(0,15), limits = c(0,500), breaks=c(5,50,500))+
  
  theme(legend.position="bottom") +
  theme(legend.justification = 0.5) +
  
  theme(legend.direction="horizontal") +
  theme(legend.box="vertical") +
  #theme(legend.box.just = 0.7) +
  
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.title.y = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+

  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())


genomes_Central

Central_fig<-plot_grid(epi_curve_Central+theme(axis.text.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.title.x = element_blank()),
                       genomes_Central,ncol=1,align='v', rel_heights = c(0.35,0.65))
Central_fig




epi_curve_North<-ggplot(data=subset(owid_data, RCC=='North'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  #geom_area(fill='tan3', alpha=0.3) +
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='tan3', alpha=0.3) +
  
  geom_line(color='tan3', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))
epi_curve_North

genomes_North<-ggplot(data=subset(df_africa, RCC=='North'))  + theme_minimal_hgrid()+
geom_count(data=subset(subset(df_africa, RCC=='North'), Nextstrain_clade %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='North'), Nextstrain_clade %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='North'), Nextstrain_clade %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='North'), Nextstrain_clade %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa, RCC=='North'), pango_lineage=='C.36.3' | pango_lineage=="C.36"),aes(fill='C.36',x=date3, y='C.36'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  
  geom_count(data=subset(subset(df_africa, RCC=='North'), !Nextstrain_clade %like% 'Beta' & !Nextstrain_clade %like% 'Alpha' & !Nextstrain_clade %like% 'Delta'  & !Nextstrain_clade %like% 'Omicron' & pango_lineage!='C.36' & pango_lineage!='C.36.3'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  #geom_count() +
  #scale_size_area()+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF",'indianred3',"#ABDDA4",'hotpink2','grey'), name='Lineages', guide='none')+
  scale_size(range=c(0,15), limits = c(0,500), breaks=c(5,50,500))+
  
  theme(legend.position="bottom") +
  theme(legend.justification = 0.5) +
  
  theme(legend.direction="horizontal") +
  theme(legend.box="vertical") +
  #theme(legend.box.just = 0.7) +
  
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.title.y = element_text(color="black", size=8, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())


genomes_North

North_fig<-plot_grid(epi_curve_North+theme(axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.title.x = element_blank()),
                     genomes_North,ncol=1,align='v', rel_heights = c(0.35,0.65))
North_fig

a<-plot_grid(Africa_fig,North_fig,SA_fig,ncol=3,labels = c('A', 'B','C'))
b<-plot_grid(West_fig,Central_fig,East_fig,ncol=3,labels = c('D', 'E','F'))


wholegfig<-plot_grid(a,b, ncol=1,align='v',rel_heights = c(0.5,0.5))

ggsave('Fig1.pdf',plot=wholegfig, width = 45, height = 28, units = "cm",limitsize = FALSE)
