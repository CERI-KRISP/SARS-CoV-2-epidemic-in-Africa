library(ggplot2)
library(lubridate)
library(readxl)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library(gridExtra)
library(data.table)
library(scales)
library(openxlsx)
library(readr)
library(zoo)
library(ggforestplot)
library(sf)
library(raster)
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
library(ggnewscale)



#####Supplementary Figure S1

##Read in genomic data

df_africa <- read_csv("africa_100k.csv")

##Format date

df_africa$date<-as.Date(df_africa$Collection.date)
df_africa_second$days<-as.Date(cut(df_africa_second$date,breaks = "day",start.on.monday = FALSE))
df_africa_second$date2<-as.Date(cut(df_africa_second$date,breaks = "2 weeks",start.on.monday = FALSE))
df_africa_second$date3<-as.Date(cut(df_africa_second$date,breaks = "1 month",start.on.monday = FALSE))
df_africa_second$date<-as.Date(cut(df_africa_second$date,breaks = "1 week",start.on.monday = FALSE))

##Harmonise vountry names between datasets

df_africa$country[df_africa$country == "Cabo Verde"] <- "Cape Verde"
df_africa$country[df_africa$country == "CÃ´te d'Ivoire"] <- "Cote d'Ivoire"
df_africa$country[df_africa$country == "Union of the Comoros"] <- "Comoros"
df_africa$country[df_africa$country == "Central African Republic"] <- "CAR"

##Create country lists and join with dataframe

country_South<-c('Angola', 'Eswatini','Lesotho', 'Malawi', 'Mozambique', 'Namibia', 'Zimbabwe')
RCC_South<-c('South','South','South','South','South','South', 'South')
country_East<-c('Comoros','Djibouti', 'Madagascar', 'Rwanda', 'Somalia', 'South Sudan', 'Sudan')
RCC_East<-c('East','East','East','East','East','East','East')
country_Islands<-c('Seychelles', 'Mauritius', 'Cape Verde', 'Sao Tome and Principe')
RCC_Islands<-c('Islands', 'Islands', 'Islands', 'Islands')
country_West<-c('Benin','Burkina Faso', "Cote d'Ivoire", 'Gambia', 'Guinea', 'Guinea-Bissau')
RCC_West<-c('West','West','West','West','West', 'West')
country_West2<-c('Liberia','Mali','Niger','Nigeria','Sierra Leone', 'Togo')
RCC_West2<-c('West2','West2','West2','West2','West2', 'West2')
country_Central<-c('Burundi','Central African Republic', 'Chad', 'Equatorial Guinea','Gabon')
RCC_Central<-c('Central','Central','Central', 'Central','Central')
country_North<-c('Algeria','Libya')
RCC_North<-c('North','North')
country<-c(country_South,country_East,country_Islands,country_West,country_West2,country_Central,country_North)
RCC<-c(RCC_South,RCC_East,RCC_Islands,RCC_West,RCC_West2,RCC_Central,RCC_North)
RCC_countries <- data.frame(country, RCC)

df_africa = df_africa %>%
  left_join(RCC_countries, by = c("country" = "country"))

##Read in epi data and format date

owid_data<-read_excel('owid-covid-data.xlsx')
owid_data$date<-as.Date(owid_data$date)

owid_data = owid_data %>%
  left_join(RCC_countries, by = c("location" = "country"))

owid_data$location[owid_data$location == "Cabo Verde"] <- "Cape Verde"
owid_data$location[owid_data$location == "Côte d'Ivoire"] <- "Cote d'Ivoire"
owid_data$location[owid_data$location == "Union of the Comoros"] <- "Comoros"
owid_data$location[owid_data$location == "Central African Republic"] <- "CAR"

##Plot figures
#West1

epi_curve_West<-ggplot(data=subset(owid_data, RCC=='West'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='black', alpha=0.3) +
  geom_line(color='black', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))
  
  theme(axis.title.y = element_text(color="black", size=8, face="bold"))+
  geom_line(color='black', size=1)+
  xlab(" ")+
  ylab('New Cases per Million\n(Smoothed)')+
  facet_grid(.~location)
epi_curve_West

genomes_West<-ggplot(data=subset(df_africa_second, RCC=='West'))  + theme_minimal_hgrid()+
  geom_count(data=subset(subset(df_africa_second, RCC=='West'), Variant %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West'), Variant %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West'), Variant %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West'), Variant %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West'), Variant %like% 'Eta'),aes(fill='Eta',x=date3, y='Eta' ), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West'), !Variant %like% 'Beta' & !Variant %like% 'Alpha' & !Variant %like% 'Delta' & !Variant %like% 'Omicron' & !Variant %like% 'Eta'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_size(range=c(0,8), limits = c(0,500), breaks=c(5,50,500))+
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
  theme(axis.text.y = element_text(color="black", size=8))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y")+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'#2B83BA','hotpink2','grey'), name='Lineages', guide='none')+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())

genomes_West

West_fig<-plot_grid(epi_curve_West+theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.title.x = element_blank()),
                    genomes_West+theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5)),ncol=1, rel_heights = c(0.35,0.65))
West_fig


###West2 figure


epi_curve_West2<-ggplot(data=subset(owid_data, RCC=='West2'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='black', alpha=0.3) +
  geom_line(color='black', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))


epi_curve_West2

genomes_West2<-ggplot(data=subset(df_africa_second, RCC=='West2'))  + theme_minimal_hgrid()+
  geom_count(data=subset(subset(df_africa_second, RCC=='West2'), Variant %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West2'), Variant %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West2'), Variant %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West2'), Variant %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West2'), Variant %like% 'Eta'),aes(fill='Eta',x=date3, y='Eta' ), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='West2'), !Variant %like% 'Beta' & !Variant %like% 'Alpha' & !Variant %like% 'Delta' & !Variant %like% 'Omicron' & !Variant %like% 'Eta'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_size(range=c(0,8), limits = c(0,500), breaks=c(5,50,500))+
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
  theme(axis.text.y = element_text(color="black", size=8))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y")+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'#2B83BA','hotpink2','grey'), name='Lineages', guide='none')+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())

genomes_West2

West_fig2<-plot_grid(epi_curve_West2+theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.title.x = element_blank()),
                    genomes_West2+theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5)),ncol=1, rel_heights = c(0.35,0.65))
West_fig2


###South figure

epi_curve_South<-ggplot(data=subset(owid_data, RCC=='South'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='grey70', alpha=0.3) +
  geom_line(color='grey70', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))


epi_curve_South

genomes_South<-ggplot(data=subset(df_africa_second, RCC=='South'))  + theme_minimal_hgrid()+
  geom_count(data=subset(subset(df_africa_second, RCC=='South'), Variant %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='South'), Variant %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='South'), Variant %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='South'), Variant %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='South'), Variant %like% 'Eta'),aes(fill='Eta',x=date3, y='Eta' ), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='South'), !Variant %like% 'Beta' & !Variant %like% 'Alpha' & !Variant %like% 'Delta' & !Variant %like% 'Omicron' & !Variant %like% 'Eta'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_size(range=c(0,8), limits = c(0,500), breaks=c(5,50,500))+
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
  theme(axis.text.y = element_text(color="black", size=8))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y")+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'#2B83BA','hotpink2','grey'), name='Lineages', guide='none')+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())

genomes_South

South_fig<-plot_grid(epi_curve_South+theme(axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.title.x = element_blank()),
                     genomes_South+theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5)),ncol=1, rel_heights = c(0.35,0.65))
South_fig


###North figure

epi_curve_North<-ggplot(data=subset(owid_data, RCC=='North'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='tan3', alpha=0.3) +
  geom_line(color='tan3', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))


epi_curve_North

genomes_North<-ggplot(data=subset(df_africa_second, RCC=='North'))  + theme_minimal_hgrid()+
  geom_count(data=subset(subset(df_africa_second, RCC=='North'), Variant %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='North'), Variant %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='North'), Variant %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='North'), !Variant %like% 'Beta' & !Variant %like% 'Alpha' & !Variant %like% 'Delta' & !Variant %like% 'Omicron' & !Variant %like% 'Eta'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_size(range=c(0,8), limits = c(0,60), breaks=c(5,25,60))+
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
  theme(axis.text.y = element_text(color="black", size=8))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y")+
  scale_fill_manual(values=c("#FDAE61","#ABDDA4",'hotpink2','grey'), name='Lineages', guide='none')+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())

genomes_North

North_fig<-plot_grid(epi_curve_North+theme(axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.title.x = element_blank()),
                     genomes_North+theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5)),ncol=1, rel_heights = c(0.35,0.65))
North_fig

###Central figure

epi_curve_Central<-ggplot(data=subset(owid_data, RCC=='Central'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='green4', alpha=0.3) +
  geom_line(color='green4', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))

epi_curve_Central

genomes_Central<-ggplot(data=subset(df_africa_second, RCC=='Central'))  + theme_minimal_hgrid()+
  geom_count(data=subset(subset(df_africa_second, RCC=='Central'), Variant %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Central'), Variant %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Central'), Variant %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Central'), Variant %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Central'), Pango.lineage=='B.1.620'),aes(fill='B.1.620',x=date3, y='B.1.620'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Central'), ! Pango.lineage %like% 'B.1.620' & !Variant %like% 'Alpha' & !Variant %like% 'Delta' & !Variant %like% 'Omicron' & !Variant %like% 'Beta'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_size(range=c(0,8), limits = c(0,50), breaks=c(5,25,50))+
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
  theme(axis.text.y = element_text(color="black", size=8))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y")+
  scale_fill_manual(values=c("#FDAE61", 'indianred3',"#FFFFBF","#ABDDA4",'hotpink2','grey'), name='Lineages', guide='none')+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())

genomes_Central

Central_fig<-plot_grid(epi_curve_Central+theme(axis.text.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.title.x = element_blank()),
                       genomes_Central+theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5)),ncol=1, rel_heights = c(0.35,0.65))
Central_fig


###East figure

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

genomes_East<-ggplot(data=subset(df_africa_second, RCC=='East'))  + theme_minimal_hgrid()+
  geom_count(data=subset(subset(df_africa_second, RCC=='East'), Variant %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='East'), Variant %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='East'), Variant %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='East'), Variant %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='East'), Pango.lineage=='A.23.1'),aes(fill='A.23.1',x=date3, y='A.23.1'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='East'), !Variant %like% 'Beta' & !Variant %like% 'Alpha' & !Variant %like% 'Delta' & !Variant %like% 'Omicron' & !Pango.lineage %like% 'A.23.1'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_size(range=c(0,8), limits = c(0,50), breaks=c(5,25,50))+
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
  theme(axis.text.y = element_text(color="black", size=8))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y")+
  scale_fill_manual(values=c('indianred2',"#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name='Lineages', guide='none')+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())

genomes_East

East_fig<-plot_grid(epi_curve_East+theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.title.x = element_blank()),
                    genomes_East+theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5)),ncol=1, rel_heights = c(0.35,0.65))
East_fig



###Islands figure

epi_curve_Islands<-ggplot(data=subset(owid_data, RCC=='Islands'), aes(x=date,y=new_cases_smoothed_per_million))+
  theme_minimal()+
  geom_ribbon(aes(ymin=0, ymax=new_cases_smoothed_per_million),fill='blue', alpha=0.3) +
  geom_line(color='blue', size=1)+
  xlab(" ")+
  ylab('Reported Daily Cases\nper Million')+
  facet_grid(.~location)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 month",limits = as.Date(c('2020/02/15', '2022/03/15')))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,1000))

epi_curve_Islands

genomes_Islands<-ggplot(data=subset(df_africa_second, RCC=='Islands'))  + theme_minimal_hgrid()+
  geom_count(data=subset(subset(df_africa_second, RCC=='Islands'), Variant %like% 'Beta'),aes(fill='Beta',x=date3, y='Beta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Islands'), Variant %like% 'Alpha'),aes(fill='Alpha',x=date3, y='Alpha'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Islands'), Variant %like% 'Delta'),aes(fill='Delta',x=date3, y='Delta'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Islands'), Variant %like% 'Omicron'),aes(fill='Omicron',x=date3, y='Omicron'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Islands'), Variant %like% 'Eta'),aes(fill='Eta',x=date3, y='Eta' ), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Islands'), Pango.lineage=='A.23.1'),aes(fill='A.23.1',x=date3, y='A.23.1'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  geom_count(data=subset(subset(df_africa_second, RCC=='Islands'), !Variant %like% 'Beta' & !Variant %like% 'Alpha' & !Variant %like% 'Delta' & !Variant %like% 'Omicron' & !Pango.lineage %like% 'A.23.1'),aes(fill='Others',x=date3, y='Others'), stroke=0.4,shape=21, col='black', alpha=0.7)+
  ylab('')+ xlab('month')+
  labs(size='Genomes\nSampled Monthly')+
  scale_size(range=c(0,8), limits = c(0,200), breaks=c(5,50,200))+
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
  theme(axis.text.y = element_text(color="black", size=8))+
  theme(axis.text.x = element_text(color="black", size=10, angle=90,hjust=1,vjust=0.5))+
  scale_x_date(date_labels = "%b-%Y")+
  scale_fill_manual(values=c('indianred2',"#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name='Lineages', guide='none')+
  xlab('Sampling Dates')+
  facet_grid(.~country)+
  theme(strip.background = element_blank(), strip.text = element_blank())

genomes_Islands

Islands_fig<-plot_grid(epi_curve_Islands+theme(axis.text.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.title.x = element_blank()),
                       genomes_Islands+theme(axis.text.x = element_text(color="black", size=8, angle=90,hjust=1,vjust=0.5)),ncol=1, rel_heights = c(0.35,0.65))
Islands_fig



#####Supplementary Figure S2

##Rename variants into factors and add levels

df_africa$Variants2 <- with(df_africa, factor(Variant, 
                                      levels = c('VOC Beta GH/501Y.V2 (B.1.351+B.1.351.2+B.1.351.3) first detected in South Africa', 'VOC Delta GK (B.1.617.2+AY.*) first detected in India', 'VOC Omicron GRA (B.1.1.529+BA.*) first detected in Botswana/Hong Kong/South Africa','VOC Alpha GRY (B.1.1.7+Q.*) first detected in the UK'), 
                                      labels = c("Beta",  "Delta", "Omicron", "Alpha")))
df_africa$Variants2 <- as.character(df_africa$Variants2)
df_africa$Variants2[is.na(df_africa$Variants2)] = "Other"

df_africa$Variants2<-factor(df_africa$Variants2,levels = c("Other","Beta","Delta", "Omicron", "Alpha"))

##Add 20 day time lag

df_africa$days<-as.Date(cut(df_africa$Collection.date,breaks = "day",start.on.monday = FALSE))
df_africa2$days20 <- df_africa$days+20

##Epi data

owid<-read_excel('owid.xlsx')
africa_owid <- subset(owid, location == 'Africa')
africa_owid$days<-as.Date(cut(africa_owid$date,breaks = "day",start.on.monday = FALSE))


##Calculate proportion of variants from genomic data and attribute deaths to specific variants

P_africa20 <- prop.table(table(df_africa$days20, df_africa$Variants2), margin=1)

temp_africa20<-as.data.frame(P_africa20)
names(temp_africa20)[1] <- 'days'
head(temp_africa20)

temp_africa20$days<-as.Date(cut(as.Date(temp_africa20$days,format='%Y-%m-%d'),
                                breaks = "day",
                                start.on.monday = FALSE))

temp2_africa20<-africa_owid[c("days","new_deaths_per_million")]
head(temp2_africa20)

temp3_africa20<-join(temp_africa20, temp2_africa20,
                     type = "left")

tail(temp3_africa20)

temp3_africa20$days<-as.Date(cut(as.Date(temp3_africa20$days,format='%Y-%m-%d'),
                                 breaks = "day",
                                 start.on.monday = FALSE))

temp3_africa20$deaths_per_variant=temp3_africa20$new_deaths_per_million*temp3_africa20$Freq


##Smoothing

temp3_africa20 <- temp3_africa20 %>%
  group_by(Var2) %>%
  dplyr::mutate(deaths_per_variant_7day = zoo::rollmean(deaths_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_africa20)

##Plotting

dateVec <- seq(from = as.Date("2020-02-01"), to = as.Date("2022-03-01"), by = "days")

deaths_africa20<- ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "3 months")+
  scale_fill_manual(values=c("grey", "#FFFFBF","#ABDDA4","hotpink2", "#FDAE61"), labels=c("Other","Beta","Delta", "Omicron","Alpha"))+
  geom_density(data=temp3_africa20, aes(x = days, y = deaths_per_variant_7day, fill = Var2),stat="identity",size=0.3, alpha = 0.6)+
  geom_line(data=africa_owid, aes(x=days, y=proportion_people_fully_vaccinated/25), size=0.8, color = "black")+
  xlab('')+
  scale_y_continuous(
    
    # Features of the first axis
    name = "Reported daily deaths/million", limits=c(0, 0.6),breaks = c(0, 0.2, 0.4, 0.6),labels = c(0, 0.2, 0.4, 0.6),
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*25, name="Percentage of population fully vaccinated")
    
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_text(size=0))+
  theme(legend.position = "top", legend.text= element_text(size=8),legend.title=element_text(size=0), legend.spacing.x = unit(0.1, 'cm')) +
  theme(axis.title.y = element_text(hjust = 0.5, size=10, face="bold"))+
  theme(legend.position = c(0.1, 0.8))
  
deaths_africa20


#####Supplementary Figure S3 B

library(MASS)
library(pscl)

##Getting GISAID data into genome counts per time interval

genomes_week <- data.frame(table(cut(africa_100k$Collection.date, breaks ="1 week",start.on.monday = TRUE)))
genomes_week$Var1<-as.Date(genomes_week$Var1)
colnames(genomes_week) <- c('week','genomes')

##Testing normality = genome data not normally distributed
hist(genomes_week$genomes)
shapiro.test(genomes_week$genomes)

##case data

owid_africa<-subset(owid_covid_data, location=='Africa',select = c('location', 'date', 'new_cases_per_million'))
owid_africa$week<-as.Date(cut(owid_africa$date,breaks = "1 week",start.on.monday = TRUE))
owid_week <- owid_africa %>%
  aggregate(cbind(new_cases_per_million) ~ week, FUN=sum)


##Testing normality = data also skewed
hist(owid_week$new_cases_per_million)
shapiro.test(owid_week$new_cases_per_million)

##Merge datasets

week<-merge(genomes_week, owid_week,
                 type = "left")


##Calculate variance and mean 
mean(week$genomes)
var(week$genomes)
mean((week$genomes - mean(week$genomes)) ^ 2)

mean(week$new_cases_per_million)
var(week$new_cases_per_million)
mean((week$new_cases_per_million - mean(week$new_cases_per_million)) ^ 2)


##Neg bin models with MASS package
summary(week_model <- glm.nb(genomes ~ new_cases_per_million, data = week))

##Check if poisson model is more appropriate
pois <- glm(genomes ~ new_cases_per_million, family = "poisson",  data = week)
summary(pois)
pchisq(2 * (logLik(week_model2) - logLik(pois)), df = 1, lower.tail = FALSE)

###estimate confidence intervals from model
(est <- cbind(Estimate = coef(week_model), confint(week_model2)))


##log-likelihoods of a negative binomial and Poisson model
odTest(week_model)

##Plotting

scatter <- ggplot(data = week)+
  theme_classic()+
  geom_point(aes(x=new_cases_per_million, y=genomes),shape=20, size = 3, col='tan')+
  geom_smooth(aes(x=new_cases_per_million, y=genomes), method = 'glm.nb',color="black")+
  xlab("Reported weekly cases per million") + ylab("Weekly number of genomes")
scatter


#####Supplementary Figure S4

options(scipen=999)
owid$days<-as.Date(cut(owid$date,breaks = "day",start.on.monday = FALSE))

#Calculate tests per million
owid$new_tests_smoothed_per_million <- (owid$new_tests_smoothed / owid$population)*1000000
owid$new_cases_per_million[is.na(owid$new_cases_per_million)] = 0
owid$new_tests_smoothed_per_million[is.na(owid$new_tests_smoothed_per_million)] = 0

#Harmonise country names
owid$location[owid$location == "Cote d'Ivoire"] <- "Ivory Coast"
owid$location[owid$location == "Guinea/Bissau"] <- "Guinea-Bissau"

#Create country lists per region
country_South<-c('Botswana', 'South Africa', 'Zambia', 'Angola', 'Eswatini','Lesotho', 'Malawi', 'Mozambique', 'Namibia', 'Zimbabwe')
RCC_South<-c('South','South','South','South','South','South', 'South', 'South','South', 'South')

country_East<-c('Ethiopia', 'Kenya', 'Uganda', 'Comoros','Djibouti', 'Madagascar', 'Rwanda', 'Somalia', 'South Sudan', 'Sudan', 'Comoros', 'Mayotte', 'Reunion', 'Seychelles', 'Mauritius', 'Tanzania', 'Eritrea')
RCC_East<-c('East','East','East','East','East','East','East','East','East','East','East','East','East','East','East', 'East', 'East')

country_West<-c('Ghana', 'Nigeria', 'Senegal', 'Benin','Burkina Faso', "Ivory Coast", 'Gambia', 'Guinea', 'Guinea-Bissau', 'Liberia','Mali','Niger','Nigeria','Sierra Leone', 'Togo',  'Mauritania', 'Sao Tome and Principe', 'Cape Verde')
RCC_West<-c('West','West','West','West','West', 'West', 'West','West','West','West','West', 'West','West','West', 'West', 'West', 'West', 'West')

country_Central<-c('Cameroon','Congo','Democratic Republic of Congo', 'Burundi','Central African Republic', 'Chad', 'Equatorial Guinea','Gabon')
RCC_Central<-c('Central','Central','Central', 'Central','Central','Central', 'Central','Central')

country_North<-c('Egypt','Morocco', 'Tunisia', 'Algeria','Libya')
RCC_North<-c('North','North','North','North','North')

country<-c(country_South,country_East,country_West,country_Central, country_North)
RCC<-c(RCC_South,RCC_East,RCC_West,RCC_Central,RCC_North)

RCC_countries <- data.frame(country, RCC)

##Create dataframes for time intervals of mean tests and cases

df1 <- owid[owid$days >= "2020-01-01" & owid$date <= "2020-04-30", ]

df1a <- df1 %>%
  group_by(location) %>%
  summarise_at(vars(new_cases_per_million), list(cases = mean))

df1b <- df1 %>%
  group_by(location) %>%
  summarise_at(vars(new_tests_smoothed_per_million), list(tests = mean))

df1c <- df1a %>% 
  left_join(df1b, by = c("location" = "location"))%>%
  left_join(RCC_countries, by = c("location" = "country"))
df1c <- df1c[df1c$tests != 0, ]


df2 <- owid[owid$days >= "2020-05-01" & owid$date <= "2020-08-31", ]
df2a <- df2 %>%
  group_by(location) %>%
  summarise_at(vars(new_cases_per_million), list(cases = mean))

df2b <- df2 %>%
  group_by(location) %>%
  summarise_at(vars(new_tests_smoothed_per_million), list(tests = mean))

df2c <- df2a %>% 
  left_join(df2b, by = c("location" = "location"))%>%
  left_join(RCC_countries, by = c("location" = "country"))
df2c <- df2c[df2c$tests != 0, ]

df3 <- owid[owid$days >= "2020-09-01" & owid$date <= "2020-12-31", ]
df3a <- df3 %>%
  group_by(location) %>%
  summarise_at(vars(new_cases_per_million), list(cases = mean))

df3b <- df3 %>%
  group_by(location) %>%
  summarise_at(vars(new_tests_smoothed_per_million), list(tests = mean))

df3c <- df3a %>% 
  left_join(df3b, by = c("location" = "location"))%>%
  left_join(RCC_countries, by = c("location" = "country"))
df3c <- df3c[df3c$tests != 0, ]

df4 <- owid[owid$days >= "2021-01-01" & owid$date <= "2021-04-30", ]
df4a <- df4 %>%
  group_by(location) %>%
  summarise_at(vars(new_cases_per_million), list(cases = mean))

df4b <- df4 %>%
  group_by(location) %>%
  summarise_at(vars(new_tests_smoothed_per_million), list(tests = mean))

df4c <- df4a %>% 
  left_join(df4b, by = c("location" = "location"))%>%
  left_join(RCC_countries, by = c("location" = "country"))
df4c <- df4c[df4c$tests != 0, ]

df5 <- owid[owid$days >= "2021-05-01" & owid$date <= "2021-08-31", ]
df5a <- df5 %>%
  group_by(location) %>%
  summarise_at(vars(new_cases_per_million), list(cases = mean))

df5b <- df5 %>%
  group_by(location) %>%
  summarise_at(vars(new_tests_smoothed_per_million), list(tests = mean))

df5c <- df5a %>% 
  left_join(df5b, by = c("location" = "location"))%>%
  left_join(RCC_countries, by = c("location" = "country"))
df5c <- df5c[df5c$tests != 0, ]

df6 <- owid[owid$days >= "2021-09-01" & owid$date <= "2021-12-31", ]
df6a <- df6 %>%
  group_by(location) %>%
  summarise_at(vars(new_cases_per_million), list(cases = mean))

df6b <- df6 %>%
  group_by(location) %>%
  summarise_at(vars(new_tests_smoothed_per_million), list(tests = mean))

df6c <- df6a %>% 
  left_join(df6b, by = c("location" = "location"))%>%
  left_join(RCC_countries, by = c("location" = "country"))
df6c <- df6c[df6c$tests != 0, ]

df7 <- owid[owid$days >= "2022-01-01" & owid$date <= "2022-03-01", ]

df7a <- df7 %>%
  group_by(location) %>%
  summarise_at(vars(new_cases_per_million), list(cases = mean))

df7b <- df7 %>%
  group_by(location) %>%
  summarise_at(vars(new_tests_smoothed_per_million), list(tests = mean))

df7c <- df7a %>% 
  left_join(df7b, by = c("location" = "location"))%>%
  left_join(RCC_countries, by = c("location" = "country"))
df7c <- df7c[df7c$tests != 0, ]

##Plotting

Fig1 <- ggplot()+
  theme_classic()+
  geom_segment(data=Fig1_coords,aes(x = x1, y = y1, xend = x2, yend = y2),colour="grey",alpha=0.5, linetype = 5)+
  geom_point(data = df1c, aes(x=cases, y = tests, fill = RCC), shape=21, size=4, stroke=0.2)+
  scale_fill_manual(values=c("#006d77","#83c5be","#edf6f9",'#ffddd2','#e29578'), name=' ')+
  ylab('Average daily tests per million')+ xlab('Average daily cases per million')+ 
  geom_text_repel(data = df1c, aes(x=cases, y = tests,label = location), size =2)+
  ggtitle('Jan - Apr 2020')+
  geom_text(data = Fig1_coords, aes(x = x2, y = y2,label = percent), vjust=-0.5, size =2, colour="grey")+
  scale_y_log10() + scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
  theme(axis.title = element_text(size = 8))+
  theme(axis.text = element_text(size = 8))

Fig1

Fig2<- ggplot()+
  theme_classic()+
  geom_segment(data=Fig2_coords,aes(x = x1, y = y1, xend = x2, yend = y2),colour="grey",alpha=0.5, linetype = 5)+
  geom_point(data = df2c, aes(x=cases, y = tests, fill = RCC), shape=21, size=4, stroke=0.2)+
  scale_fill_manual(values=c("#006d77","#83c5be","#edf6f9",'#ffddd2','#e29578'), name=' ')+
  ylab('Average daily tests per million')+ xlab('Average daily cases per million')+ 
  geom_text_repel(data = df2c, aes(x=cases, y = tests,label = location), size =2)+
  ggtitle('May - Aug 2020')+
  geom_text(data = Fig2_coords, aes(x = x2, y = y2,label = percent), vjust=-0.5, size =2, colour="grey")+
  scale_y_log10() + scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
  theme(axis.title = element_text(size = 8))+
  theme(axis.text = element_text(size = 8))

Fig2

Fig3<- ggplot()+
  theme_classic()+
  geom_segment(data=Fig3_coords,aes(x = x1, y = y1, xend = x2, yend = y2),colour="grey",alpha=0.5, linetype = 5)+
  geom_point(data = df3c, aes(x=cases, y = tests, fill = RCC), shape=21, size=4, stroke=0.2)+
  scale_fill_manual(values=c("#006d77","#83c5be","#edf6f9",'#ffddd2','#e29578'), name=' ')+
  ylab('Average daily tests per million')+ xlab('Average daily cases per million')+ 
  geom_text_repel(data = df3c, aes(x=cases, y = tests,label = location), size =2)+
  ggtitle('Sep - Dec 2020')+
  geom_text(data = Fig3_coords, aes(x = x2, y = y2,label = percent), vjust=-0.5, size =2, colour="grey")+
  scale_y_log10() + scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
  theme(axis.title = element_text(size = 8))+
  theme(axis.text = element_text(size = 8))

Fig3

Fig4<- ggplot()+
  theme_classic()+
  geom_segment(data=Fig4_coords,aes(x = x1, y = y1, xend = x2, yend = y2),colour="grey",alpha=0.5, linetype = 5)+
  geom_point(data = df4c, aes(x=cases, y = tests, fill = RCC), shape=21, size=4, stroke=0.2)+
  scale_fill_manual(values=c("#006d77","#83c5be","#edf6f9",'#ffddd2','#e29578'), name=' ')+
  ylab('Average daily tests per million')+ xlab('Average daily cases per million')+ 
  geom_text_repel(data = df4c, aes(x=cases, y = tests,label = location), size =2)+
  ggtitle('Jan - Apr 2021')+
  geom_text(data = Fig4_coords, aes(x = x2, y = y2,label = percent), vjust=-0.5, size =2, colour="grey")+
  scale_y_log10() + scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
  theme(axis.title = element_text(size = 8))+
  theme(axis.text = element_text(size = 8))

Fig4

Fig5<-  ggplot()+
  theme_classic()+
  geom_segment(data=Fig5_coords,aes(x = x1, y = y1, xend = x2, yend = y2),colour="grey",alpha=0.5, linetype = 5)+
  geom_point(data = df5c, aes(x=cases, y = tests, fill = RCC), shape=21, size=4, stroke=0.2)+
  scale_fill_manual(values=c("#006d77","#83c5be","#edf6f9",'#ffddd2','#e29578'), name=' ')+
  ylab('Average daily tests per million')+ xlab('Average daily cases per million')+ 
  geom_text_repel(data = df5c, aes(x=cases, y = tests,label = location), size =2)+
  ggtitle('May - Aug 2021')+
  geom_text(data = Fig5_coords, aes(x = x2, y = y2,label = percent), vjust=-0.5, size =2, colour="grey")+
  scale_y_log10() + scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
  theme(axis.title = element_text(size = 8))+
  theme(axis.text = element_text(size = 8))

Fig5

Fig6<- ggplot()+
  theme_classic()+
  geom_segment(data=Fig6_coords,aes(x = x1, y = y1, xend = x2, yend = y2),colour="grey",alpha=0.5, linetype = 5)+
  geom_point(data = df6c, aes(x=cases, y = tests, fill = RCC), shape=21, size=4, stroke=0.2)+
  scale_fill_manual(values=c("#006d77","#83c5be","#edf6f9",'#ffddd2','#e29578'), name=' ')+
  ylab('Average daily tests per million')+ xlab('Average daily cases per million')+ 
  geom_text_repel(data = df6c, aes(x=cases, y = tests,label = location), size =2)+
  ggtitle('Sep - Dec 2021')+
  geom_text(data = Fig6_coords, aes(x = x2, y = y2,label = percent), vjust=-0.5, size =2, colour="grey")+
  scale_y_log10() + scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
  theme(axis.title = element_text(size = 8))+
  theme(axis.text = element_text(size = 8))

Fig6

Fig7<- ggplot()+
  theme_classic()+
  geom_segment(data=Fig7_coords4,aes(x = x1, y = y1, xend = x2, yend = y2),colour="grey",alpha=0.5, linetype = 5)+
  geom_point(data = df7c, aes(x=cases, y = tests, fill = RCC), shape=21, size=4, stroke=0.2)+
  scale_fill_manual(values=c("#006d77","#83c5be","#edf6f9",'#ffddd2','#e29578'), name=' ')+
  ylab('Average daily tests per million')+ xlab('Average daily cases per million')+ 
  geom_text_repel(data = df7c, aes(x=cases, y = tests,label = location), size =2)+
  ggtitle('Jan - Mar 2022')+
  geom_text(data = Fig7_coords4, aes(x = x2, y = y2,label = percent), vjust=-0.5, size =2, colour="grey")+
  scale_y_log10() + scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
  theme(axis.title = element_text(size = 8))+
  theme(axis.text = element_text(size = 8))

Fig7

comp <- ggarrange(Fig1, Fig2, Fig3, Fig4, Fig5, Fig6, Fig7,ncol =2,nrow=4, common.legend = TRUE)
comp


######Supplementary Figure S5

#Read in cleaned data
sequencing_tech_all_cleaned <- read_excel("sequencing_tech_all_cleaned.xlsx")
sequencing_tech_all_cleaned$date2<-as.Date(cut(sequencing_tech_all_cleaned$date,breaks = "2 weeks",start.on.monday = FALSE))

#Plotting

p_seq_tech<-ggplot(data=sequencing_tech_all_cleaned, mapping = aes(x = date2, fill=seq_tech))+
theme_classic()+
geom_bar(width=10,color='black', size=0.2)+
theme(axis.title.x = element_text(color="black", size=11, face="bold"))+
theme(axis.title.y = element_text(color="black", size=11, face="bold"))+
scale_x_date(date_labels = "%b\n%Y", date_breaks = "2 month")+
scale_fill_brewer(palette = 'Spectral', name='Sequencing\nTechnology')+
ylab('Genome Count')+
xlab('Sampling Date')

p_seq_tech


#####Supplementary Figure S7

#read in golbal GISAID metadata
global<-read.csv('metadata.tsv',sep = "\t")
africa<-subset(global, region=='Africa')
africa$date<-as.Date(africa$date)

africa$days<-as.Date(cut(africa$date,breaks = "day",start.on.monday = FALSE))
africa$date2<-as.Date(cut(africa$date,breaks = "2 weeks",start.on.monday = FALSE))
africa$date_submitted<-as.Date(africa$date_submitted)
africa$days_submitted<-as.Date(cut(africa$date_submitted,breaks = "day",start.on.monday = FALSE))
africa$submission_lag=africa$days_submitted-africa$days


dateVec <- seq(from = as.Date("2020-02-01"), to = as.Date("2022-02-01"), by = "days")

Burundi<-subset(africa, country=='Burundi')

Burundilag<-ggplot(data=Burundi, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Burundi')
  

Burundilag

Cameroon<-subset(africa, country=='Cameroon')

Cameroonlag<-ggplot(data=Cameroon, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Cameroon')+
  scale_y_continuous(breaks=seq(0,500,50))

Cameroonlag

CAR<-subset(africa, country=='Central African Republic')

CARlag<-ggplot(data=CAR, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Central African Republic')+
  scale_y_continuous(breaks=seq(0,500,50))

CARlag

Chad<-subset(africa, country=='Chad')

Chadlag<-ggplot(data=Chad, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Chad')+
  scale_y_continuous(breaks=seq(0,500,50))

Chadlag

RC<-subset(africa, country=='Republic of the Congo')

RClag<-ggplot(data=RC, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Republic of the Congo')+
  scale_y_continuous(breaks=seq(0,500,50))

RClag

DRC<-subset(africa, country=='Democratic Republic of the Congo')

DRClag<-ggplot(data=DRC, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Democratic Republic of the Congo')+
  scale_y_continuous(breaks=seq(0,500,50))

DRClag

EG<-subset(africa, country=='Equatorial Guinea')

EGlag<-ggplot(data=EG, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Equatorial Guinea')+
  scale_y_continuous(breaks=seq(0,500,50))

EGlag

Gabon<-subset(africa, country=='Gabon')

Gabonlag<-ggplot(data=Gabon, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Gabon')+
  scale_y_continuous(breaks=seq(0,500,50))

Gabonlag

STP<-subset(africa, country=='Sao Tome and Principe')

STPlag<-ggplot(data=STP, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('São Tomé and Príncipe')+
  scale_y_continuous(breaks=seq(0,500,50))

STPlag


central <- ggarrange(Burundilag, Cameroonlag, CARlag, Chadlag, RClag, DRClag, EGlag, Gabonlag, STPlag, nrow=3, ncol=3, align = "v")
central
annotate_figure(central, left = text_grob('Days from specimen collection to sequence submission', color = "black", hjust = 0.5, face = "bold", rot = 90), top = text_grob('Central Africa', color = "black", hjust = 0.5, face = "bold"))


Comoros<-subset(africa, country=='Union of the Comoros')

Comoroslag<-ggplot(data=Comoros, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Comoros')
 
Comoroslag


Djibouti<-subset(africa, country=='Djibouti')

Djiboutilag<-ggplot(data=Djibouti, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Djibouti')+
  scale_y_continuous(breaks=seq(0,500,50))

Djiboutilag

Ethiopia<-subset(africa, country=='Ethiopia')

Ethiopialag<-ggplot(data=Ethiopia, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Ethiopia')+
  scale_y_continuous(breaks=seq(0,500,50))

Ethiopialag

Kenya<-subset(africa, country=='Kenya')

Kenyalag<-ggplot(data=Kenya, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Kenya')+
  scale_y_continuous(breaks=seq(0,500,50))

Kenyalag

Madagascar<-subset(africa, country=='Madagascar')

Madagascarlag<-ggplot(data=Madagascar, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Madagascar')+
  scale_y_continuous(breaks=seq(0,500,50))

Madagascarlag


Mauritius<-subset(africa, country=='Mauritius')

Mauritiuslag<-ggplot(data=Mauritius, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Mauritius')+
  scale_y_continuous(breaks=seq(0,500,50))

Mauritiuslag

Rwanda<-subset(africa, country=='Rwanda')

Rwandalag<-ggplot(data=Rwanda, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Rwanda')+
  scale_y_continuous(breaks=seq(0,500,50))

Rwandalag

Seychelles<-subset(africa, country=='Seychelles')

Seychelleslag<-ggplot(data=Seychelles, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Seychelles')+
  scale_y_continuous(breaks=seq(0,500,50))

Seychelleslag


Somalia<-subset(africa, country=='Somalia')

Somalialag<-ggplot(data=Somalia, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Somalia')+
  scale_y_continuous(breaks=seq(0,500,50))

Somalialag

SS<-subset(africa, country=='South Sudan')

SSlag<-ggplot(data=SS, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('South Sudan')+
  scale_y_continuous(breaks=seq(0,500,50))

SSlag

S<-subset(africa, country=='Sudan')

Slag<-ggplot(data=S, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Sudan')+
  scale_y_continuous(breaks=seq(0,500,50))

Slag

Tanzania<-subset(africa, country=='Tanzania')

Tanzanialag<-ggplot(data=Tanzania, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Tanzania')+
  scale_y_continuous(limits = c(min(0), max=max(50)))

Tanzanialag

Uganda<-subset(africa, country=='Uganda')

Ugandalag<-ggplot(data=Uganda, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Uganda')+
  scale_y_continuous(breaks=seq(0,500,50))

Ugandalag

eastern1 <- ggarrange(Comoroslag, Djiboutilag, Ethiopialag, Kenyalag, Madagascarlag, Mauritiuslag, Rwandalag, Seychelleslag, Somalialag,SSlag, Slag, Ugandalag, nrow=4, ncol=3, align = "v")
eastern1
annotate_figure(eastern1, left = text_grob('Days from specimen collection to sequence submission', color = "black", hjust = 0.5, face = "bold", rot = 90), top = text_grob('Eastern Africa', color = "black", hjust = 0.5, face = "bold"))

eastern2 <- ggarrange(SSlag, Slag, Ugandalag, nrow=2, ncol=3, align = "v")
eastern2
annotate_figure(eastern2, left = text_grob('Days from specimen collection to sequence submission', color = "black", hjust = 0.5, face = "bold", rot = 90), top = text_grob('Eastern Africa', color = "black", hjust = 0.5, face = "bold"))

Algeria<-subset(africa, country=='Algeria')

Algerialag<-ggplot(data=Algeria, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Algeria')+
  scale_y_continuous(breaks=seq(0,500,50))

Algerialag


Egypt<-subset(africa, country=='Egypt')

Egyptlag<-ggplot(data=Egypt, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Egypt')+
  scale_y_continuous(breaks=seq(0,500,50))

Egyptlag

Libya<-subset(africa, country=='Libya')

Libyalag<-ggplot(data=Libya, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Libya')+
  scale_y_continuous(breaks=seq(0,500,50))

Libyalag

Morocco<-subset(africa, country=='Morocco')

Moroccolag<-ggplot(data=Morocco, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Morocco')+
  scale_y_continuous(breaks=seq(0,500,50))

Moroccolag

Tunisia<-subset(africa, country=='Tunisia')

Tunisialag<-ggplot(data=Tunisia, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Tunisia')+
  scale_y_continuous(breaks=seq(0,500,50))

Tunisialag


northern <- ggarrange(Algerialag, Egyptlag, Libyalag, Moroccolag, Tunisialag, nrow=2, ncol=3, align = "v")
northern
annotate_figure(northern, left = text_grob('Days from specimen collection to sequence submission', color = "black", hjust = 0.5, face = "bold", rot = 90), top = text_grob('Northern Africa', color = "black", hjust = 0.5, face = "bold"))

Angola<-subset(africa, country=='Angola')

Angolalag<-ggplot(data=Angola, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Angola')+
  scale_y_continuous(breaks=seq(0,500,50))

Angolalag


Botswana<-subset(africa, country=='Botswana')

Botswanalag<-ggplot(data=Botswana, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Botswana')+
  scale_y_continuous(breaks=seq(0,500,50))

Botswanalag


Eswatini<-subset(africa, country=='Eswatini')

Eswatinilag<-ggplot(data=Eswatini, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Eswatini')+
  scale_y_continuous(breaks=seq(0,500,50))

Eswatinilag


Lesotho<-subset(africa, country=='Lesotho')

Lesotholag<-ggplot(data=Lesotho, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Lesotho')+
  scale_y_continuous(breaks=seq(0,150,50))

Lesotholag


Malawi<-subset(africa, country=='Malawi')

Malawilag<-ggplot(data=Malawi, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Malawi')+
  scale_y_continuous(breaks=seq(0,500,50))

Malawilag


Mozambique<-subset(africa, country=='Mozambique')

Mozambiquelag<-ggplot(data=Mozambique, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Mozambique')+
  scale_y_continuous(breaks=seq(0,500,50))

Mozambiquelag

Namibia<-subset(africa, country=='Namibia')

Namibialag<-ggplot(data=Namibia, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Namibia')+
  scale_y_continuous(breaks=seq(0,500,50))

Namibialag

SA<-subset(africa, country=='South Africa')

SAlag<-ggplot(data=SA, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('South Africa')+
  scale_y_continuous(breaks=seq(0,500,50))

SAlag

Zambia<-subset(africa, country=='Zambia')

Zambialag<-ggplot(data=Zambia, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Zambia')+
  scale_y_continuous(breaks=seq(0,500,50))

Zambialag

Zimbabwe<-subset(africa, country=='Zimbabwe')

Zimbabwelag<-ggplot(data=Zimbabwe, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Zimbabwe')+
  scale_y_continuous(breaks=seq(0,500,50))

Zimbabwelag


southern <- ggarrange(Angolalag, Botswanalag, Eswatinilag, Lesotholag, Malawilag, Mozambiquelag, Namibialag, SAlag, Zambialag, Zimbabwelag, nrow=4, ncol=3, align = "v")
southern
annotate_figure(southern, left = text_grob('Days from specimen collection to sequence submission', color = "black", hjust = 0.5, face = "bold", rot = 90), top = text_grob('Southern Africa', color = "black", hjust = 0.5, face = "bold"))


Benin<-subset(africa, country=='Benin')


Beninlag<-ggplot(data=Benin, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Benin')+
  scale_y_continuous(breaks=seq(0,500,50))

Beninlag


Cabo<-subset(africa, country=='Cape Verde')

Cabolag<-ggplot(data=Cabo, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Cape Verde')+
  scale_y_continuous(breaks=seq(0,500,50))

Cabolag


BurkinaFaso<-subset(africa, country=='Burkina Faso')

BurkinaFasolag<-ggplot(data=BurkinaFaso, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Burkina Faso')+
  scale_y_continuous(breaks=seq(0,500,50))

BurkinaFasolag

Coted_Ivoire<-subset(africa, country=='Côte d’Ivoire')

Coted_Ivoirelag<-ggplot(data=Coted_Ivoire, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_text(color="black", size=12, face = "bold"))+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Côte d’Ivoire')+
  scale_y_continuous(breaks=seq(0,500,50))

Coted_Ivoirelag

Gambia<-subset(africa, country=='Gambia')

Gambialag<-ggplot(data=Gambia, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Gambia')+
  scale_y_continuous(breaks=seq(0,500,50))

Gambialag


Ghana<-subset(africa, country=='Ghana')

Ghanalag<-ggplot(data=Ghana, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Ghana')+
  scale_y_continuous(breaks=seq(0,500,50))

Ghanalag

Guinea<-subset(africa, country=='Guinea')

Guinealag<-ggplot(data=Guinea, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Guinea')+
  scale_y_continuous(breaks=seq(0,500,50))

Guinealag

GuineaBissau<-subset(africa, country=='Guinea-Bissau')

GuineaBissaulag<-ggplot(data=GuineaBissau, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Guinea-Bissau')+
  scale_y_continuous(breaks=seq(0,500,50))

GuineaBissaulag


Liberia<-subset(africa, country=='Liberia')

Liberialag<-ggplot(data=Liberia, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Liberia')+
  scale_y_continuous(breaks=seq(0,500,50))

Liberialag

Mali<-subset(africa, country=='Mali')

Malilag<-ggplot(data=Mali, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Mali')+
  scale_y_continuous(breaks=seq(0,500,50))

Malilag

Niger<-subset(africa, country=='Niger')

Nigerlag<-ggplot(data=Niger, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Niger')+
  scale_y_continuous(breaks=seq(0,500,50))

Nigerlag

Nigeria<-subset(africa, country=='Nigeria')

Nigerialag<-ggplot(data=Nigeria, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%Y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Nigeria')+
  scale_y_continuous(breaks=seq(0,500,50))

Nigerialag

Senegal<-subset(africa, country=='Senegal')

Senegallag<-ggplot(data=Senegal, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Senegal')+
  scale_y_continuous(breaks=seq(0,500,50))

Senegallag


Sierra<-subset(africa, country=='Sierra Leone')

Sierralag<-ggplot(data=Sierra, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Sierra Leone')+
  scale_y_continuous(breaks=seq(0,500,50))

Sierralag

Togo<-subset(africa, country=='Togo')

Togolag<-ggplot(data=Togo, aes(x=date2, y=as.numeric(submission_lag)))+
  theme_minimal()+
  geom_smooth(color='dodgerblue3',size=1)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "grey30" )+
  ylab('Days from specimen collection to sequence submission')+
  xlab('')+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2.5 month")+
  theme(axis.text.x = element_text(color="black", size=9))+
  theme(axis.text.y = element_text(color="black", size=9))+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))+
  ggtitle('Togo')+
  scale_y_continuous(breaks=seq(0,500,50))

Togolag

west1<-ggarrange(Beninlag, BurkinaFasolag, Cabolag, Coted_Ivoirelag, Gambialag, Ghanalag, Guinealag, GuineaBissaulag, Liberialag, nrow=3, ncol=3, align = "v")
annotate_figure(west1, top = text_grob('Western Africa', color = "black", hjust = 0.5, face = "bold"))

west2 <- ggarrange(Malilag, Nigerlag, Nigerialag, Senegallag, Sierralag, Togolag, nrow=2, ncol=3, align = "v")
west2
annotate_figure(west2, left = text_grob('Days from specimen collection to sequence submission', color = "black", hjust = 0.5, face = "bold", rot = 90), top = text_grob('Western Africa', color = "black", hjust = 0.5, face = "bold"))









######Supplementary Figure S12

#Read in genomic data and change country names
africa_100k <- read_csv("africa_100k.csv")
df <- africa_100k
df$days<-as.Date(cut(df$Collection.date,breaks = "day",start.on.monday = FALSE))
df$country[df$country == "Cabo Verde"] <- "Cape Verde"
df$country[df$country == "Cote d'Ivoire"] <- "Ivory Coast"
df$country[df$country == "Union of the Comoros"] <- "Comoros"
df$country[df$country == "Republic of the Congo"]  <- "Republic of Congo"
df$country[df$country == "Republic of the Congo"]  <- "Republic of Congo"

##Rename variants into factors and add levels

df$Variants2 <- with(df, factor(Variant, 
                                      levels = c('VOC Alpha GRY (B.1.1.7+Q.*) first detected in the UK', 'VOC Beta GH/501Y.V2 (B.1.351+B.1.351.2+B.1.351.3) first detected in South Africa', 'VOC Gamma GR/501Y.V3 (P.1+P.1.*) first detected in Brazil/Japan','VOC Delta GK (B.1.617.2+AY.*) first detected in India', 'VOC Omicron GRA (B.1.1.529+BA.*) first detected in Botswana/Hong Kong/South Africa', 'NA'), 
                                      labels = c("Alpha","Beta", "Gamma", "Delta", "Omicron", "Other")))

df$Variants2[is.na(df$Variants2)] = "Other"

##Create country lists to append to genomic data

country_South<-c('Botswana', 'South Africa', 'Zambia', 'Angola', 'Eswatini','Lesotho', 'Malawi', 'Mozambique', 'Namibia', 'Zimbabwe')
RCC_South<-c('South','South','South','South','South','South', 'South', 'South','South', 'South')

country_East<-c('Ethiopia', 'Kenya', 'Uganda', 'Comoros','Djibouti', 'Madagascar', 'Rwanda', 'Somalia', 'South Sudan', 'Sudan', 'Comoros', 'Mayotte', 'Reunion', 'Seychelles', 'Mauritius')
RCC_East<-c('East','East','East','East','East','East','East','East','East','East','East','East','East','East','East')

country_West<-c('Ghana', 'Nigeria', 'Senegal', 'Benin','Burkina Faso', "Ivory Coast", 'Gambia', 'Guinea', 'Guinea-Bissau', 'Liberia','Mali','Niger','Nigeria','Sierra Leone', 'Togo',  'Mauritania', 'Sao Tome and Principe', 'Cape Verde')
RCC_West<-c('West','West','West','West','West', 'West', 'West','West','West','West','West', 'West','West','West', 'West', 'West', 'West', 'West')

country_Central<-c('Cameroon','Republic of Congo','Democratic Republic of Congo', 'Burundi','Central African Republic', 'Chad', 'Equatorial Guinea','Gabon')
RCC_Central<-c('Central','Central','Central', 'Central','Central','Central', 'Central','Central')

country_North<-c('Egypt','Morocco', 'Tunisia', 'Algeria','Libya')
RCC_North<-c('North','North','North','North','North')

country<-c(country_South,country_East,country_West,country_Central, country_North)
RCC<-c(RCC_South,RCC_East,RCC_West,RCC_Central,RCC_North)

RCC_countries <- data.frame(country, RCC)

df2 <- df %>% 
  left_join(RCC_countries, by = c("country" = "country"))

North <- subset(df2,RCC=='North')
West <- subset(df2,RCC=='West')
Central <- subset(df2,RCC=='Central')
East <- subset(df2,RCC=='East')
South <- subset(df2,RCC=='South')

##Plotting figure B

dateVec <- seq(from = as.Date("2020-01-01"), to = as.Date("2022-03-01"), by = "days")

North_fig <-  ggplot(data=North) +
  theme_classic()+
  geom_stripes(aes(y = country, odd = "grey90", even = "#00000000")) +
  geom_segment(aes(x=min(days), y = country, xend=max(days), yend=country, group=country), colour="grey88", size=3) +
  geom_point(data=subset(North,Variants2=='Other'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(North,Variants2=='Alpha'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(North,Variants2=='Beta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(North,Variants2=='Delta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(North,Variants2=='Omicron'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  ylab('')+ xlab('')+ 
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2 month")+
  scale_color_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 8))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  ggtitle('Northern Africa')+
  theme(plot.title = element_text(size = 9))

North_fig


West_fig <-  ggplot(data=West) +
  theme_classic()+
  geom_stripes(aes(y = country, odd = "grey90", even = "#00000000")) +
  geom_segment(aes(x=min(days), y = country, xend=max(days), yend=country, group=country), colour="grey88", size=3) +
  geom_point(data=subset(West,Variants2=='Other'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(West,Variants2=='Alpha'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(West,Variants2=='Beta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(West,Variants2=='Delta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(West,Variants2=='Omicron'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  ylab('')+ xlab('')+ 
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2 month")+
  scale_color_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  theme(axis.text.x = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  ggtitle('Western Africa')+
  theme(axis.text.y = element_text(size = 8))+
  theme(plot.title = element_text(size = 9))

West_fig


Central_fig <-  ggplot(data=Central) +
  theme_classic()+
  geom_stripes(aes(y = country, odd = "grey90", even = "#00000000")) +
  geom_segment(aes(x=min(days), y = country, xend=max(days), yend=country, group=country), colour="grey88", size=3) +
  geom_point(data=subset(Central,Variants2=='Other'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(Central,Variants2=='Alpha'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(Central,Variants2=='Beta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(Central,Variants2=='Delta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(Central,Variants2=='Omicron'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  ylab('')+ xlab('')+ 
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2 month")+
  scale_color_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  theme(axis.text.x = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  ggtitle('Central Africa')+
  theme(axis.text.y = element_text(size = 8))+
  theme(plot.title = element_text(size = 9))

Central_fig


East_fig <-  ggplot(data=East) +
  theme_classic()+
  geom_stripes(aes(y = country, odd = "grey90", even = "#00000000")) +
  geom_segment(aes(x=min(days), y = country, xend=max(days), yend=country, group=country), colour="grey88", size=3) +
  geom_point(data=subset(East,Variants2=='Other'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(East,Variants2=='Alpha'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(East,Variants2=='Beta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(East,Variants2=='Delta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(East,Variants2=='Omicron'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  ylab('')+ xlab('')+ 
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2 month")+
  scale_color_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  theme(axis.text.x = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  ggtitle('Eastern Africa')+
  theme(axis.text.y = element_text(size = 8))+
  theme(plot.title = element_text(size = 9))

East_fig


South_fig <-  ggplot(data=South) +
  theme_classic()+
  geom_stripes(aes(y = country, odd = "grey90", even = "#00000000")) +
  geom_segment(aes(x=min(days), y = country, xend=max(days), yend=country, group=country), colour="grey88", size=3) +
  geom_point(data=subset(South,Variants2=='Other'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(South,Variants2=='Alpha'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(South,Variants2=='Beta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(South,Variants2=='Delta'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  geom_point(data=subset(South,Variants2=='Omicron'),aes(x=days, y=country, fill=Variants2,color=Variants2), position = position_jitter(width=0.3, height=0.3), shape=21, size=0.8, stroke=0.2, alpha = 0.5)+
  ylab('')+ xlab('')+ 
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "2 month")+
  scale_color_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  scale_fill_manual(values=c("#FDAE61","#FFFFBF","#ABDDA4",'hotpink2','grey'), name=' ')+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  ggtitle('Southern Africa')+
  theme(axis.text.y = element_text(size = 8))+
  theme(plot.title = element_text(size = 9))

South_fig

Comp <- ggarrange(North_fig,ggplot() + theme_void(), West_fig, ggplot() + theme_void(),Central_fig, ggplot() + theme_void(),East_fig, ggplot() + theme_void(),South_fig, heights = c(1,-0.2, 2,-0.2, 1.3,-0.2, 1.65,-0.2, 1.55), ncol =1, common.legend = TRUE, align = "v")
Comp

Comp2 <- annotate_figure(Comp, top = text_grob("B.", 
                                          color = "black", face = "bold", size = 12))
Comp2

##Retrieve map data

world3 <- map_data("world")

variantprops <-as.data.frame.matrix(table(df2 %>% select("country","Variants2")))

variantprops$country[variantprops$country == "Cabo Verde"] <- "Cape Verde"
variantprops$country[variantprops$country == "Cote d'Ivoire"] <- "Ivory Coast"
variantprops$country[variantprops$country == "Union of the Comoros"] <- "Comoros"
variantprops$country[variantprops$country == "Republic of the Congo"]  <- "Republic of Congo"
variantprops$country[variantprops$country == "Eswatini"]  <- "Swaziland"


variantprops <- cbind(rownames(variantprops), variantprops)
colnames(variantprops)[1] <- "country"


variantprops <- variantprops %>% 
  left_join(world3, by = c("country" = "region"))

##remove Marion Island

variantprops2 <- variantprops[-c(9966,9967,9968,9969,9970,9971,9972,9973,9974,9975), ]

variantprops2$total <- variantprops2$Alpha + variantprops2$Beta + variantprops2$Gamma + variantprops2$Delta + variantprops2$Omicron + variantprops2$Other

variantprops2$Alpha_prop <- variantprops2$Alpha/variantprops2$total
variantprops2$Beta_prop <- variantprops2$Beta/variantprops2$total
variantprops2$Gamma_prop <- variantprops2$Gamma/variantprops2$total
variantprops2$Delta_prop <- variantprops2$Delta/variantprops2$total
variantprops2$Omi_prop <- variantprops2$Omicron/variantprops2$total


##Plot A 

map_alpha <- ggplot() +
  theme_void()+
  geom_polygon(data = subset(variantprops2,!is.na(country)), aes(x = long, y = lat, group=country, fill=as.numeric(Alpha_prop))) +
  scale_fill_distiller(palette = "Oranges", direction = 1,na.value = "white",
  breaks = c(0.1, 0.2, 0.3, 0.4, 0.8, 1.0), labels = c(0.1, 0.2, 0.3, 0.4, 0.8, 1.0), name=' ')+
  coord_fixed()+
  ggtitle('Alpha')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face = "bold"))

map_alpha 

map_beta <- ggplot() +
  theme_void()+
  geom_polygon(data = subset(variantprops2,!is.na(country)), aes(x = long, y = lat, group=country, fill=as.numeric(Beta_prop))) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,na.value = "white",
                       breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0), labels = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0), name=' ')+
  coord_fixed()+
  ggtitle('Beta')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face = "bold"))

map_beta


map_delta <- ggplot() +
  theme_void()+
  geom_polygon(data = subset(variantprops2,!is.na(country)), aes(x = long, y = lat, group=country, fill=as.numeric(Delta_prop))) +
  scale_fill_distiller(palette = "Greens", direction = 1,na.value = "white", 
                       breaks = c(0.2,0.4, 0.6, 0.8, 1.0), labels = c(0.2,0.4, 0.6, 0.8, 1.0), name=' ')+
  coord_fixed()+
  ggtitle('Delta')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face = "bold"))

map_delta


map_omi <- ggplot() +
  theme_void()+
  geom_polygon(data = subset(variantprops2,!is.na(country)), aes(x = long, y = lat, group=country, fill=as.numeric(Omi_prop))) +
  scale_fill_distiller(palette = "RdPu", direction = 1,na.value = "white",
                       breaks = c(0.1, 0.2, 0.3, 0.4, 0.5), labels = c(0.1, 0.2, 0.3, 0.4, 0.5), name=' ')+
  coord_fixed()+
  ggtitle('Omicron')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face = "bold"))

map_omi

p1 <- ggarrange(map_alpha, map_beta, map_delta, map_omi, nrow=4, ncol=1)  
p1
p2 <- annotate_figure(p1, top = text_grob("A. Proportion of genomes per country", 
                                      color = "black", face = "bold", size = 11.5))
p2

p3 <- ggarrange(p2, Comp2, ncol=2, widths = c(1,2.5))
p3

p4 <- ggarrange(p1, Comp, ncol=2, widths = c(1,2.5))
p4


