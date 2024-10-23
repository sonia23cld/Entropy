#Entropy admixture plots for alyssum
k=5

admix.matrix<- read.table(paste0("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/simpleQ", k, ".txt"))
pop.code<- read.table("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/poplabel.txt")
ind.code<- read.table("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/indlabel.txt")
ind.pop<- merge(ind.code, pop.code)

data<- read.csv('/Users/sonia.celestini/Desktop/PhD Project/Alyssum_2024/Data/All_Samples.csv', header = T, sep=';')
data<- data[1:8,]
data$code2<- paste0(substr(data$code,1,2), "0", substr(data$code,3,4))

data$lat<- as.numeric(data$lat)
data$long<- as.numeric(data$long)

#k2
colors<- c("#1B32F5", "#EB3F24")

#k3
colors<- c("#EB3F24","#1B32F5", "#F6C44C")
#k4
colors<- c( "#1B32F5", "#458933", "#EB3F24", "#F6C44C")
#k5
colors<- c("#1B32F5", "#EB3F24",  "#B047F7", "#F6C44C", "#458933")
#k6
colors<- c("#EB3F24", "#B047F7", "#1B32F5", "#56BBF9", "#F6C44C", "#458933")

#colors<- c("#EB3F24", "#B047F7", "#1B32F5", "#56BBF9", "#F6C44C", "#458933", "#B047F7", , "#EF8532", "#7FEA9C", "#FEFB53", "black")

#start with map
# Set map boundary (xmin, xmax, ymin, ymax)
library(ggplot2)
library(raster)
library(maps)
library(rworldmap)
library(rgeos)
library(ggsn)

#Austria/czech republic
#lat y, long x
#ind.coord<- data.frame("Long"=sapply(ind.pop$V2,function(x) {data$long[data$code2 %in% x]}), "Lat"=sapply(ind.pop$V2,function(x) {data$lat[data$code2 %in% x]}))
ind.coord<- data[data$country %in% c("AT", "CZ"), c(12,11)]

boundary = extent(min(ind.coord$long)-1, max(ind.coord$long)+1, min(ind.coord$lat)-1, max(ind.coord$lat)+1)
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Plot basemap
basemap1 = ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="grey",
               colour="black", linewidth=0.2) + 
  ggtitle("Austria / Czech Republic") +
  theme_void() +
  theme(plot.title=element_text(margin=margin(t=2,b=-22, l=10), size =8))

#Serbia
#lat y, long x
#ind.coord<- data.frame("Long"=sapply(ind.pop$V2,function(x) {data$long[data$code2 %in% x]}), "Lat"=sapply(ind.pop$V2,function(x) {data$lat[data$code2 %in% x]}))
ind.coord<- data[data$country %in% "SRB", c(12,11)]

boundary = extent(min(ind.coord$long)-1, max(ind.coord$long)+1, min(ind.coord$lat)-1, max(ind.coord$lat)+1)
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Plot basemap
basemap2 = ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="grey",
               colour="black", linewidth=0.2) + 
  ggtitle("Serbia") +
  theme_void() +
  theme(plot.title=element_text(margin=margin(t=2,b=-22, l=10), size =8))

library(cowplot)
basemap = plot_grid(basemap1, NULL, basemap2, labels = c('', ''), rel_widths = c(1, -0.08, 1), nrow = 1)

ggsave(basemap, filename ="~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/BaseMap.jpeg", width=8, height=6, units="cm", dpi= 600)

#plot admixture
#remove repens
ind.pop<- ind.pop[1:80,]
admix.matrix<- admix.matrix[1:80,]
#Austria/ Czechia
ind.coord<- data.frame("Long"=sapply(ind.pop$V2,function(x) {data$long[data$code2 %in% x]}), "Lat"=sapply(ind.pop$V2,function(x) {data$lat[data$code2 %in% x]}))

ind.coord2<- ind.coord[ind.pop$V2 %in% data$code2[data$country %in% c("AT", "CZ")],]
ind.coord2<- as.matrix(ind.coord2)
admix.matrix2<- admix.matrix[ind.pop$V2 %in% data$code2[data$country %in% c("AT", "CZ")],]
admix.matrix2<- as.matrix(admix.matrix2)

library(conStruct)
jpeg("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/Map_k5_ATCZ.jpeg", width=4, height=6, units="cm", res= 600)
make.admix.pie.plot(admix.matrix2, ind.coord2, add = F, radii = 1.5, layer.colors = colors, mar = c(0, 0, 0, 0)) #radii is the dimension of the pies
dev.off()

#Serbia
ind.coord2<- ind.coord[ind.pop$V2 %in% data$code2[data$country %in% "SRB"],]
ind.coord2<- as.matrix(ind.coord2)
admix.matrix2<- admix.matrix[ind.pop$V2 %in% data$code2[data$country %in% "SRB"],]
admix.matrix2<- as.matrix(admix.matrix2)

library(conStruct)
jpeg("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/Map_k5_SRB.jpeg", width=4, height=6, units="cm", res= 600)
make.admix.pie.plot(admix.matrix2, ind.coord2, add = F, radii = 1.5, layer.colors = colors, mar = c(0, 0, 0, 0)) #radii is the dimension of the pies
dev.off()

#barplot
admix.matrix<- read.table(paste0("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/simpleQ", k, ".txt"))
pop.code<- read.table("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/poplabel.txt")
ind.code<- read.table("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/indlabel.txt")
ind.pop<- merge(ind.code, pop.code)

admix.matrix$pop<- ind.pop$V2 
admix.matrix$ind<- rownames(admix.matrix)
admix.matrix$State<- sapply(admix.matrix$pop, function(x) {data$country[data$code2 %in% x]})
admix.matrix$Ploidy<- sapply(admix.matrix$pop, function(x) {data$ploidy[data$code2 %in% x]})
admix.matrix$State_Ploidy<- paste(admix.matrix$State, admix.matrix$Ploidy, sep="_")

admix.matrix<- admix.matrix[1:80,]

colnames(admix.matrix)[1:k]<- paste("K", rep(1:k), sep="")

library(tidyverse)
#reorder and create new labels
plot_data <- admix.matrix %>% 
  gather('K', 'prob', K1:paste0("K", k)) %>% 
  group_by(ind) %>% 
  mutate(likely_assignment = K[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(ind)))  %>%
  group_by(likely_assignment) %>%
  mutate(state = paste0(unique(State_Ploidy), collapse="\n"))

labs<- c(unique(plot_data$state))
names(labs)<- c(unique(plot_data$likely_assignment))

plot_data$annotation<- 1
#plot_data$Bedrock<- sapply(plot_data$pop, function(x) {data$substrate2[data$code2 %in% x]})
plot_data$Ploidy<- unlist(plot_data$Ploidy)

#plot

p1a<- ggplot(plot_data, aes(id, prob, fill = K)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = colors[1:k]) +
  facet_grid(~likely_assignment, scales = 'free_x', space = 'free_x') +  #labeller=labeller(likely_assignment=labs)) +
  theme_minimal() + 
  #labs(x = "Individuals",  y = "Ancestry") + #title = paste0("K=", k),
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.spacing.x = unit(0,"line"),
        panel.margin = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        legend.position = 'none') 

#add annotation  
a<- ggplot(plot_data, aes(x=id, y=annotation, fill=Ploidy)) +
  geom_col() +
  facet_grid(~likely_assignment, scales = 'free_x', space = 'free_x') + 
  scale_fill_manual(values=c("#1e90ff", "#ffa500")) +
  theme_minimal() + 
  theme(legend.position = 'none', 
        axis.title = element_blank(),
        axis.text=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.spacing.x = unit(0,"line"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank(), )

library(ggpubr)
jpeg("~/Desktop/PhD Project/Alyssum_2024/Results/Entropy/Barplot_k5.jpeg", width=16, height=6, units="cm", res= 600)
plot_grid(p1a, NULL, a, ncol=1, rel_heights = c(10, -0.8, 3), align = 'v', axis='lr')
dev.off()