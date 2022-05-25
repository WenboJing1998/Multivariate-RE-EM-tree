library(multiREEMtree)
#
# Code for manuscript data analysis and figures
#
data(MultiPoverty)
Datapaper <- MultiPoverty
#
# Shorten variable names for plotting purposes
#
colnames(Datapaper)[c(7,9,17,19)] <- c("Corrup","Democracy","Pop.dens","Rural.pop")
formula <- as.formula(Access.to.electricity + Access.to.drinking.water + Survival.rate + 
    School.enrollment ~ Corrup + Deaths.in.violence.per.mille + Democracy + Education.expenditure.percentage + 
    Gini.index + Health.expenditure.percentage + Logistics.performance.index + Pop.dens + Production.index + 
    Rural.pop + Temperature.change + Unemployment.rate)
#
# 1 SE-based tree using 10-fold CV and surrogate split
#
set.seed(110)
poverty.1se.fit <- multiREEMtree(formula, data = Datapaper, random = ~1 | Country, lme.algorithm="nlme",
   tree.xv="1se", tree.xval=10, tree.control = rpart.control(cp=0.01, usesurrogate=2, maxsurrogate = 5), stand.y = "marginal")
plot(poverty.1se.fit)
poverty.1se.fit$Tree
summary(poverty.1se.fit$Tree)
#
# Identifying country-years in each terminal node
#
tree_world <- poverty.1se.fit$Tree
terminal_nodes <- with(tree_world, names(table(where))[table(where) %in%  with(frame, n[var=="<leaf>"])])
table(tree_world$where)
#
# Note: plot is edited to add node labels and remove extraneous numbers
#
pdf("fignogdptree.pdf", paper="USr",width=8.5, height=11, pointsize=10)
plot(poverty.1se.fit)
dev.off()
#
# Construct the "world" data set to tie countries to nodes
#
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(geos)
library(egg)
library(pheatmap)

world <- ne_countries(scale = "medium", returnclass = "sf")
#
# Add "years in each terminal node" into the "world" dataset 
world[  ,c(paste("years in node", terminal_nodes))] <- NA
#
for (i in 1:nrow(world)) {
  if (world$name[i] %in% Datapaper$Country) {
    nodes_table <- table(tree_world$where[Datapaper$Country == world$name[i]])
    world[i ,c(paste("years in node", names(nodes_table)))] <- nodes_table
    world[i ,c(paste("years in node", setdiff(terminal_nodes, names(nodes_table))))] <- 0
  }
}
#
# Number of years each country falls in each terminal node
#
dfworld=as.data.frame(world)
na.omit(dfworld[,c(65:70)])
#
# Show the number of years in which a country falls into a terminal node on the world map, where white means 0 years and 
# gray means not being included in the dataset. 
#
# Not plotting node 4, since it is only Mauritius 2012-2017
#
gglist_treenodes <- list(ggplot(data = world) +
                           ggtitle("Years in node 1") +
                           geom_sf(aes(fill = `years in node 4`)) + 
                           scale_fill_gradient(low="white", high="red")+
                           theme(legend.title = element_blank()),
                         ggplot(data = world) +
                           ggtitle("Years in node 2") +
                           geom_sf(aes(fill = `years in node 6`)) + 
                           scale_fill_gradient(low="white", high="red")+
                           theme(legend.title = element_blank()),
                         ggplot(data = world) +
                           ggtitle("Years in node 3") +
                           geom_sf(aes(fill = `years in node 7`)) + 
                           scale_fill_gradient(low="white", high="red")+
                           theme(legend.title = element_blank()),
                         ggplot(data = world) +
                           ggtitle("Years in node 5") +
                           geom_sf(aes(fill = `years in node 10`)) + 
                           scale_fill_gradient(low="white", high="red")+
                           theme(legend.title = element_blank()),
                        ggplot(data = world) +
                           ggtitle("Years in node 6") +
                           geom_sf(aes(fill = `years in node 11`)) + 
                           scale_fill_gradient(low="white", high="red")+
                           theme(legend.title = element_blank())                         )


pdf("figmapnodesnogdp.pdf", width=11, height=8.5)   
egg::ggarrange(plots=gglist_treenodes, heights = c(.5, .5, .5))
dev.off()
