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

# Method: Multivariate RE-EM tree
#
# 1 SE-based tree using 10-fold CV and surrogate split. Note: plot is edited to add node labels and remove extraneous numbers.
#
set.seed(110)
poverty.1se.fit <- multiREEMtree(formula, data = Datapaper, random = ~1 | Country, lme.algorithm="nlme",
   tree.xv="1se", tree.xval=10, stand.y = "none", tree.control = rpart.control(cp=0.01, usesurrogate=2, maxsurrogate = 5))

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
                           theme(legend.title = element_blank())
                        )


#
# Method: MRT
#

#
# Shorten variable names for plotting purposes
#
Datapaper2 <- MultiPoverty
colnames(Datapaper2)[c(7,9,15,17,19,21)] <- c("Corrup","Democracy","Health.exp.pct","Pop.dens","Rural.pop","Unemp")
formula <- as.formula(Access.to.electricity + Access.to.drinking.water + Survival.rate +
                        School.enrollment ~ Corrup + Deaths.in.violence.per.mille + Democracy + Education.expenditure.percentage +
                        Gini.index + Health.exp.pct + Logistics.performance.index + Pop.dens + Production.index +
                        Rural.pop + Temperature.change + Unemp)
#
# 1 SE-based mvpart tree using 10-fold CV and surrogate split
#
set.seed(110)
poverty.mvpart.1se.fit <- multitree(formula, data = Datapaper2, stand.y = "none",
                                    tree.xv="1se", tree.xval=10, tree.control = rpart.control(cp=0.01, usesurrogate=2, maxsurrogate = 5))

poverty.mvpart.1se.fit
plot(poverty.mvpart.1se.fit$Tree)
text(poverty.mvpart.1se.fit$Tree)
summary(poverty.mvpart.1se.fit$Tree)

#
# Note: plot in paper is edited to add node labels and remove extraneous numbers
#

getTerminal(poverty.mvpart.1se.fit)
tnodeREEM <- getTerminal(poverty.1se.fit)
tnodemvpart <- getTerminal(poverty.mvpart.1se.fit)
table(tnodeREEM$terminal_nodes_where, tnodemvpart$terminal_nodes_where)

tree_world2 <- poverty.mvpart.1se.fit$Tree
terminal_nodes2 <- with(tree_world2, names(table(where))[table(where) %in%  with(frame, n[var=="<leaf>"])])

world2 <- ne_countries(scale = "medium", returnclass = "sf")
#
# Add "years in each terminal node" into the "world" dataset
world2[  ,c(paste("years in node", terminal_nodes2))] <- NA
#
for (i in 1:nrow(world2)) {
  if (world2$name[i] %in% Datapaper2$Country) {
    nodes_table2 <- table(tree_world2$where[Datapaper2$Country == world2$name[i]])
    world2[i ,c(paste("years in node", names(nodes_table2)))] <- nodes_table2
    world2[i ,c(paste("years in node", setdiff(terminal_nodes2, names(nodes_table2))))] <- 0
  }
}
#
# Number of years each country falls in each terminal node
#
dfworld2=as.data.frame(world2)
#
# Show the number of years in which a country falls into a terminal node on the world map, where white means 0 years and
# gray means not being included in the dataset.
#
# Not plotting node 4, since it is only Bhutan 2012-2017
#
gglist_treenodes2 <- list(ggplot(data = world2) +
                            ggtitle("Years in node 1") +
                            geom_sf(aes(fill = `years in node 5`)) +
                            scale_fill_gradient(low="white", high="red")+
                            theme(legend.title = element_blank()),
                          ggplot(data = world2) +
                            ggtitle("Years in node 2") +
                            geom_sf(aes(fill = `years in node 6`)) +
                            scale_fill_gradient(low="white", high="red")+
                            theme(legend.title = element_blank()),
                          ggplot(data = world2) +
                            ggtitle("Years in node 4") +
                            geom_sf(aes(fill = `years in node 9`)) +
                            scale_fill_gradient(low="white", high="red")+
                            theme(legend.title = element_blank()),
                          ggplot(data = world2) +
                            ggtitle("Years in node 5") +
                            geom_sf(aes(fill = `years in node 11`)) +
                            scale_fill_gradient(low="white", high="red")+
                            theme(legend.title = element_blank()),
                          ggplot(data = world2) +
                            ggtitle("Years in node 6") +
                            geom_sf(aes(fill = `years in node 12`)) +
                            scale_fill_gradient(low="white", high="red")+
                            theme(legend.title = element_blank()),
                          ggplot(data = world2) +
                            ggtitle("Years in node 7") +
                            geom_sf(aes(fill = `years in node 14`)) +
                            scale_fill_gradient(low="white", high="red")+
                            theme(legend.title = element_blank()),
                          ggplot(data = world2) +
                            ggtitle("Years in node 8") +
                            geom_sf(aes(fill = `years in node 15`)) +
                            scale_fill_gradient(low="white", high="red")+
                            theme(legend.title = element_blank())                         )

egg::ggarrange(plots=gglist_treenodes2, heights = c(.5, .5, .5, .5))

#
# Method: Separate REEM trees
#

library(REEMtree)

formula.elec <- as.formula(Access.to.electricity ~ Corrup + Deaths.in.violence.per.mille + Democracy + Education.expenditure.percentage +
                             Gini.index + Health.exp.pct + Logistics.performance.index + Pop.dens + Production.index +
                             Rural.pop + Temperature.change + Unemp)
formula.water <- as.formula(Access.to.drinking.water ~ Corrup + Deaths.in.violence.per.mille + Democracy + Education.expenditure.percentage +
                              Gini.index + Health.exp.pct + Logistics.performance.index + Pop.dens + Production.index +
                              Rural.pop + Temperature.change + Unemp)
formula.survival <- as.formula(Survival.rate ~ Corrup + Deaths.in.violence.per.mille + Democracy + Education.expenditure.percentage +
                                 Gini.index + Health.exp.pct + Logistics.performance.index + Pop.dens + Production.index +
                                 Rural.pop + Temperature.change + Unemp)
formula.school <- as.formula(School.enrollment ~ Corrup + Deaths.in.violence.per.mille + Democracy + Education.expenditure.percentage +
                               Gini.index + Health.exp.pct + Logistics.performance.index + Pop.dens + Production.index +
                               Rural.pop + Temperature.change + Unemp)

set.seed(110)
poverty.REEM.elec <- REEMtree(formula.elec, random = ~1 | Country, data=Datapaper2)
library(partykit)
plot(as.party(poverty.REEM.elec$Tree), type="simple")
poverty.REEM.water <- REEMtree(formula.water, random = ~1 | Country, data=Datapaper2)
plot(as.party(poverty.REEM.water$Tree), type="simple")
poverty.REEM.survival <- REEMtree(formula.survival, random = ~1 | Country, data=Datapaper2)
plot(as.party(poverty.REEM.survival$Tree), type="simple")
poverty.REEM.school <- REEMtree(formula.school, random = ~1 | Country, data=Datapaper2)
plot(as.party(poverty.REEM.school$Tree), type="simple")
