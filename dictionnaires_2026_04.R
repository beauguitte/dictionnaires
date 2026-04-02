# Laurent Beauguitte, Dictionnaires de géographie, données relationnelles, avril 2026

rm(list=ls())

## Propriétés des ouvrages

## Import liste de liens
ge70l <- read.delim("George1970_liens.txt") 
br92l <- read.delim("Brunetetal1992_liens.txt") 
le03l <- read.delim("Levy_Lussault2003_liens.txt")
cy20l <- read.delim("Cynorhodon2020_liens.txt")
co25l <- read.delim("GeoXXI2025_liens.txt")

# Nombre de notices
length(unique(ge70l$ENTREE))
length(unique(br92l$ENTREE))
length(unique(le03l$ENTREE))
length(unique(cy20l$ENTREE))
length(unique(co25l$ORIGINE))

## Nombre d'auteurices, hiérarchie
## Proportion de notices écrites à plusieurs
## Non pertinent pour George (1970, notices non signées) 
## et Brunet et al. (1993) (notices rarement signées) 

# Levy et Lussault 2003
le03a <- read.delim("Levy_Lussault2003_auteurices-sommets.txt")

# Nb d'auteurices
length(unique(le03a$NOM_Prenom))

# Nb de notices par auteurices
levy_1 <- as.data.frame(table(le03a$NOM_Prenom))
table(levy_1$Freq)
# proportion de notices écrites par le principal auteur
max(levy_1$Freq)/length(unique(le03l$ORIGINE))

# Pourcentage de notices écrites à plusieurs
levy_2 <- as.data.frame(table(le03a$ENTREE))
100*table(levy_2$Freq)/length(levy_2$Var1)

# Dictionnaire de l'anthropocène 2020
cy20a <- read.delim("Cynorhodon2020_auteurices-sommets.txt")
length(unique(cy20a$AUTEURICE))
anth_1 <- as.data.frame(table(cy20a$AUTEURICE))
table(anth_1$Freq)
max(anth_1$Freq)/length(unique(cy20l$ORIGINE))
anth_2 <- as.data.frame(table(cy20a$ORIGINE))
100*table(anth_2$Freq)/length(anth_2$Var1)

## Collectif GéoXXI, 2025
co25a <- read.delim("GeoXXI2025_auteurices-sommets.txt")
length(unique(co25a$AUTEURICE))
geo21_1 <- as.data.frame(table(co25a$AUTEURICE)) 
table(geo21_1$Freq)
max(geo21_1$Freq)/length(unique(co25l$ORIGINE))
geo21_2 <- as.data.frame(table(co25a$ENTREE))
table(geo21_2$Freq)/length(geo21_2$Var1)

# suppression des objets non utilisés ensuite
rm(anth_1, anth_2, co25a, cy20a, geo21_1, geo21_2, le03a, levy_1, levy_2)

# Analyse de réseau

library(igraph)

# Réseau des notices

# Notices sans le moindre renvoi (George et Brunet et al.)
sum(is.na(ge70l$DESTINATION))
sum(is.na(br92l$DESTINATION))

# Distribution des degrés entrants et sortants
# Composantes et réciprocité des liens
# Réeau des liens mutuels

#George 1970 - sélection origine - destination
d <- ge70l[,c(3,5)]

# si on souhaite examiner uniquement les notices avec degré sortant > 0
# george <- na.omit(george)

g <- graph_from_data_frame(d, directed= TRUE)

# supprimer sommet NA
g <- induced_subgraph(g, v = which(V(g)$name != "NA"))
is_simple(g)

# supprimer boucles et liens multiples
g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
is_simple(g)

# degrés entrants et sortants
V(g)$degin <- degree(g, mode = "in")
V(g)$degou <- degree(g, mode = "out")
summary(V(g)$degin)
summary(V(g)$degou)

# connexité
comp <- components(g)
# nombre d'isolés
table(comp$csize)
# proportion de sommets dans la composante principale
100* max(comp$csize) / vcount(g)

# suppression des isolés
g <- delete_vertices(g, V(g)[degree(g, mode="all") < 1])

# extraire la plus grande composante 
comp <- decompose(g)

pgcom <- comp[[1]]
plot(pgcom, 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)

# proportion de liens mutuels
dc <- dyad_census(pgcom)
100*dc$mut / dc$asym

# degré entrant
V(pgcom)$degin <- degree(pgcom, mode="in")
summary(V(pgcom)$degin)

# proportion de sommets avec degré entrant nul
100*table(V(pgcom)$degin) / vcount(pgcom)

# termes ayant le degré entrant le plus élevé
g15 <- induced_subgraph(pgcom, v = which(V(pgcom)$degin > 14))
V(g15)$name

# distribution des degrés entrants (composante principale)
plot(degree_distribution(pgcom, mode="in"), type="b")

# ego-network du terme ville
EgoNet_ville <- make_ego_graph(pgcom,
                                  nodes = V(pgcom)[name=='ville'], 
                                  order = 1,  # voisins d'ordre 1
                                  mode = c("all"))
plot(EgoNet_ville[[1]], 
     edge.arrow.size = 0.5,
     vertex.color="yellow")

# explorations complémentaires
# sélection du réseau des liens mutuels
mutg <- which_mutual(pgcom)
gmg <- delete_edges(pgcom, E(pgcom)[mutg == FALSE])

# suppression des isolés
gmg <- delete_vertices(gmg, V(gmg)[degree(gmg) < 1])
gmg

# extraire la plus grande composante connexe
comp <- decompose(gmg)

plot(comp[[3]], 
     vertex.label.cex=1, 
     vertex.size = 6, 
     edge.arrow.size = 0)

rm(list=ls())

# Brunet et al. 1993

d <- read.delim("Brunetetal1992_liens.txt") 
d <- d[,c(3,5)]

# suppression des liens vers des entrées inexistantes
d <- d[d$DESTINATION != "<entrée absente>",]
g <- graph_from_data_frame(d, directed=TRUE)
g <- induced_subgraph(g, v = which(V(g)$name != "NA"))
g<- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
is_simple(g)

# degrés entrants et sortants
V(g)$degin <- degree(g, mode = "in")
V(g)$degou <- degree(g, mode = "out")
summary(V(g)$degin)
summary(V(g)$degou)

# connexité
comp <- components(g)
# nombre d'isolés
table(comp$csize)
# proportion de sommets dans la composante principale
100* max(comp$csize) / vcount(g)

# analyse de la plus grande composante
comp <- decompose(g)
pgcom <- comp[[1]]
plot(pgcom, 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)

# proportion de liens mutuels
dc <- dyad_census(pgcom)
100*dc$mut / dc$asym

# degré entrant
V(pgcom)$degin <- degree(pgcom, mode="in")
summary(V(pgcom)$degin)

# proportion de sommets avec degré entrant nul
100*table(V(pgcom)$degin) / vcount(pgcom)

# termes ayant le degré entrant le plus élevé
g20 <- induced_subgraph(pgcom, v = which(V(pgcom)$degin > 20))
V(g20)$name

# distribution des degrés entrants (composante principale)
plot(degree_distribution(pgcom, mode="in"), type="b")

# ego-network du terme ville
EgoNet_ville <- make_ego_graph(pgcom,
                                  nodes = V(pgcom)[name=='ville'], 
                                  order = 1,  # voisins d'ordre 1
                                  mode = c("all"))
plot(EgoNet_ville[[1]], 
     edge.arrow.size = 0.5,
     vertex.color="yellow")

EgoNet_ville2 <- induced_subgraph(EgoNet_ville[[1]], 
                                 v = which(V(EgoNet_ville[[1]])$name != "ville"))
plot(EgoNet_ville2, 
     edge.arrow.size = 0.5,
     vertex.color="yellow")

# explorations complémentaires
# sélection du réseau des liens mutuels
mutg <- which_mutual(pgcom)
gmg <- delete_edges(pgcom, E(pgcom)[mutg == FALSE])

# suppression des isolés
gmg <- delete_vertices(gmg, V(gmg)[degree(gmg) < 1])
gmg

# extraire la plus grande composante connexe
comp <- decompose(gmg)

plot(comp[[3]], 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)

rm(list=ls())

# Lévy, Lussault, 2003

d <- read.delim("Levy_Lussault2003_liens.txt")

g <- graph_from_data_frame(d[,c(3,5)], directed=TRUE)
g <- induced_subgraph(g, v = which(V(g)$name != "NA"))
g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
is_simple(g)

# degrés entrants et sortants
V(g)$degin <- degree(g, mode = "in")
V(g)$degou <- degree(g, mode = "out")
summary(V(g)$degin)
summary(V(g)$degou)

# connexité
comp <- components(g)
# nombre d'isolés
table(comp$csize)
# proportion de sommets dans la composante principale
100* max(comp$csize) / vcount(g)

# analyse de la plus grande composante
comp <- decompose(g)
pgcom <- comp[[1]]
plot(pgcom, 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)

# proportion de liens mutuels
dc <- dyad_census(pgcom)
100*dc$mut / dc$asym

# degré entrant
V(pgcom)$degin <- degree(pgcom, mode="in")
summary(V(pgcom)$degin)

# proportion de sommets avec degré entrant nul
100*table(V(pgcom)$degin) / vcount(pgcom)

# termes ayant le degré entrant le plus élevé
g20 <- induced_subgraph(pgcom, v = which(V(pgcom)$degin > 20))
V(g20)$name

# distribution des degrés entrants (composante principale)
plot(degree_distribution(pgcom, mode="in"), type="b")

# ego-network du terme ville
EgoNet_ville <- make_ego_graph(pgcom,
                                  nodes = V(pgcom)[name=='Ville'], 
                                  order = 1,  # voisins d'ordre 1
                                  mode = c("all"))
plot(EgoNet_ville[[1]], 
     edge.arrow.size = 0.5,
     vertex.color="yellow",
     vertex.label=NA)

EgoNet_ville2 <- induced_subgraph(EgoNet_ville[[1]], 
                                     v = which(V(EgoNet_ville[[1]])$name != "Ville"))
plot(EgoNet_ville2, 
     edge.arrow.size = 0.2,
     vertex.color="yellow",
     vertex.size=5)

# explorations complémentaires
# sélection du réseau des liens mutuels
mutg <- which_mutual(pgcom)
gmg <- delete_edges(pgcom, E(pgcom)[mutg == FALSE])

# suppression des isolés
gmg <- delete_vertices(gmg, V(gmg)[degree(gmg) < 1])
gmg

# extraire la plus grande composante connexe
comp <- decompose(gmg)

plot(comp[[1]], 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)

rm(list=ls())

# Dictionnaire critique de l'anthropocène, 2020

d <- read.delim("Cynorhodon2020_liens.txt")

g <- graph_from_data_frame(d[,c(3,5)], directed=TRUE)
g <- induced_subgraph(g, v = which(V(g)$name != "NA"))
g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
is_simple(g)

# degrés entrants et sortants
V(g)$degin <-degree(g, mode = "in")
V(g)$degou <- degree(g, mode = "out")
summary(V(g)$degin)
summary(V(g)$degou)

# connexité
is_connected(g)

plot(g, 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)

# proportion de liens mutuels
dc <- dyad_census(g)
100*dc$mut / dc$asym

# degré entrant
V(g)$degin <- degree(g, mode="in")
summary(V(g)$degin)

# proportion de sommets avec degré entrant nul
100*table(V(g)$degin) / vcount(g)

# termes ayant le degré entrant le plus élevé
g20 <- induced_subgraph(g, v = which(V(g)$degin > 20))
V(g20)$name

# distribution des degrés entrants
plot(degree_distribution(g, mode="in"), type="b")

# ego-network des termes ville - durable et ville - intelligente
EgoNet_ville <- make_ego_graph(g,
                               nodes = V(g)[name=='ville - durable'], 
                               order = 1,  # voisins d'ordre 1
                               mode = c("all"))
plot(EgoNet_ville[[1]], 
     edge.arrow.size = 0.5,
     vertex.color="yellow")


EgoNet_ville <- make_ego_graph(g,
                               nodes = V(g)[name=='ville - intelligente'], 
                               order = 1,  # voisins d'ordre 1
                               mode = c("all"))
plot(EgoNet_ville[[1]], 
     edge.arrow.size = 0.5,
     vertex.color="yellow")

EgoNet_ville <- make_ego_graph(g,
                               nodes = V(g)[name=='ville - flottante'], 
                               order = 1,  # voisins d'ordre 1
                               mode = c("all"))
plot(EgoNet_ville[[1]], 
     edge.arrow.size = 0.5,
     vertex.color="yellow")

# explorations complémentaires
# sélection du réseau des liens mutuels
mutg <- which_mutual(g)
gmg <- delete_edges(g, E(g)[mutg == FALSE])

# suppression des isolés
gmg <- delete_vertices(gmg, V(gmg)[degree(gmg) < 1])
gmg

# extraire la plus grande composante connexe
comp <- decompose(gmg)

plot(comp[[1]], 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)

rm(list=ls())

# Géographies, un dictionnaire
d <- read.delim("GeoXXI2025_liens.txt")

g <- graph_from_data_frame(d[,c(2,4)], directed=TRUE)
g <- induced_subgraph(g, v = which(V(g)$name != "NA"))
g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
is_simple(g)

# degrés entrants et sortants
V(g)$degin <- degree(g, mode = "in")
V(g)$degou <- degree(g, mode = "out")
summary(V(g)$degin)
summary(V(g)$degou)

# connexité
is_connected(g)
comp <- components(g)

# nombre d'isolés
table(comp$csize)
# proportion de sommets dans la composante principale
100* max(comp$csize) / vcount(g)

# analyse de la plus grande composante
comp <- decompose(g)
pgcom <- comp[[1]]
plot(pgcom, 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)

# proportion de liens mutuels
dc <- dyad_census(pgcom)
100*dc$mut / dc$asym

# degré entrant
V(pgcom)$degin <- degree(pgcom, mode="in")
summary(V(pgcom)$degin)

# proportion de sommets avec degré entrant nul
100*table(V(pgcom)$degin) / vcount(pgcom)

# termes ayant le degré entrant le plus élevé
g10 <- induced_subgraph(pgcom, v = which(V(pgcom)$degin > 10))
V(g10)$name

# distribution des degrés entrants (composante principale)
plot(degree_distribution(pgcom, mode="in"), type="b")

# ego-network du terme ville
EgoNet_ville <- make_ego_graph(pgcom,
                               nodes = V(pgcom)[name=='ville'], 
                               order = 1,  # voisins d'ordre 1
                               mode = c("all"))
plot(EgoNet_ville[[1]], 
     edge.arrow.size = 0.5,
     vertex.color="yellow")

# explorations complémentaires
# sélection du réseau des liens mutuels
mutg <- which_mutual(pgcom)
gmg <- delete_edges(pgcom, E(pgcom)[mutg == FALSE])

# suppression des isolés
gmg <- delete_vertices(gmg, V(gmg)[degree(gmg) < 1])
gmg

# extraire la plus grande composante connexe
comp <- decompose(gmg)
comp

plot(comp[[1]], 
     vertex.label=NA, 
     vertex.size = 6, 
     edge.arrow.size = 0)
