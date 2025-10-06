# Install needed packages if not already installed
# Taking STRING as an example: nodes, edges can turn to network regardless where they are from
# Load libraries
library(igraph)
library(STRINGdb)
library(dplyr)

# preparation of STRING database
## threshold: low threshold: 150, meidum threshold: 400, reliable (high) threshold: 700, highest threshold: 900
set.seed(1992)
score_threshold = 700
string_db <- STRINGdb$new(version="11", species=9606, score_threshold=score_threshold, input_directory="")

# input genes: start_gene, intermediate_gene, outcome_gene
# id transition for network matching

all_genes = unique(sort(c(start_gene, intermediate_gene, outcome_gene)))
mapped <- string_db$map(data.frame(gene = all_genes), "gene", removeUnmappedRows = TRUE)
gene2id <- setNames(mapped$STRING_id, mapped$gene)
id2gene <- setNames(mapped$gene, mapped$STRING_id)

# interactions establishment for all genes (after ID transition)

full_input = gene2id[all_genes]
interactions <- string_db$get_interactions(full_input)
interactions$from <- id2gene[interactions$from]
interactions$to   <- id2gene[interactions$to]
interactions <- interactions[!is.na(interactions$from) & !is.na(interactions$to), ]

# network 1: start_gene network
network_1_interactions = interactions[((interactions$from %in% pulled_down)&(interactions$to %in% pulled_down))|
                                       ((interactions$to %in% pulled_down)&(interactions$from %in% pulled_down)),]
network_1 = graph_from_data_frame(network_1_interactions, directed = TRUE) # directed or not
network_1 = simplify(network_1, remove.multiple = TRUE, remove.loops = TRUE)

# network 2: extend start_gene network to intermediate_gene network

network_2_interactions = rbind(network_1_interactions,
                                interactions[((interactions$from %in% start_gene)&
                                                (interactions$to %in% intermediate_gene))|
                                            ((interactions$to %in% intermediate_gene)&
                                                (interactions$from %in% start_gene)),])
network_2_interactions = network_2_interactions[!duplicated(network_2_interactions),]
network_2 = graph_from_data_frame(network_2_interactions, directed = TRUE) # directed or not
network_2 = simplify(network_2, remove.multiple = TRUE, remove.loops = TRUE)

## coloring network 2, nodes, network 1

V(network_2)$color = "skyblue"
V(network_2)$color[V(network_2)$name %in% 
    unique(sort(intersect(start_gene,intermediate_gene)))] = "red"
V(network_2)$color[V(network_2)$name %in% 
    unique(sort(setdiff(intermediate_gene,start_gene)))] = "green"

# network 3: linking start_gene/intermediate_gene network further to outcome_gene

## finding nodes from intermediate_gene

part_1_interactions = interactions[((interactions$from %in% start_gene)&
                                                (interactions$to %in% intermediate_gene))|
                                            ((interactions$to %in% intermediate_gene)&
                                                (interactions$from %in% start_gene)),]
selected_oxphos = unique(sort(intersect(OXPHOS_gene,c(part_1_interactions$from,part_1_interactions$to))))
part_2_interactions = interactions[((interactions$from %in% selected_oxphos)&
                                              (interactions$to %in% hf_gene))|
                                       ((interactions$to %in% hf_gene)&
                                              (interactions$from %in% selected_oxphos)),]
network_3_interactions = rbind(network_2_interactions,part_2_interactions)
network_3_interactions = network_3_interactions[!duplicated(network_3_interactions),]
network_3 = graph_from_data_frame(network_3_interactions, directed = TRUE)
network_3 = simplify(network_3, remove.multiple = TRUE, remove.loops = TRUE)
V(network_3)$color = "skyblue"
V(network_3)$color[V(network_3)$name %in% unique(sort(intersect(OXPHOS_gene,pulled_down)))] = "red"
V(network_3)$color[V(network_3)$name %in% unique(sort(setdiff(OXPHOS_gene,pulled_down)))] = "green"
V(network_3)$color[V(network_3)$name %in% unique(sort(intersect(hf_gene,OXPHOS_gene)))] = "red"
V(network_3)$color[V(network_3)$name %in% unique(sort(setdiff(hf_gene,c(OXPHOS_gene,pulled_down))))] = "purple"

# standard visualization function 

par(mfrow = c(1,3))
plot(network_1,
     vertex.size = 10,
     vertex.label.cex = 0.5,
     vertex.label.color = "black",
     edge.arrow.size   = 0.2,
     edge.color = "gray50",
     layout = layout_with_fr, main = "Start Network")

plot(network_2,
     vertex.size = 8,
     vertex.label.cex = 0.4,
     vertex.label.color = "black",
     edge.arrow.size   = 0.4,
     edge.color = "gray50",
     layout = layout_with_fr, main = "Start-Intermediate Network")

legend("topleft", legend=c("Start","Bridging nodes","Intermediate"),
         col=colors, pch=c(rep(16,3),NA),
         lty=c(NA,NA,NA,1), pt.cex=c(2,2,2,NA), bty="n")

plot(network_3,
     vertex.size = 4,
     vertex.label.cex = 0.2,
     vertex.label.color = "black",
     edge.arrow.size   = 0.2,
     edge.color = "gray50",
     layout = layout_with_fr, main = "Start-Intermediate-Outcome Network")

legend("topleft", legend=c("Start","Bridging nodes","Intermediate","Outcome"),
         col=colors, pch=c(rep(16,4),NA),
         lty=c(NA,NA,NA,NA,1), pt.cex=c(2,2,2,2,NA), bty="n")

