# load('esets.RData')
# saveRDS(X.t, "X.t.rds")
# rm(X.t)

setwd('~/Desktop/')
X.t <- readRDS("X.t.rds")

dim(X.t)

############ Clustering ############
p<-fviz_nbclust(X.t, kmeans, method = "wss")
p
# api_create(p, filename = "opt_cluster_tss_boot1")
# 
# p<-fviz_nbclust(X.t, kmeans, method = "silhouette")
# p<-ggplotly(p)
# # api_create(p, filename = "opt_cluster_silhouette_boot1")
# 
# p<-fviz_nbclust(X.t, kmeans, method = "gap_stat", nboot = 10)
# p<-ggplotly(p)