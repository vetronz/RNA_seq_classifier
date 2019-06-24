library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms

setwd('~/Desktop/')


X.t <- readRDS("X.t.rds")
df.1 <- X.t

set.seed(123)
# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(df.1, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

wss_values <- c('2423019', '2129208', '1972591', '1867830', '1794755', '1742501', '1701956', '1668217', '1632186', '1602538', '1576787', '1549114', '1533915', '1509435', '1484579')


# saveRDS(wss_values, file = "wss_values.rds")

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


# gap_stat <- clusGap(X.t, FUN = kmeans, nstart = 25,
#                     K.max = 10, B = 20)

setwd('~/Documents/Masters/RNA_seq_classifier/Data/')
gap_full <- readRDS("gap_stat_full.rds")
class(gap_full)

gap_full$Tab[,3]
plot(1:10, gap_full$Tab[,3])
