# library(Biobase)
# library(GEOquery)
library(dplyr)
library(ggfortify)

setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
#setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
#getwd()

rm(list=setdiff(ls(), "x"))
load('esets.RData')


dim(e.set) # megaExp data
dim(status) # megaExp labels

dim(e.set.i) # iris data
dim(status.iris) # iris labels

# attributes(status)
# status$Diagnosis
# status$category

# transpose
e.set.t <- t(e.set)
e.set.i.t <- t(e.set.i)

dim(e.set.t)
dim(e.set.i.t)

# rm(targets, targets.iris, e.set, e.set.i)
# e.set.t[1:5,(ncol(e.set.t)-2):ncol(e.set.t)]

# scale
e.set.s <- scale(e.set.t)
e.set.i.s <- scale(e.set.i.t)

# e.set.i.s <- as.data.frame(scale(e.set.i.t))
# e.set.s <- as.data.frame(scale(e.set.t))
# e.set.i.s <- as.data.frame(scale(e.set.i.t))

# apply(e.set.s, 2, var)
# apply(e.set.s, 2, mean)


# labels
label.i <- rownames(e.set.i.t)
label.s.i <- status.iris$My_code
common <- match(label.i, label.s.i)
# common

e.set.l <- as.character(status$most_general)
e.set.i.l <- as.character(status.iris$most_general[common]) # pass common index to extract common iris labels

e.set.df <- as.data.frame(e.set.s)
e.set.i.df <- as.data.frame(e.set.i.s)

e.set.df$label <- e.set.l
e.set.df$exp <- 0

e.set.i.df$label <- e.set.i.l
e.set.i.df$exp <- 1

e.set.df[1:5, (ncol(e.set.df)-4):ncol(e.set.df)]
e.set.i.df[1:5, (ncol(e.set.i.df)-4):ncol(e.set.i.df)]


# combine datasets
common.genes <-intersect(colnames(e.set.df), colnames(e.set.i.df))
length(common.genes) # 42004

dim(e.set.df[,common.genes])
dim(e.set.i.df[,common.genes])

e.set.df[1:10,(ncol(e.set.df)-4):ncol(e.set.df)]
e.set.i.df[1:10,(ncol(e.set.i.df)-4):ncol(e.set.i.df)]

e.set.c <- rbind(e.set.df[,common.genes], e.set.i.df[,common.genes])

dim(e.set.c)
e.set.c[1:10,(ncol(e.set.c)-5):ncol(e.set.c)]
e.set.c[(nrow(e.set.c)-10):nrow(e.set.c),(ncol(e.set.c)-5):ncol(e.set.c)]



# View(status)
e.set.c$label == 'bacterial'
e.set.c$label == 'viral'

e.set.c$label <- ifelse(e.set.c$label == 'bacterial', 1, 0)


# PCA
gset.pca <- prcomp(e.set.c, scale = FALSE)
summary(gset.pca)
plot(gset.pca, type = 'l')
autoplot(gset.pca, data = e.set.c, colour = 'label')
autoplot(gset.pca, data = e.set.c, colour = 'exp')







## DGE
bct <- e.set.c[e.set.c$label == 'bacterial',]
vrl <- e.set.c[e.set.c$label == 'viral',]

# strip the labels for raw data
bct.r <- bct[-((ncol(bct)-1):ncol(bct))]
vrl.r <- vrl[-((ncol(vrl)-1):ncol(vrl))]
dim(bct.r)
dim(vrl.r)

# 
# dim(bct)
# bct[1:10,(ncol(bct)-5):ncol(bct)]
# bct[(nrow(bct)-10):nrow(bct),(ncol(bct)-5):ncol(bct)]
# 
# vrl[1:10,(ncol(vrl)-5):ncol(vrl)]
# vrl[(nrow(vrl)-10):nrow(vrl),(ncol(vrl)-5):ncol(vrl)]

bct.mean <- apply(bct.r, 2, mean)
vrl.mean <- apply(vrl.r, 2, mean)

head(bct.mean)
head(vrl.mean)

# Just get the maximum of all the means
limit = max(bct.mean, vrl.mean)

# Scatter plot
plot(vrl.mean ~ bct.mean, xlab = "bct", ylab = "vrl",
     main = "Bct vs Vrl - Scatter", xlim = c(0, limit), ylim = c(0, limit))

# Compute fold-change (biological significance)
fold = bct.mean - vrl.mean

# Histogram of the fold differences
hist(fold, col = "gray")

# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

i <- 1
x = as.numeric(bct.r[,i])
y = as.numeric(vrl.r[,i])
n1 <- length(x)
n2 <- length(y)
n.df <- min(n1, n2)

# point estimate is the difference between the average of the two vectors
pe <- mean(x)-mean(y)
se <- sqrt((var(x)/n1 + var(y)/n2))
t_stat <- pe/se

# df minimum of n1 or n2 -1. However as there is diff variance in expression between samples
# r will use welch correction when performing t-test
df <- n-1 
df <- 125

p_val <- (1-pt(t_stat, df))*2

# 95 % conf interval
test_stat <- qt(0.025, df)
me <- test_stat * se
ci <- c((pe + me), (pe - me))

cat(paste('Homebrew T-Test',
          '\nt_stat:', round(t_stat,3),
          '\ndf:', df,
          '\np_val:', round(p_val,10),
          '\n95% CI:', round(ci[1],3), round(ci[2],3)
))

t.test(x, y)





# saveRDS(, file = 'gset_GSE72809')
#asdf