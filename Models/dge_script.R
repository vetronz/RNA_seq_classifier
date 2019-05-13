# Read the data into R
# library (RCurl)
library(Biobase)
library(GEOquery)
library(ggplot2)
library(tidyr)

# url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
# data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Check the loaded dataset
dim(data) # Dimension of the dataset
typeof(data)
head(data) # First few rows
tail(data) # Last few rows

typeof(data) # double which is still indexable somehow

# load series and platform data from GEO
gset <- getGEO("GSE5583", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL81", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gset_df <- as.data.frame(exprs(gset))
dim(gset_df)
typeof(gset_df)
cols <- colnames(gset_df)
cols[1]
# gset_df[[cols[1]]]
typeof(gset_df[[cols[1]]])

ggplot(gather(gset_df), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

gset_l <- log(gset_df)
ggplot(gather(gset_l), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

ggplot(stack(gset_l), aes(x = ind, y = values)) +
  geom_boxplot()

# split dataframe into conditions
wt = gset_l[,1:3]
ko = gset_l[,4:6]
head(wt)
head(ko)

# wt_1 = exprs(gset[,1])
# wt_2 = exprs(gset[,2])

# return two double vectors with mean expression of the 3 genes
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)
head(wt.mean)
head(ko.mean)

# Just get the maximum of all the means
limit = max(wt.mean, ko.mean)

# Scatter plot
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO",
     main = "GSE5583 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "red")


# Compute fold-change (biological significance)
# Difference between the means of the conditions
fold = wt.mean - ko.mean

# Histogram of the fold differences
hist(fold, col = "gray")

# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics


## T TESTING ##
# construct t-test to check if there is statistically significant gene expression
# between the first two genes

# get the 3 gene expression levels for a certain gene in each of 2 conditions
x <- NULL
ncol(wt)
for (i in 1: ncol(wt)){
  val <- wt[1,][[i]]
  x[i] <- val
}

y <- NULL
for (i in 1: ncol(ko)){
  val <- ko[1,][[i]]
  y[i] <- val
}

x
y

# define number of reps for each condition which we need for df
n1 <- length(x)
n2 <- length(y)

# point estimate is the difference between the average of the two vectors
pe <- mean(x)-mean(y)

# standard error for independent samples 
se <- sqrt((var(x)/n + var(y)/n))

# t score forumla
t_stat <- pe/se


# df minimum of n1 or n2 -1. However as there is diff variance in expression between samples
# r will use welch correction when performing t-test
df <- n1-1 

p_val <- pt(t_stat, df)*2

# 95 % conf interval
test_stat <- qt(0.025, df)
me <- test_stat * se
ci <- c((pe + me), (pe - me))


cat(paste('Homebrew T-Test',
          '\nt_stat:', t_stat,
          '\ndf:', df,
          '\np_val:', p_val,
          '\n95% CI:', ci[1], ci[2]
          ))

t <- t.test(x,y)


for(i in 1 : nrow(gset_l)) { # For each gene : 
  x = wt[i,] # WT of gene number i
  y = ko[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}


# Histogram of p-values (-log10)
hist(-log10(pvalue), col = "gray")
hist(pvalue, col = "gray")

# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot(fold, -log10(pvalue), main = "GSE5583 - Volcano")

fold_cutoff = 2
pvalue_cutoff = 0.01
abline(v = fold_cutoff, col = "blue", lwd = 2)
abline(v = -fold_cutoff, col = "red", lwd = 2)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 2)

# Screen for the genes that satisfy the filtering criteria

# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(gset_l[filter_by_fold, ])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(gset_l[filter_by_pvalue, ])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue

filtered = gset_l[filter_combined,]
dim(filtered)

head(filtered)



# Let's generate the volcano plot again,
# highlighting the significantly differential expressed genes
plot(fold, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (fold[filter_combined], -log10(pvalue[filter_combined]),
	pch = 16, col = "red")

# Highlighting up-regulated in red and down-regulated in blue
plot(fold, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (fold[filter_combined & fold < 0],
	-log10(pvalue[filter_combined & fold < 0]),
	pch = 16, col = "red")
points (fold[filter_combined & fold > 0],
	-log10(pvalue[filter_combined & fold > 0]),
pch = 16, col = "blue")





# Cluster the rows (genes) & columns (samples) by correlation
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
rowv
colv
plot(rowv)
plot(colv)

mat_data <- data.matrix(data[,2:ncol(data)])


rnames <- 
# Generate a heatmap

# convert to matrix for heatmap
mat_data <- data.matrix(filtered[,1:ncol(filtered)])

rownames(mat_data)

heatmap(mat_data)


# install.packages("gplots")		# Uncomment if not already installed
# install.packages("RColorBrewer")	# Uncomment if not already installed

library(gplots)

# Enhanced heatmap
heatmap(mat_data, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row")


#############################################################################

### NOT GOT THIS TOGETHER ###


n = nrow(filtered)

cor.table = NULL
x = NULL
y = NULL
cor.val = NULL
cor.sig = NULL

cor.test(x_exps,y_exps)


for (i in 1 : (n-1)) {
	x_name = rownames(filtered)[i]
	x_exps = filtered[i, ]	

	for (j in (i+1) : n) {
		y_name = rownames(filtered)[j]
		y_exps = filtered[j, ]
		
		output = cor.test(x_exps,y_exps)
		
		x = c(x, x_name)
		y = c(y, y_name)
		cor.val = c(cor.val, output$estimate)
		cor.sig = c(cor.sig, output$p.value)
	}
}

cor.table = data.frame (x, y, cor.val, cor.sig)

dim(cor.table)
head(cor.table)

sig_cutoff = 0.001

cor.filtered = subset (cor.table, cor.sig < sig_cutoff)

dim(cor.filtered)
head(cor.filtered)













