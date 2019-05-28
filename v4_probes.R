
getwd()
setwd('/home/patrick/Downloads')
ls()

l <- readLines("HumanHT-12_V4_0_R2_15002873_B.txt")
dim(l)
length(l)
l[[11]]


my_data <- read.delim("HumanHT-12_V4_0_R2_15002873_B.txt")

my_data <- read.table("HumanHT-12_V4_0_R2_15002873_B.txt",  sep = "|", fill=TRUE)
dim(my_data)
class(my_data)

my_data[1]



