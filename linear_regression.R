### EXAMPLE ONE ###

# define our two vectors. feature and indep var
x <- c(2.1,2.5,3.6,4)
y<-c(8,10,12,14)


n<-length(x)

plot(x,y)

# calc covariance
x_spread <- x-mean(x)
y_spread <- y-mean(y)

cross_prod <- x_spread*y_spread
tot <- sum(cross_prod)

my_cov <- tot/(n-1)
my_cov
# cov(x,y)

## correlation ##
my_cor <- my_cov/(sd(x)*sd(y))
my_cor
# cor(x,y)

# define our regression coefficient and intercept terms using equations derived from 1st derivatives
m <- my_cor*sd(y)/sd(x)
c <- mean(y-(m*x))
print(paste0('Coefficient: ', m))
print(paste0('intercept: ', c))

# plot
abline(c, m)

# vertical line at x = 3
abline(v=3)
abline(h=(3*m+c))

# vertical line at x = 4
abline(v=4)
abline(h=(4*m+c))

# difference of 1 increment in x dimension adds additional m in y dimension
3*m+c+m
4*m+c


# regression model summary
model <- lm(y ~ x)
summary(model)

# check against our own values
# the standard error is the root mean squared error except we take 2 from n as a correction

pred <- function(x) x*m+c
v_pred <- pred(x)
e <- v_pred - y

sem <- sum((e)^2)/(n-2) # method 1
(sum((e)^2)/(n-1))^(1/2) # method 2
(sum((e)^2)/(n-2))^(1/2) # method 3, called the residual standard error on the model output


t_score <- m/sem
df <- 3-2

pt(t_score, df, lower.tail = FALSE)*2

# 95 % conf interval
t_stat <- qt(0.025,25)
me <- t_stat * sem
m + me
m - me



############################################################################################



### EXAMPLE TWO ###

v = 1:100
v2 = sample((v/100), 100, replace = TRUE)
plot(v, v2)

# sum(v-mean(v) * (v2-mean(v2)))/(length(v)-1)

v_spread <- v-mean(v)
v2_spread <- v2-mean(v2)

cross_prod <- v_spread*v2_spread
# cross_prod
tot <- sum(cross_prod)
tot

my_cov <- tot/(length(v)-1)
my_cov
cov(v,v2)

## CORRELATION ##
my_cor <- my_cov/(sd(v)*sd(v2))
my_cor
cor(v,v2)

m <- my_cor*sd(v2)/sd(v)
c <- mean(v2-(m*v))
abline(c, m)
print(paste0('Coefficient: ', m))
print(paste0('intercept: ', c))











