# NOTE: the current version of GPfit - version 0.2-0 is more flexible
# and more stable than version 0.1-0. 

# 1. The most important update:
#    GP_fit(x,y) # in version 0.1-0. fits a GP model with Gaussian correlation
#
#    In version 0.2-0, one can specify the correlation structure as Matern or
#    power exponential with different power. For instance, we used 
#    GP_fit(x, y,corr = list(type = "exponential", power = 2)) to produce the plots and
#    simulation results presented in this paper.
#
#    The default is set to: corr = list(type = "exponential", power = 1.95)
#
# 2. A few default plotting characters have been updated for better visuals
#    the defaults in version 0.1-0 was: 
#         "line_type = c(1, 1), legends = TRUE, cex = 2, pch = 1"
#    the defaults in version 0.2-0 are:
#         "line_type = c(1, 2), legends = FALSE, cex = 1, pch = 20"
#
# 3. The version 1.0-0 is the same as version 0.2-0 except the citation of 
#    this JSS paper 
#
# Loads the package (assuming that it is installed on the machine)


#=================================
# Likelihood plots - Fig 1
#=================================

library("GPfit", "lhs")

#-- 1d test function -------------
computer_simulator_1 <- function(x) {
  y <- log(x + 0.1) + sin(5 * pi * x)
  return(y)
}
#--------------------

n <- 10
d <- 1
set.seed(12)
# x = maximinLHS(n,d)
x <- c(0.6381151, 0.30458841, 0.51526238, 0.17747182, 0.83945864, 0.73615131, 0.242109, 
  0.43220115, 0.95521608, 0.09712655)
y <- computer_simulator_1(x)

# Figure 1(a)
theta_lb <- 10^(-15)
theta_ub <- 70
resolution <- 100

thetavector <- seq(from = theta_lb, to = theta_ub, length.out = resolution)
betavector <- matrix(log10(thetavector), resolution, 1)

likeval <- matrix(0, resolution, 1)
# for(i in 1:resolution){likeval[i] = GP_deviance(betavector[i],x,y)} # for version 0.1-0
for (i in 1:resolution) {
  likeval[i] <- GP_deviance(betavector[i], x, y, corr = list(type = "exponential", 
    power = 2))
}

pdf("fig1a.pdf", height = 5, width = 5)
plot(thetavector, likeval, type = "l", xlab = expression(theta), ylab = expression(-2 * 
  log(L[theta])))
dev.off()

# Figure 1(b)
theta_lb <- 10^(-15)
theta_ub <- 5
resolution <- 100

thetavector <- seq(from = theta_lb, to = theta_ub, length.out = resolution)
betavector <- matrix(log10(thetavector), resolution, 1)

likeval <- matrix(0, resolution, 1)
# for(i in 1:resolution){likeval[i] = GP_deviance(betavector[i],x,y)} # for version 0.1-0
for (i in 1:resolution) {
  likeval[i] <- GP_deviance(betavector[i], x, y, corr = list(type = "exponential", 
    power = 2))
}

pdf("fig1b.pdf", height = 5, width = 5)
plot(thetavector, likeval, type = "l", xlab = expression(theta), ylab = expression(-2 * 
  log(L[theta])))
dev.off()


#-- 2d test function -------------
computer_simulator_2 <- function(x) {
  x <- 4 * x - 2
  x1 <- x[, 1]
  x2 <- x[, 2]
  
  t1 <- (1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 6 * x1 * x2 + 
    3 * x2^2))
  t2 <- (30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 36 * 
    x1 * x2 + 27 * x2^2))
  y <- t1 * t2
  
  return(y)
}
#--------------------

library("lattice", "lhs", "colorspace")
#---------------------
n <- 30
d <- 2
set.seed(1)
x <- matrix(0, n, 2)

# x = maximinLHS(n,d)
x[, 1] <- c(0.7937935261, 0.4496580769, 0.2752867432, 0.6913658896, 0.8777401649, 
  0.3987284473, 0.7497571945, 0.1601039682, 0.0527092526, 0.5000438219, 0.6361983972, 
  0.1778581875, 0.3197865148, 0.4003190597, 0.7277299663, 0.0757593439, 0.8408587278, 
  0.4686633414, 0.8024329295, 0.3507487091, 0.9231584449, 0.0005985916, 0.5385867361, 
  0.2234647278, 0.6108420338, 0.5994940089, 0.9631190429, 0.1241414898, 0.9833982868, 
  0.2536672005)
x[, 2] <- c(0.17412559, 0.33371134, 0.78181996, 0.56635533, 0.66482638, 0.16330925, 
  0.09283118, 0.48925892, 0.8443556, 0.6021934, 0.23377032, 0.80397548, 0.86791042, 
  0.43871987, 0.75889476, 0.27293878, 0.93539901, 0.30493853, 0.20142953, 0.72630772, 
  0.4022192, 0.58140528, 0.123962, 0.06284088, 0.53236272, 0.90130616, 0.02741331, 
  0.97629734, 0.6809449, 0.39753711)

y <- computer_simulator_2(x)


# Figure 1(c)
theta_lb <- 10^(-5)
theta_ub <- 13
resolution <- 40
range <- c(theta_lb, theta_ub)

thetavector <- seq(from = theta_lb, to = theta_ub, length.out = resolution)
thetavec <- expand.grid(x = thetavector, y = thetavector)
betavec <- as.matrix(log10(thetavec))

likeval <- matrix(0, resolution^2, 1)
# for(i in 1:(resolution^2)){likeval[i] = GP_deviance(betavec[i,],x,y)} # for version 0.1-0
for (i in 1:(resolution^2)) {
  likeval[i] <- GP_deviance(betavec[i, ], x, y, corr = list(type = "exponential", 
    power = 2))
}
dim(likeval) <- c(length(thetavector), length(thetavector))

pdf("fig1c.pdf", height = 5, width = 5)
mypalette <- sequential_hcl(41, power = 2.2)
ff1 <- levelplot(likeval, xlab = expression(theta[1]), ylab = expression(theta[2]), 
  row.values = thetavector, column.values = thetavector, xlim = range, ylim = range, 
  contour = TRUE, at = c(seq(770, 800, by = 5), seq(805, 920, by = 10)), col.regions = mypalette)
print(ff1)
dev.off()


# Figure 1(d)
theta_lb <- 10^(-15)
theta_ub <- 1
resolution <- 20
range <- c(theta_lb, theta_ub)

thetavector <- seq(from = theta_lb, to = theta_ub, length.out = resolution)
thetavec <- expand.grid(x = thetavector, y = thetavector)
betavec <- as.matrix(log10(thetavec))

likeval <- matrix(0, resolution^2, 1)
# for(i in 1:(resolution^2)){likeval[i] = GP_deviance(betavec[i,],x,y)} # for version 0.1-0
for (i in 1:(resolution^2)) {
  likeval[i] <- GP_deviance(betavec[i, ], x, y, corr = list(type = "exponential", 
    power = 2))
}
dim(likeval) <- c(length(thetavector), length(thetavector))

pdf("fig1d.pdf", height = 5, width = 5)
mypalette <- sequential_hcl(41, power = 2.2)
ff2 <- levelplot(likeval, xlab = expression(theta[1]), ylab = expression(theta[2]), 
  row.values = thetavector, column.values = thetavector, xlim = range, ylim = range, 
  contour = TRUE, cuts = 14, col.regions = mypalette)
print(ff2)
dev.off()

#=================================
# Likelihood plots - Fig 2
#=================================

library("GPfit")
#-- 1d test function -------------
computer_simulator_1 <- function(x) {
  y <- log(x + 0.1) + sin(5 * pi * x)
  return(y)
}
#---------------------
n <- 10
d <- 1
set.seed(12)

# x = maximinLHS(n,d)
x <- c(0.6381151, 0.30458841, 0.51526238, 0.17747182, 0.83945864, 0.73615131, 0.242109, 
  0.43220115, 0.95521608, 0.09712655)
y <- computer_simulator_1(x)

# Figure 2(a)
beta_lb <- -11
beta_ub <- 5
resolution <- 100

betavector <- seq(from = beta_lb, to = beta_ub, length.out = resolution)
betavector <- matrix(betavector, resolution, 1)

likeval <- matrix(0, resolution, 1)
# for(i in 1:(resolution)){likeval[i] = GP_deviance(betavector[i],x,y)} # for version 0.1-0
for (i in 1:resolution) {
  likeval[i] <- GP_deviance(betavector[i], x, y, corr = list(type = "exponential", 
    power = 2))
}

pdf("fig2a.pdf", height = 5, width = 5)
plot(betavector, likeval, type = "l", xlab = expression(beta), ylab = expression(-2 * 
  log(L[beta])))
dev.off()

# Figure 2(b)
beta_lb <- -2 - log10(d)
beta_ub <- 2.7 - log10(d)
resolution <- 100

betavector <- seq(from = beta_lb, to = beta_ub, length.out = resolution)
betavector <- matrix(betavector, resolution, 1)

likeval <- matrix(0, resolution, 1)
# for(i in 1:(resolution)){likeval[i] = GP_deviance(betavector[i],x,y)} # for version 0.1-0
for (i in 1:resolution) {
  likeval[i] <- GP_deviance(betavector[i], x, y, corr = list(type = "exponential", 
    power = 2))
}

pdf("fig2b.pdf", height = 5, width = 5)
plot(betavector, likeval, type = "l", xlab = expression(beta), ylab = expression(-2 * 
  log(L[beta])))
dev.off()

#-- 2d test function -------------
computer_simulator_2 <- function(x) {
  x <- 4 * x - 2
  x1 <- x[, 1]
  x2 <- x[, 2]
  
  t1 <- (1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 6 * x1 * x2 + 
    3 * x2^2))
  t2 <- (30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 36 * 
    x1 * x2 + 27 * x2^2))
  y <- t1 * t2
  
  return(y)
}
#--------------------

library("lattice", "lhs", "colorspace")
#---------------------
n <- 30
d <- 2
set.seed(1)
x <- matrix(0, n, 2)

# x = maximinLHS(n,d)
x[, 1] <- c(0.7937935261, 0.4496580769, 0.2752867432, 0.6913658896, 0.8777401649, 
  0.3987284473, 0.7497571945, 0.1601039682, 0.0527092526, 0.5000438219, 0.6361983972, 
  0.1778581875, 0.3197865148, 0.4003190597, 0.7277299663, 0.0757593439, 0.8408587278, 
  0.4686633414, 0.8024329295, 0.3507487091, 0.9231584449, 0.0005985916, 0.5385867361, 
  0.2234647278, 0.6108420338, 0.5994940089, 0.9631190429, 0.1241414898, 0.9833982868, 
  0.2536672005)
x[, 2] <- c(0.17412559, 0.33371134, 0.78181996, 0.56635533, 0.66482638, 0.16330925, 
  0.09283118, 0.48925892, 0.8443556, 0.6021934, 0.23377032, 0.80397548, 0.86791042, 
  0.43871987, 0.75889476, 0.27293878, 0.93539901, 0.30493853, 0.20142953, 0.72630772, 
  0.4022192, 0.58140528, 0.123962, 0.06284088, 0.53236272, 0.90130616, 0.02741331, 
  0.97629734, 0.6809449, 0.39753711)

y <- computer_simulator_2(x)


# Figure 2(c)
beta_lb <- -5
beta_ub <- 5
resolution <- 40
range <- c(beta_lb, beta_ub)

betavector <- seq(from = beta_lb, to = beta_ub, length.out = resolution)
betavec <- expand.grid(x = betavector, y = betavector)
betavec <- as.matrix(betavec)

likeval <- matrix(0, resolution^2, 1)
# for(i in 1:(resolution^2)){likeval[i] = GP_deviance(betavec[i,],x,y)} # for version 0.1-0
for (i in 1:(resolution^2)) {
  likeval[i] <- GP_deviance(betavec[i, ], x, y, corr = list(type = "exponential", 
    power = 2))
}
dim(likeval) <- c(length(betavector), length(betavector))

pdf("fig2c.pdf", height = 5, width = 5)
mypalette <- sequential_hcl(41, power = 2.2)
ff1 <- levelplot(likeval, xlab = expression(beta[1]), ylab = expression(beta[2]), 
  row.values = betavector, column.values = betavector, xlim = range, ylim = range, 
  contour = TRUE, at = c(seq(750, 990, by = 15)), col.regions = mypalette)
print(ff1)
dev.off()


# Figure 2(d)
beta_lb <- -2 - log10(d)
beta_ub <- 2.7 - log10(d)
resolution <- 30
range <- c(beta_lb, beta_ub)

betavector <- seq(from = beta_lb, to = beta_ub, length.out = resolution)
betavec <- expand.grid(x = betavector, y = betavector)
betavec <- as.matrix(betavec)

likeval <- matrix(0, resolution^2, 1)
# for(i in 1:(resolution^2)){likeval[i] = GP_deviance(betavec[i,],x,y)} # for version 0.1-0
for (i in 1:(resolution^2)) {
  likeval[i] <- GP_deviance(betavec[i, ], x, y, corr = list(type = "exponential", 
    power = 2))
}
dim(likeval) <- c(length(betavector), length(betavector))

pdf("fig2d.pdf", height = 5, width = 5)
mypalette <- sequential_hcl(41, power = 2.2)
ff2 <- levelplot(likeval, xlab = expression(beta[1]), ylab = expression(beta[2]), 
  row.values = betavector, column.values = betavector, xlim = range, ylim = range, 
  contour = TRUE, at = c(seq(750, 990, by = 15)), col.regions = mypalette)
print(ff2)
dev.off()




#==========================
# Example 1
#==========================

computer_simulator <- function(x) {
  y <- log(x + 0.1) + sin(5 * pi * x)
  return(y)
}
#---------------------

library("GPfit", "lhs")

n <- 7
d <- 1
set.seed(2)
x <- maximinLHS(n, 1)

# It is possible that your random number generator may lead different 'x'.  The
# design (x) generated on my computer is as follows
x <- matrix(c(0.1195556, 0.4500716, 0.763896, 0.355539, 0.8784638, 0.6224375, 0.2803777), 
  ncol = 1)


y <- matrix(0, n, 1)
for (i in 1:n) {
  y[i] <- computer_simulator(x[i])
}
# GPmodel <- GP_fit(x, y) # for version 0.1-0
GPmodel <- GP_fit(x, y, corr = list(type = "exponential", power = 2))  # for version 0.2-0

print.GP(GPmodel, digits = 4)

#-------------------------------

resolution <- 100
xpred <- seq(0, 1, length.out = resolution)
ytrue <- matrix(0, 100, 1)
for (i in 1:resolution) {
  ytrue[i] <- computer_simulator(xpred[i])
}


# Figure 3
pdf("eg1_0.pdf", height = 6, width = 6)
# plot.GP(GPmodel, resolution = 100) # for version 0.1-0
plot.GP(GPmodel, resolution = 100, line_type = c(1, 1), legends = TRUE, cex = 2, 
  pch = 1)  # for version 0.2-0
lines(xpred, ytrue, col = "black", lty = 4, lwd = 2)
dev.off()

# Figure 4(a)
pdf("eg1_1.pdf", height = 6, width = 6)
# plot.GP(GPmodel, resolution = 100) # for version 0.1-0
plot.GP(GPmodel, resolution = 100, line_type = c(1, 1), legends = TRUE, cex = 2, 
  pch = 1)  # for version 0.2-0
dev.off()

# Figure 4(b)
pdf("eg1_2.pdf", height = 6, width = 6)
# plot.GP(GPmodel, resolution = 100, line_type = c(1, 2)) # for version 0.1-0
plot.GP(GPmodel, resolution = 100, legends = TRUE, cex = 2, pch = 1)  # for version 0.2-0
dev.off()

# Figure 4(c)
pdf("eg1_3.pdf", height = 6, width = 6)
# plot.GP(GPmodel, resolution = 100, cex = 3) # for version 0.1-0
plot.GP(GPmodel, resolution = 100, line_type = c(1, 1), legends = TRUE, cex = 3, 
  pch = 1)  # for version 0.2-0
dev.off()

# Figure 4(d)
pdf("eg1_4.pdf", height = 6, width = 6)
# plot.GP(GPmodel, resolution = 100, line_type = c(1, 2), pch = 2, cex = 3)# for version 0.1-0
plot.GP(GPmodel, resolution = 100, legends = TRUE, pch = 2, cex = 3)  # for version 0.2-0
dev.off()

#==========================
# Example 2
#==========================
rm(list = ls())

#-----------------------
computer_simulator <- function(x) {
  x <- 4 * x - 2
  x1 <- x[1]
  x2 <- x[2]
  
  t1 <- (1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 6 * x1 * x2 + 
    3 * x2^2))
  t2 <- (30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 36 * 
    x1 * x2 + 27 * x2^2))
  y <- t1 * t2
  
  return(y)
}
#--------------------
library("GPfit")
library("lhs")
#---------------------
n <- 20
d <- 2
#---------------------
set.seed(3)
x <- maximinLHS(n, d)

# It is possible that your random number generator may lead different 'x'.  The
# design (x) generated on my computer is as follows
x[, 1] <- c(0.01255191, 0.71814188, 0.38710637, 0.48533217, 0.13898664, 0.87006243, 
  0.94981408, 0.66655037, 0.40426038, 0.64272338, 0.07615149, 0.54985172, 0.97109437, 
  0.81143146, 0.75517806, 0.29494404, 0.31269219, 0.57066778, 0.20456861, 0.19059827)
x[, 2] <- c(0.03046171, 0.52727509, 0.89305787, 0.26445776, 0.10251103, 0.36296337, 
  0.71618854, 0.67839353, 0.79640079, 0.18882006, 0.60236048, 0.44754285, 0.99703747, 
  0.23115798, 0.4533835, 0.07648327, 0.55592569, 0.9291943, 0.31803503, 0.8159798)

y <- matrix(0, n, 1)
for (i in 1:n) {
  y[i] <- computer_simulator(x[i, ])
}
# GPmodel <- GP_fit(x, y) # for version 0.1-0
GPmodel <- GP_fit(x, y, corr = list(type = "exponential", power = 2))  # for version 0.2-0

#--------------------------
print.GP(GPmodel, digits = 4)

#-------------------------------
set.seed(5)
xnew <- matrix(runif(20), ncol = 2)
Model_pred <- predict.GP(GPmodel, xnew, corr = list(type = "exponential", power = 2))
Model_pred
#--------------------------

library("colorspace")

# Figure 5(a) --------------
pdf("eg2_1.pdf", height = 6, width = 6)
mypalette <- sequential_hcl(51, power = 2.2)
plot.GP(GPmodel, col.regions = mypalette, cuts = 50)
dev.off()
#--------------

# Figure 5(b) ---------
pdf("eg2_2.pdf", height = 6, width = 6)
mypalette <- sequential_hcl(51, power = 2.2)
plot.GP(GPmodel, response = FALSE, contour = TRUE, col.regions = mypalette, cuts = 20)
dev.off()
#---------

# Figure 5(c) ---------
pdf("eg2_3.pdf", height = 6, width = 6)
plot.GP(GPmodel, surf_check = TRUE)
dev.off()
#---------

# Figure 5(d) ---------
pdf("eg2_4.pdf", height = 6, width = 6)
plot.GP(GPmodel, response = FALSE, surf_check = TRUE)
dev.off()
#---------

#=============================================
#=============================================
library("GPfit", "lhs", "lattice", "mlegp", "colorspace")

# Real data (tidal power)
data <- read.table("tidalpower_output.txt", head = FALSE)
x1vec <- (data[, 1] - min(data[, 1]))/(max(data[, 1]) - min(data[, 1]))
x2vec <- (data[, 2] - min(data[, 2]))/(max(data[, 2]) - min(data[, 2]))
yvec <- data[, 3]
x1mat <- matrix(x1vec, nrow = 13)
x2mat <- matrix(x2vec, nrow = 13)
ymat <- matrix(yvec, nrow = 13)
Xvec <- cbind(x1vec, x2vec)

#--------------------------------------------------
# Choosing a 30-point space-filling design
n <- 30
d <- 2
set.seed(1)
bestid <- c(141, 470, 241, 394, 17, 489, 498, 285, 315, 49, 238, 54, 318, 235, 505, 
  231, 397, 474, 365, 376, 481, 109, 321, 65, 138, 145, 46, 182, 413, 174)

Xdes <- Xvec[bestid, ]
y_Xdes <- yvec[bestid]
#--------------------------------------------------
# fitting GP model
GPmodel <- GP_fit(Xdes, y_Xdes, corr = list(type = "exponential", power = 2))

GPprediction <- predict.GP(GPmodel, Xvec, 1, corr = list(type = "exponential", power = 2))
Y_hat1 <- GPprediction$Y_hat
MSE1 <- GPprediction$MSE

sRMSE1 <- sqrt(mean((Y_hat1[-bestid] - yvec[-bestid])^2))/(max(yvec[-bestid]) - min(yvec[-bestid]))
print(sRMSE1)


#-----------
# fitting mlegp model
mlegpmodel <- mlegp(Xdes, y_Xdes)

mlegp_prediction <- predict(mlegpmodel, newData = Xvec, se.fit = TRUE)
Y_hat2 <- mlegp_prediction$fit
MSE2 <- mlegp_prediction$se.fit

sRMSE2 <- sqrt(mean((Y_hat2[-bestid] - yvec[-bestid])^2))/(max(yvec[-bestid]) - min(yvec[-bestid]))
print(sRMSE2)



#--------------------------------------------------
range <- c(0, 1)
resolution <- 30

xvector <- seq(from = range[1], to = range[2], length.out = resolution)
xvec <- expand.grid(x = xvector, y = xvector)
xvec <- as.matrix(xvec)

mypalette <- sequential_hcl(51, power = 2.2)

# Figure 7(a) ---------------------
GPprediction <- predict.GP(GPmodel, xvec, 1, corr = list(type = "exponential", power = 2))
Y_hat1 <- GPprediction$Y_hat
MSE1 <- GPprediction$MSE
dim(Y_hat1) <- c(length(xvector), length(xvector))

pdf("fig7a.pdf", height = 6, width = 6)
h1 <- levelplot(Y_hat1, xlab = expression(X[1]), ylab = expression(X[2]), row.values = xvector, 
  column.values = xvector, xlim = range, ylim = range, contour = TRUE, col.regions = mypalette)
print(h1)
dev.off()

# Figure 7(b) ---------------------
mlegp_prediction <- predict(mlegpmodel, newData = xvec, se.fit = TRUE)
Y_hat2 <- mlegp_prediction$fit
MSE2 <- mlegp_prediction$se.fit
dim(Y_hat2) <- c(length(xvector), length(xvector))

pdf("fig7b.pdf", height = 6, width = 6)
h2 <- levelplot(Y_hat2, xlab = expression(X[1]), ylab = expression(X[2]), row.values = xvector, 
  column.values = xvector, xlim = range, ylim = range, contour = TRUE, col.regions = mypalette)
print(h2)
dev.off()
#---------------------

#=============================================
#--------------------------------------------------
# Choosing a 50-point space-filling design
n <- 50
d <- 2
set.seed(1)
bestid <- c(62, 69, 238, 496, 201, 512, 515, 377, 449, 48, 18, 67, 451, 301, 491, 
  129, 267, 459, 144, 468, 519, 97, 193, 78, 139, 366, 33, 195, 310, 197, 308, 
  160, 248, 15, 427, 189, 408, 121, 290, 125, 431, 331, 385, 311, 247, 374, 25, 
  244, 367, 230)

Xdes <- Xvec[bestid, ]
y_Xdes <- yvec[bestid]
#--------------------------------------------------
# fitting GP model
GPmodel <- GP_fit(Xdes, y_Xdes, corr = list(type = "exponential", power = 2))

GPprediction <- predict.GP(GPmodel, Xvec, 1, corr = list(type = "exponential", power = 2))
Y_hat1 <- GPprediction$Y_hat
MSE1 <- GPprediction$MSE

sRMSE1 <- sqrt(mean((Y_hat1[-bestid] - yvec[-bestid])^2))/(max(yvec[-bestid]) - min(yvec[-bestid]))
print(sRMSE1)


#-----------
# fitting mlegp model
mlegpmodel <- mlegp(Xdes, y_Xdes)

mlegp_prediction <- predict(mlegpmodel, newData = Xvec, se.fit = TRUE)
Y_hat2 <- mlegp_prediction$fit
MSE2 <- mlegp_prediction$se.fit

sRMSE2 <- sqrt(mean((Y_hat2[-bestid] - yvec[-bestid])^2))/(max(yvec[-bestid]) - min(yvec[-bestid]))
print(sRMSE2)


#--------------------------------------------------
range <- c(0, 1)
resolution <- 30

xvector <- seq(from = range[1], to = range[2], length.out = resolution)
xvec <- expand.grid(x = xvector, y = xvector)
xvec <- as.matrix(xvec)

mypalette <- sequential_hcl(51, power = 2.2)

# Figure 8(a) ---------------------
GPprediction <- predict.GP(GPmodel, xvec, 1, corr = list(type = "exponential", power = 2))
Y_hat1 <- GPprediction$Y_hat
MSE1 <- GPprediction$MSE
dim(Y_hat1) <- c(length(xvector), length(xvector))

pdf("fig8a.pdf", height = 6, width = 6)
h1 <- levelplot(Y_hat1, xlab = expression(X[1]), ylab = expression(X[2]), row.values = xvector, 
  column.values = xvector, xlim = range, ylim = range, contour = TRUE, col.regions = mypalette)
print(h1)
dev.off()

# Figure 8(b) ---------------------
mlegp_prediction <- predict(mlegpmodel, newData = xvec, se.fit = TRUE)
Y_hat2 <- mlegp_prediction$fit
MSE2 <- mlegp_prediction$se.fit
dim(Y_hat2) <- c(length(xvector), length(xvector))

pdf("fig8b.pdf", height = 6, width = 6)
h2 <- levelplot(Y_hat2, xlab = expression(X[1]), ylab = expression(X[2]), row.values = xvector, 
  column.values = xvector, xlim = range, ylim = range, contour = TRUE, col.regions = mypalette)
print(h2)
dev.off()
#---------------------

#=============================================
