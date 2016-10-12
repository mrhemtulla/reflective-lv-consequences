#do some simulations! 

rm(list = ls())

setwd("//tsclient/C/Users/Mijke/Documents/__STUDIES__/CausalIndicators/MediationModelSimulation")
setwd("C:/Users/Mijke/Documents/__STUDIES__/CausalIndicators/MediationModelSimulation")
setwd("C:/Users/mijke/Box Sync/MediationModelSimulation")
source("Functions.R")
windowsFonts(A=windowsFont("Verdana"))

#### vary covariance between m indicators ####
results <- matrix(0, 96, 7) 
row <- 1
for (i in seq(0, .95, by = .01)) {
  data <- generateCov(a = NA, b = .4, avec = c(.5, .3, .2, 0), covm = i, res.M = 0)
  fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000)  
  fitLV <- sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  
  results [row, 1] <- i
  results [row, 2] <- inspect(fitComp, "coef")$beta[1,3] #a path composite model
  results [row, 3] <- inspect(fitComp, "coef")$beta[2,1] #b path composite model  
  results [row, 4] <- inspect(fitLV, "coef")$beta[1,3] #a path LV model
  results [row, 5] <- inspect(fitLV, "coef")$beta[2,1] #b path LV model
  results [row, 6] <- try(inspect(fitComp, "fit")["fmin"], silent = TRUE) #fit of Composite model
  LV.converged <- inspect(fitLV, "converged")
  if (LV.converged == TRUE) {results [row, 7] <- inspect(fitLV, "fit")["fmin"]} else {results [row, 7] <- NA}  
  row <- row + 1  
}

results <- as.data.frame(results[(1:(row-1)),], stringsAsFactors = FALSE)
colnames(results) <- c("covm", "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", "Comp.fmin", "LV.fmin")
head(results)
results <- results[-(1:12),]

#setEPS()
#postscript("Sim1.eps", width=9, height=4, family = "Palatino", paper = "special")
#pdf("Sim1.pdf", width=9, height=4)  
jpeg("Sim1.jpg", width = 9, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, .5, 0, 1))
  plot(results$covm, results$CompEst.a, type = "l", xlab = "covariance between causal indicators",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1))
  lines(results$covm, results$LVEst.a, col = "darkred")
  
  par(fig = c(.5, 1, 0, 1), new = TRUE)
  plot(results$covm, results$CompEst.b, type = "l", xlab = "covariance between causal indicators",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1))
  lines(results$covm, results$LVEst.b, col = "darkred")
  
  legend(x = .2, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
dev.off()

results$RMSEA <- sqrt(results$LV.fmin/9)

jpeg("Sim1RMSEA.jpg", width = 4, height = 2.5, units = "in", res = 800)
par(mar = c(4, 4, 1, .5), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results$covm, results$RMSEA, type = "l", xlab = "covariance between causal indicators",
     ylab = "pop RMSEA", col = "darkred", ylim = c(0, .4))
dev.off()

#### vary the amount of difference in the effect of X on each m systematically ####

results <- matrix(0, 21, 7) 
row <- 1

for (i in seq(0, .2, by = .01)) {  
  a <- .4
  b1 <- .3 + i 
  b2 <- .3 + 2*i 
  b3 <- .3 - i
  b4 <- .3 - 2*i
  data <- generateCov(a = a, b = .4, b1, b2, b3, b4, res.M = 0)
  fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000)  
  fitLV <- sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  
  results [row, 1] <- i
  results [row, 2] <- (((b1 + b2 + b3 + b4)^2)/(12*a^2)) - 1/3  #covm 
  results [row, 3] <- sqrt(1/(4+12*covm)) #set causal indicator 'loading' s.t. M has total variance 1 
  results [row, 4] <- inspect(fitComp, "coef")$beta[1,3] #a path composite model
  results [row, 5] <- inspect(fitComp, "coef")$beta[2,1] #b path composite model  
  results [row, 6] <- inspect(fitLV, "coef")$beta[1,3] #a path LV model
  results [row, 7] <- inspect(fitLV, "coef")$beta[2,1] #b path LV model
  #results [row, 6] <- try(inspect(fitComp, "fit")["fmin"], silent = TRUE) #fit of LV model
  #results [row, 7] <- try(inspect(fitLV, "fit")["fmin"], silent = TRUE) #fit of LV model (df = 9)
  row <- row + 1  
}

results <- as.data.frame(results[(1:(row-1)),], stringsAsFactors = FALSE)
colnames(results) <- c("spread", "covm", "lambda", "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b")
head(results)


jpeg("Sim2.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, .5, 0, 1))
plot(results$spread*4, results$CompEst.a, type = "l", xlab = "range of X -> m path values",
     ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1))
lines(results$spread*4, results$LVEst.a, col = "darkred")
#lines(results$spread*4, results$covm, col = "darkorange") #does not vary
#lines(results$spread*4, results$lambda, col = "darkorange") #does not vary
#hm. what is the scale of these things? are they comparable? does M have variance of 1 in both models?

par(fig = c(.5, 1, 0, 1), new = TRUE)
plot(results$spread*4, results$CompEst.b, type = "l", xlab = "range of in X -> m path values",
     ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1))
lines(results$spread*4, results$LVEst.b, col = "darkred") #erm, huh? 

legend(x = .2, y = 1, c("composite model", "reflective model"), 
       fill = c("darkgreen", "darkred"), bty = "n")
dev.off()

#### vary the amount of difference in the effect of X on each m randomly ####
#NEW: rather than getting one line, I'm now generating random values of a1-a4, holding the total 
# a path constant at .4. Thus the models vary on 2 dimensions: variance in a paths, and the 
# covariance between the causal indicators (because this is determined s.t. total a = .4)

results <- matrix(0, 10000, 11) 
row <- 1

#for (i in 1:10000) { #seq(0, .04, by = .005)) {  
#note that sim2 doesn't filter out non-PD matrices
sim2 <- function(i){
  a <- .4
  b <- .4
  avec <- runif(4, 0, .5)
  #avec <- c(meana + i, meana + 2*i, meana - i, meana - 2*i)
  meana <- mean(avec)
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = 0)
  
  amat <- avec %*% t(avec)
  gamma <- sqrt(a/sum(amat))   
  covm <- (1/(12*gamma^2)) - 1/3
  
  fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000)  
  fitLV <- sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  
  #fitFull <- sem(model = FullMod, sample.cov = data, sample.nobs = 100000, sample.cov.rescale = FALSE)

  conv.LV <- try(inspect(fitLV, "converged"))
  SEs.LV <- TRUE  
  if (class(try(vcov(fitLV))) == "try-error") {SEs.LV <- FALSE}
  if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  

  return(c(meana, rangea, vara, covm, gamma,
           inspect(fitLV, "coef")$lambda[1:4,1],
           inspect(fitComp, "coef")$beta[1,3], #a path composite model
           inspect(fitComp, "coef")$beta[2,1], #b path composite model  
           inspect(fitLV, "coef")$beta[1,3], #a path LV model
           inspect(fitLV, "coef")$beta[2,1], #b path LV model
           inspect(fitComp, "fit")["fmin"], #fit of composite model
           conv.LV, SEs.LV, fit.LV))
}

X <- 1:5000
ptm <- proc.time()
out <- sapply(X, sim2)
proc.time() - ptm

save(out, file = "sim2.5000.Rdata")
#load(file = "sim2.5000.Rdata")
results <- as.data.frame(t(out))
colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", 
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "Comp.fmin", "conv.LV", "SEconv.LV", "LV.fmin")
round(results[1:25,c(1,2,4:9,12:13,16:17)], 3)

#discard rows that did not converge: 
results2 <- results[results$conv.LV == 1,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #3935 rows left

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))

results2$col.var.a <- rbPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col.covm <- ryPal(10)[as.numeric(cut(results2$covm,breaks = 10))]

# Sim2.Plot1.jpg: estimates by range in a1-a4
jpeg("Sim2.Plot1.jpg", width = 9, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, .5, 0, 1))
  plot(results2$range.a, results2$CompEst.a, type = "p", xlab = "range of X -> m path values",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = '.', cex = 2)
  points(results2$range.a, results2$LVEst.a, col = results2$col.covm, pch = '.', cex = 2)
  
  par(fig = c(.5, 1, 0, 1), new = TRUE)
  plot(results2$range.a, results2$CompEst.b, type = "p", xlab = "range of in X -> m path values",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = '.', cex = 2)
  points(results2$range.a, results2$LVEst.b, col = results2$col.covm, pch = '.', cex = 2) 
  
  legend(x = .2, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
dev.off()

# Sim2.Plot2.jpg" estimates by variance in a1-a4 
jpeg("Sim2.Plot2.jpg", width = 9, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, .5, 0, 1))
  plot(results2$var.a, results2$CompEst.a, type = "p", xlab = "variance of X -> m path values",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2)
  points(results2$var.a, results2$LVEst.a, col = results2$col.covm, pch = ".", cex = 2)
  
  par(fig = c(.5, 1, 0, 1), new = TRUE)
  plot(results2$var.a, results2$CompEst.b, type = "p", xlab = "variance of in X -> m path values",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2)
  points(results2$var.a, results2$LVEst.b, col = results2$col.covm, pch = ".", cex = 2) 
  
  legend(x = .2, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
dev.off()

# Sim2.Plot3.jpg : estimates by covariance among indicators (covm) 
jpeg("Sim2.Plot3.jpg", width = 9, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, .5, 0, 1))
  plot(results2$covm, results2$CompEst.a, type = "p", xlab = "covariance between causal indicators",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2)
  points(results2$covm, results2$LVEst.a, col = results2$col.var.a, pch = ".", cex = 2)
#abline(v = -.24) # what is this number and why? 
#abline(v = .2)

  par(fig = c(.5, 1, 0, 1), new = TRUE)
  plot(results2$covm, results2$CompEst.b, type = "p", xlab = "covariance between causal indicators",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2)
  points(results2$covm, results2$LVEst.b, col = results2$col.var.a, pch = ".", cex = 2) 
  legend(x = .2, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
#  abline(v = .2)
dev.off()

# Sim2.Plot4.jpg : Fit by covariance among indicators 
jpeg("Sim2.Plot4.Fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$LV.fmin, type = "p", xlab = "covariance between causal indicators",
     ylab = "fmin of LV model", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2)
points(results2$covm, results2$LV.fmin, col = results2$col.var.a, pch = ".", cex = 2)
abline(v = -.235) # what is this number and why? 
abline(v = .2)
dev.off()


#try to figure some stuff out#### 
LVunderest <- results2[results2$LVEst.a < .4,]
LVunderest[,c(4,8,13)]
col2rgb("#CD0000") #red dots at the very bottom -- covm ~ -.3
col2rgb("#F3C600") #yellow -- covm ~ .25
col2rgb("#F9E200") #bright yellow -- covm ~ .30

#line of red dots where estimates of a, b = 0: 
round(results2[results2$col.covm == "#CD0000",c(1:5,8,9,11)], 3)
  #gamma is >1 in all cases
  #(should record factor loadings in LV model) (added that to loop but have to re-run still)

#what happens when covm < .21 -- some weird things happening
##this all goes away when runs that failed to produce SEs are removed. 

temp <- round(results2[results2$covm < -.21, c(1:2,4:9,12:13,15:17)], 3)
dim(temp) #421 remaining cases

head(temp)

#linear function: square of mean a path perfectly correlated with covm: 
cor(results2$mean.a^2, results2$covm) #1. 
summary(lm(covm ~ I(mean.a^2), data = results2)) #covm = -1/3 + 10/3(mean.a^2) 
#so when mean.a^2 = .1, covm = 0

plot(results2$mean.a^2, results2$covm, type = "p", col = results2$col.var.a, pch = ".")
abline(h = 0)
abline(v = .1)
#so if we want covm to be positive, square of mean a path must be higher than .10
sqrt(.1) #.3162278

#plot range of a1-a4 values by covm 
cor(results2$covm, results2$range.a) #-.287 
cor(results2$mean.a, results2$range.a) #-.287 
jpeg("Sim2.Plot5.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
# plot(results2$covm, results2$range.a, pch = ".", col = "darkgreen", 
#      xlab = "covariance among causal indicators", ylab = "range of X -> m path values")
plot(results2$mean.a, results2$range.a, pch = ".", col = "darkgreen", 
     xlab = "mean of X -> m path values", ylab = "range of X -> m path values")
dev.off()

######################## SIM 3 ###################################################

sim3 <- function(i){
  a <- .6 #up this value!
  b <- .4
  res.M <- 0
  avec <- runif(4, 0, .8) #in sim2 this had max = .5
  data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = res.M)  
  while (min(eigen(data[-c(5,8),-c(5,8)])$values) < 1e-10) { #ensure that sigma is positive definite! 
    avec <- runif(4, 0, .8) #in sim2 this had max = .5
    data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = res.M)
  }
  meana <- mean(avec)
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  amat <- avec %*% t(avec)
  gamma <- a/sum(avec)
  covm <- (1 - 4*gamma^2 - res.M)/(12*gamma^2)   #NEW so total var(M) includes residual variance
  
  fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000)  
  fitLV <- sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  
  #fitFull <- sem(model = FullMod, sample.cov = data, sample.nobs = 100000, sample.cov.rescale = FALSE)
  
  conv.LV <- try(inspect(fitLV, "converged"))
  SEs.LV <- TRUE  
  if (class(try(vcov(fitLV))) == "try-error") {SEs.LV <- FALSE}
  if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  
  
  return(c(meana, rangea, vara, covm, gamma,
           inspect(fitLV, "coef")$lambda[1:4,1],
           inspect(fitComp, "coef")$beta[1,3], #a path composite model
           inspect(fitComp, "coef")$beta[2,1], #b path composite model  
           inspect(fitLV, "coef")$beta[1,3], #a path LV model
           inspect(fitLV, "coef")$beta[2,1], #b path LV model
           inspect(fitLV, "std")$lambda[1:4,1],
           inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
           inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
           inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
           inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
           inspect(fitComp, "fit")["fmin"], #fit of composite model
           conv.LV, SEs.LV, fit.LV,
           inspect(fitComp, "se")$beta[1,3], #SE a composite model
           inspect(fitComp, "se")$beta[2,1], #SE b composite model
           inspect(fitLV, "se")$beta[1,3], #SE a LV model
           inspect(fitLV, "se")$beta[2,1])) #SE b LV model
}

X <- 1:100000
ptm <- proc.time()
out <- sapply(X, sim3)
proc.time() - ptm

save(out, file = "sim3.100000.Rdata")
#load(file = "sim3.10000.Rdata")
results <- as.data.frame(t(out))

colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", 
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "std.lam1", "std.lam2", "std.lam3", "std.lam4",
                       "std.CompEst.a", "std.CompEst.b", "std.LVEst.a", "std.LVEst.b", 
                       "Comp.fmin", "conv.LV", "SEconv.LV", "LV.fmin",
                       "CompSE.a", "CompSE.b", "LVSE.a", "LVSE.b")

#discard rows that did not converge or that generated non-positive-definite sigma: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #4758 rows left

## fixed results to contain correct gamma and covm. Resave: 
#save(Sim3Results, file = "sim3.10000.Rdata")
#load(file = "sim3.10000.Rdata")
#results2 <- Sim3Results

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))

results2$col.var.a <- rbPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col.covm <- ryPal(10)[as.numeric(cut(results2$covm,breaks = 10))]

jpeg("sim3.Plot1.jpg", width = 12, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, 0, 1))
  plot(results2$var.a, results2$std.CompEst.a, type = "p", xlab = "variance of X --> m paths",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$var.a, results2$std.LVEst.a, col = results2$col.covm, pch = ".", cex = 1)
  
  par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
  plot(results2$var.a, results2$std.CompEst.b, type = "p", xlab = "variance of X --> m paths",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$var.a, results2$std.LVEst.b, col = results2$col.covm, pch = ".", cex = 1) 
  legend(x = 0, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
  
  #plot ab path estimate
  results2$std.LVEst.ab <- results2$std.LVEst.a * results2$std.LVEst.b
  par(fig = c(2/3, 1, 0, 1), new = TRUE)
  plot(results2$var.a, results2$std.CompEst.a*results2$std.CompEst.b, type = "p", 
       xlab = "variance of X --> m paths",
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$var.a, results2$std.LVEst.ab, col = results2$col.covm, pch = ".", cex = 1)
dev.off()



jpeg("sim3.Plot3.jpg", width = 12, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, 0, 1))
  plot(results2$covm, results2$std.CompEst.a, type = "p", xlab = "correlation between causal indicators",
       ylab = "estimate of a path", col = "darkgreen", xlim = c(0, 1), ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$covm, results2$std.LVEst.a, col = results2$col.var.a, pch = ".", cex = 1)
  
  par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
  plot(results2$covm, results2$std.CompEst.b, type = "p", xlab = "correlation between causal indicators",
       ylab = "estimate of b path", col = "darkgreen", xlim = c(0, 1), ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$covm, results2$std.LVEst.b, col = results2$col.var.a, pch = ".", cex = 1) 
  legend(x = 0, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
  
  #plot ab path estimate
  results2$std.LVEst.ab <- results2$std.LVEst.a * results2$std.LVEst.b
  par(fig = c(2/3, 1, 0, 1), new = TRUE)
  plot(results2$covm, results2$std.CompEst.a*results2$std.CompEst.b, type = "p", 
       xlab = "correlation between causal indicators", xlim = c(0, 1),
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$covm, results2$std.LVEst.ab, col = results2$col.var.a, pch = ".", cex = 1)
dev.off()

#plot by mean a path instead of covm
jpeg("sim3.Plot3b.jpg", width = 12, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, 0, 1))
  plot(results2$mean.a, results2$CompEst.a, type = "p", xlab = "mean a path",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 6), pch = ".", cex = 1)
  points(results2$mean.a, results2$LVEst.a, col = results2$col.var.a, pch = ".", cex = 1)
  
  par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
  plot(results2$mean.a, results2$CompEst.b, type = "p", xlab = "mean a path",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$mean.a, results2$LVEst.b, col = results2$col.var.a, pch = ".", cex = 1) 
  legend(x = 0, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
  
  #plot ab path estimate
  results2$LVEst.ab <- results2$LVEst.a * results2$LVEst.b
  #jpeg("sim3.Plot6.jpg", width = 5, height = 4, units = "in", res = 800)
  #par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(2/3, 1, 0, 1), new = TRUE)
  plot(results2$mean.a, results2$CompEst.a*results2$CompEst.b, type = "p", 
       xlab = "mean a path",
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(.24, .34), pch = ".", cex = 1)
  points(results2$mean.a, results2$LVEst.ab, col = results2$col.var.a, pch = ".", cex = 1)
dev.off()


jpeg("sim3.Plot4.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$LV.fmin, type = "p", xlab = "covariance between causal indicators",
     ylab = "min fit function LV model", col = results2$col.var.a, ylim = c(0, 1), 
     pch = ".", cex = 1)
dev.off()

#plot fit in terms of population RMSEA instead of fmin: 
results2$pop.RMSEA <- sqrt(results2$LV.fmin/9) 
jpeg("sim3.Plot4.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$pop.RMSEA, type = "p", xlab = "correlation between causal indicators",
     ylab = "Population RMSEA", col = results2$col.var.a, ylim = c(0, 1), 
     pch = ".", cex = 1)
dev.off()


#plot range of a1-a4 values by covm: 
cor(results2$covm, results2$range.a) #-.287 
cor(results2$mean.a, results2$range.a) #-.287 
jpeg("Sim3.Plot5.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$mean.a, results2$range.a, pch = ".", cex = 1, col = results2$col.var.a, 
     xlab = "covariance among causal indicators", ylab = "range of X -> m path values")
dev.off()


#plot fit as a function of the range of X -> m values, and covm (colour)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
jpeg("sim3.Plot7.jpg", width = 9, height = 4, units = "in", res = 800)
plot(results2$range.a, results2$LV.fmin, type = "p", xlab = "range in a1-a4 paths",
     ylab = "min fit function LV model", col = results2$col.covm, ylim = c(0, 1), 
     pch = ".", cex = 1)
dev.off()

#find example of same range in a1-a4 but different cov: 
head(results2)
?order

results2$range.a2 <- round(results2$range.a, 2)
which(results2$LV.fmin == min(results2[results2$range.a2 == .7,]$LV.fmin)) #2049
which(results2$LV.fmin == max(results2[results2$range.a2 == .7,]$LV.fmin)) #2970

results2[c(2049, 2970),]

cor(results2$LVEst.a, results2$LV.fmin)

#plot fit as a function of the range of X -> m values, and covm (colour)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
jpeg("sim3.Plot8.jpg", width = 9, height = 4, units = "in", res = 800)
plot(results2$LVEst.a*results2$LVEst.b, results2$LV.fmin, pch = ".", 
     col = results2$col.var.a, cex = 1, xlab = "a*b path estimate", ylab = "fmin")
dev.off()


######################## SIM 4 ###################################################

sim4 <- function(i){
  a <- .5 #increase this to increase the range of possible values.  
  b <- .5
  avec <- runif(4, 0, .8) #in sim2 this had max = .5
  meana <- mean(avec)
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = 0)

  fitLV <- try(sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE, warn = FALSE), silent = TRUE)  
  fitLV
  

  if (min(eigen(data[-5,-5])$values) > 0) {
    amat <- avec %*% t(avec)
    gamma <- a/sum(avec)
    covm <- (1/(12*gamma^2)) - 1/3
    
    fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000, warn = FALSE)  
    fitLV <- try(sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE, warn = FALSE), silent = TRUE)  
    
    conv.LV <- try(inspect(fitLV, "converged"), silent = TRUE)
    SEs.LV <- TRUE  
    if (identical(class(try(vcov(fitLV), silent = TRUE)), "try-error") == TRUE) {SEs.LV <- FALSE}
    if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  
    
    return(c(meana, rangea, vara, covm, gamma,
             inspect(fitLV, "coef")$lambda[1:4,1],
             inspect(fitComp, "coef")$beta[1,3], #a path composite model
             inspect(fitComp, "coef")$beta[2,1], #b path composite model  
             inspect(fitLV, "coef")$beta[1,3], #a path LV model
             inspect(fitLV, "coef")$beta[2,1], #b path LV model
             inspect(fitComp, "fit")["fmin"], #fit of composite model
             conv.LV, SEs.LV, fit.LV))    
  } else {return (rep(NA, 17))}
}

X <- 1:10000
ptm <- proc.time()
out <- sapply(X, sim4)
proc.time() - ptm

save(out, file = "sim4.10000.Rdata")
load(file = "sim4.10000.Rdata")
results <- as.data.frame(t(out))

#new to reflect added columns: 
colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", 
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "Comp.fmin", "conv.LV", "SEconv.LV", "LV.fmin")
round(results[1:25,c(1,2,4:9,12:13,16:17)], 3)

#discard rows that did not converge or that generated non-positive-definite sigma: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #5652 rows left

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))

results2$col.var.a <- rbPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col.covm <- ryPal(10)[as.numeric(cut(results2$covm,breaks = 10))]

jpeg("sim4.Plot3.jpg", width = 12, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, 1/3, 0, 1))
plot(results2$covm, results2$CompEst.a, type = "p", xlab = "covariance between causal indicators",
     ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 5), pch = ".", cex = 2)
points(results2$covm, results2$LVEst.a, col = results2$col.var.a, pch = ".", cex = 2)

par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
plot(results2$covm, results2$CompEst.b, type = "p", xlab = "covariance between causal indicators",
     ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2)
points(results2$covm, results2$LVEst.b, col = results2$col.var.a, pch = ".", cex = 2) 
legend(x = -.2, y = 1, c("composite model", "reflective model"), 
       fill = c("darkgreen", "darkred"), bty = "n")

#plot ab path estimate #
results2$LVEst.ab <- results2$LVEst.a * results2$LVEst.b
#jpeg("sim4.Plot6.jpg", width = 5, height = 4, units = "in", res = 800)
par(fig = c(2/3, 1, 0, 1), new = TRUE)
plot(results2$covm, results2$CompEst.a*results2$CompEst.b, type = "p", 
     xlab = "covariance between causal indicators",
     ylab = "estimate of a*b path", col = "darkgreen", ylim = c(.25, .38), pch = ".", cex = 2)
points(results2$covm, results2$LVEst.ab, col = results2$col.var.a, pch = ".", cex = 2)
dev.off()

#plot fit  #
jpeg("sim4.Plot4.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$LV.fmin, type = "p", xlab = "covariance between causal indicators",
     ylab = "min fit function LV model", col = results2$col.var.a, ylim = c(0, 1), 
     pch = ".", cex = 2)


#this calls for a 3D plot... 
library(rgl)
example(plot3d)

#fun and stingray-like plot:  
plot3d(x = results2$covm, y = results2$LVEst.ab, z = results2$range.a, col = results2$col.var.a)

#### SIM 5a ###################################
#choose fixed values of covm, get a range of DIF; vice versa. 

a <- .6
covm <- .4

#covm <- ((1/(12*a^2))*(meana*4)^2) - (1/3)
meana <- sqrt((3*a^2/4)*(covm + (1/3))) #phew

sim5a <- function(i){
  a <- .6 
  b <- .4
  covm <- .4
  meana <- sqrt((3*a^2/4)*(covm + (1/3))) #phew
  a4 <- 1
  while (a4 > .8 || a4 < 0){
    avec <- runif(3, 0, .8) #in sim2 this had max = .5
    a4 <- meana*4 - sum(avec)    
  }
  avec <- c(avec, a4)   
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = 0)
  
  if (min(eigen(data[-5,-5])$values) > 0) { #ensure that sigma is positive definite! 
    #gamma <- sqrt(a/sum(amat))  ###WRONG
    gamma <- a/sum(avec)
    covm <- (1/(12*gamma^2)) - 1/3    
    
    fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000)  
    fitLV <- sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  
    #fitFull <- sem(model = FullMod, sample.cov = data, sample.nobs = 100000, sample.cov.rescale = FALSE)
    
    conv.LV <- try(inspect(fitLV, "converged"))
    SEs.LV <- TRUE  
    if (class(try(vcov(fitLV))) == "try-error") {SEs.LV <- FALSE}
    if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  
    
    return(c(meana, rangea, vara, covm, gamma,
             inspect(fitLV, "coef")$lambda[1:4,1],
             inspect(fitComp, "coef")$beta[1,3], #a path composite model
             inspect(fitComp, "coef")$beta[2,1], #b path composite model  
             inspect(fitLV, "coef")$beta[1,3], #a path LV model
             inspect(fitLV, "coef")$beta[2,1], #b path LV model
             inspect(fitLV, "std")$lambda[1:4,1],
             inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
             inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
             inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
             inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
             inspect(fitComp, "fit")["fmin"], #fit of composite model
             conv.LV, SEs.LV, fit.LV,
             inspect(fitComp, "se")$beta[1,3], #SE a composite model
             inspect(fitComp, "se")$beta[2,1], #SE b composite model
             inspect(fitLV, "se")$beta[1,3], #SE a LV model
             inspect(fitLV, "se")$beta[2,1])) #SE b LV model))    
  } else {return (rep(NA, 29))}
}

X <- 1:100000
ptm <- proc.time()
out5a <- sapply(X, sim5a)
proc.time() - ptm

save(out5a, file = "sim5a.10000.Rdata")
#load(file = "sim5a.10000.Rdata")
results <- as.data.frame(t(out5a))

colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", 
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "std.lam1", "std.lam2", "std.lam3", "std.lam4",
                       "std.CompEst.a", "std.CompEst.b", "std.LVEst.a", "std.LVEst.b", 
                       "Comp.fmin", "conv.LV", "SEconv.LV", "LV.fmin",
                       "CompSE.a", "CompSE.b", "LVSE.a", "LVSE.b")


#discard rows that did not converge or that generated non-positive-definite sigma: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
#dim(results2) #9636 rows left


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))

results2$col.var.a <- rbPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col.range.a <- ryPal(10)[as.numeric(cut(results2$range.a,breaks = 10))]

results2$CompEst.ab <- results2$CompEst.a * results2$CompEst.b
results2$LVEst.ab <- results2$LVEst.a * results2$LVEst.b
results2$std.CompEst.ab <- results2$std.CompEst.a * results2$std.CompEst.b
results2$std.LVEst.ab <- results2$std.LVEst.a * results2$std.LVEst.b

results2$sd.a <- sqrt(results2$var.a)

jpeg("sim5a.Plot1.jpg", width = 12, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, 0, 1))
  plot(results2$range.a, results2$std.CompEst.a, type = "p", xlab = "range of X --> m paths",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$range.a, results2$std.LVEst.a, col = results2$col.var.a, pch = ".", cex = 1)
  
  par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
  plot(results2$range.a, results2$std.CompEst.b, type = "p", xlab = "range of X --> m paths",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$range.a, results2$std.LVEst.b, col = results2$col.var.a, pch = ".", cex = 1) 
  legend(x = 0, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
  
  par(fig = c(2/3, 1, 0, 1), new = TRUE)
  plot(results2$range.a, results2$std.CompEst.ab, type = "p", 
       xlab = "range of X --> m paths",
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$range.a, results2$std.LVEst.ab, col = results2$col.var.a, pch = ".", cex = 1)
dev.off()

#plot by variance instead of range: 
jpeg("sim5a.Plot2.jpg", width = 12, height = 8, units = "in", res = 500)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, .5, 1))
  plot(results2$var.a, results2$std.CompEst.a, type = "p", xlab = "variance of X --> m paths",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$var.a, results2$std.LVEst.a, col = results2$col.range.a, pch = ".", cex = 1)
  
  par(fig = c(1/3, 2/3, .5, 1), new = TRUE)
  plot(results2$var.a, results2$std.CompEst.b, type = "p", xlab = "variance of X --> m paths",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$var.a, results2$std.LVEst.b, col = results2$col.range.a, pch = ".", cex = 1) 
  legend(x = 0, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkorange"), bty = "n")
  
  par(fig = c(2/3, 1, .5, 1), new = TRUE)
  plot(results2$var.a, results2$std.CompEst.ab, type = "p", 
       xlab = "variance of X --> m paths",
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
  points(results2$var.a, results2$std.LVEst.ab, col = results2$col.range.a, pch = ".", cex = 1)

  #add pop RMSEA to the same figure: 
  results2$pop.RMSEA <- sqrt(results2$LV.fmin / 9)

  par(fig = c(0, 1, 0, .5), new = TRUE)
  plot(results2$var.a, results2$pop.RMSEA, type = "p", xlab = "variance of X --> m paths",
       ylab = "Population RMSEA", col = results2$col.range.a, ylim = c(0, 1), 
       pch = ".", cex = 1)
dev.off()


jpeg("sim5a.Plot3.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$range.a, results2$LV.fmin, type = "p", xlab = "range of X --> m paths",
     ylab = "min fit function LV model", col = results2$col.var.a, ylim = c(0, 1), 
     pch = ".", cex = 1)
dev.off()

jpeg("sim5a.Plot4.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$var.a, results2$LV.fmin, type = "p", xlab = "variance of X --> m paths",
     ylab = "min fit function LV model", col = results2$col.range.a, ylim = c(0, 1), 
     pch = ".", cex = 1)
dev.off()

#### SIM 5b ###################################
#choose fixed values of var.a, get a range of covm; vice versa. 

#given 3 values, pick 4th to achieve a particular variance! 
var.a <- .03

avec <- runif(3, 0, .8) 
avec

avec #[1] 0.01047671 0.69576863 0.47189423 0.68133162
var.a <- var(avec) #[1] 0.1022308


#quadratic formula to solve for a4 given var(ai): 
qa <- 3/4
qb <- (-1/2)*(avec[1] + avec[2] + avec[3])
qc <-  3/4*(avec[1]^2 + avec[2]^2 + avec[3]^2) - (1/2)*(avec[1]*avec[2] + avec[1]*avec[3] + avec[2]*avec[3]) - 3*var.a

a4 <- (-qb+sqrt(qb^2-(4*qa*qc)))/(2*qa)
a4 <- (-qb-sqrt(qb^2-(4*qa*qc)))/(2*qa)

sim5b <- function(i){
  a <- .6 
  b <- .4
  goal <- .03 #what we want var.a to be
  vara = 1
  while (vara > .03){ #fix range to .5
    avec <- runif(3, 0, .8) 
    qa <- 3/4
    qb <- (-1/2)*(avec[1] + avec[2] + avec[3])
    qc <-  3/4*(avec[1]^2 + avec[2]^2 + avec[3]^2) - (1/2)*(avec[1]*avec[2] + avec[1]*avec[3] + avec[2]*avec[3]) - 3*goal
    a4 <- (-qb+sqrt(qb^2-(4*qa*qc)))/(2*qa)
    if (is.na(a4) == FALSE) {
      avec <- c(avec, a4)   
      vara <- var(avec)      
    }
  }
  meana <- mean(avec)
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = 0)
  
  if (min(eigen(data[-5,-5])$values) > 0) { #ensure that sigma is positive definite! 
    amat <- avec %*% t(avec)
    #gamma <- sqrt(a/sum(amat)) #WRONG
    gamma <- a/sum(avec)    
    covm <- (1/(12*gamma^2)) - 1/3
    
    fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000)  
    fitLV <- sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE, warn = FALSE)  
    #fitFull <- sem(model = FullMod, sample.cov = data, sample.nobs = 100000, sample.cov.rescale = FALSE)
    
    conv.LV <- try(inspect(fitLV, "converged"))
    SEs.LV <- TRUE  
    if (class(try(vcov(fitLV))) == "try-error") {SEs.LV <- FALSE}
    if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  
    
    return(c(meana, rangea, vara, covm, gamma,
             inspect(fitLV, "coef")$lambda[1:4,1],
             inspect(fitComp, "coef")$beta[1,3], #a path composite model
             inspect(fitComp, "coef")$beta[2,1], #b path composite model  
             inspect(fitLV, "coef")$beta[1,3], #a path LV model
             inspect(fitLV, "coef")$beta[2,1], #b path LV model
             inspect(fitLV, "std")$lambda[1:4,1],
             inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
             inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
             inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
             inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
             inspect(fitComp, "fit")["fmin"], #fit of composite model
             conv.LV, SEs.LV, fit.LV,
             inspect(fitComp, "se")$beta[1,3], #SE a composite model
             inspect(fitComp, "se")$beta[2,1], #SE b composite model
             inspect(fitLV, "se")$beta[1,3], #SE a LV model
             inspect(fitLV, "se")$beta[2,1])) #SE b LV model
  } else {return (rep(NA, 29))}
}

X <- 1:100000
ptm <- proc.time()
out5b <- sapply(X, sim5b)
proc.time() - ptm

save(out5b, file = "sim5b.100000.Rdata")
#load(file = "sim5b.100000.Rdata")
results <- as.data.frame(t(out5b))

colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", 
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "std.lam1", "std.lam2", "std.lam3", "std.lam4",
                       "std.CompEst.a", "std.CompEst.b", "std.LVEst.a", "std.LVEst.b", 
                       "Comp.fmin", "conv.LV", "SEconv.LV", "LV.fmin",
                       "CompSE.a", "CompSE.b", "LVSE.a", "LVSE.b")

#fix gamma and covm values: 
results$gamma <- .6 /(4*results$mean.a)
results$covm <- (1/(12*results$gamma^2)) - 1/3
head(results)

#discard rows that did not converge or that generated non-positive-definite sigma: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #3879 rows left


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))

results2$col.range.a <- rbPal(10)[as.numeric(cut(results2$range.a,breaks = 10))]
results2$col.covm <- ryPal(10)[as.numeric(cut(results2$covm,breaks = 10))]

results2$std.CompEst.ab <- results2$std.CompEst.a * results2$std.CompEst.b
results2$std.LVEst.ab <- results2$std.LVEst.a * results2$std.LVEst.b

results2$pop.RMSEA <- sqrt(results2$LV.fmin/9)

#plot by covm
jpeg("sim5b.Plot1.jpg", width = 12, height = 8, units = "in", res = 500)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, .5, 1))
  plot(results2$covm, results2$std.CompEst.a, type = "p", xlab = "correlation between causal indicators",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), xlim = c(0, 1), pch = ".", cex = 1)
  points(results2$covm, results2$std.LVEst.a, col = results2$col.range.a, pch = ".", cex = 1)
  
  par(fig = c(1/3, 2/3, .5, 1), new = TRUE)
  plot(results2$covm, results2$std.CompEst.b, type = "p", xlab = "correlation between causal indicators",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), xlim = c(0, 1), pch = ".", cex = 1)
  points(results2$covm, results2$std.LVEst.b, col = results2$col.range.a, pch = ".", cex = 1) 
  legend(x = 0, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "purple4"), bty = "n")
  
  par(fig = c(2/3, 1, .5, 1), new = TRUE)
  plot(results2$covm, results2$std.CompEst.ab, type = "p", 
       xlab = "correlation between causal indicators",
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), xlim = c(0, 1), pch = ".", cex = 1)
  points(results2$covm, results2$std.LVEst.ab, col = results2$col.range.a, pch = ".", cex = 1)
  
  #add pop RMSEA to the plot: 
  par(fig = c(0, 1, 0, .5), new = TRUE)
  plot(results2$covm, results2$pop.RMSEA, type = "p", xlab = "correlation between causal indicators",
       ylab = "population RSMEA", col = results2$col.range.a, ylim = c(0, 1), xlim = c(0,1),
       pch = ".", cex = 1)
dev.off()



jpeg("sim5b.Plot2.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$LV.fmin, type = "p", xlab = "covariance between causal indicators",
     ylab = "min fit function LV model", col = results2$col.range.a, ylim = c(0, 1), xlim = c(0,1),
     pch = ".", cex = 1)
dev.off()


######################## SIM 6 ###################################################
#allow both a and gamma values to differ

sim6 <- function(i){
  a <- .6 #increase this to increase the range of possible values.  
  b <- .4

  data <- matrix(0, 8, 8)
  while (min(eigen(data[-c(5,8),-c(5,8)])$values) <=0) {
    avec <- runif(4, 0, .8) #in sim2 this had max = .5
    gen <- generateCov2(a = a, b = b, avec = avec, covm = NA, res.M = .25)
    data <- gen[[1]]
  }
    
  meana <- mean(avec)
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  gammavec <- gen[[2]]
  covm <- gen[[3]]
  
  fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000, warn = FALSE)  
  fitLV <- try(sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE, warn = FALSE), silent = TRUE)  
  
  conv.LV <- try(inspect(fitLV, "converged"), silent = TRUE)
  SEs.LV <- TRUE  
  if (identical(class(try(vcov(fitLV), silent = TRUE)), "try-error") == TRUE) {SEs.LV <- FALSE}
  if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  
  
  return(c(avec, gammavec, covm, 
           inspect(fitLV, "std")$lambda[1:4,1],
           inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
           inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
           inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
           inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
           inspect(fitComp, "fit")["fmin"], #fit of composite model
           conv.LV, SEs.LV, fit.LV))    
}

X <- 1:10000
ptm <- proc.time()
out <- sapply(X, sim6)
proc.time() - ptm

save(out, file = "sim6.10000.Rdata")
load(file = "sim6.10000.Rdata")
results <- as.data.frame(t(out))

colnames(results) <- c("a1", "a2", "a3", "a4", "gam1", "gam2", "gam3", "gam4", 
                       "covm", "lam1", "lam2", "lam3", "lam4",
                       "ComStd.a", "ComStd.b", "LVStd.a", "LVStd.b", 
                       "Com.fmin", "conv.LV", "SEconv.LV", "LV.fmin")
round(results[1:25,], 3)

#discard rows that did not converge or that generated non-positive-definite sigma: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #5652 rows left
results2$var.a <- apply(results2[,1:4], 1, var)
results2$var.gam <- apply(results2[,5:8], 1, var)


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))
gyPal <- colorRampPalette(c('darkgreen','yellow'))


results2$col.var.a <- rbPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col2.var.a <- gyPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col.var.gam <- rbPal(10)[as.numeric(cut(results2$var.gam,breaks = 10))]
results2$col2.var.gam <- gyPal(10)[as.numeric(cut(results2$var.gam,breaks = 10))]
results2$col.covm <- ryPal(10)[as.numeric(cut(results2$covm,breaks = 10))]

jpeg("sim6.Plot3.jpg", width = 12, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, 1/3, 0, 1))
plot(results2$covm, results2$ComStd.a, type = "p", xlab = "covariance between causal indicators",
     ylab = "estimate of a path", col = results2$col2.var.gam, ylim = c(0, 1), pch = ".", cex = 2)
points(results2$covm, results2$LVStd.a, col = results2$col.var.gam, pch = ".", cex = 2)
abline(h = .6)

par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
plot(results2$covm, results2$ComStd.b, type = "p", xlab = "covariance between causal indicators",
     ylab = "estimate of b path", col = results2$col2.var.gam, ylim = c(0, 1), pch = ".", cex = 2)
points(results2$covm, results2$LVStd.b, col = results2$col.var.gam, pch = ".", cex = 2) 
abline(h = .4)
legend(x = -.2, y = 1, c("composite model", "reflective model"), 
       fill = c("darkgreen", "darkred"), bty = "n")

#plot ab path estimate #
results2$LVEst.ab <- results2$LVEst.a * results2$LVEst.b
#jpeg("sim6.Plot6.jpg", width = 5, height = 4, units = "in", res = 800)
par(fig = c(2/3, 1, 0, 1), new = TRUE)
plot(results2$covm, results2$ComStd.a*results2$ComStd.b, type = "p", 
     xlab = "covariance between causal indicators",
     ylab = "estimate of a*b path", col = results2$col2.var.gam, ylim = c(0, .4), pch = ".")
points(results2$covm, results2$LVStd.a*results2$LVStd.b, col = results2$col.var.gam, pch = ".")
abline(h = .24)
dev.off()

#plot fit  #
jpeg("sim6.Plot4.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$LV.fmin, type = "p", xlab = "covariance between causal indicators",
     ylab = "min fit function LV model", col = results2$col.var.gam, ylim = c(0, 1), 
     pch = ".", cex = 2)
points(results2$covm, results2$Com.fmin, col = results2$col2.var.gam, pch = ".")
dev.off()


######################## SIM 3b ###################################################
# now indicators are unreliable so composite model is not an exact fit

sim3b <- function(i){
  a <- .6 
  b <- .4
  res.M <- 0
  avec <- runif(4, 0, .8) #in sim2 this had max = .5
  data <- generateCov3(a = a, b = b, avec = avec, covm = NA, res.M = res.M, r = .8)  
  while (min(eigen(data[-c(5,8),-c(5,8)])$values) < 1e-10) { #ensure that sigma is positive definite! 
    avec <- runif(4, 0, .8) #in sim2 this had max = .5
    data <- generateCov3(a = a, b = b, avec = avec, covm = NA, res.M = res.M)
  }
  meana <- mean(avec)
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  amat <- avec %*% t(avec)
  gamma <- a/sum(avec)
  covm <- (1 - 4*gamma^2 - res.M)/(12*gamma^2)   #NEW so total var(M) includes residual variance
  
  fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000)  
  fitLV <- sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  
  #fitFull <- sem(model = FullMod, sample.cov = data, sample.nobs = 100000, sample.cov.rescale = FALSE)
  
  conv.LV <- try(inspect(fitLV, "converged"))
  SEs.LV <- TRUE  
  if (class(try(vcov(fitLV))) == "try-error") {SEs.LV <- FALSE}
  if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  
  
  return(c(meana, rangea, vara, covm, gamma,
           inspect(fitLV, "coef")$lambda[1:4,1],
           inspect(fitComp, "coef")$beta[1,3], #a path composite model
           inspect(fitComp, "coef")$beta[2,1], #b path composite model  
           inspect(fitLV, "coef")$beta[1,3], #a path LV model
           inspect(fitLV, "coef")$beta[2,1], #b path LV model
           inspect(fitLV, "std")$lambda[1:4,1],
           inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
           inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
           inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
           inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
           inspect(fitComp, "fit")["fmin"], #fit of composite model
           conv.LV, SEs.LV, fit.LV,
           inspect(fitComp, "se")$beta[1,3], #SE a composite model
           inspect(fitComp, "se")$beta[2,1], #SE b composite model
           inspect(fitLV, "se")$beta[1,3], #SE a LV model
           inspect(fitLV, "se")$beta[2,1])) #SE b LV model))    
}

X <- 1:100000
ptm <- proc.time()
out <- sapply(X, sim3b)
time3b <- proc.time() - ptm #20 hours? 

save(out, file = "sim3b.100000.Rdata")
load(file = "sim3b.10000.Rdata")
results <- as.data.frame(t(out))

colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", 
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "std.lam1", "std.lam2", "std.lam3", "std.lam4",
                       "std.CompEst.a", "std.CompEst.b", "std.LVEst.a", "std.LVEst.b", 
                       "Comp.fmin", "conv.LV", "SEconv.LV", "LV.fmin",
                       "CompSE.a", "CompSE.b", "LVSE.a", "LVSE.b")

#discard rows that did not converge or that generated non-positive-definite sigma: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) # 54057 rows left

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))

results2$col.var.a <- rbPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col.covm <- ryPal(10)[as.numeric(cut(results2$covm,breaks = 10))]

jpeg("sim3b.Plot1.jpg", width = 12, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, 1/3, 0, 1))
plot(results2$var.a, results2$std.CompEst.a, type = "p", xlab = "variance of X --> m paths",
     ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
points(results2$var.a, results2$std.LVEst.a, col = results2$col.covm, pch = ".", cex = 1)
abline(h = .6, col = "black")

par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
plot(results2$var.a, results2$std.CompEst.b, type = "p", xlab = "variance of X --> m paths",
     ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
points(results2$var.a, results2$std.LVEst.b, col = results2$col.covm, pch = ".", cex = 1) 
abline(h = .4, col = "black")
legend(x = 0, y = 1, c("composite model", "reflective model"), 
       fill = c("darkgreen", "darkred"), bty = "n")

#plot ab path estimate
results2$std.LVEst.ab <- results2$std.LVEst.a * results2$std.LVEst.b
par(fig = c(2/3, 1, 0, 1), new = TRUE)
plot(results2$var.a, results2$std.CompEst.a*results2$std.CompEst.b, type = "p", 
     xlab = "variance of X --> m paths",
     ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
points(results2$var.a, results2$std.LVEst.ab, col = results2$col.covm, pch = ".", cex = 1)
abline(h = .24, col = "black")
dev.off()

jpeg("sim3b.Plot3.jpg", width = 12, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, 1/3, 0, 1))
plot(results2$covm, results2$std.CompEst.a, type = "p", xlab = "correlation between causal indicators",
     ylab = "estimate of a path", col = "darkgreen", xlim = c(0, 1), ylim = c(0, 1), pch = ".", cex = 1)
points(results2$covm, results2$std.LVEst.a, col = results2$col.var.a, pch = ".", cex = 1)
abline(h = .6, col = "black")

par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
plot(results2$covm, results2$std.CompEst.b, type = "p", xlab = "correlation between causal indicators",
     ylab = "estimate of b path", col = "darkgreen", xlim = c(0, 1), ylim = c(0, 1), pch = ".", cex = 1)
points(results2$covm, results2$std.LVEst.b, col = results2$col.var.a, pch = ".", cex = 1) 
abline(h = .4, col = "black")
legend(x = 0, y = 1, c("composite model", "reflective model"), 
       fill = c("darkgreen", "darkred"), bty = "n")

#plot ab path estimate
results2$std.LVEst.ab <- results2$std.LVEst.a * results2$std.LVEst.b
par(fig = c(2/3, 1, 0, 1), new = TRUE)
plot(results2$covm, results2$std.CompEst.a*results2$std.CompEst.b, type = "p", 
     xlab = "correlation between causal indicators", xlim = c(0, 1),
     ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
points(results2$covm, results2$std.LVEst.ab, col = results2$col.var.a, pch = ".", cex = 1)
abline(h = .24, col = "black")
dev.off()

#plot by mean a path instead of covm
jpeg("sim3b.Plot3b.jpg", width = 12, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, 1/3, 0, 1))
plot(results2$mean.a, results2$CompEst.a, type = "p", xlab = "mean a path",
     ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 6), pch = ".", cex = 1)
points(results2$mean.a, results2$LVEst.a, col = results2$col.var.a, pch = ".", cex = 1)
abline(h = .6, col = "black")

par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
plot(results2$mean.a, results2$CompEst.b, type = "p", xlab = "mean a path",
     ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 1)
points(results2$mean.a, results2$LVEst.b, col = results2$col.var.a, pch = ".", cex = 1) 
legend(x = 0, y = 1, c("composite model", "reflective model"), 
       fill = c("darkgreen", "darkred"), bty = "n")
abline(h = .4, col = "black")

#plot ab path estimate
results2$LVEst.ab <- results2$LVEst.a * results2$LVEst.b
#jpeg("sim3b.Plot6.jpg", width = 5, height = 4, units = "in", res = 800)
#par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(2/3, 1, 0, 1), new = TRUE)
plot(results2$mean.a, results2$CompEst.a*results2$CompEst.b, type = "p", 
     xlab = "mean a path",
     ylab = "estimate of a*b path", col = "darkgreen", ylim = c(.18, .32), pch = ".", cex = 1)
points(results2$mean.a, results2$LVEst.ab, col = results2$col.var.a, pch = ".", cex = 1)
abline(h = .24, col = "black")
dev.off()


jpeg("sim3b.Plot4.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$LV.fmin, type = "p", xlab = "covariance between causal indicators",
     ylab = "min fit function LV model", col = results2$col.var.a, ylim = c(0, 1), 
     pch = ".", cex = 1)
dev.off()

#plot fit in terms of population RMSEA instead of fmin: 
results2$pop.RMSEA <- sqrt(results2$LV.fmin/9) 
results2$pop.RMSEA.Comp <- sqrt(results2$Comp.fmin) 
jpeg("sim3b.Plot4.fit.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$pop.RMSEA, type = "p", xlab = "correlation between causal indicators",
     ylab = "Population RMSEA", col = results2$col.var.a, ylim = c(0, 1), 
     pch = ".", cex = 1)
points(results2$covm, results2$pop.RMSEA.Comp, col = "darkgreen", pch = ".", cex = 1)
dev.off()


#plot range of a1-a4 values by covm: 
cor(results2$covm, results2$range.a) #-.287 
cor(results2$mean.a, results2$range.a) #-.287 
jpeg("sim3b.Plot5.jpg", width = 9, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$mean.a, results2$range.a, pch = ".", cex = 1, col = results2$col.var.a, 
     xlab = "covariance among causal indicators", ylab = "range of X -> m path values")
dev.off()


#plot fit as a function of the range of X -> m values, and covm (colour)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
jpeg("sim3b.Plot7.jpg", width = 9, height = 4, units = "in", res = 800)
plot(results2$range.a, results2$LV.fmin, type = "p", xlab = "range in a1-a4 paths",
     ylab = "min fit function LV model", col = results2$col.covm, ylim = c(0, 1), 
     pch = ".", cex = 1)
dev.off()

#plot fit as a function of the range of X -> m values, and covm (colour)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
jpeg("sim3b.Plot8.jpg", width = 9, height = 4, units = "in", res = 800)
plot(results2$LVEst.a*results2$LVEst.b, results2$LV.fmin, pch = ".", 
     col = results2$col.var.a, cex = 1, xlab = "a*b path estimate", ylab = "fmin")
dev.off()

