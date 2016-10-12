#re: cole and preacher

#what if the true LV is not reflective but formative? 

#model: x --> m --> y, where x and y are observed and m is formative. 

rm(list = ls())
setwd("C:/Users/Mijke/Documents/__STUDIES__/CausalIndicators/MediationModelSimulation")


#regression paths a1-a4: causal indicators of M on X, vectorized:
avec = c(.5, .3, .2, 0)
amat <- avec %*% t(avec)

#res.M is the residual variance of M (i.e., if it is more than just the sum of m1-m4)
res.M <- 0

#variance of residuals for the causal indicators of M:
evec <- 1 - avec^2

#solve for gamma and total cov. between m1-m4 based on a (total x -> m path)
a <- .4
#gamma <- sqrt(a/sum(amat))   
gamma <- a/sum(avec)
covm <- (1/(12*gamma^2)) - 1/3    

#residual covariances among m1-m4: 
emat <- covm - amat  


#y on m
b <- .4

# genmod <- '
#            m1 ~ .5*X
#            m2 ~ .3*X
#            m3 ~ .2*X
#            m4 ~ 0*X
#            m1 ~~ .75*m1
#            m2 ~~ .91*m2
#            m3 ~~ .96*m3
#            m4 ~~ 1*m4
#            X ~~ 1*X
#            M ~ .3627381*m1
#            M ~ .3627381*m2
#            M ~ .3627381*m3
#            M ~ .3627381*m4
#            m1 ~~ .15*m2
#            m1 ~~ .2*m3
#            m1 ~~ .3*m4
#            m2 ~~ .24*m3
#            m2 ~~ .3*m4
#            m3 ~~ .3*m4
#            M ~~ 0*M
#            Y ~ .4*M
#            Y ~~ .84*Y'

genmod <-  paste(c("m1 ~ ", avec[1], "*X
                    m2 ~ ", avec[2], "*X 
                    m3 ~ ", avec[3], "*X
                    m4 ~ ", avec[4], "*X
                    m1 ~~ ", evec[1], "*m1
                    m2 ~~ ", evec[2], "*m2 
                    m3 ~~ ", evec[3], "*m3
                    m4 ~~ ", evec[4], "*m4  
                    X ~~ 1*X 
                    M ~ ", gamma, "*m1 
                    M ~ ", gamma, "*m2 
                    M ~ ", gamma, "*m3  
                    M ~ ", gamma, "*m4 
                    m1 ~~ ", emat[1,2], "*m2 
                    m1 ~~ ", emat[1,3], "*m3 
                    m1 ~~ ", emat[1,4], "*m4  
                    m2 ~~ ", emat[3,2], "*m3 
                    m2 ~~ ", emat[2,4], "*m4  
                    m3 ~~ ", emat[3,4], "*m4 
                    M ~~ ", res.M, "*M 
                    Y ~ ", b, "*M 
                    Y ~~ ", (1-b^2), "*Y"), sep = "")  

genData <- simulateData(model = genmod, model.type = "sem", meanstructure = "FALSE", 
                          return.type = "cov", empirical = TRUE, sample.cov.rescale = FALSE)
colnames(genData) <- c("m1", "m2", "m3", "m4", "M", "Y", "X")
round(genData, 2)


# generate cov matrix by hand: 
Beta <- matrix(0, 7, 7)
Beta[1:4,7] <- avec
Beta[5,1:4] <- gamma
Beta[6,5] <- b
I <- diag(7)
Psi <- matrix(0, 7, 7)
Psi[1:4, 1:4] <- emat
diag(Psi) <- c(evec, res.M, (1-b^2), 1)  

Sigma <- solve(I - Beta) %*% Psi %*% t(solve(I - Beta))
colnames(Sigma) <- c("m1", "m2", "m3", "m4", "M", "Y", "X")
round(Sigma, 2)

fitmod <- 'M ~ X
           Y ~ M'

fitData <- sem(model = fitmod, sample.cov = Sigma, sample.nobs = 100000)
summary(fitData, std= TRUE) #m on x is .363, we want it to be .4! 


LVmod <- 'LM =~ m1 + m2 + m3 + m4
          LM ~ X
          Y ~ LM'
fitLV <- sem(model = LVmod, sample.cov = genData, sample.nobs = 100000, std.lv = TRUE)
summary(fitLV)


generateCov(a = .4, b = .4, avec = c(.3, .3, .3, .3), covm = NA, res.M = 0)
