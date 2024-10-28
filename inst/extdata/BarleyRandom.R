#########################################
###
###  Two-stage
###
###  random at stage 1
###
###  Ari Verbyla
###
#########################################

require(asreml)
asreml.options(Cfixed=TRUE, extra=3)
require(lattice)
require(lmmtools)

data(BarleyData)

#######  Dumb gamma parameterisation causes problems
#######  vm must use variance matrix not inverse, the latter fails
#######  IMPORTANT: must use ar1v to ensure a variance is part of the residual model
#######  s1trace = TRUE will print out things that are happening along the way

s11 <- stage1(yield ~ 1 + lin(Row) + lin(Column), random= ~ Genotype + Block + Row + Column,
              residual = ~ ar1v(Column):ar1(Row), Trait="Site", n.trait=1, s1trace=TRUE,
              Genetic="Genotype", data = BarleyData, predict.options=list(ignore=c("lin(Row)", "lin(Column)")))

####  Won't work W11 <- as.matrix(s11$bdWt)  # ()
V11 <- as.matrix(s11$bdVt)

##############  diag model

s21.sv <- asreml(yield ~ Site-1, random = ~ diag(Site):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["diag"]] <- asreml(yield ~ Site-1, random = ~ diag(Site):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=200, G.param=g21, workspace="20Gb")

save.image()

##############  fa1 model

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 1):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa1"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 1):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=200, G.param=g21, workspace="20Gb")

save.image()

##############  fa2 model

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 2):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa2"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 2):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=200, G.param=g21, workspace="20Gb")

save.image()

##############  fa3 model

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 3):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa3"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 3):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=200, G.param=g21, workspace="20Gb")

save.image()

##############  fa4 model

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 4):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa4"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 4):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=200, G.param=g21, workspace="20Gb")

save.image()

##############  fa5

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 5):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa5"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 5):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=100, G.param=g21, workspace="20Gb")

save.image()

##############  fa6 model

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 6):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa6"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 6):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=100, G.param=g21, workspace="20Gb")

save.image()

##############  fa7 model

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 7):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa7"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 7):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=100, G.param=g21, workspace="20Gb")

save.image()

##############  fa8 model

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 8):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa8"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 8):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=100, G.param=g21, workspace="20Gb")

save.image()

##############  fa9 model

s21.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 9):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11$pred.df, start.values=TRUE)

g21 <- s21.sv$vparameters.table
g21
g21$Value[grep("vm\\(", g21$Component)] <- 1
g21$Constraint[grep("vm\\(", g21$Component)] <- "F"
g21

s21.asr[["fa9"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 9):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  maxit=100, G.param=g21, workspace="20Gb")

##################### Testing

loglik <- unlist(lapply(s21.asr, function(el) el$loglik))
summ.bnd <- lapply(s21.asr, function(el) summary(el)$varcomp[, "bound"])
df1  <- diff(unlist(lapply(summ.bnd, function(el) sum(el != "B" & el != "F"))))
df <- seq(24, by=-1, length.out=9)
remlrt <- 2*diff(loglik)
pval <- 1-pchisq(remlrt, df)
pval

pval1 <- 1-pchisq(remlrt, df1)
pval1

res <- icREML(s21.asr)
fasum <- faSummary(s21.asr, s11$pred.df, "Site", id="Genotype")
fasum[, 1:6]

################  Predictions

tmp <- s21.asr[["fa8"]]

s21.pred <- predict(tmp, classify="Site:Genotype",
                    workspace="10Gb", pworkspace="10Gb", maxit=1)

save.image()

gparam <- s21.asr[["fa8"]]$vparameters[-c(1, 218)]
gparam

psi <- gparam[1:24]
Lambda <- matrix(gparam[25:216], ncol=8, byrow=FALSE)
Gg <- Lambda %*% t(Lambda) + diag(psi)
round(Gg*100,2)
dimnames(Gg) <- list(levels(s11$pred.df$Site),levels(s11$pred.df$Site))
round(cov2cor(Gg),2)

Gg2 <- Gg
save(Gg2, file="Gg2.rda")

#############  According to LRT and AIC

tmp6 <- s21.asr[["fa6"]]
tmp8 <- s21.asr[["fa8"]]

s216.pred <- predict(tmp6, classify="Site:Genotype",  only = "fa(Site, 6):Genotype",
                    workspace="10Gb", pworkspace="10Gb", maxit=1)

s218.pred <- predict(tmp8, classify="Site:Genotype", only = "fa(Site, 8):Genotype",
                    workspace="10Gb", pworkspace="10Gb", maxit=1)

save.image()

length(s21.asr[["fa6"]]$vparameters)

gparam6 <- s21.asr[["fa6"]]$vparameters[-c(1,170)]
gparam6

psi6 <- gparam6[1:24]
Lambda6 <- matrix(gparam6[25:216], ncol=6, byrow=FALSE)
Gg6 <- Lambda6 %*% t(Lambda6) + diag(psi6)
round(Gg6*100,2)
dimnames(Gg6) <- list(levels(s11$pred.df$Site),levels(s11$pred.df$Site))
round(cov2cor(Gg6),2)

Gg26 <- Gg6
save(Gg26, file="Gg26.rda")

length(s21.asr[["fa8"]]$vparameters)

gparam8 <- s21.asr[["fa8"]]$vparameters[-c(1, 218)]
gparam8

psi8 <- gparam8[1:24]
Lambda8 <- matrix(gparam8[25:216], ncol=8, byrow=FALSE)
Gg8 <- Lambda8 %*% t(Lambda8) + diag(psi8)
round(Gg8*100,2)
dimnames(Gg8) <- list(levels(s11$pred.df$Site),levels(s11$pred.df$Site))
round(cov2cor(Gg8),2)

Gg28 <- Gg8
save(Gg28, file="Gg28.rda")

###################  Diagonal weights

s21d.asr <- list()

s21d.asr[["diag"]] <- asreml(yield ~ Site-1, random = ~ diag(Site):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa1"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 1):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa2"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 2):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa3"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 3):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa4"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 4):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa5"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 5):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa6"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 6):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa7"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 7):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa8"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 8):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

s21d.asr[["fa9"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 9):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11$pred.df,
                   maxit=200, workspace="20Gb")

loglikd <- unlist(lapply(s21d.asr, function(el) el$loglik))
##df <- seq(24, by=-1, length.out=9)
remlrtd <- 2*diff(loglikd)
pvald <- 1-pchisq(remlrtd, df)
round(pvald,3)

summd.bnd <- lapply(s21d.asr, function(el) summary(el)$varcomp[, "bound"])
dfd1  <- diff(unlist(lapply(summd.bnd, function(el) sum(el != "B" & el != "F"))))
pvald1 <- 1-pchisq(remlrtd, dfd1)
round(pvald1,3)

resd <- icREML(s21d.asr)
fasumd <- faSummary(s21d.asr, s11$pred.df, "Site", id="Genotype")
fasumd[, 1:6]

################  Predictions

tmpd6 <- s21d.asr[["fa6"]]
tmpd8 <- s21d.asr[["fa8"]]

#######s21d6.pred  <- s21d.pred  ##### done earlier

s21d6.pred <- predict(tmpd6, classify="Site:Genotype", only = "fa(Site, 6):Genotype",
                      workspace="10Gb", pworkspace="10Gb", maxit=1)

s21d8.pred <- predict(tmpd8, classify="Site:Genotype", only = "fa(Site, 8):Genotype",
                    workspace="10Gb", pworkspace="10Gb", maxit=1)

save.image()

gparamd6 <- s21d.asr[["fa6"]]$vparameters[-169]
gparamd6

psid6 <- gparamd6[1:24]
Lambdad6 <- matrix(gparamd6[25:168], ncol=6, byrow=FALSE)
Ggd6 <- Lambdad6 %*% t(Lambdad6) + diag(psid6)
round(Ggd6*100,2)
dimnames(Ggd6) <- list(levels(s11$pred.df$Site),levels(s11$pred.df$Site))
round(cov2cor(Ggd6),2)

Gg2d6 <- Ggd6
save(Gg2d6, file="Gg2d6.rda")

length(s21.asr[["fa8"]]$vparameters)

gparamd8 <- s21d.asr[["fa8"]]$vparameters[-217]
gparamd8

psid8 <- gparamd8[1:24]
Lambdad8 <- matrix(gparamd8[25:216], ncol=8, byrow=FALSE)
Ggd8 <- Lambdad8 %*% t(Lambdad8) + diag(psid8)
round(Ggd8*100,2)
dimnames(Ggd8) <- list(levels(s11$pred.df$Site),levels(s11$pred.df$Site))
round(cov2cor(Ggd8),2)

Gg2d8 <- Ggd8
save(Gg2d8, file="Gg2d8.rda")

#############  Put together

rresults <- cbind.data.frame(s216.pred$pvals[, 1:4], s218.pred$pvals[, 3:4],
                             s21d6.pred$pvals[, 3:4], s21d8.pred$pvals[, 3:4])

names(rresults)[3:10] <- c("RandomFull2Stage (fa6)", "RandomFull2StageSE (fa6)",
                           "RandomFull2Stage (fa8)", "RandomFull2StageSE (fa8)",
                           "RandomDiag2Stage (fa6)", "RandomDiag2StageSE (fa6)",
                           "RandomDiag2Stage (fa8)", "RandomDiag2StageSE (fa8)")

pairs(rresults[, c(3,5,7,9)])
cor(rresults[, c(3,5,7, 9)])

##############  so standard errors  ...

write.table(rresults, file="BarleyRandom2stage1.csv", quote=FALSE, row.names=FALSE, sep=",")

#################################   Variance parameters

Rg.df <- data.frame(as.vector(Gg26[upper.tri(Gg26, diag=FALSE)]), as.vector(Gg28[upper.tri(Gg28, diag=FALSE)]),
                     as.vector(Gg2d6[upper.tri(Gg2d6, diag=FALSE)]), as.vector(Gg2d8[upper.tri(Gg2d8, diag=FALSE)]))

Vg.df <- data.frame(diag(Gg26), diag(Gg28), diag(Gg2d6), diag(Gg2d8))

names(Rg.df) <- names(Vg.df) <- c("RandomFull (fa6)", "RandomFull (fa8)", "RandomDiag (fa6)", "RandomDiag (fa8)")

head(Rg.df)
pairs(Rg.df)
cor(Rg.df)

save(Rg.df, file="Rgdf.rda")

head(Vg.df)
pairs(Vg.df)
cor(Vg.df)

save(Vg.df, file="Vgdf.rda")

#############################################

