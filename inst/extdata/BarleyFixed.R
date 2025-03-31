#########################################
###
###  Two-stage
###
###  Fixed at stage 1
###
###  Ari Verbyla
###
#########################################

require(asreml)
asreml::asreml.options(Cfixed=TRUE, extra=3)
require(lattice)
require(lmmtools)

data(BarleyData)

#######  Gamma parameterisation causes problems
#######  vm must use variance matrix not inverse, the latter fails at stage 2
#######  residual must include the variance, hence ar1v
#######  s1trace = TRUE will print out things that are happening along the way

s11f <- stage1(yield ~ Genotype-1 + lin(Row) + lin(Column), random= ~ Block + Row + Column,
              residual = ~ ar1v(Column):ar1(Row), Trait="Site", n.trait=1, s1trace=TRUE,
              Genetic="Genotype", data = BarleyData, type="fixed",
              predict.options=list(ignore=c("lin(Row)", "lin(Column)")))

names(s11f)

####  vm won't work with the inverse W11 <- as.matrix(s11f$bdWt) so

#####  Block Variance matrix of the estimated Genotype means

V11 <- as.matrix(s11f$bdVt)

#########################  Stage 2 analyses

s21f.asr <- list()

##############  diag model

s21f.sv <- asreml(yield ~ Site-1, random = ~ diag(Site):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["diag"]] <- asreml(yield ~ Site-1, random = ~ diag(Site):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=20, G.param=g21f, workspace="20Gb")

##s21f.asr <- update(s21f.asr, maxit=40)

summary(s21f.asr[["diag"]])$varcomp

tmp <- s21f.asr[["diag"]]

s21f.pred[["diag"]] <- predict(tmp, classify="Site:Genotype", only=
                    workspace="10Gb", pworkspace="10Gb")

save.image()

##############  fa1 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 1):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa1"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 1):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=100, G.param=g21f, workspace="20Gb")

tmp <- s21f.asr[["fa1"]]
tmp <- update(tmp, maxit=100)
s21f.asr[["fa1"]] <- tmp

lapply(s21f.asr, function(el) el$loglik)

save.image()

##############  fa2 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 2):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa2"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 2):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=200, G.param=g21f, workspace="20Gb")

lapply(s21f.asr, function(el) el$loglik)

save.image()

##############  fa3 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 3):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa3"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 3):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=200, G.param=g21f, workspace="20Gb")

tmp <- s21f.asr[["fa3"]]
tmp <- update(tmp)
s21f.asr[["fa3"]] <- tmp

lapply(s21f.asr, function(el) el$loglik)

save.image()

##############  fa4 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 4):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa4"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 4):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=200, G.param=g21f, workspace="20Gb")

tmp <- s21f.asr[["fa4"]]
tmp <- update(tmp)
s21f.asr[["fa4"]] <- tmp

lapply(s21f.asr, function(el) el$loglik)

save.image()

##############  fa5 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 5):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa5"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 5):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=200, G.param=g21f, workspace="20Gb")

tmp <- s21f.asr[["fa5"]]
tmp <- update(tmp)
s21f.asr[["fa5"]] <- tmp

lapply(s21f.asr, function(el) el$loglik)

save.image()

##############  fa6 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 6):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa6"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 6):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=200, G.param=g21f, workspace="20Gb")

tmp <- s21f.asr[["fa6"]]
tmp <- update(tmp)
s21f.asr[["fa6"]] <- tmp

lapply(s21f.asr, function(el) el$loglik)

save.image()

##############  fa7 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 7):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa7"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 7):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=200, G.param=g21f, workspace="20Gb")

tmp <- s21f.asr[["fa7"]]
tmp <- update(tmp)
s21f.asr[["fa7"]] <- tmp

lapply(s21f.asr, function(el) el$loglik)

save.image()

##############  fa8 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 8):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa8"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 8):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=200, G.param=g21f, workspace="20Gb")

tmp <- s21f.asr[["fa8"]]
tmp <- update(tmp)
s21f.asr[["fa8"]] <- tmp

loglik <- unlist(lapply(s21f.asr, function(el) el$loglik))

save.image()

##############  fa9 model

s21f.sv <- asreml(yield ~ Site-1, random = ~ fa(Site, 9):Genotype + vm(TG, V11),
                 family = asr_gaussian(dispersion = 0.0001),
                 data=s11f$pred.df, start.values=TRUE)

g21f <- s21f.sv$vparameters.table
g21f
g21f$Value[grep("vm\\(", g21f$Component)] <- 1
g21f$Constraint[grep("vm\\(", g21f$Component)] <- "F"
g21f

s21f.asr[["fa9"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 9):Genotype + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11f$pred.df,
                  maxit=200, G.param=g21f, workspace="20Gb")

tmp <- s21f.asr[["fa9"]]
tmp <- update(tmp)
s21f.asr[["fa9"]] <- tmp

loglik <- unlist(lapply(s21f.asr, function(el) el$loglik))

save.image()

res <- icREML(s21f.asr)
fasum <- faSummary(s21f.asr, s11f$pred.df, "Site", id="Genotype")
fasum[, 1:6]

##################### Testing

loglik <- unlist(lapply(s21f.asr, function(el) el$loglik))
summ.bnd <- lapply(s21f.asr, function(el) summary(el)$varcomp[, "bound"])
df1  <- diff(unlist(lapply(summ.bnd, function(el) sum(el != "B" & el != "F"))))
df <- seq(24, by=-1, length.out=9)
remlrt <- 2*diff(loglik)
pval <- 1-pchisq(remlrt, df1)
pval
pval1 <- 1-pchisq(remlrt, df1)
pval1

################  Prediction

tmp6 <- s21f.asr[["fa6"]]
tmp8 <- s21f.asr[["fa8"]]

s21f6.pred <- predict(tmp6, classify="Site:Genotype", only = "fa(Site, 6):Genotype",
                      workspace="10Gb", pworkspace="10Gb", maxit=1)

s21f6.pred <- cs21f6.pred

s21f8.pred <- predict(tmp8, classify="Site:Genotype", only = "fa(Site, 8):Genotype",
                    workspace="10Gb", pworkspace="10Gb", maxit=1)

gparam6 <- s21f.asr[["fa6"]]$vparameters[-c(1,170)]
gparam6

psi6 <- gparam6[1:24]
Lambda6 <- matrix(gparam6[25:168], ncol=6, byrow=FALSE)
Ggf6 <- Lambda6 %*% t(Lambda6) + diag(psi6)
round(Ggf6*100,2)
dimnames(Ggf6) <- list(levels(s11f$pred.df$Site),levels(s11f$pred.df$Site))
round(cov2cor(Ggf6),2)

Gg2f6 <- Ggf6
save(Gg2f6, file="Gg2f6.rda")

gparam8 <- s21f.asr[["fa8"]]$vparameters[-c(1,218)]
gparam8

psi8 <- gparam8[1:24]
Lambda8 <- matrix(gparam8[25:216], ncol=8, byrow=FALSE)
Ggf8 <- Lambda8 %*% t(Lambda8) + diag(psi8)
round(Ggf8*100,2)
dimnames(Ggf8) <- list(levels(s11f$pred.df$Site),levels(s11f$pred.df$Site))
round(cov2cor(Ggf8),2)

Gg2f8 <- Ggf8
save(Gg2f8, file="Gg2f8.rda")


###################  Diagonal weights

s21fd.asr <- list()

s21fd.asr[["diag"]] <- asreml(yield ~ Site-1, random = ~ diag(Site):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa1"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 1):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa2"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 2):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa3"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 3):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa4"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 4):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa5"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 5):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa6"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 6):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa7"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 7):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa8"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 8):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

s21fd.asr[["fa9"]] <- asreml(yield ~ Site-1, random = ~ fa(Site, 9):Genotype,
                   weights = wt, family = asr_gaussian(dispersion = 1),
                   data=s11f$pred.df,
                   maxit=200, workspace="10Gb")

loglikd <- unlist(lapply(s21fd.asr, function(el) el$loglik))
df <- seq(24, by=-1, length.out=9)
remlrtd <- 2*diff(loglikd)
pvald <- 1-pchisq(remlrtd, df)
pvald

summd.bnd <- lapply(s21fd.asr, function(el) summary(el)$varcomp[, "bound"])
dfd1  <- diff(unlist(lapply(summd.bnd, function(el) sum(el != "B" & el != "F"))))
pvald1 <- 1-pchisq(remlrtd, dfd1)
pvald1

resd <- icREML(s21fd.asr)
fasumd <- faSummary(s21fd.asr, s11f$pred.df, "Site", id="Genotype")
fasumd[, 1:6]

tmpd6 <- s21fd.asr[["fa6"]]
tmpd8 <- s21fd.asr[["fa8"]]

s21fd6.pred <- predict(tmpd6, classify="Site:Genotype", only = "fa(Site, 6):Genotype",
                    workspace="10Gb", pworkspace="10Gb", maxit=1)

s21fd8.pred <- predict(tmpd8, classify="Site:Genotype", only = "fa(Site, 8):Genotype",
                    workspace="10Gb", pworkspace="10Gb", maxit=1)

save.image()

gparamd6 <- s21fd.asr[["fa6"]]$vparameters[-169]
gparamd6

psid6 <- gparamd6[1:24]
Lambdad6 <- matrix(gparamd6[25:168], ncol=6, byrow=FALSE)
Ggfd6 <- Lambdad6 %*% t(Lambdad6) + diag(psid6)
round(Ggfd6*100,2)
dimnames(Ggfd6) <- list(levels(s11f$pred.df$Site),levels(s11f$pred.df$Site))
round(cov2cor(Ggfd6),2)

Gg2fd6 <- Ggfd6
save(Gg2fd6, file="Gg2fd6.rda")

gparamd8 <- s21fd.asr[["fa8"]]$vparameters[-217]
gparamd8

psid8 <- gparamd8[1:24]
Lambdad8 <- matrix(gparamd8[25:216], ncol=8, byrow=FALSE)
Ggfd8 <- Lambdad8 %*% t(Lambdad8) + diag(psid8)
round(Ggfd8*100,2)
dimnames(Ggfd8) <- list(levels(s11f$pred.df$Site),levels(s11f$pred.df$Site))
round(cov2cor(Ggfd8),2)

Gg2fd8 <- Ggfd8
save(Gg2fd8, file="Gg2fd8.rda")

panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y, use="complete.obs"), digits=4)
    txt <- paste0(r)
    cex.cor <- 1.25 #/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 1)
}

pairs(Ggf.df,
      lower.panel = panel.cor,
      upper.panel = upper.panel)

write.table(Ggf.df, file="Ggfdf.csv", quote=FALSE, row.names=FALSE, sep=",")


#############  Put together

fresults <- cbind.data.frame(s21f6.pred$pvals[, 1:4], s21f8.pred$pvals[, 3:4],
                             s21fd6.pred$pvals[, 3:4], s21fd8.pred$pvals[, 3:4])

names(fresults)[3:10] <- c("FixedFull2Stage (fa6)", "FixedFull2StageSE (fa6)",
                           "FixedFull2Stage (fa8)", "FixedFull2StageSE (fa8)",
                           "FixedDiag2Stage (fa6)", "FixedDiag2StageSE (fa6)",
                           "FixedDiag2Stage (fa8)", "FixedDiag2StageSE (fa8)")

pairs(fresults[, c(3,5,7,9)])
cor(fresults[, c(3,5,7, 9)])

##############  so standard errors  ...

write.table(fresults, file="BarleyFixed2stage1.csv", quote=FALSE, row.names=FALSE, sep=",")

##############

Rgf.df <- data.frame(as.vector(Gg2f6[upper.tri(Gg2f6, diag=FALSE)]), as.vector(Gg2f8[upper.tri(Gg2f8, diag=FALSE)]),
                     as.vector(Gg2fd6[upper.tri(Gg2fd6, diag=FALSE)]), as.vector(Gg2fd8[upper.tri(Gg2fd8, diag=FALSE)]))

Vgf.df <- data.frame(diag(Gg2f6), diag(Gg2f8), diag(Gg2fd6), diag(Gg2fd8))

names(Rgf.df) <- names(Vgf.df) <- c("FixedFull (fa6)", "FixedFull (fa8)", "FixedDiag (fa6)", "FixedDiag (fa8)")

head(Rgf.df)
pairs(Rgf.df)
cor(Rgf.df)

save(Rgf.df, file="Rgfdf.rda")

head(Vgf.df)
pairs(Vgf.df)
cor(Vgf.df)

save(Vgf.df, file="Vgfdf.rda")

