###########################
###
###  Wheat data
###  Gilmour et al. (1997)
###
###  Ari Verbyla
###
###########################

require(asreml)
asreml.options(Cfixed=TRUE)

data(WheatData, package="lmmtools")

## WheatData <- asreml.read.table("wheatData.csv", header=TRUE, sep=",")

str(WheatData)
head(WheatData)

### Models

fm0 <- asreml(yield ~ 1,
              random = ~ Variety + Block,
              data = WheatData)

fm1 <- asreml(yield ~ 1,
              random = ~ Variety + Block,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm2 <- asreml(yield ~ 1,
              random = ~ Variety + Block + units,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm2 <- update(fm2)

fm3 <- asreml(yield ~ 1,
              random = ~ Variety + Block + Row + Col + units,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm3 <- update(fm3)

fm4 <- asreml(yield ~ lin(Row) + lin(Col),
              random = ~ Variety + Block + Row + Col + units,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm5 <- asreml(yield ~ lin(Row) + lin(Col),
              random = ~ Variety + Block + Row + Col +
                  spl(col) + units,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm6 <- asreml(yield ~ Rowcode + Colcode,
              random = ~ Variety + Block + units,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm6 <- update(fm6)

fm7 <- asreml(yield ~ Rowcode + Colcode,
              random = ~ Variety + Block + Row + Col + units,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm8 <- asreml(yield ~ Rowcode + Colcode + lin(Row) + lin(Col),
              random = ~ Variety + Block + Row + Col + units,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm8 <- update(fm8)

fm9 <- asreml(yield ~ Rowcode + Colcode + lin(Row) + lin(Col),
              random = ~ Variety + Block + Row + Col + spl(col) + units,
              residual = ~ ar1(Col):ar1(Row),
              data = WheatData)

fm9 <- update(fm9)

########  AIC and BIC

fm <- list(fm0=fm0, fm1=fm1, fm2=fm2, fm3=fm3, fm4=fm4, fm5=fm5, fm6=fm6, fm7=fm7,
           fm8=fm8, fm9=fm9)

library(lmmtools)

selection <- icREML(fm)
selection

###   model    loglik p q      AIC      BIC    logdet
###1    fm0 -1758.165 1 3 3524.331 3539.527  8.367356
###2    fm1 -1600.114 1 5 3212.229 3235.023  8.268164
###3    fm2 -1568.298 1 6 3150.596 3177.189  9.678423
###4    fm3 -1563.750 1 8 3145.499 3179.691  8.343257
###5    fm4 -1562.153 3 8 3146.307 3188.097 14.604338
###6    fm5 -1556.860 3 9 3137.721 3183.310  9.934153
###7    fm6 -1554.353 5 6 3130.707 3172.497 34.700807
###8    fm7 -1554.282 5 8 3134.563 3183.951 35.098493
###9    fm8 -1552.382 7 8 3134.764 3191.750 41.239778
###10   fm9 -1545.603 7 9 3123.207 3183.992 33.945625

#######  Heritabilities

h2 <- unlist(lapply(fm, function(el) genherit1(el, id="Variety")))
h2
