############## p-rep trials

require(asreml)
require(lmmtools)
require(ASExtras4)
asreml::asreml.options(Cfixed=TRUE, extra=3, ai.sing=TRUE)

load("prepData.rda")

s11 <- stage1(yield ~ 1 + lin(Row) + lin(Column), random = ~ mpid + Block + Row + Column,
             residual = ~ ar1v(Column):ar1(Row),
             data=prepData, Trait="Site", Genetic = "mpid",
             specific = list(PhenoGroup = list(type="fixed", which=3:4)),
             predGenetic=list(id = "mpid", classify = "mpid", rm = list(list(terms = "PhenoGroup", which=3:4))),
             predict.options=list(ignore = c("lin(Row)", "lin(Column"), pworkspace="5GB"),
             na.action = na.method(x="include"), keep.models=TRUE,
             s1trace=TRUE, workspace = "5GB")

V11 <- as.matrix(s11$bdVt)

##############  models for stage 2

##############  factor analytic  fitted as reduced rank plus diag for both additive and residual genetic effects

s21.sv <- asreml(yield ~ Site-1, random = ~ corgh(Site):mpid + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  start.values=TRUE)

g2 <- s21.sv$vparameters.table
g2
g2$Value[grep("vm\\(TG", g2$Component)] <- 1
g2$Constraint[grep("vm\\(TG", g2$Component)] <- "F"
g2

s21.asr <- asreml(yield ~ Site-1, random = ~ corgh(Site):mpid + vm(TG, V11),
                  family = asr_gaussian(dispersion = 0.0001),
                  data=s11$pred.df,
                  G.param=g2, workspace="10Gb", maxit=30)

summary(s21.asr)$varcomp

##                              component  std.error   z.ratio bound %ch
##vm(TG, V11)                  1.00000000         NA        NA     F   0
##Site:mpid!Site!2:!Site!1.cor 0.07594544 0.05968759  1.272382     U   0
##Site:mpid!Site!3:!Site!1.cor 0.83152884 0.03600421 23.095323     U   0
##Site:mpid!Site!3:!Site!2.cor 0.43710182 0.04881442  8.954358     U   0
##Site:mpid!Site!4:!Site!1.cor 0.18279776 0.05813619  3.144302     U   0
##Site:mpid!Site!4:!Site!2.cor 0.78479584 0.03230489 24.293408     U   0
##Site:mpid!Site!4:!Site!3.cor 0.61519378 0.04058417 15.158466     U   0
##Site:mpid!Site_1             0.13702933 0.01084824 12.631479     P   0
##Site:mpid!Site_2             0.23073430 0.01742860 13.238835     P   0
##Site:mpid!Site_3             0.23960804 0.01751181 13.682654     P   0
##Site:mpid!Site_4             0.16426024 0.01136456 14.453731     P   0
##units!R                      0.00010000         NA        NA     F   0


####################  Diagonal weights

s21d.asr <- asreml(yield ~ Site-1, random = ~ corgh(Site):mpid,
                  family = asr_gaussian(dispersion = 1), weights = wt,
                  data=s11$pred.df,
                  workspace="10Gb", maxit=30)

summary(s21d.asr)$varcomp

##                              component  std.error   z.ratio bound %ch
##Site:mpid!Site!2:!Site!1.cor 0.09676608 0.05624874  1.720324     U   0
##Site:mpid!Site!3:!Site!1.cor 0.83444938 0.03505235 23.805804     U   0
##Site:mpid!Site!3:!Site!2.cor 0.41606340 0.04618932  9.007784     U   0
##Site:mpid!Site!4:!Site!1.cor 0.13820778 0.05670395  2.437357     U   0
##Site:mpid!Site!4:!Site!2.cor 0.74962332 0.03055817 24.531030     U   0
##Site:mpid!Site!4:!Site!3.cor 0.56269483 0.04101417 13.719523     U   0
##Site:mpid!Site_1             0.13549693 0.01074795 12.606765     P   0
##Site:mpid!Site_2             0.25898522 0.01824981 14.191116     P   0
##Site:mpid!Site_3             0.24553294 0.01777902 13.810263     P   0
##Site:mpid!Site_4             0.17879251 0.01190211 15.021912     P   0
##units!R                      1.00000000         NA        NA     F   0

###################  Heritabilities

h2 <- lapply(s11$models, function(el) genherit1(el, id = "mpid"))
h2

##$`1`
##    yield
##0.6464116
##
##$`2`
##    yield
##0.6324765
##
##$`3`
##    yield
##0.7477038
##
##$`4`
##    yield
##0.7866346

