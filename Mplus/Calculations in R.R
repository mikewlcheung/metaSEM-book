###################################################
### Chapter 6: Three-level meta-analysis
###################################################

library(metaSEM)

## Select the relevant variables to export
my.df <- Bornmann07[, c(3,1,4:6)]
## Standardize "year"
my.df$Year <- scale(my.df$Year)
Fellow <- ifelse(Bornmann07$Type=="Fellowship", yes=1, no=0)

D_Phy <- ifelse(Bornmann07$Discipline=="Physical sciences", yes=1, no=0)
D_Life <- ifelse(Bornmann07$Discipline=="Life sciences/biology", yes=1, no=0)
D_Soc<- ifelse(Bornmann07$Discipline=="Social sciences/humanities", yes=1, no=0)
D_Mul <- ifelse(Bornmann07$Discipline=="Multidisciplinary", yes=1, no=0)

C_USA <- ifelse(Bornmann07$Country=="United States", yes=1, no=0)
C_Aus <- ifelse(Bornmann07$Country=="Australia", yes=1, no=0)
C_Can <- ifelse(Bornmann07$Country=="Canada", yes=1, no=0)
C_Eur <- ifelse(Bornmann07$Country=="Europe", yes=1, no=0)
C_UK <- ifelse(Bornmann07$Country=="United Kingdom", yes=1, no=0)

my.df <- cbind(my.df, Fellow, D_Phy, D_Life, D_Soc, D_Mul, C_USA, C_Aus, C_Can, C_Eur, C_UK)

## Show a few cases
head(my.df)

## Write to an external file
write.table(my.df, "Bornmann07.dat", row.names=FALSE, col.names=FALSE, na="NA", sep="\t")


###################################################
### Chapter 9: Conducting meta-analysis with Mplus
###################################################
library(metaSEM)

## Select the effect sizes
y <- wvs94a[, c("lifesat","lifecon")]

## Convert it into a column of effect sizes
y <- matrix(t(y), ncol=1)

## Prepare the design matrix
X <- matrix(rep(c(1,0,0,1), nrow(wvs94a)), ncol=2, byrow=TRUE)

## Convert the known sampling covariance matrix into 
## a block diagonal matrix
V <- matrix2bdiag(wvs94a[, c("lifesat_var","inter_cov",
                              "lifecon_var")])

## Calculate the transformation matrix
W0.5 <- chol(solve(V))

## Calculate the transformed effect size
y_new <- W0.5 %*% y

## Calculate the transformed design matrix
X_new <- W0.5 %*% X

## Center gnp and divide it by 10000 to improve numerical stability
## Prepare the gnp
gnp <- scale(wvs94a$gnp/10000, scale=FALSE)

## Convert y into one row per study
y2 <- matrix(c(t(y_new)), ncol=2, byrow=TRUE)
x2 <- matrix(c(t(X_new)), ncol=4, byrow=TRUE)

my.wide <- cbind(y2, x2, gnp)
colnames(my.wide) <- c("y1", "y2", "y1f1", "y1f2", "y2f1", "y2f2", "gnp")

## Write it as a plain text for Mplus
## y2f1 is excluded in my.wide since it contains only 0.
## Missing values are represented by *
write.table(my.wide[,-5], "wvs94a.dat", sep=" ", na="*", row.names=FALSE, 
             col.names=FALSE)

