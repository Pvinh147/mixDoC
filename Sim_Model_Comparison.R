# Philip Vinh
# DoC Model Comparisons

# Load libraries
library(OpenMx)
library(tidyverse)
set.seed(114477)

# Parameters
ax=0.7; cx=0.2; 
ay=0.3; cy=0.6;
ex=1-ax-cx; ey=1-ay-cy;
ra=0; rc=0; re=0;
mX=0.2; mY=0.1
b12=0.1; b21=0.4
raFree=FALSE; rcFree=FALSE; reFree=FALSE;

nfam_DoCXY=10000
nfam_DoCYX=100
nfam_dis= 100
nfam=nfam_DoCXY + nfam_DoCYX + nfam_dis



# Data Generating Model: DoC (X->Y) --------------------------------------------------

### Mixture Model
# Define Variables
manifests <- c("XT1", "YT1", "XT2", "YT2")
latents <- c("Axt1", "Cxt1", "Ext1", "Ayt1", "Cyt1", "Eyt1",
             "Axt2", "Cxt2", "Ext2", "Ayt2", "Cyt2", "Eyt2", "dummyra1", "dummyra2")

# Specify basemodel
basemodel <- mxModel("model", 
                     type="RAM",
                     manifestVars=manifests,
                     latentVars=latents,
                     # Causal paths
                     mxPath(from="XT1", to=c("YT1"), free=c(TRUE), value=c(b12), arrows=1, label=c("b12")),
                     mxPath(from="YT1", to=c("XT1"), free=c(TRUE), value=c(b21), arrows=1, label=c("b21")),
                     mxPath(from="XT2", to=c("YT2"), free=c(TRUE), value=c(b12), arrows=1, label=c("b12")),
                     mxPath(from="YT2", to=c("XT2"), free=c(TRUE), value=c(b21), arrows=1, label=c("b21")),
                     # a c e paths
                     mxPath(from="Axt1", to=c("XT1"), free=c(TRUE), value=c(ax), arrows=1, label=c("a11")),
                     mxPath(from="Cxt1", to=c("XT1"), free=c(TRUE), value=c(cx), arrows=1, label=c("c11")),
                     mxPath(from="Ext1", to=c("XT1"), free=c(TRUE), value=c(ex), arrows=1, label=c("e11")),
                     mxPath(from="Ayt1", to=c("YT1"), free=c(TRUE), value=c(ay), arrows=1, label=c("a22")),
                     mxPath(from="Cyt1", to=c("YT1"), free=c(TRUE), value=c(cy), arrows=1, label=c("c22")),
                     mxPath(from="Eyt1", to=c("YT1"), free=c(TRUE), value=c(ey), arrows=1, label=c("e22")),
                     mxPath(from="Axt2", to=c("XT2"), free=c(TRUE), value=c(ax), arrows=1, label=c("a11")),
                     mxPath(from="Cxt2", to=c("XT2"), free=c(TRUE), value=c(cx), arrows=1, label=c("c11")),
                     mxPath(from="Ext2", to=c("XT2"), free=c(TRUE), value=c(ex), arrows=1, label=c("e11")),
                     mxPath(from="Ayt2", to=c("YT2"), free=c(TRUE), value=c(ay), arrows=1, label=c("a22")),
                     mxPath(from="Cyt2", to=c("YT2"), free=c(TRUE), value=c(cy), arrows=1, label=c("c22")),
                     mxPath(from="Eyt2", to=c("YT2"), free=c(TRUE), value=c(ey), arrows=1, label=c("e22")),
                     mxPath(from="Cxt1", to=c("Cxt2"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Cyt1", to=c("Cyt2"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     # ra rc re
                     mxPath(from="Axt1", to=c("Axt1", "Ayt1"), free=c(FALSE, raFree), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
                     mxPath(from="Cxt1", to=c("Cxt1", "Cyt1"), free=c(FALSE, FALSE), value=c(1.0, rc), arrows=2, label=c(NA, "rc")),
                     mxPath(from="Ext1", to=c("Ext1", "Eyt1"), free=c(FALSE, FALSE), value=c(1.0, re), arrows=2, label=c(NA, "re")),
                     mxPath(from="Axt2", to=c("Axt2", "Ayt2"), free=c(FALSE, raFree), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
                     mxPath(from="Cxt2", to=c("Cxt2", "Cyt2"), free=c(FALSE, FALSE), value=c(1.0, rc), arrows=2, label=c(NA, "rc")),
                     mxPath(from="Ext2", to=c("Ext2", "Eyt2"), free=c(FALSE, FALSE), value=c(1.0, re), arrows=2, label=c(NA, "re")),
                     mxPath(from="Cxt1", to=c("Cyt2"), free=c(FALSE), value=(rc), arrows=2, label=c("rc")),
                     mxPath(from="Cyt1", to=c("Cxt2"), free=c(FALSE), value=(rc), arrows=2, label=c("rc")),
                     # variances 
                     mxPath(from="Ayt1", to=c("Ayt1"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Cyt1", to=c("Cyt1"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Eyt1", to=c("Eyt1"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Cxt2", to=c("Cxt2"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Ext2", to=c("Ext2"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Ayt2", to=c("Ayt2"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Cyt2", to=c("Cyt2"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Eyt2", to=c("Eyt2"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     mxPath(from="Axt1", to=c("Axt1"), free=c(FALSE), value=c(1.0), arrows=2, label=NA),
                     # Means
                     mxPath(from="one", to=c("XT1", "YT1", "XT2", "YT2"), 
                            free=c(TRUE, TRUE, TRUE, TRUE), value=c(0, 0, 0, 0), arrows=1, label=c("mX", "mY", "mX", "mY"))
)

# Specify MZ/DZ models
mzBaseModel <- mxModel(basemodel,
                       mxPath(from="Axt1", to=c("Axt2"), free=FALSE, value=c(1.0), arrows=2, label=NA),
                       mxPath(from="Ayt1", to=c("Ayt2"), free=FALSE, value=c(1.0), arrows=2, label=NA),
                       # Dummy path for cross twin cross trait covariance
                       mxPath(from="Axt1", to=c("dummyra1"), free=c(FALSE), value=c(ra), arrows=2, label=c("ra")),
                       mxPath(from="dummyra1", to=c("Ayt2"), free=c(FALSE), value=c(1.0), arrows=1),
                       mxPath(from="Ayt1", to=c("dummyra2"), free=c(FALSE), value=c(ra), arrows=2, label=c("ra")),
                       mxPath(from="dummyra2", to=c("Axt2"), free=c(FALSE), value=c(1.0), arrows=1),
                       mxPath(from="dummyra1", to=c("dummyra1"), free=c(FALSE), value=c(0), arrows=2, label=c(NA)),
                       mxPath(from="dummyra2", to=c("dummyra2"), free=c(FALSE), value=c(0), arrows=2, label=c(NA)),
                       mxFitFunctionML(vector=TRUE),
                       name="MZ")

dzBaseModel <- mxModel(basemodel, 
                       mxPath(from="Axt1", to=c("Axt2"), free=FALSE, value=c(0.5), arrows=2, label=NA),
                       mxPath(from="Ayt1", to=c("Ayt2"), free=FALSE, value=c(0.5), arrows=2, label=NA),
                       # Dummy path for cross twin cross trait covariance
                       mxPath(from="Axt1", to=c("dummyra1"), free=c(FALSE), value=c(ra), arrows=2, label=c("ra")),
                       mxPath(from="dummyra1", to=c("Ayt2"), free=c(FALSE), value=c(0.5), arrows=1),
                       mxPath(from="Ayt1", to=c("dummyra2"), free=c(FALSE), value=c(ra), arrows=2, label=c("ra")),
                       mxPath(from="dummyra2", to=c("Axt2"), free=c(FALSE), value=c(0.5), arrows=1),
                       mxPath(from="dummyra1", to=c("dummyra1"), free=c(FALSE), value=c(0), arrows=2, label=c(NA)),
                       mxPath(from="dummyra2", to=c("dummyra2"), free=c(FALSE), value=c(0), arrows=2, label=c(NA)),
                       mxFitFunctionML(vector=TRUE),
                       name="DZ")

# Set free status for ra/rc/re 
mzModel <- omxSetParameters(mzBaseModel, label=c("ra", "rc", "re"), 
                            free=c(raFree, rcFree, reFree))

dzModel <- omxSetParameters(dzBaseModel, label=c("ra", "rc", "re"), 
                            free=c(raFree, rcFree, reFree))



### Generate Classes
# Concordant T1:X->Y T2:X->Y (b21t1=0  b21t2=0)
MZ_XYcon <- omxSetParameters(mzModel,  label=c("b21"), free=FALSE, value=0, newlabels=c("b21zero"), name="MZ_XYcon")
MZ_XYcon <- omxSetParameters(MZ_XYcon, label=c("b12"), free=TRUE, value=b12, newlabels=c("b12"))
MZ_XYcon <- omxSetParameters(MZ_XYcon, label=c("mX", "mY"), value=c(0, mY), newlabels=c("mX0", "mY"))

DZ_XYcon <- omxSetParameters(dzModel,  label=c("b21"), free=FALSE, value=0, newlabels=c("b21zero"), name="DZ_XYcon")
DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("b12"), free=TRUE, value=b12, newlabels=c("b12"))
DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("mX", "mY"), value=c(0, mY), newlabels=c("mX0", "mY"))


# Discordant T1:X->Y T2:X<-Y (b21t1=0  b12t2=0)
MZ_XYdis <- mzModel;
MZ_XYdis$A$values[1,2]<-0; MZ_XYdis$A$free[1,2]<-F; MZ_XYdis$A$labels[1,2] <- "b21zero"
MZ_XYdis$A$values[2,1]<-b12; MZ_XYdis$A$free[2,1]<-T; MZ_XYdis$A$labels[2,1] <- "b12"
MZ_XYdis$A$values[3,4]<-b21; MZ_XYdis$A$free[3,4]<-T; MZ_XYdis$A$labels[3,4] <- "b21"
MZ_XYdis$A$values[4,3]<-0; MZ_XYdis$A$free[4,3]<-F; MZ_XYdis$A$labels[4,3] <- "b21zero"
MZ_XYdis <- mxModel(MZ_XYdis, name="MZ_XYdis")
MZ_XYdis <- omxSetParameters(MZ_XYdis, label=c("mX", "mY"), value=c(mX, mY), newlabels=c("mX", "mY"))

DZ_XYdis <- dzModel;
DZ_XYdis$A$values[1,2]<-0; DZ_XYdis$A$free[1,2]<-F; DZ_XYdis$A$labels[1,2] <- "b21zero"
DZ_XYdis$A$values[2,1]<-b12; DZ_XYdis$A$free[2,1]<-T; DZ_XYdis$A$labels[2,1] <- "b12"
DZ_XYdis$A$values[3,4]<-b21; DZ_XYdis$A$free[3,4]<-T; DZ_XYdis$A$labels[3,4] <- "b21"
DZ_XYdis$A$values[4,3]<-0; DZ_XYdis$A$free[4,3]<-F; DZ_XYdis$A$labels[4,3] <- "b12zero"
DZ_XYdis <- mxModel(DZ_XYdis, name="DZ_XYdis")
DZ_XYdis <- omxSetParameters(DZ_XYdis, label=c("mX", "mY"), value=c(mX, mY), newlabels=c("mX", "mY"))


# Discordant T1:X<-Y T2:X->Y (b12t1=0  b21t2=0)
MZ_YXdis <- mzModel; MZ_YXdis$A$values[2,1]<-0; MZ_YXdis$A$free[2,1]<-F; MZ_YXdis$A$labels[2,1] <- "b21zero"
MZ_YXdis$A$values[1,2]<-b21; MZ_YXdis$A$free[1,2]<-T; MZ_YXdis$A$labels[1,2] <- "b21"
MZ_YXdis$A$values[4,3]<-b12; MZ_YXdis$A$free[4,3]<-T; MZ_YXdis$A$labels[4,3] <- "b12"
MZ_YXdis$A$values[3,4]<-0; MZ_YXdis$A$free[3,4]<-F; MZ_YXdis$A$labels[3,4] <- "b21zero"
MZ_YXdis <- mxModel(MZ_YXdis, name="MZ_YXdis")
MZ_YXdis <- omxSetParameters(MZ_YXdis, label=c("mX", "mY"), value=c(mX, mY), newlabels=c("mX", "mY"))

DZ_YXdis <- dzModel; DZ_YXdis$A$values[2,1]<-0; DZ_YXdis$A$free[2,1]<-F; DZ_YXdis$A$labels[2,1] <- "b21zero"
DZ_YXdis$A$values[1,2]<-b21; DZ_YXdis$A$free[1,2]<-T; DZ_YXdis$A$labels[1,2] <- "b21"
DZ_YXdis$A$values[4,3]<-b12; DZ_YXdis$A$free[4,3]<-T; DZ_YXdis$A$labels[4,3] <- "b12"
DZ_YXdis$A$values[3,4]<-0; DZ_YXdis$A$free[3,4]<-F; DZ_YXdis$A$labels[3,4] <- "b21zero"
DZ_YXdis <- mxModel(DZ_YXdis, name="DZ_YXdis")
DZ_YXdis <- omxSetParameters(DZ_YXdis, label=c("mX", "mY"), value=c(mX, mY), newlabels=c("mX", "mY"))


# Concordant T1:X<-Y T2:X<-Y (b12t1=0  b12t2=0)
MZ_YXcon <- omxSetParameters(mzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="MZ_YXcon")
MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("mX", "mY"), value=c(mX, 0), newlabels=c("mX", "mY0"))

DZ_YXcon <- omxSetParameters(dzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="DZ_YXcon")
DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("mX", "mY"), value=c(mX, 0), newlabels=c("mX", "mY0"))

# Generate Data
MZ_XYcon_dat <- mxGenerateData(MZ_XYcon, nrows=(nfam_DoCXY/2), empirical=TRUE)
DZ_XYcon_dat <- mxGenerateData(DZ_XYcon, nrows=(nfam_DoCXY/2), empirical=TRUE)

mzData <- MZ_XYcon_dat
dzData <- DZ_XYcon_dat

nmz <- dim(MZ_XYcon_dat)[1]
ndz <- dim(DZ_XYcon_dat)[1]

# Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

### Generate MZ/DZ/Mix Model
modelMZ1 <- mxModel(	"MZ1", DataMZ, 
                    MZ_XYcon, MZ_XYdis, MZ_YXdis, MZ_YXcon,
                    mxMatrix(type="Full", 1, 4, free=c(F,T,T,T),
                             values=c(1), lbound=0.001, ubound=100, 
                             labels=c("p1", "p2", "p2", "p4"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("MZ_XYcon", "MZ_XYdis", "MZ_YXdis", "MZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())


modelDZ1 <- mxModel("DZ1", DataDZ, 
                   DZ_XYcon, DZ_XYdis, DZ_YXdis, DZ_YXcon,
                   mxMatrix(type="Full", 1, 4, free=c(F,T,T,T), 
                            values=c(1), 
                            lbound=0.001, ubound=100, labels=c("p5", "p6", "p6", "p8"), name="weights"),
                   mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                             name="weightT"),
                   mxAlgebra(weights/weightT, name="stdweights"),
                   mxExpectationMixture(components=c("DZ_XYcon", "DZ_XYdis", "DZ_YXdis", "DZ_YXcon"), 
                                        weights="weights", scale='sum'),
                   mxFitFunctionML())

modelMZ2 <- mxModel(	"MZ2", DataMZ, 
                    MZ_XYcon, MZ_YXcon,
                    mxMatrix(type="Full", 1, 2, free=c(F,T),
                             values=c(1), lbound=0.001, ubound=100, 
                             labels=c("p1", "p4"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("MZ_XYcon", "MZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())


modelDZ2 <- mxModel("DZ2", DataDZ, 
                   DZ_XYcon, DZ_YXcon,
                   mxMatrix(type="Full", 1, 2, free=c(F,T), 
                            values=c(1), 
                            lbound=0.001, ubound=100, labels=c("p5", "p8"), name="weights"),
                   mxAlgebra(weights[1,1]+weights[1,2],
                             name="weightT"),
                   mxAlgebra(weights/weightT, name="stdweights"),
                   mxExpectationMixture(components=c("DZ_XYcon", "DZ_YXcon"), 
                                        weights="weights", scale='sum'),
                   mxFitFunctionML())

### Build Mixture Models
mixModel_4class <- mxModel("4class", modelMZ1, modelDZ1, mxFitFunctionMultigroup(c("MZ1", "DZ1")))
mixModel_2class <- mxModel("2class", modelMZ2, modelDZ2, mxFitFunctionMultigroup(c("MZ2", "DZ2")))



### DoC Model Specification (Univariate and bidirectional)
# Mean Matrix
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")

# Path coefficient matrices 
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  label=c("a11", NA, NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  label=c("c11", NA, NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  label=c("e11", NA, NA, "e22"), lbound=0.0001, name="e")

# Matrices for Correlation Coefficients within/across Individuals
pathRa <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=ra, 
                   label=c(NA, "ra", "ra", NA), lbound=-1, ubound=1, name="Ra")
pathRc <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=rc, 
                   label=c(NA, "rc", "rc", NA), lbound=-1, ubound=1, name="Rc")
pathRe <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=re, 
                   label=c(NA, "re", "re", NA), lbound=-1, ubound=1, name="Re")

# Causal path matrix
pathB <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(F, F, F, F), 
                  labels=c("b11", "b12", "b21", 'b22'), name="b")

# Algebra to Constrain Correlation Matrices to be Positive Definite
cnstrPos <- mxMatrix(type="Full", nrow=1, ncol=6, free=FALSE, values=.0001, name="cnstrPos" )
corMin <- mxAlgebra(expression= cbind(min(eigenval(Ram)), min(eigenval(Rcm)), min(eigenval(Rem)),
                                      min(eigenval(Raf)), min(eigenval(Rcf)), min(eigenval(Ref))), name="corMin" )
corPos <- mxConstraint(expression=corMin > cnstrPos, name="corPos")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% (Ra) %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% (Rc) %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% (Re) %*% t(e), name="E")

# Create Algebra for Causal Variance Components
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invI <- mxAlgebra(expression=solve(I-b), name="IB")

# Create Algebra for Total Variances and Standard Deviations (diagonal only)
covP <- mxAlgebra(expression=A+C+E, name="V")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")
covIP <- mxAlgebra(expression=IB %&% V, name="IV")

# Create Algebra for Expected Variance/Covariance Matrices in MZ & DZ twins
covMZ <- mxAlgebra(expression=A+C, name="cMZ")
covDZ <- mxAlgebra(expression=0.5%x%(A)+C, name="cDZ")
matI2 <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I2")
expCovMZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cMZ), cbind(cMZ, V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cDZ), cbind(cDZ, V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

### Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, pathRa, pathRc, pathRe, pathB, 
             covA, covC, covE, covP, matI, invI, invSD, covIP, matI2)
modelMZ <- mxModel(pars, DataMZ, covMZ, expCovMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, DataDZ, covDZ, expCovDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Build Models
DoC <- mxModel("DoC", pars, modelMZ, modelDZ, multi)
DoC <- omxSetParameters(DoC, labels=c("ra", "rc", "re"), free=c(raFree, rcFree, reFree), values=c(ra, rc, re))
DoCXY <- omxSetParameters(DoC, labels=c("b12"), free=TRUE, values=b12, name="DoCXY")
DoCXY <- omxSetParameters(DoCXY, labels=c("mX"), free=TRUE, values=0, newlabels=c("mX0"))
DoCYX <- omxSetParameters(DoC, labels=c("b21"), free=TRUE, values=b21, name="DoCYX")
DoCYX <- omxSetParameters(DoCYX, labels=c("mY"), free=TRUE, values=0, newlabels=c("mY0"))
DoCbi <- omxSetParameters(DoC, labels=c("b12", "b21"), free=c(TRUE, TRUE), values=c(b12, b21), name="DoCbi")


### Choleksy Model Specification 
# Create Algebra for expected Mean Matrices
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")
# Create Matrices for Path Coefficients
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  labels=c("a11", "a21", NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  labels=c("c11", "c21", NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  labels=c("e11", "e21", NA, "e22"), lbound=0.0001, name="e")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% t(e), name="E")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP <- mxAlgebra(expression= A+C+E, name="V")
covMZ <- mxAlgebra(expression= A+C, name="cMZ")
covDZ <- mxAlgebra(expression= 0.5%x%A + C, name="cDZ")
expCovMZ <- mxAlgebra(expression= rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression= rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, covA, covC, covE, covP)
modelMZ <- mxModel(pars, covMZ, expCovMZ, DataMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, covDZ, expCovDZ, DataDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Create Algebra for Standardization
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")

# Calculate genetic and environmental correlations
corA <- mxAlgebra(expression=solve(sqrt(I*A))%&%A, name="rA") 
corC <- mxAlgebra(expression=solve(sqrt(I*C))%&%C, name="rC")
corE <- mxAlgebra(expression=solve(sqrt(I*E))%&%E, name="rE")

# Create Algebra for Variance Components
rowUS <- rep('US', 2)
colUS <- rep(c('A','C','E','SA','SC','SE'), each=2)
estUS <- mxAlgebra(expression=cbind(A,C,E,A/V,C/V,E/V), name="US", dimnames=list(rowUS,colUS))

# Build Model
calc <- list( matI, invSD, corA, corC, corE, estUS)
choleskymodel <- mxModel("cholesky", pars, modelMZ, modelDZ, multi, calc)



### Run and Compare Models
# Run models
mixModel_4class_run1 <- mxRun(mixModel_4class)
mixModel_4class_run1sum <- summary(mixModel_4class_run1)
mixModel_2class_run1 <- mxRun(mixModel_2class)
mixModel_2class_run1sum <- summary(mixModel_2class_run1)
DoCXY_run1 <- mxRun(DoCXY)
DoCXY_run1sum <- summary(DoCXY_run1)
DoCYX_run1 <- mxRun(DoCYX)
DoCYX_run1sum <- summary(DoCYX_run1)
DoCbi_run1 <- mxRun(DoCbi)
DoCbi_run1sum <- summary(DoCbi_run1)
cholesky_run1 <- mxRun(choleskymodel)
cholesky_run1sum <- summary(cholesky_run1)


run1 <- matrix(nrow=6, ncol=7)
colnames(run1) <- c("model", "ep", "df", "minus2LL", "AIC", "deltaAIC", "AkaikeWeight")
run1[1,1] <- "mixModel_4class" ; run1[2,1] <- "mixModel_2class" ; run1[3,1] <- "DoCXY"; 
run1[4,1] <- "DoCYX"; run1[5,1] <- "DoCbi"; run1[6,1] <- "cholesky"
run1[1,2] <- mixModel_4class_run1sum$estimatedParameters; run1[2,2] <- mixModel_2class_run1sum$estimatedParameters ; run1[3,2] <- DoCXY_run1sum$estimatedParameters; 
run1[4,2] <- DoCYX_run1sum$estimatedParameters ; run1[5,2] <- DoCbi_run1sum$estimatedParameters; run1[6,2] <- cholesky_run1sum$estimatedParameters
run1[,3] <- nfam_DoCXY*4-as.numeric(run1[,2])
run1[1,4] <- mixModel_4class_run1sum$fit; run1[2,4] <- mixModel_2class_run1sum$fit; run1[3,4] <- DoCXY_run1sum$fit; 
run1[4,4] <- DoCYX_run1sum$fit ; run1[5,4] <- DoCbi_run1sum$fit; run1[6,4] <- cholesky_run1sum$fit
run1[1,5] <- mixModel_4class_run1sum$information[3]; run1[2,5] <- mixModel_2class_run1sum$information[3]; run1[3,5] <- DoCXY_run1sum$information[3]; 
run1[4,5] <- DoCYX_run1sum$information[3] ; run1[5,5] <- DoCbi_run1sum$information[3]; run1[6,5] <- cholesky_run1sum$information[3]
run1 <- as.data.frame(run1)
run1$minus2LL <- round(as.numeric(run1$minus2LL), digits=3)
run1$AIC <- round(as.numeric(run1$AIC), digits=3)

run1$deltaAIC <- run1$AIC-min(run1$AIC)
run1$AkaikeWeight <- exp(-0.5*run1$deltaAIC)
AICsum1 <- sum(run1$AkaikeWeight)
run1$AkaikeWeight <- run1$AkaikeWeight/AICsum1
run1$AkaikeWeight <- round(as.numeric(run1$AkaikeWeight), digits=3)



# Data Generating Model: DoC (Y->X) --------------------------------------------------

# Generate Data
MZ_YXcon_dat <- mxGenerateData(MZ_YXcon, nrows=nfam_DoCYX/2, empirical=TRUE)
DZ_YXcon_dat <- mxGenerateData(DZ_YXcon, nrows=nfam_DoCYX/2, empirical=TRUE)

mzData <- MZ_YXcon_dat
dzData <- DZ_YXcon_dat

nmz <- dim(MZ_YXcon_dat)[1]
ndz <- dim(DZ_YXcon_dat)[1]

# Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

### Generate MZ/DZ/Mix Models
modelMZ1 <- mxModel(	"MZ1", DataMZ, 
                     MZ_XYcon, MZ_XYdis, MZ_YXdis, MZ_YXcon,
                     mxMatrix(type="Full", 1, 4, free=c(F,T,T,T),
                              values=c(1), lbound=0.001, ubound=100, 
                              labels=c("p1", "p2", "p2", "p4"), name="weights"),
                     mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                               name="weightT"),
                     mxAlgebra(weights/weightT, name="stdweights"),
                     mxExpectationMixture(components=c("MZ_XYcon", "MZ_XYdis", "MZ_YXdis", "MZ_YXcon"), 
                                          weights="weights", scale='sum'),
                     mxFitFunctionML())


modelDZ1 <- mxModel("DZ1", DataDZ, 
                    DZ_XYcon, DZ_XYdis, DZ_YXdis, DZ_YXcon,
                    mxMatrix(type="Full", 1, 4, free=c(F,T,T,T), 
                             values=c(1), 
                             lbound=0.001, ubound=100, labels=c("p5", "p6", "p6", "p8"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("DZ_XYcon", "DZ_XYdis", "DZ_YXdis", "DZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())

modelMZ2 <- mxModel(	"MZ2", DataMZ, 
                     MZ_XYcon, MZ_YXcon,
                     mxMatrix(type="Full", 1, 2, free=c(F,T),
                              values=c(1), lbound=0.001, ubound=100, 
                              labels=c("p1", "p4"), name="weights"),
                     mxAlgebra(weights[1,1]+weights[1,2],
                               name="weightT"),
                     mxAlgebra(weights/weightT, name="stdweights"),
                     mxExpectationMixture(components=c("MZ_XYcon", "MZ_YXcon"), 
                                          weights="weights", scale='sum'),
                     mxFitFunctionML())


modelDZ2 <- mxModel("DZ2", DataDZ, 
                    DZ_XYcon, DZ_YXcon,
                    mxMatrix(type="Full", 1, 2, free=c(F,T), 
                             values=c(1), 
                             lbound=0.001, ubound=100, labels=c("p5", "p8"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("DZ_XYcon", "DZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())

### Build Mixture Models
mixModel_4class <- mxModel("4class", modelMZ1, modelDZ1, mxFitFunctionMultigroup(c("MZ1", "DZ1")))
mixModel_2class <- mxModel("2class", modelMZ2, modelDZ2, mxFitFunctionMultigroup(c("MZ2", "DZ2")))

# DoC Models
### DoC Model Specification (Univariate and bidirectional)
# Mean Matrix
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")

# Path coefficient matrices 
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  label=c("a11", NA, NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  label=c("c11", NA, NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  label=c("e11", NA, NA, "e22"), lbound=0.0001, name="e")

# Matrices for Correlation Coefficients within/across Individuals
pathRa <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=ra, 
                   label=c(NA, "ra", "ra", NA), lbound=-1, ubound=1, name="Ra")
pathRc <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=rc, 
                   label=c(NA, "rc", "rc", NA), lbound=-1, ubound=1, name="Rc")
pathRe <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=re, 
                   label=c(NA, "re", "re", NA), lbound=-1, ubound=1, name="Re")

# Causal path matrix
pathB <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(F, F, F, F), 
                  labels=c("b11", "b12", "b21", 'b22'), name="b")

# Algebra to Constrain Correlation Matrices to be Positive Definite
cnstrPos <- mxMatrix(type="Full", nrow=1, ncol=6, free=FALSE, values=.0001, name="cnstrPos" )
corMin <- mxAlgebra(expression= cbind(min(eigenval(Ram)), min(eigenval(Rcm)), min(eigenval(Rem)),
                                      min(eigenval(Raf)), min(eigenval(Rcf)), min(eigenval(Ref))), name="corMin" )
corPos <- mxConstraint(expression=corMin > cnstrPos, name="corPos")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% (Ra) %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% (Rc) %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% (Re) %*% t(e), name="E")

# Create Algebra for Causal Variance Components
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invI <- mxAlgebra(expression=solve(I-b), name="IB")

# Create Algebra for Total Variances and Standard Deviations (diagonal only)
covP <- mxAlgebra(expression=A+C+E, name="V")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")
covIP <- mxAlgebra(expression=IB %&% V, name="IV")

# Create Algebra for Expected Variance/Covariance Matrices in MZ & DZ twins
covMZ <- mxAlgebra(expression=A+C, name="cMZ")
covDZ <- mxAlgebra(expression=0.5%x%(A)+C, name="cDZ")
matI2 <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I2")
expCovMZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cMZ), cbind(cMZ, V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cDZ), cbind(cDZ, V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

### Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, pathRa, pathRc, pathRe, pathB, 
             covA, covC, covE, covP, matI, invI, invSD, covIP, matI2)
modelMZ <- mxModel(pars, DataMZ, covMZ, expCovMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, DataDZ, covDZ, expCovDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Build Models
DoC <- mxModel("DoC", pars, modelMZ, modelDZ, multi)
DoC <- omxSetParameters(DoC, labels=c("ra", "rc", "re"), free=c(raFree, rcFree, reFree), values=c(ra, rc, re))
DoCXY <- omxSetParameters(DoC, labels=c("b12"), free=TRUE, values=b12, name="DoCXY")
DoCXY <- omxSetParameters(DoCXY, labels=c("mX"), free=TRUE, values=0, newlabels=c("mX0"))
DoCYX <- omxSetParameters(DoC, labels=c("b21"), free=TRUE, values=b21, name="DoCYX")
DoCYX <- omxSetParameters(DoCYX, labels=c("mY"), free=TRUE, values=0, newlabels=c("mY0"))
DoCbi <- omxSetParameters(DoC, labels=c("b12", "b21"), free=c(TRUE, TRUE), values=c(b12, b21), name="DoCbi")


### Choleksy Model Specification 
# Create Algebra for expected Mean Matrices
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")
# Create Matrices for Path Coefficients
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  labels=c("a11", "a21", NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  labels=c("c11", "c21", NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  labels=c("e11", "e21", NA, "e22"), lbound=0.0001, name="e")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% t(e), name="E")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP <- mxAlgebra(expression= A+C+E, name="V")
covMZ <- mxAlgebra(expression= A+C, name="cMZ")
covDZ <- mxAlgebra(expression= 0.5%x%A + C, name="cDZ")
expCovMZ <- mxAlgebra(expression= rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression= rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, covA, covC, covE, covP)
modelMZ <- mxModel(pars, covMZ, expCovMZ, DataMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, covDZ, expCovDZ, DataDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Create Algebra for Standardization
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")

# Calculate genetic and environmental correlations
corA <- mxAlgebra(expression=solve(sqrt(I*A))%&%A, name="rA") 
corC <- mxAlgebra(expression=solve(sqrt(I*C))%&%C, name="rC")
corE <- mxAlgebra(expression=solve(sqrt(I*E))%&%E, name="rE")

# Create Algebra for Variance Components
rowUS <- rep('US', 2)
colUS <- rep(c('A','C','E','SA','SC','SE'), each=2)
estUS <- mxAlgebra(expression=cbind(A,C,E,A/V,C/V,E/V), name="US", dimnames=list(rowUS,colUS))

# Build Model
calc <- list( matI, invSD, corA, corC, corE, estUS)
choleskymodel <- mxModel("cholesky", pars, modelMZ, modelDZ, multi, calc)

### Run and Compare Models
# Run models
mixModel_4class_run2 <- mxRun(mixModel_4class)
mixModel_4class_run2sum <- summary(mixModel_4class_run2)
mixModel_2class_run2 <- mxRun(mixModel_2class)
mixModel_2class_run2sum <- summary(mixModel_2class_run2)
DoCXY_run2 <- mxRun(DoCXY)
DoCXY_run2sum <- summary(DoCXY_run2)
DoCYX_run2 <- mxRun(DoCYX)
DoCYX_run2sum <- summary(DoCYX_run2)
DoCbi_run2 <- mxRun(DoCbi)
DoCbi_run2sum <- summary(DoCbi_run2)
cholesky_run2 <- mxRun(choleskymodel)
cholesky_run2sum <- summary(cholesky_run2)

# Results
run2 <- matrix(nrow=6, ncol=7)
colnames(run2) <- c("model", "ep", "df", "minus2LL", "AIC", "deltaAIC", "AkaikeWeight")
run2[1,1] <- "mixModel_4class" ; run2[2,1] <- "mixModel_2class" ; run2[3,1] <- "DoCXY"; 
run2[4,1] <- "DoCYX"; run2[5,1] <- "DoCbi"; run2[6,1] <- "cholesky"
run2[1,2] <- mixModel_4class_run2sum$estimatedParameters; run2[2,2] <- mixModel_2class_run2sum$estimatedParameters ; run2[3,2] <- DoCXY_run2sum$estimatedParameters; 
run2[4,2] <- DoCYX_run2sum$estimatedParameters ; run2[5,2] <- DoCbi_run2sum$estimatedParameters; run2[6,2] <- cholesky_run2sum$estimatedParameters
run2[,3] <- nfam_DoCXY*4-as.numeric(run2[,2])
run2[1,4] <- mixModel_4class_run2sum$fit; run2[2,4] <- mixModel_2class_run2sum$fit; run2[3,4] <- DoCXY_run2sum$fit; 
run2[4,4] <- DoCYX_run2sum$fit ; run2[5,4] <- DoCbi_run2sum$fit; run2[6,4] <- cholesky_run2sum$fit
run2[1,5] <- mixModel_4class_run2sum$information[3]; run2[2,5] <- mixModel_2class_run2sum$information[3]; run2[3,5] <- DoCXY_run2sum$information[3]; 
run2[4,5] <- DoCYX_run2sum$information[3] ; run2[5,5] <- DoCbi_run2sum$information[3]; run2[6,5] <- cholesky_run2sum$information[3]
run2 <- as.data.frame(run2)
run2$minus2LL <- round(as.numeric(run2$minus2LL), digits=3)
run2$AIC <- round(as.numeric(run2$AIC), digits=3)

run2$deltaAIC <- run2$AIC-min(run2$AIC)
run2$AkaikeWeight <- exp(-0.5*run2$deltaAIC)
AICsum2 <- sum(run2$AkaikeWeight)
run2$AkaikeWeight <- run2$AkaikeWeight/AICsum2
run2$AkaikeWeight <- round(as.numeric(run2$AkaikeWeight), digits=3)


# Data Generating Model: DoC (X->Y) + DoC (Y->X) --------------------------------------------------

# Merge data
mzData <- rbind(MZ_XYcon_dat, MZ_YXcon_dat)
dzData <- rbind(DZ_XYcon_dat, DZ_YXcon_dat)

nmz <- dim(MZ_YXcon_dat)[1]
ndz <- dim(DZ_YXcon_dat)[1]

# Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

### Generate MZ/DZ/Mix Models
modelMZ1 <- mxModel(	"MZ1", DataMZ, 
                     MZ_XYcon, MZ_XYdis, MZ_YXdis, MZ_YXcon,
                     mxMatrix(type="Full", 1, 4, free=c(F,T,T,T),
                              values=c(1), lbound=0.001, ubound=100, 
                              labels=c("p1", "p2", "p2", "p4"), name="weights"),
                     mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                               name="weightT"),
                     mxAlgebra(weights/weightT, name="stdweights"),
                     mxExpectationMixture(components=c("MZ_XYcon", "MZ_XYdis", "MZ_YXdis", "MZ_YXcon"), 
                                          weights="weights", scale='sum'),
                     mxFitFunctionML())


modelDZ1 <- mxModel("DZ1", DataDZ, 
                    DZ_XYcon, DZ_XYdis, DZ_YXdis, DZ_YXcon,
                    mxMatrix(type="Full", 1, 4, free=c(F,T,T,T), 
                             values=c(1), 
                             lbound=0.001, ubound=100, labels=c("p5", "p6", "p6", "p8"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("DZ_XYcon", "DZ_XYdis", "DZ_YXdis", "DZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())

modelMZ2 <- mxModel(	"MZ2", DataMZ, 
                     MZ_XYcon, MZ_YXcon,
                     mxMatrix(type="Full", 1, 2, free=c(F,T),
                              values=c(1), lbound=0.001, ubound=100, 
                              labels=c("p1", "p4"), name="weights"),
                     mxAlgebra(weights[1,1]+weights[1,2],
                               name="weightT"),
                     mxAlgebra(weights/weightT, name="stdweights"),
                     mxExpectationMixture(components=c("MZ_XYcon", "MZ_YXcon"), 
                                          weights="weights", scale='sum'),
                     mxFitFunctionML())


modelDZ2 <- mxModel("DZ2", DataDZ, 
                    DZ_XYcon, DZ_YXcon,
                    mxMatrix(type="Full", 1, 2, free=c(F,T), 
                             values=c(1), 
                             lbound=0.001, ubound=100, labels=c("p5", "p8"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("DZ_XYcon", "DZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())

### Build Mixture Models
mixModel_4class <- mxModel("4class", modelMZ1, modelDZ1, mxFitFunctionMultigroup(c("MZ1", "DZ1")))
mixModel_2class <- mxModel("2class", modelMZ2, modelDZ2, mxFitFunctionMultigroup(c("MZ2", "DZ2")))

# DoC Models
### DoC Model Specification (Univariate and bidirectional)
# Mean Matrix
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")

# Path coefficient matrices 
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  label=c("a11", NA, NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  label=c("c11", NA, NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  label=c("e11", NA, NA, "e22"), lbound=0.0001, name="e")

# Matrices for Correlation Coefficients within/across Individuals
pathRa <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=ra, 
                   label=c(NA, "ra", "ra", NA), lbound=-1, ubound=1, name="Ra")
pathRc <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=rc, 
                   label=c(NA, "rc", "rc", NA), lbound=-1, ubound=1, name="Rc")
pathRe <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=re, 
                   label=c(NA, "re", "re", NA), lbound=-1, ubound=1, name="Re")

# Causal path matrix
pathB <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(F, F, F, F), 
                  labels=c("b11", "b12", "b21", 'b22'), name="b")

# Algebra to Constrain Correlation Matrices to be Positive Definite
cnstrPos <- mxMatrix(type="Full", nrow=1, ncol=6, free=FALSE, values=.0001, name="cnstrPos" )
corMin <- mxAlgebra(expression= cbind(min(eigenval(Ram)), min(eigenval(Rcm)), min(eigenval(Rem)),
                                      min(eigenval(Raf)), min(eigenval(Rcf)), min(eigenval(Ref))), name="corMin" )
corPos <- mxConstraint(expression=corMin > cnstrPos, name="corPos")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% (Ra) %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% (Rc) %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% (Re) %*% t(e), name="E")

# Create Algebra for Causal Variance Components
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invI <- mxAlgebra(expression=solve(I-b), name="IB")

# Create Algebra for Total Variances and Standard Deviations (diagonal only)
covP <- mxAlgebra(expression=A+C+E, name="V")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")
covIP <- mxAlgebra(expression=IB %&% V, name="IV")

# Create Algebra for Expected Variance/Covariance Matrices in MZ & DZ twins
covMZ <- mxAlgebra(expression=A+C, name="cMZ")
covDZ <- mxAlgebra(expression=0.5%x%(A)+C, name="cDZ")
matI2 <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I2")
expCovMZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cMZ), cbind(cMZ, V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cDZ), cbind(cDZ, V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

### Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, pathRa, pathRc, pathRe, pathB, 
             covA, covC, covE, covP, matI, invI, invSD, covIP, matI2)
modelMZ <- mxModel(pars, DataMZ, covMZ, expCovMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, DataDZ, covDZ, expCovDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Build Models
DoC <- mxModel("DoC", pars, modelMZ, modelDZ, multi)
DoC <- omxSetParameters(DoC, labels=c("ra", "rc", "re"), free=c(raFree, rcFree, reFree), values=c(ra, rc, re))
DoCXY <- omxSetParameters(DoC, labels=c("b12"), free=TRUE, values=b12, name="DoCXY")
DoCXY <- omxSetParameters(DoCXY, labels=c("mX"), free=TRUE, values=0, newlabels=c("mX0"))
DoCYX <- omxSetParameters(DoC, labels=c("b21"), free=TRUE, values=b21, name="DoCYX")
DoCYX <- omxSetParameters(DoCYX, labels=c("mY"), free=TRUE, values=0, newlabels=c("mY0"))
DoCbi <- omxSetParameters(DoC, labels=c("b12", "b21"), free=c(TRUE, TRUE), values=c(b12, b21), name="DoCbi")


### Choleksy Model Specification 
# Create Algebra for expected Mean Matrices
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")
# Create Matrices for Path Coefficients
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  labels=c("a11", "a21", NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  labels=c("c11", "c21", NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  labels=c("e11", "e21", NA, "e22"), lbound=0.0001, name="e")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% t(e), name="E")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP <- mxAlgebra(expression= A+C+E, name="V")
covMZ <- mxAlgebra(expression= A+C, name="cMZ")
covDZ <- mxAlgebra(expression= 0.5%x%A + C, name="cDZ")
expCovMZ <- mxAlgebra(expression= rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression= rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, covA, covC, covE, covP)
modelMZ <- mxModel(pars, covMZ, expCovMZ, DataMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, covDZ, expCovDZ, DataDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Create Algebra for Standardization
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")

# Calculate genetic and environmental correlations
corA <- mxAlgebra(expression=solve(sqrt(I*A))%&%A, name="rA") 
corC <- mxAlgebra(expression=solve(sqrt(I*C))%&%C, name="rC")
corE <- mxAlgebra(expression=solve(sqrt(I*E))%&%E, name="rE")

# Create Algebra for Variance Components
rowUS <- rep('US', 2)
colUS <- rep(c('A','C','E','SA','SC','SE'), each=2)
estUS <- mxAlgebra(expression=cbind(A,C,E,A/V,C/V,E/V), name="US", dimnames=list(rowUS,colUS))

# Build Model
calc <- list( matI, invSD, corA, corC, corE, estUS)
choleskymodel <- mxModel("cholesky", pars, modelMZ, modelDZ, multi, calc)




### Run and Compare Models
# Run models
mixModel_4class_run3 <- mxRun(mixModel_4class)
mixModel_4class_run3sum <- summary(mixModel_4class_run3)
mixModel_2class_run3 <- mxRun(mixModel_2class)
mixModel_2class_run3sum <- summary(mixModel_2class_run3)
DoCXY_run3 <- mxRun(DoCXY)
DoCXY_run3sum <- summary(DoCXY_run3)
DoCYX_run3 <- mxRun(DoCYX)
DoCYX_run3sum <- summary(DoCYX_run3)
DoCbi_run3 <- mxRun(DoCbi)
DoCbi_run3sum <- summary(DoCbi_run3)
cholesky_run3 <- mxRun(choleskymodel)
cholesky_run3sum <- summary(cholesky_run3)

# Results
run3 <- matrix(nrow=6, ncol=7)
colnames(run3) <- c("model", "ep", "df", "minus2LL", "AIC", "deltaAIC", "AkaikeWeight")
run3[1,1] <- "mixModel_4class" ; run3[2,1] <- "mixModel_2class" ; run3[3,1] <- "DoCXY"; 
run3[4,1] <- "DoCYX"; run3[5,1] <- "DoCbi"; run3[6,1] <- "cholesky"
run3[1,2] <- mixModel_4class_run3sum$estimatedParameters; run3[2,2] <- mixModel_2class_run3sum$estimatedParameters ; run3[3,2] <- DoCXY_run3sum$estimatedParameters; 
run3[4,2] <- DoCYX_run3sum$estimatedParameters ; run3[5,2] <- DoCbi_run3sum$estimatedParameters; run3[6,2] <- cholesky_run3sum$estimatedParameters
run3[,3] <- nfam_DoCXY*4+nfam_DoCYX*4-as.numeric(run3[,2])
run3[1,4] <- mixModel_4class_run3sum$fit; run3[2,4] <- mixModel_2class_run3sum$fit; run3[3,4] <- DoCXY_run3sum$fit; 
run3[4,4] <- DoCYX_run3sum$fit ; run3[5,4] <- DoCbi_run3sum$fit; run3[6,4] <- cholesky_run3sum$fit
run3[1,5] <- mixModel_4class_run3sum$information[3]; run3[2,5] <- mixModel_2class_run3sum$information[3]; run3[3,5] <- DoCXY_run3sum$information[3]; 
run3[4,5] <- DoCYX_run3sum$information[3] ; run3[5,5] <- DoCbi_run3sum$information[3]; run3[6,5] <- cholesky_run3sum$information[3]
run3 <- as.data.frame(run3)
run3$minus2LL <- round(as.numeric(run3$minus2LL), digits=3)
run3$AIC <- round(as.numeric(run3$AIC), digits=3)

run3$deltaAIC <- run3$AIC-min(run3$AIC)
run3$AkaikeWeight <- exp(-0.5*run3$deltaAIC)
AICsum3 <- sum(run3$AkaikeWeight)
run3$AkaikeWeight <- run3$AkaikeWeight/AICsum3
run3$AkaikeWeight <- round(as.numeric(run3$AkaikeWeight), digits=3)


# Generate disconcordant data --------------------------------------------------

# Generate Data
MZ_XYdis_dat <- mxGenerateData(MZ_XYdis, nrows=(nfam_dis/4), empirical=TRUE)
MZ_YXdis_dat <- mxGenerateData(MZ_YXdis, nrows=(nfam_dis/4), empirical=TRUE)
DZ_XYcon_dat <- mxGenerateData(DZ_XYdis, nrows=(nfam_dis/4), empirical=TRUE)
DZ_YXcon_dat <- mxGenerateData(DZ_YXdis, nrows=(nfam_dis/4), empirical=TRUE)

mzData <- rbind(mzData, MZ_XYdis_dat, MZ_YXdis_dat)
dzData <- rbind(dzData, DZ_XYcon_dat, DZ_YXcon_dat)

nmz <- dim(MZ_XYcon_dat)[1]
ndz <- dim(DZ_XYcon_dat)[1]

# Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

### Generate MZ/DZ/Mix Model
modelMZ1 <- mxModel(	"MZ1", DataMZ, 
                     MZ_XYcon, MZ_XYdis, MZ_YXdis, MZ_YXcon,
                     mxMatrix(type="Full", 1, 4, free=c(F,T,T,T),
                              values=c(1), lbound=0.001, ubound=100, 
                              labels=c("p1", "p2", "p2", "p4"), name="weights"),
                     mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                               name="weightT"),
                     mxAlgebra(weights/weightT, name="stdweights"),
                     mxExpectationMixture(components=c("MZ_XYcon", "MZ_XYdis", "MZ_YXdis", "MZ_YXcon"), 
                                          weights="weights", scale='sum'),
                     mxFitFunctionML())


modelDZ1 <- mxModel("DZ1", DataDZ, 
                    DZ_XYcon, DZ_XYdis, DZ_YXdis, DZ_YXcon,
                    mxMatrix(type="Full", 1, 4, free=c(F,T,T,T), 
                             values=c(1), 
                             lbound=0.001, ubound=100, labels=c("p5", "p6", "p6", "p8"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("DZ_XYcon", "DZ_XYdis", "DZ_YXdis", "DZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())

modelMZ2 <- mxModel(	"MZ2", DataMZ, 
                     MZ_XYcon, MZ_YXcon,
                     mxMatrix(type="Full", 1, 2, free=c(F,T),
                              values=c(1), lbound=0.001, ubound=100, 
                              labels=c("p1", "p4"), name="weights"),
                     mxAlgebra(weights[1,1]+weights[1,2],
                               name="weightT"),
                     mxAlgebra(weights/weightT, name="stdweights"),
                     mxExpectationMixture(components=c("MZ_XYcon", "MZ_YXcon"), 
                                          weights="weights", scale='sum'),
                     mxFitFunctionML())


modelDZ2 <- mxModel("DZ2", DataDZ, 
                    DZ_XYcon, DZ_YXcon,
                    mxMatrix(type="Full", 1, 2, free=c(F,T), 
                             values=c(1), 
                             lbound=0.001, ubound=100, labels=c("p5", "p8"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("DZ_XYcon", "DZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())

### Build Mixture Models
mixModel_4class <- mxModel("4class", modelMZ1, modelDZ1, mxFitFunctionMultigroup(c("MZ1", "DZ1")))
mixModel_2class <- mxModel("2class", modelMZ2, modelDZ2, mxFitFunctionMultigroup(c("MZ2", "DZ2")))

# DoC Models
### DoC Model Specification (Univariate and bidirectional)
# Mean Matrix
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")

# Path coefficient matrices 
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  label=c("a11", NA, NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  label=c("c11", NA, NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  label=c("e11", NA, NA, "e22"), lbound=0.0001, name="e")

# Matrices for Correlation Coefficients within/across Individuals
pathRa <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=ra, 
                   label=c(NA, "ra", "ra", NA), lbound=-1, ubound=1, name="Ra")
pathRc <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=rc, 
                   label=c(NA, "rc", "rc", NA), lbound=-1, ubound=1, name="Rc")
pathRe <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=re, 
                   label=c(NA, "re", "re", NA), lbound=-1, ubound=1, name="Re")

# Causal path matrix
pathB <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(F, F, F, F), 
                  labels=c("b11", "b12", "b21", 'b22'), name="b")

# Algebra to Constrain Correlation Matrices to be Positive Definite
cnstrPos <- mxMatrix(type="Full", nrow=1, ncol=6, free=FALSE, values=.0001, name="cnstrPos" )
corMin <- mxAlgebra(expression= cbind(min(eigenval(Ram)), min(eigenval(Rcm)), min(eigenval(Rem)),
                                      min(eigenval(Raf)), min(eigenval(Rcf)), min(eigenval(Ref))), name="corMin" )
corPos <- mxConstraint(expression=corMin > cnstrPos, name="corPos")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% (Ra) %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% (Rc) %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% (Re) %*% t(e), name="E")

# Create Algebra for Causal Variance Components
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invI <- mxAlgebra(expression=solve(I-b), name="IB")

# Create Algebra for Total Variances and Standard Deviations (diagonal only)
covP <- mxAlgebra(expression=A+C+E, name="V")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")
covIP <- mxAlgebra(expression=IB %&% V, name="IV")

# Create Algebra for Expected Variance/Covariance Matrices in MZ & DZ twins
covMZ <- mxAlgebra(expression=A+C, name="cMZ")
covDZ <- mxAlgebra(expression=0.5%x%(A)+C, name="cDZ")
matI2 <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I2")
expCovMZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cMZ), cbind(cMZ, V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cDZ), cbind(cDZ, V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

### Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, pathRa, pathRc, pathRe, pathB, 
             covA, covC, covE, covP, matI, invI, invSD, covIP, matI2)
modelMZ <- mxModel(pars, DataMZ, covMZ, expCovMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, DataDZ, covDZ, expCovDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Build Models
DoC <- mxModel("DoC", pars, modelMZ, modelDZ, multi)
DoC <- omxSetParameters(DoC, labels=c("ra", "rc", "re"), free=c(raFree, rcFree, reFree), values=c(ra, rc, re))
DoCXY <- omxSetParameters(DoC, labels=c("b12"), free=TRUE, values=b12, name="DoCXY")
DoCXY <- omxSetParameters(DoCXY, labels=c("mX"), free=TRUE, values=0, newlabels=c("mX0"))
DoCYX <- omxSetParameters(DoC, labels=c("b21"), free=TRUE, values=b21, name="DoCYX")
DoCYX <- omxSetParameters(DoCYX, labels=c("mY"), free=TRUE, values=0, newlabels=c("mY0"))
DoCbi <- omxSetParameters(DoC, labels=c("b12", "b21"), free=c(TRUE, TRUE), values=c(b12, b21), name="DoCbi")


### Choleksy Model Specification 
# Create Algebra for expected Mean Matrices
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")
# Create Matrices for Path Coefficients
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  labels=c("a11", "a21", NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  labels=c("c11", "c21", NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  labels=c("e11", "e21", NA, "e22"), lbound=0.0001, name="e")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% t(e), name="E")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP <- mxAlgebra(expression= A+C+E, name="V")
covMZ <- mxAlgebra(expression= A+C, name="cMZ")
covDZ <- mxAlgebra(expression= 0.5%x%A + C, name="cDZ")
expCovMZ <- mxAlgebra(expression= rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression= rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, covA, covC, covE, covP)
modelMZ <- mxModel(pars, covMZ, expCovMZ, DataMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, covDZ, expCovDZ, DataDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Create Algebra for Standardization
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")

# Calculate genetic and environmental correlations
corA <- mxAlgebra(expression=solve(sqrt(I*A))%&%A, name="rA") 
corC <- mxAlgebra(expression=solve(sqrt(I*C))%&%C, name="rC")
corE <- mxAlgebra(expression=solve(sqrt(I*E))%&%E, name="rE")

# Create Algebra for Variance Components
rowUS <- rep('US', 2)
colUS <- rep(c('A','C','E','SA','SC','SE'), each=2)
estUS <- mxAlgebra(expression=cbind(A,C,E,A/V,C/V,E/V), name="US", dimnames=list(rowUS,colUS))

# Build Model
calc <- list( matI, invSD, corA, corC, corE, estUS)
choleskymodel <- mxModel("cholesky", pars, modelMZ, modelDZ, multi, calc)

### Run and Compare Models
# Run models
mixModel_4class_run4 <- mxRun(mixModel_4class)
mixModel_4class_run4sum <- summary(mixModel_4class_run4)
mixModel_2class_run4 <- mxRun(mixModel_2class)
mixModel_2class_run4sum <- summary(mixModel_2class_run4)
DoCXY_run4 <- mxRun(DoCXY)
DoCXY_run4sum <- summary(DoCXY_run4)
DoCYX_run4 <- mxRun(DoCYX)
DoCYX_run4sum <- summary(DoCYX_run4)
DoCbi_run4 <- mxRun(DoCbi)
DoCbi_run4sum <- summary(DoCbi_run4)
cholesky_run4 <- mxRun(choleskymodel)
cholesky_run4sum <- summary(cholesky_run4)

# Results
run4 <- matrix(nrow=6, ncol=7)
colnames(run4) <- c("model", "ep", "df", "minus2LL", "AIC", "deltaAIC", "AkaikeWeight")
run4[1,1] <- "mixModel_4class" ; run4[2,1] <- "mixModel_2class" ; run4[3,1] <- "DoCXY"; 
run4[4,1] <- "DoCYX"; run4[5,1] <- "DoCbi"; run4[6,1] <- "cholesky"
run4[1,2] <- mixModel_4class_run4sum$estimatedParameters; run4[2,2] <- mixModel_2class_run4sum$estimatedParameters ; run4[3,2] <- DoCXY_run4sum$estimatedParameters; 
run4[4,2] <- DoCYX_run4sum$estimatedParameters ; run4[5,2] <- DoCbi_run4sum$estimatedParameters; run4[6,2] <- cholesky_run4sum$estimatedParameters
run4[,3] <- nfam*4-as.numeric(run4[,2])
run4[1,4] <- mixModel_4class_run4sum$fit; run4[2,4] <- mixModel_2class_run4sum$fit; run4[3,4] <- DoCXY_run4sum$fit; 
run4[4,4] <- DoCYX_run4sum$fit ; run4[5,4] <- DoCbi_run4sum$fit; run4[6,4] <- cholesky_run4sum$fit
run4[1,5] <- mixModel_4class_run4sum$information[3]; run4[2,5] <- mixModel_2class_run4sum$information[3]; run4[3,5] <- DoCXY_run4sum$information[3]; 
run4[4,5] <- DoCYX_run4sum$information[3] ; run4[5,5] <- DoCbi_run4sum$information[3]; run4[6,5] <- cholesky_run4sum$information[3]
run4 <- as.data.frame(run4)
run4$minus2LL <- round(as.numeric(run4$minus2LL), digits=3)
run4$AIC <- round(as.numeric(run4$AIC), digits=3)

run4$deltaAIC <- run4$AIC-min(run4$AIC)
run4$AkaikeWeight <- exp(-0.5*run4$deltaAIC)
AICsum4 <- sum(run4$AkaikeWeight)
run4$AkaikeWeight <- run4$AkaikeWeight/AICsum4
run4$AkaikeWeight <- round(as.numeric(run4$AkaikeWeight), digits=3)




classp <- rbind(indClassProbs(mixModel_4class_run4$DZ1, mixModel_4class_run4$DZ1$stdweights$result), 
                indClassProbs(mixModel_4class_run4$MZ1, mixModel_4class_run4$MZ1$stdweights$result))
classp <-
1-entropy(classp)






# Generate bidirectional data --------------------------------------------------

# Generate Data
mzbidata <- mxGenerateData(DoCbi$MZ, nrows=(5000), empirical=TRUE)
dzbidata <- mxGenerateData(DoCbi$DZ, nrows=(5000), empirical=TRUE)


mzData <- mzbidata
dzData <- dzbidata

nmz <- dim(MZ_XYcon_dat)[1]
ndz <- dim(DZ_XYcon_dat)[1]

# Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

### Generate MZ/DZ/Mix Model
modelMZ1 <- mxModel(	"MZ1", DataMZ, 
                     MZ_XYcon, MZ_XYdis, MZ_YXdis, MZ_YXcon,
                     mxMatrix(type="Full", 1, 4, free=c(F,T,T,T),
                              values=c(1), lbound=0.001, ubound=100, 
                              labels=c("p1", "p2", "p2", "p4"), name="weights"),
                     mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                               name="weightT"),
                     mxAlgebra(weights/weightT, name="stdweights"),
                     mxExpectationMixture(components=c("MZ_XYcon", "MZ_XYdis", "MZ_YXdis", "MZ_YXcon"), 
                                          weights="weights", scale='sum'),
                     mxFitFunctionML())


modelDZ1 <- mxModel("DZ1", DataDZ, 
                    DZ_XYcon, DZ_XYdis, DZ_YXdis, DZ_YXcon,
                    mxMatrix(type="Full", 1, 4, free=c(F,T,T,T), 
                             values=c(1), 
                             lbound=0.001, ubound=100, labels=c("p5", "p6", "p6", "p8"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2]+weights[1,3]+weights[1,4],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("DZ_XYcon", "DZ_XYdis", "DZ_YXdis", "DZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())

modelMZ2 <- mxModel(	"MZ2", DataMZ, 
                     MZ_XYcon, MZ_YXcon,
                     mxMatrix(type="Full", 1, 2, free=c(F,T),
                              values=c(1), lbound=0.001, ubound=100, 
                              labels=c("p1", "p4"), name="weights"),
                     mxAlgebra(weights[1,1]+weights[1,2],
                               name="weightT"),
                     mxAlgebra(weights/weightT, name="stdweights"),
                     mxExpectationMixture(components=c("MZ_XYcon", "MZ_YXcon"), 
                                          weights="weights", scale='sum'),
                     mxFitFunctionML())


modelDZ2 <- mxModel("DZ2", DataDZ, 
                    DZ_XYcon, DZ_YXcon,
                    mxMatrix(type="Full", 1, 2, free=c(F,T), 
                             values=c(1), 
                             lbound=0.001, ubound=100, labels=c("p5", "p8"), name="weights"),
                    mxAlgebra(weights[1,1]+weights[1,2],
                              name="weightT"),
                    mxAlgebra(weights/weightT, name="stdweights"),
                    mxExpectationMixture(components=c("DZ_XYcon", "DZ_YXcon"), 
                                         weights="weights", scale='sum'),
                    mxFitFunctionML())

### Build Mixture Models
mixModel_4class <- mxModel("4class", modelMZ1, modelDZ1, mxFitFunctionMultigroup(c("MZ1", "DZ1")))
mixModel_2class <- mxModel("2class", modelMZ2, modelDZ2, mxFitFunctionMultigroup(c("MZ2", "DZ2")))

# DoC Models
### DoC Model Specification (Univariate and bidirectional)
# Mean Matrix
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")

# Path coefficient matrices 
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  label=c("a11", NA, NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  label=c("c11", NA, NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, FALSE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  label=c("e11", NA, NA, "e22"), lbound=0.0001, name="e")

# Matrices for Correlation Coefficients within/across Individuals
pathRa <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=ra, 
                   label=c(NA, "ra", "ra", NA), lbound=-1, ubound=1, name="Ra")
pathRc <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=rc, 
                   label=c(NA, "rc", "rc", NA), lbound=-1, ubound=1, name="Rc")
pathRe <- mxMatrix(type="Stand", nrow=2, ncol=2, free=TRUE, values=re, 
                   label=c(NA, "re", "re", NA), lbound=-1, ubound=1, name="Re")

# Causal path matrix
pathB <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(F, F, F, F), 
                  labels=c("b11", "b12", "b21", 'b22'), name="b")

# Algebra to Constrain Correlation Matrices to be Positive Definite
cnstrPos <- mxMatrix(type="Full", nrow=1, ncol=6, free=FALSE, values=.0001, name="cnstrPos" )
corMin <- mxAlgebra(expression= cbind(min(eigenval(Ram)), min(eigenval(Rcm)), min(eigenval(Rem)),
                                      min(eigenval(Raf)), min(eigenval(Rcf)), min(eigenval(Ref))), name="corMin" )
corPos <- mxConstraint(expression=corMin > cnstrPos, name="corPos")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% (Ra) %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% (Rc) %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% (Re) %*% t(e), name="E")

# Create Algebra for Causal Variance Components
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invI <- mxAlgebra(expression=solve(I-b), name="IB")

# Create Algebra for Total Variances and Standard Deviations (diagonal only)
covP <- mxAlgebra(expression=A+C+E, name="V")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")
covIP <- mxAlgebra(expression=IB %&% V, name="IV")

# Create Algebra for Expected Variance/Covariance Matrices in MZ & DZ twins
covMZ <- mxAlgebra(expression=A+C, name="cMZ")
covDZ <- mxAlgebra(expression=0.5%x%(A)+C, name="cDZ")
matI2 <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I2")
expCovMZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cMZ), cbind(cMZ, V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression=I2 %x% IB %&% rbind( cbind(V, cDZ), cbind(cDZ, V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

### Create mxData Objects
DataMZ <- mxData(observed=mzData, type="raw")
DataDZ <- mxData(observed=dzData, type="raw")

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, pathRa, pathRc, pathRe, pathB, 
             covA, covC, covE, covP, matI, invI, invSD, covIP, matI2)
modelMZ <- mxModel(pars, DataMZ, covMZ, expCovMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, DataDZ, covDZ, expCovDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Build Models
DoC <- mxModel("DoC", pars, modelMZ, modelDZ, multi)
DoC <- omxSetParameters(DoC, labels=c("ra", "rc", "re"), free=c(raFree, rcFree, reFree), values=c(ra, rc, re))
DoCXY <- omxSetParameters(DoC, labels=c("b12"), free=TRUE, values=b12, name="DoCXY")
DoCXY <- omxSetParameters(DoCXY, labels=c("mX"), free=TRUE, values=0, newlabels=c("mX0"))
DoCYX <- omxSetParameters(DoC, labels=c("b21"), free=TRUE, values=b21, name="DoCYX")
DoCYX <- omxSetParameters(DoCYX, labels=c("mY"), free=TRUE, values=0, newlabels=c("mY0"))
DoCbi <- omxSetParameters(DoC, labels=c("b12", "b21"), free=c(TRUE, TRUE), values=c(b12, b21), name="DoCbi")


### Choleksy Model Specification 
# Create Algebra for expected Mean Matrices
meanG <- mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(mX, mY, mX, mY), 
                  labels=c("mX", "mY", "mX", "mY"), name="meanG")
# Create Matrices for Path Coefficients
pathA <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ax, 0, 0, ay), 
                  labels=c("a11", "a21", NA, "a22"), lbound=0.0001, name="a")
pathC <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(cx, 0, 0, cy), 
                  labels=c("c11", "c21", NA, "c22"), lbound=0.0001, name="c")
pathE <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE, TRUE, FALSE, TRUE), values=c(ex, 0, 0, ey), 
                  labels=c("e11", "e21", NA, "e22"), lbound=0.0001, name="e")

# Create Algebra for Variance Components
covA <- mxAlgebra(expression=a %*% t(a), name="A")
covC <- mxAlgebra(expression=c %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% t(e), name="E")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP <- mxAlgebra(expression= A+C+E, name="V")
covMZ <- mxAlgebra(expression= A+C, name="cMZ")
covDZ <- mxAlgebra(expression= 0.5%x%A + C, name="cDZ")
expCovMZ <- mxAlgebra(expression= rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ")
expCovDZ <- mxAlgebra(expression= rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ")

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal(covariance="expCovMZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
expDZ <- mxExpectationNormal(covariance="expCovDZ", means="meanG", dimnames=c("XT1", "YT1", "XT2", "YT2"))
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars <- list(meanG, pathA, pathC, pathE, covA, covC, covE, covP)
modelMZ <- mxModel(pars, covMZ, expCovMZ, DataMZ, expMZ, funML, name="MZ")
modelDZ <- mxModel(pars, covDZ, expCovDZ, DataDZ, expDZ, funML, name="DZ")
multi <- mxFitFunctionMultigroup(c("MZ","DZ"))

# Create Algebra for Standardization
matI <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")
invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")

# Calculate genetic and environmental correlations
corA <- mxAlgebra(expression=solve(sqrt(I*A))%&%A, name="rA") 
corC <- mxAlgebra(expression=solve(sqrt(I*C))%&%C, name="rC")
corE <- mxAlgebra(expression=solve(sqrt(I*E))%&%E, name="rE")

# Create Algebra for Variance Components
rowUS <- rep('US', 2)
colUS <- rep(c('A','C','E','SA','SC','SE'), each=2)
estUS <- mxAlgebra(expression=cbind(A,C,E,A/V,C/V,E/V), name="US", dimnames=list(rowUS,colUS))

# Build Model
calc <- list( matI, invSD, corA, corC, corE, estUS)
choleskymodel <- mxModel("cholesky", pars, modelMZ, modelDZ, multi, calc)

### Run and Compare Models
# Run models
mixModel_4class_run4 <- mxRun(mixModel_4class)
mixModel_4class_run4sum <- summary(mixModel_4class_run4)
mixModel_2class_run4 <- mxRun(mixModel_2class)
mixModel_2class_run4sum <- summary(mixModel_2class_run4)
DoCXY_run4 <- mxRun(DoCXY)
DoCXY_run4sum <- summary(DoCXY_run4)
DoCYX_run4 <- mxRun(DoCYX)
DoCYX_run4sum <- summary(DoCYX_run4)
DoCbi_run4 <- mxRun(DoCbi)
DoCbi_run4sum <- summary(DoCbi_run4)
cholesky_run4 <- mxRun(choleskymodel)
cholesky_run4sum <- summary(cholesky_run4)

