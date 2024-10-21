# Philip Vinh
# Parameter Exploration (With Mean Differences)
set.seed(114477)

# Load libraries
library(OpenMx)

# Parameters
ax=0.7; cx=0.2; 
ay=0.3; cy=0.6;
ex=1-ax-cx; ey=1-ay-cy;
ra=0.0; rc=0; re=0;
mX=0.2; mY=0.1
b12=0.1; b21=0.4
raFree=FALSE; rcFree=FALSE; reFree=FALSE;
npropMZ=c(0.25, 0.25, 0.25, 0.25); npropDZ=c(0.25, 0.25, 0.25, 0.25);

dfresults <- data.frame(NULL)

### Entropy Function -----------------------------------------------------------
# Function for individual class probabilities
indClassProbs <- function(model){
  cp <- mxEval(stdweights, model)
  cp2 <- as.vector(cp)
  cps <- diag(length(cp2))
  diag(cps) <- cp2
  subs <- model@submodels
  of <- function(num) {
    return(mxEval(objective, subs[[num]]))
  }
  rl <- sapply(1:length(names(subs)), of)
  raw <- (rl %*% cps)
  tot <- 1 / apply(raw, 1, sum)
  div <- matrix(rep(tot, length(cp2)), ncol=length(cp2))
  icp <- raw * div
  return(icp)
}



# Function for entropy
entropy <- function(classProbs){
  # this function takes a matrix of class probabilities from a
  # mixture model with people in rows and classes in columns
  # and returns the entropy statistic, as given by
  # Ramaswamy, Desarbo, Reibstein, and Robinson 1993
  n <- dim(classProbs)[1]
  k <- dim(classProbs)[2]
  e <- 1-(sum(-classProbs*log(classProbs))/(n*log(k)))
  return(e)
}


# Simulations ------------------------------------------------------------------

for(i in c(5000, 10000, 20000, 50000, 100000)){
  i=i
  nMZ=i/2; nDZ=i/2
  
  # Prepare Base Model
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
                              free=c(TRUE, TRUE, TRUE, TRUE), value=c(0, 0, 0, 0), arrows=1, label=c("mXT1", "mYT1", "mXT2", "mYT2"))
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
  MZ_XYcon <- omxSetParameters(MZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  DZ_XYcon <- omxSetParameters(dzModel,  label=c("b21"), free=FALSE, value=0, newlabels=c("b21zero"), name="DZ_XYcon")
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("b12"), free=TRUE, value=b12, newlabels=c("b12"))
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  
  # Discordant T1:X->Y T2:X<-Y (b21t1=0  b12t2=0)
  MZ_XYdis <- mzModel;
  MZ_XYdis$A$values[1,2]<-0; MZ_XYdis$A$free[1,2]<-F; MZ_XYdis$A$labels[1,2] <- "b21zero"
  MZ_XYdis$A$values[2,1]<-b12; MZ_XYdis$A$free[2,1]<-T; MZ_XYdis$A$labels[2,1] <- "b12"
  MZ_XYdis$A$values[3,4]<-b21; MZ_XYdis$A$free[3,4]<-T; MZ_XYdis$A$labels[3,4] <- "b21"
  MZ_XYdis$A$values[4,3]<-0; MZ_XYdis$A$free[4,3]<-F; MZ_XYdis$A$labels[4,3] <- "b21zero"
  MZ_XYdis <- mxModel(MZ_XYdis, name="MZ_XYdis")
  MZ_XYdis <- omxSetParameters(MZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  DZ_XYdis <- dzModel;
  DZ_XYdis$A$values[1,2]<-0; DZ_XYdis$A$free[1,2]<-F; DZ_XYdis$A$labels[1,2] <- "b21zero"
  DZ_XYdis$A$values[2,1]<-b12; DZ_XYdis$A$free[2,1]<-T; DZ_XYdis$A$labels[2,1] <- "b12"
  DZ_XYdis$A$values[3,4]<-b21; DZ_XYdis$A$free[3,4]<-T; DZ_XYdis$A$labels[3,4] <- "b21"
  DZ_XYdis$A$values[4,3]<-0; DZ_XYdis$A$free[4,3]<-F; DZ_XYdis$A$labels[4,3] <- "b12zero"
  DZ_XYdis <- mxModel(DZ_XYdis, name="DZ_XYdis")
  DZ_XYdis <- omxSetParameters(DZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  
  # Discordant T1:X<-Y T2:X->Y (b12t1=0  b21t2=0)
  MZ_YXdis <- mzModel; MZ_YXdis$A$values[2,1]<-0; MZ_YXdis$A$free[2,1]<-F; MZ_YXdis$A$labels[2,1] <- "b21zero"
  MZ_YXdis$A$values[1,2]<-b21; MZ_YXdis$A$free[1,2]<-T; MZ_YXdis$A$labels[1,2] <- "b21"
  MZ_YXdis$A$values[4,3]<-b12; MZ_YXdis$A$free[4,3]<-T; MZ_YXdis$A$labels[4,3] <- "b12"
  MZ_YXdis$A$values[3,4]<-0; MZ_YXdis$A$free[3,4]<-F; MZ_YXdis$A$labels[3,4] <- "b21zero"
  MZ_YXdis <- mxModel(MZ_YXdis, name="MZ_YXdis")
  MZ_YXdis <- omxSetParameters(MZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  DZ_YXdis <- dzModel; DZ_YXdis$A$values[2,1]<-0; DZ_YXdis$A$free[2,1]<-F; DZ_YXdis$A$labels[2,1] <- "b21zero"
  DZ_YXdis$A$values[1,2]<-b21; DZ_YXdis$A$free[1,2]<-T; DZ_YXdis$A$labels[1,2] <- "b21"
  DZ_YXdis$A$values[4,3]<-b12; DZ_YXdis$A$free[4,3]<-T; DZ_YXdis$A$labels[4,3] <- "b12"
  DZ_YXdis$A$values[3,4]<-0; DZ_YXdis$A$free[3,4]<-F; DZ_YXdis$A$labels[3,4] <- "b21zero"
  DZ_YXdis <- mxModel(DZ_YXdis, name="DZ_YXdis")
  DZ_YXdis <- omxSetParameters(DZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  
  # Concordant T1:X<-Y T2:X<-Y (b12t1=0  b12t2=0)
  MZ_YXcon <- omxSetParameters(mzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="MZ_YXcon")
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0", "mX", "mY0"))
  
  DZ_YXcon <- omxSetParameters(dzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="DZ_YXcon")
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0","mX", "mY0"))
  
  ### Generate Data
  MZ_XYcon_dat <- mxGenerateData(MZ_XYcon, nrows=nMZ*npropMZ[1], empirical=TRUE)
  DZ_XYcon_dat <- mxGenerateData(DZ_XYcon, nrows=nDZ*npropDZ[1], empirical=TRUE)
  
  MZ_XYdis_dat <- mxGenerateData(MZ_XYdis, nrows=nMZ*npropMZ[2], empirical=TRUE)
  DZ_XYdis_dat <- mxGenerateData(DZ_XYdis, nrows=nDZ*npropDZ[2], empirical=TRUE)
  
  MZ_YXdis_dat <- mxGenerateData(MZ_YXdis, nrows=nMZ*npropMZ[3], empirical=TRUE)
  DZ_YXdis_dat <- mxGenerateData(DZ_YXdis, nrows=nDZ*npropDZ[3], empirical=TRUE)
  
  MZ_YXcon_dat <- mxGenerateData(MZ_YXcon, nrows=nMZ*npropMZ[4], empirical=TRUE)
  DZ_YXcon_dat <- mxGenerateData(DZ_YXcon, nrows=nDZ*npropDZ[4], empirical=TRUE)
  
  nSumMZ <- sum(c(dim(MZ_XYcon_dat)[1],dim(MZ_XYdis_dat)[1],dim(MZ_YXdis_dat)[1],dim(MZ_YXcon_dat)[1]))
  nSumDZ <- sum(c(dim(DZ_XYcon_dat)[1],dim(DZ_XYdis_dat)[1],dim(DZ_YXdis_dat)[1],dim(DZ_YXcon_dat)[1]))
  
  ### Merge data
  mzData <- rbind(MZ_XYcon_dat, MZ_XYdis_dat, MZ_YXdis_dat, MZ_YXcon_dat)
  dzData <- rbind(DZ_XYcon_dat, DZ_XYdis_dat, DZ_YXdis_dat, DZ_YXcon_dat)
  
  ### Create mxData Objects
  DataMZ <- mxData(observed=mzData, type="raw")
  DataDZ <- mxData(observed=dzData, type="raw")
  
  
  ### Generate MZ/DZ/Mix Model
  modelMZ <- mxModel(	"MZ", DataMZ, 
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
  
  
  modelDZ <- mxModel("DZ", DataDZ, 
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
  
  
  ### Generate Mixture Model
  mixModel <- mxModel("mixture", modelMZ, modelDZ, mxFitFunctionMultigroup(c("MZ", "DZ")))
  
  # Model Run
  modelfit <- mxTryHard(mixModel, extraTries=3)
  
  # modelboot <- mxBootstrap(modelfit, replications=100)

  ### Calculate Entropy
  # Calculate individual class probabilities
  MZ_indclassprops <- indClassProbs(modelfit$MZ)
  DZ_indclassprops <- indClassProbs(modelfit$DZ)
  indclassprops_modfit <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # MZ_indclassprops <- indClassProbs(modelboot$MZ, modelboot$MZ$stdweights$result)
  # DZ_indclassprops <- indClassProbs(modelboot$DZ, modelboot$DZ$stdweights$result)
  # indclassprops_modboot <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # Calculate Entropy
  S_modfit <- entropy(indclassprops_modfit)
  # S_modboot <- entropy(indclassprops_modboot)

  
  ### Store results
  results <- c(ax, cx, ex, ay, cy, ey, ra, b12, b21, nMZ+nDZ,
               mX, mY,
               npropMZ[1], npropMZ[2], npropMZ[3], npropMZ[4], npropDZ[1], npropDZ[2], npropDZ[3], npropDZ[4],
               modelfit$MZ$stdweights$result[1], modelfit$MZ$stdweights$result[2], modelfit$MZ$stdweights$result[3], modelfit$MZ$stdweights$result[4],
               modelfit$DZ$stdweights$result[1], modelfit$DZ$stdweights$result[2], modelfit$DZ$stdweights$result[3], modelfit$DZ$stdweights$result[4],
               modelfit$output$estimate["a11"], modelfit$output$estimate["c11"], modelfit$output$estimate["e11"],
               modelfit$output$estimate["a22"], modelfit$output$estimate["c22"], modelfit$output$estimate["e22"],
               modelfit$output$estimate["ra"],
               modelfit$output$estimate["b12"], modelfit$output$estimate["b21"],
               modelfit$output$estimate["mX"], modelfit$output$estimate["mY"],
               modelfit$output$estimate["mX0"], modelfit$output$estimate["mY0"],
               S_modfit)
  dfresults <- rbind(dfresults, results)
  
  colnames(dfresults) <- c("ax", "cx", "ex", "ay", "cy", "ey", "ra", "b12", "b21", "nfam",
                           "mX", "mY",
                           "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                           "estp1", "estp2", "estp3", "estp4",
                           "estp5", "estp6", "estp7", "estp8",
                           "estax", "estca", "estex",
                           "estay", "estcy", "estey",
                           "estra",
                           "estb12", "estb21",
                           "estmXc1", "estmYc1", "estmX0", "estmY0",
                           "entropyfit")
  
  
}



# Mean Difference
for(i in c(0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)){
  nMZ=2500; nDZ=2500
  mX=i; mY=0.1
  
  # Prepare Base Model
  manifests <- c("XT1", "YT1", "XT2", "YT2")
  latents <- c("Axt1", "Cxt1", "Ext1", "Ayt1", "Cyt1", "Eyt1",
               "Axt2", "Cxt2", "Ext2", "Ayt2", "Cyt2", "Eyt2", "dummyra1", "dummyra2")
  
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
                       mxPath(from="Axt1", to=c("Axt1", "Ayt1"), free=c(FALSE, FALSE), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
                       mxPath(from="Cxt1", to=c("Cxt1", "Cyt1"), free=c(FALSE, FALSE), value=c(1.0, rc), arrows=2, label=c(NA, "rc")),
                       mxPath(from="Ext1", to=c("Ext1", "Eyt1"), free=c(FALSE, FALSE), value=c(1.0, re), arrows=2, label=c(NA, "re")),
                       mxPath(from="Axt2", to=c("Axt2", "Ayt2"), free=c(FALSE, FALSE), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
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
                              free=c(TRUE, TRUE, TRUE, TRUE), value=c(0, 0, 0, 0), arrows=1, label=c("mXT1", "mYT1", "mXT2", "mYT2"))
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
  MZ_XYcon <- omxSetParameters(MZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  DZ_XYcon <- omxSetParameters(dzModel,  label=c("b21"), free=FALSE, value=0, newlabels=c("b21zero"), name="DZ_XYcon")
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("b12"), free=TRUE, value=b12, newlabels=c("b12"))
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  
  # Discordant T1:X->Y T2:X<-Y (b21t1=0  b12t2=0)
  MZ_XYdis <- mzModel;
  MZ_XYdis$A$values[1,2]<-0; MZ_XYdis$A$free[1,2]<-F; MZ_XYdis$A$labels[1,2] <- "b21zero"
  MZ_XYdis$A$values[2,1]<-b12; MZ_XYdis$A$free[2,1]<-T; MZ_XYdis$A$labels[2,1] <- "b12"
  MZ_XYdis$A$values[3,4]<-b21; MZ_XYdis$A$free[3,4]<-T; MZ_XYdis$A$labels[3,4] <- "b21"
  MZ_XYdis$A$values[4,3]<-0; MZ_XYdis$A$free[4,3]<-F; MZ_XYdis$A$labels[4,3] <- "b21zero"
  MZ_XYdis <- mxModel(MZ_XYdis, name="MZ_XYdis")
  MZ_XYdis <- omxSetParameters(MZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  DZ_XYdis <- dzModel;
  DZ_XYdis$A$values[1,2]<-0; DZ_XYdis$A$free[1,2]<-F; DZ_XYdis$A$labels[1,2] <- "b21zero"
  DZ_XYdis$A$values[2,1]<-b12; DZ_XYdis$A$free[2,1]<-T; DZ_XYdis$A$labels[2,1] <- "b12"
  DZ_XYdis$A$values[3,4]<-b21; DZ_XYdis$A$free[3,4]<-T; DZ_XYdis$A$labels[3,4] <- "b21"
  DZ_XYdis$A$values[4,3]<-0; DZ_XYdis$A$free[4,3]<-F; DZ_XYdis$A$labels[4,3] <- "b12zero"
  DZ_XYdis <- mxModel(DZ_XYdis, name="DZ_XYdis")
  DZ_XYdis <- omxSetParameters(DZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  
  # Discordant T1:X<-Y T2:X->Y (b12t1=0  b21t2=0)
  MZ_YXdis <- mzModel; MZ_YXdis$A$values[2,1]<-0; MZ_YXdis$A$free[2,1]<-F; MZ_YXdis$A$labels[2,1] <- "b21zero"
  MZ_YXdis$A$values[1,2]<-b21; MZ_YXdis$A$free[1,2]<-T; MZ_YXdis$A$labels[1,2] <- "b21"
  MZ_YXdis$A$values[4,3]<-b12; MZ_YXdis$A$free[4,3]<-T; MZ_YXdis$A$labels[4,3] <- "b12"
  MZ_YXdis$A$values[3,4]<-0; MZ_YXdis$A$free[3,4]<-F; MZ_YXdis$A$labels[3,4] <- "b21zero"
  MZ_YXdis <- mxModel(MZ_YXdis, name="MZ_YXdis")
  MZ_YXdis <- omxSetParameters(MZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  DZ_YXdis <- dzModel; DZ_YXdis$A$values[2,1]<-0; DZ_YXdis$A$free[2,1]<-F; DZ_YXdis$A$labels[2,1] <- "b21zero"
  DZ_YXdis$A$values[1,2]<-b21; DZ_YXdis$A$free[1,2]<-T; DZ_YXdis$A$labels[1,2] <- "b21"
  DZ_YXdis$A$values[4,3]<-b12; DZ_YXdis$A$free[4,3]<-T; DZ_YXdis$A$labels[4,3] <- "b12"
  DZ_YXdis$A$values[3,4]<-0; DZ_YXdis$A$free[3,4]<-F; DZ_YXdis$A$labels[3,4] <- "b21zero"
  DZ_YXdis <- mxModel(DZ_YXdis, name="DZ_YXdis")
  DZ_YXdis <- omxSetParameters(DZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  
  # Concordant T1:X<-Y T2:X<-Y (b12t1=0  b12t2=0)
  MZ_YXcon <- omxSetParameters(mzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="MZ_YXcon")
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0", "mX", "mY0"))
  
  DZ_YXcon <- omxSetParameters(dzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="DZ_YXcon")
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0","mX", "mY0"))
  
  ### Generate Data
  MZ_XYcon_dat <- mxGenerateData(MZ_XYcon, nrows=nMZ*npropMZ[1], empirical=TRUE)
  DZ_XYcon_dat <- mxGenerateData(DZ_XYcon, nrows=nDZ*npropDZ[1], empirical=TRUE)
  
  MZ_XYdis_dat <- mxGenerateData(MZ_XYdis, nrows=nMZ*npropMZ[2], empirical=TRUE)
  DZ_XYdis_dat <- mxGenerateData(DZ_XYdis, nrows=nDZ*npropDZ[2], empirical=TRUE)
  
  MZ_YXdis_dat <- mxGenerateData(MZ_YXdis, nrows=nMZ*npropMZ[3], empirical=TRUE)
  DZ_YXdis_dat <- mxGenerateData(DZ_YXdis, nrows=nDZ*npropDZ[3], empirical=TRUE)
  
  MZ_YXcon_dat <- mxGenerateData(MZ_YXcon, nrows=nMZ*npropMZ[4], empirical=TRUE)
  DZ_YXcon_dat <- mxGenerateData(DZ_YXcon, nrows=nDZ*npropDZ[4], empirical=TRUE)
  
  nSumMZ <- sum(c(dim(MZ_XYcon_dat)[1],dim(MZ_XYdis_dat)[1],dim(MZ_YXdis_dat)[1],dim(MZ_YXcon_dat)[1]))
  nSumDZ <- sum(c(dim(DZ_XYcon_dat)[1],dim(DZ_XYdis_dat)[1],dim(DZ_YXdis_dat)[1],dim(DZ_YXcon_dat)[1]))
  
  ### Merge data
  mzData <- rbind(MZ_XYcon_dat, MZ_XYdis_dat, MZ_YXdis_dat, MZ_YXcon_dat)
  dzData <- rbind(DZ_XYcon_dat, DZ_XYdis_dat, DZ_YXdis_dat, DZ_YXcon_dat)
  
  ### Create mxData Objects
  DataMZ <- mxData(observed=mzData, type="raw")
  DataDZ <- mxData(observed=dzData, type="raw")
  
  
  ### Generate MZ/DZ/Mix Model
  modelMZ <- mxModel(	"MZ", DataMZ, 
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
  
  
  modelDZ <- mxModel("DZ", DataDZ, 
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
  
  
  ### Generate Mixture Model
  mixModel <- mxModel("mixture", modelMZ, modelDZ, mxFitFunctionMultigroup(c("MZ", "DZ")))
  
  # Model Run
  modelfit <- mxTryHard(mixModel, extraTries=3)
  
  # modelboot <- mxBootstrap(modelfit, replications=20)
  
  ### Calculate Entropy
  # Calculate individual class probabilities
  MZ_indclassprops <- indClassProbs(modelfit$MZ)
  DZ_indclassprops <- indClassProbs(modelfit$DZ)
  indclassprops_modfit <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # MZ_indclassprops <- indClassProbs(modelboot$MZ, modelboot$MZ$stdweights$result)
  # DZ_indclassprops <- indClassProbs(modelboot$DZ, modelboot$DZ$stdweights$result)
  # indclassprops_modboot <- rbind(MZ_indclassprops, DZ_indclassprops)

  # Calculate Entropy
  S_modfit <- entropy(indclassprops_modfit)
  # S_modboot <- entropy(indclassprops_modboot)
  
  
  ### Store results
  results <- c(ax, cx, ex, ay, cy, ey, ra, b12, b21, nMZ+nDZ,
               mX, mY,
               npropMZ[1], npropMZ[2], npropMZ[3], npropMZ[4], npropDZ[1], npropDZ[2], npropDZ[3], npropDZ[4],
               modelfit$MZ$stdweights$result[1], modelfit$MZ$stdweights$result[2], modelfit$MZ$stdweights$result[3], modelfit$MZ$stdweights$result[4],
               modelfit$DZ$stdweights$result[1], modelfit$DZ$stdweights$result[2], modelfit$DZ$stdweights$result[3], modelfit$DZ$stdweights$result[4],
               modelfit$output$estimate["a11"], modelfit$output$estimate["c11"], modelfit$output$estimate["e11"],
               modelfit$output$estimate["a22"], modelfit$output$estimate["c22"], modelfit$output$estimate["e22"],
               modelfit$output$estimate["ra"],
               modelfit$output$estimate["b12"], modelfit$output$estimate["b21"],
               modelfit$output$estimate["mX"], modelfit$output$estimate["mY"],
               modelfit$output$estimate["mX0"], modelfit$output$estimate["mY0"],
               S_modfit)
  dfresults <- rbind(dfresults, results)
  
  colnames(dfresults) <- c("ax", "cx", "ex", "ay", "cy", "ey", "ra", "b12", "b21", "nfam",
                           "mX", "mY",
                           "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                           "estp1", "estp2", "estp3", "estp4",
                           "estp5", "estp6", "estp7", "estp8",
                           "estax", "estca", "estex",
                           "estay", "estcy", "estey",
                           "estra",
                           "estb12", "estb21",
                           "estmXc1", "estmYc1", "estmX0", "estmY0",
                           "entropyfit")
  
  
}  
 


# Casual Effect size
for(i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)){
  nMZ=2500; nDZ=2500
  mX=0.2; mY=0.1
  b12=0.1; b21=i
  
  # Prepare Base Model
  manifests <- c("XT1", "YT1", "XT2", "YT2")
  latents <- c("Axt1", "Cxt1", "Ext1", "Ayt1", "Cyt1", "Eyt1",
               "Axt2", "Cxt2", "Ext2", "Ayt2", "Cyt2", "Eyt2", "dummyra1", "dummyra2")
  
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
                       mxPath(from="Axt1", to=c("Axt1", "Ayt1"), free=c(FALSE, FALSE), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
                       mxPath(from="Cxt1", to=c("Cxt1", "Cyt1"), free=c(FALSE, FALSE), value=c(1.0, rc), arrows=2, label=c(NA, "rc")),
                       mxPath(from="Ext1", to=c("Ext1", "Eyt1"), free=c(FALSE, FALSE), value=c(1.0, re), arrows=2, label=c(NA, "re")),
                       mxPath(from="Axt2", to=c("Axt2", "Ayt2"), free=c(FALSE, FALSE), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
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
                              free=c(TRUE, TRUE, TRUE, TRUE), value=c(0, 0, 0, 0), arrows=1, label=c("mXT1", "mYT1", "mXT2", "mYT2"))
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
  MZ_XYcon <- omxSetParameters(MZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  DZ_XYcon <- omxSetParameters(dzModel,  label=c("b21"), free=FALSE, value=0, newlabels=c("b21zero"), name="DZ_XYcon")
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("b12"), free=TRUE, value=b12, newlabels=c("b12"))
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  
  # Discordant T1:X->Y T2:X<-Y (b21t1=0  b12t2=0)
  MZ_XYdis <- mzModel;
  MZ_XYdis$A$values[1,2]<-0; MZ_XYdis$A$free[1,2]<-F; MZ_XYdis$A$labels[1,2] <- "b21zero"
  MZ_XYdis$A$values[2,1]<-b12; MZ_XYdis$A$free[2,1]<-T; MZ_XYdis$A$labels[2,1] <- "b12"
  MZ_XYdis$A$values[3,4]<-b21; MZ_XYdis$A$free[3,4]<-T; MZ_XYdis$A$labels[3,4] <- "b21"
  MZ_XYdis$A$values[4,3]<-0; MZ_XYdis$A$free[4,3]<-F; MZ_XYdis$A$labels[4,3] <- "b21zero"
  MZ_XYdis <- mxModel(MZ_XYdis, name="MZ_XYdis")
  MZ_XYdis <- omxSetParameters(MZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  DZ_XYdis <- dzModel;
  DZ_XYdis$A$values[1,2]<-0; DZ_XYdis$A$free[1,2]<-F; DZ_XYdis$A$labels[1,2] <- "b21zero"
  DZ_XYdis$A$values[2,1]<-b12; DZ_XYdis$A$free[2,1]<-T; DZ_XYdis$A$labels[2,1] <- "b12"
  DZ_XYdis$A$values[3,4]<-b21; DZ_XYdis$A$free[3,4]<-T; DZ_XYdis$A$labels[3,4] <- "b21"
  DZ_XYdis$A$values[4,3]<-0; DZ_XYdis$A$free[4,3]<-F; DZ_XYdis$A$labels[4,3] <- "b12zero"
  DZ_XYdis <- mxModel(DZ_XYdis, name="DZ_XYdis")
  DZ_XYdis <- omxSetParameters(DZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  
  # Discordant T1:X<-Y T2:X->Y (b12t1=0  b21t2=0)
  MZ_YXdis <- mzModel; MZ_YXdis$A$values[2,1]<-0; MZ_YXdis$A$free[2,1]<-F; MZ_YXdis$A$labels[2,1] <- "b21zero"
  MZ_YXdis$A$values[1,2]<-b21; MZ_YXdis$A$free[1,2]<-T; MZ_YXdis$A$labels[1,2] <- "b21"
  MZ_YXdis$A$values[4,3]<-b12; MZ_YXdis$A$free[4,3]<-T; MZ_YXdis$A$labels[4,3] <- "b12"
  MZ_YXdis$A$values[3,4]<-0; MZ_YXdis$A$free[3,4]<-F; MZ_YXdis$A$labels[3,4] <- "b21zero"
  MZ_YXdis <- mxModel(MZ_YXdis, name="MZ_YXdis")
  MZ_YXdis <- omxSetParameters(MZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  DZ_YXdis <- dzModel; DZ_YXdis$A$values[2,1]<-0; DZ_YXdis$A$free[2,1]<-F; DZ_YXdis$A$labels[2,1] <- "b21zero"
  DZ_YXdis$A$values[1,2]<-b21; DZ_YXdis$A$free[1,2]<-T; DZ_YXdis$A$labels[1,2] <- "b21"
  DZ_YXdis$A$values[4,3]<-b12; DZ_YXdis$A$free[4,3]<-T; DZ_YXdis$A$labels[4,3] <- "b12"
  DZ_YXdis$A$values[3,4]<-0; DZ_YXdis$A$free[3,4]<-F; DZ_YXdis$A$labels[3,4] <- "b21zero"
  DZ_YXdis <- mxModel(DZ_YXdis, name="DZ_YXdis")
  DZ_YXdis <- omxSetParameters(DZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  
  # Concordant T1:X<-Y T2:X<-Y (b12t1=0  b12t2=0)
  MZ_YXcon <- omxSetParameters(mzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="MZ_YXcon")
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0", "mX", "mY0"))
  
  DZ_YXcon <- omxSetParameters(dzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="DZ_YXcon")
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0","mX", "mY0"))
  
  ### Generate Data
  MZ_XYcon_dat <- mxGenerateData(MZ_XYcon, nrows=nMZ*npropMZ[1], empirical=TRUE)
  DZ_XYcon_dat <- mxGenerateData(DZ_XYcon, nrows=nDZ*npropDZ[1], empirical=TRUE)
  
  MZ_XYdis_dat <- mxGenerateData(MZ_XYdis, nrows=nMZ*npropMZ[2], empirical=TRUE)
  DZ_XYdis_dat <- mxGenerateData(DZ_XYdis, nrows=nDZ*npropDZ[2], empirical=TRUE)
  
  MZ_YXdis_dat <- mxGenerateData(MZ_YXdis, nrows=nMZ*npropMZ[3], empirical=TRUE)
  DZ_YXdis_dat <- mxGenerateData(DZ_YXdis, nrows=nDZ*npropDZ[3], empirical=TRUE)
  
  MZ_YXcon_dat <- mxGenerateData(MZ_YXcon, nrows=nMZ*npropMZ[4], empirical=TRUE)
  DZ_YXcon_dat <- mxGenerateData(DZ_YXcon, nrows=nDZ*npropDZ[4], empirical=TRUE)
  
  nSumMZ <- sum(c(dim(MZ_XYcon_dat)[1],dim(MZ_XYdis_dat)[1],dim(MZ_YXdis_dat)[1],dim(MZ_YXcon_dat)[1]))
  nSumDZ <- sum(c(dim(DZ_XYcon_dat)[1],dim(DZ_XYdis_dat)[1],dim(DZ_YXdis_dat)[1],dim(DZ_YXcon_dat)[1]))
  
  ### Merge data
  mzData <- rbind(MZ_XYcon_dat, MZ_XYdis_dat, MZ_YXdis_dat, MZ_YXcon_dat)
  dzData <- rbind(DZ_XYcon_dat, DZ_XYdis_dat, DZ_YXdis_dat, DZ_YXcon_dat)
  
  ### Create mxData Objects
  DataMZ <- mxData(observed=mzData, type="raw")
  DataDZ <- mxData(observed=dzData, type="raw")
  
  
  ### Generate MZ/DZ/Mix Model
  modelMZ <- mxModel(	"MZ", DataMZ, 
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
  
  
  modelDZ <- mxModel("DZ", DataDZ, 
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
  
  
  ### Generate Mixture Model
  mixModel <- mxModel("mixture", modelMZ, modelDZ, mxFitFunctionMultigroup(c("MZ", "DZ")))
  
  # Model Run
  modelfit <- mxTryHard(mixModel, extraTries=3)
  
  # modelboot <- mxBootstrap(modelfit, replications=100)
  
  ### Calculate Entropy
  # Calculate individual class probabilities
  MZ_indclassprops <- indClassProbs(modelfit$MZ)
  DZ_indclassprops <- indClassProbs(modelfit$DZ)
  indclassprops_modfit <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # MZ_indclassprops <- indClassProbs(modelboot$MZ, modelboot$MZ$stdweights$result)
  # DZ_indclassprops <- indClassProbs(modelboot$DZ, modelboot$DZ$stdweights$result)
  # indclassprops_modboot <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # Calculate Entropy
  S_modfit <- entropy(indclassprops_modfit)
  # S_modboot <- entropy(indclassprops_modboot)
  
  
  ### Store results
  results <- c(ax, cx, ex, ay, cy, ey, ra, b12, b21, nMZ+nDZ,
               mX, mY,
               npropMZ[1], npropMZ[2], npropMZ[3], npropMZ[4], npropDZ[1], npropDZ[2], npropDZ[3], npropDZ[4],
               modelfit$MZ$stdweights$result[1], modelfit$MZ$stdweights$result[2], modelfit$MZ$stdweights$result[3], modelfit$MZ$stdweights$result[4],
               modelfit$DZ$stdweights$result[1], modelfit$DZ$stdweights$result[2], modelfit$DZ$stdweights$result[3], modelfit$DZ$stdweights$result[4],
               modelfit$output$estimate["a11"], modelfit$output$estimate["c11"], modelfit$output$estimate["e11"],
               modelfit$output$estimate["a22"], modelfit$output$estimate["c22"], modelfit$output$estimate["e22"],
               modelfit$output$estimate["ra"],
               modelfit$output$estimate["b12"], modelfit$output$estimate["b21"],
               modelfit$output$estimate["mX"], modelfit$output$estimate["mY"],
               modelfit$output$estimate["mX0"], modelfit$output$estimate["mY0"],
               S_modfit)
  dfresults <- rbind(dfresults, results)
  
  colnames(dfresults) <- c("ax", "cx", "ex", "ay", "cy", "ey", "ra", "b12", "b21", "nfam",
                           "mX", "mY",
                           "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                           "estp1", "estp2", "estp3", "estp4",
                           "estp5", "estp6", "estp7", "estp8",
                           "estax", "estca", "estex",
                           "estay", "estcy", "estey",
                           "estra",
                           "estb12", "estb21",
                           "estmXc1", "estmYc1", "estmX0", "estmY0",
                           "entropyfit")
  
  
}  

# Variance components
for(i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)){
  nMZ=2500; nDZ=2500
  mX=0.2; mY=0.1
  b12=0.1; b21=0.4
  ay=0+i ; cy=0.9-i
  
  # Prepare Base Model
  manifests <- c("XT1", "YT1", "XT2", "YT2")
  latents <- c("Axt1", "Cxt1", "Ext1", "Ayt1", "Cyt1", "Eyt1",
               "Axt2", "Cxt2", "Ext2", "Ayt2", "Cyt2", "Eyt2", "dummyra1", "dummyra2")
  
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
                       mxPath(from="Axt1", to=c("Axt1", "Ayt1"), free=c(FALSE, FALSE), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
                       mxPath(from="Cxt1", to=c("Cxt1", "Cyt1"), free=c(FALSE, FALSE), value=c(1.0, rc), arrows=2, label=c(NA, "rc")),
                       mxPath(from="Ext1", to=c("Ext1", "Eyt1"), free=c(FALSE, FALSE), value=c(1.0, re), arrows=2, label=c(NA, "re")),
                       mxPath(from="Axt2", to=c("Axt2", "Ayt2"), free=c(FALSE, FALSE), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
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
                              free=c(TRUE, TRUE, TRUE, TRUE), value=c(0, 0, 0, 0), arrows=1, label=c("mXT1", "mYT1", "mXT2", "mYT2"))
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
  MZ_XYcon <- omxSetParameters(MZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  DZ_XYcon <- omxSetParameters(dzModel,  label=c("b21"), free=FALSE, value=0, newlabels=c("b21zero"), name="DZ_XYcon")
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("b12"), free=TRUE, value=b12, newlabels=c("b12"))
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  
  # Discordant T1:X->Y T2:X<-Y (b21t1=0  b12t2=0)
  MZ_XYdis <- mzModel;
  MZ_XYdis$A$values[1,2]<-0; MZ_XYdis$A$free[1,2]<-F; MZ_XYdis$A$labels[1,2] <- "b21zero"
  MZ_XYdis$A$values[2,1]<-b12; MZ_XYdis$A$free[2,1]<-T; MZ_XYdis$A$labels[2,1] <- "b12"
  MZ_XYdis$A$values[3,4]<-b21; MZ_XYdis$A$free[3,4]<-T; MZ_XYdis$A$labels[3,4] <- "b21"
  MZ_XYdis$A$values[4,3]<-0; MZ_XYdis$A$free[4,3]<-F; MZ_XYdis$A$labels[4,3] <- "b21zero"
  MZ_XYdis <- mxModel(MZ_XYdis, name="MZ_XYdis")
  MZ_XYdis <- omxSetParameters(MZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  DZ_XYdis <- dzModel;
  DZ_XYdis$A$values[1,2]<-0; DZ_XYdis$A$free[1,2]<-F; DZ_XYdis$A$labels[1,2] <- "b21zero"
  DZ_XYdis$A$values[2,1]<-b12; DZ_XYdis$A$free[2,1]<-T; DZ_XYdis$A$labels[2,1] <- "b12"
  DZ_XYdis$A$values[3,4]<-b21; DZ_XYdis$A$free[3,4]<-T; DZ_XYdis$A$labels[3,4] <- "b21"
  DZ_XYdis$A$values[4,3]<-0; DZ_XYdis$A$free[4,3]<-F; DZ_XYdis$A$labels[4,3] <- "b12zero"
  DZ_XYdis <- mxModel(DZ_XYdis, name="DZ_XYdis")
  DZ_XYdis <- omxSetParameters(DZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  
  # Discordant T1:X<-Y T2:X->Y (b12t1=0  b21t2=0)
  MZ_YXdis <- mzModel; MZ_YXdis$A$values[2,1]<-0; MZ_YXdis$A$free[2,1]<-F; MZ_YXdis$A$labels[2,1] <- "b21zero"
  MZ_YXdis$A$values[1,2]<-b21; MZ_YXdis$A$free[1,2]<-T; MZ_YXdis$A$labels[1,2] <- "b21"
  MZ_YXdis$A$values[4,3]<-b12; MZ_YXdis$A$free[4,3]<-T; MZ_YXdis$A$labels[4,3] <- "b12"
  MZ_YXdis$A$values[3,4]<-0; MZ_YXdis$A$free[3,4]<-F; MZ_YXdis$A$labels[3,4] <- "b21zero"
  MZ_YXdis <- mxModel(MZ_YXdis, name="MZ_YXdis")
  MZ_YXdis <- omxSetParameters(MZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  DZ_YXdis <- dzModel; DZ_YXdis$A$values[2,1]<-0; DZ_YXdis$A$free[2,1]<-F; DZ_YXdis$A$labels[2,1] <- "b21zero"
  DZ_YXdis$A$values[1,2]<-b21; DZ_YXdis$A$free[1,2]<-T; DZ_YXdis$A$labels[1,2] <- "b21"
  DZ_YXdis$A$values[4,3]<-b12; DZ_YXdis$A$free[4,3]<-T; DZ_YXdis$A$labels[4,3] <- "b12"
  DZ_YXdis$A$values[3,4]<-0; DZ_YXdis$A$free[3,4]<-F; DZ_YXdis$A$labels[3,4] <- "b21zero"
  DZ_YXdis <- mxModel(DZ_YXdis, name="DZ_YXdis")
  DZ_YXdis <- omxSetParameters(DZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  
  # Concordant T1:X<-Y T2:X<-Y (b12t1=0  b12t2=0)
  MZ_YXcon <- omxSetParameters(mzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="MZ_YXcon")
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0", "mX", "mY0"))
  
  DZ_YXcon <- omxSetParameters(dzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="DZ_YXcon")
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0","mX", "mY0"))
  
  ### Generate Data
  MZ_XYcon_dat <- mxGenerateData(MZ_XYcon, nrows=nMZ*npropMZ[1], empirical=TRUE)
  DZ_XYcon_dat <- mxGenerateData(DZ_XYcon, nrows=nDZ*npropDZ[1], empirical=TRUE)
  
  MZ_XYdis_dat <- mxGenerateData(MZ_XYdis, nrows=nMZ*npropMZ[2], empirical=TRUE)
  DZ_XYdis_dat <- mxGenerateData(DZ_XYdis, nrows=nDZ*npropDZ[2], empirical=TRUE)
  
  MZ_YXdis_dat <- mxGenerateData(MZ_YXdis, nrows=nMZ*npropMZ[3], empirical=TRUE)
  DZ_YXdis_dat <- mxGenerateData(DZ_YXdis, nrows=nDZ*npropDZ[3], empirical=TRUE)
  
  MZ_YXcon_dat <- mxGenerateData(MZ_YXcon, nrows=nMZ*npropMZ[4], empirical=TRUE)
  DZ_YXcon_dat <- mxGenerateData(DZ_YXcon, nrows=nDZ*npropDZ[4], empirical=TRUE)
  
  nSumMZ <- sum(c(dim(MZ_XYcon_dat)[1],dim(MZ_XYdis_dat)[1],dim(MZ_YXdis_dat)[1],dim(MZ_YXcon_dat)[1]))
  nSumDZ <- sum(c(dim(DZ_XYcon_dat)[1],dim(DZ_XYdis_dat)[1],dim(DZ_YXdis_dat)[1],dim(DZ_YXcon_dat)[1]))
  
  ### Merge data
  mzData <- rbind(MZ_XYcon_dat, MZ_XYdis_dat, MZ_YXdis_dat, MZ_YXcon_dat)
  dzData <- rbind(DZ_XYcon_dat, DZ_XYdis_dat, DZ_YXdis_dat, DZ_YXcon_dat)
  
  ### Create mxData Objects
  DataMZ <- mxData(observed=mzData, type="raw")
  DataDZ <- mxData(observed=dzData, type="raw")
  
  
  ### Generate MZ/DZ/Mix Model
  modelMZ <- mxModel(	"MZ", DataMZ, 
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
  
  
  modelDZ <- mxModel("DZ", DataDZ, 
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
  
  
  ### Generate Mixture Model
  mixModel <- mxModel("mixture", modelMZ, modelDZ, mxFitFunctionMultigroup(c("MZ", "DZ")))
  
  # Model Run
  modelfit <- mxTryHard(mixModel, extraTries=3)
  
  # modelboot <- mxBootstrap(modelfit, replications=100)
  
  ### Calculate Entropy
  # Calculate individual class probabilities
  MZ_indclassprops <- indClassProbs(modelfit$MZ)
  DZ_indclassprops <- indClassProbs(modelfit$DZ)
  indclassprops_modfit <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # MZ_indclassprops <- indClassProbs(modelboot$MZ, modelboot$MZ$stdweights$result)
  # DZ_indclassprops <- indClassProbs(modelboot$DZ, modelboot$DZ$stdweights$result)
  # indclassprops_modboot <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # Calculate Entropy
  S_modfit <- entropy(indclassprops_modfit)
  # S_modboot <- entropy(indclassprops_modboot)
  
  
  ### Store results
  results <- c(ax, cx, ex, ay, cy, ey, ra, b12, b21, nMZ+nDZ,
               mX, mY,
               npropMZ[1], npropMZ[2], npropMZ[3], npropMZ[4], npropDZ[1], npropDZ[2], npropDZ[3], npropDZ[4],
               modelfit$MZ$stdweights$result[1], modelfit$MZ$stdweights$result[2], modelfit$MZ$stdweights$result[3], modelfit$MZ$stdweights$result[4],
               modelfit$DZ$stdweights$result[1], modelfit$DZ$stdweights$result[2], modelfit$DZ$stdweights$result[3], modelfit$DZ$stdweights$result[4],
               modelfit$output$estimate["a11"], modelfit$output$estimate["c11"], modelfit$output$estimate["e11"],
               modelfit$output$estimate["a22"], modelfit$output$estimate["c22"], modelfit$output$estimate["e22"],
               modelfit$output$estimate["ra"],
               modelfit$output$estimate["b12"], modelfit$output$estimate["b21"],
               modelfit$output$estimate["mX"], modelfit$output$estimate["mY"],
               modelfit$output$estimate["mX0"], modelfit$output$estimate["mY0"],
               S_modfit)
  dfresults <- rbind(dfresults, results)
  
  colnames(dfresults) <- c("ax", "cx", "ex", "ay", "cy", "ey", "ra", "b12", "b21", "nfam",
                           "mX", "mY",
                           "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                           "estp1", "estp2", "estp3", "estp4",
                           "estp5", "estp6", "estp7", "estp8",
                           "estax", "estca", "estex",
                           "estay", "estcy", "estey",
                           "estra",
                           "estb12", "estb21",
                           "estmXc1", "estmYc1", "estmX0", "estmY0",
                           "entropyfit")
  
  
}  

# rA
for(i in c(0.1, 0.15, 0.20, 0.25, 0.3)){
  nMZ=2500; nDZ=2500
  mX=0.2; mY=0.1
  b12=0.1; b21=0.4
  raFree=TRUE; ra=i
  
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
                       mxPath(from="Axt1", to=c("Axt1", "Ayt1"), free=c(FALSE, FALSE), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
                       mxPath(from="Cxt1", to=c("Cxt1", "Cyt1"), free=c(FALSE, FALSE), value=c(1.0, rc), arrows=2, label=c(NA, "rc")),
                       mxPath(from="Ext1", to=c("Ext1", "Eyt1"), free=c(FALSE, FALSE), value=c(1.0, re), arrows=2, label=c(NA, "re")),
                       mxPath(from="Axt2", to=c("Axt2", "Ayt2"), free=c(FALSE, FALSE), value=c(1.0, ra), arrows=2, label=c(NA, "ra")),
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
                              free=c(TRUE, TRUE, TRUE, TRUE), value=c(0, 0, 0, 0), arrows=1, label=c("mXT1", "mYT1", "mXT2", "mYT2"))
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
  MZ_XYcon <- omxSetParameters(MZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  DZ_XYcon <- omxSetParameters(dzModel,  label=c("b21"), free=FALSE, value=0, newlabels=c("b21zero"), name="DZ_XYcon")
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("b12"), free=TRUE, value=b12, newlabels=c("b12"))
  DZ_XYcon <- omxSetParameters(DZ_XYcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, 0, mY), newlabels=c("mX0", "mY", "mX0", "mY"))
  
  
  # Discordant T1:X->Y T2:X<-Y (b21t1=0  b12t2=0)
  MZ_XYdis <- mzModel;
  MZ_XYdis$A$values[1,2]<-0; MZ_XYdis$A$free[1,2]<-F; MZ_XYdis$A$labels[1,2] <- "b21zero"
  MZ_XYdis$A$values[2,1]<-b12; MZ_XYdis$A$free[2,1]<-T; MZ_XYdis$A$labels[2,1] <- "b12"
  MZ_XYdis$A$values[3,4]<-b21; MZ_XYdis$A$free[3,4]<-T; MZ_XYdis$A$labels[3,4] <- "b21"
  MZ_XYdis$A$values[4,3]<-0; MZ_XYdis$A$free[4,3]<-F; MZ_XYdis$A$labels[4,3] <- "b21zero"
  MZ_XYdis <- mxModel(MZ_XYdis, name="MZ_XYdis")
  MZ_XYdis <- omxSetParameters(MZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  DZ_XYdis <- dzModel;
  DZ_XYdis$A$values[1,2]<-0; DZ_XYdis$A$free[1,2]<-F; DZ_XYdis$A$labels[1,2] <- "b21zero"
  DZ_XYdis$A$values[2,1]<-b12; DZ_XYdis$A$free[2,1]<-T; DZ_XYdis$A$labels[2,1] <- "b12"
  DZ_XYdis$A$values[3,4]<-b21; DZ_XYdis$A$free[3,4]<-T; DZ_XYdis$A$labels[3,4] <- "b21"
  DZ_XYdis$A$values[4,3]<-0; DZ_XYdis$A$free[4,3]<-F; DZ_XYdis$A$labels[4,3] <- "b12zero"
  DZ_XYdis <- mxModel(DZ_XYdis, name="DZ_XYdis")
  DZ_XYdis <- omxSetParameters(DZ_XYdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(0, mY, mX, 0), newlabels=c("mX0", "mY", "mX", "mY0"))
  
  
  # Discordant T1:X<-Y T2:X->Y (b12t1=0  b21t2=0)
  MZ_YXdis <- mzModel; MZ_YXdis$A$values[2,1]<-0; MZ_YXdis$A$free[2,1]<-F; MZ_YXdis$A$labels[2,1] <- "b21zero"
  MZ_YXdis$A$values[1,2]<-b21; MZ_YXdis$A$free[1,2]<-T; MZ_YXdis$A$labels[1,2] <- "b21"
  MZ_YXdis$A$values[4,3]<-b12; MZ_YXdis$A$free[4,3]<-T; MZ_YXdis$A$labels[4,3] <- "b12"
  MZ_YXdis$A$values[3,4]<-0; MZ_YXdis$A$free[3,4]<-F; MZ_YXdis$A$labels[3,4] <- "b21zero"
  MZ_YXdis <- mxModel(MZ_YXdis, name="MZ_YXdis")
  MZ_YXdis <- omxSetParameters(MZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  DZ_YXdis <- dzModel; DZ_YXdis$A$values[2,1]<-0; DZ_YXdis$A$free[2,1]<-F; DZ_YXdis$A$labels[2,1] <- "b21zero"
  DZ_YXdis$A$values[1,2]<-b21; DZ_YXdis$A$free[1,2]<-T; DZ_YXdis$A$labels[1,2] <- "b21"
  DZ_YXdis$A$values[4,3]<-b12; DZ_YXdis$A$free[4,3]<-T; DZ_YXdis$A$labels[4,3] <- "b12"
  DZ_YXdis$A$values[3,4]<-0; DZ_YXdis$A$free[3,4]<-F; DZ_YXdis$A$labels[3,4] <- "b21zero"
  DZ_YXdis <- mxModel(DZ_YXdis, name="DZ_YXdis")
  DZ_YXdis <- omxSetParameters(DZ_YXdis, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, 0, mY), newlabels=c("mX", "mY0",  "mX0", "mY"))
  
  
  # Concordant T1:X<-Y T2:X<-Y (b12t1=0  b12t2=0)
  MZ_YXcon <- omxSetParameters(mzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="MZ_YXcon")
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  MZ_YXcon <- omxSetParameters(MZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0", "mX", "mY0"))
  
  DZ_YXcon <- omxSetParameters(dzModel,  label=c("b12"), free=FALSE, value=0, newlabels=c("b12zero"), name="DZ_YXcon")
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("b21"), free=TRUE, value=b21, newlabels=c("b21"))
  DZ_YXcon <- omxSetParameters(DZ_YXcon, label=c("mXT1", "mYT1", "mXT2", "mYT2"), value=c(mX, 0, mX, 0), newlabels=c("mX", "mY0","mX", "mY0"))
  
  ### Generate Data
  MZ_XYcon_dat <- mxGenerateData(MZ_XYcon, nrows=nMZ*npropMZ[1], empirical=TRUE)
  DZ_XYcon_dat <- mxGenerateData(DZ_XYcon, nrows=nDZ*npropDZ[1], empirical=TRUE)
  
  MZ_XYdis_dat <- mxGenerateData(MZ_XYdis, nrows=nMZ*npropMZ[2], empirical=TRUE)
  DZ_XYdis_dat <- mxGenerateData(DZ_XYdis, nrows=nDZ*npropDZ[2], empirical=TRUE)
  
  MZ_YXdis_dat <- mxGenerateData(MZ_YXdis, nrows=nMZ*npropMZ[3], empirical=TRUE)
  DZ_YXdis_dat <- mxGenerateData(DZ_YXdis, nrows=nDZ*npropDZ[3], empirical=TRUE)
  
  MZ_YXcon_dat <- mxGenerateData(MZ_YXcon, nrows=nMZ*npropMZ[4], empirical=TRUE)
  DZ_YXcon_dat <- mxGenerateData(DZ_YXcon, nrows=nDZ*npropDZ[4], empirical=TRUE)
  
  nSumMZ <- sum(c(dim(MZ_XYcon_dat)[1],dim(MZ_XYdis_dat)[1],dim(MZ_YXdis_dat)[1],dim(MZ_YXcon_dat)[1]))
  nSumDZ <- sum(c(dim(DZ_XYcon_dat)[1],dim(DZ_XYdis_dat)[1],dim(DZ_YXdis_dat)[1],dim(DZ_YXcon_dat)[1]))
  
  ### Merge data
  mzData <- rbind(MZ_XYcon_dat, MZ_XYdis_dat, MZ_YXdis_dat, MZ_YXcon_dat)
  dzData <- rbind(DZ_XYcon_dat, DZ_XYdis_dat, DZ_YXdis_dat, DZ_YXcon_dat)
  
  ### Create mxData Objects
  DataMZ <- mxData(observed=mzData, type="raw")
  DataDZ <- mxData(observed=dzData, type="raw")
  
  
  ### Generate MZ/DZ/Mix Model
  modelMZ <- mxModel(	"MZ", DataMZ, 
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
  
  
  modelDZ <- mxModel("DZ", DataDZ, 
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
  
  
  ### Generate Mixture Model
  mixModel <- mxModel("mixture", modelMZ, modelDZ, mxFitFunctionMultigroup(c("MZ", "DZ")))
  
  # Model Run
  modelfit <- mxTryHard(mixModel, extraTries=3)
  
  # modelboot <- mxBootstrap(modelfit, replications=100)
  
  ### Calculate Entropy
  # Calculate individual class probabilities
  MZ_indclassprops <- indClassProbs(modelfit$MZ)
  DZ_indclassprops <- indClassProbs(modelfit$DZ)
  indclassprops_modfit <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # MZ_indclassprops <- indClassProbs(modelboot$MZ, modelboot$MZ$stdweights$result)
  # DZ_indclassprops <- indClassProbs(modelboot$DZ, modelboot$DZ$stdweights$result)
  # indclassprops_modboot <- rbind(MZ_indclassprops, DZ_indclassprops)
  
  # Calculate Entropy
  S_modfit <- entropy(indclassprops_modfit)
  # S_modboot <- entropy(indclassprops_modboot)
  
  
  ### Store results
  results <- c(ax, cx, ex, ay, cy, ey, ra, b12, b21, nMZ+nDZ,
               mX, mY,
               npropMZ[1], npropMZ[2], npropMZ[3], npropMZ[4], npropDZ[1], npropDZ[2], npropDZ[3], npropDZ[4],
               modelfit$MZ$stdweights$result[1], modelfit$MZ$stdweights$result[2], modelfit$MZ$stdweights$result[3], modelfit$MZ$stdweights$result[4],
               modelfit$DZ$stdweights$result[1], modelfit$DZ$stdweights$result[2], modelfit$DZ$stdweights$result[3], modelfit$DZ$stdweights$result[4],
               modelfit$output$estimate["a11"], modelfit$output$estimate["c11"], modelfit$output$estimate["e11"],
               modelfit$output$estimate["a22"], modelfit$output$estimate["c22"], modelfit$output$estimate["e22"],
               modelfit$output$estimate["ra"],
               modelfit$output$estimate["b12"], modelfit$output$estimate["b21"],
               modelfit$output$estimate["mX"], modelfit$output$estimate["mY"],
               modelfit$output$estimate["mX0"], modelfit$output$estimate["mY0"],
               S_modfit)
  dfresults <- rbind(dfresults, results)
  
  colnames(dfresults) <- c("ax", "cx", "ex", "ay", "cy", "ey", "ra", "b12", "b21", "nfam",
                           "mX", "mY",
                           "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                           "estp1", "estp2", "estp3", "estp4",
                           "estp5", "estp6", "estp7", "estp8",
                           "estax", "estca", "estex",
                           "estay", "estcy", "estey",
                           "estra",
                           "estb12", "estb21",
                           "estmXc1", "estmYc1", "estmX0", "estmY0",
                           "entropyfit")
  
  
}  

