#JCR Replication Files - Networks of Cooperation:
#Rebel Alliances in Fragmented Civil Wars
# Emily Kalah Gade - ekgade@uw.edu
# 3 Jan 2019

rm(list = ls())
set.seed(32314)
options(digits = 3)

library(statnet)
library(network)
library(tidyverse)
library(amen)

# Main document replication (all supplementary analysis below)


################################
# AME Regression Analysis - Table 5 #
################################

# Loading DV matrices for AME regression analysis for 31 groups 
DV_Colab_31<-read.csv("jo_All_Nov2017", stringsAsFactors = F, header = T)

# Making a var with the names of groups in order, to use to match the order of all other variables 
# (This is left as a column because R adds "X" before numeric indicators when it imports, and we have some groups that are abbreviated by number)
colnames<-DV_Colab_31$X # Taken from row names in this first dataset
rownames(DV_Colab_31)<-colnames # Make the list rownames for this data frame
DV_Colab_31<-DV_Colab_31[,-1] # Subtract duplicate row
colnames(DV_Colab_31)<-colnames # Also rename colnames (as some will be listed as "X1st" rather than "1st")

colnames<-colnames(DV_Colab_31)

## Setting up dyadic level and node level controls 

ideol<-read.csv("IdeologyVars_JPR.csv", header = T, stringsAsFactors = F)

# Select covariate data for only 30 infighting groups we are currently analyzing  
ideol2<-filter(ideol, ideol$GroupCode %in% colnames) 

ideolVars<-ideol2[,-1] # Get rid of long form names 

# Cut down to only vars of intrest for this analysis 
ideolVars<-ideolVars[,-10] # Drop power broker score - alternative measure of power
ideolVars<-ideolVars[,-6:-7] # Drop location and state sponsor, will add in later in matrix form
ideolVars$GroupSize <- as.numeric(gsub(",","",ideolVars$GroupSize)) # Remove commas
ideolVars<-ideolVars[,-7] # Remove group size "full"
idv<-ideolVars[match(colnames, ideolVars$GroupCode),]
idv<-na.omit(idv) # AME won't take NAs (that's why we removed state sponsor etc for now)

## now state sponsor node var (just does a group have a sponsor, Y/N)

Nodespons<-as.data.frame(cbind(ideol2$GroupCode, ideol2$State.Sponsors), stringsAsFactors = F) # Create two column dataframe of state sponsors
Nodespons$V2[Nodespons$V2==""]<-NA # Make sure all blanks are read as NA
Nodespons$V2<- ifelse(is.na(Nodespons$V2), 0, 1) # Replace groups that are NA with Zero (no sponsor) or 1 (sponsor)
names(Nodespons)<-c("GroupCode", "spons")
NodeS<-Nodespons[match(colnames, Nodespons$GroupCode),]
a<-array(NodeS[,2]) # Convert to array 

## Add a node level var for ISIL/ASIM (in case these are confounding)

NodeS$V3[NodeS$GroupCode=="ISIL"]<-"1" # Give ISIL a 1
NodeS$V3<-as.numeric(ifelse(is.na(NodeS$V3), 0, NodeS$V3)) # Give all non-ISIL groups  zero
ISIL<-array(NodeS[,3]) # Convert to array

NodeS$V3[NodeS$GroupCode=="ASIM"]<-"1" # Give ASIM a 1
NodeS$V3<-as.numeric(ifelse(is.na(NodeS$V3), 0, NodeS$V3)) # Give all non-ASIM groups  zero
ASIM<-array(NodeS[,3]) # Convert to array


## Creating Xn and Xd for AME analysis 
b<-array(idv[,5]) # Average ideology, could use a compoent of ideology instead
c<-array(idv[,6]) # Size
a_a<-c/1000 # Putting them on similar scale for rough interp. 

Xn<-as.array(cbind(b, a_a, a, ISIL)) # AME takes arrays
dimnames(Xn)[[1]]<-colnames # Labling rows (groups), already ordered/matched above
dimnames(Xn)[[2]]<-c("averageId", "size", "spons", "ISIL") # Lable columns

Xn_c<-as.array(cbind(b, a_a, a, ASIM)) # AME takes arrays
dimnames(Xn_c)[[1]]<-colnames # Labling rows (groups), already ordered/matched above
dimnames(Xn_c)[[2]]<-c("averageId", "size", "spons", "ASIM") # Lable columns



## Creating Xn and Xd for AME analysis  WITHOUT ISIS

Xn_b<-as.array(cbind(b, a_a, a)) # AME takes arrays
dimnames(Xn_b)[[1]]<-colnames # Labling rows (groups), already ordered/matched above
dimnames(Xn_b)[[2]]<-c("averageId", "size", "spons") # Lable columns


## Constructing dyadic vars
ad<-as.data.frame(Xn) # Generate a dataframe from your array
ideol_diff<- outer(ad$averageId, ad$averageId, "-") # Make a matrix where each square represents the difference in ideology scores of a given dyad
ideol_diff<-abs(ideol_diff) # Has to be postive - so take the absolute value 

power_diff<- outer(ad$size, ad$size, "-") # Replicates the above for power
power_diff<-abs(power_diff)

# Load State Sponsor and Location Matrcies 
loc<-read.csv("LocationDiff_CollabPaper_Oct2018_31Groups.csv", stringsAsFactors = F, header = T) # Load shared location matrix 

x<-loc$X # Deal with row and column names as above
loc<-loc[,-1]
colnames(loc)<-x
rownames(loc)<-x
dim(loc) # Check dimentions 
loc<-as.matrix(loc) # Convert to matrix format
dimnames(loc)[[1]]<-x
dimnames(loc)[[2]]<-x

locFinal<-loc[colnames,] # Make order of "loc" match order of ideology and everything else
locFinal<-locFinal[,colnames] # And match on the other axis... 

## state sponsorship

spons<-read.csv("stateSponsonershipCOLLAB_Oct2018-31Groups.csv", stringsAsFactors = F, header = T)
x<-spons$X 
spons<-spons[,-1]
rownames(spons)<-x
colnames(spons)<-x
spons<-as.matrix(spons)
dimnames(spons)[[1]]<-x
dimnames(spons)[[2]]<-x
spons<-as.matrix(spons)
sponsFinal<-spons[colnames,]
sponsFinal<-sponsFinal[,colnames]


# Create an empty array to hold your matrcies (dyad vars) for AME
Xd<- array(dim=c(nrow(DV_Colab_31), nrow(DV_Colab_31), 4)) # Must be same dimentions as number of groups (rows and columns) and have enough "slices" to fit all dyad vars  
Xd[,,1]<- ideol_diff # Average ideology dyad
Xd[,,2]<-power_diff  # Size dyad
Xd[,,3]<-as.matrix(locFinal) # location shared dyad
Xd[,,4]<-as.matrix(spons) #sponsFinal # sponsorship shared dyad 

dimnames(Xd)[[1]]<-colnames
dimnames(Xd)[[2]]<-colnames
dimnames(Xd)[[3]]<-c("ideol_diff", "powerdiff", "loc", "spons")

# Modeling 
diag(DV_Colab_31)<-0 #AMEN won't take NAs
Yreal<-as.matrix(DV_Colab_31) # For raw count 
Yrealsqr<-sqrt(Yreal) # For square root - best practice as AME doesn't take a count model, 
#this should better approximate a normal distirubiton

### for sqrt transformed 
fit_collab_Bivariate_power_sqrt<-ame(Yrealsqr, power_diff, R=1, model="nrm",
                                     symmetric=TRUE,burn=10000,nscan=10000,odens=10)
fit_collab_Bivariate_ideol_sqrt<-ame(Yrealsqr, ideol_diff, R=1, model="nrm",
                                     symmetric=TRUE,burn=10000,nscan=10000,odens=10)
fit_collab_Bivariate_spons_sqrt<-ame(Yrealsqr, sponsFinal,  R=1, model="nrm",
                                     symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# sqrt - no ISIS
fit_collab_nodesdyads_noIsis_sqrt<-ame(Yrealsqr, Xd, Xn_b, R=1, model="nrm",
                                       symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# sqrt - WITH ISIS
fit_collab_nodesdyads_Isis_sqrt<-ame(Yrealsqr, Xd, Xn, R=1, model="nrm",
                                     symmetric=TRUE,burn=10000,nscan=10000,odens=10)


# sqrt - WITH ASIM
fit_collab_nodesdyads_ASIM_sqrt<-ame(Yrealsqr, Xd, Xn_c, R=1, model="nrm",
                                     symmetric=TRUE,burn=10000,nscan=10000,odens=10)


### doing some quick predicted prob plus marginal effects stuff
#YPM<-fit_collab_nodesdyads_ASIM_sqrt$YPM

# Create an empty array to hold your matrcies (dyad vars) for AME
Xd_noidea<- array(dim=c(nrow(DV_Colab_31), nrow(DV_Colab_31), 3)) # Must be same dimentions as number of groups (rows and columns) and have enough "slices" to fit all dyad vars  
#Xd[,,1]<- ideol_diff # Average ideology dyad
Xd_noidea[,,1]<-power_diff  # Size dyad
Xd_noidea[,,2]<-as.matrix(locFinal) # location shared dyad
Xd_noidea[,,3]<-as.matrix(sponsFinal) # sponsorship shared dyad 

dimnames(Xd_noidea)[[1]]<-colnames
dimnames(Xd_noidea)[[2]]<-colnames
dimnames(Xd_noidea)[[3]]<-c("powerdiff", "loc", "spons")

# sqrt - WITH ASIM
fit_collab_nodesdyads_ASIM_sqrt_noidea<-ame(Yrealsqr, Xd_noidea, Xn_c, R=1, model="nrm",
                                            symmetric=TRUE,burn=10000,nscan=10000,odens=10)


ybin <-Yreal
ybin[ybin>0] <-1 # binary DV


### for count transformed 
fit_collab_Bivariate_power<-ame(Yreal, power_diff, R=1, model="nrm",
                                symmetric=TRUE,burn=10000,nscan=10000,odens=10)
fit_collab_Bivariate_ideol<-ame(Yreal, ideol_diff, R=1, model="nrm",
                                symmetric=TRUE,burn=10000,nscan=10000,odens=10)
fit_collab_Bivariate_spons<-ame(Yreal, sponsFinal,  R=1, model="nrm",
                                symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# count - no ISIS
fit_collab_nodesdyads_noIsis<-ame(Yreal, Xd, Xn_b, R=1, model="nrm",
                                  symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# count - WITH ISIS
fit_collab_nodesdyads_Isis<-ame(Yreal, Xd, Xn, R=1, model="nrm",
                                symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# count - WITH ASIM
fit_collab_nodesdyads_ASIM<-ame(Yreal, Xd, Xn_c, R=1, model="nrm",
                                symmetric=TRUE,burn=10000,nscan=10000,odens=10)

YPM_nonSqrt<-fit_collab_nodesdyads_ASIM$YPM


# ordinal
fit_nodesdyads_ord <-ame(Yreal, Xd, Xn_b, R=1, model="ord", symmetric=TRUE,burn=10000,nscan=10000,odens=10)



summary(fit_collab_Bivariate_power_sqrt)
summary(fit_collab_Bivariate_ideol_sqrt)
summary(fit_collab_Bivariate_spons_sqrt)
summary(fit_collab_nodesdyads_noIsis_sqrt)
summary(fit_collab_nodesdyads_Isis_sqrt)

# for non-square root transformed, re-run models with "Yreal"


summary(fit_collab_Bivariate_power)
summary(fit_collab_Bivariate_ideol)
summary(fit_collab_Bivariate_spons)
summary(fit_collab_nodesdyads_noIsis)
summary(fit_collab_nodesdyads_Isis)
summary(fit_collab_nodesdyads_ASIM)
summary(fit_nodesdyads_ord)


# use summary(fit_xxx) to see a printout of findings

# To see AME diagnostic plots: plot(name of the fit object) - these are in Supplementary Materials


##### plotting the latent space


fit_yreal<-ame(Yreal, R=1, model="nrm",
               symmetric=TRUE,burn=10000,nscan=10000,odens=10)

fit_yrealsqrt<-ame(Yrealsqr, R=1, model="nrm",
                   symmetric=TRUE,burn=10000,nscan=10000,odens=10)

fit_collab_Bivariate_ideol<-ame(Yreal, ideol_diff2, R=1, model="nrm",
                                symmetric=TRUE,burn=10000,nscan=10000,odens=10)
plot(fit_yreal$U, b)
plot(fit_yrealsqrt$U, b)
plot(fit_yreal$U, fit_collab_Bivariate_ideol$U)
plot(fit_yrealsqrt$U, fit_collab_Bivariate_ideol_sqrt$U)

#pdf(file="LatentSpace_compare_3", paper="letter",width = 7,height = 5)
par(mfrow=c(2,1))

#circplot(Yreal)
u<-fit_nodesdyads_noIsis_Yrealsqr2$U
cnames<-rownames(u)
urank<-rank(u) 

plot(urank,u,type="n",xlab="rank order of u")
abline(h=0,col="gray") 
addlines(Yrealsqr>0 , cbind(urank,u),col="green")
addlines(Yrealsqr<0 , cbind(urank,u),col="pink")
text(urank,u,cnames,srt=-45,cex=1.0) 


u2<-rawYreal$U
cnames2<-rownames(u2)
urank2<-rank(u2) 
plot(urank2,u2,type="n",xlab="rank order of u2")
abline(h=0,col="gray") 
addlines(Yrealsqr>0 , cbind(urank2,u2),col="green")
addlines(Yrealsqr<0 , cbind(urank2,u2),col="pink")
text(urank,u,cnames,srt=-45,cex=1.0) 

circplot(Yrealsqr)

u3<-fit_Bivariate_ideol$U
cnames2<-rownames(u3)
urank3<-rank(u3) 
plot(urank3,u3,type="n",xlab="rank order of u3")
abline(h=0,col="gray") 
addlines(Yrealsqr>0 , cbind(urank3,u3),col="green")
addlines(Yrealsqr<0 , cbind(urank3,u3),col="pink")
text(urank,u,cnames,srt=-45,cex=1.0) 

dev.off()






