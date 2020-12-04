library(devtools)
library(dplyr)
library(tibble)
library(mice)
library(data.table)
library(bc3)
#library(bc3, lib.loc ="/home/sandovall2/R/R-3.6.1/")

############ READ PREPROCESSED DATA----------------------------------------------------------------
getwd()
setwd("C:/Users/sandovall2/Box/BCAST Risk Factor/Datasets/Input Data")
origin_data<- data.frame(read.csv("analyse_cleaned_reduced_20200403.csv"))

#testing
# ############ Add tslb column where time since last birth = refage-lastchildage
# origin_data$tslb<- NA
# idx<-which(origin_data$parity!=0  & !is.na(origin_data$lastChildAge))
# origin_data$tslb[idx]<- origin_data$refage[idx]-origin_data$lastChildAge[idx]
# summary(origin_data$tslb)

############ check missings in lastchildage for parous women
print("------------ Pre-imputation: parous women------------")
idx<- which(origin_data$parous==1)
nrow(origin_data[idx,])
print("------------ Pre-imputation: nulliparous women------------")
idx<- which(origin_data$parous==0)
nrow(origin_data[idx,])
print("------------ Pre-imputation: Missing lastchildage among parous women------------")
idx<- which(origin_data$parous==1 & is.na(origin_data$lastChildAge))
nrow(origin_data[idx,])

############ make list of variables to be used and excluded in imputations
non_imp_cols<-colnames(select(origin_data,-c("ageMenarche","parous","parity","mensAgeLast","postmeno","ageFFTP","breastMos","lastChildAge","study","refage","StudyCountry")))
imp_cols<- colnames(select(origin_data,c("ageMenarche","parous","parity","mensAgeLast","postmeno","ageFFTP","breastMos","lastChildAge","study","refage","StudyCountry")))

print("----------------- Summary of data pre-imputation ---------------")
summary(select(origin_data,imp_cols))

############################ IMPUTATION 1: PAROUS AND AGEMENARCHE  ################################################
print("BEGIN 1ST IMPUTATION:  PAROUS AND AGEMENARCHE")
# initialize imputation parameters
init <- mice(origin_data, maxit = 0)
# remove variables as a predictors 
predM = init$predictorMatrix
predM[, non_imp_cols]=0
# skip variables from imputation
meth = init$method
meth[non_imp_cols]=""


first_imp <- mice(origin_data,m=1,seed=1,method=meth,predictorMatrix = predM, visitSequence = "monotone", maxit=3)  #imputation 1, removed seed
first_imp$method
imp1=complete(first_imp,1)
print("summary: parity, parous (post-imp1), and parous (pre-imp1)")
summary(imp1$parity)
summary(imp1$parous)
summary(origin_data$parous)

print("----reset values of parity,tslb,lastChildAge,ageFFTP,breastMos,mensAgeLast,and ageMenarche-----")
imp1$parity = origin_data$parity
imp1$lastChildAge = origin_data$lastChildAge
imp1$ageFFTP = origin_data$ageFFTP
imp1$breastMos = origin_data$breastMos
imp1$mensAgeLast = origin_data$mensAgeLast

print("-----------Samples where parous=0 and parity is missing-----------")
idx<- which(imp1$parous==0 & is.na(imp1$parity)) 
nrow(imp1[idx,])
print("----------- Fixed: if parous=0, then parity should be zero---------")
imp1$parity[idx]<-0
print("----------- Fix: if parous=0, then breastMos should be zero---------")
idx<- which(imp1$parous==0 & is.na(imp1$breastMos))
imp1$breastMos[idx]<-0

print("-----------Samples where parous=1 and parity is missing-----------")
idx<- which(imp1$parous==1 & is.na(imp1$parity))
nrow(select(imp1,imp_cols)[idx,])

# checks for parous and parity
print("----Number of nulliparous women")
idx= which(imp1$parous==0)
nrow(imp1[idx,])
print("----Number of women with parity=0")
idx1= which(imp1$parity==0)
nrow(imp1[idx1,])
print("----Number of parous women")
idx= which(imp1$parous==1)
nrow(imp1[idx,])
print("----Number of women with parity>0 or parity=NA")
idx1= which(is.na(imp1$parity) | imp1$parity>0)
nrow(imp1[idx1,])

######################## IMPUTATION 2: PARITY (PAROUS SAMPLES) ##########################################
# Subset Parous and nulliparous
Parous= imp1[imp1$parous==1,] #parous
Nulliparous= imp1[imp1$parous==0,] #nulliparous
summary(select(Parous,imp_cols))


# initialize imputation parameters
init <- mice(Parous, maxit = 0)
# remove variables as a predictors 
predM = init$predictorMatrix
predM[, non_imp_cols]=0
# skip variables from imputation
meth = init$method
meth[non_imp_cols]=""
# set post processing paramters
post <- init$post
# ensure imputed parity is between 1 and 16
post["parity"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(1,16))"

second_imp <- mice(Parous,m=1,method=meth, predictorMatrix=predM, post= post,seed=2, visitSequence = "monotone", maxit=5,print=FALSE)  #imputation 1, removed seed
second_imp$method
imp2=complete(second_imp,1)

############## Return some imputed columns to original with missings to feed into third imputation
print("----reset values of mensAgeLast,lastChildAge,ageFFTP,breastMos-----")
imp2$mensAgeLast <- Parous$mensAgeLast
imp2$lastChildAge <- Parous$lastChildAge
imp2$ageFFTP <- Parous$ageFFTP
imp2$breastMos <- Parous$breastMos

print("------------------summary of imp2 ----------------")
summary(select(imp2,imp_cols))



######################## IMPUTATION 3: AGEFFTP (PAROUS SAMPLES) ##########################################

# initialize imputation parameters
init <- mice(imp2, maxit = 0)
# remove variables as a predictors 
predM = init$predictorMatrix
predM[, non_imp_cols]=0
# skip variables from imputation
meth = init$method
meth[non_imp_cols]=""
# set post processing paramters
post <- init$post
# ensure imputed ageFFTP is equal to or less than refage
post["ageFFTP"] <- "imp[[j]][, i] <- pmin(imp[[j]][, i], init$data$refage[!r[,j]])"


third_imp <- mice(imp2,m=1,method=meth, predictorMatrix=predM, post= post,seed=3, visitSequence = "monotone", maxit=2,print=FALSE)  #imputation 1, removed seed
third_imp$method
imp3=complete(third_imp,1)

############## Return some imputed columns to original with missings to feed into fourth imputation
print("----reset values of mensAgeLast,lastChildAge,breastMos-----")
imp3$mensAgeLast <- imp2$mensAgeLast
imp3$lastChildAge <- imp2$lastChildAge
imp3$breastMos <- imp2$breastMos

print("---------- summary of imp3-------------")
summary(select(imp3,imp_cols))

############## Check results
# ageFFTP<= refage
idx = which(imp3$ageFFTP>imp3$refage)
ageFFTP_errors_imp3 = imp3[idx,imp_cols]
print("-------Check results of imp3, # samples where agefftp>refage:-----------")
nrow(ageFFTP_errors_imp3[idx,])

############## manual post-processing: set lastchildage to ageFFTP when parity=1 
idx = which(imp3$parity==1 & is.na(imp3$lastChildAge))
imp3$lastChildAge[idx] = imp3$ageFFTP[idx]

########################## IMPUTATION 4: lastChildAge,breastMos  #########################################

print("BEGIN 4th IMPUTATION: lastChildAge,breastMos")
# initialize imputation parameters
init <- mice(imp3, maxit = 0)

#pred<-quickpred(imp3, minpuc=0.15, mincor = 0.15)

# remove variables as a predictors 
predM = init$predictorMatrix
predM[, non_imp_cols]=0

# skip variables from imputation
meth = init$method
meth[non_imp_cols]=""

# set post processing paramters
post <- init$post
#post["lastChildAge"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i],c(init$data$ageFFTP[!r[,j]], init$data$refage[!r[,j]]))"
post["lastChildAge"] <- "imp[[j]][, i] <-pmin(pmax(imp[[j]][, i], init$data$ageFFTP[!r[,j]]), init$data$refage[!r[,j]])"

fourth_imp <- mice(imp3,maxit=5,m=1,seed=4,post=post, method=meth,predictorMatrix=predM)
print("Methods in fourth imputation")
fourth_imp$method
#fourth_imp$loggedEvents
imp4= complete(fourth_imp,1)
hist(imp4$lastChildAge)
print("-------summary imp4$lastchildage-------------")
summary(imp4$lastChildAge)

print("-----------rows where lastchildage is erroneoulsy greater than refage or less than agefftp-----------")
idx<- which(imp4$lastChildAge<imp2$ageFFTP | imp4$lastChildAge>imp2$refage)
lastch <- select(imp4[idx,],imp_cols)
nrow(lastch) 

# check distribution of lastchildage pre and post imputation
hist1<- hist(Parous$lastChildAge)
hist2<- hist(imp4$lastChildAge)

############## Reset mensagelast and tslb values back
imp4$mensAgeLast <- imp3$mensAgeLast

print("---------- summary of imp4-------------")
summary(select(imp4,imp_cols))

############# Combine rows of parous and nulliparous samples 
imp4=rbind(imp4,Nulliparous)

############# Add tslb data to newly imputed parous samples
imp4$tslb = NA
idx = which(imp4$parous==1)
imp4$tslb[idx]<- imp4$refage[idx]-imp4$lastChildAge[idx]


print("---------Negative tslb values------")
idx<- which(imp4$tslb<0)
neg_tslb <- imp4[idx,]
nrow(neg_tslb )

###################################### FIFTH IMPUTATION OF MENSAGELAST ###################################### 

############# Subset data based on pre-menopausal and post-menopausal
print("--------------BEGIN FIFTH IMPUTATION-------------")
PostMeno= imp4[imp4$postmeno==1,] #postmenopausal should not have missing mensagelast
PreMeno= imp4[imp4$postmeno==0,] #premeno should have missing mensagelast

# initialize imputation parameters
init <- mice(PostMeno, maxit = 0)
# remove variables as a predictors 
predM = init$predictorMatrix
predM[, non_imp_cols]=0
# skip variables from imputation
meth = init$method
meth[non_imp_cols]=""

fifth_imp <- mice(PostMeno,predictorMatrix=predM, method=meth,visitSequence="monotone", maxit=2,m=1,print=FALSE) # 3rd imp of post menopausal women
print("5th imputation complete: mensagelast")
fifth_imp$method
imp5= complete(fifth_imp,1)

############# Reset lastChildAge, ageFFTP, breastMos and tslb except mensagelast
imp5$lastChildAge<-PostMeno$lastChildAge
imp5$ageFFTP<-PostMeno$ageFFTP
imp5$breastMos<-PostMeno$breastMos
imp5$tslb<-PostMeno$tslb
summary(imp4$breastMos)
print("----------Summary of imp5-------------")
summary(select(imp5,imp_cols))

############# Combine imputed and non-imputed rows
imp5<-rbind(imp5,PreMeno)

print("----- Return imputed some imputed columns to original with missings----")
summary(imp5)
table(imp5$parous,imp5$parity)

print("---- Summary of lastChildAge in imp5-------")
summary(imp5$lastChildAge)
print("--------Summary of tslb in imp5----------")
summary(imp5$tslb)

summary(select(imp5,imp_cols))
print("--------Table of parous in imp5----------")
table(imp5$parous)

#save subset and preprocessed file
#save.image(file="three_imputations_workspace_incl_parous.Rdata")
write.csv(imp3, paste0("5imp_local_120320.csv"), row.names = FALSE)

#################################### CATEGORIZATION OF IMPUTED DATA  ###################################### 

data<-imp5
print("---------------------------Total samples ---------------------------")
print(nrow(data))

############### FOCUS ON CONTROLS AND INVASIVE BREAST CANCER CASES

print("--------- Only keep samples of controls and invasive breast cancer cases (o,1) ---------")
idx<- which(data$status==0|data$status==1)
data <- data[idx,]
print(nrow(data))

print("--------------------------- Table of status ---------------------------")
table(data$status)
print("--------------------------- Table of study ---------------------------")
table(data$study)
summary(data$study)

############## TUMOR MARKERS

print("--------------------------- Table of study by ER_status1 ---------------------------")
table(data$study,data$ER_status1)
summary(data$ER_status1)
print("--------------------------- Table of study by PR_status1 ---------------------------")
table(data$study,data$PR_status1)
print("--------------------------- Table of study by HER2_status1 ---------------------------")
table(data$study,data$HER2_status1)
print("--------------------------- Table of study by Grade1 ---------------------------")
table(data$study,data$Grade1)


############## CATEGORIZE  IMPUTED DATA
#---------------------------------------------------------------------------------------
# First round decimals down breastMos, ageMenarche, ageFFTP and mensAgeLast variables
data$tslb<- ceiling(data$tslb)
data$lastChildAge<- ceiling(data$lastChildAge)
data$breastMos<- ceiling(data$breastMos)
data$ageMenarche<- ceiling(data$ageMenarche)
data$ageFFTP<- ceiling(data$ageFFTP)
data$mensAgeLast<- ceiling(data$mensAgeLast)

# * subgroups-----------------------------------------------------------------------------
data$subgroup<-NA
data$subgroup<- 9
data$subgroup[data$status==3]<- 9
data$subgroup[data$status==2]<- 3
data$subgroup[data$ER_status1==0]<- 2
data$subgroup[data$ER_status1==1]<- 1
data$subgroup[data$status==0]<- 0

# * 5 molecular subtypes;(Audrey)------------------------------------------------------------------

data$molgroup=NA
data$molgroup[(data$ER_status1==1 | data$PR_status1==1) & (data$HER2_status1==0) & (data$Grade1==1 | data$Grade1==2) ]<- 1
data$molgroup[(data$ER_status1==1 | data$PR_status1==1) & (data$HER2_status1==0) & (data$Grade1==3) ]<- 2
data$molgroup[(data$ER_status1==1 | data$PR_status1==1) & (data$HER2_status1==1)]<- 3
data$molgroup[(data$ER_status1==0 & data$PR_status1==0) & (data$HER2_status1==1)]<- 4
data$molgroup[(data$ER_status1==0 & data$PR_status1==0) & (data$HER2_status1==0)]<- 5
# from above: if status=0 then do; subgroup=0; molgroup=0; /* er_status1=0; pr_status1=0; her2_status1=0; grade1=0; */ end;
data$molgroup[(data$status==0)]<- 0

################### * AGE AT MENARCHE;--------------------------------------------------

data$agemenarche_cat<-NA
data$agemenarche_cat[data$ageMenarche<8 | data$ageMenarche>26 ]<- 9 # 0:7 & >26
data$agemenarche_cat[data$ageMenarche >= 15 & data$ageMenarche <= 25]<- 0 #15:25
data$agemenarche_cat[data$ageMenarche==14]<- 1            #14
data$agemenarche_cat[data$ageMenarche==13]<- 2            #13
data$agemenarche_cat[data$ageMenarche >=8 & data$ageMenarche <=12 ]<- 3   #8:12


print("----------------- Table and summary of data$agemenarche_cat -----------------")
table(data$agemenarche_cat)
summary(data$agemenarche_cat)
print("------- Table and summary of data$agemenarche_cat with collapsed cat9(8> ageMenarche >12) into cat3 (8<= ageMenarche <=12 ) -------")
data$agemenarche_cat[data$agemenarche_cat==9] <- 3
table(data$agemenarche_cat)
summary(data$agemenarche_cat)

############## * PARITY CATEGORIES

data$parity_cat=NA
# if parity<=2 then parity_cat=parity;
data$parity_cat[data$parity<=2 & !is.na(data$parity)]<-data$parity[data$parity <= 2 & !is.na(data$parity)]
data$parity_cat[data$parity >= 3 ]<- 3  # >=3
data$parity_cat[is.na(data$parity) | data$parity >25]<- 9 # there shouldn't be parity missings or greater than 16
data$parity_cat[data$parity ==0]<- 0

print("----------------- Table of data$parity_cat -----------------")
table(data$parity_cat)

print("----------------- table of data$parity -----------------")
table(data$parity)
print("----------------- Number of samples of parous women -----------------")
idx<-which(data$parity>=1)
length(idx)
print("----------------- Number of samples of nulliparous women -----------------")
idx<-which(data$parity==0)
length(idx)

############## * AGE AT FIRST FULL-TERM PREGNANCY;----------------------------------------------------

data$agefftp_cat<- NA
data$agefftp_cat[data$ageFFTP >= 30 & data$ageFFTP <= 65]<- 4  # RANGE: 30:65
data$agefftp_cat[data$ageFFTP >= 25 & data$ageFFTP <= 29]<- 3  # RANGE: 25:29
data$agefftp_cat[data$ageFFTP >= 20 & data$ageFFTP <= 24]<- 2 # RANGE: 20-24
data$agefftp_cat[data$ageFFTP >= 10 & data$ageFFTP <= 19]<- 1 # RANGE: 10:19
data$agefftp_cat[(data$ageFFTP >= 1 & data$ageFFTP <= 9) | data$ageFFTP > 65]<- 9 # RANGE: 1-9 & >65
data$agefftp_cat[data$parity==0 ]<- 0 # RANGE: 0

print("----------------- Table and summary of data$agefftp_cat -----------------")
table(data$agefftp_cat)
summary(data$agefftp_cat)


############## * TSLB = NA PARITY=0, CONVERT TO 0s

print("----------------- TSLB=NA PARITY=0 CONVERSION TO TSLB=0 -------------")
idx<- which(is.na(data$tslb) & data$parity==0)
print("----------------- Number of missings in tslb where parity=0 -------------")
length(idx)
print("----------------- Summary tslb continuous -----------------")
table(data$tslb)
summary(data$tslb)

data$tslb[idx]<-0 
print("------------- Summary tslb after missings replaced with 0 ------------")
table(data$tslb)
summary(data$tslb)

############## * TSLB = 0 PARITY!= 0, CONVERT TO 0.5s

print("----------------- TSLB=NA PARITY!=0 CONVERSION TO TSLB=0.5 -------------")
idx<- which(data$tslb ==0& data$parity!=0)
print("----------------- Number of tslb=0 where parity!=0 -----------------")
length(idx) 
print(" ----------------- Table tslb by parity -----------------")
table(data$tslb, data$parity)
print(" ----------------- New table tslb by parity with 0.5 category-----------------")
data$tslb[idx]<-0.5 
table(data$tslb, data$parity)

############## * LAST CHILD AGE -------------------------------------------------------------------
print(" ----------------- lastChildAge cat -----------------")
data$lastchildage_cat=NA
data$lastchildage_cat<- 4                               # RANGE: 35-65
data$lastchildage_cat[floor(data$lastChildAge) <35 ]<- 3 # RANGE: 30-34
data$lastchildage_cat[floor(data$lastChildAge) <30 ]<- 2 # RANGE: 25-29
data$lastchildage_cat[floor(data$lastChildAge) <25 ]<- 1  # RANGE: 10-24
data$lastchildage_cat[floor(data$lastChildAge) <10 | floor(data$lastChildAge) >65]<- 9 # RANGE: 0-9 & >65

############## * time since last birth tslb; (Audrey) ------------------------------------------------
print(" ----------------- tslb cat -----------------")
data$tslb_cat=NA # Add tslb_cat column
data$tslb_cat[data$tslb >= 55 & data$parity!=0]<- 12
data$tslb_cat[data$tslb >= 50 & data$tslb < 55 & data$parity!=0]<- 11
data$tslb_cat[data$tslb >= 45 & data$tslb < 50 & data$parity!=0]<- 10
data$tslb_cat[data$tslb >= 40 & data$tslb < 45 & data$parity!=0]<- 9
data$tslb_cat[data$tslb >= 35 & data$tslb < 40 & data$parity!=0]<- 8
data$tslb_cat[data$tslb >= 30 & data$tslb < 35 & data$parity!=0]<- 7
data$tslb_cat[data$tslb >= 25 & data$tslb < 30 & data$parity!=0]<- 6
data$tslb_cat[data$tslb >= 20 & data$tslb < 25 & data$parity!=0]<- 5
data$tslb_cat[data$tslb >= 15 & data$tslb < 20 & data$parity!=0]<- 4
data$tslb_cat[data$tslb >= 10 & data$tslb < 15 & data$parity!=0]<- 3
data$tslb_cat[data$tslb >= 5 & data$tslb < 10 & data$parity!=0]<- 2 #5:9
data$tslb_cat[data$tslb >= 0 & data$tslb < 5 & data$parity!=0]<- 1 #0:4
data$tslb_cat[data$parity==0]<- 0 #nulliparous
print("---------------Table and summary of data$tslb_cat---------------")
table(data$tslb_cat)
summary(data$tslb_cat)

############## Check for invalid negative values in tslb
# idx<-which(data$tslb_cat==12)
# print("Total tslb_cat==12 (greater or equal to 55)")
# length(idx)
# 
# print("Negative values in tslb:")
# idx<-which(data$tslb < 0)
# check_tslb<- data$tslb[idx]
# length(check_tslb)


#---------------------------------------------------------------------------------------
############ * AGE AT LAST MENSTRUATION; 
print(" ----------------- mensAgeLast cat -----------------")
print("----------- Has NAs for premenopausal women ---------")
idx<- which(data$postmeno==0)
print("----------- Total premenopausal women ---------")
length(idx)

data$mensagelast_cat <-NA
data$mensagelast_cat[data$postmeno==0 ]<- 1 # premeno
data$mensagelast_cat[data$mensAgeLast >=10 & data$mensAgeLast <= 49]<- 0 #10:49
data$mensagelast_cat[data$mensAgeLast >= 50 & data$mensAgeLast <= 53 ]<- 2 #50:53
data$mensagelast_cat[data$mensAgeLast >= 54 & data$mensAgeLast <= 65]<- 3 #55:65
data$mensagelast_cat[(data$mensAgeLast >= 1 & data$mensAgeLast <= 9) | data$mensAgeLast > 65]<- 9 # outliers 1:9, >65


print("--------------- Table and summary of data$mensagelast_cat ---------------")
table(data$mensagelast_cat)
summary(data$mensagelast_cat)

print("------ Table and summary of data$mensagelast_cat after collapsing outlier category (data$mensAgeLast >= 1 & data$mensAgeLast <= 9) into cat 3 (<=65)-------")
data$mensagelast_cat[data$mensagelast_cat==9] <- 3
table(data$mensagelast_cat)
summary(data$mensagelast_cat)

############ * AGE AT MENOPAUSE MEDIAN CATEGORIES;---------------------------------------
# data$mensagelast_cat_med <-NA
# idx<-which(data$postmeno==0) 
# data$mensagelast_cat_med[idx]<-1 # Pre-menopausal women
# idx<-which(data$mensAgeLast >=10 & data$mensAgeLast <= 49) 
# data$mensagelast_cat_med[idx]<-median(data$mensAgeLast[idx]) # AGE AT MENOPAUSE RANGE: 10-49
# idx<-which(data$mensAgeLast >= 50 & data$mensAgeLast <= 53) 
# data$mensagelast_cat_med[idx] <- median(data$mensAgeLast[idx]) # AGE AT MENOPAUSE RANGE: 50-53
# idx<-which(data$mensAgeLast >= 54 & data$mensAgeLast <= 65) 
# data$mensagelast_cat_med[idx] <- median(data$mensAgeLast[idx])  # AGE AT MENOPAUSE RANGE: 54-65
# idx<-which(data$mensAgeLast >= 1 & data$mensAgeLast <= 9 | data$mensAgeLast > 65 ) # RANGE: 1-9 & >65
# data$mensagelast_cat_med[idx] <- 99

#################### * BREASTFEEDING DURATION;------------------------------------------
print(" ----------------- breastMos cat -----------------")
data$breastmos_cat=NA
data$breastmos_cat[data$breastMos>300]<- 9
data$breastmos_cat[data$breastMos>=25 & data$parity<=300 ]<- 5 # 25:299
data$breastmos_cat[data$breastMos>=13 & data$breastMos<=24 ]<- 4
data$breastmos_cat[data$breastMos>=7 & data$breastMos <= 12 ]<- 3
data$breastmos_cat[data$breastMos>=1 & data$breastMos<=6]<- 2
data$breastmos_cat[data$breastMos == 0 & data$parity!=0 ]<- 1
data$breastmos_cat[data$parity ==0]<- 0 #nulliparous

print("--------------- Table and summary of data$breastmos_cat ---------------")
table(data$breastmos_cat)
summary(data$breastmos_cat)

print("---------------Table and summary of data$breastMos ---------------")
table(data$breastMos)
summary(data$breastMos)

#################### * BREASTFEEDING DURATION MEDIAN CATEGORIES; ------------------------------------------
# data$breastmos_cat_med <- NA
# data$breastmos_cat_med[data$parity==0]<-1 #nulliparous coded as one TO DO
# data$breastmos_cat_med[data$breastMos == 0 & data$parity!=0 ]<-0 # MEDIAN FOR BREAST FEEDING DURATION RANGE: 0
# idx<-which(data$breastMos>=1 & data$breastMos<=6) 
# data$breastmos_cat_med[idx]<- median(data$breastMos[idx]) # MEDIAN FOR BREAST FEEDING DURATION RANGE: 1-6
# idx<-which(data$breastMos>=7 & data$breastMos<=12) # MEDIAN FOR BREAST FEEDING DURATION RANGE: 7-12
# data$breastmos_cat_med[idx]<- median(data$breastMos[idx]) 
# idx<-which(data$breastMos>=13 & data$breastMos<=24) # MEDIAN FOR BREAST FEEDING DURATION RANGE: 13-24
# data$breastmos_cat_med[idx] <- median(data$breastMos[idx]) 
# idx<-which(data$breastMos>300)                      # MEDIAN FOR BREAST FEEDING DURATION RANGE: >300
# data$breastmos_cat_med[idx] <- 99

#################### * Cross tabulation parity_tslb------------------------------------------
# 
data$parity_tslb=NA # add parity_tslb column
# if parity_cat = 0 then parity_tslb = 0;
idx <- which(data$parity_cat==0)
length(idx)
data$parity_tslb[idx] <- 0
# else if parity_cat = 1 then do; %create_cross(tslb_cat, 1 2 3 4 5 6 7 8 9 10 11 12, parity_tslb,  1  2  3  4  5  6  7  8  9 10 11 12); end;
idx <- which(data$parity_cat==1)
length(idx)
data$parity_tslb[idx] <- data$tslb_cat[idx]
print("table(data$tslb_cat 1:12, parity_cat=1)")
table(data$tslb_cat[idx])
# else if parity_cat = 2 then do; %create_cross(tslb_cat, 1 2 3 4 5 6 7 8 9 10 11 12, parity_tslb, 13 14 15 16 17 18 19 20 21 22 23 24); end;
idx <- which(data$parity_cat==2)
length(idx)
data$parity_tslb[idx] <- data$tslb_cat[idx]+12
print("table(data$tslb_cat 13:24, parity_cat=2)")
table(data$tslb_cat[idx])
# else if parity_cat = 3 then do; %create_cross(tslb_cat, 1 2 3 4 5 6 7 8 9 10 11 12, parity_tslb, 25 26 27 28 29 30 31 32 33 34 35 36); end;
idx <- which(data$parity_cat==3)
length(idx)
data$parity_tslb[idx] <- data$tslb_cat[idx]+24
print("table(data$tslb_cat 25:36, parity_cat=3)")
table(data$tslb_cat[idx])
print("-------------------table of data$parity_tslb-------------------")
table(data$parity_tslb)
print("-------------------table of data$parity_tslb by study-------------------")
table(data$parity_tslb, data$study)

################################## Begin two-stage models ##################################################
############################################################################################################
############################################################################################################

##############collapse breastmost_cat 0 agefftp 0 and lastchildage 0 as 1, otherwise there will be collinearity issue with parity. Since women with no children will also fall into these categories
##############subset data? instead ofchange 0 to 1 etc. (FIXED: was removig cat 0)
# data$breastmos_cat[data$breastmos_cat==0] <- 1
# data$agefftp_cat[data$agefftp_cat==0] <- 1
# data$lastchildage_cat[data$lastchildage_cat==0] <- 1

##############only focus on population based study
data1 = data %>% filter(design_cat==0)

#############put the missing tumor characteristics as 888
idx.ER.mis <- which(data1$status==1&is.na(data1$ER_status1))
data1$ER_status1[idx.ER.mis] <- 888
idx.PR.mis <- which(data1$status==1&is.na(data1$PR_status1))
data1$PR_status1[idx.PR.mis] <- 888
idx.HER2.mis <- which(data1$status==1&is.na(data1$HER2_status1))
data1$HER2_status1[idx.HER2.mis] <- 888
idx.grade.mis <- which(data1$status==1&is.na(data1$Grade1))
data1$Grade1[idx.grade.mis] <- 888
#############put the subject BCAC-16687664 HER2 status as 1
idx <- which(data1$HER2_status1==2)
data1$HER2_status1[idx] <- 1
#############check the result
table(data1$status,data1$ER_status1)
table(data1$status,data1$PR_status1)
table(data1$status,data1$HER2_status1)
table(data1$status,data1$Grade1)
library(nnet)


########putting missing variable as NA, comment for now
# agemenarche_cat <- data1$agemenarche_cat
# idx <- which(agemenarche_cat==9)
# agemenarche_cat[idx] <- NA
# parity_cat <- data1$parity_cat
# idx <- which(parity_cat==9)
# parity_cat[idx] <- NA
# mensagelast_cat <- data1$mensagelast_cat
# idx <- which(mensagelast_cat==9)
# mensagelast_cat[idx] <- NA
# agefftp_cat <- data1$agefftp_cat
# idx <- which(agefftp_cat==9)
# agefftp_cat[idx] <- NA
# breastmos_cat <- data1$breastmos_cat
# idx <- which(breastmos_cat==9)
# breastmos_cat[idx] <- NA
# lastchildage_cat <- data1$lastchildage_cat
# idx <- which(lastchildage_cat==9)
# lastchildage_cat[idx] <- NA
# study <- data1$study
# refage <- data1$refage

all.covariates <- select(data1,c(
  
  "agemenarche_cat",
  "parity_cat",
  "tslb",
  "parity_tslb",
  #"parity_tslb_3",
  "mensagelast_cat",
  "agefftp_cat",
  "breastmos_cat",
  #"lastchildage_cat",
  "study",
  "StudyCountry",
  "refage"))


new.data1 <- data.frame(data1$molgroup,all.covariates)

colnames(new.data1) <- c("molgroup",
                         "agemenarche_cat",
                         "parity_cat",
                         "tslb",
                         "parity_tslb",
                         #"parity_tslb_3",
                         "mensagelast_cat",
                         "agefftp_cat",
                         "breastmos_cat",
                         #"lastchildage_cat",
                         "study",
                         "StudyCountry",
                         "refage")

sapply(new.data1, function(x) sum(is.na(x)))

####Categorize on new data after imputation

###########create the dummy variable for agemenarche_cat
###########create the dummy variable for parity_cat #parity_tslb?
###########create the dummy variable for mensagelat_cat
###########create the dummy variable for agefftp_cat
###########create the dummy variable for breastmos_cat
###########create the dummy variable for study

parity_tslb_mat <- model.matrix(~as.factor(new.data1$parity_tslb)-1)[,-1] 
colnames(parity_tslb_mat) <- paste0("parity_tslb",c(1:36)) 
print("parity_tslb_mat complete")

print("-------------------Parity_tslb_mat column names-------------------")
colnames(parity_tslb_mat)
print("-------------------table of data$parity_tslb-------------------")
table(data$parity_tslb)
agemenarche_mat <- model.matrix(~as.factor(new.data1$agemenarche_cat)-1)[,-1]
colnames(agemenarche_mat) <- paste0("agemenarche_cat",c(1:3))
print("agemenarche_mat (1:3) complete")

parity_mat <- model.matrix(~as.factor(new.data1$parity_cat)-1)[,-1:-2]
colnames(parity_mat) <- paste0("parity_cat",c(2:3))
colnames(parity_mat)

#################### * Cross tabulation parity_cat 1 and 2 and tslb continuous
print("----------Begin cross class of parity_cat 2 and 3 with tslb continuous for cont paritycat2 and 3-----")
parity_mat_2<-as.integer(parity_mat[1:119672,1])
parity_cat_2_tslb=as.integer(parity_mat_2*data1$tslb) # add parity_1_tslb column

parity_mat_3<-as.integer(parity_mat[1:119672,2])
parity_cat_3_tslb=as.integer(parity_mat_3*data1$tslb)  # add parity_2_tslb column
##############################################################################

mensagelast_mat <- model.matrix(~as.factor(new.data1$mensagelast_cat)-1)[,-1:-2]
colnames(mensagelast_mat) <- paste0("mensagelast_cat",c(2:3))
print("mensagelast_mat (2:3) complete")

agefftp_mat <- model.matrix(~as.factor(new.data1$agefftp_cat)-1)[,-1:-2]
colnames(agefftp_mat) <- paste0("agefftp",c(2:4)) # full= fill in TO DO
print("agefftp_mat (2:4) complete")

breastmos_mat <- model.matrix(~as.factor(new.data1$breastmos_cat)-1)[,-1:-2]
colnames(breastmos_mat) <- paste0("breastmos_cat",c(2:5)) # full= fill in TO DO
print("breastmos_mat (2:5) complete")

study_mat <- model.matrix(~as.factor(new.data1$study)-1)[,-1]

print(colnames(study_mat))
refage <- new.data1$refage
print("study_mat and refage complete")
study<-new.data1$study

###########create the phenotype file
#four different tumor characteristics were included,
#ER (positive vs negative),
#PR (positive vs negative),
#HER2 (positive vs negative)
#grade (oridinal 1,2,3)
#the phenotype file
y <- cbind(data1$status,data1$ER_status1,data1$PR_status1,
           data1$HER2_status1,data1$Grade1)
colnames(y) <- c("casecontrol",
                 "ER",
                 "PR",
                 "HER2",
                 "Grade")
idx <- which(data1$status==0)  #this was commented out. I added it back in but there are still issues
y[idx,2:5] <- NA

#by default, we remove all the subtypes with less than 10 cases
z.standard <- GenerateZstandard(y)


############ Generate the combinations of all the subtypes---------------------------------------------------
############ two-stage model with intrinsic subtypes
#
# # Check newer version zdesign bc3
z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
rowSums(z.design)

colnames(z.design) <- c("Luminal A-like","Luminal B,HER2-negative-like",
                        "Luminal B-like",
                        "HER2 enriched-like",
                        "TN")

z.ER <- cbind(1,z.standard[,1])

############# Define model paramteers---------------------------------------------------

parity_tslb = data1$parity_tslb #numeric
tslb = new.data1$tslb #numeric including 0.5
StudyCountry = data1$StudyCountry  #factor
refage = as.numeric(new.data1$refage)

#setwd("C:/Users/sandovall2/Box/BCAST Risk Factor/Datasets/Model Results/imputed datasets/imputation and model 10 29") #Set working directory

############# MODEL 1 & 2: We want a global test for heterogeneity between the variables and the tumor markers
############# Hoayu's 2nd model

############# generate the additive design matrix for two-stage model
z.standard <- GenerateZstandard(y)
z.additive <- cbind(1,z.standard)
z.baselineonly = z.additive[,1,drop=F]

#setwd("C:/Users/sandovall2/Box/BCAST Risk Factor/Datasets/Model Results/Model_1_and_2_parity_cross_tab/models_parity_cross_tab/Hoayu_Model_11_06/Model_1")
print("running model.1 <- TwoStageModel Case-case ORs")
model.1 <- EMmvpolySelfDesignnew(y=y,
                                 z.design = z.additive,
                                 x.self.design = cbind(
                                   parity_tslb_mat,
                                   #parity_cat_2_tslb,
                                   #parity_cat_3_tslb,
                                   agemenarche_mat,
                                   breastmos_mat,
                                   agefftp_mat,
                                   mensagelast_mat,
                                   #tslb,
                                   refage
                                 ),
                                 z.design2 = z.baselineonly,
                                 x.self.design2 = study_mat,
                                 missingTumorIndicator = 888)
#save(model.1,file=paste0("/data/sandovall2/Imputation7_m1/imputation7_m1_model1_",i1,".Rdata"))

#setwd("C:/Users/sandovall2/Box/BCAST Risk Factor/Datasets/Model Results/Model_1_and_2_parity_cross_tab/models_parity_cross_tab/Hoayu_Model_11_06/Model_2")
print("running model.2 <- EMmvpolySelfDesignnew Case-control ORs")
model.2 <- EMmvpolySelfDesignnew(y=y,
                                 z.design = z.design,
                                 x.self.design = cbind(
                                   parity_tslb_mat,
                                   #parity_cat_2_tslb,
                                   #parity_cat_3_tslb,
                                   agemenarche_mat, 
                                   breastmos_mat, 
                                   agefftp_mat,
                                   mensagelast_mat,
                                   #tslb,
                                   refage
                                 ),
                                 z.design2 = z.baselineonly,
                                 x.self.design2 = study_mat,
                                 missingTumorIndicator = 888)

#save(model.2,file=paste0("/data/sandovall2/Imputation7_m1/imputation7_m1_model2_",i1,".Rdata"))
