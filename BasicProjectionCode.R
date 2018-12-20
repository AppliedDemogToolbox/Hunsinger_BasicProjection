##############################################################################################################################
##############################################################################################################################
##R CODE FOR A BASIC POPULATION PROJECTION/FORECAST, APPLIED TO ALASKA
##
##EDDIE HUNSINGER, FEBRUARY 2010 (LAST UPDATED DECEMBER 2018)
##http://www.demog.berkeley.edu/~eddieh/
##edyhsgr@gmail.com
##
##THIS WORK IS BASED ON CODE FOR A STOCHASTIC POPULATION FORECAST FOR ALASKA (EDDIE HUNSINGER, 2010), WHICH IS AVAILABLE AT: 
##https://applieddemogtoolbox.github.io/Toolbox/#StochasticForecast
##TO SEE THE RELATED STOCHASTIC POPULATION FORECAST PAPER (WITH REFERENCES), GO TO http://www.demog.berkeley.edu/~eddieh/documents/ExpertForecastPaper.pdf
##
##IF YOU WOULD LIKE TO USE, SHARE OR REPRODUCE ANY INFORMATION OR IDEAS FROM THIS WORK, BE SURE TO CITE THE SOURCE
##
##MUCH OF THE STOCHASTIC POPULATION FORECAST CODE IS BASED ON CODE THAT I DEVELOPED WHILE WORKING FOR THE ALASKA DEPARTMENT OF LABOR AND WORKFORCE DEVELOPMENT
##I ALSO BENEFITED FROM FACULTY, STAFF AND FELLOW STUDENTS AT THE UC BERKELEY DEPARTMENT OF DEMOGRAPHY
##
##THERE IS NO WARRANTY FOR THIS CODE
##THIS CODE HAS NOT BEEN CAREFULLY REVIEWED
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
##############################################################################################################################
##SELECT THE INPUT DATA
##############################################################################################################################
##############################################################################################################################

##DIMENSIONS
##SIZE OF PROJECTION MATRIX
SIZE<-21
##NUMBER OF PROJECTION STEPS
STEPS<-6
BASEANDSTEPS<-STEPS+1

##SURVIVAL PARAMETERS
##
##YEAR 2000 US LIFE TABLE lx CURVES FROM THE NATIONAL CENTER FOR HEALTH STATISTICS
Survival<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_BasicProjection/raw/master/lx2000_US_NCHS.csv",header=TRUE,sep=",")
lxF<-Survival$F_2000
lxM<-Survival$M_2000
##"BA" IS THE BRASS RELATIONAL LOGIT MODEL ALPHA. IT IS CALIBRATED FOR THE JUMP-OFF PERIOD AND SET TO FOLLOW A PATH OF INCREASE
##OF .03 FOR EACH FIVE-YEAR STEP
BA_baseF<-.03
BA_baseM<-.08
BA_consF<-.03
BA_consM<-.03

##FERTILITY PARAMETERS
##
##YEAR 2005 US FERTILITY RATES FROM THE NATIONAL CENTER FOR HEALTH STATISTICS, SUMMED TO 1
Fertility<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_BasicProjection/raw/master/Fx2005_US_NCHS.csv",header=TRUE,sep=",")
PropFx<-c(Fertility$PropFx)
##FRACTION FEMALE AT BIRTH
ffab<-.4886
##"TFR" IS THE TOTAL FERTILITY RATE. IT IS SET TO 2.3 FOR EACH FORECAST STEP
TFR_base<-2.37
TFR_cons<-c(TFR_base,2.3,2.3,2.3,2.3,2.3,2.3)

##NETMIGRATION PARAMETERS
##
##1995-2000 MIGRATION PROFILES FROM THE 2000 CENSUS, SUMMED TO 1
Migration<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_BasicProjection/raw/master/MigProf2000_AKCACO_Census2000.csv",header=TRUE,sep=",")
PropInM<-c(Migration$AKIn)
PropInF<-c(Migration$AKIn)
PropOutM<-c(Migration$AKOut)
PropOutF<-c(Migration$AKOut)
##THE OUT MIGRATION NUMBER, WHICH IS SET TO 140,000 FOR EACH FIVE YEAR FORECAST STEP 
OutRate<-140000
##"TIR" IS THE NET NUMBER. IT IS SET TO 0 FOR EACH FORECAST STEP
TIR_base<-1601
TIR_cons<-c(TIR_base,0,0,0,0,0,0)

##BASE POPULATION
##
##VINTAGE 2009 POPULATION ESTIMATES FOR 2005, FROM THE ALASKA DEPARTMENT OF LABOR AND WORKFORCE DEVELOPMENT (http://almis.labor.state.ak.us/?PAGEID=67&SUBID=115)
K05<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_BasicProjection/raw/master/AgeSex2005_AK_AKDOLv2009.csv",header=TRUE,sep=",")
KF05<-K05$AK_F_2005
KM05<-K05$AK_M_2005

##############################################################################################################################
##############################################################################################################################
##RUN THE FORECAST CODE
##############################################################################################################################
##############################################################################################################################

##SURVIVAL
BAF<-array(0,c(BASEANDSTEPS))
BAF[1]<-BA_baseF
for(i in 2:BASEANDSTEPS){BAF[i]<-BAF[i-1]+BA_consF}
BAF<-t(BAF)
BrassF00<-data.frame(Alpha=BAF[1],Beta=1)
BrassF05<-data.frame(Alpha=BAF[2],Beta=1)
BrassF10<-data.frame(Alpha=BAF[3],Beta=1)
BrassF15<-data.frame(Alpha=BAF[4],Beta=1)
BrassF20<-data.frame(Alpha=BAF[5],Beta=1)
BrassF25<-data.frame(Alpha=BAF[6],Beta=1)
BrassF30<-data.frame(Alpha=BAF[7],Beta=1)
BAM<-array(0,c(BASEANDSTEPS))
BAM[1]<-BA_baseM
for(i in 2:BASEANDSTEPS){BAM[i]<-BAM[i-1]+BA_consM}
BAM<-t(BAM)
BrassM00<-data.frame(Alpha=BAM[1],Beta=1)
BrassM05<-data.frame(Alpha=BAM[2],Beta=1)
BrassM10<-data.frame(Alpha=BAM[3],Beta=1)
BrassM15<-data.frame(Alpha=BAM[4],Beta=1)
BrassM20<-data.frame(Alpha=BAM[5],Beta=1)
BrassM25<-data.frame(Alpha=BAM[6],Beta=1)
BrassM30<-data.frame(Alpha=BAM[7],Beta=1)

##FERTILITY
fmab <- 1-ffab
Fx<-array(0,c(SIZE))
Fx[1:SIZE]<-PropFx
TFRFORE<-array(0,c(BASEANDSTEPS))
TFRFORE[1]<-TFR_base
for(i in 2:BASEANDSTEPS){TFRFORE[i]<-TFR_cons[i]}
TFRFORE<-t(TFRFORE)

##NETMIGRATION
TIRFORE<-array(0,c(BASEANDSTEPS))
TIRFORE[1]<-TIR_base
for(i in 2:BASEANDSTEPS){TIRFORE[i]<-TIR_cons[i]}
TIRFORE<-t(TIRFORE)

##CALCULATE THE Yx FOR THE lx'S
YxM<-YxF<-NULL
for (i in 1:length(lxF)){YxF[i]<-.5*log(lxF[i]/(1-lxF[i]))}
for (i in 1:length(lxM)){YxM[i]<-.5*log(lxM[i]/(1-lxM[i]))}

##Improve Survival and make lx's for each period
lxF30<-lxF25<-lxF20<-lxF15<-lxF10<-lxF05<-lxF00<-array(0,c(SIZE+1))
lxM30<-lxM25<-lxM20<-lxM15<-lxM10<-lxM05<-lxM00<-array(0,c(SIZE+1))
for (i in 1:length(lxF)){lxF00[i]<-1/(1+exp(-2*BrassF00$Alpha-2*BrassF00$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM00[i]<-1/(1+exp(-2*BrassM00$Alpha-2*BrassM00$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF05[i]<-1/(1+exp(-2*BrassF05$Alpha-2*BrassF05$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM05[i]<-1/(1+exp(-2*BrassM05$Alpha-2*BrassM05$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF10[i]<-1/(1+exp(-2*BrassF10$Alpha-2*BrassF10$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM10[i]<-1/(1+exp(-2*BrassM10$Alpha-2*BrassM10$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF15[i]<-1/(1+exp(-2*BrassF15$Alpha-2*BrassF15$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM15[i]<-1/(1+exp(-2*BrassM15$Alpha-2*BrassM15$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF20[i]<-1/(1+exp(-2*BrassF20$Alpha-2*BrassF20$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM20[i]<-1/(1+exp(-2*BrassM20$Alpha-2*BrassM20$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF25[i]<-1/(1+exp(-2*BrassF25$Alpha-2*BrassF25$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM25[i]<-1/(1+exp(-2*BrassM25$Alpha-2*BrassM25$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF30[i]<-1/(1+exp(-2*BrassF30$Alpha-2*BrassF30$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM30[i]<-1/(1+exp(-2*BrassM30$Alpha-2*BrassM30$Beta*YxM[i]))}

##MAKE nLx's FOR EACH PERIOD
LxF30<-LxF25<-LxF20<-LxF25<-LxF15<-LxF10<-LxF05<-array(0,c(SIZE))
LxM30<-LxM25<-LxM20<-LxM15<-LxM15<-LxM10<-LxM05<-array(0,c(SIZE))
##**THIS IS A LITTLE OFF FOR THE FIRST AGE GROUP**
for (i in 1:SIZE){LxF05[i]<-.5*(lxF05[i]+lxF05[i+1])}
for (i in 1:SIZE){LxM05[i]<-.5*(lxM05[i]+lxM05[i+1])}
for (i in 1:SIZE){LxF10[i]<-.5*(lxF10[i]+lxF10[i+1])}
for (i in 1:SIZE){LxM10[i]<-.5*(lxM10[i]+lxM10[i+1])}
for (i in 1:SIZE){LxF15[i]<-.5*(lxF15[i]+lxF15[i+1])}
for (i in 1:SIZE){LxM15[i]<-.5*(lxM15[i]+lxM15[i+1])}
for (i in 1:SIZE){LxF20[i]<-.5*(lxF20[i]+lxF20[i+1])}
for (i in 1:SIZE){LxM20[i]<-.5*(lxM20[i]+lxM20[i+1])}
for (i in 1:SIZE){LxF25[i]<-.5*(lxF25[i]+lxF25[i+1])}
for (i in 1:SIZE){LxM25[i]<-.5*(lxM25[i]+lxM25[i+1])}
for (i in 1:SIZE){LxF30[i]<-.5*(lxF30[i]+lxF30[i+1])}
for (i in 1:SIZE){LxM30[i]<-.5*(lxM30[i]+lxM30[i+1])}

##TABLE e0
e0MFORE<-array(0,c(BASEANDSTEPS))
e0MFORE[2]<-sum(LxM05*5)
e0MFORE[3]<-sum(LxM10*5)
e0MFORE[4]<-sum(LxM15*5)
e0MFORE[5]<-sum(LxM20*5)
e0MFORE[6]<-sum(LxM25*5)
e0MFORE[7]<-sum(LxM30*5)
e0FFORE<-array(0,c(BASEANDSTEPS))
e0FFORE[2]<-sum(LxF05*5)
e0FFORE[3]<-sum(LxF10*5)
e0FFORE[4]<-sum(LxF15*5)
e0FFORE[5]<-sum(LxF20*5)
e0FFORE[6]<-sum(LxF25*5)
e0FFORE[7]<-sum(LxF30*5)

##MAKE nSx's FOR EACH PERIOD
SxF30<-SxF25<-SxF20<-SxF25<-SxF15<-SxF10<-SxF05<-array(0,c(SIZE-1))
SxM30<-SxM25<-SxM20<-SxM15<-SxM15<-SxM10<-SxM05<-array(0,c(SIZE-1))
for (i in 1:SIZE-1){SxF05[i]<-(LxF05[i+1]/LxF05[i])}
for (i in 1:SIZE-1){SxM05[i]<-(LxM05[i+1]/LxM05[i])}
for (i in 1:SIZE-1){SxF10[i]<-(LxF10[i+1]/LxF10[i])}
for (i in 1:SIZE-1){SxM10[i]<-(LxM10[i+1]/LxM10[i])}
for (i in 1:SIZE-1){SxF15[i]<-(LxF15[i+1]/LxF15[i])}
for (i in 1:SIZE-1){SxM15[i]<-(LxM15[i+1]/LxM15[i])}
for (i in 1:SIZE-1){SxF20[i]<-(LxF20[i+1]/LxF20[i])}
for (i in 1:SIZE-1){SxM20[i]<-(LxM20[i+1]/LxM20[i])}
for (i in 1:SIZE-1){SxF25[i]<-(LxF25[i+1]/LxF25[i])}
for (i in 1:SIZE-1){SxM25[i]<-(LxM25[i+1]/LxM25[i])}
for (i in 1:SIZE-1){SxF30[i]<-(LxF30[i+1]/LxF30[i])}
for (i in 1:SIZE-1){SxM30[i]<-(LxM30[i+1]/LxM30[i])}

##PUT THE Sx DATA INTO THE SUBDIAGONAL OF WHAT WILL BE THE LESLIE MATRICES
SF30<-SF25<-SF20<-SF15<-SF10<-SF05<-array(0,c(SIZE,SIZE))
SM30<-SM25<-SM20<-SM15<-SM10<-SM05<-array(0,c(SIZE,SIZE))
SF05<-rbind(0,cbind(diag(SxF05),0))
SF10<-rbind(0,cbind(diag(SxF10),0))
SF15<-rbind(0,cbind(diag(SxF15),0))
SF20<-rbind(0,cbind(diag(SxF20),0))
SF25<-rbind(0,cbind(diag(SxF25),0))
SF30<-rbind(0,cbind(diag(SxF30),0))
SM05<-rbind(0,cbind(diag(SxM05),0))
SM10<-rbind(0,cbind(diag(SxM10),0))
SM15<-rbind(0,cbind(diag(SxM15),0))
SM20<-rbind(0,cbind(diag(SxM20),0))
SM25<-rbind(0,cbind(diag(SxM25),0))
SM30<-rbind(0,cbind(diag(SxM30),0))

##PUT FERTILITY INTO AGE PROFILES
TFR2005<-TFRFORE[2]
TFR2010<-TFRFORE[3]
TFR2015<-TFRFORE[4]
TFR2020<-TFRFORE[5]
TFR2025<-TFRFORE[6]
TFR2030<-TFRFORE[7]

Fert2005<-TFR2005*Fx
Fert2010<-TFR2010*Fx
Fert2015<-TFR2015*Fx
Fert2020<-TFR2020*Fx
Fert2025<-TFR2025*Fx
Fert2030<-TFR2030*Fx

##MAKE MIGRATION AGE PROFILES
IxM<-array(0,c(SIZE))
IxM[1:SIZE]<-PropInM
IxF<-array(0,c(SIZE))
IxF[1:SIZE]<-PropInF

OxM<-array(0,c(SIZE))
OxM[1:SIZE]<-PropOutM
OxF<-array(0,c(SIZE))
OxF[1:SIZE]<-PropOutF

##MAKE THE LESLIE MATRICES FOR FEMALES
BF30<-BF25<-BF20<-BF15<-BF10<-BF05<-0*SF05

for(j in 1:SIZE-1)
  {BF05[1,j]<-(LxF05[1]/2)*(Fert2005[j]+Fert2005[j+1]*(SxF05[j]))*ffab}
AF05 = SF05 + BF05
for(j in 1:SIZE-1)
  {BF10[1,j]<-(LxF10[1]/2)*(Fert2010[j]+Fert2010[j+1]*(SxF10[j]))*ffab}
AF10 = SF10 + BF10
for(j in 1:SIZE-1)
  {BF15[1,j]<-(LxF15[1]/2)*(Fert2015[j]+Fert2015[j+1]*(SxF15[j]))*ffab}
AF15 = SF15 + BF15
for(j in 1:SIZE-1)
  {BF20[1,j]<-(LxF20[1]/2)*(Fert2020[j]+Fert2020[j+1]*(SxF20[j]))*ffab}
AF20 = SF20 + BF20
for(j in 1:SIZE-1)
  {BF25[1,j]<-(LxF25[1]/2)*(Fert2025[j]+Fert2025[j+1]*(SxF25[j]))*ffab}
AF25 = SF25 + BF25
for(j in 1:SIZE-1)
  {BF30[1,j]<-(LxF30[1]/2)*(Fert2030[j]+Fert2030[j+1]*(SxF30[j]))*ffab}
AF30 = SF30 + BF30

##MAKE ARRAYS TO HOLD THE DATA
KF05<-array(KF05,c(SIZE,1))
KF10<-array(0,c(SIZE,1))
KF35<-KF30<-KF25<-KF20<-KF15<-KF10

##PROJECT THE FEMALE POPULATION (NATURAL INCREASE, LESS OUT MIGRATION, PLUS IN MIGRATION 
##(IN MIGRATION IS OUT MIGRATION SUM PLUS NET MIGRATION SUM)
Out2005<-array(0,c(SIZE,1))
Out2005<-OxF*OutRate/2
In2005<-array(0,c(SIZE,1))
In2005<-(sum(Out2005)+TIRFORE[2]/2)*IxF
KF10<-(AF05%*%KF05)+t(t(In2005))-t(t(Out2005))

Out2010<-array(0,c(SIZE,1))
Out2010<-OxF*OutRate/2
In2010<-array(0,c(SIZE,1))
In2010<-(sum(Out2010)+TIRFORE[3]/2)*IxF
KF15<-(AF10%*%KF10)+t(t(In2010))-t(t(Out2010))

Out2015<-array(0,c(SIZE,1))
Out2015<-OxF*OutRate/2
In2015<-array(0,c(SIZE,1))
In2015<-(sum(Out2015)+TIRFORE[4]/2)*IxF
KF20<-(AF15%*%KF15)+t(t(In2015))-t(t(Out2015))

Out2020<-array(0,c(SIZE,1))
Out2020<-OxF*OutRate/2
In2020<-array(0,c(SIZE,1))
In2020<-(sum(Out2020)+TIRFORE[5]/2)*IxF
KF25<-(AF20%*%KF20)+t(t(In2020))-t(t(Out2020))

Out2025<-array(0,c(SIZE,1))
Out2025<-OxF*OutRate/2
In2025<-array(0,c(SIZE,1))
In2025<-(sum(Out2025)+TIRFORE[6]/2)*IxF
KF30<-(AF25%*%KF25)+t(t(In2025))-t(t(Out2025))

Out2030<-array(0,c(SIZE,1))
Out2030<-OxF*OutRate/2
In2030<-array(0,c(SIZE,1))
In2030<-(sum(Out2030)+TIRFORE[7]/2)*IxF
KF35<-(AF30%*%KF30)+t(t(In2030))-t(t(Out2030))

##MAKE THE LESLIE MATRICES FOR MALES
BM30<-BM25<-BM20<-BM15<-BM10<-BM05<-0*SF05

for(j in 1:SIZE-1)
  {BM05[1,j]<-(LxM05[1]/2)*(Fert2005[j]+Fert2005[j+1]*(SxF05[j]))*fmab}
AM05 = SM05 + BM05
for(j in 1:SIZE-1)
  {BM10[1,j]<-(LxM10[1]/2)*(Fert2010[j]+Fert2010[j+1]*(SxF10[j]))*fmab}
AM10 = SM10 + BM10
for(j in 1:SIZE-1)
  {BM15[1,j]<-(LxM15[1]/2)*(Fert2015[j]+Fert2015[j+1]*(SxF15[j]))*fmab}
AM15 = SM15 + BM15
for(j in 1:SIZE-1)
  {BM20[1,j]<-(LxM20[1]/2)*(Fert2020[j]+Fert2020[j+1]*(SxF20[j]))*fmab}
AM20 = SM20 + BM20
for(j in 1:SIZE-1)
  {BM25[1,j]<-(LxM25[1]/2)*(Fert2025[j]+Fert2025[j+1]*(SxF25[j]))*fmab}
AM25 = SM25 + BM25
for(j in 1:SIZE-1)
  {BM30[1,j]<-(LxM30[1]/2)*(Fert2030[j]+Fert2030[j+1]*(SxF30[j]))*fmab}
AM30 = SM30 + BM30

##MAKE ARRAYS TO HOLD THE DATA
KM05<-array(KM05,c(SIZE,1))
KM10<-array(0,c(SIZE,1))
KBM10<-array(0,c(SIZE,1))
KSM10<-array(0,c(SIZE,1))
KM35<-KM30<-KM25<-KM20<-KM15<-KM10
KBM35<-KBM30<-KBM25<-KBM20<-KBM15<-KBM10
KSM35<-KSM30<-KSM25<-KSM20<-KSM15<-KSM10

##PROJECT THE MALE POPULATION (NATURAL INCREASE, LESS OUT MIGRATION, PLUS IN MIGRATION 
##(IN MIGRATION IS OUT MIGRATION SUM PLUS NET MIGRATION SUM)
OutM2005<-array(0,c(SIZE))
OutM2005<-OxM*OutRate/2
InM2005<-array(0,c(SIZE))
InM2005<-(sum(OutM2005)+TIRFORE[2]/2)*IxM
KBM10<-(BM05%*%KF05)
KSM10<-(SM05%*%KM05)+t(t(InM2005))-t(t(OutM2005))
KM10<-(KBM10+KSM10)

OutM2010<-array(0,c(SIZE))
OutM2010<-OxM*OutRate/2
InM2010<-array(0,c(SIZE))
InM2010<-(sum(OutM2010)+TIRFORE[3]/2)*IxM
KBM15<-(BM10%*%KF10)
KSM15<-(SM10%*%KM10)+t(t(InM2010))-t(t(OutM2010))
KM15<-(KBM15+KSM15)

OutM2015<-array(0,c(SIZE))
OutM2015<-OxM*OutRate/2
InM2015<-array(0,c(SIZE))
InM2015<-(sum(OutM2015)+TIRFORE[4]/2)*IxM
KBM20<-(BM15%*%KF15)
KSM20<-(SM15%*%KM15)+t(t(InM2015))-t(t(OutM2015))
KM20<-(KBM20+KSM20)

OutM2020<-array(0,c(SIZE))
OutM2020<-OxM*OutRate/2
InM2020<-array(0,c(SIZE))
InM2020<-(sum(OutM2020)+TIRFORE[5]/2)*IxM
KBM25<-(BM20%*%KF20)
KSM25<-(SM20%*%KM20)+t(t(InM2020))-t(t(OutM2020))
KM25<-(KBM25+KSM25)

OutM2025<-array(0,c(SIZE))
OutM2025<-OxM*OutRate/2
InM2025<-array(0,c(SIZE))
InM2025<-(sum(OutM2025)+TIRFORE[6]/2)*IxM
KBM30<-(BM25%*%KF25)
KSM30<-(SM25%*%KM25)+t(t(InM2025))-t(t(OutM2025))
KM30<-(KBM30+KSM30)

OutM2030<-array(0,c(SIZE))
OutM2030<-OxM*OutRate/2
InM2030<-array(0,c(SIZE))
InM2030<-(sum(OutM2030)+TIRFORE[7]/2)*IxM
KBM35<-(BM30%*%KF30)
KSM35<-(SM30%*%KM30)+t(t(InM2030))-t(t(OutM2030))
KM35<-(KBM35+KSM35)

##MAKE TABLES OF DATA OF INTEREST
KT05<-sum(KF05)+sum(KM05)
KT10<-sum(KF10)+sum(KM10)
KT15<-sum(KF15)+sum(KM15)
KT20<-sum(KF20)+sum(KM20)
KT25<-sum(KF25)+sum(KM25)
KT30<-sum(KF30)+sum(KM30)
KT35<-sum(KF35)+sum(KM35)
KT<-c(KT05,KT10,KT15,KT20,KT25,KT30,KT35)

##############################################################################################################################
##############################################################################################################################
##FIGURES AND TABLES
##############################################################################################################################
##############################################################################################################################

#FIGURE 1
TFR<-array(0,c(13))
TFR<-c(2.58,2.16,2.40,2.35,2.58,2.48,TFRFORE)
plot(TFR,type="l",ylim=c(1,3),xlim=c(0,13),col="black",xlab="Time Period",ylab="Total Fertility Rate",axes=F,lwd=8)
axis(side=1,at=1:13,labels=c("1970-75","1975-80","1980-85","1985-90","1990-95","1995-00","2000-05","2005-10","2010-15","2015-20","2020-25","2025-30","2030-35"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
title("TOTAL FERTILITY RATE: HISTORICAL AND FORECAST",cex.main=1)
Sys.sleep(5)

#FIGURE 2
TNRFORE<-array(0,c(13))
TNRFORE<-c(47469,881,75984,-39444,3942,-10674,1601,TIRFORE)
plot(TNRFORE,type="l",ylim=c(-100000,100000),xlim=c(0,13),col="black",xlab="Time Period",ylab="Net Migration",axes=F,lwd=8)
axis(side=1,at=1:13,labels=c("1970-75","1975-80","1980-85","1985-90","1990-95","1995-00","2000-05","2005-10","2010-15","2015-20","2020-25","2025-30","2030-35"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
title("NET MIGRATION: HISTORICAL AND FORECAST",cex.main=1)
Sys.sleep(5)

#FIGURE 4
e0F<-array(0,c(13))
e0F<-c(74.6,75.8,77.0,78.1,78.9,79.4,79.9,e0FFORE[2:7])
plot(e0F,type="l",ylim=c(65,90),xlim=c(0,13),col="black",xlab="Time Period",ylab="Life Expectancy at Birth",axes=F,lwd=8)
axis(side=1,at=1:13,labels=c("1970-75","1975-80","1980-85","1985-90","1990-95","1995-00","2000-05","2005-10","2010-15","2015-20","2020-25","2025-30","2030-35"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
e0M<-array(0,c(13))
e0M<-c(66.8,68.1,69.5,70.9,72.4,74.0,75.2,e0MFORE[2:7])
points(e0M,type="l",ylim=c(65,80),xlim=c(0,13),col="black",lwd=8)
title("LIFE EXPECTANCY AT BIRTH: HISTORICAL AND FORECAST",cex.main=1)
mtext(side=3,line=.75,text="(Historical estimates are interpolated by the author from estimates for 1970, 1980, 1990 and 2000.)",font=1,cex=.75)
mtext(side=1,line=-16.5,text="Female",font=2,cex=1)
mtext(side=1,line=-10,adj=.6,text="Male",font=2,cex=1)
Sys.sleep(5)

#FIGURE 5
KTH<-array(0,c(14))
KTH<-c(308500,384100,419800,543900,553171,601581,627533,664334,KT[2:7])
plot(KTH,type="l",ylim=c(0,KT35*1.5),xlim=c(0,14),col="black",xlab="Year",ylab="Population",axes=F,lwd=8)
axis(side=1,at=1:14,labels=c("1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020","2025","2030","2035"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
title("TOTAL POPULATION: HISTORICAL AND FORECAST",cex.main=1)
Sys.sleep(5)

#FIGURE 6 (CARL MASON'S GREAT PYRAMID FUNCTION)
poppyr3<-function(male,female,cat){
split.screen(figs=rbind(c(0,.58,0,1),c(.43,1,0,1)))
screen(1)
barplot(male,horiz=T,names=cat,space=0,
xlim=c(50000,0),col="dodgerblue")
title("Male",line=-3,cex.main=1)
screen(2)
barplot(female,horiz=T,names=F,space=0,
xlim=c(0,50000),col="gold")
title("Female",line=-3,cex.main=1)
close.screen(all=T)}
male<-array(0,c(SIZE))
female<-array(0,c(SIZE))
for (i in 1:SIZE){male[i]<-KM05[i,]}
for (i in 1:SIZE){female[i]<-KF05[i,]}
cat<-seq(0,100,5)
poppyr3(male,female,cat)
mtext(side=3,line=1,text="ALASKA 2005 POPULATION FORECAST BY AGE AND SEX",font=2,cex=1)
mtext(side=2,line=3,text="Age",font=1,cex=1)
mtext(side=1,line=3,text="Population",font=1,cex=1)
#dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure13.eps")
Sys.sleep(5)

#FIGURE 7 (CARL MASON'S GREAT PYRAMID FUNCTION)
poppyr3<-function(male,female,cat){
split.screen(figs=rbind(c(0,.58,0,1),c(.43,1,0,1)))
screen(1)
barplot(male,horiz=T,names=cat,space=0,
xlim=c(50000,0),col="dodgerblue")
title("Male",line=-3,cex.main=1)
screen(2)
barplot(female,horiz=T,names=F,space=0,
xlim=c(0,50000),col="gold")
title("Female",line=-3,cex.main=1)
close.screen(all=T)}
male<-array(0,c(SIZE))
female<-array(0,c(SIZE))
for (i in 1:SIZE){male[i]<-KM35[i,]}
for (i in 1:SIZE){female[i]<-KF35[i,]}
cat<-seq(0,100,5)
poppyr3(male,female,cat)
mtext(side=3,line=1,text="ALASKA 2035 POPULATION FORECAST BY AGE AND SEX",font=2,cex=1)
mtext(side=2,line=3,text="Age",font=1,cex=1)
mtext(side=1,line=3,text="Population",font=1,cex=1)
#dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure13.eps")
Sys.sleep(5)

KTH

#write.table(###, file="G:/###/###.csv", sep=",")

