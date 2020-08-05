#' Runs NicheMapR model
#'
#' @description Wrapper function for running NicheMapR to derive soil moistures
#' @param weather a data.frame of hourly weather variables in same format as [microclimc::weather()].
#' Column obs_time must give times in GMT/UTC.
#' @param prec a vector of daily or hourly recipitation values (mm).
#' @param lat latitude (decimal degrees, positive in northern hemisphere).
#' @param long longitude (decimal degrees, negative west of Greenwich meridian).
#' @param Usrhyt Local height (m) at which air temperature, wind speed and humidity are to be computed (see details).
#' @param Veghyt At of vegetation canopy (see details).
#' @param Refhyt Reference height (m) at which air input climate variables are measured.
#' @param PAI single value or vector of hourly or daily values of total plant area per unit ground area of canopy.
#' @param LOR Campbell leaf angle distrubution coefficient
#' @param pLAI fraction of PAI that is green vegetation.
#' @param clump clumpiness factor for canopy (0 to 1, 0 = even)
#' @param REFL A single numeric value of ground reflectivity of shortwave radiation.
#' @param LREFL A single numeric value of average canopy reflectivity of shortwave radiation.
#' @param SLE Thermal emissvity of ground
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values
#' starting from 0, and more closely spaced near the surface.
#' @param ALTT elevation (m) of location.
#' @param SLOPE slope of location (decimal degrees).
#' @param ASPECT aspect of location (decimal degrees, 0 = north).
#' @param ERR Integrator error tolerance for soil temperature calculations.
#' @param soiltype soil type as for [microclimc::soilparams()]. Set to NA to use soil parameters
#' specified below.
#' @param PE Air entry potentials (J/kg) (19 values descending through soil for specified soil
#' nodes in parameter DEP and points half way between). Ignored if soil type is given.
#' @param KS Saturated conductivities (kg s/m^3) (19 values descending through soil for
#' specified soil nodes in parameter DEP and points half way between). Ignored if soil type
#' is given.
#' @param BB  Campbell's soil 'b' parameter (19 values descending through soil for specified
#' soil nodes in parameter DEP and points half way between). Ignored if soil type is given.
#' @param BD Soil bulk density (Mg/m^3) (19 values descending through soil for specified
#' soil nodes in parameter DEP and points half way between). Ignored if soil type is given.
#' @param DD Soil density (Mg/m^3) (19 values descending through soil for specified
#' soil nodes in parameter DEP and points half way between). Ignored if soil type is given.
#' @param cap Is organic cap present on soil surface? (1 = Yes, 0 = No).
#' @param hori Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 10 degree intervals.
#' @param maxpool Max depth for water pooling on the surface (mm), to account for runoff.
#' @param rainmult Rain multiplier for surface soil moisture (used to induce runon).
#' @param SoilMoist_Init Initial volumetric soil water content at each soil node (m^3/m^3)
#' @return a list with the following components:
#' @return `metout` a data.frame of microclimate variables for each hour at height `Usrhyt`.
#' @return `soiltemps` a data.frame of hourly soil temperatures for each depth node.
#' @return `soilmoist` a data.frame of hourly soil volumetric water content for each depth node.
#' @return `snowtemp` if snow present, a data,frame of snow temperature (°C), at each of the potential 8 snow layers (see details).
#' 0 if snow not present.
#' @return `plant` a data.frame of plant transpiration rates (g/m^2/hr), leaf water potentials
#' (J/kg) and root water potential (J/kg) at each of the 10 specified depths.
#' @details Requires NicheMapR: devtools::install_github('mrke/NicheMapR'). NicheMapR is an integrated
#' soil moisture and temperature model that treats the vegetation as a single layer (https://mrke.github.io/).
#' This wrapper function, runs the NicheMapR::microclimate function with reduced parameter
#' inputs, by specifying sensible default values for other parameters where it is unlikely that these
#' would be known or explicitely measured at the site. The degree of canopy shading is worked out
#' explicitely from `PAI` at each hourly time interval by estimating canopy tranmission of direct
#' and diffuse radiation, and by adjusting the amount of radiation that would be absorbed
#' at the ground surface for a given slope and aspect. Roughness lengths and zero plane displacement
#' heights are limited to <2m in NicheMapR, so if `Veghyt` > 2, it is set to 2m.
#' This is adequate for computing soil moisture and total evapotranspiration, but will give false
#' air temperatures for heights below canopy (see [runwithNMR()]. If `soiltype` is given, the subsequent
#' soil parameters are computed from soil type.  If `soiltype` is NA, soil paramaters must be specified.
#'
#' @import microctools
#' @import NicheMapR
#' @export
#' @examples
#' # Run NicheMapR with default parameters and inbuilt weather datasets
#' library(NicheMapR)
#' microout<-runNMR(weather, dailyprecip, 50.2178, -5.32656, 0.05, 0.02, PAI = 1)
#' # Plot air temperatures
#' metout <- microout$metout
#' tmn <- min(metout$TALOC,metout$TAREF)
#' tmx <- max(metout$TALOC,metout$TAREF)
#' dday <- metout$DOY+metout$TIME/1440 # Decimal day
#' plot(metout$TALOC~dday, type="l", col = "red", ylim=c(tmn,tmx), ylab = "Temperature")
#' par(new = T)
#' plot(metout$TAREF~dday, type="l", col = "blue", ylim=c(tmn,tmx), ylab = "", xlab = "")
#' # Plot soil temperatures
#' soiltemp<-microout$soiltemps
#' tmn <- min(soiltemp$D0cm,soiltemp$D30cm)
#' tmx <- max(soiltemp$D0cm,soiltemp$D30cm)
#' plot(soiltemp$D0cm~dday, type="l", col = "red", ylim=c(tmn,tmx), ylab = "Temperature")
#' par(new = T)
#' plot(soiltemp$D30cm~dday, type="l", col = "blue", ylim=c(tmn,tmx), ylab = "", xlab = "")
#' # Plot soil moistures
#' soilm <- microout$soilmoist
#' mmn <- min(soilm$WC2.5cm,soilm$WC30cm)
#' mmx <- max(soilm$WC2.5cm,soilm$WC30cm)
#' plot(soilm$WC2.5cm~dday, type="l", col = "red", ylim=c(mmn,mmx), ylab = "Soil moisture")
#' par(new = T)
#' plot(soilm$WC30cm~dday, type="l", col = "blue", ylim=c(mmn,mmx), ylab = "", xlab = "")
runNMR <- function(weather, prec, lat, long, Usrhyt, Veghyt, Refhyt = 2, PAI = 3, LOR = 1,
                   pLAI = 0.8, clump = 0, REFL = 0.15, LREFL = 0.4, SLE = 0.95,
                   DEP = c(0,2.5,5,10,15,20,30,50,100,200), ALTT = 0, SLOPE = 0, ASPECT = 0,
                   ERR = 1.5, soiltype = "Loam", PE=rep(1.1,19), KS=rep(0.0037,19), BB=rep(4.5,19),
                   BD=rep(1.3,19), DD=rep(2.65,19), cap = 1, hori = rep(0,36), maxpool = 1000,
                   rainmult = 1, SoilMoist_Init = c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3)) {
  if (Veghyt > 2) Veghyt<-2
  loc<-c(long,lat)
  tmehr<-as.POSIXlt(weather$obs_time,tz="UTC")
  nyears<-length(unique(tmehr$year))
  fail<-nyears*24*365
  ystart<-tmehr$year[1]+1900
  yfinish<-tmehr$year[length(tmehr)]+1900
  yearlist<-seq(ystart,(ystart+(nyears-1)),1)
  tme<-seq(tmehr[1],tmehr[length(tmehr)],"days")
  doy <- as.numeric(strftime(tme, format = "%j"))
  ndays<-length(doy)
  doynum<-ndays
  ida<-ndays
  microdaily<-1
  daystart<-1
  # Set root properties
  L<-c(0,0,8.2,8,7.8,7.4,7.1,6.4,5.8,4.8,4,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000
  R1<-0.001
  RW<-2.5e+10
  RL<-2e+06
  PC<- -1500
  SP<-10
  IM<-1e-06
  # Set snow properties
  snowtemp<-1.5
  snowdens<-0.375
  densfun<-c(0.5979,0.2178,0.001,0.0038)
  snowmelt<-1
  undercatch<-1
  rainmelt<-0.0125
  grasshade<-ifelse(Veghyt<0.5,0,1)
  ### LAI etc
  if (length(PAI)==1) {
    PAI<-rep(PAI,ndays)
  } else if (length(PAI)==ndays*24) {
    PAI<-matrix(PAI,ncol=24,byrow=T)
    PAI<-apply(PAI,1,mean)
  }
  if (length(PAI) != ndays) stop("PAI must be a single value or hourly/daily values")
  MAXSHADES<-rep(100,ndays)
  # Work out canopy shading
  diftr<-cantransdif(PAI,LREFL,clump)
  lt<-tmehr$hour+tmehr$min/60+tmehr$sec/3600
  jd<-jday(tme=tmehr)
  sa<-solalt(lt,lat,long,jd,0)
  dirtr<-cantransdir(PAI,LOR,sa,LREFL,clump)
  si<-solarcoef(SLOPE,ASPECT,lt,lat,long,jd,merid=0)
  rad_dir<-weather$swrad-weather$difrad
  rad_ground<-si*dirtr*rad_dir+weather$difrad*diftr
  MINSHADES<-(1-rad_ground/weather$swrad)*100
  MINSHADES[is.na(MINSHADES)]<-mean(MINSHADES,na.rm=T)
  MINSHADES[MINSHADES>99.9]<-99.9
  MINSHADES<-matrix(MINSHADES,ncol=24,byrow=T)
  MINSHADES<-apply(MINSHADES,1,mean)
  intercept<-mean(MINSHADES)/100*0.3 # snow interception
  x<-t(as.matrix(as.numeric(c(loc[1],loc[2]))))
  ALREF<-abs(trunc(x[1]))
  HEMIS<-ifelse(x[2]<0,2,1)
  ALAT<-abs(trunc(x[2]))
  AMINUT<-(abs(x[2])-ALAT)*60
  ALONG<- abs(trunc(x[1]))
  ALMINT<-(abs(x[1])-ALONG)*60
  azmuth<-ASPECT
  lat<-as.numeric(loc[2])
  long<-as.numeric(loc[1])
  Density<-2.56
  Thcond<-2.5
  SpecHeat<-870
  if (is.na(soiltype) == FALSE) {
    sel<-which(microclimc::soilparams$Soil.type == soiltype)
    if (length(sel) == 0) stop("Erroneous soil type specified")
    PE<-rep(microclimc::soilparams$psi_e[sel],19)
    BB<-rep(microclimc::soilparams$b[sel],19)
    BD<-rep(microclimc::soilparams$rho[sel],19)
    KS<-rep(CampNormTbl9_1$Ks[sel],19)
  }
  BulkDensity<-BD[seq(1,19,2)]
  VIEWF<-1-sum(sin(as.data.frame(hori)*pi/180))/length(hori)
  #########################
  lt<-tmehr$hour+tmehr$min/60+tmehr$sec/3600
  jd<-jday(tme=tmehr)
  sa<-solalt(lt,lat,long,jd,0)
  ZENhr<-90-sa
  ZENhr[ZENhr>90]<-90
  TAIRhr<-weather$temp
  SOLRhr<-weather$swrad*VIEWF
  SOLRhr[SOLRhr<0]<-0
  sb<-5.67*10^-8
  IRDhr<-weather$skyem*sb*(weather$temp+273.15)^4
  RHhr<-weather$relhum
  RHhr[RHhr>100]<-100
  RHhr[RHhr<0]<-0
  e0<-satvap(TAIRhr,ice = TRUE)
  ea<-e0*(RHhr/100)
  eo<-1.24*(10*ea/(TAIRhr+273.15))^(1/7)
  CLDhr<-((weather$skyem-eo)/(1-eo))*100
  CLDhr[CLDhr<0]<-0
  CLDhr[CLDhr>100]<-100
  WNhr<-weather$windspeed
  WNhr[is.na(WNhr)]<-0.1
  PRESShr<-weather$pressure*1000
  RAINFALL<-prec
  RAINFALL[RAINFALL<0.1]<-0
  ZENhr2<-ZENhr
  ZENhr2[ZENhr2!=90]<-0
  dmaxmin<-function(x,fun) {
    dx <- t(matrix(x, nrow = 24))
    apply(dx, 1, fun)
  }
  TMAXX<-dmaxmin(TAIRhr,max)
  TMINN<-dmaxmin(TAIRhr,min)
  CCMAXX<-dmaxmin(CLDhr,max)
  CCMINN<-dmaxmin(CLDhr,min)
  RHMAXX<-dmaxmin(RHhr,max)
  RHMINN<-dmaxmin(RHhr,min)
  WNMAXX<-dmaxmin(WNhr,max)
  WNMINN<-dmaxmin(WNhr,min)
  PRESS<-dmaxmin(PRESShr,min)
  ###
  slope <- 0
  azmuth <- 0
  relhum <- 1
  optdep.summer<-as.data.frame(rungads(loc[2],loc[1],relhum, 0))
  optdep.winter<-as.data.frame(rungads(loc[2],loc[1],relhum, 1))
  optdep<-cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
  optdep<-as.data.frame(optdep)
  colnames(optdep)<-c("LAMBDA","OPTDEPTH")
  a<-lm(OPTDEPTH~poly(LAMBDA,6,raw=TRUE),data=optdep)
  LAMBDA<-c(290,295,300,305,310,315,320,330,340,350,360,370,380,390,400,420,440,460,480,
            500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,
            880,900,920,940,960,980,1000,1020,1080,1100,1120,1140,1160,1180,1200,1220,
            1240,1260,1280,1300,1320,1380,1400,1420,1440,1460,1480,1500,1540,1580,1600,
            1620,1640,1660,1700,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,
            2150,2200,2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,
            2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000)
  TAI<-predict(a,data.frame(LAMBDA))
  RAINFALL<-RAINFALL
  ALLMINTEMPS<-TMINN
  ALLMAXTEMPS<-TMAXX
  ALLTEMPS<-cbind(ALLMAXTEMPS,ALLMINTEMPS)
  WNMAXX<-WNMAXX
  WNMINN<-WNMINN
  WNhr<-WNhr
  REFLS<-rep(REFL,ndays)
  PCTWET <-rep(0,ndays)
  soilwet<-RAINFALL
  soilwet[soilwet<=1.5]<-0
  soilwet[soilwet>0]<-90
  if (ndays < 1) PCTWET<-pmax(soilwet,PCTWET)
  Intrvls<-rep(0,ndays)
  Intrvls[1]<-1
  Numtyps<-10
  Nodes<-matrix(data=0,nrow=10,ncol=ndays)
  Nodes[1:10,]<-c(1:10)
  ALREF<-abs(trunc(x[1]))
  HEMIS<-ifelse(x[2]<0,2,1)
  ALAT<-abs(trunc(x[2]))
  AMINUT<-(abs(x[2])-ALAT)*60
  ALONG<-abs(trunc(x[1]))
  ALMINT<-(abs(x[1])-ALONG)*60
  avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
  soilinit<-rep(avetemp,20)
  tannul<-mean(unlist(ALLTEMPS))
  deepsoil<-rep(mean(TAIRhr),ndays)
  SLES<-matrix(nrow=ndays,data=0)
  SLES<-SLES+SLE
  moists2<-matrix(nrow=10,ncol=ndays,data=0)
  moists2[1:10,]<-SoilMoist_Init
  moists<-moists2
  soilprops<-matrix(data=0,nrow=10,ncol=5)
  soilprops[,1]<-BulkDensity
  soilprops[,2]<-min(0.26,1-BulkDensity/Density)
  soilprops[,3]<-Thcond
  soilprops[,4]<-SpecHeat
  soilprops[,5]<-Density
  if (cap==1) {
    soilprops[1:2,3] <- 0.2
    soilprops[1:2,4] <- 1920
  }
  hourly<-1
  if (length(prec) == length(TAIRhr)) {
    rainhourly<-1
    RAINhr<-prec
    RAINhr[RAINhr<0.1]<-0
    raintest<-RAINhr
  } else if (length(prec) == length(TAIRhr)/24) {
    rainhourly<-0
    RAINhr<-rep(0,24*ndays)
    raintest<-rep(prec/24,each=24)
  } else stop("Rainfall must be daily or hourly")
  # Decide whether to run snowmodel
  snowtest<-ifelse(TAIRhr>0,0,-TAIRhr)*raintest
  snowmodel <- 0
  if (max(snowtest)>0) snowmodel<-1
  RUF<-roughlength(Veghyt,mean(PAI),0.0003)
  D0<-zeroplanedis(Veghyt, mean(PAI))
  microinput<-c(ndays,RUF,ERR,Usrhyt,Refhyt,Numtyps,0,0,0,0,1,ida,
                HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,1,
                microdaily,tannul,0.0167238,VIEWF,snowtemp,snowdens,snowmelt,undercatch,
                rainmult,0,1,maxpool,0,snowmodel,rainmelt,
                0,densfun,hourly,rainhourly,0,0,RW,PC,RL,SP,R1,IM,
                500,0,0,fail,0,intercept,grasshade,0,0,D0)
  doy1<-matrix(data=0,nrow=ndays,ncol=1)
  SLES1<-matrix(data=0,nrow=ndays,ncol=1)
  MAXSHADES1<-matrix(data=0,nrow=ndays,ncol=1)
  MINSHADES1<-matrix(data=0,nrow=ndays,ncol=1)
  TMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  TMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  CCMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  CCMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  RHMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  RHMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  WNMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  WNMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  REFLS1<-matrix(data= 0,nrow=ndays,ncol=1)
  PCTWET1<-matrix(data=0,nrow=ndays, ncol=1)
  RAINFALL1<-matrix(data=0,nrow=ndays,ncol=1)
  tannul1<-matrix(data=0,nrow=ndays,ncol=1)
  moists1<-matrix(data=0,nrow=10,ncol=ndays)
  doy1[1:ndays]<-doy
  SLES1[1:ndays]<-SLES
  MAXSHADES1[1:ndays]<-MAXSHADES
  MINSHADES1[1:ndays]<-MINSHADES
  TMAXX1[1:ndays]<-TMAXX
  TMINN1[1:ndays]<-TMINN
  CCMAXX1[1:ndays]<-CCMAXX
  CCMINN1[1:ndays]<-CCMINN
  RHMAXX1[1:ndays]<-RHMAXX
  RHMINN1[1:ndays]<-RHMINN
  WNMAXX1[1:ndays]<-WNMAXX
  WNMINN1[1:ndays]<-WNMINN
  REFLS1[1:ndays]<-REFLS
  PCTWET1[1:ndays]<-PCTWET
  RAINFALL1[1:ndays]<-RAINFALL
  tannul1[1:ndays]<-tannul
  moists1[1:10, 1:ndays] <- moists
  tides<-matrix(data=0,nrow=24*ndays,ncol=3)
  TIMAXS<-c(1,1,0,0)
  TIMINS<-c(0,0,1,1)
  LAI<-pLAI*PAI
  micro<-list(tides=tides,microinput=microinput,doy=doy,SLES=SLES1,DEP=DEP,Nodes=Nodes,
              MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX1,
              TMINN=TMINN1,RHMAXX=RHMAXX1,RHMINN=RHMINN1,CCMAXX=CCMAXX1,CCMINN=CCMINN1,
              WNMAXX=WNMAXX1,WNMINN=WNMINN1,TAIRhr=TAIRhr,RHhr=RHhr,WNhr=WNhr,CLDhr=CLDhr,
              SOLRhr=SOLRhr,RAINhr=RAINhr,ZENhr=ZENhr,IRDhr=IRDhr,REFLS=REFLS1,PCTWET=PCTWET1,
              soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists1,
              RAINFALL=RAINFALL1,tannulrun=deepsoil,PE=PE,KS=KS,BB=BB,BD=BD,DD=DD,L=L,LAI=LAI)
  microut<-microclimate(micro)
  metout<-as.data.frame(microut$metout)
  soil<-as.data.frame(microut$soil)
  soilmoist<-as.data.frame(microut$soilmoist)
  plant<-as.data.frame(microut$plant)
  if (snowmodel == 1) {
    snow <- microut$sunsnow
  } else snow <- 0
  if (max(metout[,1]==0)) stop("ERROR: the model crashed - try a different error tolerance spacing in DEP")
  return(list(metout=metout,soiltemps=soil,soilmoist=soilmoist,snowtemp=snow,plant=plant))
}

