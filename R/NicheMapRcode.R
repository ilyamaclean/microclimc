#' Runs NicheMapR model
#'
#' @description Wrapper function for running NicheMapR to derive soil moistures
#' @param climdata a data.frame of hourly weather variables in same format as [microclimc::weather()].
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
#' @param SLE Thermal emissivity of ground
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
runNMR <- function(climdata, prec, lat, long, Usrhyt, Veghyt, Refhyt = 2, PAI = 3, LOR = 1,
                   pLAI = 0.8, clump = 0, REFL = 0.15, LREFL = 0.4, SLE = 0.95,
                   DEP = c(0,2.5,5,10,15,20,30,50,100,200), ALTT = 0, SLOPE = 0, ASPECT = 0,
                   ERR = 1.5, soiltype = "Loam", PE=rep(1.1,19), KS=rep(0.0037,19), BB=rep(4.5,19),
                   BD=rep(1.3,19), DD=rep(2.65,19), cap = 1, hori = rep(0,36), maxpool = 1000,
                   rainmult = 1, SoilMoist_Init = c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3)) {
  if (Veghyt > 2) Veghyt<-2
  loc<-c(long,lat)
  tmehr<-as.POSIXlt(climdata$obs_time,tz="UTC")
  nyears<-length(unique(tmehr$year))
  fail<-nyears*24*365
  ystart<-tmehr$year[1]+1900
  yfinish<-tmehr$year[length(tmehr)]+1900
  yearlist<-seq(ystart,(ystart+(nyears-1)),1)
  doy <- unique(as.numeric(strftime(tmehr, format = "%j")))
  ndays <- unique(paste(as.numeric(strftime(tmehr, format = "%j")),
                        as.numeric(strftime(tmehr, format = "%y"))))
  ndays <- length(ndays)
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
  rad_dir<-climdata$swrad-climdata$difrad
  rad_ground<-si*dirtr*rad_dir+climdata$difrad*diftr
  MINSHADES<-(1-rad_ground/climdata$swrad)*100
  MINSHADES[MINSHADES>99.9]<-99.9
  MINSHADES[MINSHADES<0]<-0
  MINSHADES[is.na(MINSHADES)]<-mean(MINSHADES,na.rm=T)
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
  TAIRhr<-climdata$temp
  SOLRhr<-climdata$swrad*VIEWF
  sb<-5.67*10^-8
  IRDhr<-climdata$skyem*sb*(climdata$temp+273.15)^4
  RHhr<-climdata$relhum
  RHhr[RHhr>100]<-100
  RHhr[RHhr<0]<-0
  e0<-satvap(TAIRhr)
  ea<-e0*(RHhr/100)
  eo<-1.24*(10*ea/(TAIRhr+273.15))^(1/7)
  CLDhr<-((climdata$skyem-eo)/(1-eo))*100
  CLDhr[CLDhr<0]<-0
  CLDhr[CLDhr>100]<-100
  WNhr<-climdata$windspeed
  WNhr[is.na(WNhr)]<-0.1
  PRESShr<-climdata$pres*1000
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
    RAINhr<-rep(prec/24,each=24)
  } else stop("Rainfall must be daily or hourly")
  # Decide whether to run snowmodel
  snowtest<-ifelse(TAIRhr>0,0,-TAIRhr)*RAINhr
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
  raind<-matrix(RAINhr,ncol=24,byrow=TRUE)
  raind<-apply(raind,1,sum)
  RAINFALL1[1:ndays]<-raind
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
    snow <- as.data.frame(microut$sunsnow)
  } else snow <- 0
  return(list(metout=metout,soiltemps=soil,soilmoist=soilmoist,snowtemp=snow,plant=plant))
}
#' Internal function for calculating lead absorbed radiation on vector
.leafabs2 <-function(Rsw, tme, tair, tground, lat, long, PAIt, PAIu, pLAI, x, refls, refw, refg, vegem, skyem, dp,
                     merid = round(long/15, 0) * 15, dst = 0, clump = 0) {
  jd<-jday(tme = tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-solalt(lt,lat,long,jd,merid=merid,dst=dst)
  sa2<-ifelse(sa<5,5,sa)
  ref<-pLAI*refls+(1-pLAI)*refw
  sunl<-psunlit(PAIu,x,sa,clump)
  mul<-radmult(x,sa2)
  aRsw <- (1-ref) * cansw(Rsw,dp,jd,lt,lat,long,PAIu,x,ref,merid=merid,dst=dst,clump=clump)
  aGround <- refg * cansw(Rsw,dp,jd,lt,lat,long,PAIt,x,ref,merid=merid,dst=dst,clump=clump)
  aGround <- (1-ref) * cansw(aGround,dp,jd,lt,lat,long,PAIt,x,ref,merid=merid,dst=dst,clump=clump)
  aRsw <- sunl * mul * aRsw + (1-sunl) * aRsw + aGround
  aRlw <- canlw(tair, PAIu, 1-vegem, skyem = skyem, clump = clump)$lwabs
  return(aRsw+aRlw)
}
#' Internal function for calculating lead absorbed radiation on vector
.windprofile <- function(ui, zi, zo, a = 2, hgt, PAI, psi_m = 0, hgtg = 0.05 * hgt, zm0 = 0.004) {
  d<-zeroplanedis(hgt,PAI)
  zm<-roughlength(hgt,PAI,zm0)
  ln2<-suppressWarnings(log((zi-d)/zm)+psi_m)
  ln2[ln2<0.55]<-0.55
  uf<-(0.4*ui)/ln2
  if (zo >= hgt) {
    ln2<-suppressWarnings(log((zo-d)/zm)+psi_m)
    ln2[ln2<0.55]<-0.55
    uo<-(uf/0.4)*ln2
  } else {
    ln2<-suppressWarnings(log((hgt-d)/zm)+psi_m)
    ln2[ln2<0.55]<-0.55
    uh<-(uf/0.4)*ln2
    if (zi >= (0.1*hgt)) {
      uo<-uh*exp(a*(zo/hgt)-1)
    } else {
      uo<-uh*exp(a*((0.1*hgt)/hgt)-1)
      zmg<-0.1*hgtg
      ln1<-log((0.1*hgt)/zmg)+psi_m
      uf<-0.4*(uo/ln1)
      ln2<-log(zo/zmg)+psi_m
      uo<-(uf/0.4)*ln2
    }
  }
  uo
}
#' Steady-state leaf and air temperature
#' @description Rapid method for calculating steady-state leaf and air temperature at one height
#' @param tair air temperature at reference height (deg C)
#' @param tground ground surface temperature (as returned by e.g. [runNMR()])
#' @param relhum relative humidity at reference height (percentage)
#' @param pk atmospheric pressure (kPa)
#' @param theta volumetric water content of upper most soil layer in current time step (m^3 / m^3)
#' as returned by e.g. [runNMR()]
#' @param gtt molar conductivity from leaf height to reference height as returned by [gturb()]
#' and [gcanopy()]
#' @param gt0 molar conductivity from leaf height to ground as returned by [gcanopy()]
#' @param gha boundary layer conductance of leaf as returned by [gforcedfree()]
#' @param gL combined boundary layer and leaf-air conductance given by 1/(1/gha+1/(uz*ph))
#' @param gv combined boundary layer and stomatal conductance of leaf
#' @param Rabs radiation absorbed by leaf as returned by e.g. [leafabs()]
#' @param soilb Shape parameter for Campbell soil model (dimensionless, > 1) as returned by
#' [soilinit()]
#' @param Psie Soil matric potential (J / m^3) as returned by [soilinit()]
#' @param Smax Volumetric water content at saturation [m3 / m3]
#' @param surfwet proportion of leaf surface acting as free water surface
#' @param leafdens Total one sided leaf area per m^3 at desired height
#' @return a list of the following:
#' @return `tleaf` leaf temperature
#' @return `tn` air temperature
#' @return `rh` relative humidity
#' @import microctools
#' @export
tleafS <- function(tair, tground, relhum, pk, theta, gtt, gt0, gha, gv, gL, Rabs, vegem, soilb,
                   Psie, Smax, surfwet, leafdens) {
  ###
  cp<-cpair(tair)
  # Air temperature expressed as leaf temperature
  aL<-(gtt*(tair+273.15)+gt0*(tground+273.15))/(gtt+gt0)
  bL<-(leafdens*gL)/(gtt+gt0)
  # Vapour pressures
  es<-satvap(tair)
  eref <- (relhum/100)*es
  rhsoil<-soilrh(theta,soilb,Psie,Smax,tground)
  esoil<-rhsoil*satvap(tground)
  # add a small correction to improve delta estimate
  sb<-5.67*10^-8
  Rnet<-Rabs-0.97*sb*(tair+273.15)^4
  tle<-tair+(0.5*Rnet)/(cp*leafdens*gha)
  tapprox<-(tle+tair)/2
  delta <- 4098*(0.6108*exp(17.27*tapprox/(tapprox+237.3)))/(tapprox+237.3)^2
  ae<-(gtt*eref+gt0*esoil+gv*es)/(gtt+gt0+gv)
  be<-(gv*delta)/(gtt+gt0+gv)
  # Sensible heat
  bH<-gha*cp
  # Latent heat
  lambda <- (-42.575*tair+44994)
  aX<-((lambda*gv)/pk)*(surfwet*es-ae)
  bX<-((lambda*gv)/pk)*(surfwet*delta-be)
  aX[aX<0]<-0
  bX[bX<0]<-0
  # Emmited radiation
  aR<-sb*vegem*aL^4
  bR<-4*vegem*sb*(aL^3*bL+(tair+273.15)^3)
  # Leaf temperature
  dTL <- (Rabs-aR-aX)/(1+bR+bX+bH)
  # tz pass 1
  tn<-aL-273.15+bL*dTL
  tleaf<-tn+dTL
  # new vapour pressure
  eanew<-ae+be*dTL
  eanew[eanew<0.01]<-0.01
  tmin<-dewpoint(eanew,tn,ice = TRUE)
  esnew<-satvap(tn)
  eanew<-ifelse(eanew>esnew,esnew,eanew)
  rh<-(eanew/esnew)*100
  # Set both tair and tleaf so as not to drop below dewpoint
  tleaf<-ifelse(tleaf<tmin,tmin,tleaf)
  tn<-ifelse(tn<tmin,tmin,tn)
  tmax<-ifelse(tn+20<80,tn+20,80)
  tleaf<-ifelse(tleaf>tmax,tmax,tleaf)
  # cap upper limits of both tair and tleaf
  tmx<-pmax(tground+5,tair+5)
  tn<-ifelse(tn>tmx,tmx,tn)
  tmx<-pmax(tground+30,tair+30)
  tleaf<-ifelse(tleaf>tmx,tmx,tleaf)
  # cap lower limits
  tmn<-tair-7
  tn<-ifelse(tn<tmn,tmn,tn)
  tmn<-tair-20
  tleaf<-ifelse(tleaf<tmn,tmn,tleaf)
  return(list(tleaf=tleaf,tn=tn,rh=rh))
}
#' Internal function for running model with snow
.runmodelsnow <- function(climdata, vegp, soilp, nmrout, reqhgt, lat, long, metopen = TRUE, windhgt = 2) {
  # Snow
  snow<-nmrout$snow
  if (class(snow) == "data.frame") {
    snowtemp<-snow$SN1
  } else snowtemp<-rep(0,length(L))
  metout<-nmrout$metout
  snowdep<-metout$SNOWDEP
  # (1) Unpack variables
  tair<-climdata$temp
  relhum<-climdata$relhum
  pk<-climdata$pres
  u<-climdata$windspeed
  hgt<-vegp$hgt
  # Esimate PAIu
  PAIt<-apply(vegp$PAI,2,sum)
  PAIt[PAIt<0.001]<-0.001
  LAI<-(vegp$pLAI*vegp$PAI)
  LAI<-apply(LAI,2,mean)
  pLAI<-LAI/PAIt
  pLAI[is.na(pLAI)]<-0.8
  m<-length(vegp$iw)
  z<-c((1:m)-0.5)/m*hgt
  if (reqhgt < hgt) {
    sel<-which(z>reqhgt)
    wgt1<-abs(z[sel[1]]-reqhgt)
    wgt2<-abs(z[sel[1]-1]-reqhgt)
    sel2<-c(sel[1]-1,sel)
    PAIu1<-vegp$PAI[sel,]
    PAIu1<-apply(PAIu1,2,sum)
    PAIu2<-vegp$PAI[sel2,]
    PAIu2<-apply(PAIu2,2,sum)
    PAIu<-PAIu1+(wgt1/(wgt1+wgt2))*(PAIu2-PAIu1)
    dif<-abs(z-reqhgt)
    sel<-which(dif==min(dif))
    leafdens<-vegp$PAI[sel,]/(z[2]-z[1])
  } else PAIu<-rep(0,length(PAIt))
  # Calculate wind speed 2 m above canopy
  if (metopen) {
    if (windhgt != 2) u <- u*4.87/log(67.8*windhgt-5.42)
    u2<-u*log(67.8*(hgt+2)-5.42)/log(67.8*2-5.42)
  } else {
    u2 <- u
    if (windhgt != (2+hgt)) u2 <- windprofile(u, windhgt, hgt+2, a = 2, PAIt, hgt)
  }
  u2[u2<0.5]<-0.5
  cp<-cpair(tair)
  ph<-phair(tair,pk)
  sb<-5.67*10^-8
  Rsw<-0.06*climdata$swrad
  Rlw<-(1-climdata$skyem)*sb*0.85*(snowtemp+273.15)^4
  Rnet<-Rsw-Rlw
  H<-0.65*Rnet
  # Calculate diabatic correction factor
  zm<-ifelse(hgt>snowdep,roughlength(hgt, PAIt),0.003)
  d<-ifelse(hgt>snowdep,zeroplanedis(hgt, PAIt),snowdep-0.003)
  hgt2<-ifelse(hgt>snowdep,hgt,snowdep)
  uf <- (0.4*u2)/log((2+hgt2-d)/zm)
  dba <- diabatic_cor(tair,pk,H,uf,hgt+2,d)
  dba$psi_m<-ifelse(dba$psi_m>2.5,2.5,dba$psi_m)
  dba$psi_m<-ifelse(dba$psi_h>2.5,2.5,dba$psi_h)
  # Wind speed at user height
  a1<-attencoef(hgt,PAIt,vegp$x,vegp$lw,vegp$cd,mean(vegp$iw))
  a2<-attencoef(hgt,0,vegp$x,vegp$lw,vegp$cd,mean(vegp$iw))
  uz1<-.windprofile(u2,hgt+2,reqhgt,a1,hgt,PAIt,dba$psi_m)
  uz2<-.windprofile(u2,snowdep+2,reqhgt,a2,0,0,dba$psi_m)
  uz<-ifelse(hgt>snowdep,uz1,uz2)
  uf<-(0.4*u2)/(log((hgt2+2-d)/zm)+dba$psi_m)
  uf[uf<0.065]<-0.065
  sas <- which(reqhgt>=snowdep)
  sbs <- which(reqhgt<snowdep)
  # Below snow
  if (length(sbs)>0) {
    lyr<-round((reqhgt/snowdep[sbs])*9,0)
    lyr[lyr<1]<-1
    lyr[lyr>9]<-9
    tz2<-0
    for(i in 1:length(lyr)) tz2[i]<-snow[sbs[i],lyr[i]+2]
    rh2<-rep(100,length(tz2))
    Rsw2<-rep(0,length(tz2))
    Rlw2<-5.67*10^-8*0.85*(snowtemp+273.15)^4
    Rlw2<-Rlw2[sbs]
    if (reqhgt > hgt) {
      tleaf2<-rep(-999,length(sbs))
    } else tleaf2<-snowtemp[sbs]
    ws2<-rep(0,length(sbs))
  }
  if (reqhgt >= hgt) {
    # Above snow
    if (length(sas)>0) {
      xx<-(H[sas]/(0.4*ph[sas]*cp[sas]*uf[sas]))
      T0<-snowtemp[sas]
      psihe<-(T0-tair[sas])/xx-log((hgt+2-d[sas])/(0.2*zm[sas]))
      rat<-log((reqhgt-d[sas])/(0.2*zm[sas]))/log((hgt+2-d[sas])/(0.2*zm[sas]))
      psihe<-rat*psihe
      tz1<-T0-xx*(log((reqhgt-d[sas])/(0.2*zm[sas]))+psihe)
      tmx<-pmax(tair[sas],T0)
      tmn<-pmin(tair[sas],T0)
      tz1<-ifelse(tz1>tmx,tmx,tz1)
      tz1<-ifelse(tz1<tmn,tmn,tz1)
      ea<-satvap(tair[sas])*(relhum[sas]/100)
      rh1<-(ea/satvap(tz1))*100
      Rsw<-climdata$swrad
      Rlw<-5.67*10^-8*climdata$skyem*(tair+273.15)^4
      Rsw<-Rsw[sas]
      Rlw<-Rlw[sas]
      tleaf1<-rep(-999,length(sas))
    }
    tz<-rep(0,length(uf))
    rh<-rep(0,length(uf))
    ws<-rep(0,length(uf))
    tz[sas]<-tz1
    tz[sbs]<-tz2
    rh[sas]<-rh1
    rh[sbs]<-rh2
    tleaf<-rep(-999,length(tz))
    Rswa<-rep(0,length(uf))
    Rlwa<-rep(0,length(uf))
    Rswa[sas]<-Rsw
    Rswa[sbs]<-Rsw2
    Rlwa[sas]<-Rlw
    Rlwa[sbs]<-Rlw2
    ws[sas]<-uz[sas]
    ws[sbs]<-ws2
  } else {
    # Above snow
    if (length(sas)>0) {
      ln2 <- suppressWarnings(log((hgt-d[sas])/zm[sas])+dba$psi_m[sas])
      ln2[ln2<0.55]<-0.55
      uh <- (uf[sas]/0.4)*ln2
      # Conductivities
      gta <- gturb(u2[sas],hgt+2,hgt+2,hgt,hgt,PAIt[sas],tair[sas],dba$psi_m[sas],dba$psi_h[sas],0.004,pk[sas])
      gtc <- gcanopy(uh,hgt,0.004,tair[sas],tair[sas],hgt,PAIt[sas],vegp$x,vegp$lw*2,vegp$cd,mean(vegp$iw),1,pk[sas])
      gt0 <- gcanopy(uh,0.004,0,tair[sas],tair[sas],hgt,PAIt[sas],vegp$x,vegp$lw,vegp$cd,mean(vegp$iw),1,pk[sas])
      gtt <- 1/(1/gta+1/gtc)
      gha <- 1.41*gforcedfree(vegp$lw*0.71,uz[sas],tair[sas],5,pk[sas],5)
      # Radiation
      jd<-jday(tme=tme)
      lt <- tme$hour+tme$min/60+tme$sec/3600
      lt<-lt
      dp<-climdata$difrad/climdata$swrad
      dp[is.na(dp)]<-0.5
      dp[is.infinite(dp)]<-0.5
      dp[dp>1]<-1
      Rsw <- cansw(climdata$swrad,dp,tme=tme,lat=lat,long=long,x=vegp$x,l=PAIu,ref=vegp$refls)
      Rsw<-Rsw[sas]
      Rlw <- canlw(tair, PAIu, 1-vegp$vegem, climdata$skyem, vegp$clump)$lwin
      Rlw <-Rlw[sas]
      gv <- layercond(Rsw, vegp$gsmax, vegp$q50)
      gv <-1/(1/gv+1/gha)
      # Leaf absorbed radiation
      Rabs<-.leafabs2(climdata$swrad,tme,tair,tground,lat,long,PAIt,PAIu,pLAI,vegp$x,0.95,0.95,
                      0.95,0.85,climdata$skyem,dp,vegp$clump)
      Rabs<-Rabs[sas]
      soilm<-nmrout$soilmoist
      theta<-soilm$WC0cm[sas]
      tln<-tleafS(tair[sas],snowtemp[sas],relhum[sas],pk[sas],theta,gtt,gt0,gha,gv,Rabs,
                  0.5,soilp$b,soilp$psi_e,soilp$Smax,1,leafdens[sas])
      tleaf1<-pmax(tln$tleaf,snowtemp[sas])
      tz1<-tln$tn
      rh1<-tln$rh
    }
    tz<-rep(0,length(uf))
    rh<-rep(0,length(uf))
    tleaf<-rep(0,length(uf))
    ws<-rep(0,length(uf))
    tz[sas]<-tz1
    tz[sbs]<-tz2
    tleaf[sas]<-tleaf1
    tleaf[sbs]<-tleaf2
    rh[sas]<-rh1
    rh[sbs]<-rh2
    Rswa<-rep(0,length(uf))
    Rlwa<-rep(0,length(uf))
    Rswa[sas]<-Rsw
    Rswa[sbs]<-Rsw2
    Rlwa[sas]<-Rlw
    Rlwa[sbs]<-Rlw2
    ws[sas]<-uz[sas]
    ws[sbs]<-ws2
  }
  metout<-data.frame(obs_time=climdata$obs_time,Tref=climdata$temp,Tloc=tz,tleaf=tleaf,
                     RHref=relhum,RHloc=rh,RSWloc=Rswa,RLWloc=Rlwa,windspeed=ws)
  return(metout)
}
#' Run model under steady state conditions
#' @description Rapid method for calculating below or above canopy or below ground microclimate
#' under steady-state at one user specified height.
#' @param climdata  data.frame of climate variables needed to run the run the model (dataset should follow format of [weather()])
#' @param vegp a list of vegetation parameters as returned by [microctools::habitatvars()].
#' @param soilp a list of soil parameters as returned by [soilinit()]
#' @param nmrout a list of putputs from NicheMapR as returned by [runNMR()].
#' @param reqhgt height (m) for which microclimate is needed.
#' @param lat latitude of location (decimal degrees).
#' @param long longitude of location (decimal degrees).
#' @param metopen optional logical indicating whether the wind measurement used as an input to
#' the model is from a nearby weather station located in open ground (TRUE) or above the canopy
#' for which temperatures are modelled (FALSE - see details)
#' @param windhgt height above ground of wind measurement. If `metopen` is FALSE, must be above
#' canopy.
#' @param surfwet proportion of leaf surface acting as free water surface
#' @param groundem thermal emissivity of ground layer
#' @return a data.frame with the following columns:
#' @return `obs_time` time of observation. Same as in climdata.
#' @return `Tref` Temperature (deg C) at reference height as in climdata
#' @return `Tloc` Temperature (deg C) at height `reqhgt`
#' @return `tleaf` Leaf temperature (deg C) at height `reqhgt`. -999 if `reqhgt` above canopy
#' or below ground
#' @return `RHref` Relative humidity (percentage) at reference height as in climdata.
#' @return `RHloc` Relative humidity (percentage) at height `reqhgt`
#' @return `RSWloc` Total incoming shortwave radiation at height `reghgt` (W/m^2)
#' @return `RLWloc` Total downward longwave radiation at height `reqhgt` (W/m^2)
#' @import microctools
#' @export
#'
#' @details This is a rapid implementation of model when time increments are hourly such that
#' transient heat fluxes and heat storage in canopy can be ignored. Computations are performed
#' simultaniously on all data, negating need to run in timesteps. Also includes implementation of
#' model with snow present. Should generally be run using wrapper function [runwithNMR()] but
#' provided as a standalone function in case data for multiple heights are needed, in whihc case
#' it can be run multiple times without also running NicheMapR.
runmodelS <- function(climdata, vegp, soilp, nmrout, reqhgt,  lat, long, metopen = TRUE, windhgt = 2,
                      surfwet = 1, groundem = 0.95) {
  # (1) Unpack variables
  tme<-as.POSIXlt(climdata$obs_time)
  tair<-climdata$temp
  relhum<-climdata$relhum
  pk<-climdata$pres
  u<-climdata$windspeed
  hgt<-vegp$hgt
  # Ensure at least some PAI
  vegp$PAI[vegp$PAI<0.0001]<-0.0001
  # Esimate total PAI and proportion LAI
  PAIt<-apply(vegp$PAI,2,sum)
  LAI<-(vegp$pLAI*vegp$PAI)
  LAI<-apply(LAI,2,mean)
  pLAI<-LAI/PAIt
  m<-length(vegp$iw)
  z<-c((1:m)-0.5)/m*hgt
  # Esimate PAI above point and PAI of layer
  if (reqhgt < hgt) {
    sel<-which(z>reqhgt)
    if (length(sel) > 1) {
      wgt1<-abs(z[sel[1]]-reqhgt)
      wgt2<-abs(z[sel[1]-1]-reqhgt)
      sel2<-c(sel[1]-1,sel)
      PAIu1<-vegp$PAI[sel,]
      if (length(sel) > 1) PAIu1<-apply(PAIu1,2,sum)
      PAIu2<-vegp$PAI[sel2,]
      if (length(sel2) > 1) PAIu2<-apply(PAIu2,2,sum)
      if (length(wgt2) > 1) {
        PAIu<-PAIu1+(wgt1/(wgt1+wgt2))*(PAIu2-PAIu1)
      } else PAIu<-PAIu1
    } else {
      zu<-z[length(z)]
      wgt<-(hgt-reqhgt)/(hgt-zu)
      PAIu<-wgt*vegp$PAI[length(z),]
    }
    dif<-abs(z-reqhgt)
    sel<-which(dif==min(dif))[1]
    leafdens<-vegp$PAI[sel,]/(z[2]-z[1])
  } else PAIu<-rep(0,length(PAIt))
  # (2) Estimate Sensible Heat flux
  # (2a) Latent heat flux
  plant<-nmrout$plant
  lambda <- (-42.575*tair+44994)
  L<-(lambda*plant$TRANS)/(3600*18.01528)
  metout<-nmrout$metout
  snowdep<-metout$SNOWDEP
  selsnow<-which(snowdep > 0)
  # (2b) Ground heat flux
  # Calculate wind speed 2 m above canopy
  if (metopen) {
    if (windhgt != 2) u <- u*4.87/log(67.8*windhgt-5.42)
    u2<-u*log(67.8*(hgt+2)-5.42)/log(67.8*2-5.42)
  } else {
    u2 <- u
    if (windhgt != (2+hgt)) u2 <- windprofile(u, windhgt, hgt+2, a = 2, PAIt, hgt)
  }
  cp<-cpair(tair)
  ph<-phair(tair,pk)
  u2[u2<0.5]<-0.5
  zm<-roughlength(hgt, PAIt)
  d<-zeroplanedis(hgt, PAIt)
  uf <- (0.4*u2)/log((2+hgt-d)/zm)
  # Wind speed at user height
  a<-attencoef(hgt,PAIt,vegp$x,vegp$lw,vegp$cd,mean(vegp$iw))
  uz<-.windprofile(u2,hgt+2,reqhgt,a,hgt,PAIt)
  ln2 <- suppressWarnings(log((hgt - d) / zm))
  ln2[ln2<0.55]<-0.55
  uh <- (uf/0.4)*ln2
  # Conductivities
  hgt2<-ifelse(hgt>2,2,hgt)
  gta <- gturb(u,2,2,NA,hgt2,PAIt,tair,0,0,0.004,pk)
  # Extract soil temperature
  soilt<-nmrout$soiltemps
  tground<- soilt$D0cm
  mult<-cp*gta
  # Calculate max conductivity
  maxflux<-(cp*ph*(hgt+2))/3600
  mult<-ifelse(mult>maxflux,maxflux,mult)
  G<-mult*(tair-tground)
  sb<-sb<-5.67*10^-8
  Rsw<-(1-vegp$refg)*climdata$swrad
  Rlw<-(1-climdata$skyem)*sb*groundem*(tground+273.15)^4
  Rnet<-Rsw-Rlw
  H<-Rnet-G-L
  # Recalculate with diabatic lapse rate
  dba <- diabatic_cor(tair,pk,H,uf,hgt+2,d)
  dba$psi_m<-ifelse(dba$psi_m>2.5,2.5,dba$psi_m)
  dba$psi_m<-ifelse(dba$psi_h>2.5,2.5,dba$psi_h)
  # Wind speed at user height
  a<-attencoef(hgt,PAIt,vegp$x,vegp$lw,vegp$cd,mean(vegp$iw))
  uz<-.windprofile(u2,hgt+2,vegp$hgtg,a,hgt,PAIt,dba$psi_m)
  uf<-(0.4*u2)/(log((hgt+2-d)/zm)+dba$psi_m)
  uf[uf<0.065]<-0.065
  xx<-(H/(0.4*ph*cp*uf))
  T0<-tair+xx*(log((hgt+2-d)/(0.2*zm))+dba$psi_h)
  xx<-(H/(0.4*ph*cp*uf))
  T0<-tair+xx*(log((hgt+2-d)/(0.2*zm))+dba$psi_h)
  ea<-satvap(tair)*(relhum/100)
  tmn<-pmax(dewpoint(ea,tair),tair-5)
  T0<-ifelse(T0<tmn,tmn,T0)
  tmx<-pmax(tair+20,tground)
  T0<-ifelse(T0>tmx,tmx,T0)
  T0<-ifelse(T0>80,80,T0)
  if (reqhgt >= hgt) {
    psihe<-(T0-tair)/xx-log((hgt+2-d)/(0.2*zm))
    rat<-log((reqhgt-d)/(0.2*zm))/log((hgt+2-d)/(0.2*zm))
    psihe<-rat*psihe
    tz<-T0-xx*(log((reqhgt-d)/(0.2*zm))+psihe)
    tmx<-pmax(tair,T0)
    tmn<-pmin(tair,T0)
    tz<-ifelse(tz>tmx,tmx,tz)
    tz<-ifelse(tz<tmn,tmn,tz)
    rh<-(ea/satvap(tz))*100
    tleaf<-rep(-999,length(tz))
    Rsw<-climdata$swrad
    Rlw<-5.67*10^-8*climdata$skyem*(tair+273.15)^4
  } else {
    ln2 <- suppressWarnings(log((hgt - d) / zm) + dba$psi_m)
    ln2[ln2<0.55]<-0.55
    uh <- (uf/0.4)*ln2
    # Conductivities
    gta <- gturb(u2,hgt+2,hgt+2,hgt,hgt,PAIt,tair,dba$psi_m,dba$psi_h,0.004,pk)
    gtc <- gcanopy(uh,hgt,reqhgt,tair,tair,hgt,PAIt,vegp$x,vegp$lw,vegp$cd,mean(vegp$iw),1,pk)
    gtc[gtc<1]<-1
    gt0 <- gcanopy(uh,reqhgt,0,tair,tair,hgt,PAIt,vegp$x,vegp$lw,vegp$cd,mean(vegp$iw),1,pk)
    gt0[gt0<1]<-1
    gtt <- 1/(1/gta+1/gtc)
    gha <- 1.41*gforcedfree(vegp$lw*0.71,uz,tair,5,pk,5)
    # Radiation
    jd<-jday(tme=tme)
    lt <- tme$hour+tme$min/60+tme$sec/3600
    dp<-climdata$difrad/climdata$swrad
    dp[is.na(dp)]<-0.5
    dp[is.infinite(dp)]<-0.5
    dp[dp>1]<-1
    Rsw <- cansw(climdata$swrad,dp,tme=tme,lat=lat,long=long,x=vegp$x,l=PAIu,ref=vegp$refls)
    Rlw <- canlw(tair, PAIu, 1-vegp$vegem, climdata$skyem, vegp$clump)$lwin
    gv <- layercond(Rsw, vegp$gsmax, vegp$q50)
    gv <-1/(1/gv+1/gha)
    # Leaf absorbed radiation
    Rabs<-.leafabs2(climdata$swrad,tme,tair,tground,lat,long,PAIt,PAIu,pLAI,vegp$x,vegp$refls,vegp$refw,
                    vegp$refg,vegp$vegem,climdata$skyem,dp,vegp$clump)
    soilm<-nmrout$soilmoist
    theta<-soilm$WC0cm
    tln<-tleafS(tair,tground,relhum,pk,theta,gtt,gt0,gha,gv,Rabs,vegp$vegem,soilp$b,soilp$psi_e,
                soilp$Smax,surfwet,leafdens)
    tleaf<-tln$tleaf
    tz<-tln$tn
    tleaf<-pmin(tleaf,T0)
    tz<-pmin(tz,T0)
    rh<-tln$rh
  }
  metout<-data.frame(obs_time=climdata$obs_time,Tref=climdata$temp,Tloc=tz,tleaf=tleaf,
                     RHref=relhum,RHloc=rh,RSWloc=Rsw,RLWloc=Rlw,windspeed=uz)
  # Cap at theoretical upper and lower limits
  mx<-pmax(tground,T0)
  mn<-pmin(tground,tair-7)
  metout$Tloc[metout$Tloc>mx]<-mx[metout$Tloc>mx]
  metout$tleaf[metout$tleaf>mx]<-mx[metout$tleaf>mx]
  metout$Tloc[metout$Tloc<mn]<-mn[metout$Tloc<mn]
  metout$tleaf[metout$tleaf<mn]<-mn[metout$tleaf<mn]
  # Consider snow
  snow<-nmrout$snow
  if (class(snow) == "data.frame") {
    snowtemp<-nmrout$SN1
  } else snowtemp<-rep(0,length(L))
  mo<-nmrout$metout
  snowdep<-mo$SNOWDEP
  if (max(snowdep) > 0) {
    mos <- .runmodelsnow(climdata,vegp,soilp,nmrout,reqhgt,lat,long,metopen,windhgt)
    sel<-which(snowdep>0)
    metout$Tloc[sel]<-mos$Tloc[sel]
    metout$tleaf[sel]<-mos$tleaf[sel]
    metout$RHloc[sel]<-mos$RHloc[sel]
    metout$RSWloc[sel]<-mos$RSWloc[sel]
    metout$RLWloc[sel]<-mos$RLWloc[sel]
    metout$windspeed[sel]<-mos$windspeed[sel]
  }
  return(metout)
}
#' Runs microclimate model in hourly timesteps with NicheMapR
#'
#' @description This function is a wrapper function for [runNMR()] and [runmodelS()] enabling
#' full microclimate model at any height above or below ground rapidly. NicheMapR is to
#' compute snowdepths, soil moistures and below-ground soil temperatures. The function
#' runmodelS is used to compute below or above canopy air temperatures.
#'
#' @param climdata data.frame of climate variables needed to run the run the model. The dataset should follow format
#' of [weather()]). Times must be in UTC.
#' @param prec vector of hourly or daily precipitation (mm) (see details).
#' @param vegp a list of vegetation parameters as returned by [microctools::habitatvars()].
#' @param soilp a list of soil parameters as returned by [soilinit()].
#' @param reqhgt height (m) for which microclimate is needed.
#' @param lat latitude of location (decimal degrees).
#' @param long longitude of location (decimal degrees).
#' @param altt elevation of location (m).
#' @param slope slope of ground surface (decimal degrees).
#' @param aspect aspect of ground surface (decimal degrees, 0 = north).
#' @param metopen optional logical indicating whether the wind measurement used as an input to
#' the model is from a nearby weather station located in open ground (TRUE) or above the canopy
#' for which temperatures are modelled (FALSE - see details)
#' @param windhgt height above ground of wind measurement. If `metopen` is FALSE, must be above
#' canopy.
#' @param surfwet proportion of leaf surface acting as free water surface.
#' @param groundem thermal emissivity of ground layer.
#' @param ERR Integrator error tolerance for soil temperature calculations.
#' @param cap Is organic cap present on soil surface? (1 = Yes, 0 = No).
#' @param hori Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 10 degree intervals.
#' @param maxpool Max depth for water pooling on the surface (mm), to account for runoff.
#' @param rainmult Rain multiplier for surface soil moisture (used to induce runon).
#' @param SoilMoist_Init Initial volumetric soil water content at each soil node (m^3/m^3)
#' @return a list of two objects:
#' @return a data.frame of microclimate variables for height `reqhgt` asreturned by [runmodelS()]
#' @return a list of NicheMapR outputs as returned by [runNMR()]
#'
#' @import microctools
#' @import NicheMapR
#' @export
#' @seealso [runmodelS()], [runNMR()]
#'
#' @details This is a wrapper function for [runmodelS()] and [runNMR()] that allows fairly rapid
#' calaulation of microclimatic conditions below ground and/or above canopy. It uses the NicheMapR package
#' to calculate soil moisture-dependent soil temperatures, treating the vegetation canopy as a
#' single layer and then calculates air and leaf temperatures,air humidity and radiation using
#' runmodelS. Precipitation values are needed for computation of soil moisture and can either be
#' supplied as a vector of hourly values (must have identical number of values as in `climdata`) or
#' daily values (`climdata` must contain climate for entire days). It compares the length of the vector
#' to `climdata` to determine whether values are daily or hourly. In this function constant soil
#' poperties are are assumed throughout the soil profile. For more flexible implementation, [runNMR()]
#' and [runmodelS()] can be run seperately.
#'
#' @examples
#' require(NicheMapR)
#' # ====================================================================== #
#' # ~~~~~~~~~~~~~~~ Run model using inbuilt weather datasets ~~~~~~~~~~~~~ #
#' # ====================================================================== #
#' tme<-as.POSIXlt(weather$obs_time,tz="UTC")
#' vegp <- habitatvars(4, 50.2178, -5.32656, tme, m = 10) # Decidious broadleaf forest
#' soilp<- soilinit("Loam")
#' climrun <- runwithNMR(weather,dailyprecip,vegp,soilp,1,50.2178,-5.32656)
#' # ====================================================================== #
#' # ~~~~~~~~~~~~~~~ Extract and view data from outputs ~~~~~~~~~~~~~~~~~~~ #
#' # ====================================================================== #
#' metout<-climrun$metout
#' nmrout<-climrun$nmrout
#' soiltemps<-nmrout$soiltemps
#' head(metout)
#' head(soiltemps)
#' head(nmrout$soilmoist)
#' head(nmrout$snowtemp)
#' head(nmrout$plant)
#' # ====================================================================== #
#' # ~~~~~~~~~~~~~~~ Plot soil and air temperatures ~~~~~~~~~~~~~~~~~~~~~~~ #
#' # ====================================================================== #
#' tsoil1<-soiltemps$D2.5cm
#' tsoil2<-soiltemps$D50cm
#' tair<-metout$Tloc
#' tref<-metout$Tref
#' tleaf<-metout$tleaf
#' Month<-as.POSIXct(metout$obs_time,tz="UTC")
#' tmn<-min(tsoil1,tsoil2,tair,tref,tleaf)
#' tmx<-max(tsoil1,tsoil2,tair,tref,tleaf)
#' par(mfrow=c(2,1))
#' # Air
#' plot(tref~Month,type="l",col=rgb(0,0,0,0.3),ylab="Temperature",ylim=c(tmn,tmx))
#' par(new=T)
#' plot(tleaf~Month,type="l",col=rgb(0,1,0,0.3),xlab="",ylab="",ylim=c(tmn,tmx))
#' par(new=T)
#' plot(tair~Month,type="l",col=rgb(1,0,0,0.3),xlab="",ylab="",ylim=c(tmn,tmx))
#' # Soil
#' plot(tref~Month,type="l",col=rgb(0,0,0,0.5),ylab="Temperature",ylim=c(tmn,tmx))
#' par(new=T)
#' plot(tsoil1~Month,type="l",col="red",xlab="",ylab="",ylim=c(tmn,tmx))
#' par(new=T)
#' plot(tsoil2~Month,type="l",col="blue",lwd=2,xlab="",ylab="",ylim=c(tmn,tmx))
#' # ====================================================================== #
#' # ~~~~~~~~~~~~~~~ Plot data by monthly mean ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#' # ====================================================================== #
#' tme<-as.POSIXlt(metout$obs_time,tz="UTC")
#' tsoil1m<-aggregate(tsoil1,by=list(tme$mon),mean)$x
#' tsoil2m<-aggregate(tsoil2,by=list(tme$mon),mean)$x
#' tairm<-aggregate(tair,by=list(tme$mon),mean)$x
#' trefm<-aggregate(tref,by=list(tme$mon),mean)$x
#' tleafm<-aggregate(tleaf,by=list(tme$mon),mean)$x
#' tmn<-min(tsoil1m,tsoil2m,tairm,trefm,tleafm)
#' tmx<-max(tsoil1m,tsoil2m,tairm,trefm,tleafm)
#' Month<-c(1:12)
#' # Air
#' par(mfrow=c(2,1))
#' plot(trefm~Month,type="l",lwd=2,col="darkgray",ylab="Temperature",ylim=c(tmn,tmx))
#' par(new=T)
#' plot(tleafm~Month,type="l",lwd=2,col="darkgreen",xlab="",ylab="",ylim=c(tmn,tmx))
#' par(new=T)
#' plot(tairm~Month,type="l",lwd=2,col="red",xlab="",ylab="",ylim=c(tmn,tmx))
# Soil
#' plot(trefm~Month,type="l",lwd=2,col="darkgray",ylab="Temperature",ylim=c(tmn,tmx))
#' par(new=T)
#' plot(tsoil1m~Month,type="l",lwd=2,col="red",xlab="",ylab="",ylim=c(tmn,tmx))
#' par(new=T)
#' plot(tsoil2m~Month,type="l",lwd=2,col="blue",xlab="",ylab="",ylim=c(tmn,tmx))
#' # ====================================================================== #
runwithNMR <- function(climdata, prec, vegp, soilp, reqhgt, lat, long, altt = 0, slope = 0,
                       aspect = 0,  metopen = TRUE, windhgt = 2, surfwet = 1, groundem = 0.95,
                       ERR = 1.5, cap = 1, hori = rep(0,36), maxpool =1000, rainmult = 1,
                       SoilMoist_Init = c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3)) {
  # (1) Unpack variables
  tair<-climdata$temp
  relhum<-climdata$relhum
  pk<-climdata$pres
  u<-climdata$windspeed
  hgt<-vegp$hgt
  # (1) Run NicheMapR
  PAIt<-apply(vegp$PAI,2,sum)
  LAI<-(vegp$pLAI*vegp$PAI)
  LAI<-apply(LAI,2,mean)
  pLAI<-LAI/PAIt
  LREFL<-mean(pLAI*vegp$refls+(1-pLAI)*vegp$reflp)
  soiltype<-as.character(soilp$Soil.type)
  DEP<-c(0,2.5,5,10,15,20,30,50,100,200)
  PE = rep(1.1, 19)
  KS = rep(0.0037, 19)
  BB = rep(4.5, 19)
  BD = rep(1.3, 19)
  DD = rep(2.65, 19)
  nmrout<-runNMR(climdata,prec,lat,long,1.95,hgt,2,PAIt,vegp$x,pLAI,
                 vegp$clump,vegp$refg,LREFL,0.95,DEP,altt,slope,aspect,
                 ERR,soiltype,PE,KS,BB,BD,DD,cap,hori,maxpool,rainmult,SoilMoist_Init)
  if (reqhgt < 0) {
    dep<-c(0,2.5,5,10,15,20,30,50,100,200)
    soilz<- -reqhgt*100
    dif<-abs(soilz-dep)
    o<-order(dif)
    dif1<-dif[o[1]]
    dif2<-dif[o[2]]
    tz1<-soilt[,o[1]+1]
    tz2<-soilt[,o[2]+1]
    tz<-tz1*(dif2/(dif1+dif2))+tz2*(dif1/(dif1+dif2))
    metout <- data.frame(obs_time=climdata$obs_time,Tref=tair,Tloc=tz,
                         tleaf=-999,RHref=relhum,RHloc=-999)
  } else if (reqhgt < (hgt+2)) {
    metout <- runmodelS(climdata,vegp,soilp,nmrout,reqhgt,lat,long,metopen,windhgt,surfwet,0.95)
  } else {
    metout <- data.frame(obs_time=climdata$obs_time,Tref=tair,Tloc=tair,
                         tleaf=-999,RHref=relhum,RHloc=relhum)
    warning("Height out of range. Output climate identical to input")
  }
  return(list(metout = metout, nmrout = nmrout))
}

