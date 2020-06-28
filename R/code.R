globalVariables(c("soilparams", "climvars", "weather"))
#' Calculates radiation absorbed by leaf
#'
#' @description Calculates the flux density of radiation absorbed by leaf.
#'
#' @param Rsw Incoming shortwave radiation (W / m2)
#' @param tme POSIXlt object of time(s)
#' @param tair air temperature (deg C)
#' @param tground ground temperature (deg C)
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees)
#' @param PAI Vector of plant area indices for each canopy layer
#' @param pLAI Proportion of plant area that is green vegetation
#' @param x the ratio of vertical to horizontal projections of leaf foliage
#' @param refls reflectivity of green vegetation to shortwave radiation
#' @param refw reflectivity of woody vegetation to shortwave radiation
#' @param refg reflectivity of ground to shortwave radiation
#' @param vegem emissivity of vegetation
#' @param skyem sky emissivity
#' @param dp proportion of `Rsw` that is diffuse radiation. If not provided, then calculated using [microctools::difprop()]
#' @param merid optional longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst optional value representing the time difference from the timezone meridian
#' @param clump clumpiness factor for canopy (0-1)
#'
#' @return a list of with the following components:
#' @return `aRsw` absorbed Shortwave radiation (W / m2)
#' @return `aRlw` absorbed longwave radiation (W / m2)
#' @return `ref` mean area-wighted reflectivity of vegetation (green and woody)
#' @import microctools
#' @export
leafabs <-function(Rsw, tme, tair, tground, lat, long, PAI, pLAI, x, refls, refw, refg, vegem, skyem, dp = NA,
                   merid = round(long/15, 0) * 15, dst = 0, clump = 0) {
  PAIc <- rev(cumsum(PAI))
  PAIg <- cumsum(PAI)
  jd<-jday(tme = tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  if (is.na(dp)) dp<-difprop(Rsw,jd,lt,lat,long,merid=merid,dst=dst)
  sa<-solalt(lt,lat,long,jd,merid=merid,dst=dst)
  sa2<-ifelse(sa<5,5,sa)
  ref<-pLAI*refls+(1-pLAI)*refw
  sunl<-psunlit(PAIc,x,sa,clump)
  mul<-radmult(x,sa2)
  aRsw <- (1-ref) * cansw(Rsw,dp,jd,lt,lat,long,PAIc,x,ref,merid=merid,dst=dst,clump=clump)
  aGround <- refg * cansw(Rsw,dp,jd,lt,lat,long,sum(PAI),x,mean(ref),merid=merid,dst=dst,clump=clump)
  aGround <- (1-ref) * cansw(aGround,dp,jd,lt,lat,long,PAIg,x,ref,merid=merid,dst=dst,clump=clump)
  aRsw <- sunl * mul * aRsw + (1-sunl) * aRsw + aGround
  aRlw <- canlw(tair, PAIc, 1-vegem, skyem = skyem, clump = clump)$lwabs
  return(list(aRsw=aRsw, aRlw=aRlw, ref=ref))
}
#' Calculates radiation emitted by leaf
#'
#' @description Calculates the flux density of radiation emitted by the leaf.
#'
#' @param tc temperature (deg C)
#' @param vegem emissivity of leaf
#'
#' @return Flux density of radiation emitted by leaf (W / m2)
#' @export
leafem <- function(tc, vegem) {
  eRlw <- vegem * 5.67*10^-8 * (tc + 273.15)^4
  eRlw
}
#' Thomas algorithm for solving simultanious heat fluxes
#'
#' @description `Thomas` implements the Thomas algorithm for solving simultanious heat
#' fluxes between soil / air layers.
#'
#' @param tc vector of soil and air temperatures (deg C) from previous timestep (see details)
#' @param tmsoil temperature (deg C) of deepest soil layer. Typically mean annual temperature
#' @param tair air temperature at reference height 2 m above canopy in current time step (deg C)
#' @param k vector of thermal conductances between layers (W / m^2 / K) (see details)
#' @param cd thermal heat capacity of layers (W / m^2 / K)
#' @param f forward / backward weighting of algorithm (see details)
#' @param X vector of temperatures to be added resulting from e.g. leaf heat fluxes or radiation
#' absorbed by top soil layer
#' @return a vector of temperatures (deg C) for each soil / air layer for current time step. The first value
#' is `tair` and the last `tmsoil`
#' @export
#' @details The vector `tc` must be ordered with reference air temperature first and the soil temperature
#' of the  deepest layer last. I.e. the length of the vector `tc` is the number of nodes + 2.
#' The vector `k` is the conductivity between each node and that diectly below it, the first value
#' representing conductivity between reference height and the top canopy node. I.e. the length
#' of the vector `k` is the number of nodes + 1. The vector `cd` is the heat storage at each node.
#' I.e. the length of the vector `cd` is the same as the number of nodes. The  weighting factor `f`  may
#' range from 0 to 1. If `f` = 0, the flux is determined by the temperature difference at the beginning
#' of the time step. If `f` = 0.5, the average of the old and new temperatures is used to compute heat flux.
#' If `f` = 1, fluxes are computed using only the new temperatures. The best value to use
#' for `f` is determined by considerations of numerical stability and accuracy and experimentation
#' may be required. If `f` = 0  more heat transfer between nodes is predicted than would actually
#' occur, and can therefore become unstable if time steps are too large. When `f` > 0.5,
#' stable solutions will always be obtained, but heat flux will be underestimated. The
#' best accuracy is obtained with `f` around 0.4, while best stability is at `f` = 1.
#' A typical compromise is `f` = 0.6.
Thomas <- function(tc, tmsoil, tair, k, cd, f = 0.6, X = 0) {
  m <- length(tc) - 2
  tn<-rep(0,m+2)
  tn[m+2]<-tmsoil
  tn[1]<-tair
  g<-1-f
  a <- c(0,0); b <- 0; cc <-0; d <- 0
  xx<-(2:(m+1))
  tc[xx]<-tc[xx]+g*X
  cc[xx]<- -k[xx]*f
  a[xx+1]<-cc[xx]
  b[xx]<-f*(k[xx]+k[xx-1])+cd
  d[xx]<-g*k[xx-1]*tc[xx-1]+(cd-g*(k[xx]+k[xx-1]))*tc[xx]+g*k[xx]*tc[xx+1]
  d[2]<-d[2]+k[1]*tn[1]*f
  d[m+1]<-d[m+1] + k[m+1] * f * tn[m+2]
  for (i in 2:m) {
    cc[i]<-cc[i]/b[i]
    d[i]<-d[i]/b[i]
    b[i+1]<-b[i+1]-a[i+1]*cc[i]
    d[i+1]<-d[i+1]-a[i+1]*d[i]
  }
  tn[m+1] <- d[m+1] / b[m+1]
  for (i in m:2) {
    tn[i]<-d[i]-cc[i]*tn[i+1]
  }
  xmn<-pmin(tc[xx],tc[xx-1],tc[xx+1])
  xmx<-pmax(tc[xx],tc[xx-1],tc[xx+1])
  tn[xx]<-ifelse(tn[xx]<xmn,xmn,tn[xx])
  tn[xx]<-ifelse(tn[xx]>xmx,xmx,tn[xx])
  tn[xx]<-tn[xx]+f*X
  tn
}
#' Thomas algorithm for solving simultanious vapour fluxes
#'
#' @description `ThomasV` implements the Thomas algorithm for solving simultanious vapour
#' fluxes between air layers.
#'
#' @param Vo a vector of air vapour concentrations for each canopy node in the previos timestep (mol fraction)
#' @param tn vector of air temperatures (deg C) for each canopy node in the current timestep (deg C)
#' @param pk atmospheric pressure (kPa)
#' @param theta Volumetric water content of the upper most soil layer in the current time step (m3 / m3)
#' @param thetap Volumetric water content of the upper most soil layer in the previous time step (m3 / m3)
#' @param relhum relative humidity (percentage) at reference height 2 m above canopy in current time step (percentage)
#' @param tair air temperature at reference height 2 m above canopy in current time step (deg C)
#' @param tsoil temperature of upper soil layer in current time step (deg C)
#' @param zth heightdifference between each canopy node and that directly below it. the first value is
#' the height difference between the lowest canopy node and the ground
#' @param gt vector of molar conductances between each canopy node at that directly below it (mol / m^2 / sec).
#' The first value is the conductivity between the ground and the lowest node, and the last value the
#' conductivity between the highest node and reference height.
#' @param Vflux Total vapour flux from leaves to air (mol /m^3)
#' @param f forward / backward weighting of algorithm (as for [Thomas()])
#' @param previn a list of model outputs form the previous timestep
#' @param soilp a list of soil parameters as returned by [soilinit()]
#' @return a vector of vapour concentrations expressed as mole fractions for each canopy node in the
#' current time step. The first value is that for the ground and the last value that at reference height
#' @import microctools
#' @export
#' @seealso [Thomas()]
ThomasV <- function(Vo, tn, pk, theta, thetap, relhum, tair, tsoil, zth, gt, Vflux, f = 0.6, previn, soilp) {
  m<-length(zth)
  ph<-phair(tn,pk)
  ea<-0.6108*exp(17.27*tair/(tair+237.3))*(relhum/100)
  eap<-0.6108*exp(17.27*previn$tair/(previn$tair+237.3))*(previn$relhum/100)
  Vair<-ea/pk
  rhsoil<-soilrh(theta,soilp$b,soilp$psi_e,soilp$Smax,tsoil)
  rhsoilp<-soilrh(thetap,soilp$b,soilp$psi_e,soilp$Smax,previn$soiltc[1])
  rhsoil[rhsoil>1]<-1
  rhsoilp[rhsoilp>1]<-1
  eas<-0.6108*exp(17.27*tsoil/(tsoil+237.3))*rhsoil
  easp<-0.6108*exp(17.27*previn$tsoil/(previn$tsoil+237.3))*rhsoilp
  Vsoil<-eas/pk
  Vo<-c(easp/previn$pk,Vo,eap/previn$pk)
  Vn<-Thomas(rev(Vo), Vsoil, Vair, rev(gt), rev(ph), f, Vflux)
  Vn<-rev(Vn)
}
#' Calculates wind profile for entire canopy
#'
#' @description calculates wind speed at any given point above or below canopy
#' from wind speed at reference height
#'
#' @param ui wind at reference height (m / s)
#' @param zi height of wind measurement (m)
#' @param zo height for which wind is required (m)
#' @param a attenuation coefficient as returned by [microctools::attencoef()]
#' @param PAI plant area index (m / m)
#' @param hgt height of canopy (m)
#' @param psi_m diabatic correction factor as return by [microctools::diabatic_cor()]
#' @param hgtg height of ground vegetation layer below canopy (m)
#' @param zm0 roughness length (m) of vegetation layer below canopy as returned by [microctools::roughlength()]
#' @return wind speed (m /s) at height `zo`
#' @import microctools
#' @export
#' @examples
#' zo<- c(0:200) / 100
#' uz <- windprofile(2, 2, zo, 2.3, 3, 2)
#' plot(zo ~ uz, type = "l", xlab = "wind speed", ylab = "height")
windprofile <- function(ui, zi, zo, a = 2, PAI, hgt, psi_m = 0, hgtg = 0.05 * hgt, zm0 = 0.004) {
  if (length(as.vector(hgt)) == 1) hgt <- hgt + zo * 0
  if (length(as.vector(psi_m)) == 1) psi_m <- psi_m + zo * 0
  if (length(as.vector(a)) == 1) a <- a + zo * 0
  if (length(as.vector(hgtg)) == 1) hgtg <- hgtg + zo * 0
  d <- zeroplanedis(hgt, PAI)
  zm <- roughlength(hgt, PAI, zm0)
  dg <- 0.65 * hgtg
  zmg <- 0.1 * hgtg
  ln1 <- log((zi - d) / zm) + psi_m
  uf <- 0.4 * (ui / ln1)
  # zo above canopy
  ln2 <- suppressWarnings(log((zo - d) / zm) + psi_m)
  ln2[ln2<0.55]<-0.55
  uo <- (uf / 0.4) * ln2
  # zo below canopy
  sel <- which(zo < hgt)
  ln2 <- log((hgt[sel] - d[sel]) / zm[sel]) + psi_m[sel]
  uo[sel] <- (uf[sel] / 0.4) * ln2
  # zo above 10% of canopy hgt, but below top of canopy
  sel <- which(zo > (0.1 * hgt) & zo < hgt)
  uo[sel] <- uo[sel] * exp(a[sel] * ((zo[sel] / hgt[sel]) - 1))
  # zo below 10% of canopy hgt
  sel <- which(zo <= (0.1 * hgt))
  uo[sel] <- uo[sel] * exp(a[sel] * (((0.1 * hgt[sel]) / hgt[sel]) - 1))
  ln1 <- log((0.1 * hgt[sel]) / zmg[sel]) + psi_m[sel]
  uf <- 0.4 * (uo[sel] / ln1)
  ln2 <- log(zo[sel] / zmg[sel]) + psi_m[sel]
  uo[sel] <- (uf / 0.4) * ln2
  uo
}
#' Calculates wind profile for individual canopy layers
#'
#' @description Used for calculating wind speed profile in canopies with variable
#' `PAI`
#' @param uh wind speed at top of canopy layer (m /s)
#' @param z height of canopy layer node (m)
#' @param hgt height to top of canopy layer (m)
#' @param PAI Plant Area Index of canopy layer (m / m)
#' @param x the ratio of vertical to horizontal projections of leaf foliage
#' @param lw mean leaf width (m)
#' @param cd drag coefficient
#' @param iw turbulence intensity
#' @param phi_m diabatic correction factor
#' @param edgedist optional numeric value indicating distance (m) to open
#' habitat (see details)
#' @param uref wind at reference height (m) (see details)
#' @param zref height (m) of `uref`
#' @details if 'edgedist' is not NA, then horizontal wind component is also added. Here,
#' the wind profile of a reference grass surface is calaculate, and at any given height `z`
#' an attenuated wind speed inside thew canopy calculated, with the degree of attenuation
#' determined by foliage density, distance from edge and the ratio of horizontal to vertical
#' wind  speed. If the the attenuated horizontal wind speed exceeds the wind speed due
#' to vertical attenuation, the horizontal wind speed is returned. The parameter `uref`
#' represents the wind speed at height `zref` above the reference grass surface. If `edgedist`
#' is set to NA (the default) the horizontal wind component and hence `uref` and `zref` are
#' ignored.
#' @return wind speed at height of canopy node (m / s)
#' @import microctools
#' @export
#' @examples
#' # ==== Generate plant area index values
#' m <- 100
#' hgt <- 5
#' z<-c(1:m) * (hgt / m)
#' PAI <- microctools::PAIgeometry(m, 3, 7, 70)
#' plot(z~PAI, type = "l")
#' cPAI <- cumsum(PAI)
#' # ==== Calculate wind speed at top of canopy
#' uref <- 2
#' a <- microctools::attencoef(hgt, 3, 1)
#' uh <- windprofile(uref, hgt + 2, hgt, a, 3, hgt)
#' # === Calculate canopy profile (near edge)
#' uz1 <- 0
#' for (i in m:1) {
#'   uz1[i] <- windcanopy(uh, z[i], z[i] + 0.05, cPAI[i], edgedist = 20, uref = uref)
#'   uh <- windcanopy(uh, z[i] - 0.05, z[i] + 0.05, cPAI[i], edgedist = 20, uref = uref)
#' }
#' # === Calculate canopy profile (far from edge)
#' uh <- windprofile(uref, hgt + 2, hgt, a, 3, hgt)
#' uz2 <- 0
#' for (i in m:1) {
#'   uz2[i] <- windcanopy(uh, z[i], z[i] + 0.05, cPAI[i], edgedist = 500, uref = uref)
#'  uh <- windcanopy(uh, z[i] - 0.05, z[i] + 0.05, cPAI[i], edgedist = 500, uref = uref)
#' }
#' plot(z ~ uz1, type = "l", xlab = "wind speed", ylab = "height", xlim = c(0,1.5),
#'      col = rgb(1,0,0,0.5), lwd = 2)
#' par(new = TRUE)
#' plot(z ~ uz2, type = "l", xlab = "", ylab = "", xlim = c(0,1.5),
#'      col = rgb(0,0,1,0.5), lwd = 2)
windcanopy <- function(uh, z, hgt, PAI = 3, x = 1, lw = 0.05, cd = 0.2,
                       iw = 0.5, phi_m  = 1, edgedist = NA, uref, zref = hgt + 2) {
  a <- attencoef(hgt, PAI, x, lw, cd, iw, phi_m)
  uz <- uh * exp(a * ((z / hgt) - 1))
  # horizontal wind component
  if (is.na(edgedist) == F) {
    a2 <- attencoef(hgt, PAI, 1/x, lw, cd, iw, phi_m)
    ur <- windprofile(uref, zref, z, 0.7989537, 1.146, 0.12)
    uwr<-suppressWarnings(2.5*log((z-0.08)/0.01476))
    uwr[uwr < 1] <- 1; uwr[is.na(uwr)] <- 1
    edr <- (hgt - edgedist/uwr) / hgt
    uhr <- ur * exp(a2 * (edr - 1))
    uz <- pmax(uz,uhr,na.rm=T)
  }
  uz
}
#' Calculates temperature above canopy
#'
#' @description calaculates the temperature at the top of the canopy based
#' on the standard logarithmic height profile
#' @param tz temperature at height `zu` above the canopy (deg C)
#' @param uz wind speed at height `zu` above the canopy (m / s)
#' @param zu height of `tz` and `uz` (m)
#' @param zo height (m) for which temperature is required (see details)
#' @param H sensible heat flux density (W / m^2). See details.
#' @param hgt height of the canopy (m)
#' @param PAI total plant area index of the canopy. Used to calculate roughness lengths.
#' @param zm0 roughness length governing momentum transfer of ground vegetation
#' @param pk atmospheric pressure (kPa)
#' @param psi_h diabatic correction factor for heat transfer
#' @details Estimation of `H` requires estimation of temperature,
#' so most either be derived by iteration of taken form the previous timestep.
#' `H` is given by the net energy balance equation: `H` = `Rabs` - `Rem` - `L` - `G`
#' where `Rabs` is absorbed radiation, `Rem` emitted radiation, `L` Latent heat
#' exchange and `G` ground heat flux. This function is not valid for temperatures
#' below canopy so `hgt` must be lower than `zi`.
#' @export
#' @examples
#' abovecanopytemp(11, 2, 2, 1.5, 500, 1, 3)
#' abovecanopytemp(11, 2, 2, 0.5, 500, 0.25, 3)
abovecanopytemp <- function(tz, uz, zu, zo, H, hgt, PAI, zm0 = 0.004, pk = 101.3, psi_h = 0) {
  if (min(zo) < hgt) stop("zi must be greater or equal to hgt")
  d <- zeroplanedis(hgt, PAI)
  zm <- roughlength(hgt, PAI, zm0)
  zh <- 0.2 * zm
  uf <- (0.4 * uz) /  (log((zu - d) / zm) + psi_h)
  ph <- phair(tz, pk)
  cp <- cpair(tz)
  m <- H / (0.4 * ph * cp * uf)
  tref <- tz + m * (log((zu - d) / zh) + psi_h)
  tc <- tref - m * (log((zo - d) / zh) + psi_h)
  tc
}
#' Calculates leaf temperature of canopy layers
#'
#' @description Calculates leaf temperature from radiation fluxes and reference air
#' and ground tempertaure
#'
#' @param tair air temperature at two metres above canopy (deg C)
#' @param relhum relative humidity at two metres above canopy (Percentage)
#' @param pk air pressure at two metres above canopy (kPa)
#' @param timestep duration of time step (s)
#' @param gt heat conductance by turbulent convection (mol / m^2 / s)
#' @param gha leaf-air heat conductance (mol / m^2 / s)
#' @param gv leaf-air vapour conductance (mol / m^2 / s)
#' @param Rabs Flux density of absorbed radiation as returned by [leafabs()] (W/m2)
#' @param previn list of values from previous timestep
#' @param vegp list of vegetation paramaters (see e.g. `vegparams` dataset)
#' @param soilp list Soil paramaters (see e.g. `soilparams` dataset)
#' @param theta volumetric soil moisture fraction of top soil layer (m^3 / m^3)
#' @param zlafact numeric value indicating how close to leaves air temperatures are
#' needed for (1 - average leaf-air distance, 0.5 = half average leaf-air distance etc.).
#' Must be greater than 1.
#' @param surfwet Fraction of surface area acting like a free‐water surface
#' @return a list with the following elements:
#' @return `tn` air temperature of each layer (deg C)
#' @return `tleaf` leaf temperature of each layer (deg C)
#' @return `ea` vapour pressure (kPa)
#' @return `gtt` heat conductance to reference height two m above canopy
#' @return `Rem` emitted radiation (W / m^2)
#' @return `H` Sensible heat flux from leaf to air (W / m^2)
#' @return `L` Latent heat of vapourisation from leaf to air (W / m^2)
#' @import microctools
#' @export
#' @details `leaftemp` computes the average leaf and air temperature of each canopy layer based on
#' radiation and evaorative fluxes and reference air and ground temperature. The function
#' automatically determines whether heat storage should be considered, based the specific heat
#' capacity of the vegetation layer and the time step of the model.
#' @examples
#' # Generate paramaters for function:
#' tme <- as.POSIXlt(0, origin = "2020-05-04 12:00", tz = "GMT")
#' previn <- paraminit(20, 10, 10, 15, 2, 80, 11, 500)
#' vegp <- microctools::habitatvars(3, 50, -5, tme, m = 20)
#' soilp <- soilinit("Loam")
#' z<-c((1:20) - 0.5) / 20 * vegp$hgt
#' # run function (setting conductances in current time step to same as in previous):
#' ltemp <- leaftemp(11, 80, 101.3, 60, previn$gt, previn$gha, previn$gv, previn$Rabs, previn,
#'                   vegp, soilp, 0.3)
#' plot(z ~ ltemp$tleaf, type = "l", xlab = "Leaf temperature", ylab = "Height")
leaftemp <- function(tair, relhum, pk, timestep, gt, gha, gv, Rabs, previn, vegp, soilp,
                     theta, zlafact = 1, surfwet = 1) {
  edf<-function(ea1,ea2,tc) {
    tk<-tc+273.15
    es<-0.6108*exp(17.27*tc/(tc+237.3))
    y<-ifelse(ea1>ea2,ea1-ea2,0)
    y<-ifelse(ea2>es,0,y)
    y
  }
  z<-previn$z
  lambda <- (-42.575*tair+44994)*surfwet
  # Sort out thicknesses
  m<-length(gt)-1
  zth<-c(z[2:m]-z[1:(m-1)],vegp$hgt-(z[m]+z[m-1])/2)
  zla<-mixinglength(zth,vegp$PAI,vegp$x,vegp$lw)*0.5*zlafact
  # Sort out conductivitites
  gt<-0.5*gt+0.5*previn$gt
  gv1 <- 1 / ((1/previn$gv) + (1/previn$gha))
  gv2 <- 1 / ((1/gv) + (1/gha))
  gv<-0.5*gv1+0.5*gv2
  gha<-0.5*gha+0.5*previn$gha
  mtref<-0.5*tair+0.5*previn$tair
  mrh<-0.5*relhum+0.5*previn$relhum
  igtt<-rep(1/gt[m+1],m+1)
  igtt2<-1/gt[1]
  for (i in m:1) igtt[i]<-igtt[i+1]+1/gt[i]
  for(i in 2:m) igtt2[i]<-igtt2[i-1]+1/gt[i]
  gtt<-1/igtt[-1]
  gtt2<-1/igtt2
  zref<-(vegp$hgt+2)-z
  mpk<-0.5*previn$pk+0.5*pk
  # Vapour pressure
  esj<-0.6108*exp(17.27*previn$tc/(previn$tc+237.3))
  eaj<-(previn$rh/100)*esj
  estl<-satvap(previn$tleaf,ice = TRUE)
  esref<-satvap(mtref,ice=TRUE)
  eref<-(mrh/100)*esref
  delta <- 4098*(0.6108*exp(17.27*previn$tleaf/(previn$tleaf+237.3)))/(previn$tleaf+237.3)^2
  rhsoil<-soilrh(theta,soilp$b,-soilp$psi_e,soilp$Smax, previn$soiltc[1])
  esoil<-rhsoil*satvap(previn$soiltc[1],ice = TRUE)
  # Test whether steady state
  PAIm<-2*vegp$PAI/zth
  gv2<-(PAIm/zth)*gv
  test<-pmax(timestep*gtt/zref,timestep*gv/zla,timestep*gtt2/z)
  sel<-which(test>1)
  btm<-(1/timestep)+0.5*(gtt/zref+gtt2/z+gv/zla)
  ae<-eaj+0.5*((gtt/zref)*edf(eref,eaj,previn$tc)+(gtt2/z)*edf(esoil,eaj,previn$tc)+
                 (gv/zla)*edf(estl,eaj,previn$tc))/btm
  be<-(0.25*gv*delta)/btm
  ae2<-(gtt*eref+gtt2*esoil+gv2*estl)/(gtt+gtt2+gv2)
  be2<-(0.5*delta)/(gtt+gtt2+gv2)
  ae[sel]<-ae2[sel]
  be[sel]<-be2[sel]
  # Air temperature
  # Test whether steady state
  test<-pmax(timestep*gtt/zref,timestep*gha/zla,timestep*gtt2/z)
  # ~Transient
  ae2<-ifelse(ae>estl,estl,ae)
  ae3<-ifelse(ae>esoil,esoil,ae)
  sel<-which(test>1)
  ph<-phair(previn$tc,previn$pk)
  cp<-cpair(previn$tc)
  vden<-vegp$thickw*vegp$PAI
  ma<-(timestep*PAIm*2)/(cp*ph*(1-vden))
  K1<-gtt*cp/zref; K2<-gtt2*cp/z; K3<-gha*cp*PAIm/zla
  btm<-1+0.5*ma*(K1+K2+K3)
  aL<-(previn$tc+0.5*ma*(K1*mtref+K2*previn$soiltc[1]+K3*previn$tleaf))/btm
  bL<-(0.25*ma*K3)/btm
  # ~Steady state
  K1<-gtt[sel]*cp[sel]; K2<-gtt2[sel]*cp[sel]; K3<-gha[sel]*cp[sel]
  aL[sel]<-(K1*mtref+K2*previn$soiltc[1]+K3*previn$tleaf[sel])/(K1+K2+K3)
  bL[sel]<-(0.5*K3)/(K1+K2+K3)
  bL<-ifelse(bL<0,0,bL)
  # Sensible Heat
  aH<-cp*gha*(previn$tleaf-aL)
  bH<-cp*gha*(0.5-0.5*bL)
  bH[bH<0]<-0
  # Latent Heat
  aX<-(lambda*gv/mpk)*edf(estl,ae,previn$tc)
  bX<-(lambda*gv/mpk)*(delta-be)
  bX[bX<0]<-0
  # Radiation
  Rabs<-0.5*Rabs+0.5*previn$Rabs
  sb<-5.67*10^-8
  aR<-vegp$vegem*sb*(previn$tleaf+273.15)^4
  bR<-vegp$vegem*sb*2*(previn$tleaf+273.15)^3
  # Leaf temperature
  Ch<-vden*vegp$cpw*vegp$phw
  ml<-(timestep*PAIm)/Ch
  dTL<-(ml*(Rabs-aR-aX-aH))/(zla+ml*(bR+bX+bH))
  # check whether steady state
  dTL2<-(Rabs-aR-aX-aH)/(bR+bX+bH)
  sel<-which(abs(dTL2)<abs(dTL))
  dTL[sel]<-dTL2[sel]
  # set limits to air temperature
  tn<-aL+bL*dTL
  tleaf<-previn$tleaf+dTL
  sel<-which(tleaf>tn)
  tmx<-rep(tair,length(tn[sel]))+15
  tmx<-pmax(tmx,tleaf[sel])
  tn[sel]<-ifelse(tn[sel]>tmx,tmx,tn[sel])
  sel<-which(tleaf<tn)
  tmn<-rep(tair,length(tn[sel]))-5
  tmn<-pmin(tmn,tleaf[sel])
  tn[sel]<-ifelse(tn[sel]<tmn,tmn,tn[sel])
  # check whether saturated or not
  eam<-ae+be*dTL
  ean<-eaj+2*(eam-eaj)
  es<-0.6108*exp(17.27*tn/(tn+237.3))
  ean<-ifelse(ean>es,es,ean)
  # Calculate minimum potential vapour pressure
  eamn <- (relhum/100)*0.6108*exp(17.27*tair/(tair+237.3))
  ean<-ifelse(ean<eamn,eamn,ean)
  # Dew point temperature
  a<-log(ean/0.6108)
  Tdew<- -(237.3*a)/(a-17.27)
  dT1 <- -(previn$tleaf+dTL - Tdew)
  dT1 <- ifelse(dT1<0,0,dT1)
  estl <- 0.6108*exp(17.27*tleaf/(tleaf+237.3))
  Lc <- lambda * gv * (estl - ean) / mpk
  dT2 <- -(ml/zla)*Lc
  dT2 <- ifelse(dT2 < 0, 0, dT2)
  dT <- pmin(dT1,dT2)
  tleaf<-tleaf+dT
  dTL<-tleaf-previn$tleaf
  # Set limits to leaf temperature
  dTmx <- previn$tc + 20.2 - previn$tleaf
  dTmn <- previn$tc - 4.5 - previn$tleaf
  dTL[dTL>dTmx] <- dTmx[dTL>dTmx]
  dTL[dTL<dTmn] <- dTmn[dTL<dTmn]
  tn<-aL+bL*dTL
  tleaf<-previn$tleaf+dTL
  return(list(tn=tn, tleaf=previn$tleaf+dTL, ea=ean, gtt=gtt, Rem=aR+bR*dTL,
              H=aH+bH*dTL, L=aX+bX*dTL, esoil = esoil))
}
#' Initialise paramaters for first time step of model
#'
#' @description Generates a set of climate and conductivity parameters for running the
#' first time step of the model
#'
#' @param m number of canopy layer nodes
#' @param sm number of soil layers
#' @param hgt height of canopy (m)
#' @param tair air temperature at reference height (deg C)
#' @param u wind speed at reference height (m/s)
#' @param relhum relative humidity at reference height (percentage)
#' @param tsoil stable temperature below deepest soil layer. Usually ~mean annual temperature
#' (deg C). See details.
#' @param Rsw total incoming shortwave radiation (W / m^2)
#' @return a list with the following elements:
#' @return `tc` a vector of air temperatures for each canopy layer (deg C)
#' @return `soiltc` a vector of airsoil temperatures for each soil layer (deg C)
#' @return `tleaf` a vector of leaf temperatures for each canopy layer (deg C)
#' @return `tabove` initially set above canopy temperature, here set as `tair`
#' @return `uz` a vector of wind speeds for each canopy layer (m/s)
#' @return `z` a vector of canopy node heights (m)
#' @return `sz` a vector of soil node depths (m)
#' @return `zabove` height of `tabove`. Here set to `hgt`
#' @return `rh` a vector of relative humidities
#' @return `relhum` relative humidity at 2 m above canopy (percentage)
#' @return `tair` air temperature at 2 m above canopy (deg C)
#' @return `tsoil` temperature below deepest soil layer. Same as input
#' @return `pk` pressure at 2 m above canopu (kPA)
#' @return `Rabs` Absorbed radiation (W / m^2)
#' @return `gt` Conductivity in air of each canopy layer node (mol/m^2/sec)
#' @return `gv` Leaf conductivity to vapour loss for each canopy layer node  (mol/m^2/sec)
#' @return `gha` Conductivity between air and leaf for each canopy layer node (mol/m^2/sec)
#' @return `H` Sensible heat flux (W / m^2). Here set to zero
#' @return `L` Latent heat flux (W / m^2). Here set to zero
#' @return `G` Ground heat flux (W / m^2). Here set to zero
#' @return `Rswin` Incoming shortwave radiation. Here set to `Rabs`
#' @return `Rlwin` Incoming longwave radiation. Here set to 0.2 * `Rabs`
#' @return `psi_m` Diabatic correction factor. Here set to zero
#' @importFrom stats spline
#' @import microctools
#' @export
#' @examples
#' paraminit(20, 10, 10, 15, 2, 80, 11, 500)
#'
#' @details All values are approximate. Values for `tc` and `tsoil` are derived by
#' linear intepolation between `tair` and `tsoil`. Values for `Rabs` are derived from
#' `Rsw` but attenuate through the canopy. Values for `tleaf` are derived
#' from `tc` and `Rabs`. Values for `rh` are the same as `relhum`. Values for `gt`
#' `gv` and `gha` are typical for decidious woodland with wind above canopy at 2 m/s.
#' `gt` is scaled by canopy height and `m` (and hence distance between nodes). The first
#' value represents conductivity between the ground and the lowest canopy node. The last
#' value represents conductivity between the air at 2 m above canopy and the highest
#' canopy node.
paraminit <- function(m, sm, hgt, tair, u, relhum, tsoil, Rsw) {
  tcb <- spline(c(1,2), c(tsoil, tair), n = m + sm)
  tc <- tcb$y[(sm + 1): (m + sm)]
  soiltc <- rev(tcb$y[1:sm])
  rh <- rep(relhum, m)
  Rabs <- spline(c(1, 2), c(0.2 * Rsw + 20, Rsw + 20), n = m)$y
  tleaf <- tc + 0.01 * Rabs
  gt <- rep(50, m + 1) * (m / hgt) * (2 / m)
  gt[1] <- 2 * gt[1]
  gt[m+1] <- gt[m+1] * (hgt / m) * 0.5
  gv <- spline(c(1, 2), c(0.25, 0.32), n = m)$y
  gha <- spline(c(1, 2), c(0.13, 0.19), n = m)$y
  z<-c((1:m)-0.5)/m*hgt
  sz<-2/sm^2.42*c(1:sm)^2.42
  uz <- (c(1:m)/m)*u
  return(list(tc = tc, soiltc = soiltc, tleaf = tleaf, tabove = tair, uz = uz, z = z, sz = sz,
              zabove = hgt, rh = rh, relhum = relhum, tair = tair, tsoil = tsoil,
              pk = 101.3, Rabs = Rabs, gt = gt, gv = gv, gha = gha, H = 0, L = rep(0, m),
              G = 0, Rswin = Rabs, Rlwin = 0.2 * Rabs, psi_m = 0))
}
#' Returns soil parameters for a given soil type
#'
#' @description Returns soil parameters needed to run the microclimate model
#' for a given soil type
#'
#' @param soiltype one of `Sand`, `Loamy sand`, `Sandy loam`, `Loam`, `Silt`,
#' `Silt loam`, `Sandy clay loam`, `Clay loam`, `Silty clay loam`, `Sandy clay`,
#' `Silty clay` or `Clay`.
#' @param m number of soil nodes (see details)
#' @param sdepth depth of deepest soil node (m)
#' @param reqdepth optional depth (m) for which temperature is required (m)
#' @return a list with the following items:
#' @return `Soil.type` description of soil type
#' @return `Smax` Volumetric water content at saturation (m^3 / m^3)
#' @return `Smin` Residual water content (m^3 / m^3)
#' @return `alpha` Shape parameter of the van Genuchten model (cm^-1)
#' @return `n` Pore size distribution parameter (dimensionless, > 1)
#' @return `Ksat` Saturated hydraulic conductivity (cm / day)
#' @return `Vq` Volumetric quartz content of soil
#' @return `Vm` Volumetric mineral content of soil
#' @return `Vo` Volumetric organic content of soil
#' @return `Mc` Mass fraction of clay
#' @return `rho` Soil bulk density (Mg / m^3)
#' @return `b` Shape parameter for Campbell model (dimensionless, > 1)
#' @return `psi_e` Matric potential (J / m^3)
#' @return `z` depth (m) of soil nodes (see details)
#' @details the depth of the soil nodes are automatically calculated so
#' as to ensure to ensure progressively thicker soil layers. If `reqdepth` is not NA, the soil node
#' nearest to `reqdepth` is set to `reqdepth`.
#' @import microctools
#' @export
#' @examples
#' soilinit("Loam")
soilinit <- function(soiltype, m = 10, sdepth = 2, reqdepth = NA) {
  sel<-which(soilparams$Soil.type==soiltype)
  soilp<-soilparams[sel,]
  z<-(sdepth/m^1.2)*c(1:m)^1.2
  if (is.na(reqdepth) == F) z[abs(z-reqdepth)==min(abs(z-reqdepth))][1]<-reqdepth
  soilp<-as.list(soilp)
  soilp$z<-z
  return(soilp)
}
#' Run canopy model for a single time step
#' @description run the above or below-canopy model for a single timestep
#'
#' @param climvars a of climate variables needed to run the run the model for one timestep (see dataset [climvars()])
#' @param previn a list of model outputs form the previous timestep as returned initially by [paraminit()]
#' @param vegp a list of vegetation parameters as returned by [microctools::habitatvars()]
#' @param soilp a list of soil parameters as returned by [soilinit()]
#' @param timestep length of model timestep (s)
#' @param tme POSIXlt object of the date and time of the current time step
#' @param lat latitude of location (decimal degrees)
#' @param long lonitude of location (decimal degrees)
#' @param edgedist distance to open ground (m)
#' @param reqhgt optional height for which temperature is required (see details)
#' @param sdepth depth of deepest soil node (m)
#' @param zu height above ground of reference climate measurements (m)
#' @param theta volumetric water content of upper most soil layer in current time step (m^3 / m^3)
#' @param thetap volumetric water content of upper most soil layer in previous time step (m^3 / m^3)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param n forward / backward weighting for Thomas algorithm (see [Thomas()])
#' @param metopen optional logical indicating whether the wind measurement used as an input to
#' the model is from a nearby weather station located in open ground (TRUE) or above the canopy
#' for which temperatures are modelled (FALSE - see details)
#' @param windhgt height above ground of wind measurement. If `metopen` is FALSE, must be above
#' canopy.
#' @param zlafact numeric value indicating how close to leaves air temperatures are
#' needed for (1 - average leaf-air distance, 0.5 = half average leaf-air distance etc.).
#' Must be greater than 1.
#' @param surfwet fraction of vegetation surface acting like a free-water surface
#' @return a list of of model outputs for the current timestep with the same format as `previn`
#' @import microctools
#' @export
#' @details model outputs are returned for each canopy node, with the number of canopy nodes (m)
#' determined by `previn`. Canopy nodes are spaced at equal heights throughout the canopy as in the
#' example. If `reqhgt` is set, and below the height of the canopy, the node nearest to that
#' height is set at the value specified. The returned value `tabove` is then the temperature
#' at te top of the canopy. If `reqhgt` is above canopy, nodes are calculated automatically,
#' but `tabove` is the temperature at height `reqhgt`. If temperatures below ground are
#' needed, the depth can be set using [soilinit()]. The wind profile of the canopy depends on
#' the nature of the canopy itself, and often available wind measurements are for a nearby
#' weather station located in open ground where it is possible that the height of the wind
#' measurement is below the height of the canopy being studied. When `metopen` is TRUE, the
#' wind profile of reference grass surface is used to derive estimates for two metres above the
#' canopy of interest. When `metopen` is FALSE, `windhgt` must be above canopy and the profile above
#' the canopy being studied is used.
#'
#' @examples
#' # Create initail parameters
#' tme <- as.POSIXlt(0, origin = "2020-05-04 12:00", tz = "GMT")
#' previn <- paraminit(20, 10, 10, 15, 2, 80, 11, 500)
#' vegp <- microctools::habitatvars(4, 50, -5, tme, m = 20)
#' z<-c((1:20)-0.5)/20*vegp$hgt
#' soilp<- soilinit("Loam")
#' climvars <- list(tair=16,relhum=90,pk=101.3,u=2.1,tsoil=11,skyem=0.9,Rsw=500,dp=NA)
#' # Run model 100 times for current time step
#' for (i in 1:100) {
#'   plot(z ~ previn$tc, type = "l", xlab = "Temperature", ylab = "Height", main = i)
#'   previn <- runonestep(climvars, previn, vegp, soilp, 60, tme, 50, -5)
#' }
runonestep <- function(climvars, previn, vegp, soilp, timestep, tme, lat, long, edgedist = 1000,
                       sdepth = 2, reqhgt = NA, zu = 2, theta = 0.3, thetap = 0.3, merid = 0,
                       dst = 0, n = 0.6, metopen = TRUE, windhgt = 2, zlafact = 1, surfwet = 1) {
  # =============   Unpack climate variables ========== #
  m <- length(previn$tc)
  tair<-climvars$tair; relhum<-climvars$relhum; pk<-climvars$pk; u<-climvars$u
  tsoil<-climvars$tsoil; skyem<-climvars$skyem; Rsw<-climvars$Rsw; dp<-climvars$dp
  H <- previn$H
  # ========== Calculate baseline variables ============== #
  tc<-previn$tc; ppk<-previn$pk; hgt<-vegp$hgt
  ph<-phair(tc,ppk); pha<-phair(tair,pk) # molar density of air
  cp<-cpair(tc) # specific heat of air
  lambda <- -42.575*tc+44994 # Latent heat of vapourisation (J / mol)
  # Adjust wind to 2 m above canopy
  if (metopen) {
    if (windhgt != 2) u <- u*4.87/log(67.8*windhgt-5.42)
    u2<-u*log(67.8*(hgt+2)-5.42)/log(67.8*2-5.42)
  } else {
    u2 <- u
    if (windhgt != (2+hgt)) u2 <- windprofile(u, windhgt, hgt+2, a = 2, sum(vegp$PAI), hgt, previn$psi_m)
  }
  u2[u2 < 0.5] <- 0.5
  # Generate heights of nodes
  z<-c((1:m)-0.5)/m*hgt
  if (is.na(reqhgt) == F) z[abs(z-reqhgt)==min(abs(z-reqhgt))][1]<-reqhgt
  zt<-z[2:(m)]-z[1:(m-1)] # difference in height between layers
  #  Set z above
  zabove <- vegp$hgt
  if (is.na(reqhgt) == F) {
    if (reqhgt > hgt) zabove <- reqhgt
  }
  # ========== Calculate diabatic correction factors ============== #
  d<-zeroplanedis(hgt,sum(vegp$PAI))
  zm<-roughlength(hgt, sum(vegp$PAI), vegp$zm0)
  uf<-(0.4*u2)/log(((hgt+2)-d)/zm)
  abod<-diabatic_cor(tair, pk, H, uf, (hgt+2), d)
  cand<-diabatic_cor_can(tc, previn$uz, z, vegp$PAI, vegp$x, vegp$lw)
  psi_m<-abod$psi_m; psi_h<-abod$psi_h; phi_m<-cand$phi_m; phi_h<-cand$phi_h
  # Calculate temperatures and relative humidities for top of canopy
  tcan <- abovecanopytemp(tair,u2,zu+hgt,zabove,H,hgt,sum(vegp$PAI),vegp$zm0,pk,psi_h)
  # Set limits to tcan
  tcan<-ifelse(tcan>(tair+15),(tair+15),tcan)
  tcan<-ifelse(tcan<(tair-4.5),(tair-4.5),tcan)
  # Adjust relative humidity
  ea<-0.6108*exp(17.27*tair/(tair+237.3))*(relhum/100)
  eas<-0.6108*exp(17.27*tcan/(tcan+237.3))
  relhum<-(ea/eas)*100; relhum[relhum>100]<-100
  # ========== Calculate wind speed and turbulent conductances ======== #
  ac<- attencoef(hgt,sum(vegp$PAI),vegp$x,vegp$lw,vegp$cd,vegp$iw,phi_m)
  uz<-rep(0,m); gt<-rep(0,m)
  uh<-windprofile(u2,hgt+2,hgt,ac[m],sum(vegp$PAI),hgt,psi_m,vegp$hgtg,vegp$zm0)
  uz[m]<-windcanopy(uh,z[m],hgt,sum(vegp$PAI),vegp$x,vegp$lw,vegp$cd,vegp$iw[m],
                    phi_m,edgedist,u,zu)
  gt[m]<-gcanopy(uh,z[m],z[m-1],tc[m],tc[m-1],hgt,sum(vegp$PAI),vegp$x,vegp$lw,
                 vegp$cd,vegp$iw[m],phi_m,pk)
  tc<-c(previn$soiltc[1],tc)
  # ========== Calculate canopy turbulences ========== #
  for (i in (m-1):1) {
    zi<-z[i+1]+0.5*zt[i]; zo<-z[i+1]-0.5*zt[i]
    uh<-windcanopy(uh,zo,zi,sum(vegp$PAI[1:i]),vegp$x,vegp$lw,vegp$cd,vegp$iw[i],
                   phi_m,edgedist,u,zu)
    if (i>1) {
      z0<-z[i-1]
    } else z0<-0
    uz[i]<-windcanopy(uh,zi,zo,sum(vegp$PAI[1:i]),vegp$x,vegp$lw,vegp$cd,
                      vegp$iw[i],phi_m,edgedist,u,zu)
    gt[i]<-gcanopy(uh,z[i],z0,tc[i+1],tc[i],zi,sum(vegp$PAI[1:i]),vegp$x,vegp$lw,
                   vegp$cd,vegp$iw[i],phi_m,pk)
  }
  # ========== Calculate conductivity to top of canopy and merge ============ #
  gt[m+1]<-gcanopy(uh,hgt,z[m],tcan,tc[i+1],hgt,vegp$PAI[m]*0.5,vegp$x,vegp$lw,
                   vegp$cd,vegp$iw[m],phi_m,pk)
  # Turbulent air conductivity and layer merge
  lmm<-layermerge(z,gt,hgt,timestep,ph)
  mrge<-lmm$mrge; gtx<-lmm$gtx; zth<-lmm$zth
  tc <- tc[-1]
  # ========== Calculate absorbed radiation =========== #
  Rabss<-leafabs(Rsw,tme,tair,previn$soiltc[1],lat,long,vegp$PAI,vegp$pLAI,vegp$x,
                 vegp$refls,vegp$refw,vegp$refg,vegp$vegem,skyem,dp,merid,dst,vegp$clump)
  Rabs<-Rabss$aRsw+Rabss$aRlw
  # ============= Conductivities =============== #
  # Vapour conductivity
  gv<-layercond(Rabss$aRsw/(1-Rabss$ref),vegp$gsmax,vegp$q50)
  # Leaf conductivity
  dtc<-previn$tleaf-previn$tc
  gha<-1.41*gforcedfree(vegp$lw*0.71,uz,tc,dtc,pk)
  tln<-leaftemp(tcan,relhum,pk,timestep,gt,gha,gv,Rabs,previn,vegp,soilp,theta,zlafact = zlafact, surfwet = surfwet)
  eaj<-0.6108*exp(17.27*tc/(tc+237.3))*(previn$rh/100)
  Vo<-eaj/previn$pk
  # =============== Soil conductivity =========== #
  sm<-length(previn$soiltc)
  cdk<-soilk(timestep,theta,soilp)
  sz<-soilp$z
  # conductivity and specific heat
  vden<-(vegp$PAI*vegp$thickw)
  mult<-1-vden
  zla <- mixinglength(vegp$hgt, vegp$PAI, vegp$x, vegp$lw)*0.5*zlafact
  X<-tln$tn-tc # Heat to add
  TT<-cumsum((ph/gt[1:m])*(z-c(0,z[1:(m-1)])))
  # Soil heat
  vdef<-tln$esoil-tln$ea[1]
  if(vdef>0) {
    mm<-(lambda[1]*ph[1]*z[1])/pk
    delta<-4098*(0.6108*exp(17.27*previn$soiltc[1]/(previn$soiltc[1]+237.3)))/(previn$soiltc[1]+237.3)^2
    btm<-1+0.5*(mm/cdk$cd[1])*delta
    Xs<-((1-vegp$refg)*(Rabss$aRsw[1]/cdk$cd[1])-(mm/cdk$cd[1])*(tln$esoil-tln$ea[1]))/btm
  } else {
    Xs<-(1-vegp$refg)*(Rabss$aRsw[1]/cdk$cd[1])
  }
  mm<-ifelse(Xs<0,-1,1)
  Xsmx<-((1-vegp$refg)*Rabss$aRsw[1])/(gha[1]*cp[1])+(tln$tn[1]-previn$soiltc[1])
  Xs<-ifelse(abs(Xs)>abs(Xsmx),Xsmx,Xs)
  Xs<- abs(Xs)*mm
  # Canopy air layer not in equilibrium with above canopy
  if (length(unique(lmm$u))>1) {
    lav<-layeraverage(lmm,tc,vegp$hgt,gha,gt,zla,z,Vo,tln$ea,X,vden,ppk,vegp$PAI,TT)
    cda<-lav$cp*lav$ph*(1-lav$vden)*(lav$z[3:(lav$m+2)]-lav$z[1:(lav$m)])/2*timestep
    ka<-lav$gt*c(lav$cp,cpair(previn$tabove))
    ka[1:lav$m]<-ifelse(ka[1:lav$m]>cda,cda,ka[1:lav$m])
    k<-c(rev(ka),cdk$k)
    cd<-c(rev(cda),cdk$cd)
    # Heat to add / loose
    X<-c(rev(lav$X),Xs,rep(0,(sm-1)))
    tc2<-c(previn$tabove,rev(lav$tc),previn$soiltc, previn$tsoil)
    tn2<-Thomas(tc2, tsoil, tcan, k, cd, n, X)
    tnair<-rev(tn2[1:(lav$m+2)])
    tnsoil<-tn2[(lav$m+2):(length(tn2)-1)]
    # vapour exchange
    zth2<-lav$z[2:(lav$m+1)]-lav$z[1:lav$m]
    Vmflux<-(lav$ea/pk)-lav$Vo
    Vmflux[Vmflux<0]<-0
    gtmx<-c(lav$ph/zth2*(2*timestep),Inf)
    gt2<-ifelse(lav$gt>gtmx,gtmx,lav$gt)
    tn2<-tnair[-1]; tn2<-tn2[-length(tn2)]
    if (length(lav$Vo)>1) {
      Vn<-ThomasV(lav$Vo,tn2,pk,theta,thetap,relhum,tcan,tnsoil[1],zth2,gt2,Vmflux,n,previn,soilp)
    } else {
      Vair<-(0.6108*exp(17.27*tcan/(tcan+237.3))*(relhum/100))/pk
      rhsoil<-soilrh(theta,soilp$b,soilp$psi_e,soilp$Smax,tnsoil[1])
      rhsoil[rhsoil>1]<-1
      Vsoil<-(0.6108*exp(17.27*tnsoil[1]/(tnsoil[1]+237.3))*rhsoil)/pk
      Vn<-c(Vsoil,lav$ea/pk,Vair)
    }
    # Interpolate
    if (length(lav$Vo) < m) {
      TX<-TT[length(TT)]+(pha/gt[m+1])*2
      T2<-c(0,lav$TT,TX)
      tn<-layerinterp(T2, TT, tnair)
      vn<-layerinterp(T2, TT, Vn)
    } else {
      tn <- tnair[2:(m+1)]
      vn <- Vn[2:(m+1)]
    }
    ea<-vn*pk
    # Canopy air layer in equilibrium with above canopy
  } else {
    X<-c(Xs,rep(0,(sm-1)))
    gt2<-1/cumsum(1/gt)[round(length(gha)/2,0)]
    ka<-gt2*cpair(previn$tabove)
    k<-c(ka,cdk$k)
    tc2<-c(previn$tabove,previn$soiltc,previn$tsoil)
    tn2<-Thomas(tc2, tsoil, tcan, k, cdk$cd, n, X)
    tnsoil<-tn2[2:(length(tn2)-1)]
    TX<-TT[length(TT)]+(pha/gt[m+1])*2
    wgt<-TT/TX
    tn<-(1-wgt)*tnsoil[1]+wgt*tcan
    eair<-(relhum/100)*0.6108*exp(17.27*tcan/(tcan+237.3))
    ea<-(1-wgt)*tln$esoil+wgt*eair
  }
  es<-0.6108*exp(17.27*tn/(tn+237.3))
  rh<-(ea/es)*100
  rh[rh>100]<-100
  rh[rh<0]<-0
  # Calculate Heat flux
  shade<-(1- exp(-Rabss$ref[m]*sum(vegp$PAI)))
  alb<-shade*Rabss$ref[m]+(1-shade)*vegp$refg
  Rlw<-5.67*10^-8*0.97*(0.5*tnsoil[1]+0.5*previn$soiltc[1]+273.15)^4
  # Latent heat
  L<-tln$L*vegp$PAI; L[L<-0]<-0
  Lt<-sum(L)
  tdif<-0.5*previn$soiltc[1]+0.5*tnsoil[1]-0.5*previn$tair-0.5*tair
  zdif<-abs(z-d+0.2*zm)
  sel<-which(zdif==min(zdif))[1]
  gt2<-gcanopy(uh,d+0.2*zm,0,tn[sel],tnsoil[1],hgt,sum(vegp$PAI),vegp$x,vegp$lw,
               vegp$cd,mean(vegp$iw),phi_m,pk)
  G<-gt2*cp[1]*(tn[sel]-tnsoil[1])
  tstG <- (1-alb)*Rsw-Rlw-Lt
  G <- ifelse(G > abs(tstG), abs(tstG), G)
  G <- ifelse(G < -abs(tstG), -abs(tstG), G)
  H<-(1-alb)*Rsw-Rlw-G-Lt
  # Incoming radiation
  Rswin<-Rabss$aRsw/(1-Rabss$ref)
  Rlwin<-canlw(tn, sum(vegp$PAI), 1-vegp$vegem, skyem = climvars$skyem)$lwin
  dataout<-list(tc=tn,soiltc=tnsoil,tleaf=tln$tleaf,tabove=tcan,uz=uz,z=z,sz=sz,zabove=zabove,
                rh=rh,relhum=relhum,tair=tair,tsoil=tsoil,pk=pk,Rabs=Rabs,
                gt=gt,gv=gv,gha=gha,H=H,L=L,G=G,Rswin=Rswin,Rlwin=Rlwin,psi_m=psi_m)
  return(dataout)
}
#' internal function to sort out vegetation parameters
.vegpsort <- function(vegp, i) {
  if (class(vegp$PAI) == "matrix") PAI <- as.vector(vegp$PAI[,i])
  if (class(vegp$pLAI) == "matrix") pLAI <- as.vector(vegp$pLAI[,i])
  if (length(vegp$zm0) > 1) zm0 <- vegp$zm0[i]
  vegp$PAI <- PAI
  vegp$pLAI <- pLAI
  vegp$zm0 <- zm0
  return(vegp)
}
#' Model spin-up for first time-step
#'
#' @description `spinup` runs the model repeatedly using data form the first time-step
#' for a set number of steps to ensure initial conditions are stable and
#' appropriate soil temperatures are set.
#'
#' @param climdata a data.frame of climate variables needed to run the run the model (dataset should follow format of [weather()])
#' @param vegp a list of vegetation parameters as returned by [microctools::habitatvars()]
#' @param soilp a list of soil parameters as returned by [soilinit()]
#' @param lat Latitude (decimal degrees)
#' @param long Longitude (decimal degrees, negative west of Greenwich meridion)
#' @param edgedist distance to open ground (m)
#' @param reqhgt optional height for which temperature is required (see details)
#' @param sdepth depth of deepest soil node (m)
#' @param zu height above ground of reference climate measurements (m)
#' @param theta volumetric water content of upper most soil layer in current time step (m^3 / m^3)
#' @param thetap volumetric water content of upper most soil layer in previous time step (m^3 / m^3)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param n forward / backward weighting for Thomas algorithm (see [Thomas()])
#' @param plotout optional logical indicating whether to a plot a profile of temperatures
#' upon completion.
#' @param steps number of iterations for which which to run spin-up
#' @param metopen optional logical indicating whether the wind measurement used as an input to
#' the model is from a nearby weather station located in open ground (TRUE) or above the canopy
#' for which temperatures are modelled (FALSE - see details)
#' @param windhgt height above ground of wind measurement. If `metopen` is FALSE, must be above
#' canopy.
#' @param zlafact numeric value indicating how close to leaves air temperatures are
#' needed for (1 - average leaf-air distance, 0.5 = half average leaf-air distance etc.).
#' Must be greater than 1.
#' @param surfwet fraction of vegetation surface acting like a free-water surface
#'
#' @return a list of model outputs as for [paraminit()] or [runonestep()]
#'
#' @details If `reqhgt` is set, and below the height of the canopy, the canopy node nearest
#' to that height is set at the value specified. The returned value `tabove` is then the
#' temperature at the top of the canopy. If `reqhgt` is above canopy, nodes are calculated
#' automatically, but `tabove` is the temperature at height `reqhgt`. If `reqhgt` is
#' negative, the soil node nearest to that height is set at the value specified.
#' @export
#' @examples
#' tme<-as.POSIXlt(weather$obs_time, format = "%Y-%m-%d %H:%M", tz = "UTC")
#' vegp <- habitatvars(4, 50, -5, tme, m = 20)
#' soilp<- soilinit("Loam")
#' # run spinup and produce profile plot
#' modelout<-spinup(weather, vegp, soilp, lat = 50, long = -5)
spinup <- function(climdata, vegp, soilp, lat, long, edgedist = 100, reqhgt = NA,
                   sdepth = 2, zu = 2, theta = 0.3, thetap = 0.3, merid = 0,
                   dst = 0, n = 0.6, plotout = TRUE, steps = 200, metopen = TRUE,
                   windhgt = 2, zlafact = 1, surfwet = 1) {
  tme<-as.POSIXlt(climdata$obs_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  timestep<-round(as.numeric(tme[2])-as.numeric(tme[1]),0)
  reqdepth <- NA
  if (is.na(reqhgt) == F) {
    if (reqhgt < 0) reqdepth <- reqhgt
  }
  tsoil<-mean(climdata$temp, na.rm = T)
  m <- length(vegp$thickw)
  sm <- length(soilp$z)
  previn <- paraminit(m, sm, vegp$hgt, climdata$temp[1], climdata$windspeed[1], climdata$relhum[1],
                      tsoil, climdata$swrad[1])
  dp <- climdata$difrad[1] / climdata$swrad[1]
  dp[!is.finite(dp)] <- 0.5
  climvars <- list(tair = climdata$temp[1], relhum = climdata$relhum[1], pk = climdata$pres[1],
                   u2 = climdata$windspeed[1], tsoil = tsoil, skyem = climdata$skyem[1],
                   Rsw = climdata$swrad[1], dp = dp)
  H<-0
  vegp2<-.vegpsort(vegp,1)
  for (i in 1:steps) {
    previn  <- runonestep(climvars, previn, vegp2, soilp, timestep, tme[1], lat,
                          long, edgedist, sdepth, reqhgt, zu, theta, thetap,
                          merid, dst, n, metopen, windhgt, zlafact, surfwet)
    H[i]<-previn$H
    previn$H<-mean(H)
  }
  if (plotout) plotresults(previn, vegp, climvars)
  previn
}
#' Runs the microclimate model over time
#'
#' @description `runmodel` is used to run the full model over time
#'
#' @param climdata a data.frame of climate variables (see e.g. `weather`)
#' @param vegp a list of vegetation parameters as returned by [microctools::habitatvars()]
#' @param soilp a list of soil parameters as returned by [soilinit()]
#' @param lat Latitude (decimal degrees)
#' @param long Longitude (decimal degrees, negative west of Greenwich meridion)
#' @param edgedist distance to open ground (m)
#' @param reqhgt optional height for which temperature is required (see details)
#' @param sdepth depth of deepest soil node (m)
#' @param zu height above ground of reference climate measurements (m)
#' @param theta volumetric water content of upper most soil layer in current time step (m^3 / m^3)
#' @param thetap volumetric water content of upper most soil layer in previous time step (m^3 / m^3)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param n forward / backward weighting for Thomas algorithm (see [Thomas()])
#' @param steps number of iterations over which to run [spinup()]
#' @param plotout optional logical indicating whether to a plot a profile of temperatures
#' upon completion.
#' @param plotsteps number of iterations run before resuls plotted if `plotout` set to TRUE
#' @param tsoil optional stable temperature of the deepest soil layer (see details)
#' @param metopen optional logical indicating whether the wind measurement used as an input to
#' the model is from a nearby weather station located in open ground (TRUE) or above the canopy
#' for which temperatures are modelled (FALSE - see details)
#' @param windhgt height above ground of wind measurement. If `metopen` is FALSE, must be above
#' canopy.
#' @param zlafact numeric value indicating how close to leaves air temperatures are
#' needed for (1 - average leaf-air distance, 0.5 = half average leaf-air distance etc.).
#' Must be greater than 1.
#' @param surfwet fraction of vegetation surface acting like free-water surface
#' @return a data.frame with the following elements:
#' @return `obs_time` POSIXlt object of times associated wiht eahc output
#' @return `reftemp` air temperature (deg C) at reference height - i.e. `climdata$temp`
#' @return `tout` air or soil temperature (deg C) (see details)
#' @return `tleaf` Leaf temperature (deg C) (see details)
#' @return `relhum` Relative humidity (Percentage) (see details)
#' @return `SWin` Incoming shortwave radiation (W / m^2) (see details)
#' @return `LWin` Incoming longwave radiation (W / m^2) (see details)
#' @return `H` Total sensible heat flux to/from canopy (W / m^2). Positive values indicate
#' sensible heat release from canopy to air.
#' @return `L` Sensible heat flux (W / m^2) from each canopy node or for the entire canopy (see details)
#' @return `G` Total heat flux to/from  ground (W / m^2). Positive values indicate
#' flux from canopy to ground. Negative values indicate flux from ground to canopy.
#' @import microctools
#' @export
#'
#' @details If `reqhgt` is set, and below the height of the canopy, the canopy node nearest
#' to that height is set at the value specified. The returned values `tout`, `tleaf`, `relhum`,
#' `L` and `Swin` and `Lwin` are then the values at that height. If `reqhgt` is above canopy, nodes are
#' calculated automatically, and `tout` and `relhum` are the temperature and relative humidity
#' at height `reqhgt`, `SWin` and `LWin` the values above canopy, `tleaf` the mean leaf
#' temperature and `L` the summed latent heat exchange for the entire canopy. If `reqhgt` is
#' negative, the soil node nearest to that height is set at the value specified, `tout` is
#' soil temperature at that node, `tleaf`, `relhum`, `Swin` and `Lwin` are the mean values for
#' the whole canopy. If `reqhgt` is not set, `tout`, `tleaf`, `relhum`, `SWin`
#' and `LWin` are mean values for the whole canopy. The parameter `tsoil` is the temperature
#' of soil below `sdepth`, which is assumed constant. If `tsoil` is not provided, it is
#' assigned a value equivelent to mean of `climdata$temp`.
#'
#' The wind profile of the canopy depends on the nature of the canopy itself, and often
#' available wind measurements are for a nearby weather station located in open ground where
#' it is possible that the height of the wind measurement is below the height of the canopy
#' being studied. When `metopen` is TRUE, the wind profile of reference grass surface is
#' used to derive estimates for two metres above the canopy of interest. When `metopen` is
#' FALSE, `windhgt` must be above canopy and the profile above the canopy being studied is used.
#'
#' @examples
#' tme<-as.POSIXlt(weather$obs_time, format = "%Y-%m-%d %H:%M", tz = "UTC")
#' vegp <- microctools::habitatvars(4, lat = 50, long = -5, tme, m = 20)
#' soilp<- soilinit("Loam")
#' dataout <- runmodel(weather, vegp, soilp, lat = 50, long = -5)
#' par(mfrow=c(2,1))
#' plot(tout~as.POSIXct(obs_time), data = dataout, type = "l", col = "red",
#'      xlab = "Month", ylab = "Temperature", ylim = c(-8.5, 27.5))
#' par(new=TRUE)
#' plot(reftemp~as.POSIXct(obs_time), data = dataout, type = "l", col = rgb(0,0,0,0.5),
#'      xlab = "", ylab = "Temperature", ylim = c(-8.5, 27.5), main = "Air temperature")
#' plot(tleaf~as.POSIXct(obs_time), data = dataout, type = "l", col = "darkgreen",
#'      xlab = "Month", ylab = "Temperature", ylim = c(-8.5, 27.5))
#' par(new=TRUE)
#' plot(reftemp~as.POSIXct(obs_time), data = dataout, type = "l", col = rgb(0,0,0,0.5),
#'     xlab = "", ylab = "", ylim = c(-8.5, 27.5), main = "Leaf temperature")
runmodel <- function(climdata, vegp, soilp, lat, long, edgedist = 100, reqhgt = NA,
                     sdepth = 2, zu = 2, theta = 0.3, thetap = 0.3, merid = 0,
                     dst = 0, n = 0.6, steps = 200, plotout = TRUE, plotsteps = 100,
                     tsoil = NA, metopen = TRUE, windhgt = 2, zlafact = 1, surfwet = 1) {
  tme<-as.POSIXlt(climdata$obs_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  previn <- spinup(climdata,vegp,soilp,lat,long,edgedist,reqhgt,sdepth,zu,theta,
                   thetap,merid,dst,n,plotout,steps,metopen,windhgt)
  timestep<-round(as.numeric(tme[2])-as.numeric(tme[1]),0)
  reqdepth <- NA
  if (is.na(reqhgt) == F) {
    if (reqhgt < 0) reqdepth <- reqhgt
  }
  if (is.na(tsoil)) tsoil<-mean(climdata$temp, na.rm=T)
  dp <- climdata$difrad / climdata$swrad
  dp[!is.finite(dp)] <- 0.5
  tout <- 0
  tleaf <- 0
  rh <- 0
  H <- 0
  L <- 0
  G <- 0
  Rswin <- 0
  Rlwin <- 0
  for (i in 1:length(tme)) {
    climvars<-list(tair = climdata$temp[i], relhum = climdata$relhum[i], pk = climdata$pres[i],
                   u2 = climdata$windspeed[i], tsoil = tsoil, skyem = climdata$skyem[i],
                   Rsw = climdata$swrad[i], dp = dp[i])
    vegp2<-.vegpsort(vegp, i)
    previn <- runonestep(climvars,previn,vegp2,soilp,timestep,tme[i],lat,long,
                         edgedist,sdepth,reqhgt,zu,theta,thetap,merid,dst,n,metopen,windhgt,
                         zlafact, surfwet)
    if (i%%plotsteps == 0 & plotout) plotresults(previn, vegp, climvars, i)
    if (is.na(reqhgt)) {
      tout[i] <- mean(previn$tc)
      tleaf[i] <- mean(previn$tleaf)
      rh[i] <- mean(previn$rh)
      Rswin[i] <- mean(previn$Rswin)
      Rlwin[i] <- mean(previn$Rlwin)
      L[i] <- sum(previn$L)
    } else {
      if (reqhgt >= vegp$hgt) {
        tout[i] <- previn$tabove
        tleaf[i] <- mean(previn$tleaf)
        es1 <- 0.6108*exp(17.27*climdata$temp[i]/(climdata$temp[i]+237.3))
        es2 <- 0.6108*exp(17.27*previn$tabove/(previn$tabove+237.3))
        ea<-(climdata$relhum/100)*es1
        relhum<-(ea/es2)*100
        relhum[relhum>100]<-100
        rh[i] <- relhum
        Rswin[i] <- climdata$swrad
        Rlwin[i] <- climdata$skyem*5.67*10^-8*(previn$tabove+273.15)^4
        L[i] <- sum(previn$L)
      }
      if (reqhgt >= 0 & reqhgt < vegp$hgt) {
        sel <- which(previn$z == reqhgt)
        tout[i] <- previn$tc[sel]
        tleaf[i] <- previn$tleaf[sel]
        rh[i] <- previn$rh[sel]
        Rswin[i] <- previn$Rswin[sel]
        Rlwin[i] <- previn$Rlwin[sel]
        L[i] <- previn$L[sel]
      }
      if (reqhgt < 0) {
        sel <- which(previn$sz == -reqhgt)
        tout[i] <- previn$soiltc[sel]
        tleaf[i] <- mean(previn$tleaf)
        rh[i] <- mean(previn$rh)
        Rswin[i] <- mean(previn$Rswin)
        Rlwin[i] <- mean(previn$Rlwin)
        L[i] <- sum(previn$L)
      }
    }
    H[i] <- previn$H
    G[i] <- previn$G
  }
  dataout<-data.frame(obs_time=tme,reftemp=climdata$temp,tout=tout,tleaf=tleaf,relhum=rh,
                      SWin=Rswin,LWin=Rlwin,H=H,L=L,G=G)
  return(dataout)
}

