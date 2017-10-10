#' Run the Linear Theory model. THis is the function that calls the others.
#'
#' @export
#' @param h matrix of elevation
#' @param dx horizontal resolution (m)
#' @param dy vertical resolution (m)
#' @param Pinf background precip (mm/hr)
#' @param tauc condensation time delay (s)
#' @param tauf fallout time delay (s)
#' @param windspeed windspeed (m/s)
#' @param winddir wind direction (m/s)
#' @param trunc should the result truncate negative values (default=TRUE)
#' @param method default = 'smith' (as the method is described in Smith and Barstad 2004); decides whether truncation should happen before or after adding background precip. If 'smith', truncation happens after adding background precip. If otherwise, truncation happens before adding background precip. For 'smith', you will have a lot of zeros if background precip is low whereas with the other method those grid cells will have the value of the background precip
#' @param trunc_value value below which modeled precip rates should go to zero. This will only make a difference if `method=LEM`
#' @details returns a matrix the same size as h with precip mm/hr

LTmodel <- function(h,dx,dy,Pinf,tauc,tauf,windspeed, winddir, trunc=TRUE, method='smith',trunc_value=0){

  # Parameters
  # Pinf=5 # backgroud precipitation rate (mm hr^-1)
  # windspeed=15  #m/s
  # winddir=270
  T0=7 #ground temperature, deg C
  # -5 is avg moist adiabatic lapse rate from wikipedia. 9.8 deg/km is dry adiabatic lapse rate. 6.49 deg/km is environmental lapse rate
  gamma = -6.5e-3# -5.8e-3#-6.49 #-5.8 #-6.49 #environmental lapse rate deg/km
  Gamma_m = -5e-3#-6.5e-3#-5
  # tauc=1000
  # tauf=1000
  # Nm = 0.005
  lat_coriolis = 0

  # Approximate the effective moist static stability (Fraser et al. 1973)
  # Smith and Barstad, 2004 - see appendix
  # (Gamma_m/gamma) > 1     - atmosphere is statically stable, (Cloud formation in stable air is unlikely, https://en.wikipedia.org/wiki/Lapse_rate).
  # Gamma_m = gamma         - moist neutral conditions
  # Gamma_m = gamma - ((Nm^2*T0)/g); # average moist adiabatic lapse rate (degC km^-1) - (Eq A13)

  #### It seems strange to me that Gamma_m is always largr than gamma if you specify Nm 0-0.01 as suggested in Smith and Barstad. the Internet suggests that Gamma_m can be between gamma and 9.8 deg/km also, or even >9.8 in some cases


  # Constants
  L = 2.5e6 # latent heat (J kg^-1)
  Rv = 461; #gas constant for vapor (J kg^-1 K^-1)
  Tref = T0 + 273.15 # ground temperature (K)
  g = 9.81; # gravitational acceleration (m s^-2)
  rhosref = 7.4e-3 # saturation water vapor density (kg m^-3)
  eps=.Machine$double.eps

  # Computed Constants
  Hw = -(Rv*Tref^2)/(L*(gamma)) # water vapor scale height (m) - (Eq A4)
  #esref = rhosref*Rv*Tref; % ground saturation vapor pressure (Pa) - (Eq A1)
  Cw = rhosref*Gamma_m/gamma # uplift sensitivity factor (kg m^-3) - (Eq A9)
  Nm=sqrt(abs(g/Tref*(gamma - Gamma_m)))

  # Compute x and y components of windspeed
  U <- -sin(winddir * 2 * pi / 360) * windspeed
  V <- cos(winddir * 2 * pi / 360) * windspeed


  # Pad DEM with 0 to remove edge effects from FFT ----
  alist <- pad_dem(h,0)
  hpad=alist$dem
  pad=alist$pad

  nxpad=dim(hpad)[1]
  nypad=dim(hpad)[2]

  # image.matrix(hpad)

  #' Prepare Fourier transform of the topography ----
  # alist <- fftfreq_kim(dx,dy,nxpad, nypad) <------ DON'T USE. doesn't seem to work when  compared with gaussian hill dem
  alist <- fftfreq_qgis(dx, dy, nxpad, nypad)
  k=alist$k
  l=alist$l

  #' Do 2D FFT of Topography ----
  hhat=fft(hpad)
  # print(hhat[1:2,1:2])

  #' Prepare precipitation transfer function relating the Fourier transforms of the terrain hhat(k,l) and the precipitiation field Phat(k,l) (Eq 49)

  ## Compute vertical wave number. Ignoring Coriolis effect, which might be important for scales >100 km (per cell? or mtn range?)
  alist=vertical_wavenumber(dim(h),Nm,U,V,k,l,eps,method='qgis') #<---- There does not appear to be a big difference if at all betwee kim and qgis  method.
  m=alist$m
  sigma=alist$sigma

  #' Transfer function relating the Fourier transforms of the terrain hhat(k,l) and the precipitiation field Phat(k,l) (Eq 49)
  Phat = (Cw*1i*sigma*hhat)/((1-1i*m*Hw)*(1+1i*sigma*tauc)*(1+1i*sigma*tauf));

  #' Perform the inverse Fourier transform on Phat
  Ppad = fft(Phat,inverse=TRUE)/length(Phat)
  Ppad = as.numeric(Ppad) #remove imaginary values

  #' Should background precip be added before or after truncation (if truncation happens)?
  if(method=='smith'){
    #' Add the background precipitation rate
    Ppad = Pinf + as.numeric(Ppad)*3600#convert mm/sec to mm/hr

    #' Apply the positive cutoff
    if(trunc){
      Ppad[Ppad < 0] = 0
    }

  } else if(method=='LEM') {
    #' Apply the positive cutoff
    if(trunc){
      Ppad[Ppad < trunc_value] = 0
    }
    #' Add the background precipitation rate
    Ppad = Pinf + as.numeric(Ppad) * 3600#convert mm/sec to mm/hr
    Ppad[Ppad < 0] = 0 #if background precip is less than trunc_value you need to zero remaining negative values
  }

  #' library(tidyverse)
  #' tbl <- tibble(x=as.numeric(dem_list$x),y=as.numeric(dem_list$y),p=Ppad)
  #' Return to matrix
  Ppad=matrix(Ppad,ncol=nypad,byrow=F)

  # image.matrix(Ppad)
  #
  # fields::image.plot(dem_list$xcoords,dem_list$ycoords,t(Ppad[nrow(Ppad):1,]))

  # Extract original region
  P=Ppad[(pad+1):(nrow(h)+pad),(pad+1):(ncol(h)+pad)]

  return(P)

}
