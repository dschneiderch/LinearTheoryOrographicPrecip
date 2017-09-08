#' Run the Linear Theory model. THis is the function that calls the others.
#'
#' @export
#' @param h matrix of elevation
#' @param dx horizontal resolution
#' @param dy vertical resolution
#' @details returns a matrix the same size as h with precip values

LTmodel <- function(h,dx,dy){

# Parameters
Pinf=5 # backgroud precipitation rate (mm hr^-1)
windspeed=15  #m/s
winddir=270
T0=7 #ground temperature, deg C
# -5 is avg moist adiabatic lapse rate from wikipedia. 9.8 deg/km is dry adiabatic lapse rate. 6.49 deg/km is environmental lapse rate
gamma = -5.8e-3#-6.49 #-5.8 #-6.49 #environmental lapse rate deg/km
Gamma_m = -6.5e-3#-5
tauc=1000
tauf=1000
Nm=0.005
lat_coriolis = 0

# Approximate the effective moist static stability (Fraser et al. 1973)
# Smith and Barstad, 2004 - see appendix
# (Gamma_m/gamma) > 1     - atmosphere is statically stable, (Cloud formation in stable air is unlikely, https://en.wikipedia.org/wiki/Lapse_rate).
# Gamma_m = gamma         - moist neutral conditions
# Gamma_m = gamma - ((Nm^2*T0)/g); # average moist adiabatic lapse rate (degC km^-1) - (Eq A13)
# Nm=sqrt(abs(g/Tref*(gamma - Gamma_m)))
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


# Compute x and y components of windspeed
U <- -sin(winddir * 2 * pi / 360) * windspeed
V <- cos(winddir * 2 * pi / 360) * windspeed



# Pad DEM with 0 to remove edge effects from FFT ----
alist <- pad_dem(h)
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

#' Prepare precipitation transfer function relating the Fourier transforms of the terrain hhat(k,l) and the precipitiation field Phat(k,l) (Eq 49)

## Compute vertical wave number. Ignoring Coriolis effect, which might be important for scales >100 km (per cell? or mtn range?)
alist=vertical_wavenumber(ny,k,l,method='qgis') #<---- There does not appear to be a big difference if at all betwee kim and qgis  method.
m=alist$m
sigma=alist$sigma

#' Transfer function relating the Fourier transforms of the terrain hhat(k,l) and the precipitiation field Phat(k,l) (Eq 49)
Phat = (Cw*1i*sigma*hhat)/((1-1i*m*Hw)*(1+1i*sigma*tauc)*(1+1i*sigma*tauf));

#' Perform the inverse Fourier transform on Phat
Ppad = fft(Phat,inverse=TRUE)/length(Phat)

#' Add the background precipitation rate
Ppad = Pinf + as.numeric(Ppad)*3600#convert mm/sec to mm/hr

#' Apply the positive cutoff
Ppad[Ppad<0]=0

#' library(tidyverse)
#' tbl <- tibble(x=as.numeric(dem_list$x),y=as.numeric(dem_list$y),p=Ppad)
#' #' Return to matrix
Ppad=matrix(Ppad,ncol=nypad,byrow=F)

# image.matrix(Ppad)
#
# fields::image.plot(dem_list$xcoords,dem_list$ycoords,t(Ppad[nrow(Ppad):1,]))

# Extract original region
P=Ppad[(pad+1):(nx+pad),(pad+1):(ny+pad)]

return(P)

}
