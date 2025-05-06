---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(one_ebm)=

# 14. One-dimensional Energy Balance Model

```{note}
The Python scripts used below and some materials are modified from Prof. Brian E. J. Rose's climlab [website](https://brian-rose.github.io/ClimateLaboratoryBook/courseware/heat-transport.html).
```

## Spatial structure of insolation and surface temperature

```{code-cell} ipython3
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import climlab
from climlab import constants as const

#  Calculate daily average insolation as function of latitude and time of year
lat = np.linspace( -90., 90., 500 )
days = np.linspace(0, const.days_per_year, 365 )
Q = climlab.solar.insolation.daily_insolation( lat, days )

##  daily surface temperature from  NCEP reanalysis
#ncep_url = "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/"
ncep_url = "/Users/yuchiaol_ntuas/Desktop/ebooks/data/"
ncep_temp = xr.open_dataset( ncep_url + "skt.sfc.day.1981-2010.ltm.nc", decode_times=False)
#url = 'http://apdrc.soest.hawaii.edu:80/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/'
#skt_path = 'surface_gauss/skt'
#ncep_temp = xr.open_dataset(url+skt_path)
ncep_temp_zon = ncep_temp.skt.mean(dim='lon')

fig = plt.figure(figsize=(12,6))

ax1 = fig.add_subplot(121)
CS = ax1.contour( days, lat, Q , levels = np.arange(0., 600., 50.) )
ax1.clabel(CS, CS.levels, inline=True, fmt='%1.0f', fontsize=10)
ax1.set_title('Daily average insolation', fontsize=18 )
ax1.contourf ( days, lat, Q, levels=[-100., 0.], colors='k' )

ax2 = fig.add_subplot(122)
CS = ax2.contour( (ncep_temp.time - ncep_temp.time[0])/const.hours_per_day, ncep_temp.lat, 
                 ncep_temp_zon.T, levels=np.arange(210., 310., 10. ) )
ax2.clabel(CS, CS.levels, inline=True, fmt='%1.0f', fontsize=10)
ax2.set_title('Observed zonal average surface temperature', fontsize=18 )

for ax in [ax1,ax2]:
    ax.set_xlabel('Days since January 1', fontsize=16 )
    ax.set_ylabel('Latitude', fontsize=16 )
    ax.set_yticks([-90,-60,-30,0,30,60,90])
    ax.grid()
```

```{note}
Dr Zack Labe has nice animations to demonstration solar radiation: [here](https://zacklabe.com/arctic-climate-seasonality-and-variability/)
```

## A RCE model with latitudinal structure

```{code-cell} ipython3
# A two-dimensional domain
state = climlab.column_state(num_lev=30, num_lat=40, water_depth=10.)
#  Specified relative humidity distribution
h2o = climlab.radiation.ManabeWaterVapor(name='Fixed Relative Humidity', state=state)
#  Hard convective adjustment
conv = climlab.convection.ConvectiveAdjustment(name='Convective Adjustment', state=state, adj_lapse_rate=6.5)
#  Daily insolation as a function of latitude and time of year
sun = climlab.radiation.DailyInsolation(name='Insolation', domains=state['Ts'].domain)
#  Couple the radiation to insolation and water vapor processes
rad = climlab.radiation.RRTMG(name='Radiation',
                             state=state, 
                             specific_humidity=h2o.q, 
                             albedo=0.125,
                             coszen=sun.coszen,
                             irradiance_factor=sun.insolation)
model = climlab.couple([rad,sun,h2o,conv], name='RCM')
print( model)

model.compute_diagnostics()

fig, ax = plt.subplots()
ax.plot(model.lat, model.insolation)
ax.set_xlabel('Latitude')
ax.set_ylabel('Insolation (W/m2)');

#climlab.radiation.DailyInsolation
#climlab.radiation.AnnualMeanInsolation

#  model is initialized on Jan. 1

#  integrate forward just under 1/4 year... should get about to the NH spring equinox
model.integrate_days(31+28+22)
Q_spring = model.insolation.copy()
#  Then forward to NH summer solstice
model.integrate_days(31+30+31)
Q_summer = model.insolation.copy()
#  and on to autumnal equinox
model.integrate_days(30+31+33)
Q_fall = model.insolation.copy()
#  and finally to NH winter solstice
model.integrate_days(30+31+30)
Q_winter = model.insolation.copy()

fig, ax = plt.subplots()
ax.plot(model.lat, Q_spring, label='Spring')
ax.plot(model.lat, Q_summer, label='Summer')
ax.plot(model.lat, Q_fall, label='Fall')
ax.plot(model.lat, Q_winter, label='Winter')
ax.legend()
ax.set_xlabel('Latitude')
ax.set_ylabel('Insolation (W/m2)');

```

```{code-cell} ipython3
# time integration

model.integrate_years(4.)
#model.integrate_years(1.)

model.timeave.keys()

fig, ax = plt.subplots()
ax.plot(model.lat, model.timeave['insolation'])
ax.set_xlabel('Latitude')
ax.set_ylabel('Annual-mean Insolation (W/m2)')

# Plot annual mean surface temperature in the model,
#   compare to observed annual mean surface temperatures
fig, ax = plt.subplots()
ax.plot(model.lat, model.timeave['Ts'], label='RCE')
ax.plot(ncep_temp_zon.lat, ncep_temp_zon.mean(dim='time'), label='obs')
ax.set_xticks(range(-90,100,30))
ax.grid(); ax.legend();

```

- The tropical regions warm too much in RCE.
- The polar regions cool too much in RCE.

What process is lacking here?

We now turn to look at the vertical structure.

```{code-cell} ipython3
#  Observed air temperature from NCEP reanalysis
## The NOAA ESRL server is shutdown! January 2019
## may shutdown again in 2025
ncep_url = "/Users/yuchiaol_ntuas/Desktop/ebooks/data/"
ncep_air = xr.open_dataset( ncep_url + "air.mon.1981-2010.ltm.nc", decode_times=False)
#air = xr.open_dataset(url+'pressure/air')
#ncep_air = air.rename({'lev':'level'})
level_ncep_air = ncep_air.level
lat_ncep_air = ncep_air.lat
Tzon = ncep_air.air.mean(dim=('time','lon'))

#  Compare temperature profiles in RCE and observations
contours = np.arange(180., 350., 15.)

fig = plt.figure(figsize=(14,6))
ax1 = fig.add_subplot(1,2,1)
cax1 = ax1.contourf(lat_ncep_air, level_ncep_air, Tzon+const.tempCtoK, levels=contours)
fig.colorbar(cax1)
ax1.set_title('Observered temperature (K)')

ax2 = fig.add_subplot(1,2,2)
field = model.timeave['Tatm'].transpose()
cax2 = ax2.contourf(model.lat, model.lev, field, levels=contours)
fig.colorbar(cax2)
ax2.set_title('RCE temperature (K)')

for ax in [ax1, ax2]:
    ax.invert_yaxis()
    ax.set_xlim(-90,90)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])

```

The TOA energy budget is nicely closed:

```{code-cell} ipython3
fig, ax = plt.subplots()
ax.plot(model.lat, model.timeave['ASR'], label='ASR')
ax.plot(model.lat, model.timeave['OLR'], label='OLR')
ax.set_xlabel('Latitude')
ax.set_ylabel('W/m2')
ax.legend(); ax.grid()
```

We now look at the reanalysis data:
```{code-cell} ipython
# Get TOA radiative flux data from NCEP reanalysis
# downwelling SW
ncep_url = "/Users/yuchiaol_ntuas/Desktop/ebooks/data/"
dswrf = xr.open_dataset(ncep_url + 'dswrf.ntat.mon.1981-2010.ltm.nc', decode_times=False)
#dswrf = xr.open_dataset(url + 'other_gauss/dswrf')
#  upwelling SW
uswrf = xr.open_dataset(ncep_url + 'uswrf.ntat.mon.1981-2010.ltm.nc', decode_times=False)
#uswrf = xr.open_dataset(url + 'other_gauss/uswrf')
#  upwelling LW
ulwrf = xr.open_dataset(ncep_url + 'ulwrf.ntat.mon.1981-2010.ltm.nc', decode_times=False)
#ulwrf = xr.open_dataset(url + 'other_gauss/ulwrf')

ASR = dswrf.dswrf - uswrf.uswrf
OLR = ulwrf.ulwrf

ASRzon = ASR.mean(dim=('time','lon'))
OLRzon = OLR.mean(dim=('time','lon'))

ticks = [-90, -60, -30, 0, 30, 60, 90]
fig, ax = plt.subplots()
ax.plot(ASRzon.lat, ASRzon, label='ASR')
ax.plot(OLRzon.lat, OLRzon, label='OLR')
ax.set_ylabel('W m$^{-2}$')
ax.set_xlabel('Latitude')
ax.set_xlim(-90,90); ax.set_ylim(50,310)
ax.set_xticks(ticks);
ax.set_title('Observed annual mean radiation at TOA')
ax.legend(); ax.grid();

```

We quickly look at the CESM simulations:

```{code-cell} ipython3
casenames = {'cpl_control': 'cpl_1850_f19',
             'cpl_CO2ramp': 'cpl_CO2ramp_f19',
             'som_control': 'som_1850_f19',
             'som_2xCO2':   'som_1850_2xCO2',
            }
# The path to the THREDDS server, should work from anywhere
#basepath = 'http://thredds.atmos.albany.edu:8080/thredds/dodsC/CESMA/'
#basepath = "/Users/yuchiaol_ntuas/Desktop/ebooks/data/"
# For better performance if you can access the roselab_rit filesystem (e.g. from JupyterHub)
#basepath = '/roselab_rit/cesm_archive/'
casepaths = {}
for name in casenames:
#    casepaths[name] = basepath + casenames[name] + '/concatenated/'
    casepaths[name] = "/Users/yuchiaol_ntuas/Desktop/ebooks/data/"

# make a dictionary of all the CAM atmosphere output
atm = {}
for name in casenames:
    path = casepaths[name] + casenames[name] + '.cam.h0.nc'
    print('Attempting to open the dataset ', path)
    atm[name] = xr.open_dataset(path)

lat_cesm = atm['cpl_control'].lat
ASR_cesm = atm['cpl_control'].FSNT
OLR_cesm = atm['cpl_control'].FLNT

# extract the last 10 years from the slab ocean control simulation
# and the last 20 years from the coupled control
nyears_slab = 10
nyears_cpl = 20
clim_slice_slab = slice(-(nyears_slab*12),None)
clim_slice_cpl = slice(-(nyears_cpl*12),None)

#  For now we're just working with the coupled control simulation
#  Take the time and zonal average
ASR_cesm_zon = ASR_cesm.isel(time=clim_slice_slab).mean(dim=('lon','time'))
OLR_cesm_zon = OLR_cesm.isel(time=clim_slice_slab).mean(dim=('lon','time'))

fig, ax = plt.subplots()
ax.plot(lat_cesm, ASR_cesm_zon, label='ASR')
ax.plot(lat_cesm, OLR_cesm_zon, label='OLR')
ax.set_ylabel('W m$^{-2}$')
ax.set_xlabel('Latitude')
ax.set_xlim(-90,90); ax.set_ylim(50,310)
ax.set_xticks([-90, -60, -30, 0, 30, 60, 90]);
ax.set_title('CESM control simulation: Annual mean radiation at TOA')
ax.legend(); ax.grid();

```

- Note the energy balance latitude!!!

## Zonal mean energy budget and heat transport
For a regional energy balance model, we may need to consider the energy transport.

```{figure} /_static/lecture_specific/lecture1_figures/heat_transport_tmp1.png
:scale: 50%
```

```{math}
:label: my_label122
\begin{eqnarray}
\frac{\partial E}{\partial t} &=& ASR - OLR + \frac{1}{2\pi a^{2}\delta \phi \cos(\phi)}(H(\phi)-H(\phi + \delta \phi)) \\
 &=& ASR - OLR - \frac{1}{2\pi a^{2}\cos(\phi)}(\frac{H(\phi + \delta \phi)-H(\phi)}{\delta \phi})
\end{eqnarray}
```

if we take $\delta \phi \to 0$, we can get:
```{math}
:label: my_label123
\frac{\partial E}{\partial t} = ASR - OLR - \frac{1}{2\pi a^{2}\cos(\phi)}\frac{\partial H(\phi)}{\partial \phi}
```

We can define the dynamical heating rate (the convergence of energy transport):
```{math}
:label: my_label124
h = - \frac{1}{2\pi a^{2}\cos(\phi)}\frac{\partial H(\phi)}{\partial \phi}
```

If the climate system reaches an equilibrium state, that is $\frac{\partial E}{\partial t}=0$, then the divergence of heat transport must balance the net radiative heating/cooling at every latitude.

```{math}
:label: my_label125
R_{TOA} = ASR - OLR = \frac{1}{2\pi a^{2}\cos(\phi)}\frac{\partial H(\phi)}{\partial \phi}
```

```{math}
:label: my_label126
\frac{\partial H(\phi)}{\partial \phi} = 2\pi a^{2}\cos(\phi)R_{TOA}
```

We can take integral along latitude from the South Pole:
```{math}
:label: my_label127
\begin{eqnarray}
\int_{-\pi /2}^{\phi} \frac{\partial H(\phi)}{\partial \phi'} d\phi' = 2\pi a^{2}\int_{-\pi /2}^{\phi}\cos(\phi')R_{TOA}d\phi' \\
\rightarrow H(\phi) - H(-\pi /2) = 2\pi a^{2}\int_{-\pi /2}^{\phi}\cos(\phi')R_{TOA}d\phi' \\
\rightarrow H(\phi) = 2\pi a^{2}\int_{-\pi /2}^{\phi}\cos(\phi')R_{TOA}d\phi'
\end{eqnarray}
```

We assume $H(-\pi /2) = 0$ (make sense?) and integral all the way to the North Pole:
```{math}
:label: my_label128
H(\pi /2) = 2\pi a^{2}\int_{-\pi /2}^{\pi /2}\cos(\phi')R_{TOA}d\phi' = 0 (?)
```

If $H(\pi /2) \ne 0$, what does this mean?

Let's look at the poleward heat transport in reanalysis data and CESM1 simulations.

```{code-cell} ipython3
def inferred_heat_transport(energy_in, lat=None, latax=None):
    '''Compute heat transport as integral of local energy imbalance.
    Required input:
        energy_in: energy imbalance in W/m2, positive in to domain
    As either numpy array or xarray.DataArray
    If using plain numpy, need to supply these arguments:
        lat: latitude in degrees
        latax: axis number corresponding to latitude in the data
            (axis over which to integrate)
    returns the heat transport in PW.
    Will attempt to return data in xarray.DataArray if possible.
    '''
    from scipy import integrate
    from climlab import constants as const
    if lat is None:
        try: lat = energy_in.lat
        except:
            raise InputError('Need to supply latitude array if input data is not self-describing.')
    lat_rad = np.deg2rad(lat)
    coslat = np.cos(lat_rad)
    field = coslat*energy_in
    if latax is None:
        try: latax = field.get_axis_num('lat')
        except:
            raise ValueError('Need to supply axis number for integral over latitude.')
    #  result as plain numpy array
    integral = integrate.cumulative_trapezoid(field, x=lat_rad, initial=0., axis=latax)
    result = (1E-15 * 2 * np.pi * const.a**2 * integral)
    if isinstance(field, xr.DataArray):
        result_xarray = field.copy()
        result_xarray.values = result
        return result_xarray
    else:
        return result

fig, ax = plt.subplots()
ax.plot(lat_cesm, inferred_heat_transport(ASR_cesm_zon - OLR_cesm_zon))
ax.set_ylabel('PW')
ax.set_xlabel('Latitude')
ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
ax.grid()
ax.set_title('Total northward heat transport inferred from CESM control simulation')
```

```{code-cell} ipython3
#  Need to flip the arrays because we want to start from the south pole
Rtoa_ncep = ASRzon-OLRzon
lat_ncep = ASRzon.lat
fig, ax = plt.subplots()
ax.plot(lat_ncep, inferred_heat_transport(Rtoa_ncep))
ax.set_ylabel('PW')
ax.set_xlabel('Latitude')
ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
ax.grid()
ax.set_title('Total northward heat transport inferred from NCEP reanalysis')

```

We need to make a correction to adjust.

```{code-cell} ipython3
#  global average of TOA radiation in reanalysis data
weight_ncep = np.cos(np.deg2rad(lat_ncep)) / np.cos(np.deg2rad(lat_ncep)).mean(dim='lat')
imbal_ncep = (Rtoa_ncep).weighted(weight_ncep).mean(dim='lat')
print( 'The net downward TOA radiation flux in NCEP renalysis data is %0.1f W/m2.' %imbal_ncep)

Rtoa_ncep_balanced = Rtoa_ncep - imbal_ncep
newimbalance = float(Rtoa_ncep_balanced.weighted(weight_ncep).mean(dim='lat'))
print( 'The net downward TOA radiation flux after balancing the data is %0.2e W/m2.' %newimbalance)

fig, ax = plt.subplots()
ax.plot(lat_ncep, inferred_heat_transport(Rtoa_ncep_balanced))
ax.set_ylabel('PW')
ax.set_xlabel('Latitude')
ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
ax.grid()
ax.set_title('Total northward heat transport inferred from NCEP reanalysis (after global balancing)')
```

Does this correction make sense?



```{figure} /_static/lecture_specific/lecture1_figures/transport_tmp1.jpg
---
scale: 90%
---
The ERBE period zonal mean annual cycle of the meridional energy transport in PW by (a) the atmosphere and ocean as inferred from ERBE RT, NRA δAE/δt, and GODAS δOE/δt; (b) the atmosphere based on NRA; and (c) by the ocean as implied by ERBE + NRA FS and GODAS δOE/δt. Stippling and hatching in (a)–(c) represent regions and times of year in which the standard deviation of the monthly mean values among estimates, some of which include the CERES period (see text), exceeds 0.5 and 1.0 PW, respectively. (d) The median annual mean transport by latitude for the total (gray), atmosphere (red), and ocean (blue) accompanied with the associated ±2σ range (shaded). Source: [Fasullo and Trenberth (2008)](https://journals.ametsoc.org/view/journals/clim/21/10/2007jcli1936.1.xml)

```

## Separate energy budget for the atmosphere and ocean
We can separate the energy budget for atmospheric and oceanic components. We need surface heat fluxes (SHF) between the interface.

SHF usually consistes of four components:
- shortwave radiation
- longwave radiation
- sensible heat flux
- latent heat flux

What we might be missing?

```{code-cell} ipython3
# monthly climatologies for surface flux data from reanalysis
#  all defined as positive UP
ncep_url = "/Users/yuchiaol_ntuas/Desktop/ebooks/data/"
ncep_nswrs = xr.open_dataset(ncep_url + "nswrs.sfc.mon.1981-2010.ltm.nc", decode_times=False)
ncep_nlwrs = xr.open_dataset(ncep_url + "nlwrs.sfc.mon.1981-2010.ltm.nc", decode_times=False)
ncep_shtfl = xr.open_dataset(ncep_url + "shtfl.sfc.mon.1981-2010.ltm.nc", decode_times=False)
ncep_lhtfl = xr.open_dataset(ncep_url + "lhtfl.sfc.mon.1981-2010.ltm.nc", decode_times=False)

#  Calculate ANNUAL AVERAGE net upward surface flux
ncep_net_surface_up = (ncep_nlwrs.nlwrs
                     + ncep_nswrs.nswrs
                     + ncep_shtfl.shtfl
                     + ncep_lhtfl.lhtfl
                      ).mean(dim='time')

lon_ncep = ncep_net_surface_up.lon
fig, ax = plt.subplots()
cax = ax.pcolormesh(lon_ncep, lat_ncep, ncep_net_surface_up, 
               cmap=plt.cm.seismic, vmin=-200., vmax=200. )
fig.colorbar(cax, ax=ax)
ax.set_title('Net upward surface energy flux in NCEP Reanalysis data')

```

```{note}
It is always important to check the sign of these fluxes. Different dataset may have difference definition of the sign.
```

We can now separate the whole system to upper atmosphere and lower ocean with the SHF exchanges energy between them. If we define the sign of SHF as positive upward, that is from the ocean to the atmosphere.

```{figure} /_static/lecture_specific/lecture1_figures/colum_atm_ocn_tmp1.png
---
scale: 50%
---
```

For the ocean:
```{math}
:label: my_label129
\frac{\partial E_o}{\partial t} = -SHF + OHT
 = -SHF - \frac{1}{2\pi a^{2}\cos(\phi)}\frac{\partial H_o}{\partial \phi}
```

If we consider an equilibrium state:
```{math}
:label: my_label130
H_{o}(\phi) = -2\pi a^{2} \int_{-\pi/2}^{\phi} SHF \cos(\phi')d\phi'
```

For the atmosphere:
```{math}
:label: my_label131
\frac{\partial E_a}{\partial t} = R_{TOA} -SHF + AHT
 = R_{TOA} -SHF - \frac{1}{2\pi a^{2}\cos(\phi)}\frac{\partial H_a}{\partial \phi}
```

Similarly if we consider an equilibrium state:
```{math}
:label: my_label132
H_{a}(\phi) = 2\pi a^{2} \int_{-\pi/2}^{\phi} (R_{TOA} + SHF) \cos(\phi')d\phi'
```

We can further separate the AHT into wet and dry component, that is $AHT = AHT_{dry} + AHT_{wet}$.
$AHT_{wet}$ is also called latent heat transport. Considering the water budget in the atmosphere:
```{math}
:label: my_label133
L_{v}\frac{\partial Q}{\partial t} = L_{v}(E-P) - AHT_{wet} = L_{v}(E-P) - \frac{1}{2\pi a^{2}\cos(\phi)}\frac{\partial H_{wet,LH}}{\partial \phi}

```
where $Q$ is the depth-integrated water vapor in kg/m$^2$, $E$ and $P$ are evaporation and precipitation in kg/m$^2$/s (=mm/s, why?), $L_{v}=2.5\times 10^{-6}$ is the latend heat of vaporization.

```{math}
:label: my_label134
H_{wet, LH}(\phi) = 2\pi a^{2} \int_{-\pi/2}^{\phi} L_{v}(E-P) \cos(\phi')d\phi'
```

So the dry component can be calculated as:
```{math}
:label: my_label135
H_{dry} = H_{a} - H_{wet, LH}
```

Let's look at CESM1 results:


```{code-cell} ipython3
def CESM_heat_transport(run, timeslice=clim_slice_cpl):
    #  Take zonal and time averages of the necessary input fields
    fieldlist = ['FLNT','FSNT','LHFLX','SHFLX','FLNS','FSNS','PRECSC','PRECSL','QFLX','PRECC','PRECL']
    zon = run[fieldlist].isel(time=timeslice).mean(dim=('lon','time'))
    OLR = zon.FLNT
    ASR = zon.FSNT
    Rtoa = ASR - OLR  # net downwelling radiation
    #  surface energy budget terms, all defined as POSITIVE UP
    #    (from ocean to atmosphere)
    LHF = zon.LHFLX
    SHF = zon.SHFLX
    LWsfc = zon.FLNS
    SWsfc = -zon.FSNS
    SnowFlux =  ((zon.PRECSC + zon.PRECSL) *
                      const.rho_w * const.Lhfus)
    # net upward radiation from surface
    SurfaceRadiation = LWsfc + SWsfc
    # net upward surface heat flux
    SurfaceHeatFlux = SurfaceRadiation + LHF + SHF + SnowFlux
    # net heat flux into atmosphere
    Fatmin = Rtoa + SurfaceHeatFlux
    #  hydrological cycle, all terms in  kg/m2/s or mm/s
    Evap = zon.QFLX
    Precip = (zon.PRECC + zon.PRECL) * const.rho_w
    EminusP = Evap - Precip
        
    # heat transport terms
    HT = {}
    HT['total'] = inferred_heat_transport(Rtoa)
    HT['atm'] = inferred_heat_transport(Fatmin)
    HT['ocean'] = inferred_heat_transport(-SurfaceHeatFlux)
    HT['latent'] = inferred_heat_transport(EminusP*const.Lhvap) # atm. latent heat transport from moisture imbal.
    HT['dse'] = HT['atm'] - HT['latent']  # dry static energy transport as residual

    return HT

#  Compute heat transport partition for both control and 2xCO2 simulations
HT_control = CESM_heat_transport(atm['cpl_control'])
HT_2xCO2 = CESM_heat_transport(atm['cpl_CO2ramp'])

fig = plt.figure(figsize=(16,6))
runs = [HT_control, HT_2xCO2]
N = len(runs)

for n, HT in enumerate([HT_control, HT_2xCO2]):
    ax = fig.add_subplot(1, N, n+1)
    ax.plot(lat_cesm, HT['total'], 'k-', label='total', linewidth=2)
    ax.plot(lat_cesm, HT['atm'], 'r-', label='atm', linewidth=2)
    ax.plot(lat_cesm, HT['dse'], 'r--', label='dry')
    ax.plot(lat_cesm, HT['latent'], 'r:', label='latent')
    ax.plot(lat_cesm, HT['ocean'], 'b-', label='ocean', linewidth=2)

    ax.set_xlim(-90,90)
    ax.set_ylim(-5,6)
    ax.set_xticks(ticks)
    ax.legend(loc='upper left')
    ax.grid()
    ax.set_xlabel('Latitude')
    ax.set_ylabel('PW')

```

A couple of notes below:

- the energy transport can be also calculated as (for ocean?):
```{math}
:label: my_label136
H = \int_{0}^{2\pi}\int_{bottom}^{top} \rho \cdot v \cdot e \cdot a \cos(\phi) dzd\lambda
```
- for the ocean, $e_{o} \approx c_{w}T$, where $c_{w} = 4.2\times 10^{3}$ J/kg/K is the specific heat of seawater.
- for the atmosphere, we usually use moist static energy (MSE):
```{math}
:label: my_label137
e_{a} \approx MSE = c_{p}T + gz + L_{v}q
```
- We can decomposite MSE into dry and wet/latent heat component. How?
- We usually neglect the kinetic energy $e_{k} = \frac{|\vec{v}|^2}{2}$, why?

Professor Yen-Ting Hwang's GRL paper:
```{figure} /_static/lecture_specific/lecture1_figures/hwang_2011_grl_tmp1.png
---
scale: 140%
---
Changes in northward energy transports in PW from 2001∼2020 to 2081 ∼2100 in the A2 Scenario: (a) atmospheric energy transport, (b) moisture (solid) and DSE (dashed) transport, and (c) oceanic energy transport. Source: [Hwang et al. (2011)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011GL048546)
```

## Parameterization for heat transport
If the climate system is not in equilibrium, how do you calculate the AHT and eventually the temperature? I think this is a very difficult problem, rooted in the turbulent nonclosure (a lot to discuss...). One way to do it is that we parameterize the processes involve, or make it "simpler".

We start with the energy budget for a thin zonal band at latitude $\phi$:

```{math}
:label: my_label138
\frac{\partial E(\phi)}{\partial t} = ASR(\phi) - OLR(\phi) - \frac{1}{2\pi a^{2}\cos(\phi)}\frac{\partial H}{\partial \phi}
``` 

We approximate or "parameterize" the heat transport as a down-gradient diffusion process:
```{math}
:label: my_label139
H(\phi) = -2\pi a^{2}\cos(\phi)D\frac{\partial T_{s}}{\partial \phi}
```
where $D$ is the diffusivity or thermal conductivity in unit of W/m$^2$/K.
It is noted that we consider surface air temperature in the gradient. Why can we do this?

Let's put paramertized $H$ into the budget equation:
```{math}
:label: my_label140
\frac{\partial E(\phi)}{\partial t} = ASR(\phi) - OLR(\phi) - \frac{1}{2\pi a^{2}\cos(\phi)}\frac{\partial }{\partial \phi}(-2\pi a^{2}\cos(\phi)D\frac{\partial T_{s}}{\partial \phi})
```

Assume $D$ is a constant:
```{math}
:label: my_label141
\frac{\partial E(\phi)}{\partial t} = ASR(\phi) - OLR(\phi) + \frac{D}{\cos(\phi)}\frac{\partial }{\partial \phi}(\cos(\phi)\frac{\partial T_{s}}{\partial \phi})
```

Next we parameterize the internal energy as:
```{math}
:label: my_label142
E(\phi) = C(\phi)T_{s}(\phi)
```
where $C$ is the effective heat capacity in units of J/m$^2$/K.

So the energy buget equation becomes:
```{math}
:label: my_label143
C(\phi)\frac{\partial T_{s}}{\partial t} = ASR(\phi) - OLR(\phi) + \frac{D}{\cos(\phi)}\frac{\partial }{\partial \phi}(\cos(\phi)\frac{\partial T_{s}}{\partial \phi})
``` 

- looks like heat equation?
- take global average?

Let's exercise to solve a temperature diffusion equation with climlab:
```{math}
:label: my_label144
C(\phi)\frac{\partial T_{s}}{\partial t} = \frac{D}{\cos(\phi)}\frac{\partial }{\partial \phi}(\cos(\phi)\frac{\partial T_{s}}{\partial \phi})
```

```{code-cell} ipython3
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import climlab
from climlab import constants as const

#  First define an initial temperature field
#   that is warm at the equator and cold at the poles
#   and varies smoothly with latitude in between

from climlab.utils import legendre
sfc = climlab.domain.zonal_mean_surface(num_lat=90, water_depth=10.)
lat = sfc.lat.points
initial = 12. - 40. * legendre.P2(np.sin(np.deg2rad(lat)))

fig, ax = plt.subplots()
ax.plot(lat, initial)
ax.set_xlabel('Latitude')
ax.set_ylabel('Temperature (deg C)');

##  Set up the climlab diffusion process

# make a copy of initial so that it remains unmodified
Ts = climlab.Field(np.array(initial), domain=sfc)
# thermal diffusivity in W/m**2/degC
D = 0.55
# create the climlab diffusion process
#  setting the diffusivity and a timestep of ONE MONTH
d = climlab.dynamics.MeridionalHeatDiffusion(name='Diffusion', 
            state=Ts, D=D, timestep=const.seconds_per_month)
print(d)

#  We are going to step forward one month at a time
#  and store the temperature each time
n_iter = 5
temp = np.zeros((Ts.size, n_iter+1))
temp[:, 0] = np.squeeze(Ts)
for n in range(n_iter):
    d.step_forward()
    temp[:, n+1] = np.squeeze(Ts)

#  Now plot the temperatures
fig,ax = plt.subplots()
ax.plot(lat, temp)
ax.set_xlabel('Latitude')
ax.set_ylabel('Temperature (deg C)')
ax.legend(range(n_iter+1)); ax.grid();

```

Here we use the 2nd Legendre polynomial:
```{math}
:label: my_label145
P_{2}(x) = \frac{1}{2}(3x^{2}-1),
```
where $x=\sin(\phi)$.

```{code-cell} ipython3
x = np.linspace(-1,1)
fig,ax = plt.subplots()
ax.plot(x, legendre.P2(x))
ax.set_title('$P_2(x)$')

```

## Parameterization for the radiative processes
For shortwave, we assume the planetary albedo is fixed:
```{math}
:label: my_label146
ASR(\phi) = (1-\alpha)Q(\phi)
```

For longwave radiation, we simplifiy it with a linear form:
```{math}
:label: my_label147
OLR(\phi) = A+B\times T_{s}(\phi)
```

- A in units of W/m$^2$. What is this? (forcing?)
- B in units of W/m$^2$/K. What is this? (feedback parameter?)

We need to find A and B. We now fit their values using reanalysis data.

```{code-cell} ipython3
import xarray as xr
## The NOAA ESRL server is shutdown! January 2019
ncep_url = "/Users/yuchiaol_ntuas/Desktop/ebooks/data/"
#ncep_url = "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/"
ncep_Ts = xr.open_dataset( ncep_url + "skt.sfc.mon.1981-2010.ltm.nc", decode_times=False)
#url = 'http://apdrc.soest.hawaii.edu:80/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/'
#ncep_Ts = xr.open_dataset(url + 'surface_gauss/skt')
lat_ncep = ncep_Ts.lat; lon_ncep = ncep_Ts.lon
print( ncep_Ts)

#  Take the annual and zonal average!
Ts_ncep_annual = ncep_Ts.skt.mean(dim=('lon','time'))

# TOA radiation data
ncep_ulwrf = xr.open_dataset( ncep_url + "ulwrf.ntat.mon.1981-2010.ltm.nc", decode_times=False)
ncep_dswrf = xr.open_dataset( ncep_url + "dswrf.ntat.mon.1981-2010.ltm.nc", decode_times=False)
ncep_uswrf = xr.open_dataset( ncep_url + "uswrf.ntat.mon.1981-2010.ltm.nc", decode_times=False)
OLR_ncep_annual = ncep_ulwrf.ulwrf.mean(dim=('lon','time'))
ASR_ncep_annual = (ncep_dswrf.dswrf - ncep_uswrf.uswrf).mean(dim=('lon','time'))

#  Use a linear regression package to compute best fit for the slope and intercept
from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(Ts_ncep_annual, OLR_ncep_annual)

print( 'Best fit is A = %0.0f W/m2 and B = %0.1f W/m2/degC' %(intercept, slope))

#  More standard values
A = 210.
B = 2.

fig, ax1 = plt.subplots(figsize=(8,6))
ax1.plot( Ts_ncep_annual, OLR_ncep_annual, 'o' , label='data')
ax1.plot( Ts_ncep_annual, intercept + slope * Ts_ncep_annual, 'k--', label='best fit')
ax1.plot( Ts_ncep_annual, A + B * Ts_ncep_annual, 'r--', label='B=2')
ax1.set_xlabel('Surface temperature (C)', fontsize=16)
ax1.set_ylabel('OLR (W m$^{-2}$)', fontsize=16)
ax1.set_title('OLR versus surface temperature from NCEP reanalysis', fontsize=18)
ax1.legend(loc='upper left')
ax1.grid()

```

```{note}
- Try take global average before regression. What A and B do you get? 

- And if global average temperature is 288 K, what global-averaged OLR do you get?

- $B=2$, is this consistent to the total feedback parameter we obtain before?

```

## Specification for the albedo
```{code-cell} ipython3
days = np.linspace(1.,50.)/50 * const.days_per_year
Qann_ncep = climlab.solar.insolation.daily_insolation(lat_ncep, days ).mean(dim='day')
albedo_ncep = 1 - ASR_ncep_annual / Qann_ncep

albedo_ncep_global = np.average(albedo_ncep, weights=np.cos(np.deg2rad(lat_ncep)))

print( 'The annual, global mean planetary albedo is %0.3f' %albedo_ncep_global)
fig,ax = plt.subplots()
ax.plot(lat_ncep, albedo_ncep)
ax.grid();
ax.set_xlabel('Latitude')
ax.set_ylabel('Albedo');

# Add a new curve to the previous figure
a0 = albedo_ncep_global
a2 = 0.25
ax.plot(lat_ncep, a0 + a2 * legendre.P2(np.sin(np.deg2rad(lat_ncep))))

```

We use the 2nd Legendre polynomial:
```{math}
:label: my_label148
\alpha(\phi) = \alpha_{0} + \alpha_{2}P_{2}(\sin(\phi))
```

```{math}
:label: my_label149
\alpha_{0} = 0.354, \mbox{    } \alpha_{2}=0.25
```

## Annual-mean EBM
Let's put everything together and consider annual-mean, that is $Q(\phi,t)=Q(\phi)$ (no seasonal cycle).

```{math}
:label: my_label150
C(\phi)\frac{\partial T_{s}}{\partial t} = (1-\alpha(\phi))Q(\phi) - (A+BT_{s}) + \frac{D}{\cos(\phi)}\frac{\partial }{\partial \phi}(\cos(\phi)\frac{\partial T_{s}}{\partial \phi})
```

In equilibrium state, we can drop tendency term:
```{math}
:label: my_label151
0 = (1-\alpha(\phi))Q(\phi) - (A+BT_{s}) + \frac{D}{\cos(\phi)}\frac{\partial }{\partial \phi}(\cos(\phi)\frac{\partial T_{s}}{\partial \phi})
```

Let $x=\sin(\phi)$, we can get:
```{math}
:label: my_label152
\frac{D}{B}\frac{d}{dx}((1-x^{2})\frac{dT_{s}}{dx}) - T_{s} = -\frac{(1-\alpha(x))Q(x)-A}{B}
```

```{note}
D/B is a very important parameter for the efficiency of heat transport.
```

```{code-cell} ipython3
#  Some imports needed to make and display animations
from IPython.display import HTML
from matplotlib import animation

def setup_figure():
    templimits = -20,32
    radlimits = -340, 340
    htlimits = -6,6
    latlimits = -90,90
    lat_ticks = np.arange(-90,90,30)

    fig, axes = plt.subplots(3,1,figsize=(8,10))
    axes[0].set_ylabel('Temperature (deg C)')
    axes[0].set_ylim(templimits)
    axes[1].set_ylabel('Energy budget (W m$^{-2}$)')
    axes[1].set_ylim(radlimits)
    axes[2].set_ylabel('Heat transport (PW)')
    axes[2].set_ylim(htlimits)
    axes[2].set_xlabel('Latitude')
    for ax in axes: ax.set_xlim(latlimits); ax.set_xticks(lat_ticks); ax.grid()
    fig.suptitle('Diffusive energy balance model with annual-mean insolation', fontsize=14)
    return fig, axes

def initial_figure(model):
    #  Make figure and axes
    fig, axes = setup_figure()
    # plot initial data
    lines = []
    lines.append(axes[0].plot(model.lat, model.Ts)[0])
    lines.append(axes[1].plot(model.lat, model.ASR, 'k--', label='SW')[0])
    lines.append(axes[1].plot(model.lat, -model.OLR, 'r--', label='LW')[0])
    lines.append(axes[1].plot(model.lat, model.net_radiation, 'c-', label='net rad')[0])
    lines.append(axes[1].plot(model.lat, model.heat_transport_convergence, 'g--', label='dyn')[0])
    lines.append(axes[1].plot(model.lat, 
            model.net_radiation+model.heat_transport_convergence, 'b-', label='total')[0])
    axes[1].legend(loc='upper right')
    lines.append(axes[2].plot(model.lat_bounds, model.heat_transport)[0])
    lines.append(axes[0].text(60, 25, 'Day 0'))
    return fig, axes, lines

def animate(day, model, lines):
    model.step_forward()
    #  The rest of this is just updating the plot
    lines[0].set_ydata(model.Ts)
    lines[1].set_ydata(model.ASR)
    lines[2].set_ydata(-model.OLR)
    lines[3].set_ydata(model.net_radiation)
    lines[4].set_ydata(model.heat_transport_convergence)
    lines[5].set_ydata(model.net_radiation+model.heat_transport_convergence)
    lines[6].set_ydata(model.heat_transport)
    lines[-1].set_text('Day {}'.format(int(model.time['days_elapsed'])))
    return lines   

#  A model starting from isothermal initial conditions
e = climlab.EBM_annual()
e.Ts[:] = 15.  # in degrees Celsius
e.compute_diagnostics()

#  Plot initial data
fig, axes, lines = initial_figure(e)

ani = animation.FuncAnimation(fig, animate, frames=np.arange(1, 100), fargs=(e, lines))

HTML(ani.to_html5_video())

```

We try an example using climlab:
```{code-cell} ipython3
D = 0.1
model = climlab.EBM_annual(name='EBM', A=210, B=2, D=D, a0=0.354, a2=0.25)
print(model)

# The model object stores a dictionary of important parameters
print(model.param)

# Run it out long enough to reach equilibrium
model.integrate_years(10)

fig, axes = plt.subplots(1,2, figsize=(12,4))
ax = axes[0]
ax.plot(model.lat, model.Ts, label=('D = %0.1f' %D))
ax.plot(lat_ncep, Ts_ncep_annual, label='obs')
ax.set_ylabel('Temperature (degC)')
ax = axes[1]
energy_in = np.squeeze(model.ASR - model.OLR)
ax.plot(model.lat, energy_in, label=('D = %0.1f' %D))
ax.plot(lat_ncep, ASR_ncep_annual - OLR_ncep_annual, label='obs')
ax.set_ylabel('Net downwelling radiation at TOA (W m$^{-2}$)')
for ax in axes:
    ax.set_xlabel('Latitude'); ax.legend(); ax.grid();

def inferred_heat_transport( energy_in, lat_deg ):
    '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
    from scipy import integrate
    from climlab import constants as const
    lat_rad = np.deg2rad( lat_deg )
    return ( 1E-15 * 2 * np.pi * const.a**2 * 
            integrate.cumulative_trapezoid( np.cos(lat_rad)*energy_in,
            x=lat_rad, initial=0. ) )

fig, ax = plt.subplots()
ax.plot(model.lat, inferred_heat_transport(energy_in, model.lat), label=('D = %0.1f' %D))
ax.set_ylabel('Heat transport (PW)')
ax.legend(); ax.grid()
ax.set_xlabel('Latitude')

```

What do we see?


## Finding optimal diffusivity $D$
We would like to find a $D$ such that temperature gradient between the pole and equator $\Delta T = 45$ K, and peak heat transport $H = 5.5$ PW.

```{code-cell} ipython3
Darray = np.arange(0., 2.05, 0.05)

model_list = []
Tmean_list = []
deltaT_list = []
Hmax_list = []

for D in Darray:
    ebm = climlab.EBM_annual(A=210, B=2, a0=0.354, a2=0.25, D=D)
    ebm.integrate_years(20., verbose=False)
    Tmean = ebm.global_mean_temperature()
    deltaT = np.max(ebm.Ts) - np.min(ebm.Ts)
    energy_in = np.squeeze(ebm.ASR - ebm.OLR)
    Htrans = ebm.heat_transport
    Hmax = np.max(Htrans)
    model_list.append(ebm)
    Tmean_list.append(Tmean)
    deltaT_list.append(deltaT)
    Hmax_list.append(Hmax)

color1 = 'b'
color2 = 'r'

fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(111)
ax1.plot(Darray, deltaT_list, color=color1)
ax1.plot(Darray, Tmean_list, 'b--')
ax1.set_xlabel(r'D (W m$^{-2}$ K$^{-1}$)', fontsize=14)
ax1.set_xticks(np.arange(Darray[0], Darray[-1], 0.2))
ax1.set_ylabel(r'$\Delta T$ (equator to pole)', fontsize=14,  color=color1)
for tl in ax1.get_yticklabels():
    tl.set_color(color1)
ax2 = ax1.twinx()
ax2.plot(Darray, Hmax_list, color=color2)
ax2.set_ylabel('Maximum poleward heat transport (PW)', fontsize=14, color=color2)
for tl in ax2.get_yticklabels():
    tl.set_color(color2)
ax1.set_title('Effect of diffusivity on temperature gradient and heat transport in the EBM', fontsize=16)
ax1.grid()

ax1.plot([0.6, 0.6], [0, 140], 'k-');

```

## State-dependent diffusivity
```{figure} /_static/lecture_specific/lecture1_figures/chang_merlis_2023_tmp1.png
---
scale: 60%
---
EBM solutions with climate invariant diffusivity. (a) Change in surface temperature DT vs latitude for
numerical solutions [Eq. (1); solid lines] and analytic theory (dashed lines) for different values of control diffusivity
D indicated in the legend. (b) Change in the second-order Legendre polynomial component of temperature DT2 vs
D for numerical solutions [Eq. (1)] and the analytic theory [Eq. (7) scaled by the global-mean temperature change
DT0]. (c) As in (a), but for the numerical solutions for the EBM with the linearized MSE approximation, Eq. (5).
(d) As in (b), but for the change in the second-order Legendre polynomial component of MSE Dh2 and the analytic
theory [Eq. (10) scaled by DT0]. See section 2b of the text for the detailed calculations of the theoretical
estimates for DT in (a) and (c). Source: [Change and Merlis (2023)](https://journals.ametsoc.org/view/journals/clim/36/22/JCLI-D-23-0121.1.xml?tab_body=pdf)

```



