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

(polar_amplification)=

# A1. Polar Amplification in simple models

```{note}
The materials and scripts are obtained from climlab tutorial: [here](https://climlab.readthedocs.io/en/latest/courseware/PolarAmplification.html)
```

```{code-cell} ipython3
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import climlab
from climlab import constants as const
```

## EBM with surface and atmosphere layers

```{code-cell} ipython3
ebm = climlab.GreyRadiationModel(num_lev=1, num_lat=90)
insolation = climlab.radiation.AnnualMeanInsolation(domains=ebm.Ts.domain)
ebm.add_subprocess('insolation', insolation)
ebm.subprocess.SW.flux_from_space = ebm.subprocess.insolation.insolation
print(ebm)

#  add a fixed relative humidity process
#  (will only affect surface evaporation)
h2o = climlab.radiation.ManabeWaterVapor(state=ebm.state, **ebm.param)
ebm.add_subprocess('H2O', h2o)

#  Add surface heat fluxes
shf = climlab.surface.SensibleHeatFlux(state=ebm.state, Cd=3E-4)
lhf = climlab.surface.LatentHeatFlux(state=ebm.state, Cd=3E-4)
# couple water vapor to latent heat flux process
lhf.q = h2o.q
ebm.add_subprocess('SHF', shf)
ebm.add_subprocess('LHF', lhf)

ebm.integrate_years(1)

plt.plot(ebm.lat, ebm.Ts)
plt.plot(ebm.lat, ebm.Tatm)
```

```{code-cell} ipython3
co2ebm = climlab.process_like(ebm)
co2ebm.subprocess['LW'].absorptivity = ebm.subprocess['LW'].absorptivity*1.1

co2ebm.integrate_years(3.)

#  no heat transport but with evaporation -- no polar amplification
plt.plot(ebm.lat, co2ebm.Ts - ebm.Ts)
plt.plot(ebm.lat, co2ebm.Tatm - ebm.Tatm)

```

## Include meridional heat transport
```{code-cell} ipython3
diffebm = climlab.process_like(ebm)
# thermal diffusivity in W/m**2/degC
D = 0.6
# meridional diffusivity in m**2/s
K = D / diffebm.Tatm.domain.heat_capacity * const.a**2
d = climlab.dynamics.MeridionalDiffusion(K=K, state={'Tatm': diffebm.Tatm}, **diffebm.param)
diffebm.add_subprocess('diffusion', d)
print(diffebm)

diffebm.integrate_years(3)

plt.plot(diffebm.lat, diffebm.Ts)
plt.plot(diffebm.lat, diffebm.Tatm)

```

```{code-cell} ipython3
def inferred_heat_transport( energy_in, lat_deg ):
    '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
    from scipy import integrate
    from climlab import constants as const
    lat_rad = np.deg2rad( lat_deg )
    return ( 1E-15 * 2 * np.math.pi * const.a**2 * integrate.cumtrapz( np.cos(lat_rad)*energy_in,
            x=lat_rad, initial=0. ) )

#  Plot the northward heat transport in this model
Rtoa = np.squeeze(diffebm.timeave['ASR'] - diffebm.timeave['OLR'])
plt.plot(diffebm.lat, inferred_heat_transport(Rtoa, diffebm.lat))

```

```{code-cell} ipython3
##  Now warm it up!
co2diffebm = climlab.process_like(diffebm)
co2diffebm.subprocess['LW'].absorptivity = diffebm.subprocess['LW'].absorptivity*1.1

co2diffebm.integrate_years(5)

#  with heat transport and evaporation
#   Get some modest polar amplifcation of surface warming
#    but larger equatorial amplification of atmospheric warming
#    Increased atmospheric gradient = increased poleward flux.
plt.plot(diffebm.lat, co2diffebm.Ts - diffebm.Ts, label='Ts')
plt.plot(diffebm.lat, co2diffebm.Tatm - diffebm.Tatm, label='Tatm')
plt.legend()

```

```{code-cell} ipython3
Rtoa = np.squeeze(diffebm.timeave['ASR'] - diffebm.timeave['OLR'])
Rtoa_co2 = np.squeeze(co2diffebm.timeave['ASR'] - co2diffebm.timeave['OLR'])
plt.plot(diffebm.lat, inferred_heat_transport(Rtoa, diffebm.lat), label='1xCO2')
plt.plot(diffebm.lat, inferred_heat_transport(Rtoa_co2, diffebm.lat), label='2xCO2')
plt.legend()

```

## Exclude evaporation
```{code-cell} ipython3
diffebm2 = climlab.process_like(diffebm)
diffebm2.remove_subprocess('LHF')
diffebm2.integrate_years(3)
co2diffebm2 = climlab.process_like(co2diffebm)
co2diffebm2.remove_subprocess('LHF')
co2diffebm2.integrate_years(3)

#  With transport and no evaporation...
#  No polar amplification, either of surface or air temperature!
plt.plot(diffebm2.lat, co2diffebm2.Ts - diffebm2.Ts, label='Ts')
plt.plot(diffebm2.lat, co2diffebm2.Tatm[:,0] - diffebm2.Tatm[:,0], label='Tatm')
plt.legend()
plt.figure()

#  And in this case, the lack of polar amplification is DESPITE an increase in the poleward heat transport.
Rtoa = np.squeeze(diffebm2.timeave['ASR'] - diffebm2.timeave['OLR'])
Rtoa_co2 = np.squeeze(co2diffebm2.timeave['ASR'] - co2diffebm2.timeave['OLR'])
plt.plot(diffebm2.lat, inferred_heat_transport(Rtoa, diffebm2.lat), label='1xCO2')
plt.plot(diffebm2.lat, inferred_heat_transport(Rtoa_co2, diffebm2.lat), label='2xCO2')
plt.legend()

```

## Column model

```{code-cell} ipython3
model = climlab.GreyRadiationModel(num_lev=30, num_lat=90, abs_coeff=1.6E-4)
insolation = climlab.radiation.AnnualMeanInsolation(domains=model.Ts.domain)
model.add_subprocess('insolation', insolation)
model.subprocess.SW.flux_from_space = model.subprocess.insolation.insolation
print(model)

#  Convective adjustment for atmosphere only
conv = climlab.convection.ConvectiveAdjustment(state={'Tatm':model.Tatm}, adj_lapse_rate=6.5,
                                                       **model.param)
model.add_subprocess('convective adjustment', conv)

#  add a fixed relative humidity process
#  (will only affect surface evaporation)
h2o = climlab.radiation.water_vapor.ManabeWaterVapor(state=model.state, **model.param)
model.add_subprocess('H2O', h2o)

#  Add surface heat fluxes
shf = climlab.surface.SensibleHeatFlux(state=model.state, Cd=1E-3)
lhf = climlab.surface.LatentHeatFlux(state=model.state, Cd=1E-3)
lhf.q = model.subprocess.H2O.q
model.add_subprocess('SHF', shf)
model.add_subprocess('LHF', lhf)

model.integrate_years(3.)

def plot_temp_section(model, timeave=True):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if timeave:
        field = model.timeave['Tatm'].transpose()
    else:
        field = model.Tatm.transpose()
    cax = ax.contourf(model.lat, model.lev, field)
    ax.invert_yaxis()
    ax.set_xlim(-90,90)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    fig.colorbar(cax)

plot_temp_section(model, timeave=False)

```

```{code-cell} ipython3
co2model = climlab.process_like(model)
co2model.subprocess['LW'].absorptivity = model.subprocess['LW'].absorptivity*1.1

co2model.integrate_years(3)

plot_temp_section(co2model, timeave=False)

#  Without transport, get equatorial amplification
plt.figure()
plt.plot(model.lat, co2model.Ts - model.Ts, label='Ts')
plt.plot(model.lat, co2model.Tatm[:,0] - model.Tatm[:,0], label='Tatm')
plt.legend()

```

## Include meridional heat transport

```{code-cell} ipython3
diffmodel = climlab.process_like(model)

# thermal diffusivity in W/m**2/degC
D = 0.05
# meridional diffusivity in m**2/s
K = D / diffmodel.Tatm.domain.heat_capacity[0] * const.a**2
print(K)

d = climlab.dynamics.MeridionalDiffusion(K=K, state={'Tatm':diffmodel.Tatm}, **diffmodel.param)
diffmodel.add_subprocess('diffusion', d)
print(diffmodel)

diffmodel.integrate_years(3)

plot_temp_section(diffmodel)

#  Plot the northward heat transport in this model
plt.figure()
Rtoa = np.squeeze(diffmodel.timeave['ASR'] - diffmodel.timeave['OLR'])
plt.plot(diffmodel.lat, inferred_heat_transport(Rtoa, diffmodel.lat))

```

```{code-cell} ipython3
##  Now warm it up!
co2diffmodel = climlab.process_like(diffmodel)
co2diffmodel.subprocess['LW'].absorptivity = diffmodel.subprocess['LW'].absorptivity*1.1

co2diffmodel.integrate_years(3)

#  With transport, get polar amplification...
#   of surface temperature, but not of air temperature!
plt.plot(diffmodel.lat, co2diffmodel.Ts - diffmodel.Ts, label='Ts')
plt.plot(diffmodel.lat, co2diffmodel.Tatm[:,0] - diffmodel.Tatm[:,0], label='Tatm')
plt.legend()


Rtoa = np.squeeze(diffmodel.timeave['ASR'] - diffmodel.timeave['OLR'])
plt.figure()
Rtoa_co2 = np.squeeze(co2diffmodel.timeave['ASR'] - co2diffmodel.timeave['OLR'])
plt.plot(diffmodel.lat, inferred_heat_transport(Rtoa, diffmodel.lat), label='1xCO2')
plt.plot(diffmodel.lat, inferred_heat_transport(Rtoa_co2, diffmodel.lat), label='2xCO2')

```

## Exclude evaporation
```{code-cell} ipython3
diffmodel2 = climlab.process_like(diffmodel)
diffmodel2.remove_subprocess('LHF')
print(diffmodel2)

diffmodel2.integrate_years(3)
```

```{code-cell} ipython3
co2diffmodel2 = climlab.process_like(co2diffmodel)
co2diffmodel2.remove_subprocess('LHF')
co2diffmodel2.integrate_years(3)

```

```{code-cell} ipython3
#  With transport and no evaporation...
#  No polar amplification, either of surface or air temperature!
plt.plot(diffmodel2.lat, co2diffmodel2.Ts - diffmodel2.Ts, label='Ts')
plt.plot(diffmodel2.lat, co2diffmodel2.Tatm[:,0] - diffmodel2.Tatm[:,0], label='Tatm')
plt.legend()


Rtoa = np.squeeze(diffmodel2.timeave['ASR'] - diffmodel2.timeave['OLR'])
Rtoa_co2 = np.squeeze(co2diffmodel2.timeave['ASR'] - co2diffmodel2.timeave['OLR'])
plt.figure()
plt.plot(diffmodel2.lat, inferred_heat_transport(Rtoa, diffmodel2.lat), label='1xCO2')
plt.plot(diffmodel2.lat, inferred_heat_transport(Rtoa_co2, diffmodel2.lat), label='2xCO2')

```



