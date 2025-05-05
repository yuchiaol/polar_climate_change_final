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

(ice_albedo_ebm)=

# Ice Albedo Feedback in EBM

```{note}
The Python scripts used below and some materials are modified from Prof. Brian E. J. Rose's climlab [website](https://brian-rose.github.io/ClimateLaboratoryBook/courseware/seasonal-cycle.html).
```

## Interactive ice/snow line (step function)
The idea is to change the albedo as a function of latitude and a threshold temperatre $T_{f} = -10 ^{o}C$. This is historically inheritated from Budyko's model.

```{math}
:label: my_label171
\begin{eqnarray}
\alpha(\phi,T(\phi)) &=& \alpha_{0} + \alpha_{2}P_{2}(\sin(\phi)), &\mbox{  when }& T(\phi) > T_{f}\\
&=& \alpha_{i}, &\mbox{  when }& T(\phi) \le T_{f}
\end{eqnarray}
```
where $P_{2}(sin(\phi))=\frac{1}{2}(3\sin(\phi)^{2}-1)$. 


```{code-cell} ipython3
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import climlab
from climlab import constants as const

#  for convenience, set up a dictionary with our reference parameters
param = {'D':0.55, 'A':210, 'B':2, 'a0':0.3, 'a2':0.078, 'ai':0.62, 'Tf':-10.}
model1 = climlab.EBM_annual(name='EBM with interactive ice line',
                            num_lat=180, 
                            D=0.55, 
                            A=210., 
                            B=2., 
                            Tf=-10., 
                            a0=0.3, 
                            a2=0.078, 
                            ai=0.62)
print(model1)

print(model1.param)

#  A python shortcut... we can use the dictionary to pass lots of input arguments simultaneously:

#  same thing as before, but written differently:
model1 = climlab.EBM_annual(name='EBM with interactive ice line',
                            num_lat=180,
                            **param)
print(model1)

def ebm_plot(e, return_fig=False):    
    templimits = -60,32
    radlimits = -340, 340
    htlimits = -6,6
    latlimits = -90,90
    lat_ticks = np.arange(-90,90,30)
    
    fig = plt.figure(figsize=(8,12))

    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(e.lat, e.Ts)
    ax1.set_ylim(templimits)
    ax1.set_ylabel('Temperature (°C)')
    
    ax2 = fig.add_subplot(3,1,2)
    ax2.plot(e.lat, e.ASR, 'k--', label='SW' )
    ax2.plot(e.lat, -e.OLR, 'r--', label='LW' )
    ax2.plot(e.lat, e.net_radiation, 'c-', label='net rad' )
    ax2.plot(e.lat, e.heat_transport_convergence, 'g--', label='dyn' )
    ax2.plot(e.lat, e.net_radiation + e.heat_transport_convergence, 'b-', label='total' )
    ax2.set_ylim(radlimits)
    ax2.set_ylabel('Energy budget (W m$^{-2}$)')
    ax2.legend()
    
    ax3 = fig.add_subplot(3,1,3)
    ax3.plot(e.lat_bounds, e.heat_transport )
    ax3.set_ylim(htlimits)
    ax3.set_ylabel('Heat transport (PW)')
    
    for ax in [ax1, ax2, ax3]:
        ax.set_xlabel('Latitude')
        ax.set_xlim(latlimits)
        ax.set_xticks(lat_ticks)
        ax.grid()
    
    if return_fig:
        return fig

model1.integrate_years(5)

model1.ASR.to_xarray()

climlab.global_mean(model1.ASR - model1.OLR)

# Integrate out to equilibrium.
model1.integrate_years(5)
#  Check for energy balance
print(climlab.global_mean(model1.net_radiation))
f = ebm_plot(model1)

#  There is a diagnostic that tells us the current location of the ice edge:
model1.icelat

```

## Polar-amplified warming in EBM (Yeah!!!)
If we consider a scenario of doubling CO2, that is to change paraemter $A$ to $A-\delta A$, where $\delta A= 4$ W/m$^2$.

```{code-cell} ipython3
model1.subprocess['LW'].A

deltaA = 4.

#  This is a very handy way to "clone" an existing model:
model2 = climlab.process_like(model1)

#  Now change the longwave parameter:
model2.subprocess['LW'].A = param['A'] - deltaA
#  and integrate out to equilibrium again
model2.integrate_years(5, verbose=False)

plt.plot(model1.lat, model1.Ts, label='model1')
plt.plot(model2.lat, model2.Ts, label='model2')
plt.legend(); plt.grid()

model2.icelat

```

```{code-cell} ipython3

model3 = climlab.process_like(model1)
model3.subprocess['LW'].A = param['A'] - 2*deltaA
model3.integrate_years(5, verbose=False)

plt.plot(model1.lat, model1.Ts, label='model1')
plt.plot(model2.lat, model2.Ts, label='model2')
plt.plot(model3.lat, model3.Ts, label='model3')
plt.xlim(-90, 90)
plt.grid()
plt.legend()

```

For me, this is the minimal model to produce polar amplification. Starting from this model, I want to increase the model hierachy to study the drivers of polar mechanism.

## Solar forcing

```{code-cell} ipython3

m = climlab.EBM_annual(num_lat=180, **param)
#  The current (default) solar constant, corresponding to present-day conditions:
m.subprocess.insolation.S0

#  First, get to equilibrium
m.integrate_years(5.)
#  Check for energy balance
climlab.global_mean(m.net_radiation)

m.icelat

#  Now make the solar constant smaller:
m.subprocess.insolation.S0 = 1300.

#  Integrate to new equilibrium
m.integrate_years(10.)
#  Check for energy balance
climlab.global_mean(m.net_radiation)

m.icelat

ebm_plot(m)

```

## The effect of diffusivity with albedo feedback
We repeat the results of [Stone (1978)](https://www.sciencedirect.com/science/article/pii/0377026578900064?via%3Dihub)

```{code-cell} ipython3
param = {'A':210, 'B':2, 'a0':0.3, 'a2':0.078, 'ai':0.62, 'Tf':-10.}
print( param)

Darray = np.arange(0., 2.05, 0.05)


model_list = []
Tmean_list = []
deltaT_list = []
Hmax_list = []
for D in Darray:
    ebm = climlab.EBM_annual(num_lat=360, D=D, **param )
    ebm.integrate_years(5., verbose=False)
    Tmean = ebm.global_mean_temperature()
    deltaT = np.max(ebm.Ts) - np.min(ebm.Ts)
    HT = np.squeeze(ebm.heat_transport)
    ind = np.where(ebm.lat_bounds==35.5)[0]
    Hmax = HT[ind]
    model_list.append(ebm)
    Tmean_list.append(Tmean)
    deltaT_list.append(deltaT)
    Hmax_list.append(Hmax)

color1 = 'b'
color2 = 'r'

fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(111)
ax1.plot(Darray, deltaT_list, color=color1, label=r'$\Delta T$')
ax1.plot(Darray, Tmean_list, '--', color=color1, label=r'$\overline{T}$')
ax1.set_xlabel('D (W m$^{-2}$ °C$^{-1}$)', fontsize=14)
ax1.set_xticks(np.arange(Darray[0], Darray[-1], 0.2))
ax1.set_ylabel(r'Temperature (°C)', fontsize=14,  color=color1)
for tl in ax1.get_yticklabels():
    tl.set_color(color1)
ax1.legend(loc='center right')
ax2 = ax1.twinx()
ax2.plot(Darray, Hmax_list, color=color2)
ax2.set_ylabel('Poleward heat transport across 35.5° (PW)', fontsize=14, color=color2)
for tl in ax2.get_yticklabels():
    tl.set_color(color2)
ax1.set_title('Effect of diffusivity on EBM with albedo feedback', fontsize=16)
ax1.grid()

```

## Diffusive response to a point source of energy at 45$^o$N

- without albedo feedback

```{code-cell} ipython3
param_noalb = {'A': 210, 'B': 2, 'D': 0.55, 'Tf': -10.0, 'a0': 0.3, 'a2': 0.078}
m1 = climlab.EBM_annual(num_lat=180, **param_noalb)
print(m1)

m1.integrate_years(5.)

m2 = climlab.process_like(m1)

point_source = climlab.process.energy_budget.ExternalEnergySource(state=m2.state, timestep=m2.timestep)
ind = np.where(m2.lat == 45.5)
point_source.heating_rate['Ts'][ind] = 100.

m2.add_subprocess('point source', point_source)
print( m2)

m2.integrate_years(5.)

plt.plot(m2.lat, m2.Ts - m1.Ts)
plt.xlim(-90,90)
plt.xlabel('Latitude')
plt.ylabel(r'$\Delta T_s$')
plt.grid()

```

```{note}
The length scale of the point warming is suggested to be proportional to $\sqrt{\frac{D}{B}}$

```

- with albedo
```{code-cell} ipython3
m3 = climlab.EBM_annual(num_lat=180, **param)
m3.integrate_years(5.)
m4 = climlab.process_like(m3)
point_source = climlab.process.energy_budget.ExternalEnergySource(state=m4.state, timestep=m4.timestep)
point_source.heating_rate['Ts'][ind] = 100.
m4.add_subprocess('point source', point_source)
m4.integrate_years(5.)

plt.plot(m4.lat, m4.Ts - m3.Ts)
plt.xlim(-90,90)
plt.xlabel('Latitude')
plt.ylabel(r'$\Delta T_s$')
plt.grid()

```



## Homework assignment X (due xxx)
- Solar forcing. Please change the solar forcing and discuss the results.

```{code-cell} ipython3

m = climlab.EBM_annual(num_lat=180, **param)
#  The current (default) solar constant, corresponding to present-day conditions:
m.subprocess.insolation.S0

#  First, get to equilibrium
m.integrate_years(5.)
#  Check for energy balance
climlab.global_mean(m.net_radiation)

m.icelat

#  Now make the solar constant smaller:
m.subprocess.insolation.S0 = 1300.

#  Integrate to new equilibrium
m.integrate_years(10.)
#  Check for energy balance
climlab.global_mean(m.net_radiation)

m.icelat

ebm_plot(m)

```




