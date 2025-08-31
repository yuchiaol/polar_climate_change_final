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

(syllabus)=

# 1. Syllabus

2025 fall semester; 
Credit: 4; 
Class Number: AtmScixxxx

## Instructor
Professor Yu-Chiao Liang (梁禹喬)<br>
Contacting email: yuchiaoliang@ntu.edu.tw<br>
Office phone: 02-3366-3907<br>

### Teaching Assistants
Ya-Fan Chung (鍾雅帆) <br>
b10209012@ntu.edu.tw


### Office Hours
by appointment

### Location and Time
Monday 3:30-5:10 pm and Wednesday 10:20 am-12:10

### Course website
[https://yuchiaol.github.io/polar_climate_change_final/docs/index.html](https://yuchiaol.github.io/polar_climate_change_final/docs/index.html)


### Grading
- Homework assignments (60%)
- Mid-term presentation (15%) <br>
-- review one key paper <br>
-- describe your final project topic and your hypothesis <br>
-- design numerical experiments or analysis approach <br>
-- you should be able to run your model! <br>
- Final presentation (15%+5% report)<br>
-- test your hypothesis<br>
-- present your findings based on your experiments<br>
- Class interaction (5%)<br>

### Course Description

(in Chinese) 本課程探討極區的氣候變化。我們不僅介紹極區冰圈以及大氣圈過去幾十年的快速變化，包括海冰，永凍土，降雪，冰河，極地渦旋，也會介紹北極暖
化增強現象的成因與影響，以及探討高緯度氣候變化與中低緯度大氣海洋環流的交互作用。我們會閱讀IPCC最新關於冰圈的特別報告以及最新的研究文獻，來幫助我>們掌握極區最新的研究成果與發展趨勢，同時也會使用不同複雜度的大氣模式，來加深了解極區與中低緯度大氣環流交互作用背後的物理動力機制。

Rapid polar climate change in the past decades was the dominant signature of anthropogenic global warming. This course aims to understand the polar climate change and its interaction with regional and global atmospheric circulations. The first part gives an overview of the fast-changing polar cryosphere and atmosphere, including sea ice, permafrost, snow, glaciers, and stratospheric polar vortex. Students will read IPCC’s AR6 WG1 and Special Report on the Ocean and Cryosphere in a Changing Climate and present main conclusions. The second part discusses the two-way interactions of the polar climate change and lower-latitude (including mid-latitude and tropical) atmospheric and oceanic circulations at various spatial and temporal scales, with an emphasis on the cause and effect of Arctic Amplification. We will use a hierarchical modelling approach with different complexity in attempt to understand the potential atmospheric circulation changes in response to polar warming. The anticipated results will be compared to the observations and the state-of-the-art global climate model simulations, for example CMIP5/6 and Polar Amplification Model Intercomparison Project.

### Course Objectives
This course aims to 1) understand the polar climate change and its interaction with regional and global atmospheric circulations, 2) cover some materials of Chapter 12 - Middle Atmosphere Dynamics - in Holton’s “An Introduction to Dynamic Meteorology” (I will try!), and 3) train students to use simple atmospheric models to investigate Arctic-midlatitude connections.

### Course Requirements
Willingness to lead discussion for reading materials and participate in group cooperation, and basic FORTRAN programing and plotting skills (Python preferred because we will use a machine learning package written in Python!).

### Tentative Topics
- Overview of Arctic
  - Geography, climatology and meteorology
  - A case study for Siberian record-breaking warming
  - Explorations

- Snow and Permafrost

- Glacier and Sea-level Rise

- Sea Ice
  - Dynamics and thermodynamics
  - Modelling and prediction
  - Past, recent, and future changes

- Cause and Effect of Arctic Amplification
  - Local vs remote impacts
  - Mechanisms: climate forcings, climate feedbacks, and poleward energy transport
  - Debates on Arctic-midlatitude linkages

- Polar Stratospheric circulation
  - Polar vortex and stratospheric sudden warming
  - Stratosphere-troposphere coupling

### Tentative Schedule
- Week1 (9/1; 9/3)      – introduction; visitors
- Week2 (9/8; 9/10)     – overview of polar climate change
- Week3 (9/15; 9/17)    – zero-dim energy balance model
- Week4 (9/22; 9/24)    – layered model
- Week5 (10/1)          – RCE
- Week6 (10/8)          – climate sensitivity and feedback
- Week7 (10/13; 10/15)  – midterm presentation
- Week8 (10/20; 10/22)  - RAE
- Week9 (10/27; 10/29)  - one-dim EBM
- Week10 (11/3; 11/5)   – seasonal cycle
- Week11 (11/10; 11/12) – ice-albedo feedback in EBM
- Week12 (10/17; 10/19) – diffusive response; snowball earth
- Week13 (11/24; 11/26) - causes of polar amplification
- Week14 (12/1; 12/3)   – Arctic-midlatitude linkage
- Week15 (12/8; 12/10)  – final presentation
- Week16 (12/15; 12/17) – AGU week

### Relevant Texts and References
- [IPCC AR6 WG1 Report](https://www.ipcc.ch/report/ar6/wg1/#FullReport)
- [IPCC Special Report on the Ocean and Cryosphere in a Changing Climate](https://www.ipcc.ch/srocc)
  - [Chapter 3 Polar Regions](https://www.ipcc.ch/srocc/chapter/chapter-3-2/)

- [An Introduction to Dynamic Meteorology (2012), J. R. Holton and G. J. Hakim](https://www.amazon.com/Introduction-Dynamic-Meteorology-International-Geophysics/dp/0123848660/ref=asc_df_0123848660/?tag=hyprod-20&linkCode=df0&hvadid=312091458201&hvpos=&hvnetw=g&hvrand=14614331955549249595&hvpone=&hvptwo=&hvqmt=&hvdev=c&hvdvcmdl=&hvlocint=&hvlocphy=9003483&hvtargid=pla-465623449605&psc=1&tag=&ref=&adgrpid=63669393113&hvpone=&hvptwo=&hvadid=312091458201&hvpos=&hvnetw=g&hvrand=14614331955549249595&hvqmt=&hvdev=c&hvdvcmdl=&hvlocint=&hvlocphy=9003483&hvtargid=pla-465623449605)

- Arctic Amplification and Feedbacks
  - Review article: [Goosse et al. (2018)](https://www.nature.com/articles/s41467-018-04173-0)
  - Review article: [Previdi et al. (2021)](https://iopscience.iop.org/article/10.1088/1748-9326/ac1c29)
  - Review article: [Taylor et al. (2022)](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2021.758361/full)
  - Temperature feedback: [Pithan and Mauritsen (2014)](https://www.nature.com/articles/ngeo2071)
  - Forcing latitudes: [Stuecker et al. (2018)](https://www.nature.com/articles/s41558-018-0339-y)

- Arctic Amplification and Midlatitude Weather and Climate
  - Review article: [Cohen et al. (2014)](https://www.nature.com/articles/ngeo2234)
  - Key paper: [Francis and Vavrus 2012](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012gl051000)
  - Weakened evidence: [Blackport and Screen (2020)](https://www.nature.com/articles/s41558-020-00954-y)
  - Model consistency and discrepancy: [Screen et al. (2018)](https://www.nature.com/articles/s41561-018-0059-y)
  - Divergent consensus [Cohen et al. (2019)](https://www.nature.com/articles/s41558-019-0662-y)
  - Robust but weak response: [Smith et al. (2022)](https://www.nature.com/articles/s41467-022-28283-y)

- Debates on [Mori et al. (2019)](https://www.nature.com/articles/s41558-018-0379-3)
  - First debate: [Screen and Blackport (2019)](https://www.nature.com/articles/s41558-019-0635-1) and [Mori et al. (2019)](https://www.nature.com/articles/s41558-019-0636-0)
  - Second debate: [Zappa et al. (2021)](https://www.nature.com/articles/s41558-020-00982-8) and [Mori et al. (2021)](https://www.nature.com/articles/s41558-020-00983-7)

- [The Climate Laboratory](https://brian-rose.github.io/ClimateLaboratoryBook/home.html)

```{note}
I would like to especially thank Professor Brian E. J. Rose (University at Albany) to kindly share his climate modeling course materials and allow us to use the Python scripts and data. Some materials are obtained from [The Climate Laboratory](https://brian-rose.github.io/ClimateLaboratoryBook/home.html).
```

- Idealized models


