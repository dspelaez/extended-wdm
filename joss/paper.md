---
title: 'EWDM: A wavelet-based method for estimating directional spectra of ocean waves'
tags:
  - python
  - physical oceanography
  - ocean surface waves
  - directional wave spectra
  - continous wavelet transform
authors:
  - name:
        firstname: Daniel
        surname:  Peláez-Zapata
    orcid: 0000-0001-5862-6194
    affiliation: 1
  - name:
        firstname: Frédéric
        surname: Dias
    orcid: 0000-0002-5123-4929
    affiliation: "1, 2"
affiliations:
 - name: Centre Borelli, École Normale Supérieure Paris-Saclay, France
   index: 1
   ror: 02hcn4061
 - name: School of Mathematics and Statistics, University College Dublin, Ireland
   index: 2
date: 5 April 2025
bibliography: paper.bib

---

# Summary

The accurate estimation of directional wave spectra is crucial for understanding ocean wave dynamics, air-sea interactions, and coastal processes [@Barstow2005]. The directional wave spectrum is often computed by resolving the distribution of wave energy as a function of frequency and direction using mathematical methods applied to time series from measuring instruments (e.g., bouys). These time series typically consist of either triplets of variables at a single point (e.g., velocity or acceleration components) or arrays of measurements distributed across multiple spatial locations (e.g., wave staffs). Traditional methods, based on the conventional Fourier cross-spectral analysis, such as the Truncated Fourier Series (TFS), Maximum Likelihood (MLM), or Maximum Entropy (MEM) methods, often suffer from limitations in spectral resolution, inaccurate estimations, and spurious spectral peaks. This disadvantage is generally attributed to certain simplifications on the shape of the directional distribution and assumptions such as stationarity, which may not always capture the complexity of ocean waves [@Benoit1997]. Wavelet-based methods provide a more flexible approach, thanks to their time-frequency decomposition capabilities [@Donelan1996; @Donelan2015; @PelaezZapata2024a; @Krogstad2006]. This paper presents the EWDM (Extended Wavelet Directional Method), a Python toolkit developed to estimate the directional spectrum of ocean waves based on a wavelet-based technique. EWDM aims to address the limitations of conventional techniques by providing a robust estimation of the directional wave spectrum from both single-point triplet data [@PelaezZapata2024a; @Krogstad2006] and spatially-distributed arrays [@Donelan1996; @Donelan2015]. Consequently, EWDM can be used on diverse sources of data, including GPS buoys, pitch-roll buoys, arrays of wave staffs, acoustic Doppler current profilers (ADCP) and sampled points from stereo-videos of the sea surface. Key features of the EWDM include the implementation of wavelet-based algorithms for extracting directional information from wave time series, improved estimation of wave directional distribution using Kernel Density Estimation (KDE), tools for processing and visualizing directional wave data, and compatibility with popular data sources, including Spotter buoys and CDIP (Coastal Data Information Program) database. The package is powered by [`xarray`](https://github.com/pydata/xarray)'s [@Hoyer2017] labelled multi-dimensional arrays, which enhances its efficiency, scalability and compatibility with scientific tools and workflows.

# Statement of need

Accurate knowledge of how wave energy is distributed across frequencies and directions is critical for understanding key processes in physical oceanography, such as air-sea interactions, long-term wave climate, and upper ocean dynamics, as well as for applications in coastal and ocean engineering, including wave forecasting, safe navigation, assessment of coastal erosion, design of maritime structures and operation of wave energy converters and offshore wind turbines. The directional wave spectrum offers essential insight into this energy distribution [@Barstow2005]. Despite its critical importance, the dynamic and complex nature of ocean waves, combined with limitations in current measurement technologies and analysis techniques, continues to challenge the accurate estimation of directional wave information.

Available wave sensors typically measure at single-point locations (e.g. buoys) or at a few spatially distributed points (e.g., wave staffs). This constrains the information available for inferring wave directionality. The estimation of the directional wave spectrum from single-point data generally requires three perpendicular wave quantities (e.g., velocity or acceleration components), often called triplets. Similarly, when spatial arrays of surface elevation data is available, the directional wave spectrum can be inferred. However, the accuracy of this estimation might be significantly influenced by the sensor distribution, occasionally resulting in unreliable estimations [@Young1994]. The difficulties in obtaining the directional wave spectrum have led to the development of various methods. However, these methods often rely on truncated Fourier series, which struggle to represent complex or multi-modal sea states, or on statistical techniques that impose assumptions that may not hold under real ocean conditions [@PelaezZapata2024a; @Benoit1993; @Benoit1993a; @Benoit1995; @Benoit1997; @Ochoa1990].

@Donelan1996 introduced the WDM (Wavelet Directional Method) as an alternative technique, based on the Continuous Wavelet Transform (CWT), to overcome some of the limitations of conventional methods. However, this method was only applicable to spatial arrays of wave gauges. Building on this, @PelaezZapata2024a and @Krogstad2006 demonstrated the extension of this wavelet-based approach to single-point measurements, offering a promising advancement in wavelet-based directional wave spectrum estimation. This package is therefore an implementation of these algorithms. It has been successfully employed in the recent research by @PelaezZapata2024a, delivering robust and accurate results. This application further demonstrates its effectiveness in addressing the challenges posed by real-world complex sea-state conditions. It is expected to be a valuable tool for oceanographers, engineers, students, and practitioners working in the field of ocean wave analysis.


# State of the field

Numerous software tools are available for analysing ocean wave data, each designed to specific applications and needs. Well-established tools, such as [`WAFO`](https://www.maths.lth.se/matstat/wafo/) [Wave Analysis for Fatigue and Oceanography, @Brodtkorb2000wafo], focus on the statistical analysis of random waves, fatigue analysis, wave-induced loads, and reliability assessments in ocean engineering. Alternative, recent software packager, such as [`FOWD`](https://github.com/dionhaefner/FOWD) [Free Ocean Wave Dataset for Data Mining and Machine Learning, @Hafner2021fowd], provide methods for processing and analysis wave elevation time series, including the estimation of zero-crossing parameters, such as crest-to-trough height, skewness and kurtosis; and spectral parameters, such as significant wave height, peak and mean periods, directional spreading and spectral bandwidths. [`FOWD`](https://github.com/dionhaefner/FOWD) focuses on providing an extensive dataset of wave parameters optimised for machine learning and data mining applications. While these packages facilitate the estimation of certain directional parameters, they do not have the functionality to calculate the directional spectrum from wave observations.

[`DIWASP`](https://github.com/metocean/diwasp) [DIrectional WAve SPectrum analysis, @Johnson2002diwasp], is a widely adopted MATLAB-based toolbox specifically designed for this task. [`DIWASP`](https://github.com/metocean/diwasp) implements a variety of routines for directional wave spectrum estimation, such as the Iterative Maximum Likelihood Method (IMLM), the Extended Maximum Entropy Principle (EMEP), the Bayesian Directional Method (BDM), among others. However, all these methods rely on conventional Fourier cross-spectral analysis, and therefore, are subject to the limitations previously discussed.

Modern Python-based tools, such as [`wavespectra`](https://github.com/wavespectra/wavespectra) and [`oceanwaves-python`](https://github.com/openearth/oceanwaves-python), provide an abstraction layer for typical directional wave spectra manipulation. These packages provide functionalities for common operations, such as reading and writing spectral data from numerical models (e.g., SWAN, WAVEWATCH III) and wave instruments (e.g., Spotter buoys, TRIAXYS, Datawell), spectral partitioning, estimation of spectral parameters, standard conversions and advanced plotting. They do not incorporate, however, algorithms for estimating the directional wave spectrum from raw wave measurements. EWDM complements the ecosystem, providing routines for estimating the directional wave spectrum directly from measurements of wave time series. The outputs obtained from EWDM can potentially be integrated with these tools to broaden its functionalities, opening up new possibilities to the scientific community.


# Acknowledgements

This work was funded by the European Research Council (ERC) under the EU Horizon 2020 research and innovation program (Grant Agreement 833125- HIGHWAVE) and by the Centre National d'Études Spatiales (CNES) under the project ARANSAT. DPZ would like to thank the maintainers and developers of the open-source scientific Python ecosystem for making data processing and visualization more efficient. DPZ is grateful to the Navier (ENS Paris-Saclay) and Stokes (University College Dublin) teams for the interesting discussions.


# References
