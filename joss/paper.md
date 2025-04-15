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

The accurate estimation of directional wave spectra is crucial for understanding ocean wave dynamics, air-sea interactions, and coastal processes [@Barstow_etal_2005]. Traditional methods, based on the conventional Fourier transform, such as the Truncated Fourier Series (TFS), Maximum Likelihood (MLM), or Maximum Entropy (MEM) methods, often suffer from limitations in spectral resolution, inaccurate estimations, and spurious spectral peaks. This disadvantage is generally attributed to certain simplifications on the shape of the directional distribution and assumptions such as stationarity, which may not accurately capture the complexity of ocean waves [@Benoit_etal_1997]. Wavelet-based methods provide a more flexible approach, thanks to their time-frequency decomposition capabilities [@Donelan_etal_1996; @PelaezZapata_2024a; @Krogstad_etal_2006]. This paper presents the EWDM (Extended Wavelet Directional Method), a Python toolkit developed to estimate the directional spectra of ocean waves using the Continuous Wavelet Transform (CWT). EWDM aims to address the limitations of conventional techniques by providing a robust estimation of the directional wave spectrum from both single-point triplet data [@PelaezZapata_2024a; @Krogstad_etal_2006] and spatially-distributed arrays, extending the capabilities of the original WDM [@Donelan_etal_1996], which was initially proposed to work only on arrays. Consequently, EWDM can be used on diverse sources of data, including GPS buoys, pitch-roll buoys, arrays of wave staffs, stereo-images and acoustic Doppler current profilers (ADCP). Key features of the EWDM include the implementation of wavelet-based algorithms for extracting directional information from wave time series, improved estimation of wave directional distribution using Kernel Density Estimation (KDE), tools for processing and visualizing directional wave data, and compatibility with popular data sources such as the so-called Spotter buoys and the directional waveriders (DWR) from the CDIP (Coastal Data Information Program) database. The package is powered by `xarray`'s [@Hoyer_2017] labelled multi-dimensional arrays, which enhances its efficiency and scalability.


# Statement of need

Various aspects of physical oceanography, including air-sea interactions, wave-induced mixing, and upper layer dynamics, as well as several applications in coastal engineering, such as wave forecasting, prediction of rogue waves, quantification of coastal erosion, and the design of offshore structures and renewable energy devices, require accurate and precise knowledge of how wave energy is distributed across frequencies and directions [@Barstow_etal_2005]. The directional wave spectrum provides essential information about this energy distribution. However, the dynamic and unpredictable nature of ocean waves, coupled with the constraints of current measurement technologies and analysis methods, poses challenges in attaining accurate directional information.

The directional wave spectra can be estimated from measurements using either spatial arrays of wave gauges or single-point triplets, where three perpendicular wave quantities are measured at the same point [@Barstow_etal_2005]. This latter is, for instance, the case of wave buoys, which are the most widespread means for obtaining _in-situ_ wave data. However, current methods rely either on truncated Fourier series that are unable to capture complex sea states or on statistical fitting of the distribution of wave directions and often involve assumptions that may not hold in real-world conditions [@PelaezZapata_2024a; @Benoit_1993; @Benoit_1993a; @Benoit_and_Teisson_1995; @Benoit_etal_1997; @Ochoa_and_Delgado-Gonzalez_199].

@Donelan_etal_1996 introduced an alternative method based on the Continuous Wavelet Transform (CWT) to overcome some of the limitations of classic methods, though specifically applicable to spatial arrays of wave gauges. Building on this, @PelaezZapata_2024a and @Krogstad_etal_2006 demonstrated the extension of this wavelet-based approach to single-point measurements, offering a promising advancement in directional wave spectrum estimation. This package is therefore an implementation of these algorithms. It has been successfully employed in the recent research by @PelaezZapata_2024a, delivering robust and accurate results. This application further demonstrates its effectiveness in addressing the challenges posed by real-world complex sea-state conditions. It is expected to be a valuable tool for oceanographers, engineers, students, and practitioners working in the field of ocean wave analysis.


# State of the field

Numerous software tools are available for analysing ocean wave data, each designed to specific applications and needs. Well-established tools, such as [`WAFO`](https://www.maths.lth.se/matstat/wafo/) (Wave Analysis for Fatigue and Oceanography; @Brodtkorb_2000_wafo) focuses on the statistical analysis of random waves, fatigue analysis, wave-induced loads, and reliability assessments in ocean engineering. Recent software tools, such as [`FOWD`](https://github.com/dionhaefner/FOWD) (Free Ocean Wave Dataset for Data Mining and Machine Learning; @Hafner_2021_fowd) provides methods for processing and analysis wave elevation time-series, including the estimation of zero-crossing parameters, such as crest-to-trough height, skewness and kurtosis; and spectral parameters, such as significant wave height, peak and mean periods, directional spreading and spectral bandwidths. FOWD focuses on providing an extensive dataset of wave parameters optimised for machine learning and data mining applications. While these packages facilitate the estimation of certain directional parameters, they do not have the functionality to calculate the directional spectrum from wave observations.

[`DIWASP`](https://github.com/metocean/diwasp) (DIrectional WAve SPectrum analysis; @Johnson_2002_diwasp), is a widely adopted MATLAB-based toolbox specifically designed for this task. DIWASP implements a variety of algorithms for directional wave spectrum estimation, such as the Iterative Maximum Likelihood Method (IMLM), the Extended Maximum Entropy Principle (EMEP), the Bayesian Directional Method (BDM), among others. However, all these methods rely on conventional Fourier cross-spectral analysis, and therefore, are subject to the limitations previously discussed.

Modern Python-based tools, such as [`wavespectra`](https://github.com/wavespectra/wavespectra) and [`OceanWaves`](https://github.com/openearth/oceanwaves-python), provide an abstraction layer for typical directional wave spectra manipulation. These packages provide functionalities for common operations, such as reading and writing spectral data from numerical models (e.g., SWAN, WAVEWATCH III, etc.) and observational platforms (e.g., buoys, ADCP, etc.), spectral partitioning, estimation of spectral parameters, standard conversions and advanced plotting. They do not incorporate, however, algorithms for estimating the directional wave spectrum from raw wave measurements. EWDM complemenents the ecosystems of tools, providing routines for estimating the directional wave spectrum directly from measururements of wave time series. The results obtained from EWDM can potentially
be integrated with [`wavespectra`](https://github.com/wavespectra/wavespectra) or [`OceanWaves`](https://github.com/openearth/oceanwaves-python) to broaden their functionalities.


# Acknowledgements

This work was funded by the European Research Council (ERC) under the EU Horizon 2020 research and innovation program (Grant Agreement 833125- HIGHWAVE) and by the Centre National d'Études Spatiales (CNES) under the project ARANSAT. DPZ would like to thank the maintainers and developers of the open-source scientific Python ecosystem for making data processing and visualization more efficient. DPZ is grateful to the Navier (ENS Paris-Saclay) and Stokes (University College Dublin) teams for the interesting discussions.


# References
