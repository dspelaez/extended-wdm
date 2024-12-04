---
title: 'EWDM: A wavelet-based method for estimating directional spectra of ocean
waves'
tags:
  - python
  - physical oceanography
  - ocean surface waves
  - directional wave spectra
  - wavelets
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
date: 1 December 2024
bibliography: paper.bib

---

# Summary

The EWDM (Extended Wavelet Directional Method) is a Python toolkit developed to accurately estimate the directional spectra of ocean waves using wavelet-based techniques. It aims to provide a robust estimation of the directional wave spectrum from diverse sources of data, including GPS buoys, pitch-roll buoys, arrays of wave staffs, and ADCPs. This research tool addresses the limitations of conventional techniques in capturing accurate directional wave information. Furthermore, it extends the capabilities of the original WDM [@Donelan_etal_1996], which was previously only able to resolve spatial arrays of wave staffs, to include other common sources of single-point data such as ADCPs and wave buoys [@PelaezZapata_2024a; @Krogstad_etal_2006]. Key features of the EWDM include the implementation of wavelet-based algorithms for extracting directional information from wave time series, improved estimation of wave directional distribution using Kernel Density Estimation (KDE), tools for processing and visualizing directional wave data, and compatibility with popular data sources such as SOFAR Spotter buoys and the CDIP database. The package is powered by `xarray` [@Hoyer_2017] labelled multi-dimensional arrays, enhancing its efficiency and scalability.


# Statement of need

Various aspects of physical oceanography, including air-sea interactions, wave-induced mixing and upper layer dynamics, as well as several applications in coastal engineering, such as, wave forecasting, prediction of rogue waves, quantification of coastal erosion, and the design of offshore structures and renewable energy devices, require accurate and precise knowledge of the directional distribution of ocean waves [@Barstow_etal_2005]. The directional wave spectrum provide essential information about the distribution of wave energy with respect to frequency and direction. However, the dynamic and unpredictable nature of ocean waves, coupled with the constraints of current measurement technologies and analysis methods, pose challenges in attaining accurate directional information.

The directional wave spectra can be estimated from measurements using either spatial arrays of wave gauges or single-point triplets, where three perpendicular wave quantities are measured at the same point [@Barstow_etal_2005]. This latter is for instance the case of wave buoys, which are the most wide-spread means for obtaining _in-situ_ wave data. However, current methods rely either on truncated Fourier series -that are unable to capture complex sea states- or on statistical fitting of the distribution of wave directions and often involve assumptions that may not hold in real-world conditions. See @PelaezZapata_2024a, @Benoit_1993, @Benoit_1993a, @Benoit_and_Teisson_1995 and @Benoit_etal_1997 for a comprehensive review of various methods.

@Donelan_etal_1996 introduced an alternative method based on the continuous wavelet transform (CWT) to overcome some of the limitations of classic methods, though specifically applicable to spatial arrays of wave gauges. Building on this, @PelaezZapata_2024a and @Krogstad_etal_2006 demonstrated the extension of this wavelet-based approach to single-point measurements, offering a promising advancement in directional wave spectrum estimation. This package is therefore an implementation of these algorithms. It has been successfully employed in the recent research by @PelaezZapata_2024a, delivering robust and accurate results. This application further demonstrates its effectiveness in addressing the challenges posed by real-world complex sea-state conditions. It is expected to be a valuable tool for oceanographers, engineers, students, and practitioners working in the field of ocean wave analysis.


# Acknowledgements

This work was funded by the European Research Council (ERC) under the EU Horizon 2020 research and innovation program (Grant Agreement 833125- HIGHWAVE) and by the Centre National d'Études Spatiales (CNES) under the project ARANSAT. DPZ would like to thank the maintainers and developers of open-source scientific Python ecosystem for making data processing and visualization more efficient. DPZ is grateful to the Navier (ENS Paris-Saclay) and Stokes (University College Dublin) teams for the fruitful discussions during the analysis of the data.


# References
