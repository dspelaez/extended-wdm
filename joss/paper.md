---
title: 'EWDM: Extended Wavelet Directional Method for estimating ocean waves directional spectra'
tags:
  - Python
  - oceanography
authors:
  - name:
        firstname: Daniel
        surname:  Peláez-Zapata
    orcid: 0000-0001-5862-6194
    affiliation: 1
  - name:
        firstname: Frédéric
        surname: Dias
  - orcid: 0000-0002-5123-4929
    affiliation: "1, 2"
affiliations:
 - name: École Normale Supérieure Paris-Saclay, France
   index: 1
   ror: 02hcn4061
 - name: University College Dublin, Ireland
   index: 2
date: 1 December 2024
bibliography: paper.bib

---

# Summary

The research purpose of the Extended Wavelet Directional Method (EWDM) software is to address the limitations of conventional Fourier-based techniques in accurately capturing directional wave information in physical oceanography. It aims to provide a robust estimation of the directional wave spectrum from diverse sources of data such as GPS buoys, pitch-roll buoys, arrays of wave staffs, and ADCPs. The software implements wavelet-based algorithms, Kernel Density Estimation (KDE), and tools for processing and visualizing directional wave data, making it suitable for researchers, students, and engineers in physical oceanography.

In the context of related work, the EWDM builds upon the traditional Wavelet Directional Method (WDM) by extending its capabilities to incorporate a wide range of data sources and configurations, thus offering an alternative and improved methodology for estimating directional wave spectra. Furthermore, the collaborative and open approach of the EWDM welcomes contributions, feedback, and collaboration from the community, aligning with the principles of transparency, reproducibility, and accessibility within the physical oceanography research community.

Considering the significance of improved directional wave spectrum estimation and the unique features offered by the EWDM, its publication as an open-source software package in the Journal of Open Source Software (JOSS) would not only make the method more accessible but also promote transparency, reproducibility, and collaboration within the field of physical oceanography research.


# Statement of need

The estimation of directional wave spectra is paramount in the field of physical oceanography for understanding wave dynamics and coastal processes. Conventional Fourier-based techniques have inherent limitations in accurately capturing directional wave information, prompting the need for alternative methodologies. The wavelet-based approach has emerged as a promising alternative, offering improved directional wave spectrum estimation.

The Extended Wavelet Directional Method (EWDM) is a Python toolkit designed to address the shortcomings of conventional Fourier-based techniques and provide a robust estimation of the directional wave spectrum from diverse sources of data, including GPS buoys, pitch-roll buoys, arrays of wave staffs, and ADCPs. With specific implementations for spatial arrays of wave staffs inspired by Donelan's WDM, as well as methods for single-point triplet data drawn from [@PelaezZapata_2024a] and Krogstad et al. (2006), the EWDM extends the capabilities of the original WDM to incorporate a wide range of data sources and configurations.

Key features of the EWDM include the implementation of wavelet-based algorithms for extracting directional information from wave time series, improved estimation of wave directional distribution using Kernel Density Estimation (KDE), tools for processing and visualizing directional wave data, and compatibility with popular data sources such as SOFAR Spotter buoys and the CDIP database. The package is powered by xarray labelled multi-dimensional arrays, enhancing its efficiency and scalability.

The utility of the EWDM extends to researchers, students, and engineers in physical oceanography, offering a powerful, user-friendly toolkit for estimating directional wave spectra. By fostering an open and collaborative approach, the EWDM welcomes contributions, feedback, and collaboration from the community to further enhance its capabilities and facilitate advancements in the field of directional wave spectrum estimation.

In light of the significance of improved directional wave spectrum estimation and the unique features offered by the EWDM, its publication as an open-source software package in JOSS would not only enhance the accessibility of the method but also promote transparency, reproducibility, and collaboration within the physical oceanography research community.


# Acknowledgements


# References
