.. extended-wdm documentation master file, created by
   sphinx-quickstart on Thu Jun 27 21:31:51 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/logo.png
   :width: 200 px
   :align: right

Welcome to extended-wdm package documentation!
==============================================

*EWDM (Extended Wavelet Directional Method) is a Python toolkit for a
wavelet-based estimation of the directional wave spectra*

This algorithm enhances the capabilities of `Donelan's`_ WDM (Wavelet
Directional Method), which is specifically designed to calculate directional
spectra from spatial arrays, extending it to single-point triplet data from
sources like GPS wave buoys, pitch-roll buoys, and Acoustic Doppler Current
Profilers.

The wavelet-based methods for estimating the directional wave spectra have
emerged as a practical alternative to the conventional Fourier-based techniques,
particularly well-suited for the analysis of data from `SOFAR Spotter buoys`_
as demonstrated by `Pelaez-Zapata et al. (2024)`_ and from other triplet data
such as wave slopes as shown by `Krogstad et al. (2006)`_.

Key features of the EWDM include:

* Implementation of the wavelet-based algorithms for extracting directional
  information from wave time series.

* Improve estimation of wave directional distribution using KDE (Kernel Density
  Estimation)

* Tools for processing and visualising directional wave data.

* Powered by `xarray` labelled multi-dimensional arrays.

* Helper functions to handle commonly used data sources such as `SOFAR Spotter
  buoys`_ and `CDIP database`_.

* Documentation, examples, and comparison with conventional methods.

.. _Donelan's: https://doi.org/10.1175/1520-0485(1996)026<1901:naotdp>2.0.co;2
.. _Pelaez-Zapata et al. (2024): https://doi.org/10.1175/JTECH-D-23-0058.1
.. _Krogstad et al. (2006): https://onepetro.org/IJOPE/article-abstract/28936/Wavelet-And-Local-Directional-Analysis-of-Ocean?redirectedFrom=fulltext
.. _SOFAR Spotter buoys: https://www.sofarocean.com/products/spotter
.. _CDIP database: https://cdip.ucsd.edu/


Whether you are a researcher, student, or engineer in physical oceanography, the
Extended Wavelet Directional Method provides a powerful, user-friendly toolkit
for in-depth analysis of directional ocean wave spectra. Join us in exploring
the fascinating world of directional wave analysis and making meaningful
contributions to the understanding of ocean wave dynamics. We welcome
contributions, feedback, and collaboration from the community to further enhance
the capabilities of the Extended Wavelet Directional Method.


.. note::
   This is an experimental version.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   maths
   examples


API Documentation
=================

.. toctree::
   :maxdepth: 4

   ewdm


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
