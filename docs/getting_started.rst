Getting started
===============

Overview
--------

The current estimation techniques for the directional wave spectrum, such as the Maximum Likelihood (MLM) and Maximum Entropy Methods (MEM), rely on the assumption of stationary data due to their Fourier-based nature. `Donelan et al. (1996)`_ introduced the Wavelet Directional Method (WDM), which is a wavelet-based time-frequency decomposition approach to estimate the directional wave spectrum from surface elevation data. The WDM provides temporal information at the expense of reduced frequency resolution, but allows for detailed information on wave direction at a given frequency, even in non-stationary sea-states.

The deployment of spatial arrays in the field for wave measurements can be challenging. While wave buoys represent the vast majority of world-wide wave observations, they only provide information at a single point and are limited to three variables. `Krogstad et al. (2006)`_ proposed using analog wavelet-based first-five directional moments from an interpolated dataset of surface elevation and wave slopes to estimate the directional spectrum. Their results were compared with array-based WDM, showing good agreement. Furthermore, `Pelaez-Zapata et al. (2024)`_ proposed a more intuitive method using cross-wavelet coefficients of wave-induced velocities from a GPS wave buoy to estimate local wave direction and therefore directional wave spectrum. This approach is also applicable to wave-induced displacements and accelerations.

This pacakge can be considered as an extension of the original WDM (hence the name Extendend Wavelet Directional Method) as it can produce wavelet-based directional wave spectra, not only from spatial arrays, but also from single point triplets. Consequently, this package offers two main classes: :py:class:`ewdm.Arrays` and :py:class:`ewdm.Triplets` to handle the different kind of input data.


.. _Donelan et al. (1996): https://doi.org/10.1175/1520-0485(1996)026<1901:naotdp>2.0.co;2
.. _Pelaez-Zapata et al. (2024): https://doi.org/10.1175/JTECH-D-23-0058.1
.. _Krogstad et al. (2006): https://onepetro.org/IJOPE/article-abstract/28936/Wavelet-And-Local-Directional-Analysis-of-Ocean?redirectedFrom=fulltext


Installation
------------


Stable release
^^^^^^^^^^^^^^

To install extended-wdm, run this command in your terminal:

.. code-block:: console

    $ pip install ewdm

This is the preferred method to install extended-wdm, as it will always install
the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
^^^^^^^^^^^^

The sources for extended-wdm can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/dspelaez/extended-wdm

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/dspelaez/extended-wdm/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/dspelaez/extended-wdm
.. _tarball: https://github.com/dspelaez/extended-wdm/tarball/master


Conventions
-----------
