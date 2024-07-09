Mathematical framework
======================

Continuous Wavelet Transform
----------------------------

The Continuous Wavelet Transform (CWT) of say a single-point sea surface
elevation signal :math:`\eta(t)`, is defined as:

.. math:: W(\tau, \sigma) = \int_t \eta(t) \psi(\tau-t,\sigma) dt

where :math:`\psi^*(\tau, \sigma)` is the core of the Wavelet transform and is known as the *Mother wavelet*, and :math:`\sigma` and :math:`\tau` are a scale and time factors, respectively. Using the convolution theorem, it is possible to obtain the complex wavelet coefficients in terms of the surface elevation spectrum as:

.. math:: W(\tau, \sigma) = \int_t \hat{\eta}(\omega) \hat{\psi}(\sigma \omega) e^{i\omega t} dt

In this case, :math:`\hat{\psi}(\sigma\omega)` represents the Fourier transform of the *Mother wavelet*. There are several functions used as Mother Wavelet, and in principle, the method can be used with any of them. Donelan proposed to use the Morlet wavelet because of its characterists and its Fourier transform can be expresed as:

.. math:: \hat{\psi}(\sigma\omega) = \left(\frac{2\pi\sigma}{\Delta t}\right)^{1/2} \pi^{-1/4} H(\omega) e^{(-\sigma\omega-\omega_0)/2}

where :math:`H(\omega)` is the step-function and :math:`\omega_0=6` to the admisiblity condition could be satisfied. For practical purposes, the wavelet scale and time factors, :math:`\sigma` and :math:`\tau`, can be replaced by the wave frequency and time scale, respectively.


Estimation of local wave direction
----------------------------------

Arrays of wave gauges
^^^^^^^^^^^^^^^^^^^^^

Assume we have a spatial array of sensors. In each point with coordinates :math:`\mathrm{x}_j=(x_j, y_j)` we have a time series of sea surface elevation :math:`\eta_j(t)`. The first step is to compute the Continuous Wavelet Transform (CWT) to each one of these time series. Let's also assume that we can express the surface elevation at each point as a sum of plane waves with different frequencies and directions, as:

.. math:: \eta_j(t) = \sum a \cos(\mathbf{k} \cdot \mathbf{x} - \omega t + \phi)

If the assumption this is valid, we can take a pair of points, located at :math:`\mathbf{x}_m` and :math:`\mathbf{x}_n`, in a spatial array, and express the surface elevation as:

.. math:: |W_m| e^{i(\mathbf{k}\cdot\mathbf{x}_m-\omega t + \phi_m)} = |W_n| e^{i(\mathbf{k}\cdot\mathbf{x}_n-\omega t + \phi_n)},

where :math:`|W_{m}|` and :math:`|W_{n}|` are the wave amplitude at the points :math:`m` and :math:`n`, respectively; :math:`\mathbf{k}=(k_x, k_y)` is the wavenumber; :math:`\omega` is the wave frequency, and :math:`\phi` is the wave phase at each point. If we consider that :math:`|W_m| = |W_n|` for a given wave component, then the phase difference between the two points is:

.. math:: \Delta\phi_{mn} = \mathbf{k}\cdot\mathbf{x}_m - \mathbf{k}\cdot\mathbf{x}_n

which can be expressed in a compact way, considering all the possible pairs in
the array, as:

.. math:: X \mathbf{k} = \Delta\phi

where,

.. math::
   X = \left(\begin{array}{cc}
        \Delta x_{12} & \Delta y_{12} \\
        \Delta x_{13} & \Delta y_{13} \\
        \vdots        & \vdots        \\
        \Delta x_{mn} & \Delta y_{mn} \\
        \end{array}\right),

is the pair separation matrix and, the corresponding phase difference is:

.. math::
   \Delta\phi = \left(\begin{array}{c}
        \Delta \phi_{12}  \\
        \Delta \phi_{13}  \\
        \vdots         \\
        \Delta \phi_{mn}  \\
        \end{array}\right)

The system can be solved using a least-squares estimator as:

.. math::
   \mathbf{k}^{\mathrm{LS}} = (X^{\mathrm{T}} X)^{-1} X^{\mathrm{T}} \Delta\phi

Using the least-square estimations, the root mean squared error of the wavenumber estimation can be expressed as :math:`\epsilon = ||X \mathbf{k}^{\mathrm{LS}} - \Delta\phi||^2`. As we seek pairs in the array to be approximately perpendicular to wave traveling direction, a restriction in the phase difference :math:`\Delta\phi` can be imposed, like :math:`\Delta\phi \le \alpha 2\pi`. A value of :math:`\alpha \approx 0.5` is set by default.

Finally, applying the wavelet transform to each time series in the array, allow one to compute the power :math:`|W_{mn}(f,t)|^2` and phase :math:`\phi_{ij}(f,t)` as function of the frequency (wavelet scale) and time. Then, if the linear system is solved, one can compute the wavenumber components :math:`k_x(f,t)` and :math:`k_y(f,t)`, and likewise the wave direction as,

.. math:: \theta(f,t) = \tan^{-1}\left( \frac{k_y}{k_x} \right).

One advantage arrays is that one can compute the wavenumber spectrum or even the full three-dimensional wave spectrum.


Single-point triplets
^^^^^^^^^^^^^^^^^^^^^

In the context of single-point triplet, such as in the case of GPS wave buoys
that deliver sea surface elevation and the two ortogoonal components of the
horizonal wave-induced velocities (or displacements), the procedure to obtain
the local wave direction is simplet. First, the cross-wavelet transform between
sea surface elevation and each velocity component is computed:

.. math:: W_{u\eta}(f,\theta) = W_u(f,\theta) W_\eta^*(f,\theta)

and

.. math:: W_{v\eta}(f,\theta) = W_v(f,\theta) W_\eta^*(f,\theta)

where the start represents complex conjugate.

Considering that surface elevation and wave orbital velocities are in phase, i.e., when water surface goes up, the horizonal velocity is positive, and vice-verse, so the cross-wavelet coefficients provide an idea on the amount of wave energy travelling eastwards and northward, resectively. Therefore, we compute the local wave direction simply as:

.. math:: \theta(f,t) = \tan^{-1} \left[ \frac{\mathcal{R} \{W_{v\eta}(f,\theta) \}}{\mathcal{R} \{W_{u\eta}(f,\theta) \}} \right].


Directional distribution and directional spectrum
-------------------------------------------------

It is a common practice to express the directional spectrum as a product the frequency spectrum and the directional distribution function, i.e.:

.. math:: E(f, \theta) = S(f) D(f, \theta)

where :math:`S(f) = \int E(f,\theta) d\theta` is the azimutally-integrated wave spectra. In the context of CWT, it is simply the time-integrated wavelet power:

.. math:: S(f) = \frac{1}{T} \int |W(f,t)| ^ 2 dt

where :math:`T` is the time-series length. The directional distribution function must satisfy:

.. math:: \int D(f, \theta) d\theta = 1

This can be seen as a probability distribution of wave directions, hence it can be estimated counting the occurrences of a given direction along the time, per each frequency. Alternatively, it can be estimated -- as implemented here -- using KDE (Kernel Density Estimation) to achieve a better balance between accuracy and smootness. As we are dealing here with angular data, we use a kernel function based on von Karman distribution. See Appendix in `Pelaez-Zapata et al. (2024)`_ for further details.

.. _Pelaez-Zapata et al. (2024): https://doi.org/10.1175/JTECH-D-23-0058.1
