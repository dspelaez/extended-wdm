#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

# physical parameters
GRAV = 9.8

# variables metadata
VARIABLE_NAMES = {
     'time': {
         'standard_name': 'time',
         'long_name': 'Time',
         'units': 'seconds since 1970-01-01 00:00:00'
     },
     'frequency': {
         'standard_name': 'sea_surface_wave_frequency',
         'long_name': 'Wave frequency',
         'units': 'Hz'
     },
     'direction': {
         'standard_name': 'sea_surface_wave_direction',
         'long_name': 'Wave direction',
         'units': 'degrees'
     },
     'frequency_spectrum': {
         'standard_name': 'sea_surface_frequency_wave_spectrum',
         'long_name': 'Azimuthally-integrated wave energy density spectrum',
         'units': 'm^2/Hz'
     },
     'directional_spectrum': {
         'standard_name': 'sea_surface_directional_wave_spectrum',
         'long_name': 'Directional wave energy density spectrum',
         'units': 'm^2/Hz/degrees'
     },
     'directional_distribution': {
         'standard_name': 'sea_surface_wave_directional_distribution_function',
         'long_name': 'Directional distribution function of wave energy',
         'units': 'dimensionless'
     },
     'wavenumber': {
         'standard_name': 'sea_surface_wave_wavenumber',
         'long_name': 'Wavenumber magnitude',
         'units': 'rad/m'
     },
     'nu': {
         'standard_name': 'sea_surface_wave_inverse_phase_speed',
         'long_name': 'Inverse phase speed (nu = k / omega)',
         'units': 's/m'
     },
     'wavenumber_spectrum': {
         'standard_name': 'sea_surface_wavenumber_wave_spectrum',
         'long_name': 'Omnidirectional wavenumber wave spectrum F(k)',
         'units': 'm^3'
     },
     'nu_spectrum': {
         'standard_name': 'sea_surface_inverse_phase_speed_wave_spectrum',
         'long_name': 'Omnidirectional inverse-phase-speed wave spectrum Q(nu)',
         'units': 'm^3/s'
     },
     'directional_spectrum_f': {
         'standard_name': 'sea_surface_directional_wave_spectrum',
         'long_name': 'Directional wave energy density spectrum F(f, theta)',
         'units': 'm^2/Hz/degrees'
     },
     'directional_spectrum_k': {
         'standard_name': 'sea_surface_wavenumber_directional_wave_spectrum',
         'long_name': 'Directional wavenumber wave spectrum Psi(k, theta)',
         'units': 'm^4/rad'
     },
     'directional_spectrum_nu': {
         'standard_name': 'sea_surface_inverse_phase_speed_directional_wave_spectrum',
         'long_name': 'Directional inverse-phase-speed wave spectrum Q(nu, theta)',
         'units': 'm^4/(s^2 rad)'
     },
     'mean_direction': {
         'standard_name': 'sea_surface_wave_mean_direction',
         'long_name': 'Energy-weighted mean wave direction',
         'units': 'degrees'
     },
     'directional_spread': {
         'standard_name': 'sea_surface_wave_directional_spread',
         'long_name': 'Energy-weighted directional spread',
         'units': 'degrees'
     },
     'lh_direction': {
         'standard_name': 'sea_surface_wave_longuet_higgins_direction',
         'long_name': 'Longuet-Higgins first-moment mean direction',
         'units': 'degrees'
     },
     'sign_reference': {
         'standard_name': 'sea_surface_wave_sign_reference_direction',
         'long_name': 'Unwrapped sign-reference direction',
         'units': 'degrees'
     },
     'var_eta': {
         'standard_name': 'sea_surface_elevation_variance',
         'long_name': 'Sea surface elevation variance (m0)',
         'units': 'm^2'
     },
    'surface_elevation': {
        'standard_name': 'sea_surface_wave_elevation',
        'long_name': 'Sea surface wave elevation',
        'units': 'm'
    },
    'eastward_displacement': {
        'standard_name': 'sea_surface_wave_eastward_displacement',
        'long_name': 'Wave-induced eastward displacement',
        'units': 'm'
    },
    'northward_displacement': {
        'standard_name': 'sea_surface_wave_northward_displacement',
        'long_name': 'Wave-induced buoy displacement',
        'units': 'm'
    },
    'eastward_velocity': {
        'standard_name': 'sea_surface_wave_eastward_velocity',
        'long_name': 'Wave-induced eastward velocity',
        'units': 'm/s'
    },
    'northward_velocity': {
        'standard_name': 'sea_surface_wave_northward_velocity',
        'long_name': 'Wave-induced northward velocity',
        'units': 'm/s'
    },
    'eastward_acceleration': {
        'standard_name': 'sea_surface_wave_eastward_acceleration',
        'long_name': 'Wave-induced eastward acceleration (surge)',
        'units': 'm/s^2'
    },
    'northward_acceleration': {
        'standard_name': 'sea_surface_wave_northward_acceleration',
        'long_name': 'Wave-induced northward acceleration (sway)',
        'units': 'm/s^2'
    },
    'upward_acceleration': {
        'standard_name': 'sea_surface_wave_upward_acceleration',
        'long_name': 'Wave-induced upward acceleration (heave)',
        'units': 'm/s^2'
    },
    'eastward_slope': {
        'standard_name': 'sea_surface_wave_eastward_slope',
        'long_name': 'Wave-induced eastward slope (roll)',
        'units': 'dimensionless'
    },
    'northward_slope': {
        'standard_name': 'sea_surface_wave_northward_slope',
        'long_name': 'Wave-induced northward slope (pitch)',
        'units': 'dimensionless'
    },
    'eastward_gyro': {
        'standard_name': 'sea_surface_wave_eastward_gyro',
        'long_name': 'Wave-induced eastward gyro',
        'units': 's^-1'
    },
    'northward_gyro': {
        'standard_name': 'sea_surface_wave_northward_gyro',
        'long_name': 'Wave-induced northward gyro',
        'units': 's^-1'
    },
    'upward_gyro': {
        'standard_name': 'sea_surface_wave_up_gyro',
        'long_name': 'Wave-induced upward gyro',
        'units': 's^-1'
    },
    'position_x': {
        'standard_name': 'position_x_wave_array',
        'long_name': 'Coordinate x of wave array element (eastward)',
        'units': 'm'
    },
    'position_y': {
        'standard_name': 'position_y_wave_array',
        'long_name': 'Coordinate y of wave array element (northward)',
        'units': 'm'
    }
}
