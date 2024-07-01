#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import xarray as xr
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.ticker as mticker

from typing import Union, Tuple


# plot directional wave spectrum {{{
def _smooth(E, ws=(5, 2)):
    """Applies a simple circular smoothing to 2D array.

    Args:
        E (ndarray): Input function.
        ws (tuple): Window size. For example, for a directional spectrum
            size E(360,64) and ws=(10,2), the filter acts averaging
            10 directions and 2 frequencies.

    Returns:
        ndarray: Function smoothed.
    """

    # define window
    nd, nf = ws
    if nf == nd:
        frqwin = signal.hamming(nf)
    else:
        frqwin = np.ones(nf)

    dirwin = signal.hamming(nd)
    window = frqwin[None,:] * dirwin[:,None]
    window = window / window.sum()

    # permorm convolution and return output
    return signal.convolve2d(E, window, mode='same', boundary='wrap')

def _get_axes(
        ax=None, rmin=0.1, rmax=0.5, rstep=0.1, angle=-135,
        color="0.8", as_period=False
    ):
    """Draw polar grid on specific axes.

    Args:
        ax (matplotlib.axes, optional): Axes object to draw the grid on.
            If None, a new figure and axes are created.
        rmin (float, optional): Minimum radius for the grid. Defaults to 0.1.
        rmax (float, optional): Maximum radius for the grid. Defaults to 0.5.
        rstep (float, optional): Step size for the radius. Defaults to 0.1.
        angle (int, optional): Angle for the radius labels. Defaults to -135.
        color (str, optional): Color of the grid lines and labels. Defaults to "0.8".
        as_period (bool, optional): If True, labels are formatted as wave periods.
            Defaults to False.

    Returns:
        tuple: A tuple containing the figure and axes objects if a new figure is created. Otherwise, returns None.
    """

    if ax is None:
        fig, ax = plt.subplots(1, figsize=(5,5))

    ax.set_aspect("equal")

    for radii in np.arange(rmin, rmax+rstep, rstep):
        circle = plt.Circle(
            (0,0), radii, color=color, alpha=0.5,
            linestyle="dashed", fill=False, zorder=2
        )
        ax.add_artist(circle)
        if radii <= rmax:
            if as_period:
                radii_label = f"{1/radii:.1f}"
            else:
                radii_label = f"{radii:.2f}"
            ax.text(
                radii*np.cos(np.radians(angle)),
                radii*np.sin(np.radians(angle)),
                radii_label, fontsize="small",
                ha="center", va="center", zorder=3
            )

    ax.axhline(0, color=color, ls="dashed", alpha=0.5)
    ax.axvline(0, color=color, ls="dashed", alpha=0.5)
    ax.plot([0,1], [0,1], "--", c=color, alpha=0.5, transform=ax.transAxes)
    ax.plot([0,1], [1,0], "--", c=color, alpha=0.5, transform=ax.transAxes)

    _label_args = {"fontsize": "small", "ha": "center", "va": "center"}
    ax.text(0.50, 0.95, "N", transform=ax.transAxes, **_label_args)
    ax.text(0.95, 0.50, "E", transform=ax.transAxes, **_label_args)
    ax.text(0.50, 0.05, "S", transform=ax.transAxes, **_label_args)
    ax.text(0.05, 0.50, "W", transform=ax.transAxes, **_label_args)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if as_period:
        ax.set_ylabel("Wave period [s]")
    else:
        ax.set_ylabel("$f$ [Hz]")

    ax.set_xlim([-rmax+rstep, rmax-rstep])
    ax.set_ylim([-rmax+rstep, rmax-rstep])

    try:
        return fig, ax
    except NameError:
        pass


def _get_cmap(colors=None, N=256):
    """Return colormap for a given color list.

    Args:
        colors (list): List of colours.
        N (int): Number of colours

    Returns:
        Colormap object. If not colors passed, it returns viridis with white bottom.

    """

    if colors is None:
        colors = [(1, 1, 1), *plt.cm.viridis(np.linspace(0, 1, 8))]
    return mcolors.LinearSegmentedColormap.from_list('cmap', colors, N=N)

def _add_cbar(
        pc, ax, cax=None, style="outside", ticks=None, orientation="horizontal",
        label="$\\log_{10} E \\; \\mathrm{[m^2 Hz^{-1} deg^{-1}]}$"
    ):
    """Return a colorbar object"""

    if cax is None:
        if style == "inside":
            ticks = mticker.LinearLocator(ticks)
            orientation = "horizontal"
            cax = ax.inset_axes(
                [0.07, 0.92, 0.3, 0.035], transform=ax.transAxes
            )
        if style == "outside":
            ticks = mticker.AutoLocator()
            orientation = "vertical"
            cax = ax.inset_axes(
                [1.02, 0.0, 0.04, 1.0], transform=ax.transAxes
            )
    else:
        pass

    return plt.colorbar(
        pc, cax=cax, orientation=orientation, ticks=ticks, label=label
    )

def _add_wind_info(ax, wspd, wdir, color="k", wind_sea_radius=True):
    """Add some info relative to wind speed"""

    fwind = 9.8 / (2*np.pi * wspd)
    uwnd, vwnd = (
        fwind * np.cos(np.radians(wdir)), fwind * np.sin(np.radians(wdir))
    )

    ax.arrow(
        0, 0, uwnd, vwnd, color=color, head_width=0.01,
        length_includes_head=True
    )

    if wind_sea_radius:
        circle = plt.Circle(
            (0,0), fwind, color=color, alpha=0.5,
            linestyle="dashed", fill=False
        )
        ax.add_artist(circle)


#def plot_directional_spectrum(
#        da: xr.DataArray, frqs="frequency", dirs="direction", ax=None,
#        smooth=None, cmap=None, levels=30, vmin=None, vmax=None, contours=None,
#        colorbar=False, wspd=None, wdir=None, wind_sea_radius=None,
#        curspd=None, curdir=None, cbar_kw={}, axes_kw={}
#    ):
#    """Make a simple plot of a direcitonal wave spectrum"""

def plot_directional_spectrum(
        da: Union[xr.DataArray, np.ndarray],
        frqs: Union[str, np.ndarray] = "frequency",
        dirs: Union[str, np.ndarray] = "direction",
        ax = None,
        smooth = None,
        cmap = None,
        levels: int = 30,
        vmin: Union[float, int] = None,
        vmax: Union[float, int] = None,
        contours: Union[float, int] = None,
        colorbar: bool = False,
        wspd: Union[float, int] =None,
        wdir: Union[float, int] =None,
        wind_sea_radius: Union[float, int] = None,
        curspd: Union[float, int] = None,
        curdir: Union[float, int] = None,
        cbar_kw={},
        axes_kw={}
    ) -> None:
    """
    Make a simple plot of a directional wave spectrum.

    Args:
        da: Directional spectrum data.
        frqs: Frequency label name or numpy array.
        dirs: Direction label name or numpy array.
        ax (optional): Matplotlib axis object.
        smooth (tuple): Smoothing factor for visualisation. Defaults to None.
        cmap (optional): Colormap for the plot. Defaults to None.
        levels (int): Number of contour levels. Defaults to 30. If None, pseudo-color plot is made.
        vmin (int or float): Minimum value for colormap.
        vmax (int or float]): Maximum value for colormap.
        contours (int or float): Specific contour levels to plot on top.
        colorbar (bool): Whether to display colorbar. Defaults to False.
        wspd (int or float): Wind speed value to draw arrow. Defaults to None.
        wdir (int or float): Wind direction value in cartesian convention.
        wind_sea_radius (int or float):Whether to display wind sea separation radius.
        curspd (int or float): Current speed value.
        curdir (int or float): Current direction value in cartesian convention.
        cbar_kw (dict): Additional arguments for colorbar.
        axes_kw (dict): Additional arguments for plot axes.

    Returns:
        fig, ax
    """
    # get axis if not given
    if ax is None:
        fig, ax = _get_axes(**axes_kw)
    else:
        _get_axes(ax=ax, **axes_kw)

    # wrap spectra around the circle
    if da[dirs][0] != da[dirs][-1]:
        padded = da.pad({dirs: (1,0)}, mode="wrap")
    else:
        padded = da.copy()

    # calculate cartesian fx,fy to mimic polar coordinates
    _deg2rad = np.pi/180
    frqx = padded[frqs].data[:,None]*np.cos(padded[dirs].data[None,:]*_deg2rad)
    frqy = padded[frqs].data[:,None]*np.sin(padded[dirs].data[None,:]*_deg2rad)

    # smooth spectra
    if smooth:
        smoothed = xr.apply_ufunc(lambda x: _smooth(x, ws=smooth), padded)
    else:
        smoothed = padded.copy()

    # get colormap
    if cmap is None:
        cmap = _get_cmap()
    else:
        if isinstance(cmap, list):
            cmap = _get_cmap(cmap)

    # plot pcolormesh if levels is None, otherwise go for contourf
    if levels is None:
        pc = ax.pcolormesh(
            frqx, frqy, smoothed[1:,1:], cmap=cmap,
            shading="flat", vmin=vmin, vmax=vmax
        )
    else:
        pc = ax.contourf(
            frqx, frqy, smoothed, levels=levels, cmap=cmap,
            vmin=vmin, vmax=vmax
        )

    # plot contour lines if contours is not None
    if contours is not None:
        cts = ax.contour(
            frqx, frqy, smoothed, levels=contours, colors="k"
        )

    # add colobar if True
    if colorbar:
        cbar = _add_cbar(pc, ax, **cbar_kw)

    # add wind data
    if (wspd is None) or (wspd <= 0):
        wspd = 10
    if wdir is not None:
        _add_wind_info(
            ax, wspd=wspd, wdir=wdir, color="k",
            wind_sea_radius=wind_sea_radius
        )

    # add current data
    if (curspd is None) or (curspd <= 0):
        curspd = 100
    if curdir is not None:
        _add_wind_info(
            ax, wspd=curspd, wdir=curdir, color="red",
            wind_sea_radius=None
        )

    try:
        return fig, ax, pc
    except NameError:
        pc

# }}}
