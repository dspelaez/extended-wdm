#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import numpy as np
import pandas as pd
import xarray as xr
import os

from pathlib import Path
from typing import Union, Tuple

from .helpers import get_sampling_frequency
from .parameters import VARIABLE_NAMES


class SpotterBuoysDataSource(object):
    """Class to handle Spotter buoy data parsed from SD Card"""

    def __init__(
            self,
            fname: str ="displacement.csv",
        ):
        if os.path.isfile(fname):
            self.fname = Path(fname)
        else:
            raise FileNotFoundError(
                "File not found. The name of the file parsed from the "
                "SD Card is likely `displacement.csv`. Please provide "
                "the full path to this file."
            )


    def read_dataset(self):
        """Read raw data file and return xarray object."""

        # define data columns
        variables = {
            k: VARIABLE_NAMES[k]
            for k in [
                "eastward_displacement",
                "northward_displacement",
                "surface_elevation"
                ]
            }
        data_cols = [col for col in variables]
        columns = ['year','month','day','hour','min','sec','msec'] + data_cols
        
        # load data (deprecated)
        # fmt = "%Y %m %d %H %M %S %f"
        # parser = lambda t: pd.Timestamp(
            # *map(int,t.split()[:-1]), int(t.split()[-1])*1000
        # )
        # data = pd.read_csv(
            # self.fname, names=columns, header=None, skiprows=1,
            # parse_dates={"time": columns[:7]}, date_parser=parser
        # )
        # read the CSV file without parsing date
        data = pd.read_csv(
            self.fname, names=columns, header=None, skiprows=1
        )
        # combine the desired date-time columns and convert them to datetime
        data["time"] = pd.to_datetime(
            data[columns[:7]].astype(str).agg(','.join, axis=1),
            format="%Y,%m,%d,%H,%M,%S,%f"
        )
        data = data.drop(columns=columns[:7])

        # convert to xarray and add metadata
        ds = data.set_index("time").to_xarray()
        for name, attrs in variables.items():
            ds[name].attrs = attrs

        ds["time"].encoding = {
            'units': 'days since 1970-01-01 00:00:00',
            'calendar': 'proleptic_gregorian'
        }
        ds.attrs["sampling_rate"] = get_sampling_frequency(
            ds.time, precision=3
        )

        return ds


class CDIPDataSourceRealTime(object):
    """Class to handle CDIP real-time buoy data"""

    def __init__(self, station_id: Union[str, int]):
        self.station_id = station_id
        self.baseurl = (
            "https://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime"
        )

    def read_dataset(
            self,
            time_start: Union[str, np.datetime64] = None,
            time_end: Union[str, np.datetime64] = None,
        ) -> xr.Dataset:
        """Read dataset from a specified time range and location.

           Parameters:
            - time_start (str or np.datetime64): Start time of the range
              to extract from the dataset. Format: 'YYYY-MM-DD HH:MM:SS'.
            - time_end (str or np.datetime64): End time of the range to 
              extract from the dataset. Format 'YYYY-MM-DD HH:MM:SS'.

         Returns:
         - ds (xr.Dataset): Dataset containing the extracted data.
     """

        url = f"{self.baseurl}/{self.station_id}p1_xy.nc"
        ds = xr.open_dataset(url)
    
        # sampling frequency
        fs = np.round(ds.xyzSampleRate.item(), 4)

        # create time array
        # filter_delay = np.timedelta64(int(1E9 * ds.xyzFilterDelay.item()), 'ns')
        time_delta = np.timedelta64(int(1E9 / fs), 'ns')
        miliseconds = np.arange(len(ds.xyzCount)) * time_delta
        time = ds.xyzStartTime.values + miliseconds
        self.time = time
        self.time_coverage_start = ds.time_coverage_start
        self.time_coverage_end = ds.time_coverage_end

        # check if time_start is passed otherwise use first time 
        if time_start is None:
            time_start = time[0]
        else:
            try:
                time_start = np.datetime64(time_start)
            except TypeError:
                raise Exception(
                    "`time_end` should be either string yyy-mm-dd HH:MM:SS "
                    "or datetime64 instance."
                )

        # if time_end is not defined, add 30 minutes to start time
        if time_end is None:
            time_end = time_start + np.timedelta64(30, "m")
        else:
            try:
                time_end = np.datetime64(time_end)
            except TypeError:
                raise Exception(
                    "`time_end` should be either string yyy-mm-dd HH:MM:SS "
                    "or datetime64 instance."
                )

        # compute indices
        try:
            ij = slice(
                np.where(time >= time_start)[0][0],
                np.where(time >= time_end)[0][0]
            )
        except IndexError:
            raise Exception(
                f"`time_start` and `time_end` couldn't be found in the dataset. "
                f"This dataset covers from {self.time_coverage_start} to "
                f"{self.time_coverage_end}."
            )


        # originally X is northward and Y westward, in my convection
        # X should be eastward and Y northward
        _x = -ds.xyzYDisplacement.isel(xyzCount=ij).data
        _y =  ds.xyzXDisplacement.isel(xyzCount=ij).data
        _z =  ds.xyzZDisplacement.isel(xyzCount=ij).data
        
        ds = xr.Dataset(
                coords={"time": time[ij]},
                data_vars={
                    "eastward_displacement": ("time", _x),
                    "northward_displacement": ("time", _y),
                    "surface_elevation": ("time", _z)
                },
                attrs={
                    "sampling_rate": fs,
                    **ds.attrs
                }
        )
        variables = {
            k: VARIABLE_NAMES[k]
            for k in [
                "eastward_displacement",
                "northward_displacement",
                "surface_elevation"
                ]
            }
        for name, attrs in variables.items():
            ds[name].attrs = attrs

        ds["time"].encoding = {
            'units': 'days since 1970-01-01 00:00:00',
            'calendar': 'proleptic_gregorian'
        }
        return ds






