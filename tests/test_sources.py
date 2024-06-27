#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
@author: Daniel Santiago
@github: http://github.com/dspelaez
@created: 2024-06-27
"""

import numpy as np
import xarray as xr
import unittest

from datetime import datetime

from .sources import SpotterBuoysDataSource, CDIPDataSourceRealTime

class TestSpotterBuoysDataSource(unittest.TestCase):

    def test_init(self):
        # testing file existence
        with self.assertRaises(FileNotFoundError):
            SpotterBuoysDataSource(fname="nonexistent.csv")

    def test_read_dataset(self):
        # testing return of xarray object
        spotter = SpotterBuoysDataSource(fname="displacement.csv")
        ds = spotter.read_dataset()
        self.assertIsNotNone(ds)


class TestCDIPDataSourceRealTime(unittest.TestCase):

    def test_read_dataset(self):
        # testing return of xr.Dataset
        cdip = CDIPDataSourceRealTime(station_id="189")
        dataset = cdip.read_dataset(
            time_start="2022-01-01 00:00:00", time_end="2022-01-01 01:00:00"
        )
        self.assertIsNotNone(dataset)

# if __name__ == '__main__':
    # unittest.main()
