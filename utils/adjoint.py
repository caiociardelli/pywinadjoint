#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
 PyWinAdjoint

 Author: Caio Ciardelli, University of São Paulo, May 2021

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

-----------------------------------------------------------------------------------------------

 ADJOINT

 USAGE
   ./utils/adjoint.py FIRST_EVENT LAST_EVENT INPUT_DIRECTORY MIN_PERIOD MAX_PERIOD OUTPUT_DIRECTORY

 EXAMPLE
   ./utils/adjoint.py 1 4 events/ 30 60 adjoint/

 COMMAND LINE ARGUMENTS
   FIRST_EVENT            - index of the first event
   LAST_EVENT             - index of the last event
   INPUT_DIRECTORY        - directory containing the input files
   MIN_PERIOD             - minimum period of the bandpass filter
   MIN_PERIOD             - maximum period of the bandpass filter
   OUTPUT_DIRECTORY       - directory where the routine will write the output files

 DESCRIPTION
   Reads observed and synthetic seismograms from the INPUT_DIRECTORY directory and computes the
   exponentiated phase adjoint sources, saved to the OUTPUT_DIRECTORY.

-----------------------------------------------------------------------------------------------
"""
import os
import sys
import ctypes as cp
import numpy as np
import math

from glob import glob
from obspy.core import read
from obspy.signal.filter import bandpass
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import ticker
from multiprocessing import Process, Lock


libsignal = cp.cdll.LoadLibrary (os.path.dirname (os.path.realpath (__file__)) + '/../lib/libsignal.so')
libdistance = cp.cdll.LoadLibrary (os.path.dirname (os.path.realpath (__file__)) + '/../lib/libdistance.so')
libweights = cp.cdll.LoadLibrary (os.path.dirname (os.path.realpath (__file__)) + '/../lib/libweights.so')


libsignal.gaussianSmooth.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                     cp.c_double,
                                     cp.c_double,
                                     cp.c_int]
libsignal.gaussianSmooth.restype = None

libsignal.getTracesAmplitude.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                         np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                         cp.c_int]
libsignal.getTracesAmplitude.restype = cp.c_double

libsignal.getAdjSourceAmplitude.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                            np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                            np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                            cp.c_int]
libsignal.getAdjSourceAmplitude.restype = cp.c_double

libsignal.taper.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                            cp.c_double,
                            cp.c_double,
                            cp.c_double,
                            cp.c_double,
                            cp.c_int]
libsignal.taper.restype = None

libsignal.rotateNE2RT.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                  np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                  cp.c_double,
                                  cp.c_int]
libsignal.rotateNE2RT.restype = None

libsignal.rotateRT2NE.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                  np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                  cp.c_double,
                                  cp.c_int]
libsignal.rotateRT2NE.restype = None

libsignal.adjointSourceWF.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      cp.c_double,
                                      cp.c_double,
                                      cp.POINTER (cp.c_double),
                                      cp.POINTER (cp.c_double),
                                      cp.c_int]
libsignal.adjointSourceWF.restype = None

libsignal.adjointSourceEP.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      cp.c_double,
                                      cp.c_double,
                                      cp.c_double,
                                      cp.POINTER (cp.c_double),
                                      cp.POINTER (cp.c_double),
                                      cp.c_int]
libsignal.adjointSourceEP.restype = None

libsignal.adjointSourceEV.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                      cp.c_double,
                                      cp.c_double,
                                      cp.c_double,
                                      cp.POINTER (cp.c_double),
                                      cp.POINTER (cp.c_double),
                                      cp.c_int]
libsignal.adjointSourceEV.restype = None

libsignal.writeBinary.argtypes = [cp.c_char_p,
                                  cp.c_char,
                                  np.ctypeslib.ndpointer (cp.c_float, flags = ['C', 'A', 'W', 'O']),
                                  cp.c_double,
                                  cp.c_double,
                                  cp.c_int]
libsignal.writeBinary.restype = None

libdistance.vincenty.argtypes = [cp.c_double,
                                 cp.c_double,
                                 cp.c_double,
                                 cp.c_double]
libdistance.vincenty.restype = cp.c_double

libweights.weights.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                               np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                               cp.c_double,
                               cp.c_double,
                               cp.c_uint,
                               cp.c_uint]
libweights.weights.restype = None


class Window (object):

    def __init__ (self, left = 0, right = 0, sampling_rate = None):

        self.left = left
        self.right = right
        self.sampling_rate = sampling_rate

    def __str__ (self):

        return 'Window object: ({}, {})'.format (self.left, self.right)

    def size (self):

        return self.right - self.left

    def width (self):

        if self.sampling_rate:

            return (self.right - self.left) / self.sampling_rate

        else:

            sys.exit ("Error: 'sampling_rate' not defined.")


class NoiseBar (object):

    def __init__ (self, left = 0, right = 0, ratio = 0, color = 'white', sampling_rate = None):

        self.left = left
        self.right = right
        self.ratio = ratio
        self.color = color
        self.sampling_rate = sampling_rate

    def __str__ (self):

        return 'NoiseBar object: ({}, {})'.format (self.left, self.right)

    def width (self):

        if self.sampling_rate:

            return (self.right - self.left) / self.sampling_rate

        else:

            sys.exit ("Error: 'sampling_rate' not defined.")


class Trace (object):

    def __init__ (self, trace, min_period, max_period):

        self.min_period = min_period
        self.max_period = max_period
        self.delta = trace.stats.sac.delta
        self.back_azimuth = trace.stats.sac.baz
        self.distance = trace.stats.sac.dist
        self.begin = trace.stats.sac.b
        self.end = trace.stats.sac.e
        self.npts = trace.stats.npts
        self.sampling_rate = trace.stats.sampling_rate
        self.station_latitude = trace.stats.sac.stla
        self.station_longitude = trace.stats.sac.stlo
        self.station_elevation = trace.stats.sac.stel
        self.data = np.require (trace.data, dtype = np.float64,
                                requirements = ['C', 'A', 'W', 'O', 'E'])


def getFilesNames (input_directory_name, event, kind):

    name = '{}/event{}/*.{}.sac'.format (input_directory_name, event, kind)

    return zip (*[sorted (glob (name))[i::3] for i in range (3)])


def getTracesAmplitude (obs, syn):

    return libsignal.getTracesAmplitude (obs.data, syn.data, obs.data.size)


def getAdjSourceAmplitude (adjoint):

    c1_id, c2_id, c3_id = adjoint.keys ()

    return libsignal.getAdjSourceAmplitude (adjoint[c1_id]['source'], adjoint[c2_id]['source'],
                                            adjoint[c3_id]['source'], adjoint[c1_id]['source'].size)


def taper (trace, window = 'hanning', width = 0.05):

    trace = np.require (trace, dtype = np.float64, requirements = ['C', 'A', 'W', 'O', 'E'])

    if window == 'hanning':

        c1 = 0.5
        c2 = 0.5

        omega = np.pi / (trace.size * width)

    elif window == 'hamming':

        c1 = 0.54
        c2 = 0.46

        omega = np.pi / (trace.size * width)

    elif window == 'cosine':

        c1 = 1
        c2 = 1

        omega = np.pi / (2 * trace.size * width)

    else:

        sys.exit ("Error: unknown window type '{}'.".format (window))

    libsignal.taper (trace, width, c1, c2, omega, trace.size)

    return trace


def rotateSystem (C1, C2, conversion, back_azimuth):

    if conversion == 'NE->RT':

        libsignal.rotateNE2RT (C1, C2, back_azimuth, C1.size)

    elif conversion == 'RT->NE':

        libsignal.rotateRT2NE (C1, C2, back_azimuth, C1.size)

    else:

        sys.exit ("Error: unknown type of coordinates system conversion '{}'.".format (conversion))


def adjointSource (obs, syn, tapered_windows, config):

    adjoint_type = config['adjoint_type']

    filter_adj_src = True if config['filter_adj_src'] == 'True' else False

    truncate    = float (config['truncate'])
    taper_kind  = config['taper_kind']
    taper_width = float (config['taper_width'])

    water_level = float (config['water_level'])

    tr1 = np.require (obs.data, dtype = np.float64, requirements = ['C', 'A', 'W', 'O', 'E'])
    tr2 = np.require (syn.data, dtype = np.float64, requirements = ['C', 'A', 'W', 'O', 'E'])

    windows = np.require (tapered_windows, dtype = np.float64, requirements = ['C', 'A', 'W', 'O', 'E'])

    adj_src = np.empty (tr1.size, dtype = np.float64)

    Xi   = cp.c_double ()
    Xi_w = cp.c_double ()

    sigma = obs.min_period * obs.sampling_rate / (2 * np.pi)

    libsignal.gaussianSmooth (tr1, sigma, truncate, tr1.size)
    libsignal.gaussianSmooth (tr2, sigma, truncate, tr2.size)

    if adjoint_type == 'waveform':

        libsignal.adjointSourceWF (tr1, tr2, windows, adj_src, obs.begin, obs.end,
                                   cp.byref (Xi), cp.byref (Xi_w), tr1.size)

    elif adjoint_type == 'exponentiated_phase':

        libsignal.adjointSourceEP (tr1, tr2, windows, adj_src, obs.begin, obs.end,
                                   water_level, cp.byref (Xi), cp.byref (Xi_w), tr1.size)

    elif adjoint_type == 'envelope':

        libsignal.adjointSourceEV (tr1, tr2, windows, adj_src, obs.begin, obs.end,
                                   water_level, cp.byref (Xi), cp.byref (Xi_w), tr1.size)

    else:

        sys.exit ("Error: unknown adjoint source type '{}'.".format (adjoint_type))

    if filter_adj_src:

        adj_src = np.require (bandpass (adj_src, freqmin = 1.0 / obs.max_period, freqmax = 1.0 / obs.min_period,
                              df = obs.sampling_rate, corners = 2, zerophase = True), dtype = np.float64,
                              requirements = ['C', 'A', 'W', 'O', 'E'])

        libsignal.gaussianSmooth (adj_src, sigma, truncate, obs.data.size)

    return taper (adj_src, taper_kind, taper_width), Xi.value, Xi_w.value


def readWindows (min_period, max_period, event):

    windows = dict ()

    with open ('windows_{:.0f}_{:.0f}/windows_{}.txt'.format (min_period, max_period, event), 'r') as WINDOWS:

        lines = WINDOWS.readlines ()

        sampling_rate = float (lines[0])

        for index in range (1, len (lines)):

            data = lines[index].split ('.')

            if len (data) == 4:

                network, station, channel = data[0], data[1], data[3][:-1]

                if not network in windows: windows[network] = dict ()
                if not station in windows[network]: windows[network][station] = dict ()
                if not channel in windows[network][station]: windows[network][station][channel] = list ()

                n = int (lines[index + 1])

                if n > 0:

                    for line in lines[index + 2:index + n + 2]:

                        data = line.split ()

                        left = int (data[0])
                        right = int (data[1])

                        windows[network][station][channel] += [Window (left, right, sampling_rate)]

    return windows


def taperedWindows (windows, length, config):

    taper_kind  = config['taper_kind']
    taper_width = float (config['taper_width'])

    tapered = dict ()

    for channel in windows:

        component = channel[-1]

        tapered[component] = np.zeros (length, dtype = np.float64)

        for window in windows[channel]:

            left = window.left
            right = window.right

            tapered[component][left:right] = taper (np.ones (window.size (), dtype = np.float64),
                                                    taper_kind, taper_width)

    return tapered


def writeBinary (name, axis, array):

    if sys.byteorder == 'little':

        endianess = '<'

    else:

        endianess = '>'

    libsignal.writeBinary (name.encode ('ascii'), endianess.encode ('ascii'),
                           array, axis[0], axis[-1], axis.size)


def writeAdjointSources (adjoint, axis, output_directory_name, network, station, index,
                         rotate, back_azimuth, config):

    save_binary = True if config['save_binary'] == 'True' else False

    if rotate:

        rotateSystem (adjoint['R']['source'], adjoint['T']['source'], 'RT->NE', back_azimuth)

    for component in adjoint:

        array = np.require (adjoint[component]['source'], dtype = 'float32', requirements = ['C', 'A', 'W', 'O', 'E'])

        if component == 'R':

            component = 'N'

        elif component == 'T':

            component = 'E'

        if save_binary:

            name = '{}/event{}/{}.{}.MX{}_adj.bin'.format (output_directory_name, index, network, station, component)

            writeBinary (name, axis, array)

        else:

            name = '{}/event{}/{}.{}.MX{}.adj'.format (output_directory_name, index, network, station, component)

            np.savetxt (name, np.c_[axis, array], fmt = '%1.6e %1.6e')


def receiversCoordinates ():

    receivers = dict ()

    with open ('STATIONS.txt', 'r') as STATIONS:

        for line in STATIONS:

            if len (line) > 5 and line[0:5] != 'Event':

                data = line.split ()

                name = data[0]
                latitude  = float (data[2])
                longitude = float (data[3])

                if name in receivers: continue

                receivers[name] = {'latitude' : latitude,
                                   'longitude' : longitude}

    return receivers


def vincenty (latitude1, longitude1, latitude2, longitude2):

    return libdistance.vincenty (latitude1, longitude1, latitude2, longitude2)


def computePairs (receivers):

    codes = sorted (receivers.keys ())
    n = len (codes)
    pairs = dict ()

    for i in range (n):

        for j in range (i + 1, n):

            key = '{}-{}'.format (codes[i], codes[j])

            pairs[key] = vincenty (receivers[codes[i]]['latitude'],
                                   receivers[codes[i]]['longitude'],
                                   receivers[codes[j]]['latitude'],
                                   receivers[codes[j]]['longitude'])

    return pairs


def getCodes (input_directory_name, event):

    codes = dict ()
    stations = list ()

    with open ('{}/event{}/output_selector.txt'.format (input_directory_name, event), 'r') as FILE:

        next (FILE)

        for line in FILE:

            data = line.split ()

            network, station, component, selected = data[0], data[1], data[2][-1], float (data[-1])

            if len (data[2]) != 3 or data[2][1] != 'H':

                component = data[3][-1]

            if not component in codes: codes[component] = list ()

            if selected > 0: codes[component] += [station]

            if not station in stations: stations += [station]

    for component in codes:

        codes[component] = sorted (codes[component])

    return codes, stations


def computeDistances (codes, pairs):

    n = len (codes)
    distances = np.zeros ((n, n))

    for i in range (n):

        for j in range (i + 1, n):

            key = '{}-{}'.format (codes[i], codes[j])

            distances[i, j] = pairs[key]
            distances[j, i] = distances[i, j]

    return distances


def computeWeights (codes, component, pairs):

    distances = computeDistances (codes[component], pairs)
    weights = np.zeros (len (codes[component]))
    libweights.weights (distances.flatten (), weights, 0.5, 0.35, 3, weights.size)

    return weights


def createAdjointSources (c1, c2, network, station, channel, component, tapered_windows,
                          weights, adjoint, output_adjoint, config):

    adj_src, Xi, Xi_w = adjointSource (c1, c2, tapered_windows[component], config)

    adj_src *= weights[station][component]
    Xi *= weights[station][component]
    Xi_w *= weights[station][component]

    adjoint[component] = {'source' : adj_src, 'misfit' : Xi_w}

    channel = channel[:-1] + component

    output_adjoint += '{:^7} {:^7} {:^8} {:>11.2f} {:>16.2f}\n'.format (network,
                                                                        station,
                                                                        channel,
                                                                        Xi,
                                                                        Xi_w)

    return adjoint, output_adjoint


def plotMisfits (ax, c1, c2, axis, name, windows, adjoint, component, adjoint_type):

    amplitude = getTracesAmplitude (c1, c2)

    ax.plot (axis, c1.data, color = 'black', label = 'Observed')
    ax.plot (axis, c2.data, color = 'red', label = 'Synthetic')
    ax.set_ylabel ('Amplitude (m)', fontsize = 15)
    ax.set_xlim (axis[0], axis[-1])
    ax.set_ylim (-amplitude, amplitude)

    network, station, location, channel = name.split ('.')

    channel += component

    for window in windows[channel]:

        left = axis[window.left]
        right = axis[window.right]

        rectangle = Rectangle ((left, -amplitude), window.width (),
                                2 * amplitude,
                                color = 'blue',
                                alpha = 0.25)
        ax.add_patch (rectangle)

    phase_misfit = adjoint[component]['misfit']

    if adjoint_type == 'waveform':

        title = '{}{}   Waveform Misfit: {:.2E}'.format (name,
                                                         component,
                                                         phase_misfit)

    elif adjoint_type == 'exponentiated_phase':

        title = '{}{}    Phase Misfit: {:.2E}'.format (name,
                                                       component,
                                                       phase_misfit)

    elif adjoint_type == 'envelope':

        title = '{}{}   Envelope Misfit: {:.2E}'.format (name,
                                                         component,
                                                         phase_misfit)

    ax.set_title (title, fontsize = 16, y = 1.01)

    ax.legend (loc = 1)

    ax.tick_params (axis = 'x', labelsize = 15)
    ax.tick_params (axis = 'y', labelsize = 15)
    ax.yaxis.set_major_locator (ticker.AutoLocator ())
    ax.yaxis.set_minor_locator (ticker.AutoMinorLocator ())
    ax.yaxis.set_major_formatter (ticker.ScalarFormatter (useMathText = True))
    ax.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))

    if component != 'Z': ax.label_outer ()
    else: ax.set_xlabel ('Time (s)', fontsize = 15)


def plotAdjointSources (ax, axis, amplitude, adjoint, component, adjoint_type):

    ax.plot (axis, adjoint[component]['source'], color = 'black', linewidth = 1.5)

    ax.set_xlim (axis[0], axis[-1])
    ax.set_ylim (-amplitude, amplitude)
    ax.grid ()

    if adjoint_type == 'waveform':

        title = 'Waveform Adjoint Source'

    elif adjoint_type == 'exponentiated_phase':

        title = 'Exponentiated Phase Adjoint Source'

    elif adjoint_type == 'envelope':

        title = 'Envelope Adjoint Source'

    ax.set_title (title, fontsize = 16, y = 1.01)

    ax.tick_params (axis = 'x', labelsize = 15)
    ax.tick_params (axis = 'y', labelsize = 15)
    ax.yaxis.set_major_locator (ticker.AutoLocator ())
    ax.yaxis.set_minor_locator (ticker.AutoMinorLocator ())
    ax.yaxis.set_major_formatter (ticker.ScalarFormatter (useMathText = True))
    ax.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))

    if component == 'Z': ax.set_xlabel ('Time (s)', fontsize = 15)
    else: plt.setp (ax.get_xticklabels (), visible = False)


def readConfig ():

    config = dict ()

    with open (os.path.dirname (os.path.realpath (__file__)) + '/../Config.cfg', 'r') as FILE:

        for line in FILE:

            if line[0] == '#':

                continue

            data = line.split ()

            config[data[0]] = data[2]

    return config


def main (lock, first, last, receivers, pairs, input_directory_name, output_directory_name, config):

    min_period = float (config['min_period'])
    max_period = float (config['max_period'])

    plot                  = True if config['plot'] == 'True' else False
    write_adjoint_sources = True if config['write_adjoint_sources'] == 'True' else False
    rotate                = True if config['rotate'] == 'True' else False
    receiver_weights      = True if config['receiver_weights'] == 'True' else False

    adjoint_type = config['adjoint_type']

    for event in range (first, last + 1, 4):

        Windows = readWindows (min_period, max_period, event)

        obs = getFilesNames (input_directory_name, event, 'd')
        syn = getFilesNames (input_directory_name, event, 's')

        stations_adjoint = str ()
        output_adjoint = 'Network Station Channel  Full Misfit  Windowed Misfit\n'

        codes, stations = getCodes (input_directory_name, event)

        weights = dict ()

        for station in stations:

            if rotate:

                weights[station] = {'R' : 1,
                                    'T' : 1,
                                    'Z' : 1}

            else:

                weights[station] = {'E' : 1,
                                    'N' : 1,
                                    'Z' : 1}

        if receiver_weights:

            for component in codes:

                wt = computeWeights (codes, component, pairs)

                for index, station in enumerate (codes[component]):

                    weights[station][component] *= wt[index]

        for obs, syn in zip (obs, syn):

            e1 = Trace (read (obs[0])[0], min_period = min_period, max_period = max_period)
            e2 = Trace (read (syn[0])[0], min_period = min_period, max_period = max_period)

            n1 = Trace (read (obs[1])[0], min_period = min_period, max_period = max_period)
            n2 = Trace (read (syn[1])[0], min_period = min_period, max_period = max_period)

            z1 = Trace (read (obs[2])[0], min_period = min_period, max_period = max_period)
            z2 = Trace (read (syn[2])[0], min_period = min_period, max_period = max_period)

            if rotate:

                rotateSystem (n1.data, e1.data, 'NE->RT', e1.back_azimuth)
                rotateSystem (n2.data, e2.data, 'NE->RT', e2.back_azimuth)

            string = obs[0].split ('/')[-1].split ('.')

            network, station, location, channel = string[0], string[1], string[2], string[3][:-1]

            name = '{}.{}.{}.{}'.format (network, station, location, channel)

            axis = np.linspace (z1.begin, z1.end, z1.npts)

            windows = Windows[network][station]

            tapered_windows = taperedWindows (windows, axis.size, config)

            adjoint = dict ()

            component = 'T' if rotate else 'E'

            adjoint, output_adjoint = createAdjointSources (e1, e2, network, station, channel, component,
                                                                    tapered_windows, weights, adjoint,
                                                                    output_adjoint, config)

            component = 'R' if rotate else 'N'

            adjoint, output_adjoint = createAdjointSources (n1, n2, network, station, channel, component,
                                                                    tapered_windows, weights, adjoint,
                                                                    output_adjoint, config)

            adjoint, output_adjoint = createAdjointSources (z1, z2, network, station, channel, 'Z',
                                                                    tapered_windows, weights, adjoint,
                                                                    output_adjoint, config)

            if plot:

                fig = plt.figure (figsize = (17, 8.5))

                amplitude = getAdjSourceAmplitude (adjoint)

                component = 'T' if rotate else 'E'

                ax1 = fig.add_subplot (321)

                plotMisfits (ax1, e1, e2, axis, name, windows, adjoint, component, adjoint_type)

                ax2 = fig.add_subplot (322)

                plotAdjointSources (ax2, axis, amplitude, adjoint, component, adjoint_type)

                component = 'R' if rotate else 'N'

                ax3 = fig.add_subplot (323)

                plotMisfits (ax3, n1, n2, axis, name, windows, adjoint, component, adjoint_type)

                ax4 = fig.add_subplot (324)

                plotAdjointSources (ax4, axis, amplitude, adjoint, component, adjoint_type)

                ax5 = fig.add_subplot (325)

                plotMisfits (ax5, z1, z2, axis, name, windows, adjoint, 'Z', adjoint_type)

                ax6 = fig.add_subplot (326)

                plotAdjointSources (ax6, axis, amplitude, adjoint, 'Z', adjoint_type)

                fig.subplots_adjust (left = 0.07, bottom = 0.08, top = 0.95, right = 0.97, wspace = 0.1, hspace = 0.18)

                lock.acquire ()

                print ('Saving {}_Adjoint_Source.pdf (event {})'.format (name, event))

                lock.release ()

                plt.savefig ('{}/event{}/{}_Adjoint_Source.pdf'.format (input_directory_name, event, name))

                plt.close (fig)

            if write_adjoint_sources:

                writeAdjointSources (adjoint, axis, output_directory_name, network,
                                     station, event, rotate, z1.back_azimuth, config)

            stations_adjoint += '{:5} {:2} {:10.5f} {:10.5f} {:7.1f} {:6.1f}\n'.format (station,
                                                                                        network,
                                                                                        z1.station_latitude,
                                                                                        z1.station_longitude,
                                                                                        z1.station_elevation,
                                                                                        0)

        with open ('{}/event{}/output_adjoint.txt'.format (input_directory_name, event), 'w') as FILE:

            print ("Saving 'output_adjoint.txt' (event {})".format (event))

            FILE.write (output_adjoint)

        with open ('{}/event{}/STATIONS_ADJOINT'.format (output_directory_name, event), 'w') as FILE:

            FILE.write (stations_adjoint)


def helpMenu ():

  help = """\n Error: wrong number of parameters on the comand line...

 ADJOINT

 USAGE
   ./utils/adjoint.py FIRST_EVENT LAST_EVENT INPUT_DIRECTORY MIN_PERIOD MAX_PERIOD OUTPUT_DIRECTORY

 EXAMPLE
   ./utils/adjoint.py 1 4 events/ 30 60 adjoint/

 COMMAND LINE ARGUMENTS
   FIRST_EVENT            - index of the first event
   LAST_EVENT             - index of the last event
   INPUT_DIRECTORY        - directory containing the input files
   MIN_PERIOD             - minimum period of the bandpass filter
   MIN_PERIOD             - maximum period of the bandpass filter
   OUTPUT_DIRECTORY       - directory where the routine will write the adjoint sources

 DESCRIPTION
   Reads observed and synthetic seismograms from the INPUT_DIRECTORY directory and computes the
   exponentiated phase adjoint sources, saved to the OUTPUT_DIRECTORY.\n"""

  print (help)


if __name__ == '__main__':

    if len (sys.argv) != 7:

        sys.exit (helpMenu ())

    config = readConfig ()

    lock = Lock ()

    receivers = receiversCoordinates ()
    pairs = computePairs (receivers)

    first = int (sys.argv[1])
    last  = int (sys.argv[2])

    input_directory_name = sys.argv[3]

    config['min_period'] = sys.argv[4]
    config['max_period'] = sys.argv[5]

    output_directory_name = sys.argv[6]

    Process (target = main, args = (lock, first + 0, last, receivers, pairs,
                                    input_directory_name, output_directory_name, config)).start ()
    Process (target = main, args = (lock, first + 1, last, receivers, pairs,
                                    input_directory_name, output_directory_name, config)).start ()
    Process (target = main, args = (lock, first + 2, last, receivers, pairs,
                                    input_directory_name, output_directory_name, config)).start ()
    Process (target = main, args = (lock, first + 3, last, receivers, pairs,
                                    input_directory_name, output_directory_name, config)).start ()

