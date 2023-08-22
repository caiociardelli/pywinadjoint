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

 SELECTOR

 USAGE
   ./utils/selector.py FIRST_EVENT LAST_EVENT MIN_PERIOD MAX_PERIOD

 EXAMPLE
   ./utils/selector.py 1 4 events/ 30 60 adjoint/

 COMMAND LINE ARGUMENTS
   FIRST_EVENT            - index of the first event
   LAST_EVENT             - index of the last event
   MIN_PERIOD             - minimum period of the bandpass filter
   MIN_PERIOD             - maximum period of the bandpass filter

 DESCRIPTION
   Reads observed and synthetic seismograms from the 'events/' directory and computes the time-
   widows, saved to the 'windows/' directory.

-----------------------------------------------------------------------------------------------
"""
import os
import sys
import ctypes as cp
import numpy as np

from glob import glob
from obspy.core import read
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import ticker
from multiprocessing import Process, Lock


libsignal = cp.cdll.LoadLibrary (os.path.dirname (os.path.realpath (__file__)) + '/../lib/libsignal.so')
libdistance = cp.cdll.LoadLibrary (os.path.dirname (os.path.realpath (__file__)) + '/../lib/libdistance.so')

libsignal.gaussianSmooth.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                     cp.c_double,
                                     cp.c_double,
                                     cp.c_int]
libsignal.gaussianSmooth.restype = None

libsignal.getTracesAmplitude.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                         np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                         cp.c_int]
libsignal.getTracesAmplitude.restype = cp.c_double

libsignal.booleanZeros.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                   np.ctypeslib.ndpointer (cp.c_bool, flags = ['C', 'A', 'W', 'O']),
                                   cp.c_int,
                                   cp.c_int,
                                   cp.c_int]
libsignal.booleanZeros.restype = cp.c_int

libsignal.findEdges.argtypes = [np.ctypeslib.ndpointer (cp.c_bool, flags = ['C', 'A', 'W', 'O']),
                                np.ctypeslib.ndpointer (cp.c_int, flags = ['C', 'A', 'W', 'O']),
                                cp.c_int]
libsignal.findEdges.restype = None

libsignal.correlate.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                cp.c_int,
                                cp.c_int]
libsignal.correlate.restype = cp.c_double

libsignal.misfit.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                             np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                             cp.c_int]
libsignal.misfit.restype = cp.c_double

libsignal.similar.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                              np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                              cp.c_int,
                              cp.c_int,
                              cp.c_double,
                              cp.c_double]
libsignal.similar.restypes = cp.c_bool

libsignal.signalToNoiseRatio.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                         np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                         cp.c_int,
                                         cp.c_double,
                                         cp.c_double,
                                         cp.c_double,
                                         cp.POINTER (cp.c_int)]
libsignal.signalToNoiseRatio.restype = cp.c_int

libsignal.noiseLevels.argtypes = [np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                  np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                  cp.c_int,
                                  np.ctypeslib.ndpointer (cp.c_double, flags = ['C', 'A', 'W', 'O']),
                                  cp.c_int,
                                  cp.c_double,
                                  cp.c_double]
libsignal.noiseLevels.restype = cp.c_bool

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

libdistance.vincenty.argtypes = [cp.c_double,
                                 cp.c_double,
                                 cp.c_double,
                                 cp.c_double]
libdistance.vincenty.restype = cp.c_double


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

            return abs (self.right - self.left) / self.sampling_rate

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


def getFilesNames (min_period, max_period, event, kind):

    name = 'events_{:.0f}_{:.0f}/event{}/*.{}.sac'.format (min_period, max_period, event, kind)

    return zip (*[sorted (glob (name))[i::3] for i in range (3)])


def signalToNoiseRatio (obs, syn, config):

    min_initial_signal_to_noise_ratio = float (config['min_initial_signal_to_noise_ratio'])
    min_signal_to_noise_ratio         = float (config['min_signal_to_noise_ratio'])
    min_quiet_window_length_ratio     = float (config['min_quiet_window_length_ratio'])

    p_arrival = cp.c_int ()

    test = libsignal.signalToNoiseRatio (obs.data, syn.data, obs.data.size,
                                         min_initial_signal_to_noise_ratio, min_signal_to_noise_ratio,
                                         min_quiet_window_length_ratio, cp.byref (p_arrival))

    return test, int (p_arrival.value)


def noiseLevels (obs, syn, config):

    min_peak_ratio    = float (config['min_peak_ratio'])
    min_ratio         = float (config['min_ratio'])
    floor_ratio       = float (config['floor_ratio'])
    number_of_windows = int (config['number_of_windows'])
    max_bad_windows   = int (config['max_bad_windows'])

    ratios = np.zeros (number_of_windows, dtype = np.float64)

    compromised = libsignal.noiseLevels (obs.data, syn.data, obs.data.size, ratios, number_of_windows,
                                         min_peak_ratio, floor_ratio)

    noise_bars = list ()

    test, p_arrival = signalToNoiseRatio (obs, syn, config)

    if compromised:

        status = 'compromised'

    elif test == -1:

        noise_bars += [NoiseBar (0, p_arrival, 2, 'magenta', obs.sampling_rate)]

        return 'low signal-to-noise ratio', noise_bars

    elif test == 1:

        return 'acceptable', noise_bars

    else:

        bad_windows = 0

        total_length = obs.data.size
        box_length = total_length // number_of_windows

        l_edges = np.arange (0, total_length, box_length)
        r_edges = np.arange (box_length - 1, total_length, box_length)

        if r_edges[-1] != total_length - 1:

            r_edges[-1] = total_length - 1
            l_edges = l_edges[:-1]

        boxes = zip (l_edges, r_edges, ratios)

        for left, right, ratio in boxes:

            if ratio == -1:

                noise_bars += [NoiseBar (left, right, -1, 'green', obs.sampling_rate)]

            elif ratio < min_ratio:

                noise_bars += [NoiseBar (left, right, 1, 'black', obs.sampling_rate)]
                noise_bars += [NoiseBar (left, right, ratio, 'red', obs.sampling_rate)]

                bad_windows += 1

            elif ratio <= 1:

                noise_bars += [NoiseBar (left, right, 1, 'black', obs.sampling_rate)]
                noise_bars += [NoiseBar (left, right, ratio, 'red', obs.sampling_rate)]
                noise_bars += [NoiseBar (left, right, -1, 'green', obs.sampling_rate)]

            elif ratio <= 1.0 / min_ratio:

                noise_bars += [NoiseBar (left, right, 1, 'red', obs.sampling_rate)]
                noise_bars += [NoiseBar (left, right, 1.0 / ratio, 'black', obs.sampling_rate)]
                noise_bars += [NoiseBar (left, right, -1, 'green', obs.sampling_rate)]

            else:

                noise_bars += [NoiseBar (left, right, 1, 'red', obs.sampling_rate)]
                noise_bars += [NoiseBar (left, right, 1.0 / ratio, 'black', obs.sampling_rate)]

                bad_windows += 1

        if bad_windows > max_bad_windows:

            status = 'noisy'

        else:

            status = 'acceptable'

    return status, noise_bars


def getTracesAmplitude (obs, syn):

    return libsignal.getTracesAmplitude (obs.data, syn.data, obs.data.size)


def findEdges (data, min_allowed_index, max_allowed_index):

    boolean = np.zeros (data.size, dtype = bool)

    n = libsignal.booleanZeros (data, boolean, min_allowed_index, max_allowed_index, data.size)

    edges = np.empty (n, dtype = np.int32, order = 'C')

    libsignal.findEdges (boolean, edges, boolean.size)

    return edges


def separationIndex (trace, max_surface_waves_velocity):

    return int (trace.sampling_rate * (trace.distance / max_surface_waves_velocity - trace.begin))


def maxAllowedIndex (trace, max_surface_waves_velocity, max_period):

    return int (trace.sampling_rate * (trace.distance / max_surface_waves_velocity + 2 * max_period - trace.begin))


def getWindows (obs, syn, config):

    min_cross_correlation            = float (config['min_cross_correlation'])
    min_amplitude_ratio              = float (config['min_amplitude_ratio'])

    max_lag                          = None if config['max_lag'] == 'None' else config['max_lag']
    merge                            = True if config['merge'] == 'True' else False

    min_window_gap_in_wavelengths    = float (config['min_window_gap_in_wavelengths'])
    min_window_length_in_wavelengths = float (config['min_window_length_in_wavelengths'])
    max_surface_waves_velocity       = float (config['max_surface_waves_velocity'])
    truncate                         = float (config['truncate'])

    body_waves_only                  = True if config['body_waves_only'] == 'True' else False
    surface_waves_only               = True if config['surface_waves_only'] == 'True' else False

    min_period = int (obs.sampling_rate * obs.min_period)
    min_gap    = int (min_window_gap_in_wavelengths * min_period)
    min_length = int (min_window_length_in_wavelengths * min_period)

    min_allowed_index = 0
    max_allowed_index = maxAllowedIndex (syn, max_surface_waves_velocity, obs.max_period)
    max_allowed_index = obs.data.size - 1

    if max_allowed_index >= obs.data.size: max_allowed_index = obs.data.size - 1

    if body_waves_only:

        if surface_waves_only:

            sys.exit ("Error: cannot have 'body_waves_only' and 'surface_waves_only' flags simultaneously set to 'True'.")

        max_allowed_index = separationIndex (syn, max_surface_waves_velocity)

    elif surface_waves_only:

        min_allowed_index = separationIndex (syn, max_surface_waves_velocity)

    if not max_lag:

        max_lag = int (0.2 * min_period)

    sigma = obs.min_period / (2 * np.pi * obs.delta)

    libsignal.gaussianSmooth (obs.data, sigma, truncate, obs.data.size)
    libsignal.gaussianSmooth (syn.data, sigma, truncate, syn.data.size)

    edges = findEdges (obs.data, min_allowed_index, max_allowed_index)

    valid = np.zeros (edges.size - 1, dtype = bool)

    for index in range (1, edges.size):

        left = edges[index - 1]
        right = edges[index]

        tr1 = np.require (obs.data[left - max_lag:right + max_lag], dtype = np.float64,
                          requirements = ['C', 'A', 'W', 'O', 'E'])

        tr2 = np.require (syn.data[left:right], dtype = np.float64,
                          requirements = ['C', 'A', 'W', 'O', 'E'])

        valid[index - 1] = libsignal.similar (tr1, tr2, tr2.size, max_lag, min_cross_correlation, min_amplitude_ratio)

    if merge:

        merged = list ()

        for index in range (1, valid.size):

            if valid[index - 1] == valid[index]:

                merged += [index]

        edges = np.delete (edges, merged)
        valid = np.delete (valid, merged)

    windows = list ()

    for index in range (valid.size):

        left = edges[index]
        right = edges[index + 1]

        if valid[index] and (right - left) >= min_length:

            if windows and (left - windows[-1].right) < min_gap:

                windows[-1].right = right

            else:

                windows += [Window (left, right, obs.sampling_rate)]

    axis = np.linspace (obs.begin, obs.end, obs.npts)
    selection_region = axis[max_allowed_index] - axis[edges[0]]

    return windows, selection_region


def correlate (obs, syn, max_lag = 0):

    return libsignal.correlate (obs.data, syn.data, obs.data.size, max_lag)


def misfit (obs, syn):

    return libsignal.misfit (obs.data, syn.data, obs.data.size)


def rotateSystem (c1, c2, conversion, back_azimuth):

    if conversion == 'NE->RT':

        libsignal.rotateNE2RT (c1, c2, back_azimuth, c1.size)

    elif conversion == 'RT->NE':

        libsignal.rotateRT2NE (c1, c2, back_azimuth, c1.size)

    else:

        sys.exit ("Error: unknown type of coordinates system conversion '{}'.".format (conversion))


def plotWindows (ax, axis, c1, c2, componente, name, rotate, output, windows_file, config):

    if rotate:

        if componente == 'E':

            componente = 'T'

        elif componente == 'N':

            componente = 'R'

    amplitude = getTracesAmplitude (c1, c2)

    status, noise_bars = noiseLevels (c1, c2, config)

    win_list = list ()

    if status == 'compromised':

        result = 'Rejected'

        selected = 0

        ax.patch.set_facecolor ('gold')

    elif status == 'low signal-to-noise ratio':

        result = 'Rejected'

        selected = 0

        bar = noise_bars[0]

        left = axis[bar.left]
        right = axis[bar.right]

        rectangle = Rectangle ((left, -amplitude), bar.width (),
                                bar.ratio * amplitude,
                                color = bar.color,
                                alpha = 0.25)
        ax.add_patch (rectangle)

    elif status == 'noisy':

        result = 'Rejected'

        selected = 0

        for bar in noise_bars:

            left = axis[bar.left]
            right = axis[bar.right]

            rectangle = Rectangle ((left, 0), bar.width (),
                                    bar.ratio * amplitude,
                                    color = bar.color,
                                    alpha = 0.25)
            ax.add_patch (rectangle)

    else:

        result = 'Approved'

        if componente == 'T':

            config['max_surface_waves_velocity'] = 4.0

        windows, selection_region = getWindows (c1, c2, config)

        amplitude = getTracesAmplitude (c1, c2)

        total_width = 0

        for window in windows:

            left = axis[window.left]
            right = axis[window.right]

            rectangle = Rectangle ((left, -amplitude), window.width (),
                                    2 * amplitude,
                                    color = 'blue',
                                    alpha = 0.25)
            ax.add_patch (rectangle)

            total_width += window.width ()

        selected = 100 * total_width / selection_region

        win_list = windows

    ax.plot (axis, c1.data, color = 'black', label = 'Observed')
    ax.plot (axis, c2.data, color = 'red', label = 'Synthetic')
    ax.set_ylabel ('Amplitude (m)', fontsize = 15)
    ax.set_xlim (axis[0], axis[-1])
    ax.set_ylim (-amplitude, amplitude)

    network, station, location, channel = name.split ('.')

    channel += componente

    correlation = correlate (c1, c2)
    residual = misfit (c1, c2)

    title = '{}{} Correlation: {:.2f}  Waveform Misfit: {:.2f}  {} ({})  Selected: {:.1f}%'.format (name,
                                                                                                    componente,
                                                                                                    correlation,
                                                                                                    residual,
                                                                                                    result,
                                                                                                    status,
                                                                                                    selected)

    ax.set_title (title, fontsize = 16, y = 1)

    ax.legend (loc = 1)

    ax.tick_params (axis = 'x', labelsize = 14)
    ax.tick_params (axis = 'y', labelsize = 14)
    ax.yaxis.set_major_locator (ticker.AutoLocator ())
    ax.yaxis.set_minor_locator (ticker.AutoMinorLocator ())
    ax.yaxis.set_major_formatter (ticker.ScalarFormatter (useMathText = True))
    ax.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))

    if componente != 'Z': ax.label_outer ()
    else: ax.set_xlabel ('Time (s)', fontsize = 15)

    windows_file += '{}.{}.{}.{}\n{}\n'.format (network,
                                                station,
                                                location,
                                                channel,
                                                len (win_list))

    for window in win_list:

        windows_file += '{} {}\n'.format (window.left, window.right)

    output += '{:^8} {:^7} {:^11} {:^6} {:^+14.2f} {:^15.2f} {:^9} {:^26} {:^13.1f}\n'.format (network,
                                                                                               station,
                                                                                               location,
                                                                                               channel,
                                                                                               correlation,
                                                                                               residual,
                                                                                               result,
                                                                                               status,
                                                                                               selected)

    return output, windows_file


def readConfig ():

    config = dict ()

    with open (os.path.dirname (os.path.realpath (__file__)) + '/../Config.cfg', 'r') as FILE:

        for line in FILE:

            if line[0] == '#':

                continue

            data = line.split ()

            config[data[0]] = data[2]

    return config


def main (lock, first, last, config):

    min_period = float (config['min_period'])
    max_period = float (config['max_period'])

    rotate     = True if config['rotate'] == 'True' else False

    for event in range (first, last + 1, 4):

        observed  = getFilesNames (min_period, max_period, event, 'd')
        synthetic = getFilesNames (min_period, max_period, event, 's')

        windows = ''
        output = '{} {:>8} {:>9} {:>8} {:>12} {:>16} {:>8} {:>18} {:>23}\n'.format ('Network', 'Station', 'Location',
                                                                                    'Channel', 'Correlation',
                                                                                    'Waveform Misfit', 'Result',
                                                                                    'Status', 'Selected (%)')

        for obs, syn in zip (observed, synthetic):

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

            name = '{}.{}.{}.{}'.format (string[0], string[1], string[2], string[3][:-1])

            if not windows: windows += '{}\n'.format (z1.sampling_rate)

            fig = plt.figure (figsize = (17, 8.5))

            axis = np.linspace (e1.begin, e1.end, e1.npts)

            ax1 = fig.add_subplot (311)
            output, windows = plotWindows (ax1, axis, e1, e2, 'E', name, rotate, output, windows, config)

            ax2 = fig.add_subplot (312, sharex = ax1)
            output, windows = plotWindows (ax2, axis, n1, n2, 'N', name, rotate, output, windows, config)

            ax3 = fig.add_subplot (313, sharex = ax1)
            output, windows = plotWindows (ax3, axis, z1, z2, 'Z', name, rotate, output, windows, config)

            fig.subplots_adjust (left = 0.07, bottom = 0.07, top = 0.93, right = 0.96, wspace = 0.1, hspace = 0.18)

            lock.acquire ()

            print ('Saving {}_Windows.pdf (event {})'.format (name, event))

            lock.release ()

            plt.savefig ('events_{:.0f}_{:.0f}/event{}/{}_Windows.pdf'.format (min_period, max_period, event, name))

            plt.close (fig)

        with open ('windows_{:.0f}_{:.0f}/windows_{}.txt'.format (min_period, max_period, event), 'w') as FILE:

            FILE.write (windows[:-1])

        with open ('events_{:.0f}_{:.0f}/event{}/output_selector.txt'.format (min_period, max_period, event), 'w') as FILE:

            FILE.write (output)


def helpMenu ():

  help = """\n Error: wrong number of parameters on the comand line...

 SELECTOR

 USAGE
   ./utils/selector.py FIRST_EVENT LAST_EVENT MIN_PERIOD MAX_PERIOD

 EXAMPLE
   ./utils/selector.py 1 4 events/ 30 60 adjoint/

 COMMAND LINE ARGUMENTS
   FIRST_EVENT            - index of the first event
   LAST_EVENT             - index of the last event
   MIN_PERIOD             - minimum period of the bandpass filter
   MIN_PERIOD             - maximum period of the bandpass filter

 DESCRIPTION
   Reads observed and synthetic seismograms from the 'events/' directory and computes the time-
   widows, saved to the 'windows/' directory.\n"""

  print (help)


if __name__ == '__main__':

    if len (sys.argv) != 5:

        sys.exit (helpMenu ())

    config = readConfig ()

    lock = Lock ()

    first = int (sys.argv[1])
    last  = int (sys.argv[2])

    config['min_period'] = sys.argv[3]
    config['max_period'] = sys.argv[4]

    Process (target = main, args = (lock, first + 0, last, config)).start ()
    Process (target = main, args = (lock, first + 1, last, config)).start ()
    Process (target = main, args = (lock, first + 2, last, config)).start ()
    Process (target = main, args = (lock, first + 3, last, config)).start ()

