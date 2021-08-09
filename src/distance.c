/*
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

----------------------------------------------------------------------------------------------- */

#include <math.h>
#include "event.h"

static const double R = 6378.137;
static const double TO_RAD = 3.14159265358979323846 / 180;
static const double TO_DEG = 180 / 3.14159265358979323846;

static inline double square (double v)
{
    return v * v;
}

static inline double degrees2Radians (double angle)
{
    return angle * TO_RAD;
}

static inline double radians2Degrees (double angle)
{
    return angle * TO_DEG;
}

double haversine (double lat1, double lon1, double lat2, double lon2)
{
    double dlat = degrees2Radians (lat2 - lat1);
    double dlon = degrees2Radians (lon2 - lon1);

    return 2 * radians2Degrees (asin (sqrt (square (sin (0.5 * dlat)) + cos (degrees2Radians (lat2))
                                                                      * cos (degrees2Radians (lat1))
                                                                      * square (sin (0.5 * dlon)))));
}

double vincenty (double lat1, double lon1, double lat2, double lon2)
{
    lat1 = degrees2Radians (lat1);
    lat2 = degrees2Radians (lat2);
    lon1 = degrees2Radians (lon1);
    lon2 = degrees2Radians (lon2);

    double sin_lat1 = sin (lat1);
    double cos_lat1 = cos (lat1);
    double sin_lat2 = sin (lat2);
    double cos_lat2 = cos (lat2);

    double sin_dlon =  sin (lon2 - lon1);
    double cos_dlon =  cos (lon2 - lon1);

    return radians2Degrees (atan2 (sqrt (square (cos_lat2 * sin_dlon)
                                   + square (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_dlon)),
                                   sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_dlon));
}

double distance (struct Event *event1, struct Event *event2)
{
    double r1 = R - event1->depth;
    double r2 = R - event2->depth;

    return sqrt (square (r1) + square (r2) - 2 * r1 * r2 * cos (degrees2Radians (vincenty (event1->latitude,
                                                                                           event1->longitude,
                                                                                           event2->latitude,
                                                                                           event2->longitude))));
}

