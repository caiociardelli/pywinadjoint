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

#include <stdio.h>
#include <math.h>

static const unsigned CNMAX   = 100;
static const double THRESHOLD = 1.001;


static inline unsigned index (unsigned i, unsigned j, unsigned n)
{
    return n * i + j;
}

static inline double square (double x)
{
    return x * x;
}

static double max (double D[], double Cn[], unsigned i)
{
    double x0 = D[i - 1];
    double x1 = D[i];
    double x2 = D[i + 1];
    double y0 = Cn[i - 1];
    double y1 = Cn[i];
    double y2 = Cn[i + 1];

    double f0 = y0;
    double f01 = (y1 - y0) / (x1 - x0);
    double f12 = (y2 - y1) / (x2 - x1);
    double f02 = (f12 - f01) / (x2 - x0);

    double a = f02;
    double b = f01 - f02 * (x1 + x0);
    double c = f0 + (f02 * x1 - f01) * x0;

    return (4 * a * c - square (b)) / (4 * a);
}

static double interpolate (double D[], double Cn[], unsigned i, double cn)
{
    double x0 = D[i - 1];
    double x1 = D[i];
    double x2 = D[i + 1];
    double y0 = Cn[i - 1];
    double y1 = Cn[i];
    double y2 = Cn[i + 1];

    double m = (y2 - y0) / (x2 - x0);
    double n = y0 - m * x0 - cn;

    double r = -n / m;

    double f0 = y0;
    double f01 = (y1 - y0) / (x1 - x0);
    double f12 = (y2 - y1) / (x2 - x1);
    double f02 = (f12 - f01) / (x2 - x0);

    double a = f02;
    double b = f01 - f02 * (x1 + x0);
    double c = f0 + (f02 * x1 - f01) * x0 - cn;

    double delta = square (b) - 4 * a * c;

    if (delta < 0) return r;

    double r1 = (-b + sqrt (delta)) / (2 * a);
    double r2 = (-b - sqrt (delta)) / (2 * a);

    return (fabs (r1 - r) < fabs (r2 - r)) ? r1 : r2;
}

void weights (double *pairs, double *w, double step, double cN, unsigned nC, unsigned n)
{
    unsigned imax = 0;
    double D[CNMAX];
    double Cn[CNMAX];

    D[0]  = 0;
    Cn[0] = 1;

    for (unsigned i = 1; i < CNMAX; i++)
    {
        double c;
        D[i] = step * i;

        double min_w =  1e100;
        double max_w = -1e100;

        for (unsigned j = 0; j < n; j++)
        {
            c = 0;

            for (unsigned k = 0; k < n; k++)

                c += exp (-square (pairs[index (j, k, n)] / D[i]));

            w[j] = 1 / c;

            if (w[j] < min_w) min_w = w[j];
            if (w[j] > max_w) max_w = w[j];
        }

        Cn[i] = max_w / min_w;

        if (i > 1 && Cn[i - 1] > Cn[i] && Cn[i - 1] > Cn[i - 2])
        {
            imax = i - 1;
            break;
        }
    }

    double max_cn = max (D, Cn, imax);
    double cn = (cN * max_cn > THRESHOLD) ? cN * max_cn : THRESHOLD;

    unsigned i = 1;

    while (Cn[i + 1] < cn) i++;

    if (i > 1 && Cn[i + 1] - cn > cn - Cn[i]) i--;

    double c, sw = 0;
    double delta = interpolate (D, Cn, i, cn);

    for (unsigned i = 0; i < n; i++)
    {
        c = 0;

        for (unsigned j = 0; j < n; j++)

            c += exp (-square (pairs[index (i, j, n)] / delta));

        w[i] = 1 / c;
        sw += w[i];
    }

    double Nf = 1 / (sw * (double) nC);

    for (unsigned i = 0; i < n; i++)

        w[i] *= Nf;
}

