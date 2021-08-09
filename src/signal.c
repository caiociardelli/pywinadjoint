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
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>

static const double TO_RAD = 3.14159265358979323846 / 180;


static char machineEndianess (void)
{
    int e = 1;

    if (*(char*) &e == 1)

        return '<';

    return '>';
}

static int intSwap (int input)
{
    int output;

    char *in = (char*) &input;
    char *out = (char*) &output;

    out[0] = in[3];
    out[1] = in[2];
    out[2] = in[1];
    out[3] = in[0];

    return output;
}


static float floatSwap (float input)
{
    float output;

    char *in = (char*) &input;
    char *out = (char*) &output;

    out[0] = in[3];
    out[1] = in[2];
    out[2] = in[1];
    out[3] = in[0];

    return output;
}

static double doubleSwap (double input)
{
    double output;

    char *in = (char*) &input;
    char *out = (char*) &output;

    out[0] = in[7];
    out[1] = in[6];
    out[2] = in[5];
    out[3] = in[4];
    out[4] = in[3];
    out[5] = in[2];
    out[6] = in[1];
    out[7] = in[0];

    return output;
}

static inline double square (double value)
{
    return value * value;
}

static inline double degrees2Radians (double angle)
{
    return angle * TO_RAD;
}

static inline bool isLocalMax (double *v, int i)
{
    return v[i - 1] < v[i] && v[i] > v[i + 1];
}

static inline bool isLocalMin (double *v, int i)
{
    return v[i - 1] > v[i] && v[i] < v[i + 1];
}

static inline bool noSignalReverse (double *v, int i)
{
    return (v[i] > 0 && v[i + 1] > 0) || (v[i] < 0 && v[i + 1] < 0) ? true : false;
}

static double max (double *x, int n)
{
    double max = x[0];

    for (int i = 1; i < n; i++)

        if (x[i] > max) max = x[i];

    return max;
}

static double min (double *x, int n)
{
    double min = x[0];

    for (int i = 1; i < n; i++)

        if (x[i] < min) min = x[i];

    return min;
}

static double threshold (double *v, int n)
{
    double max_v = max (v, n);
    double min_v = min (v, n);

    return max_v > -min_v ? 0.0075 * max_v : -0.0075 * min_v;
}

static bool isLargeMax (double *v, int i, double tol)
{
    int l = i - 2;

    while (l > 0 && (v[l - 1] < v[l] || v[l] > v[l + 1]))

        l--;

    if (v[i] - v[l] > tol)

        return true;

    int r = i + 2;

    while (v[r - 1] < v[r] || v[r] > v[r + 1])

        r++;

    if (v[i] - v[r] > tol)

        return true;

    return false;
}

static bool isLargeMin (double *v, int i, double tol)
{
    int l = i - 2;

    while (l > 0 && (v[l - 1] > v[l] || v[l] < v[l + 1]))

        l--;

    if (v[l] - v[i] > tol)

        return true;

    int r = i + 2;

    while (v[r - 1] > v[r] || v[r] < v[r + 1])

        r++;

    if (v[r] - v[i] > tol)

        return true;

    return false;
}

static double amplitudeRatio (double *d, double *s, int n, int lag)
{
    double n1 = 0, n2 = 0;

    for (int i = 0; i < n; i++)
    {
        n1 += fabs (d[i + lag]);
        n2 += fabs (s[i]);
    }

    if (n1 > 0 && n2 > 0)

        if (n1 > n2)

            return n2 / n1;

        else if (n1 < n2)

            return n1 / n2;

        else

            return 1;

    else

        return 0;
}

static double simpson13 (double a, double b, double *f, int n)
{
    double sum1 = f[1];
    double sum2 = 0.0;

    for (int i = 2; i < n - 2; i += 2)
    {
        sum1 += f[i + 1];
        sum2 += f[i];
    }

    return (b - a) / (3.0 * (n - 1)) * (f[0] + f[n - 1] + 4.0 * sum1 + 2.0 * sum2);
}

static inline double simpson38 (double a, double b, double *f, int n)
{
    return (b - a) * (f[n - 4] + 3 * (f[n - 3] + f[n - 2]) + f[n - 1]) / 8.0;
}

static double integrate (double a, double b, double *f, int n)
{
    if (n == 1)

        return (b - a) * f[0];

    else if (n == 2)

        return (b - a) * (f[1] - f[0]) / 2.0;

    else if (n % 2 == 1)

        return simpson13 (a, b, f, n);

    else if (n == 4)

        return simpson38 (a, b, f, n);

    else
    {
        double m = a + (n - 4) * (b - a) / (n - 1);

        return simpson13 (a, m, f, n - 3) + simpson38 (m, b, &f[n - 4], 4);
    }
}

static void slidingCorrelation (int m, double weights[], double *x, int n)
{
    int hw = m / 2, k;

    double c[n];

    memcpy (c, x, n * sizeof (double));

    for (int i = 0; i < n; i++)
    {
        x[i] = 0;

        for (int j = 0; j < m; j++)
        {
            k = i + j - hw;

            x[i] += (k >= 0 && k < n) ? c[k] * weights[j] : 0;
        }
    }
}

void gaussianSmooth (double *x, double sigma, double truncate, int n)
{
    int lw = (int) (truncate * sigma + 0.5), m = 2 * lw + 1;

    double sum = 1.0, sigma2 = square (sigma), weight, weights[m];

    weights[lw] = 1.0;

    for (int i = 1; i < lw + 1; i++)
    {
        weight = exp (-0.5 * square (i) / sigma2);

        weights[lw + i] = weight;
        weights[lw - i] = weight;

        sum += 2.0 * weight;
    }

    for (int i = 0; i < m; i++)

        weights[i] /= sum;

    slidingCorrelation (m, weights, x, n);
}

double getTracesAmplitude (double *d, double *s, int n)
{
    double absmax = -1e15;

    for (int i = 0; i < n; i++)
    {
        if (fabs (d[i]) > absmax)

            absmax = fabs (d[i]);

        if (fabs (s[i]) > absmax)

            absmax = fabs (s[i]);
    }

    return 1.05 * absmax;
}

double getAdjSourceAmplitude (double *adjE, double *adjN, double *adjZ, int n)
{
    double absmax = -1e15;

    for (int i = 0; i < n; i++)
    {
        if (fabs (adjE[i]) > absmax)

            absmax = fabs (adjE[i]);

        if (fabs (adjN[i]) > absmax)

            absmax = fabs (adjN[i]);

        if (fabs (adjZ[i]) > absmax)

            absmax = fabs (adjZ[i]);
    }

    return absmax > 0 ? 1.05 * absmax : 1.05;
}


int booleanZeros (double *v, bool *boolean, int min_i, int max_i, int n)
{
    int i = 0, c = 1;
    double tol = threshold (v, n);

    if (max_i >= n)

        max_i = n - 1;

    while (i < min_i || fabs (v[i]) < tol)

        boolean[i++] = false;

    boolean[i++] = true;

    while (i < max_i)

        if (isLocalMax (v, i) && isLargeMax (v, i, 10 * tol))
        {
            while (i < max_i && noSignalReverse (v, i))

                boolean[i++] = false;

            boolean[i++] = true;

            c++;
        }

        else if (isLocalMin (v, i) && isLargeMin (v, i, 10 * tol))
        {
            while (i < max_i && noSignalReverse (v, i))

                boolean[i++] = false;

            boolean[i++] = true;

            c++;
        }

        else

            boolean[i++] = false;

    while (i < n)

        boolean[i++] = false;

    return c;
}

void findEdges (bool *boolean, int *edges, int n)
{
    int j = 0;

    for (int i = 0; i < n; i++)

        if (boolean[i])

            edges[j++] = i;
}

double misfit (double *d, double *s, int n)
{
    double sum = 0, n1 = 0, n2 = 0;

    for (int i = 0; i < n; i++)
    {
        n1 += d[i] * d[i];
        n2 += s[i] * s[i];
    }

    if (n1 > 0 && n2 > 0)
    {
        n1 = sqrt (n1);
        n2 = sqrt (n2);

        for (int i = 0; i < n; i++)

            sum += square (n2 * d[i] - n1 * s[i]);

        return sqrt (sum) / (2 * n1 * n2);
    }

    else

        return 0;
}

double correlate (double *d, double *s, int n, int lag)
{
    double sum = 0, n1 = 0, n2 = 0;

    for (int i = 0; i < n; i++)
    {
        n1 += d[i + lag] * d[i + lag];
        n2 += s[i] * s[i];

        sum += d[i + lag] * s[i];
    }

    if (n1 > 0 && n2 > 0)

        return sum / sqrt (n1 * n2);

    else

        return 0;
}

bool similar (double *d, double *s, int n, int lag, double min_cc, double min_rt)
{
    double max = -1;

    for (int l = 0; l <= 2 * lag; l++)
    {
        double cr = correlate (d, s, n, l);

        if (cr > max)

            max = cr;

        if (max >= min_cc && amplitudeRatio (d, s, n, l) > min_rt)

            return true;
    }

    return false;
}

int signalToNoiseRatio (double *d, double *s, int n, double min_instnr, double min_stnr, double min_qwlr, int *p_arrival)
{
    int min_qwl = (int) (min_qwlr * n);

    double max1 = 0, max2 = 0;

    for (int i = 0; i < n; i++)
    {
        if (fabs (d[i]) > max1)

            max1 = fabs (d[i]);

        if (fabs (s[i]) > max2)

            max2 = fabs (s[i]);
    }

    *p_arrival = 0;

    while (fabs (s[*p_arrival]) < 0.01 * max2)

        *p_arrival += 1;

    *p_arrival -= (int) (0.02 * n);

    if (*p_arrival > min_qwl)
    {
        double sum = 0;

        for (int i = 0; i < *p_arrival; i++)
        {
            double value = fabs (d[i]);

            if (value * 0.35 * min_instnr > max1)

                return -1;

            sum += value;
        }

        double mean = sum / (*p_arrival);

        if (mean * min_instnr > max1)

            return -1;
    }

    int qwl = 0;

    for (int i = 0; i < n; i++)
    {
        if (fabs (d[i]) * min_stnr <= max1)

            qwl++;

        else

            qwl = 0;

        if (qwl == min_qwl)

            return 1;
    }

    return 0;
}

bool noiseLevels (double *d, double *s, int size, double *ratios, int n, double min_pkr, double flr_rt)
{
    double max1 = 0, max2 = 0;

    for (int i = 0; i < size; i++)
    {
        if (fabs (d[i]) > max1)

            max1 = fabs (d[i]);

        if (fabs (s[i]) > max2)

            max2 = fabs (s[i]);
    }

    if (min_pkr * max2 > max1 || min_pkr * max1 > max2)

        return true;

    double floor;

    if (max2 > max1)

        floor = flr_rt * max2;

    else

        floor = flr_rt * max1;

    int window = size / n;

    for (int i = 0, l = 0, r = window; i < n; i++, l = r, r += window)
    {
        double sum1 = 0, sum2 = 0;

        while (l < r)
        {
            sum1 += fabs (d[l++]);
            sum2 += fabs (s[l++]);
        }

        double mean1 = sum1 / window;
        double mean2 = sum2 / window;

        if (mean1 > floor)

            ratios[i] = mean2 / mean1;

        else

            ratios[i] = -1;
    }

    return false;
}

void analyticSignal (double *x, double complex *z, int n)
{
    double complex X[n], Z[n];
    fftw_plan plan;

    plan = fftw_plan_dft_r2c_1d (n, x, X, FFTW_ESTIMATE);
    fftw_execute (plan);
    fftw_destroy_plan (plan);

    Z[0] = X[0];

    for (int i = 1; i < n / 2; i++)

        Z[i] = 2 * X[i];

    if (n % 2 == 0)

        Z[n / 2] = X[n / 2];

    else

        Z[n / 2] = 2 * X[n / 2];

    for (int i = n / 2 + 1; i < n; i++)

        Z[i] = 0;

    plan = fftw_plan_dft_1d (n, Z, z, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute (plan);
    fftw_destroy_plan (plan);

    for (int i = 0; i < n; i++)

        z[i] /= n;
}

void hilbert (double *x, double *h, int n)
{
    double complex z[n];

    analyticSignal (x, z, n);

    for (int i = 0; i < n; i++)

        h[i] = cimag (z[i]);
}

void taper (double *x, double width, double c1, double c2, double omega, int n)
{
    int length = width * n;

    for (int i = 0; i < length; i++)
    {
        x[i] *= (c1 - c2 * cos (omega * i));
        x[n - i - 1] *= (c1 - c2 * cos (omega * i));
    }
}

void rotateNE2RT (double *N, double *E, double ba, int n)
{
    ba = degrees2Radians (ba);

    for (int i = 0; i < n; i++)
    {
        double n = N[i];
        double e = E[i];

        N[i] = -e * sin (ba) - n * cos (ba);
        E[i] = -e * cos (ba) + n * sin (ba);
    }
}

void rotateRT2NE (double *R, double *T, double ba, int n)
{
    ba = degrees2Radians (ba);

    for (int i = 0; i < n; i++)
    {
        double r = R[i];
        double t = T[i];

        R[i] = t * sin (ba) - r * cos (ba);
        T[i] = -t * cos (ba) - r * sin (ba);
    }
}

void adjointSourceWF (double *d, double *s, double *w, double *adjs, double a, double b,
                      double *Xi, double *Xi_w, int n)
{
    double dwf[n], dwfw[n];

    for (int i = 0; i < n; i++)
    {
        dwf[i]  = s[i] - d[i];
        dwfw[i] = w[i] * dwf[i];

        adjs[i] = dwfw[i];

        dwf[i]  = square (dwf[i]);
        dwfw[i] = square (dwfw[i]);
    }

    *Xi   = 0.5 * integrate (a, b, dwf, n);
    *Xi_w = 0.5 * integrate (a, b, dwfw, n);
}

void adjointSourceEP (double *d, double *s, double *w, double *adjs, double a, double b,
                      double water_level, double *Xi, double *Xi_w, int n)
{
    double complex ans_d[n], ans_s[n];

    analyticSignal (d, ans_d, n);
    analyticSignal (s, ans_s, n);

    double E_s[n], k1[n], k2[n], dphi[n], dphiw[n];

    for (int i = 0; i < n; i++)

        E_s[i] = cabs (ans_s[i]);

    double wtl = water_level * max (E_s, n);

    for (int i = 0; i < n; i++)
    {
        double H_d = cimag (ans_d[i]);
        double H_s = cimag (ans_s[i]);

        double E_2 = square (E_s[i] + wtl);

        dphi[i]  = atan2 (s[i] * H_d - d[i] * H_s, d[i] * s[i] + H_d * H_s);
        dphiw[i] = w[i] * dphi[i];

        k1[i] = sin (dphiw[i]) * H_s / E_2;
        k2[i] = sin (dphiw[i]) * s[i] / E_2;

        dphi[i]  = 2 - 2 * cos (dphi[i]);
        dphiw[i] = 2 - 2 * cos (dphiw[i]);
    }

    double K2[n];

    hilbert (k2, K2, n);

    for (int i = 0; i < n; i++)

        adjs[i] = k1[i] + K2[i];

    *Xi   = 0.5 * integrate (a, b, dphi, n);
    *Xi_w = 0.5 * integrate (a, b, dphiw, n);
}

void adjointSourceEV (double *d, double *s, double *w, double *adjs, double a, double b,
                      double water_level, double *Xi, double *Xi_w, int n)
{
    double complex ans_d[n], ans_s[n];

    analyticSignal (d, ans_d, n);
    analyticSignal (s, ans_s, n);

    double E_s[n], E_d[n], k1[n], k2[n], denv[n], denvw[n];

    for (int i = 0; i < n; i++)
    {
        E_s[i] = cabs (ans_s[i]);
        E_d[i] = cabs (ans_d[i]);
    }

    double wtl = water_level * max (E_s, n);

    for (int i = 0; i < n; i++)
    {
        double H_s = cimag (ans_s[i]);

        double E_2 = square (E_s[i] + wtl);

        denv[i]  = log (E_d[i] / E_s[i]);
        denvw[i] = w[i] * denv[i];

        k1[i] = denvw[i] * H_s / E_2;
        k2[i] = denvw[i] * s[i] / E_2;

        denv[i]  = square (denv[i]);
        denvw[i] = square (denvw[i]);
    }

    double K2[n];

    hilbert (k2, K2, n);

    for (int i = 0; i < n; i++)

        adjs[i] = k1[i] + K2[i];

    *Xi   = 0.5 * integrate (a, b, denv, n);
    *Xi_w = 0.5 * integrate (a, b, denvw, n);
}

int readBinaryLength (char *name, int *n)
{
    FILE *input = fopen (name, "rb");

    if (input == NULL)
    {
        return 1;
    }

    int rbytes; char endianess;

    rbytes = fread (&endianess, sizeof (char), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 2;
    }

    bool swap = false;

    if (endianess != machineEndianess ())

        swap = true;

    double begin;

    rbytes = fread (&begin, sizeof (double), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 3;
    }

    if (swap) begin = doubleSwap (begin);

    double end;

    rbytes = fread (&end, sizeof (double), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 4;
    }

    if (swap) end = doubleSwap (end);

    rbytes = fread (n, sizeof (int), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 5;
    }

    if (swap) *n = intSwap (*n);

    char type;

    rbytes = fread (&type, sizeof (char), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 6;
    }

    fclose (input);

    return 0;
}

int readBinary (char *name, float *axis, float *array)
{
    FILE *input = fopen (name, "rb");

    if (input == NULL)
    {
        return 1;
    }

    int rbytes; char endianess;

    rbytes = fread (&endianess, sizeof (char), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 2;
    }

    bool swap = false;

    if (endianess != machineEndianess ())

        swap = true;

    double begin;

    rbytes = fread (&begin, sizeof (double), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 3;
    }

    if (swap) begin = doubleSwap (begin);

    double end;

    rbytes = fread (&end, sizeof (double), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 4;
    }

    if (swap) end = doubleSwap (end);

    int n;

    rbytes = fread (&n, sizeof (int), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 5;
    }

    if (swap) n = intSwap (n);

    char type;

    rbytes = fread (&type, sizeof (char), 1, input);

    if (rbytes != 1)
    {
        fclose (input); return 6;
    }

    if (type == 'z')

        for (int i = 0; i < n; i++)

            array[i] = 0;

    else
    {
        rbytes = fread (array, sizeof (float), n, input);

        if (rbytes != n)
        {
            fclose (input); return 7;
        }

        if (swap)

            for (int i = 0; i < n; i++)

                array[i] = floatSwap (array[i]);
    }

    double dt = (end - begin) / (n - 1);

    for (int i = 0; i < n; i++)

        axis[i] = begin + i * dt;

    fclose (input);

    return 0;
}

void writeBinary (char *name, char endianess, float *x, double begin, double end, int n)
{
    FILE *output = fopen (name, "wb");

    fwrite (&endianess, sizeof (char), 1, output);
    fwrite (&begin, sizeof (double), 1, output);
    fwrite (&end, sizeof (double), 1, output);
    fwrite (&n, sizeof (int), 1, output);

    char type = 'z';

    for (int i = 0; i < n; i++)

        if (fabs (x[i]) > 1e-15)
        {
            type = 'v';

            break;
        }

    fwrite (&type, sizeof (char), 1, output);

    if (type == 'v')

        fwrite (x, sizeof (float), n, output);

    fclose (output);
}

