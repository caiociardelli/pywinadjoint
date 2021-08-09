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

int main (int argc, char **argv)
{
    if (argc != 2)
    {
        fprintf (stderr, "Error: Wrong number of parameters at the command line...\n");

        return 1;
    }

    FILE *input = fopen (argv[1], "rb");

    if (input == NULL)
    {
        fprintf (stderr, "Error: Could not open %s file...", argv[1]);

        return 1;
    }

    char endianess, type;
    bool swap = false;
    double begin, end;
    int rbytes, n;

    rbytes = fread (&endianess, sizeof (char), 1, input);

    if (rbytes != 1)
    {
        fprintf (stderr, "Error: Could not read in the endianess of the file...\n");
        fclose (input);

        return 1;
    }

    if (endianess != machineEndianess ())

        swap = true;

    rbytes = fread (&begin, sizeof (double), 1, input);

    if (rbytes != 1)
    {
        fprintf (stderr, "Error: Could not read in the starttime...\n");
        fclose (input);

        return 1;
    }

    if (swap)

        begin = doubleSwap (begin);

    rbytes = fread (&end, sizeof (double), 1, input);

    if (rbytes != 1)
    {
        fprintf (stderr, "Error: Could not read in the endtime...\n");
        fclose (input);

        return 1;
    }

    if (swap)

        end = doubleSwap (end);

    rbytes = fread (&n, sizeof (int), 1, input);

    if (rbytes != 1)
    {
        fprintf (stderr, "Error: Could not read in the length of the file...\n");
        fclose (input);

        return 1;
    }

    if (swap)

        n = intSwap (n);

    float array[n];

    rbytes = fread (&type, sizeof (char), 1, input);

    if (rbytes != 1)
    {
        fprintf (stderr, "Error: Could not read in the type of data values...\n");
        fclose (input);

        return 1;
    }

    if (type == 'z')

        for (int i = 0; i < n; i++)

            array[i] = 0;

    else
    {
        rbytes = fread (array, sizeof (float), n, input);

        if (rbytes != n)
        {
            fprintf (stderr, "Error: Could not read in the right number of non-zero values...\n");
            fclose (input);

            return 1;
        }

        if (swap)

            for (int i = 0; i < n; i++)

                array[i] = floatSwap (array[i]);
    }

    double dt = (end - begin) / (n - 1);

    for (int i = 0; i < n; i++)

        fprintf (stdout, "%1.6e %1.6e\n", begin + i * dt, array[i]);

    fclose (input);

    return 0;
}

