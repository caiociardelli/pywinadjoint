#!/usr/bin/env bash

# PyWinEPAdjoint

# Author: Caio Ciardelli, University of São Paulo, May 2021

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#-----------------------------------------------------------------------------------------------

# RUN_EXAMPLE

# USAGE
#   ./run_example.bash

# DESCRIPTION
#   Runs the example.

#-----------------------------------------------------------------------------------------------

FIRST_EVENT=1
LAST_EVENT=4

MIN_PERIOD=30
MAX_PERIOD=60

INPUT_DIR='events'
OUTPUT_DIR='adjoint'

echo 'Cleaning previous run...'

rm -rv events/event?/*.txt
rm -rv events/event?/*.pdf
rm -rfv windows/
rm -rfv adjoint/

echo 'Creating directory to store windows'

mkdir windows

echo 'Running window selector...'

./utils/selector.py $FIRST_EVENT $LAST_EVENT $MIN_PERIOD $MAX_PERIOD

echo 'Creating directory to store the adjoint sources...'

for index in $(seq $FIRST_EVENT $LAST_EVENT);
do
  mkdir -p adjoint/event$index;
done

echo 'Creating the adjoint sources...'

./utils/adjoint.py $FIRST_EVENT $LAST_EVENT $INPUT_DIR $MIN_PERIOD $MAX_PERIOD $OUTPUT_DIR

echo 'Done!'

