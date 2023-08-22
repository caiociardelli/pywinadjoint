#!/usr/bin/env bash

# PyWinAdjoint

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

# EXPAND

# USAGE
#   ./utils/expand.bash FIRST_EVENT LAST_EVENT

# DESCRIPTION
#   Convert the adjoint sources from binary to ASCII from the FIRST_EVENT to the LAST_EVENT.

#-----------------------------------------------------------------------------------------------

echo 'Erasing previous directory...'

rm -rfv adjoint_ascii/

echo 'Creating directory to store the adjoint sources in ASCII format'

mkdir adjoint_ascii

for index in $(seq $1 $2)
do
  mkdir adjoint_ascii/event$index

  for file in $(ls adjoint/event$index/*.bin)
  do
    name="adjoint_ascii/event$index/$(echo $file | cut -d'/' -f3 | cut -d'_' -f1)".adj
    ./bin/expand $file > $name
    echo $name
  done

  cp -v adjoint/event$index/STATIONS_ADJOINT adjoint_ascii/event$index/
done

echo 'Done!'
