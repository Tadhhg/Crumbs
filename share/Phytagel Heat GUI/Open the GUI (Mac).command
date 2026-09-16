#!/bin/sh
#  Opens the Phytagel heat GUI.  Nothing to install, nothing to type.
#  Needs MATLAB R2022b or newer already on this Mac.
#
#  If double-clicking does nothing, macOS has not marked this file runnable.
#  Open Terminal, type  chmod +x   then drag this file in, and press Return.
DIR="$(cd "$(dirname "$0")" && pwd)"
CODE="$DIR/code/Heatmodel/matlab"

if [ ! -f "$CODE/heat_gui.m" ]; then
  echo "Cannot find the model code next to this file."
  echo "Keep this file inside the 'Phytagel Heat GUI' folder."
  read -r _ ; exit 1
fi

MATLAB="$(command -v matlab 2>/dev/null)"
[ -z "$MATLAB" ] && MATLAB="$(ls -d /Applications/MATLAB_R20*.app/bin/matlab 2>/dev/null | sort -r | head -1)"

if [ -z "$MATLAB" ]; then
  echo "Could not find MATLAB on this Mac."
  echo "Install MATLAB first, then double-click this file again."
  read -r _ ; exit 1
fi

echo "Starting MATLAB and opening the GUI -- this takes about 30 seconds..."
exec "$MATLAB" -sd "$CODE" -r "heat_gui"
