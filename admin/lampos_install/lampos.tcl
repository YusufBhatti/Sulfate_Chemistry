#!/usr/bin/env bash
# the next line restarts wish \
exec wish "$0" "$@"

# This is the top level script for LAMPOS
#
# Operation is fully described in the help documents.
#

set script [ file tail $argv0 ]
puts stdout $script

# Load the Tcl scripts.

set tcl_dir "[ file dirname $argv0 ]/source/Tcl"
cd $tcl_dir
set source_files [glob *.tcl]
foreach source_file $source_files {
  source $source_file
}

# Initialise variables for operations and the GUI. Load libraries.
#initialise

# Enters the event loop.
main


