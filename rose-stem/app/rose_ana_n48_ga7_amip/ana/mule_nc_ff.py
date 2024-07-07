#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
import os
import re
import mule
import numpy as np
from netCDF4 import Dataset
import mule.stashmaster
from rose.apps.rose_ana import AnalysisTask

# For rose_ana we want to use the suite's STASHmaster instead of the centrally
# installed one (the user may have made required changes to it)
um_dir = "UM_INSTALL_DIR"
if um_dir not in os.environ:
    raise ValueError(um_dir + " not set")
mule.stashmaster.STASHMASTER_PATH_PATTERN = (
    re.sub(r"UMDIR", um_dir, mule.stashmaster.STASHMASTER_PATH_PATTERN))


class CompareFFandNC(AnalysisTask):
    """Compare a FieldsFile and NetCDF file."""
    def run_analysis(self):
        """Main analysis routine called from rose_ana."""
        self.process_opt_files()
        # Deal with any unexpected options
        self.process_opt_unhandled()

        # Currently this analysis class can only handle comparing two files
        if len(self.files) != 2:
            raise ValueError("Must specify exactly two files for comparison.")

        self.perform_comparison()

    def perform_comparison(self):
        """Compare the files using mule-cumf."""
        # Turn the filenames into their absolute equivalents
        file1 = os.path.realpath(self.files[0])
        file2 = os.path.realpath(self.files[1])

        # Load the first file using Mule
        self.ff = mule.FieldsFile.from_file(file1, remove_empty_lookups=True)

        # And the second one using NetCDF
        self.nc = Dataset(file2, endian="little")

        # The relevant domain and time names from the app
        mlev_domname = "D10TH"
        mlev_timname = "T12H"
        slev_timname = "TONCE"

        # Lists of the level and time numbers for multi-level/time fields
        nc_levels = self.nc.variables[mlev_domname + "_model_level_number"]
        nc_times = self.nc.variables[mlev_timname]

        # Counters
        passed, failed = 0, 0

        for ifield, field in enumerate(self.ff.fields):
            if field.lbrel in (2, 3):
                # Convert the Field's stash code into a NetCDF variable name
                # and find the associated variable
                nc_varname = (
                    "STASH_m01s{0:02d}i{1:03d}"
                    .format(field.lbuser4 / 1000, field.lbuser4 % 1000))
                nc_var = self.nc.variables[nc_varname]

                # Extract and save the field data
                ff_data = field.get_data()

                # Now find the level and forecast hour of the field
                fld_level_value = field.lblev
                fld_time_value = field.lbft

                # The name of the level dimension will depend on what
                # type of multi-level field this is; use the STASH
                # entry to find out
                if field.stash.levelT == 2:
                    level_dimname = mlev_domname + "_eta_theta"
                elif field.stash.levelT == 1:
                    level_dimname = mlev_domname + "_eta_rho"

                # Detect if this is one of the multi-level parameters
                if (nc_var.ndim == 4
                        and level_dimname in nc_var.dimensions
                        and mlev_timname in nc_var.dimensions):

                    # And the equivalent for the netCDF file, which will be
                    # a given index in the dimension
                    nc_level_pos = (
                        list(nc_levels[:]).index(fld_level_value))
                    nc_time_pos = (
                        list(nc_times[:]).index(fld_time_value))

                    # The indices of the time and level dimensions in the
                    # netCDF variable ( which are needed to index it)
                    level_index = nc_var.dimensions.index(level_dimname)
                    time_index = nc_var.dimensions.index(mlev_timname)

                    # Setup the indexing array, and populate the correct
                    # elements with the correct indices from above
                    indices = [
                        slice(None), slice(None), slice(None), slice(None)]
                    indices[level_index] = nc_level_pos
                    indices[time_index] = nc_time_pos

                    # Now indexing the variable should give the correct slice
                    nc_data = nc_var[indices]
                elif (nc_var.ndim == 3
                      and slev_timname in nc_var.dimensions):
                    # If this is one of the single level fields, the slice
                    # just needs to omit the single time profile
                    time_index = nc_var.dimensions.index(slev_timname)
                    indices = [
                        slice(None), slice(None), slice(None)]
                    indices[time_index] = 0
                    # Now indexing the variable should give the correct slice
                    nc_data = nc_var[indices]
                else:
                    prefix = "[FAIL] "
                    self.parent.reporter(
                        "Unable to match Field {0} ({1}) to an equivalent "
                        "in the NetCDF file".format(ifield, field.stash.name),
                        prefix=prefix)
                    failed += 1

            # Should now have the 2 data arrays to compare
            if np.sum(nc_data - ff_data.astype(nc_data.dtype)) == 0:
                passed += 1
            else:
                prefix = "[FAIL] "
                self.parent.reporter(
                    "Field {0} ({1}) does not compare between the files"
                    .format(ifield, field.stash.name), prefix=prefix)
                failed += 1

        # Check if everything passed
        if failed == 0:
            self.passed = True
            prefix = "[INFO] "
            self.parent.reporter(
                "{0} compared fields matched".format(passed), prefix=prefix)
        else:
            prefix = "[FAIL] "
            self.parent.reporter(
                "{0} compared fields matched, but {1} did not"
                .format(passed, failed))

    def process_opt_files(self):
        """Process the files option; a list of one or more filenames."""
        # Get the file list from the options dictionary
        files = self.options.pop("files", None)
        # Make sure it appears as a sensible list
        if files is None:
            files = []
        elif isinstance(files, str):
            files = [files]
        self.files = files

        # Report the filenames
        for ifile, fname in enumerate(files):
            self.parent.reporter(
                "File {0}: {1}".format(ifile + 1,
                                       os.path.realpath(files[ifile])))
