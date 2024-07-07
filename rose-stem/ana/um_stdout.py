#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

import os
import re
import sys
import itertools
from rose.apps.rose_ana import AnalysisTask

UNKNOWN_SECTION = -1
TIMESTEP_SECTION = 1
NORM_SECTION = 2
SEARCH_PATTERN = '\*\s+\d+\s+\d+\s+\d+\s+(\S+)\s+\*'
TIMESTEP_HEADER_STRING = (
    'Model time:[ ]+([0-9]{2,4}[-]{0,1}){3}[ ]+([0-9]{2,4}[:]{0,1}){3}')
NORM_HEADER_STRING = 'Linear solve for Helmholtz problem'
NORM_VALUE_STRING = '\\*[ ]*[0-9]+[ ]+[0-9]+[ ]+[0-9]+[ ]+.*\\*'
TIMESTEP_NUMBER_PATTERN = 'Atm_Step: Timestep[ ]+[0-9]+'


class Timestep(object):
    """
    Represent a timestep in a model pe_output file. Stores the timestamp and
    output norms for the timestep
    """
    def __init__(self, year, month, day, hour, minute, second, number):
        self.Year = year
        self.Month = month
        self.Day = day
        self.Hour = hour
        self.Minute = minute
        self.Second = second
        self.Number = number
        self.NormList = []

    def __str__(self):
        retVal = (
            '{year:04}/{month:02}/{day:02} {hour:02}:{minute:02}:{second:02}')
        retVal = retVal.format(year=self.Year,
                               month=self.Month,
                               day=self.Day,
                               hour=self.Hour,
                               minute=self.Minute,
                               second=self.Second)
        return retVal

    def addNorm(self, norm1):
        self.NormList += [norm1]

    def compareTimes(self, other):
        if self.Year != other.Year:
            return False
        if self.Month != other.Month:
            return False
        if self.Day != other.Day:
            return False
        if self.Hour != other.Hour:
            return False
        if self.Minute != other.Minute:
            return False
        if self.Second != other.Second:
            return False
        return True

    def compareNorms(self, other):
        if len(self.NormList) != len(other.NormList):
            return False
        for n1, n2 in itertools.izip(self.NormList, other.NormList):
            if n1.Norm != n2.Norm:
                return False
        return True


class Norm(object):
    """
    Represents a norm value output by an UM EndGame run.
    """
    TOLERANCE = 1.0e-8

    def __init__(self, outer, inner, iterations, norm):
        self.Outer = outer
        self.Inner = inner
        self.Iterations = iterations
        self.Norm = norm

    def __eq__(self, other):
        if self.Outer != other.Outer:
            return False
        if self.Inner != other.Inner:
            return False
        if self.Iterations != other.Iterations:
            return False
        if abs(self.Norm - other.Norm) < self.TOLERANCE:
            return False
        return True


class CompareEGNorms(AnalysisTask):
    """Compare the norms from two EG stdout files."""
    def run_analysis(self):
        """Main analysis routine called from rose_ana."""
        self.process_opt_files()
        self.process_opt_kgo()
        self.process_opt_unmatched()
        # Deal with any unexpected options
        self.process_opt_unhandled()

        # Check that the files to be compared are present
        if len(self.files) != 2:
            raise ValueError("Must specify exactly two files")

        self.perform_comparison()
        self.update_kgo()

    def perform_comparison(self):
        """Compare the norms in the two files."""
        file1 = os.path.realpath(self.files[0])
        file2 = os.path.realpath(self.files[1])

        # Identify which of the two files is the KGO file
        if self.kgo is not None:
            kgo_file = [file1, file2][self.kgo]
            # If this file is missing, no comparison can be performed; it
            # could be that this task is brand new
            if not os.path.exists(kgo_file):
                self.parent.reporter(
                    "KGO File (file {0}) appears to be missing"
                    .format(self.kgo + 1), prefix="[FAIL] ")
                # Note that by exiting early this task counts as failed
                return

        tsList1 = self.extractNorms(file1)
        tsList2 = self.extractNorms(file2)

        prefix = "[INFO] "
        self.parent.reporter(
            '{0} time steps found in input 1'.format(len(tsList1)),
            prefix=prefix)
        self.parent.reporter(
            '{0} time steps found in input 2'.format(len(tsList2)),
            prefix=prefix)

        if (self.allow_unmatched.lower() == "false"
                and len(tsList1) != len(tsList2)):
            prefix = "[FAIL] "
            self.parent.reporter("Number of timesteps different in each file "
                                 "and \"allow_unmatched\" is false",
                                 prefix=prefix)
        else:
            misMatches, num_comps = self.compareTimestepNorms(tsList1, tsList2)

            self.parent.reporter(
                'Compared {0} timesteps'
                .format(num_comps), prefix=prefix)
            if len(misMatches) > 0:
                prefix = "[FAIL] "
                self.parent.reporter(
                    'The following timesteps have different norms:',
                    prefix=prefix)
                for ix1, ix2 in misMatches:
                    self.parent.reporter(
                        'Model time: {0}'.format(str(tsList1[ix1])),
                        prefix=prefix)
            else:
                self.passed = True
                self.parent.reporter(
                    'All matching timesteps have equal norms.', prefix=prefix)

    def process_opt_unmatched(self):
        """Process the files option; a flag to determine if an inconsistent
        number of timesteps are allowed or not."""
        self.allow_unmatched = self.options.pop("allow_unmatched", "false")

        if self.allow_unmatched.lower() not in ("true", "false"):
            raise ValueError("Allow unmatched switch must be true or false")

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

        # Expand environment variables in the filenames and report them
        for ifile, fname in enumerate(files):
            self.parent.reporter(
                "File {0}: {1}".format(ifile + 1,
                                       os.path.realpath(files[ifile])))

    def process_opt_kgo(self):
        """
        Process the KGO option; an index indicating which file (if any) is
        the KGO (Known Good Output) - this may be needed later to assist in
        updating of test results.

        """
        # Get the kgo index from the options dictionary
        kgo = self.options.pop("kgo_file", None)
        # Parse the kgo index
        if kgo is not None:
            if kgo.strip() == "":
                kgo = None
            elif kgo.isdigit():
                kgo = int(kgo)
                if kgo > len(self.files) - 1:
                    msg = "KGO index cannot be greater than number of files"
                    raise ValueError(msg)
            else:
                msg = "KGO index not recognised; must be a digit or blank"
                raise ValueError(msg)
        if kgo is not None:
            self.parent.reporter("KGO is file {0}".format(kgo + 1))
        self.kgo = kgo

    def update_kgo(self):
        """
        Update the KGO database with the status of any files marked by the
        kgo_file option (i.e. whether they have passed/failed the test.)

        """
        if self.kgo is not None and self.parent.kgo_db is not None:
            self.parent.reporter(
                "Adding entry to KGO database (File {0} is KGO)"
                .format(self.kgo + 1), prefix="[INFO] ")
            # Take the KGO file
            kgo_file = self.files[self.kgo]
            # The other file is the suite file
            suite_file = list(self.files)
            suite_file.remove(kgo_file)

            # Set the comparison status
            status = ["FAIL", " OK "][self.passed]

            # Update the database
            self.parent.kgo_db.enter_comparison(
                self.options["full_task_name"],
                os.path.realpath(kgo_file),
                os.path.realpath(suite_file[0]),
                status, "Compared with CompareEGNorms")

    def processTimestepString(self, line, timeStepHeaderPattern,
                              timeStepNumberPattern):
        """
        Process a line of a pe_output file containg a timestep header

        timeStep = processTimestepString(line, timeStepHeaderPattern)
        line: string containing the line to be processed
        timeStepHeaderPattern: string with the pattern of the timestamp
        to be matched
        timeStep: A Timestep object extract from line input argument
        """
        regExOutput = re.finditer(timeStepHeaderPattern, line)
        timeStepStr = ''.join([x1.group() for x1 in regExOutput])
        tsList1 = timeStepStr[11:].lstrip().rstrip().split(' ')
        dateList1 = [int(x1) for x1 in tsList1[0].split('-')]
        timeList1 = [int(x1) for x1 in tsList1[1].split(':')]
        regExOutput2 = re.finditer(timeStepNumberPattern, line)
        tsNum = -1
        try:
            tsNumMatch = [x for x in re.finditer(timeStepNumberPattern,
                                                 line)][0]
            tsNum = int(tsNumMatch.group()[18:].strip(' '))
        except:
            tsNum = -1

        currentTimestep = Timestep(year=dateList1[0],
                                   month=dateList1[1],
                                   day=dateList1[2],
                                   hour=timeList1[0],
                                   minute=timeList1[1],
                                   second=timeList1[2],
                                   number=tsNum)
        return currentTimestep

    def processNormString(self, line, normValuePattern):
        """
        Process a line of a pe_output file containing output model norms

        newNorm = processNormString(line, normValuePattern)
        line: string containing the line to processed
        normValuePattern: string containing the pattern to be matched
        newNorm: a Norm object containing the extracted norm value
        """
        regExOutput = re.finditer(normValuePattern, line)
        normStrRaw = ''.join([x1.group() for x1 in regExOutput])
        normStrRaw = normStrRaw[1:-1].lstrip().rstrip()
        normSet1 = [s1 for s1 in normStrRaw.split(' ') if len(s1) > 0]
        outerVal = int(normSet1[0])
        innerVal = int(normSet1[1])
        iterationsVal = int(normSet1[2])
        normVal = float(normSet1[3])
        newNorm = Norm(outer=outerVal,
                       inner=innerVal,
                       iterations=iterationsVal,
                       norm=normVal)
        return newNorm

    def extractNorms(self, filename):
        """
        Process a file at the specified location

        timeStepList = extractNorms(filename)
        filename: string containing path to the pe_output file
        timeStepList: A list fo Timestep objects extracted from the file
        """
        timeStepHeaderPattern = re.compile(TIMESTEP_HEADER_STRING)
        normHeaderPattern = re.compile(NORM_HEADER_STRING)
        normValuePattern = re.compile(NORM_VALUE_STRING)
        timeStepNumberPattern = re.compile(TIMESTEP_NUMBER_PATTERN)

        timeStepList = []
        with open(filename) as resultFile:
            status = UNKNOWN_SECTION
            currentTimestep = None
            for i, line in enumerate(resultFile):
                if re.findall(timeStepHeaderPattern, line):
                    # check that there is a current timestep, and that it has
                    # at least one norm associated with it. First timestep
                    # in a CRUN output file has a non-zero timestep number,
                    # but no norms because it is an initialisation timestep,
                    # not a calculation step so should be ignored.
                    if (currentTimestep is not None and
                            currentTimestep.NormList):
                        timeStepList += [currentTimestep]

                    currentTimestep = (
                        self.processTimestepString(line,
                                                   timeStepHeaderPattern,
                                                   timeStepNumberPattern))
                    status = TIMESTEP_SECTION
                elif (re.findall(normHeaderPattern, line)
                      and status == TIMESTEP_SECTION):
                    status = NORM_SECTION
                elif (re.findall(normValuePattern, line)
                      and status == NORM_SECTION):
                    newNorm = self.processNormString(line, normValuePattern)
                    currentTimestep.addNorm(newNorm)

        return timeStepList

    def compareTimestepNorms(self, tsList1, tsList2):
        """
        Compare 2 lists of Timestep objects
        The function considers all possible pairs from the 2 lists. If a pair
        of Timestep objects have the same timestamp, the norms for that pair
        are compared and if they differ, the indices of the each Timestep
        object is stored. The function returns a list of Tuples containing 2
        integers, which refer to a Timestep in the first and second input
        arguments that have equal timestamps and unequal norms

        misMatches = compareTimestepNorms(tsList1, tsList2)
        tsList1: A list of Timestep objects
        tsList2: A list of Timestep objects
        misMatches: a list of Tuples containing 2 integers, which refer to a
                    Timestep in the first and second input arguments that have
                    equal timestamps and unequal norms
        """
        misMatches = []
        num_comps = 0
        iter1 = ((ix1, ts1, ix2, ts2)
                 for ix1, ts1 in enumerate(tsList1)
                 for ix2, ts2 in enumerate(tsList2))
        for ix1, ts1, ix2, ts2 in iter1:
            # We are not comparing norms for timestep 0, because none are
            # output as no calculation has been done!
            if ts1.Number > 0 and ts2.Number > 0 and ts1.compareTimes(ts2):
                num_comps += 1
                if not ts1.compareNorms(ts2):
                    misMatches += [(ix1, ix2)]
        return misMatches, num_comps
