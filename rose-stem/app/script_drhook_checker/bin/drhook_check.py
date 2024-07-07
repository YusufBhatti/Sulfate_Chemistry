#!/usr/bin/env python
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
'''drhook_check.py <source> [<source>, ...]

   Given a file or directory, update all DrHook calls in the
   resulting tree to the recommended standard of
   [<ModuleName>//':'//]<RoutineName>

   Input sources can be files or directories but should be file system
   paths, not Subversion repository locations.

   Method:
    Each Fortran file is read in turn.
    Any #includes are replaced by the contents of the relevant include file.
    For each file that contains a call to dr_hook:
    - if it contains a module:
        Check that any zhook* variables have the correct attributes
        Find a call that contains a <modulename> variable;
        Check all calls contain this variable;
        Check the declaration of <modulename> itself.
    - for each routine in the file that calls dr_hook:
        Find the names of all <zhook_handle> variables
        Check that all zhook* variables have the correct attributes
        Find a call that contains a <routinename> variable;
        Check all calls contain this variable;
        Check the declaration of <routinename> itself;
        Check the calls themselves are sensible wrt in/out and no. of calls.'''

import os
import re
import sys


def get_args():
    '''Parse arguments'''
    if len(sys.argv) < 2:
        sys.exit("Provide at least one argument: a file or directory path")
    else:
        return sys.argv[1:]


def get_matching_files(targets, file_pattern):
    '''Return a list "files" of all files in "targets"
       matching pattern "file_pattern".'''

    files = []  # List of all files to process
    for source in targets:

        if os.path.isfile(source) and file_pattern.search(source):
            # Append the single file and move on to the next source:
            files.append(source)
            continue

        # Process directories:
        for root, _, fnames in os.walk(source):
            for fname in fnames:
                if file_pattern.search(fname):
                    files.append(os.path.join(root, fname))

    return files


def get_all_files(targets):
    '''Find all files in the provided targets, as two lists
       (Fortran files and include files).'''

    # Fortran files:
    ffile_pattern = re.compile("\.(F|f)(\d{2})?$")
    ffiles = get_matching_files(targets, ffile_pattern)

    # Include files:
    ifile_pattern = re.compile("\.h$")
    ifiles = get_matching_files(targets, ifile_pattern)

    return ffiles, ifiles


def read_file(ffile):
    '''Read file and return contents'''
    with open(ffile, 'r') as fhandle:
        lines = fhandle.readlines()
    return lines


def insert_includes(lines, ifiles):
    '''Insert any includes files into the read Fortran file.
       Each new include is in turn checked for further includes.'''
    warnings = []

    # A while loop is needed as we modify the list being iterated over.
    # The counter is only incremented when:
    # a) the line does not contain an include
    # b) we have failed to replace an include.
    # When an include is replaced the counter is NOT incremented so that
    # the first line of the include is also checked for includes.
    i = 0
    while i < len(lines):
        # Check for an include file:
        result = re.match(r'\s*#include\s*(["\'])(\w+\.h)\1', lines[i])
        if result:
            include_file = result.group(2)
            # Was the indicated include file given to the checker?
            for ifile in ifiles:
                if os.path.basename(ifile) == include_file:
                    # Read the include file:
                    ilines = read_file(ifile)
                    # Use array slicing to replace element i (the line
                    # containing the #include) with the include file contents:
                    lines[i:i+1] = ilines
                    break

            # Did the file get replaced?
            if re.match(r'\s*#include\s*(["\']){0}\1'.format(include_file),
                        lines[i]):
                # The (same) #include is still there.
                # Report an error and move on to the next line.
                warnings += ['Include file {0}'.format(include_file)
                             + ' could not be found']
                i += 1
        else:
            # Move on to the next line.
            i += 1

    return lines, warnings


def test_for_drhook(lines):
    '''Determine if the provided lines contain a call to Dr Hook'''
    for line in lines:
        if re.match(r'(?i)((\s*IF\s*\(\s*lhook\s*\))?\s*CALL\s+dr_hook)',
                    line):
            return True
    return False


def concatenate_lines(lines):
    '''Roll multi-line statements into a single line, after stripping out
       OpenMP sentinels and comments.'''
    # Convert the file into one long string
    text = ''.join(lines)
    # Remove C pre-processing directives
    text = re.sub(r'(?m)^\s*#.*$', '', text)
    # Strip OpenMP sentinels:
    text = re.sub(r'(?im)^\s*!\$(OMP)?\s*&?', '', text)
    # Strip comments. This may find a '!' inside a string, but if so
    # it shouldn't be anywhere near any Dr Hook-related code.
    text = re.sub(r'(?m)!.*$', r'\n', text)
    # Strip blank lines:
    text = re.sub(r'(?m)^\s*\n', '', text)
    # Merge continuation lines:
    text = re.sub(r'(?m)\s*&\s*$\s*', '', text)
    # Turn the string back into separate lines:
    lines = text.split('\n')
    return lines


def get_mod_name(lines):
    '''Return the name of the parent module'''
    for line in lines:
        result = re.match(r'\s*(?i)MODULE\s+(\w+)', line)
        if result:
            return result.group(1).upper()


def get_routine_names(lines):
    '''Return the names (uppercase) of all subroutines or functions declared'''
    routines = []
    for line in lines:
        # This has to be flexible enough to match things like
        # PURE RECURSIVE SUBROUTINE or
        # ELEMENTAL RECURSIVE CHARACTER(LEN=10) FUNCTION
        if ((re.match(r'\s*(?i)(PROGRAM|(\w+\s+){0,2}SUBROUTINE)\s+\w+', line)
           or re.match(r'\s*(?i)((\w+(\(.*\))?\s+){0,3}FUNCTION)\s+\w+', line))
           and not re.match(r'\s*END\s*', line)):
            result = re.search(
                r'\s*(?i)(?:PROGRAM|SUBROUTINE|FUNCTION)\s+(\w+)', line)
            routine = result.group(1).upper()
            routines.append(routine)
    return routines


def check_module(module, lines):
    '''Check a module for Dr Hook-related errors'''
    warnings = []
    modulevar = None
    string_literal = False

    warnings += check_zhook_declarations(lines, module)

    # Search for a call in the <name>//':'//<name> pattern.
    for line in lines:
        if test_for_drhook([line]):
            # Check for raw strings first
            result = re.search(r'(?i)(dr_hook\s*\(\s*(\'|"))', line)
            if result:
                string_literal = True
            else:
                result = re.search(
                    r'(?i)(dr_hook\s*\(\s*(\w+)\s*//\s*'
                    '(\'|")\s*::?\s*(\'|")\s*//\w+)', line)
                if result:
                    # We've found a variable; store it and stop looking:
                    modulevar = result.group(2)
                    break

    if modulevar:
        # If we found a <modulename> variable, do further checks...
        warnings += check_module_variable(module, modulevar, lines)
        warnings += check_module_calls(module, modulevar, lines)
    elif string_literal:
        # ...otherwise if we only found literals issue a warning...
        warnings += ['Calls to dr_hook contains strings in place of variables '
                     'for module name']
    else:
        # ...else we found nothing suitable at all.
        warnings += ['Calls to dr_hook do not include name of module']

    return warnings


def check_module_calls(module, modulevar, lines):
    ''' Check all calls to Dr Hook include a common <modulename> variable.'''
    warnings = []
    different_variable = False
    no_modulename = False
    string_literal = False
    for line in lines:
        if test_for_drhook([line]):
            # Check for raw strings first:
            if ((re.search(r'(?i)(dr_hook\s*\(\s*(\'|"))', line)
               and modulevar not in line)
               or re.search(r'(?i)(\'|")({0})'.format(module), line)):
                # Starts with a literal and doesn't contain the variable;
                # or contains a literal starting with the actual name.
                string_literal = True
            else:
                result = re.search(r'(?i)(dr_hook\s*\(\s*(\w+)//\s*(\'|")\s*'
                                   '::?\s*(\'|")\s*//\w+)', line)
                if result:
                    var = result.group(2)
                    if var.lower() != modulevar.lower():
                        different_variable = True
                else:
                    no_modulename = True

    if string_literal:
        warnings += ['Call to dr_hook contains a string literal '
                     'in place of module name']
    if different_variable:
        warnings += ['Multiple variables used to store name of module '
                     'in calls to dr_hook']
    if no_modulename:
        warnings += ['Call to dr_hook does not include name of module']

    return warnings


def check_module_variable(module, modulevar, lines):
    '''Check that <modulevar> is correctly declared'''
    # Check a) defined in module (not subroutine), b) PRIVATE,
    # c) assumed length parameter, d) value matches name of module,
    # e) correct case.
    warnings = []
    module_header = True
    is_public = True
    result = None
    for line in lines:

        if re.match(r'\s*(?i)(CONTAINS)(\W|$)', line):
            module_header = False
        if re.match(r'\s*(?i)(PRIVATE)(\W|$)', line):
            is_public = False

        result = re.match(
            r'\s*(?i)(CHARACTER(.*)::\s*{0}\s*=\s*(\'|")(\w+)(\'|"))'
            .format(modulevar), line)
        if result:

            if not re.search(r'(?i)(LEN\s*=\s*)?\*', result.group(2)):
                warnings += [modulevar + ' is not an assumed length parameter']

            if not result.group(4).upper() == module.upper():
                warnings += ['Value of variable {0} '.format(modulevar)
                             + 'does not match module name {0}'.format(module)]

            if not result.group(4).isupper():
                warnings += ['Value of variable {0} '.format(modulevar)
                             + 'is not all uppercase']

            if module_header:
                if re.search(r'(?i)\WPRIVATE\W', result.group(2)):
                    is_public = False
                elif re.search(r'(?i)\WPUBLIC\W', result.group(2)):
                    warnings += [modulevar + ' should not have '
                                 'the PUBLIC attribute']
            else:
                warnings += [modulevar + ' should be declared '
                             'in the module header, not a routine']
                is_public = False  # Don't check routine variables for PRIVATE

            break    # We only need to reach here once; then stop looking.

    if is_public and result is not None:
        warnings += ['Variable {0} should be declared PRIVATE'
                     .format(modulevar)]

    return warnings


def get_routine_limits(lines, routine, start):
    '''Return the start and end point of a routine'''
    # Start looking the line after the previous routine's start point so that
    # we can cope with multiple routines of the same name in a file
    # (e.g. because of ifdefs).
    # While it'd be fractionally more efficient to begin at the end of the
    # previous routine, it's much safer to do it this way.
    start = start + 1
    end = -1
    drhook_enabled = False
    for num, line in enumerate(lines[start:]):
        # Trailing horribleness in matches is to avoid multiple matches in
        # a file containing e.g. mysub and mysub1. See also get_routine_names.
        if ((re.match(r'\s*(?i)(PROGRAM|(\w+\s+){{0,2}}SUBROUTINE)\s+({0})'
                      '(\s|\(|&|$)'.format(routine), line) or
             re.match(r'\s*(?i)(\w+(\(.*\))?\s+){{0,3}}FUNCTION\s+({0})'
                      '(\s|\(|&|$)'.format(routine), line))
           and not re.match(r'\s*END\s*', line)):
            start = start + num
            break

    for num, line in enumerate(lines[start:]):
        # The CONTAINS is to match subroutines that contain a function
        # before their END SUBROUTINE line (which is why we need two loops).
        # We also need to stop searching when Dr Hook is disabled, making sure
        # it has been enabled first (some programs may ensure it has been
        # disabled before enabling it later.)
        if (re.match(r'\s*CALL\s*drhook_control_enable\W', line)):
            drhook_enabled = True

        match_contains = re.match(r'\s*(?i)CONTAINS(\W|$)', line)
        match_end = re.match(
            r'\s*END\s+(PROGRAM|SUBROUTINE|FUNCTION)\s+(?i)({0})(\s|$)'
            .format(routine), line)
        match_disable = re.match(r'\s*CALL\s*drhook_control_disable\W', line)

        if (match_contains
           or match_end
           or (drhook_enabled and match_disable)):
            end = start + num
            break

    if end == -1:
        # Routine ends with e.g. END SUBROUTINE and not END SUBROUTINE <NAME>
        for num, line in enumerate(lines[start:]):
            if re.match(r'\s*(?i)CONTAINS(\W|$)', line) or re.match(
               r'\s*END\s+(PROGRAM|SUBROUTINE|FUNCTION)(?i)(\s|$)',
               line):
                end = start + num
                break

    return start, end


def check_routine(routine, lines, module):
    '''Check a subroutine/function for Dr Hook-related errors.'''
    warnings = []
    string_literal = False
    routinevar = None

    warnings += check_zhook_declarations(lines)

    for line in lines:
        if test_for_drhook([line]):
            # Search for a routinename variable, stopping on success.

            # First, check if last part of the name is a raw string:
            result = re.search(r'(?i)(dr_hook\s*\(.*?(\'|")\w+\2\s*,)', line)
            if result:
                string_literal = True
                # ...and keep checking the remaining calls.
            elif module:
                # First, try looking for the second name...
                result = re.search(r'(?i)dr_hook\s*\(.*?:\s*(\'|")\s*//'
                                   '\s*(\w+)', line)
                if result:
                    routinevar = result.group(2)
                    break
                else:
                    # ...and if the module name is missing, try the first:
                    result = re.search(r'(?i)dr_hook\s*\(\s*(\w+)', line)
                    if result:
                        routinevar = result.group(1)
                        break
            else:
                # Take the first word:
                result = re.search(r'(?i)(dr_hook\s*\(\s*(\w+))', line)
                if result:
                    routinevar = result.group(2)
                    break

    if routinevar:
        # If we found a <routinename> variable, do further checks...
        warnings += check_routine_declarations(routine, routinevar, lines)
        warnings += check_routine_calls(routinevar, lines, module)
    elif string_literal:
        # ...otherwise if we only found literals issue a warning...
        warnings += ['Calls to dr_hook contains strings in place of variables '
                     'for routine name']
    else:
        # ...else we found nothing suitable at all.
        warnings += ['Calls to dr_hook do not match any recognised pattern']

    return warnings


def check_routine_declarations(routine, routinevar, lines):
    '''Check the declaration of <routinename> is valid.'''
    # Check a) assumed length parameter, b) value matches name of routine,
    # c) correct case.
    warnings = []
    result = None
    split_declaration = False
    for line in lines:
        result = re.match(r'\s*(?i)(CHARACTER(.*)::\s*{0}\s*='
                          '\s*(\'|")(\w+)(\'|"))'.format(routinevar), line)

        if result:
            if not re.search(r'(?i)(LEN\s*=\s*)?\*', result.group(2)):
                warnings += ['Variable {0} is not an '
                             'assumed length parameter'.format(routinevar)]

            if not result.group(4).upper() == routine.upper():
                warnings += ['Value of variable {0} does not match '
                             'routine name {1}'.format(routinevar, routine)]

            if not result.group(4).isupper():
                warnings += ['Value of variable {0} is not all uppercase'
                             .format(routinevar)]
            break  # Stop looking

        else:
            # Declaration may be split over two lines, with a separate
            # "PARAMETER (var = value)" statement.
            result = re.match(r'\s*(?i)(CHARACTER(.*){0})'
                              .format(routinevar), line)
            if result:
                split_declaration = True

                if not (re.search(r'(?i)(LEN=\*)', result.group(2))
                        or re.search(r'\(\*\)', result.group(2))):
                    warnings += ['Variable {0} is not an '
                                 'assumed length parameter'.format(routinevar)]

                break  # Stop looking.

    if not result:
        # Didn't find anything. Unlikely to compile, but anyway...
        warnings += ['Variable {0} is not declared inside routine {1}'
                     .format(routinevar, routine)]

    if split_declaration:
        # We have to scan the rest of the file to find the initialisation.
        for line in lines:
            result = re.match(
                r'\s*(?i)(PARAMETER\s*\(.*{0}\s*=\s*(\'|")(\w+)(\'|"))'
                .format(routinevar), line)
            if result:
                if not result.group(3).upper() == routine.upper():
                    warnings += ['Value of variable {0} does not match routine'
                                 ' name {1}'.format(routinevar, routine)]
                if not result.group(3).isupper():
                    warnings += ['Value of variable {0} is not all uppercase'
                                 .format(routinevar)]

    return warnings


def check_routine_calls(routinevar, lines, module):
    '''Check all the dr_hook calls within a routine.'''
    # Check a) all calls include the same <routinename> variable,
    # b) first call is 'in', c) last call is 'out',
    # d) all calls prefaced with "IF (lhook)",
    # e) all RETURN statements after the first dr_hook call
    #    are preceded by a dr_hook (out) call,
    # f) the subroutine ends with a dr_hook call.
    warnings = []
    different_variable = False
    no_routinename = False
    string_literal = False
    drhook_active = False   # Have we called Dr Hook in this routine yet?
    no_of_calls = 0
    first_call_error = False
    no_lhook_test = False
    colons = False
    whitespace = False
    return_without_call = False
    return_on_zhook_in = False
    end_without_call = False
    for num, line in enumerate(lines):
        var = None
        if test_for_drhook([line]):
            no_of_calls += 1
            lastline = line  # Save a copy of the most recent call...
            lastline_num = num  # ...and remember which line it was

            # Check the first call:
            if not drhook_active:
                drhook_active = True  # First call has been reached
                if re.search('(?i)(zhook_out)', line):
                    first_call_error = True

            # Check for missing IF (lhook) tests
            if re.match(r'(?i)\s*CALL\s+dr_hook', line):
                no_lhook_test = True

            # Check for wrong number of colons:
            if re.search(r'//\s*(\'|")\s*::\s*(\'|")\s*//', line):
                colons = True

            # Check for unnecessary whitespace around the separator:
            if (re.search(r'//\s*(\'|")\s+:', line)
               or re.search(r':\s+(\'|")\s*//', line)):
                whitespace = True

            # Check for raw strings
            if (re.search(r'(?i)(dr_hook\s*\(.*(\'|")\s*,)', line)
               and routinevar not in line):
                # Ends with a literal & doesn't contain the expected variable.
                # Chances of using a different variable AND adding an
                # extra sub-identifier to the string are low, so:
                string_literal = True
            else:
                result = re.search(r'(?i)(dr_hook\s*\((.*?)\s*,)', line)
                if result:
                    names = re.split('\W+', result.group(2))  # List of words
                    if module:
                        # Take the second word if it exists, otherwise
                        # take the first. (Accounts for 3-part strings such as
                        # ModuleName//':'//RoutineName//':subclass'.)
                        if len(names) > 1:
                            var = names[1]
                        else:
                            var = names[0]
                    else:
                        # Take the first matching word
                        var = names[0]

                    if var.lower() != routinevar.lower():
                        different_variable = True

        elif re.match(r'(?i)\s*RETURN(\W|$)', line) and drhook_active:
            # Found a RETURN after Dr Hook was activated in this routine.
            # Check the preceding line for Dr Hook:
            if not test_for_drhook([lines[num-1]]):
                return_without_call = True
            # Check it's an 'out' call:
            elif re.search('(?i)zhook_in', lines[num-1]):
                faulty_return_num = num-1
                return_on_zhook_in = True

        elif (re.match(r'(?i)\s*IF\s*\(.*\)\s*RETURN(\W|$)', line) and
              drhook_active):
            # A one line IF (something) RETURN statement (after Dr Hook has
            # been activated) automatically fails:
            return_without_call = True

    # Test the last line for either a RETURN (which is then covered by
    # the tests above) or a call to dr_hook:
    if not (re.match(r'(?i)\s*RETURN(\W|$)', lines[-1])
       or test_for_drhook([lines[-1]])):
        end_without_call = True

    # If only one call then say that, rather than having a slightly odd message
    # about the first or last call being invalid.
    if no_of_calls == 1:
        warnings += ['Routine calls dr_hook once; '
                     'at least two calls (or none) are required']
    else:
        if first_call_error:
            warnings += ["First call to dr_hook passes 'zhook_out', "
                         "not 'zhook_in'"]
        # Test the last call:
        if re.search('(?i)(zhook_in)', lastline):
            warnings += ["Last call to dr_hook passes 'zhook_in', "
                         "not 'zhook_out'"]

    if string_literal:
        warnings += ['Call to dr_hook contains a string literal '
                     'in place of routine name']
    if different_variable:
        warnings += ['Multiple variables used to store name of routine '
                     'in calls to dr_hook']
    if no_routinename:
        warnings += ['Call to dr_hook does not include name of routine']
    if no_lhook_test:
        warnings += ["Call to dr_hook is not protected by 'IF (lhook)'"]

    if colons:
        warnings += ['Call to dr_hook contains too many colons as a separator']
    if whitespace:
        warnings += ['Call to dr_hook contains unnecessary whitespace around '
                     'colon separator']

    # These tests may require special care in the top-most routine that uses
    # Dr Hook, e.g. checking against drhook_control_disable or catering for
    # other program shutdown functions such as appterminate.
    if return_without_call:
        warnings += ['A RETURN statement is not preceded by a call to '
                     'dr_hook']
    if end_without_call:
        warnings += ['Routine ends without a call to dr_hook as the '
                     'last action']

    # Extra check on line number here to ensure the last call (which may or may
    # not come with a RETURN statement) doesn't trigger this warning too:
    if return_on_zhook_in and (lastline_num != faulty_return_num):
        warnings += ["A RETURN statement is preceded by a call to dr_hook "
                     "which passes 'zhook_in', not 'zhook_out'"]

    return warnings


def check_zhook_declarations(lines, module=None):
    '''Check the zhook variables are declared correctly'''

    warnings = []
    is_public = True
    result = None
    handles = dict()  # A dict. so we have mark each one as declared or not.
    for line in lines:

        if module:
            if re.match(r'\s*(?i)(PRIVATE)(\W|$)', line):
                is_public = False
            if re.match(r'\s*(?i)(CONTAINS)(\W|$)', line):
                break  # Stop at the first CONTAINS for modules.

        # Get the names of the "handle" variables (usually zhook_handle)
        # We only need to look for these inside the routine where they're used.
        else:
            result = re.search(r'(?i)dr_hook\s*\((.*?,){2}(\w+)', line)
            if result:
                handles[result.group(2)] = False  # False until known declared

        zhook_result = re.match(r'(?i)\s*INTEGER(.*)::\s*(zhook_(in|out))'
                                '\s*=\s*(\d)', line)
        if zhook_result:
            attributes = zhook_result.group(1)
            zhook = zhook_result.group(2)
            inout = zhook_result.group(3).lower()
            value = zhook_result.group(4)

            if inout == 'in' and value != '0':
                warnings += ['Value of zhook_in is incorrect']
            if inout == 'out' and value != '1':
                warnings += ['Value of zhook_out is incorrect']
            if not re.search(r'(?i)\WPARAMETER\W', attributes):
                warnings += [zhook + ' is not a parameter']
            if not re.search(r'(?i)KIND\s*=\s*jpim', attributes):
                warnings += [zhook + " is not of kind 'jpim'"]

    # Second loop to check:
    # 1) handle variables are declared correctly;
    # 2) zhook_in/out are PRIVATE if declared at module level.
    for line in lines:

        for handle in handles:
            result = re.match(r'(?i)\s*REAL(.*){0}'.format(handle), line)
            if result:
                attributes = result.group(1)
                handles[handle] = True
                if not re.search(r'(?i)KIND\s*=\s*jprb', result.group(1)):
                    warnings += [handle + "is not of kind 'jprb'"]

        if module and is_public:
            if re.match(r'\s*(?i)(CONTAINS)(\W|$)', line):
                break  # Stop at the first CONTAINS for modules.

            zhook_result = re.match(r'(?i)\s*INTEGER(.*)::\s*(zhook_(in|out))'
                                    '\s*=\s*\d', line)
            if zhook_result:
                attributes = zhook_result.group(1)
                zhook = zhook_result.group(2)

                if not re.search(r'(?i)\WPRIVATE\W', attributes):
                    warnings += [zhook + ' should be PRIVATE '
                                 'or declared inside each routine']
                if re.search(r'(?i)\WPUBLIC\W', attributes):
                    warnings += [zhook + ' should not have '
                                 'the PUBLIC attribute']

    for handle, is_declared in handles.iteritems():
        if not is_declared:
            warnings += [handle + ' is not declared inside the routine']

    return warnings


def print_warnings(ffile, warn_file, warn_module, warn_routine):
    '''Print the collection of warnings for this file.'''
    print
    print 'In file {0}:'.format(ffile)
    if warn_file:
        for ffile, warnings in warn_file.iteritems():
            for warning in warnings:
                print '  ' + warning
    if warn_module:
        for module, warnings in warn_module.iteritems():
            print '  In module {0}:'.format(module.lower())
            for warning in warnings:
                print '    ' + warning
    if warn_routine:
        for routine, warnings in warn_routine.iteritems():
            print '  In routine {0}:'.format(routine.lower())
            for warning in warnings:
                print '    ' + warning
    return


def main():
    '''Main program - open a source tree, check files and report problems.'''

    success = True
    targets = get_args()
    ffiles, ifiles = get_all_files(targets)

    if len(ffiles) == 0:
        sys.exit("No Fortran files found in {0}".format(targets))

    for ffile in ffiles:
        warn_file = dict()
        warn_module = dict()
        warn_routine = dict()
        lines = read_file(ffile)

        lines, warnings = insert_includes(lines, ifiles)
        # Only populate the dictionary if there's actually a warning:
        if warnings:
            warn_file[ffile] = warnings

        lines = concatenate_lines(lines)
        file_calls_drhook = test_for_drhook(lines)

        if file_calls_drhook:
            module = get_mod_name(lines)         # Module name or None.
            routines = get_routine_names(lines)  # Subs & functions.

            if module:
                # Only populate the dictionary if there's actually a warning.
                warnings = check_module(module, lines)
                if warnings:
                    warn_module[module] = warnings

            start = -1
            for routine in routines:
                start, end = get_routine_limits(lines, routine, start)
                routine_calls_drhook = test_for_drhook(lines[start:end])

                if routine_calls_drhook:
                    # Only populate the dict. if there's a warning.
                    warnings = check_routine(routine, lines[start:end], module)
                    if warnings:
                        warn_routine[routine] = warnings

        if warn_file or warn_module or warn_routine:
            print_warnings(ffile, warn_file, warn_module, warn_routine)
            success = False

    if success:
        print "No errors were detected."
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
