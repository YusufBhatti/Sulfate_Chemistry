#!/usr/bin/env python
# *********************************COPYRIGHT************************************
# (C) Crown copyright Met Office. All rights reserved. 
# For further details please refer to the file COPYRIGHT.txt 
# which you should have received as part of this distribution. 
# *********************************COPYRIGHT************************************
"""
Script which tests the declarations of namelist-read broadcasting in
the UM source code
"""
import re
import os
import sys
import itertools

from optparse import OptionParser
from textwrap import wrap

# Low-limit for primes (i.e. the first prime returned must be
# above this number)  I think this ought to be set *just* higher
# than the largest likely scalar value you might expect a user
# to put into the expressions, but realistically let's just put
# it as something decently high like 1000 so that it definitely
# won't get confused with any scalar values
_MINIMUM_PRIME_LIMIT = 1000

# Global debugging output switch
_DEBUG = False

# Name mappings which will be used later
_VARTYPES = [("int",   "INTEGER"),
             ("real",  "REAL"),
             ("log",   "LOGICAL"),
             ("chars", "CHARACTER")]

# Desired maximum column width for output - we make an exception
# for filenames, which are always printed on a single line to aid
# ease of selection by the user
_OUTPUT_LINE_WIDTH = 80

# Pattern which will match the intended input files
_FILE_PATTERN = r".*\.(F|f)90"

#-------------------------------------------------------------------------------
def banner_print(message, maxwidth=_OUTPUT_LINE_WIDTH, char="%"):
    """
    Simple routine which prints a banner message
    """
    wrap_message = [char + " " + elt.ljust(maxwidth-4) + " " + char
                    for elt in wrap(message, width=maxwidth-4)]
    print "\n{0:s}\n{1:s}\n{0:s}".format(char*maxwidth, "\n".join(wrap_message))
    
#-------------------------------------------------------------------------------
def gen_primes():
    """
    Generator which returns prime numbers (once initialised the generator will
    return the next prime number when it is passed to the "next" intrinsic)

    The algorithm in use is based on the "Sieve of Eratosthenes"; modified so
    that instead of returning the first N primes it is unbounded, and it will
    only return primes above a particular starting number
    """
    primedict = {}  
    prime = 2  
    while True:
        if prime not in primedict:
            if prime > _MINIMUM_PRIME_LIMIT:
                yield prime
            primedict[prime * prime] = [prime]
        else:
            for found in primedict[prime]:
                primedict.setdefault(found + prime, []).append(found)
            del primedict[prime]
        prime += 1

#-------------------------------------------------------------------------------
def create_prime_dict(expr):
    """
    Given an expression, finds variable names and generates a replacement
    dictionary which maps each unique name onto a prime number
    """
    # Setup generator
    prime_gen = gen_primes()
    
    # Find all variable names - to do this get all unique words in the
    # expression, and then filter out digits and sort the result
    words = re.findall("\w+", expr)
    variable_names = sorted([word for word in words if not word.isdigit()])

    # Create an equivalent list of primes (strings), one for each variable name
    prime_list = [str(next(prime_gen)) for name in variable_names]

    # Turn into a dictionary and report its contents
    replace_dict = dict(zip(variable_names, prime_list))
    if len(replace_dict) > 0:
        output  = "Associating the following primes with variable names:"
        for name, prime in replace_dict.items():
            output += "\n   {0:s}: {1:s}".format(name, prime) 
    else:
        output = "No variables used, no prime substitution required"

    # Construct a dictionary out of the two lists - the variable names are
    # the keys and the primes are the values
    return replace_dict, output

#-------------------------------------------------------------------------------
def expr_eval(expr, replace_dict):
    """
    Given an expression and a replacement dictionary performs replacement
    of variable names with the number representing them, and then evaluates
    them to a numeric result
    """

    # If the expression is blank, return
    if expr == "":
        return 0, "    0 = 0"

    # Perform the replacements
    new_expr = expr
    for name, prime in replace_dict.items():
        new_expr = new_expr.replace(name, prime)

    # If this has worked we should be able to evaluate the expression to get
    # a numerical result, which will be returned
    try:
        result = eval(new_expr)
    except:
        raise ValueError("Failure to evaluate expression: {0:s}".format(expr))

    if len(replace_dict) > 0:
        output = "   {0:s} = {1:s} = {2:d}".format(expr, new_expr, result)
    else:
        output = "   {0:s} = {1:d}".format(expr, result)

    return result, output

#-------------------------------------------------------------------------------
def read_and_preprocess_file(filename):
    """
    Reads in the source code from a given file, and pre-processes it
    to remove all comments and strings, and collapse continuation lines;
    to make subsequent parsing of the contents easier
    """
    with open(filename, "r") as input_file:
        # Read in the entire file source as one big string
        code = input_file.read()

        # Cut the content of all quoted strings of both types, since these
        # could confuse any further parsing
        code = re.sub(r"((?P<q1>[\"']).*?(?P=q1))", "\"\"", code)

        # Remove all comments (note that this regex will not match past a
        # newline character, so this serves to strip all comments and
        # trailing comments)
        code = re.sub(re.compile(r"!.*$", re.MULTILINE), "", code)

        # Next get rid of all lines which are now completely blank, to
        # make the final result compact but also to deal with cases where
        # a comment line appears before a continuation line
        code = re.sub(re.compile(r"^\s*\n", re.MULTILINE), " ", code)

        # Collapse the continuation lines (The regex will match whitespace
        # at the end of a line followed by a continuation character, possibly
        # some more space, the newline, and then any leading whitespace on the
        # following line
        code = re.sub(re.compile(r"\s*&\s*$\s*", re.MULTILINE), " ", code)

    return code

#-------------------------------------------------------------------------------
def get_type_counts(read_nml_code, n_types_name, count_names):
    """
    Extract the expressions giving the count of the number of types from
    the setup part of the namelist broadcast routine
    """
    # Get the number of types counter
    search = re.search(
            r"""
            \s*::\s*            # Double-colon, possibly with spaces
            {0:s}\s*=\s*        # Type-count variable (substituted in)
            (?P<expr>.*)        # Expression setting it (captured as "expr")
            """.format(n_types_name),
            read_nml_code, flags=re.IGNORECASE|re.VERBOSE)

    # This should exist, abort for this routine if it doesn't
    if search:
        num_types = int(search.group("expr"))
    else:
        msg = "ERROR: Unable to find type count {0:s}".format(n_types_name)
        return False, msg

    # Get the different variable count expressions
    variable_counts = {}
    for vartype, _ in _VARTYPES:

        # If no name exists, save a zero and continue
        if not count_names.has_key(vartype):
            variable_counts[vartype] = "0"
            continue
        
        # Otherwise look for the variable count line, saving the expression
        # on the RHS (if the search matches)
        search = re.search(
            r"""
            \s*::\s*        # Double-colon, possibly with spaces
            {0:s}\s*=\s*    # Type instance count variable (substituted)
            (?P<expr>.*)    # Expression setting it (captured as "expr")
            """.format(count_names[vartype]),
            read_nml_code, flags=re.IGNORECASE|re.VERBOSE)

        # If an expression was found save it against the variable type
        # in the dictionary, otherwise this is an error
        if search:
            variable_counts[vartype] = search.group("expr")
        else:
            msg = ("ERROR: Unable to find count expression "
                   "{0:s}".format(count_names[vartype]))
            return False, msg

    # We can now check that the number of types matches
    num_types_found = sum([count != "0"
                           for count in variable_counts.values()])
    if num_types_found > num_types:
        msg = ("ERROR: More \"n_<type>\" variables found than "
               "indicated by \"no_of_types\" count")
        return False, msg
    elif num_types_found < num_types:
        msg = ("ERROR: Not enough \"n_<type>\" variables found "
               "based on \"no_of_types\" count")
        return False, msg
    
    return True, variable_counts

#-------------------------------------------------------------------------------
def get_definitions(broadcast_code):
    """
    Extract the expressions giving the variable definition information of
    each variable type from the broadcast type definition
    """
    # Initialise dictionary to hold definition lists, these are
    # initialised to empty lists as we will be appending to them
    variable_definitions = {}
    variable_names = []
    for vartype, _ in _VARTYPES:
        variable_definitions[vartype] = []

    # Go through the code to find the dimensions and name of each
    # variable for each of the types
    for vartype, varname in _VARTYPES:
        if vartype != "chars":

            # For non-character variables, this statement will find
            # all variable names of a particular type, including any
            # bracketed dimensioning
            defs = re.findall(
                r"""
                (?:^|\s+)  # Start of line or at least one space
                           # (non-capturing)
                {0:s}.*::  # Variable name with double-colon
                \s*(.*)    # Possible space followed by variable
                           # definition content in-full (captured)
                """.format(varname),
                broadcast_code, flags=re.IGNORECASE|re.VERBOSE)

            for val in defs:
                # Now extract the dimensions if this is an array
                search = re.search(
                    r"""
                    (?P<name>\w+)\s*    # The variable name (captured) possibly
                                        # followed by some spaces
                    (\((?P<dims>.*)\)|) # An optional bracketed expression
                                        # containing dimensions (captured)
                    """, val, flags=re.VERBOSE)
                if search:
                    # Add the name of the variable to the list
                    variable_names.append(search.group("name").lower())
                    if search.group("dims") is not None:
                        # Build up an expression with these using
                        # multiplication (if multiple dims are found)
                        dim_expr = search.group("dims").split(",")
                        for iexp, expr in enumerate(dim_expr):
                            # If the dimension has bounds specified then the
                            # result should give the upper bound minus the
                            # lower bound (+1 for fortran indexing)
                            if ":" in expr:
                                dim_expr_split = expr.split(":")
                                dim_expr[iexp] = (
                                    "("+dim_expr_split[1] + "-"+
                                    dim_expr_split[0] + "+1)")
                            # The expression ought to be bracketed to make
                            # sure any arithmetic is handled correctly
                            dim_expr[iexp] = "(" + dim_expr[iexp] + ")"

                        variable_definitions[vartype].append((val,
                            "*".join(dim_expr)))
                    else:
                        # If this isn't an array the dimension is just 1
                        variable_definitions[vartype].append((val, "1"))
        else:
            # Character variables are treated differently because they
            # have the LEN= dimension to consider as well.  This
            # statement will return a list of tuples with the LEN value
            # as well as the variable name (including bracketed
            # dimensioning)
            defs = re.findall(
                r"""
                (?:^|\s+)       # Start of line or at least one space
                                # (non-capturing)
                {0:s}           # The variable name
                .*              # Match anything here, until...
                LEN\s*=\s*(\w+) # The length specification (capturing
                                # the expression used for the length)
                .*::\s*         # Match anything again until the 
                                # double-colon
                (.*)            # And finally the variable definition
                                # content in-full (captured)
                """.format(varname),
                broadcast_code, flags=re.IGNORECASE|re.VERBOSE)
            for val_len, val in defs:
                # As above, extract the dimensions
                search = re.search(
                    r"""
                    (?P<name>\w+)\s*    # Variable name (captured)
                    (\((?P<dims>.*)\)|) # Optional bracketed expression
                    """, val, flags=re.VERBOSE)
                if search:
                    # Add the name of the variable to the list
                    variable_names.append(search.group("name").lower())
                    if search.group("dims") is not None:
                        # But here also include the LEN as if it was
                        # another dimension
                        dim_expr = search.group("dims").split(",") + [val_len]
                        for iexp, expr in enumerate(dim_expr):
                            dim_expr[iexp] = "("+dim_expr[iexp]+")"

                        variable_definitions[vartype].append(
                            (val, "*".join(dim_expr)))
                    else:
                        # If it isn't an array just use the LEN
                        variable_definitions[vartype].append((val, val_len))

    return variable_names, variable_definitions

#-------------------------------------------------------------------------------
def get_setup_nml_type_calls(code):
    """
    Given the (pre-processed) code making up a UM routine, return all
    occurrences of calls to "setup_nml_type", which indicates the
    presence of a namelist broadcast call
    """
    # Look for subroutines containing calls to "setup_nml_type", retrieving
    # the information which will be needed to parse the contents
    setup_nml_type_calls = re.findall(
        r"""
        (?:^|\s+)                      # Start of line or at least 1 space
                                       # (non-capturing)
        SUBROUTINE\s+                  # Subroutine definition followed by at
                                       # least 1 space
        (?P<subname>\w+)               # Full name of routine (captured)
        (?:.(?!SUBROUTINE))*           # Any code following this (but not going
                                       # past another subroutine statement)
        CALL\s*setup_nml_type\s*\(\s*  # The start of the call
        (?P<n_types>\w+)\s*,\s*        # Argument 1 (number of types, captured)
        (?P<id>\w+)\s*                 # Argument 2 (nml type id, captured)
        ((?:,\s*\w+\s*=\s*\w+\s*)+)    # A bunch of keyword arguments
                                       # (captured as a single string)
        \s*\)
        """,
        code, flags=re.IGNORECASE|re.VERBOSE|re.MULTILINE|re.DOTALL)

    for call in setup_nml_type_calls:
        yield call

#-------------------------------------------------------------------------------
def check_code(file_to_check, failed_files):
    """
    A function to read a file and control the code checking
    """
    failures = 0

    # Read in the code
    code = read_and_preprocess_file(file_to_check)

    # Find the calls to "setup_nml_type"
    for subname, n_types_name, bcast_id_name, type_keywords in (
                                            get_setup_nml_type_calls(code)):

        if _DEBUG:
            banner_print("Processing \"{0:s}\"".format(subname))
            print "(File: {0:s})".format(file_to_check)
            print "\n(No of types variable: {0:s})".format(n_types_name)
            print "(Broadcast id variable: {0:s})".format(bcast_id_name)

        # Extract the type count variable names from the keywords
        # (these describe how many of each type there are)
        count_names = {}
        for varname, _ in _VARTYPES:
            search = re.search(r"n_{0:s}_in\s*=\s*(?P<var>\w+)".format(varname),
                               type_keywords, flags=re.IGNORECASE)
            if search:
                count_names[varname] = search.group("var")

        if _DEBUG:
            for var, name in count_names.items():
                print "(Count variable ({0:s}): {1:s})".format(var, name)

        # Extract the code contents of the subroutine containing this call
        # (to restrict the searching which follows in case of multiple
        # broadcasting calls in the same file)
        search = re.search(
            r"""
            (?:^|\s+)        # Start of line or at least 1 space (no capture)
            SUBROUTINE\s+    # Subroutine definition and at least 1 space
            {0:s}\W          # Will be substituted with subroutine name above
            (?P<code>.*)     # Entire contents of subroutine code (captured)
            (?:^|\s+)        # Start of line or at least 1 space (no capture)
            END\s*SUBROUTINE # End subroutine statement
            \s+{0:s}\W       # Followed by a space and the same name again
            """.format(subname), code,
            flags=re.MULTILINE|re.DOTALL|re.IGNORECASE|re.VERBOSE)
        
        subcode = search.group("code")

        # Extract the information about the number of each type
        success, output = get_type_counts(subcode, n_types_name, count_names)
        if success:
            variable_counts = output
        else:
            failures += 1
            print "\nERROR: {0:s} ({0:1})".format(file_to_check, subname)
            print output
            continue

        # Next, we need to find the call to mpl_bcast which references the
        # id found above
        search = re.search(
            r"""
            (?:^|\s+)                     # Start of line or at least 1 space
            CALL\s*mpl_bcast\s*\(\s*      # The start of the call
            (?P<btype>\w+)\s*,\s*         # Argument 1 (name of type, captured)
            \w+,\s*                       # Argument 2 (throw away)
            {0:s},\s*                     # Argument 3 (nml type id, must
                                          #             match our id)
            .*\)                          # Anything else (throw away)
            """.format(bcast_id_name), subcode, flags=re.VERBOSE|re.IGNORECASE)

        # This should exist
        if search:
            mpl_bcast_type_name = search.group("btype")
            if _DEBUG:
                print "(Broadcast type name: {0:s})".format(mpl_bcast_type_name)
        else:
            print "\nERROR: {0:s} ({1:s})".format(file_to_check, subname)
            print ("ERROR: Unable to find broadcast type name "
                   "in \"{0:s}\"".format(subname))
            failures += 1
            continue

        # From this, we can find out the name of the instance holding the type
        search = re.search(
            r"""
            (?:^|\s+)                    # Start of line or at least 1 space
            TYPE\s*                      # Type definition
            \(\s*(?P<name>\w+)\s*\)      # Name of type (captured)
            \s*::\s*                     # Declaration double-colon
            {0:s}                        # Name (substituted from above)
            """.format(mpl_bcast_type_name),
            subcode, flags=re.VERBOSE|re.IGNORECASE)

        # This should exist
        if search:
            mpl_bcast_instance_name = search.group("name")
            if _DEBUG:
                print ("(Broadcast type variable: {0:s})"
                       .format(mpl_bcast_instance_name))
        else:
            print "\nERROR: {0:s} ({1:s})".format(file_to_check, subname)
            print ("ERROR: Unable to find broadcast type instance "
                   "in \"{0:s}\"".format(subname))
            failures += 1
            continue

        # Now find the code which defines the broadcast type
        broadcast_type = re.search(
            r"""
            (?:^|\s+)         # Start of line or at least one space
                              # (non-capturing)
            TYPE\s+           # Opening of type definition plus at least
                              # one space
            {0:s}             # The type name (from above)
            (?P<code>.*)      # Entirety of the code inside the type
                              # definition (captured as "code") until..
            (?:^|\s+)         # Another start of line + space
            END\s*TYPE        # The end type statement
            \s+{0:s}          # Followed by at least one space and the
                              # same name as above
            """.format(mpl_bcast_instance_name),
            subcode, flags=re.MULTILINE|re.DOTALL|re.IGNORECASE|re.VERBOSE)

        # This should exist, bail out otherwise
        if broadcast_type:
            bcastcode = broadcast_type.group("code")
        else:
            print "\nERROR: {0:s} ({1:s})".format(file_to_check, subname)
            print ("ERROR: Unable to find broadcast namelist "
                   "type definition in \"{0:s}\"".format(subname))
            failures += 1
            continue

        # Get the definitions of each variable from the broadcast type, as well
        # as the name of each variable contained in the definition
        variable_names, variable_definitions = get_definitions(bcastcode)

        # Check for the pairs of transfer statements for the variables used in
        # each broadcast call
        in_search = re.findall(
            r"""
            {0:s}\s*%\s*(?P<name>\w+) # Type attribute name
            (?:\s*\(.*?\)|)           # Possible dimension information
            \s*=\s*                   # Equals
            (?P<local>\w+             # Local variable name (captured)
            (?:\s*%\s*\w+\s*)*)\s*    # possibly including attributes
            """.format(mpl_bcast_type_name),
            subcode.lower(), flags=re.IGNORECASE|re.VERBOSE)

        out_search = re.findall(
            r"""
            (?P<local>\w+             # Local variable name (captured)
            (?:\s*%\s*\w+\s*)*)\s*    # possible including attributes
            (?:\s*\(.*?\)|)           # Possible dimension information
            \s*=\s*                   # Equals                
            {0:s}\s*%\s*(?P<name>\w+) # Type attribute name
            """.format(mpl_bcast_type_name),
            subcode.lower(), flags=re.IGNORECASE|re.VERBOSE)

        # Form dictionaries from the results above, but flip around the
        # second one so that the keys are always the variable name taken
        # from the broadcast namelist type
        in_search  = dict(in_search)
        out_search = dict([reversed(elt) for elt in out_search])

        # These can now be checked
        for name in variable_names:
            if name not in in_search.keys():
                print "\nERROR: {0:s} ({1:s})".format(file_to_check, subname)
                print ("ERROR: Unable to find transfer INTO broadcast type: "
                       "{0:s} % {1:s}".format(mpl_bcast_type_name, name))
                failures += 1
                continue

            if name not in out_search.keys():
                print "\nERROR: {0:s} ({1:s})".format(file_to_check, subname)
                print ("ERROR: Unable to find transfer FROM broadcast type: "
                       "{0:s} % {1:s}"
                       .format(mpl_bcast_type_name, name))
                failures += 1
                continue

            local_in  = in_search[name].replace(" ","").strip()
            local_out = out_search[name].replace(" ","").strip()
            
            if  local_in != local_out :
                print "\nERROR: {0:s} ({1:s})".format(file_to_check, subname)
                print ("ERROR: Local variable names in broadcast transfer "
                       "do not agree: {0:s} != {1:s}"
                       .format(local_in, local_out))                
                failures += 1

        # Now we are ready to compare the defined dimensions to the
        # indicated counts and flag up any errors
        for vartype, varname in _VARTYPES:

            # Construct an equivalent expression for the total size of the
            # defined variables by adding their individual expressions
            definitions = "+".join(
                [dim for _, dim in variable_definitions[vartype]])

            # Create a prime-dictionary which maps a prime against each
            # variable name in the variable count statement or definitions
            pdict, pdict_desc = create_prime_dict(
                variable_counts[vartype] + " " + definitions)

            # Use this to evaluate the variable count statement
            count, count_desc = expr_eval(variable_counts[vartype], pdict)

            if definitions == "":
                defs = 0
            else:
                defs, def_desc = expr_eval(definitions, pdict)

            if _DEBUG and (count != 0 or defs != 0):
                print "\nType: {0:s}".format(varname)
                print pdict_desc
                print "Variable counts: "
                print "\n".join(wrap(count_desc, width=_OUTPUT_LINE_WIDTH,
                                     subsequent_indent="   "))
                print "Variable definitions: "
                print "\n".join(wrap(def_desc, width=_OUTPUT_LINE_WIDTH,
                                     subsequent_indent="   "))

            # These two measures should be equivalent if the routine has
            # been setup correctly, but now try to determine why they
            # aren't
            if defs != count:
                print "\nERROR: {0:s} ({1:s})".format(file_to_check, subname)
                print ("ERROR: {0:s} count and definition do not agree"
                       .format(varname))
                    
                failures += 1
                diff = defs - count

                abs_diff = abs(diff)
                diff_accounted_for = 0

                # We can try to work out if it was a named dimension by
                # substituting combinations of the variables with zeros,
                # and examine what effect this has on the difference

                # But first, we need to account for any differences caused
                # by scalar variables, which we can do by setting all of
                # the variables to zero
                test_dict = pdict.copy()
                for key in test_dict.keys():
                    test_dict[key] = "0"
                    
                # Re-calculate the two results using the above dict
                count, count_desc = expr_eval(variable_counts[vartype],
                                              test_dict)
                defs, def_desc = expr_eval(definitions, test_dict)

                if _DEBUG:
                    print "\nChecking for scalar diffs (set variables to zero):"
                    print "Scalar Count:"
                    print "\n".join(wrap(count_desc, width=_OUTPUT_LINE_WIDTH,
                                         subsequent_indent="   "))
                    print "Scalar Defs:"
                    print "\n".join(wrap(def_desc, width=_OUTPUT_LINE_WIDTH,
                                         subsequent_indent="   "))

                # If these are different it means some scalars are missing
                scalar_diff = count - defs
                if scalar_diff < 0:
                    print ("ERROR: Scalar variable/s missing from "
                           "'n_{0:s}' count totalling size: {1:d}"
                           .format(vartype, -scalar_diff))
                elif scalar_diff > 0:
                    print ("ERROR: {0:s} scalar variable/s missing "
                           "from definitions totalling size: {1:d}"
                           .format(varname, scalar_diff))
                # And add this difference to the total
                diff_accounted_for += abs(scalar_diff)

                # Now build up all possible combinations of the variables
                combos = []
                for i in range(len(pdict)):
                    combos.extend([j for j in
                                   itertools.combinations(
                                       pdict.keys(), i+1)])

                for combo in combos:
                    # A running total of the product of the variables
                    # not set to zero in this combination
                    combo_total = 1

                    # Create a test dictionary based on the original
                    test_dict = pdict.copy()
                    # This will be the name of the combo in any printed
                    # output message
                    var = "*".join(combo)
                    for key in test_dict.keys():
                        if key in combo:
                            # For keys being included build up the
                            # product by multiplying them by it
                            combo_total *= int(test_dict[key])
                        else:
                            # For those being omitted set their
                            # value to zero
                            test_dict[key] = "0"


                    # Re-calculate the two results using the above dict
                    count, count_desc = expr_eval(variable_counts[vartype],
                                                  test_dict)
                    defs, def_desc = expr_eval(definitions, test_dict)
                    if _DEBUG:
                        print "\nChecking for differences with:"
                        for key, val in test_dict.items():
                            print "   {0:s}: {1:s}".format(key, val)
                        print "Count:"
                        print "\n".join(wrap(count_desc,
                                             width=_OUTPUT_LINE_WIDTH,
                                             subsequent_indent="   "))
                        print "Defs:"
                        print "\n".join(wrap(def_desc,
                                             width=_OUTPUT_LINE_WIDTH,
                                             subsequent_indent="   "))      

                    # Remove the effect of any scalar difference found above
                    subtr_total = count - defs - scalar_diff

                    if (subtr_total != 0 and
                        abs(subtr_total) % combo_total >= 0):
                        num = abs(subtr_total) / combo_total
                        if num < 1:
                            # If this matches it's because a difference was 
                            # found due to the effect of only one component of
                            # a non-linear combination being tested
                            continue
                        elif num == 1:
                            # Output is just the variable name
                            out = var
                        elif num > 1:
                            # Or an indication of how many multiples
                            out = "{0:d}*{1:s}".format(num, var)

                        if subtr_total < 0:
                            print ("ERROR: Variable/s missing from "
                                   "'n_{0:s}' count with size: {1:s}"
                                   .format(vartype, out))
                        else:
                            print ("ERROR: {0:s} variable/s missing from "
                                   "definitions with size: {1:s}"
                                   .format(varname, out))

                        # Subtract the difference from the total
                        diff_accounted_for += num*combo_total

                # After the above if the difference has anything left
                # over it is due to a mistake in how this algorithm works
                if abs_diff - diff_accounted_for > 0:
                    print ("ERROR: Un-handled problem/s with "
                           "{0:s} count/definitions".format(varname))

    if failures > 0:
        failed_files.append(file_to_check)
                            
    return failed_files

#-------------------------------------------------------------------------------
def files_to_process(filepath, ignore_list):
    """
    Generate list of files in given filepath.  Ignore any files matching
    the patterns in the ignore list
    """
    files = []  
    for root, _, filenames in os.walk(filepath):
        for filename in filenames:
            if re.match(_FILE_PATTERN, filename):
                path_to_file = os.path.join(root, filename)
                if any([ignore in path_to_file for ignore in ignore_list]):
                    print "WARNING: Ignoring file: {0:s}".format(path_to_file)
                    continue
                files.append(path_to_file) 

    return files  

#-------------------------------------------------------------------------------
def main(inputs, ignore_list):
    """ main program block """

    banner_print("Running namelist declaration checker")
    files_to_check = []
    for file_input in inputs:
        if os.path.isfile(file_input):
            print "Source (file): {0:s}".format(file_input)
            files_to_check.append(file_input)
        elif os.path.isdir(file_input):
            print "Source (dir) : {0:s}".format(file_input)
            files_to_check.extend(files_to_process(file_input, ignore_list))
        else:
            print "ERROR: Input sources must be files/directories"
            print "ERROR: {0:s} is neither".format(file_input)
            sys.exit(1)

    print "\nFound {0:d} files to check".format(len(files_to_check))

    failed_files = []
    for item in files_to_check:
        failed_files = check_code(item, failed_files)

    failures = len(failed_files)
    banner_print("Checks completed with {0:d} failure{1:s}"
                 .format(failures, "s" if failures != 1 else ""))
    print ""

    if len(failed_files) > 0:
        print "ERROR: Failed files:"
        for failure in failed_files:
            print " * {0:s}".format(failure)
        print ""
        print "\n".join(wrap("To obtain a more detailed breakdown of what "
                             "went wrong, you could re-run with the debug "
                             "output turned on, supplying only the name/s "
                             "of the affected file/s above",
                             width=_OUTPUT_LINE_WIDTH))
        print ""
        sys.exit(1)

#-------------------------------------------------------------------------------
def parse_options():
    """Parse command line code options."""
    usage = "usage: %prog [options] directory/file1 [[directory/file2] ...]"
    description = ("This script will scan UM Fortran source for "
                   "namelist-broadcasting statements and check "
                   "for common errors in their definition.  Arguments may "
                   "be any combination of files and/or directories; "
                   "individual files will be scanned, and all files within "
                   "directories and their sub-directories will be scanned.")
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-d", "--debug",
                      action="store_true", dest="debug", default=False,
                      help="debugging output")
    parser.add_option("--ignore", action="store", dest="ignore", default=None,
                      help=("ignore filename/s containing "
                            "(comma separated list of patterns)"))
    (opts, args) = parser.parse_args()
     
    if len(args) >= 1:
        opts.inputs = args[0:]
    else:
        parser.error("Please supply directory/files to work on")

    return opts

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Parse command line options
    OPTS = parse_options()

    # Toggle global debugging on/off
    if OPTS.debug:
        _DEBUG = True

    # Change the ignore input into a list
    if OPTS.ignore is None:
        OPTS.ignore = []
    else:
        OPTS.ignore = OPTS.ignore.split(",")
    
    main(OPTS.inputs, OPTS.ignore)
