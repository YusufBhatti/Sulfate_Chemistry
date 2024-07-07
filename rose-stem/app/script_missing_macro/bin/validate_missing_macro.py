#!/usr/bin/env python
# *********************************COPYRIGHT************************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *********************************COPYRIGHT************************************

import os
import re
import sys

AFTER_TAGS = {}
INITIAL_REVISIONS = {
    'um-createbc': 'vn10.0_t2',
    'um-fcm-make': 'vn8.6',
    'um-atmos': 'vn8.6',
    }


def read_file(fname):
    '''Return contents of a file given filename.'''
    with open(fname, 'r') as fh:
        lines = fh.readlines()
    return lines


def index_after_tags(meta_path):
    '''Find the correct AFTER_TAG for each type of app'''
    meta_types = os.listdir(meta_path)

    # Loop over metadata types
    for meta_type in meta_types:
        # Ignore metadata without valid upgrade macros
        if meta_type not in INITIAL_REVISIONS:
            break

        # Each upgrade macro tag should appear twice, once as BEFORE and once
        # as AFTER, except for the initial and final revisions
        tag_index = {}
        macro_files = []
        location = os.path.join(meta_path, meta_type)
        
        # Find upgrade macro files
        for item in os.listdir(location):
            result = re.search(r'\.py$', item)
            if result:
                macro_files.append(item)

        # Search each macro file for tags
        for macro_file in macro_files:
            lines = read_file(os.path.join(location, macro_file))

            for line in lines:
                before_result = re.search(r'^\s*BEFORE_TAG\s*=\s*"(\S+)"',
                                          line)
                if before_result:
                    before_tag = before_result.group(1)
                    if before_tag in tag_index:
                        del tag_index[before_tag]
                    else:
                        tag_index[before_tag] = 1

                after_result = re.search(r'^\s*AFTER_TAG\s*=\s*"(\S+)"',
                                         line)
                if after_result:
                    after_tag = after_result.group(1)
                    if after_tag in tag_index:
                        del tag_index[after_tag]
                    else:
                        tag_index[after_tag] = 1

        # Remove the initial revision for this metadata type
        del tag_index[INITIAL_REVISIONS[meta_type]]

        if len(tag_index) != 1:
            sys.exit("Unable to process tag {0} - multiple AFTER_TAGS?".format(
                     meta_type))

        # The remaining value is the last AFTER_TAG
        for remaining_value in tag_index:
            # Process _tXXXX values
            remaining_value = re.sub(r'_tXXXX', r'', remaining_value)
            AFTER_TAGS[meta_type] = remaining_value

    return


def test_app(app_path, app):
    '''Test an app has the correct metadata revision'''

    # Read the rose-app.conf file
    app_lines = read_file(os.path.join(app_path, app, 'rose-app.conf'))

    success = True
    valid_test = False

    # Find the meta= line (if any)
    for line in app_lines:
        result = re.search(r'meta=(\S+)/(\S+)', line)
        if result:
            meta_type = result.group(1)
            tag = result.group(2)

            # Check if this is a app type with upgrade macros
            if meta_type in AFTER_TAGS:
                valid_test = True

                # Check if this compares with the correct AFTER_TAG
                if tag != AFTER_TAGS[meta_type]:
                    print ("[FAIL] {0} has meta-data revision {1}".format(app,
                           tag) + " but the most recent AFTER_TAG is " +
                           "{0}".format(AFTER_TAGS[meta_type]))
                    success = False

    if success and valid_test:
        print "[INFO] {0} has the correct metadata revision".format(app) 
    return success


def main():
    '''Main program.'''
    result = 0

    # Locate the meta-data and app directories
    work_dir = os.environ['CYLC_TASK_WORK_DIR']
    meta_path = os.path.join(work_dir, '..',
                             'script_source', 'um', 'rose-meta')
    app_path = os.path.join(work_dir, '..',
                            'script_source', 'um', 'rose-stem', 'app')

    # Find the correct after tags for each app type
    index_after_tags(meta_path)

    # Loop over apps to check they have the correct metadata revision
    for app in os.listdir(app_path):
        success = test_app(app_path, app)
        if not success:
            result += 1

    if result == 0:
        print "[INFO] All apps have the correct AFTER_TAG value"
    else:
        print "[FAIL] {0} apps do not have the correct AFTER_TAG value".format(
            result)

    sys.exit(result)


if __name__ == '__main__':
    main()
