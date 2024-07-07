# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------                                                                    #
# (C) Crown copyright 1990-2018 Met Office. All rights reserved.               #
#                                                                              #
# Use, duplication or disclosure of this code is subject to the restrictions   #
# as set forth in the licence. If no licence has been raised with this copy    #
# of the code, the use, duplication or disclosure of it is strictly            #
# prohibited. Permission to do so must first be obtained in writing from the   #
# IPR Manager at the following address:                                        #
#                                                                              #
# Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB, United Kingdom             #
#                                                                              #
#------------------------------------------------------------------------------#

import os

import rose.config_editor.plugin.um.widget.stash
import stash_parse

META_DIR = os.path.dirname(  # lib
              os.path.dirname(  # python
                 os.path.dirname(  # widget
                    os.path.dirname(os.path.abspath(__file__)))))


class StashSummaryDataPanelv1(
      rose.config_editor.plugin.um.widget.stash.BaseStashSummaryDataPanelv1):

    """This is a class for displaying and editing STASH requests.

    It should be referenced via a:
    widget[rose-config-edit:sub-ns]=stash.StashSummaryDataPanelv1
    entry in the metadata. The default STASHmaster file used in the
    widget is the one in the metadata directory - supplying a path
    as an argument to the widget metadata will override it, e.g.:
    widget[rose-config-edit:sub-ns]=stash.StashSummaryDataPanelv1 /some/path/to/STASHmaster/files/

    For more information, see the base class within the Rose code.

    """

    STASHMASTER_PATH = os.path.join(META_DIR, "etc", "stash",
                                    "STASHmaster")
    STASHMASTER_META_PATH = STASHMASTER_PATH
    STASH_PACKAGE_PATH = os.path.join(META_DIR, "etc", "stash",
                                      "package")

    def get_stashmaster_lookup_dict(self):
        """Provide the necessary rose.config.ConfigNode object."""
        parser = stash_parse.StashMasterParserv1(
                       self.stashmaster_directory_path)
        return parser.get_lookup_dict()
