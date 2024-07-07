#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved. 
# For further details please refer to the file COPYRIGHT.txt 
# which you should have received as part of this distribution. 
# *****************************COPYRIGHT******************************* 

# Owner: UM System Development Team
# Purpose: Display graphs created by monitoring.py

import cgitb
cgitb.enable()
import glob
import re
import sys

DIRECTORY='/home/h01/frzz/public_html/monitoring'

graphs = glob.glob("%s/*.png"%(DIRECTORY))


print 'Content-Type: text/html\n\n'
print '<h1>Unified Model Continuous Trunk Monitoring</h1>'

for plot in sorted(graphs):
    if 'scm' in plot:
        continue
    name = re.sub(r'.*atmos', r'atmos', plot)
    name = re.sub(r'\.png', r'', name)
    if name.endswith('_progr'):
        continue
    print "<h2>%s</h2>"%(name)
    webplot = re.sub(r'/home/h01/frzz/public_html', r'www-nwp/~frzz', plot)
    print '<img src="http://%s">'%(webplot)
