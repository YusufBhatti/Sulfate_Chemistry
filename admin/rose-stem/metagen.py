#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved. 
# For further details please refer to the file COPYRIGHT.txt 
# which you should have received as part of this distribution. 
# *****************************COPYRIGHT******************************* 
"""Generate rose-stem metadata for RUN_NAMES.

Author/Owner: UM System Development Team
"""

import glob
import os
import re
import sys

def read_file(fname):
    '''Return contents of a file given filename.'''
    with open(fname, 'r') as fh:
        lines = fh.readlines()
    return lines


def write_file(fname, lines, newline=False):
    '''Write a file given names and contents. The optional newline argument
    adds a newline at the end of each element of the list.'''
    with open(fname, 'w') as fh:
        for line in lines:
            if newline:
                fh.write("%s\n"%(line))
            else:
                fh.write(line)


def run_command(command, shell=False):
    '''Given a command as a string, run it and return the exit code, standard
    out and standard error. The optional shell argument allows a shell to
    be spawned to allow multiple commands to be run.'''

    import subprocess

    if shell:
        # Create the Popen object and connect out and err to pipes using
        # the shell=True option.
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE, shell=True)
    else:
        # Turn command into a list
        command_list = command.split()

        # Create the Popen object and connect out and err to pipes
        p = subprocess.Popen(command_list, stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                                       
    # Do the communicate and wait to get the results of the command
    stdout, stderr = p.communicate()
    rc = p.wait()

    # Reformat stdout
    stdout = ''.join(stdout)
    stdout = stdout.split("\n")
    
    return rc, stdout, stderr


def parse_meta_and_delete_old_run_names(metafile):
    '''Read in the current rose-stem metadata and delete the metadata for 
    RUN_NAMES so we can replace it later.'''
    contents = read_file(metafile)
    newlines = []
    in_section = False
    for line in contents:
        if re.search(r'RUN_NAMES', line):
            in_section = True
        if in_section and re.search(r"^\[.*\]\s*$" , line):
            in_section = False
        if in_section:
            continue
        newlines.append(line)
    return newlines


def new_run_names(choices, jobs):
    '''Generate the new RUN_NAMES metadata.'''
    choices_string = ','.join(choices)
    content = [
      '',
      '[jinja2:suite.rc=RUN_NAMES]',
      'compulsory=true',
      'description=Task groups to run',
      'ns=Tasks',
      'widget[rose-config-edit]=rose.config_editor.valuewidget.choice.ChoicesValueWidget',
      '                        =--all-group=all --editable --format=python --guess-groups',
      '                  =--choices=%s'%(choices_string),
              ]

    for job in sorted(jobs):
        content.append('                  = %s'%(job))

    return content

def find_choices(stemdir):
    '''Examine the graph-group file for each site (if it exists) and generate
    the list of valid choices.'''
    sitedir = os.path.join(stemdir, 'site')
    sites = glob.glob("%s/*"%(sitedir))
    choices = []
    for site in sites:
      groupfile = os.path.join(sitedir, site, 'graph-group.rc')
      if not os.path.isfile(groupfile):
          continue
      lines = read_file(groupfile)
      for line in lines:
          result = re.search(r'^\s*"(.*)"\s*:', line)
          if result:
              choices.append(result.group(1))
    return choices

def find_jobs(stemdir):
    '''Examine all runtime files for all sites and generate a list of valid
    jobs.'''
    sitedir = os.path.join(stemdir, 'site')
    sites = glob.glob("%s/*"%(sitedir))
    job_list = set([])
    for site in sites:
      runtime_files = glob.glob(os.path.join(sitedir, site, "runtime-*"))
      runtime_files.append('runtime-common.rc')
      for fname in runtime_files:
          if not os.path.isfile(fname):
              continue
          lines = read_file(fname)
          for line in lines:
              result = re.search(r'"(.*?)"\s*in\s*name_graphs_out', line)
              if result:
                  name = result.group(1)
                  if re.search(r'~MACHINE~', name):
                      spice_name = re.sub(r'"~MACHINE~"', 'spice', name)
                      desktop_name = re.sub(r'"~MACHINE~"', 'desktop', name)
                      job_list.add(spice_name)
                      job_list.add(desktop_name)
                  else:
                      job_list.add(name)
                  
    # Manually append "fcm_make" as that's specified differently
    job_list.add("fcm_make")
    
    # Manually ignore the AQUM compiler checking job as it's currently not
    # working
    job_list.remove("metohpc_aqum_nd_comp_check")

    return list(job_list)

def divide_jobs(jobs):
    linux = []
    hpc = []
    xc40 = []
    hpc_recon = []
    xc40_recon = []
    spice = []
    spice_recon = []
    desktop = []
    desktop_recon = []
    
    for job in jobs:
        if re.search(r'metohpc', job):
            hpc.append(job)
            if re.search(r'recon', job):
                hpc_recon.append(job)
        if re.search(r'metolinux', job):
            linux.append(job)
            if re.search(r'recon', job):
                linux_recon.append(job)
        if re.search(r'meto_xc40', job):
            xc40.append(job)
            if re.search(r'recon', job):
                xc40_recon.append(job)
        if re.search(r'spice', job):
            spice.append(job)
            if re.search(r'recon', job):
                spice_recon.append(job)
        if re.search(r'desktop', job):
            desktop.append(job)
            if re.search(r'recon', job):
                desktop_recon.append(job)

    return linux, hpc, xc40, spice, desktop, hpc_recon, xc40_recon, spice_recon, desktop_recon

def convert_to_graph_group(jobs, name):
    lines = [ '    "%s" : [\n'%name,]
    for job in jobs:
      lines.append('                   "%s",\n'%(job))
    lines.append('           ],\n\n')
    return lines
    
        
def remove_from_file(fname, groups):
    lines = read_file(fname)
    newlines = []
    deleting = False
    for line in lines:
        for group in groups:
            if re.search(r'\s*"%s"\s*:'%(group), line):
                deleting = True
        if not deleting:
            if re.search(r'}', line):
                continue
            newlines.append(line)
        if deleting and re.search(r']', line):
            deleting = False
            
    return newlines


def define_meto_groups(jobs):
    linux, hpc, xc40, spice, desktop, hpc_recon, xc40_recon, spice_recon, desktop_recon = divide_jobs(jobs)
    meto_group_file = 'site/meto/graph-group.rc'
    lines = remove_from_file(meto_group_file, ['hpc', 'linux', 'hpc_recon', 'linux_recon', 'xc40', 'xc40_recon', 'spice', 'desktop', 'spice_recon', 'desktop_recon'])
    lines = lines + convert_to_graph_group(hpc, 'hpc')
    lines = lines + convert_to_graph_group(xc40, 'xc40')
    lines = lines + convert_to_graph_group(linux, 'linux')
    lines = lines + convert_to_graph_group(spice, 'spice')
    lines = lines + convert_to_graph_group(desktop, 'desktop')
    lines = lines + convert_to_graph_group(hpc_recon, 'hpc_recon')
    lines = lines + convert_to_graph_group(xc40_recon, 'xc40_recon')
    lines = lines + convert_to_graph_group(spice_recon, 'spice_recon')
    lines = lines + convert_to_graph_group(desktop_recon, 'desktop_recon')
    lines = lines + [ "    }\n", "%}\n"]
    write_file(meto_group_file, lines, newline=False)


if __name__== '__main__':
    
    cwd = os.getcwd()
    if re.search(r'rose-stem$', cwd):
        stemdir = cwd
    else:
        sys.exit("Please run in rose-stem subdirectory.")    

    metafile = os.path.join(stemdir, 'meta', 'rose-meta.conf')
    newlines = parse_meta_and_delete_old_run_names(metafile)

    choices = find_choices(stemdir)
    jobs = find_jobs(stemdir)
    newlines = newlines + new_run_names(choices, jobs)
    write_file(metafile, newlines, newline=True)

    # Tidy file
    rc, stdout, stderr = run_command('rose config-dump --file=%s'%(metafile))

    # Define the meto automatically-generated groups
    define_meto_groups(jobs)
