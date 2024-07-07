#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# *****************************COPYRIGHT*******************************
"""Parse the trunk log into TrunkCommit and TrunkTicket objects.

Author/Owner: Stuart Whitehouse
"""

from collections import defaultdict
from datetime import datetime, timedelta
from optparse import OptionParser
import re
import subprocess
import sys

FCM = 'fcm'


# Log-interpretation mappings

# 'commits' contains a list of properties set in the log message
# which are colon-separated. Note that ticket number is an implicit
# first item, so the first item in the list is the second
# colon-separated property in the log message

# 'tickets' contains the mapping of properties in tickets to those
# in commits; the key is the name of the ticket property (e.g.
# 'initial_log'). The value is a list of the corresponding property
# in the commit (e.g. 'description', together with the number in the
# list of commits from which to take that item. Useful values of the
# latter are 0 (first commit in the list) and -1 (final commit in the
# list).
UM = {
    'commits': ['author', 'description', 'ticket_type', 'code_area',
                'regression', 'severity'],
    'sep': ':',
    'sep_sub': None,
    'tickets': {
        'author': ['author', 0],
        'initial_log': ['description', 0],
        'code_area': ['code_area', 0],
        'ticket_type': ['ticket_type', -1],
        'regression': ['regression', -1],
        'severity': ['severity', -1]
    }
}

VAR = {
    'commits': ['author', 'site', 'description', 'ticket_type', 'code_area',
                'regression', 'severity'],
    'sep': ' : ',
    'sep_sub': ' ; ',
    'tickets': {
        'author': ['author', 0],
        'site' : [ 'site', 0],
        'initial_log': ['description', 0],
        'code_area': ['code_area', 0],
        'ticket_type': ['ticket_type', -1],
        'regression': ['regression', -1],
        'severity': ['severity', -1]
    }
}

class TrunkLogParser(object):
    """Class to parse the trunk log and return commit objects. Requires
       a list as a mandatory argument, with the first value as the first
       revision of trunk to search, and the second as the last revision
       to search. Optionally takes a method (defaults to UM), the
       trunk name, and whether to abort if it fails to parse a
       changeset."""
    def __init__(self, revs, method=None, trunk='fcm:um.x_tr',
                 abort_on_fail=False):
        if method is None:
            method = UM
        self.method = method
        self.first_rev = revs[0]
        self.last_rev = revs[1]
        self.trunk = trunk
        self.abort_on_fail = abort_on_fail
        self.log_message = self._get_trunk_log()
        self.commits = self._parse_log(method)
        if len(self.commits) == 0 and self.abort_on_fail:
            raise NoTrunkRevisionsError(self.first_rev, self.last_rev,
                                        self.trunk)
        self.tickets = self._convert_to_tickets(method)

    def _get_trunk_log(self):
        """Runs "fcm log" and returns stdout."""
        command = "%s log -r %s:%s %s" % (FCM, self.first_rev,
                                          self.last_rev, self.trunk)
        retc, stdout, stderr = _run_command(command)

        if retc is not 0 and self.abort_on_fail:
            raise TrunkParseFailed(command, retc, stdout, stderr)

        return stdout

    def _parse_log(self, method):
        """Parse the log file and create objects."""
        commits = []
        current_message = []
        for line in self.log_message:
            if re.search(r'^-+$', line):
                if not current_message:
                    continue
                commit = TrunkCommit(current_message,
                                     self.abort_on_fail, method)
                if commit.revision is not 0:
                    commits.append(commit)
                current_message = []
            else:
                current_message.append(line)
        return commits

    def _convert_to_tickets(self, method):
        """Convert commits to tickets."""
        tickets = []
        commits_by_ticket = self.group('ticket', 'commits')
        for i in commits_by_ticket:
            tickets.append(TrunkTicket(commits_by_ticket[i], method))
        return tickets

    def group(self, attribute, item, search=None, invert=False):
        """Return a dictionary of items grouped by a given attribute.
        Optionally provide a dictionary containing additional properties
        and values which must match, or, if invert=True, must not match,
        for them to be part of the returned dictionary. Note that if
        you have multiple things to search in the dictionary the
        invert option may not work as expected."""
        if not hasattr(self, item):
            raise AttributeNotFoundError(item, self)

        items = getattr(self, item)
        dictionary = defaultdict(list)
        for i in items:
            if not hasattr(i, attribute):
                raise AttributeNotFoundError(attribute, i)
            attrib = getattr(i, attribute)
            if isinstance(attrib, list):
                dictionary[str(len(attrib))].append(i)
            else:
                dictionary[attrib].append(i)

        if search:
            matching = {}
            for i, listitems in dictionary.iteritems():
                for j in listitems:
                    add_to_group = True
                    for prop, value in search.iteritems():
                        if invert:
                            if re.search(value, getattr(j, prop)):
                                add_to_group = False
                        else:
                            if not re.search(value, getattr(j, prop)):
                                add_to_group = False
                    if add_to_group:
                        if i in matching:
                            matching[i].append(j)
                        else:
                            matching[i] = [j]
            dictionary = matching
        return dictionary

    def most_frequent(self, attribute, item, search=None, invert=False):
        """Return a list of the items with the largest frequency with the
        names of the matching categories."""
        groups = self.group(attribute, item, search, invert)
        freq = None
        for i in groups:
            if len(groups[i]) > freq:
                freq = len(groups[i])
        things = []
        props = []
        for j in groups:
            if len(groups[j]) == freq:
                props.append(j)
                [things.append(k) for k in groups[j]]
        return props, things

    def least_frequent(self, attribute, item, search=None, invert=False):
        """Return a list of the items with the smaller frequency with the
        names of the matching categories. Note that the frequency must
        exist in at least one commit."""
        groups = self.group(attribute, item, search, invert)
        freq = 9999999999999
        for i in groups:
            if len(groups[i]) < freq:
                freq = len(groups[i])
        things = []
        props = []
        for j in groups:
            if len(groups[j]) == freq:
                props.append(j)
                [things.append(k) for k in groups[j]]
        return props, things

    def largest(self, attribute, item, search=None, invert=False):
        """Return a list of the items with the largest value."""
        groups = self.group(attribute, item, search, invert)
        freq = None
        for i in groups:
            if int(i) > freq:
                freq = int(i)
        things = []
        props = []
        for j in groups:
            if int(j) == freq:
                props.append(j)
                [things.append(k) for k in groups[j]]
        return props, things

    def ticket_by_number(self, number):
        """Return a ticket given ticket number."""
        value = None
        for i in self.tickets:
            if i.ticket_number == number:
                value = i
                break
        return value

    def commit_by_number(self, number):
        """Return a commit given commit number."""
        value = None
        for i in self.commits:
            if i.revision == number:
                value = i
                break
        return value

    def trac_format_log(self):
        """Return the entire log in trac wiki formatting."""
        commits = self.commits
        firstline = " || '''Revision''' || '''Ticket ''' || "
        for i in self.method['commits']:
            firstline += "'''"+i.title()+"'''"
            firstline += " || "
        print firstline
        for commit in commits:
            line = " || [%s] || #%s || "%(commit.revision, commit.ticket)
            for i in self.method['commits']:
              line += getattr(commit, i)
              line += " || "
            print line            

class TrunkTicket(object):
    """Class representing a ticket committed to the trunk."""
    def __init__(self, commits, method):
        # Properties defined by this module
        self.commits = commits
        self.num_commits = len(self.commits)

        # Properties automatically set by FCM
        self.initial_datetime = self.commits[0].datetime
        self.final_datetime = self.commits[-1].datetime
        self.reviewer = self.commits[0].reviewer

        # Properties set by universal convention
        self.ticket_number = self.commits[0].ticket

        # System-defined settings
        for ticket_prop, commit_prop in method['tickets'].iteritems():
            commit_property_to_set = getattr(self.commits[commit_prop[1]],
                                             commit_prop[0])
            setattr(self, ticket_prop, commit_property_to_set)

    def __str__(self):
        return "#%s" % (self.ticket_number)

    __repr__ = __str__


class TrunkCommit(object):
    """Class which describes a trunk commit."""
    def __init__(self, log_message, abort_on_fail, method):
        self.original_log = log_message
        self.revision = 0
        message = str()
        for line in self.original_log:
            if re.search(r'^Reversed', line):
                continue
            if re.search(r'^Custom merge', line):
                continue
            elif re.search(r'^r\d+', line):
                elements = line.split(" | ")
                self.revision = int(re.sub(r'^r', r'', elements[0]))
                self.reviewer = elements[1].strip()
                self.text_datetime = elements[2].strip()
                elements[2] = re.sub(r'(\d+-\d+-\d+\s*\d+:\d+:\d+).*', r'\1',
                                     elements[2].strip())
                self.datetime = datetime.strptime(elements[2],
                                                  "%Y-%m-%d %H:%M:%S")
            elif re.search(r'^Merged into', line) or re.search(r'^\s*$', line):
                continue
            else:
                message = message + line.rstrip()
        if message:
            elements = message.split(method['sep'])
            if len(elements) < len(method['commits']) + 1:
                if abort_on_fail:
                    raise LogParseError(self.revision, message)
                else:
                    print "[WARN] Failed to parse revision r%s" % (
                          self.revision) + ' "%s" - ignoring changeset' % (
                              message)
                    self.revision = 0
                    return
            self.ticket = re.sub(r'^#', r'', elements[0]).strip()
            i = 0
            for prop in method['commits']:
                i += 1
                if method['sep_sub'] is not None:
                    if method['sep_sub'] in elements[i]:
                        element_subs = elements[i].split(method['sep_sub'])
                        first = True
                        count = 0
                        for element_sub in element_subs:
                            count += 1
                            if first:
                                first = False
                                elements[i] = " '''" + str(count) + ".''' " + \
                                    element_sub.strip()
                            else:
                                elements[i] = elements[i] + " [[BR]] '''" + \
                                    str(count) + ".''' " + element_sub.strip()
                setattr(self, prop, elements[i].strip())
        else:
            if abort_on_fail:
                raise LogParseError(self.revision, self.original_log)
            else:
                self.revision = 0
                print "[WARN] Failed to parse revision r%s" % (
                      self.revision) + ' "%s" - ignoring changeset' % (message)

    def __str__(self):
        return 'r%s' % (self.revision)

    __repr__ = __str__


class LogParseError(Exception):
    """Exception when parsing the log fails."""
    def __init__(self, revision, message):
        self.revision = revision
        self.message = message
        super(LogParseError, self).__init__()

    def __str__(self):
        return "Parse failed for revision %s: " % (self.revision) + \
               "failed to parse log message '%s'" % (self.message)

    __repr__ = __str__


class AttributeNotFoundError(Exception):
    """Exception when trying to access an attribute which doesn't exist."""
    def __init__(self, attrib, item):
        self.attrib = attrib
        self.item = item
        super(AttributeNotFoundError, self).__init__()

    def __str__(self):
        return "Attempted to access parser attribute '%s' " % (self.attrib) +\
               "which does not exist in '%s'" % (self.item)

    __repr__ = __str__


class NoTrunkRevisionsError(Exception):
    """Exception when no trunk revisions are found."""
    def __init__(self, first, last, trunk):
        self.first = first
        self.last = last
        self.trunk = trunk
        super(NoTrunkRevisionsError, self).__init__()

    def __str__(self):
        return "No trunk revisions found on " + \
            "%s between %s and %s" % (self.trunk, self.first, self.last)

    __repr__ = __str__


class TrunkParseFailed(Exception):
    """Exception when fcm log on specified trunk fails."""
    def __init__(self, cmd, retc, out, err):
        self.cmd = cmd
        self.retc = retc
        self.out = out
        self.err = err
        super(TrunkParseFailed, self).__init__()

    def __str__(self):
        return 'Problem running "%s":\nReturn Code: %s\n' % (self.cmd,
                 self.retc) + 'Stdout: %s\nStderr: %s' % (self.out, self.err)

    __repr__ = __str__


def _run_command(command):
    '''Given a command as a string, run it and return the exit code,
    standard out and standard error. The optional shell argument allows
    a shell to be spawned to allow multiple commands to be run.'''

    # Turn command into a list
    command_list = command.split()

    # Create the Popen object and connect out and err to pipes
    proc = subprocess.Popen(command_list, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # Do the communicate and wait to get the results of the command
    stdout, stderr = proc.communicate()
    retcode = proc.wait()

    # Reformat stdout
    stdout = ''.join(stdout)
    stdout = stdout.split("\n")

    return retcode, stdout, stderr

############################# END OF LIBRARY ##################################

# Constants and subroutines associated with program execution from hereon.



CALENDAR = {
    '1': 'January',
    '2': 'February',
    '3': 'March',
    '4': 'April',
    '5': 'May',
    '6': 'June',
    '7': 'July',
    '8': 'August',
    '9': 'September',
    '10': 'October',
    '11': 'November',
    '12': 'December',
}


TIMES_OF_DAY = {
    'Small Hours (Midnight - 7am)': [0, 7],
    'Early Morning (7am - 9am)': [7, 9],
    'Morning (9am - Noon)': [9, 12],
    'Lunchtime (Noon - 2pm)': [12, 14],
    'Afternoon (2pm - 5pm)': [14, 17],
    'Evening (5pm - 7pm)': [17, 19],
    'Night (7pm - Midnight)': [19, 24],
}


GROUPS = {
    'UM sys': ['Glenn Greed', 'Stuart Whitehouse', 'Sean Swarbrick',
               'Richard Barnes', 'Paul Cresswell', 'Joe Mancell',
               'Steve Wardle', 'Ricky Wong', 'Roddy Sharp', 'Sam Cusworth',
               'Stuart Whitehouse and Joe Mancell', 'paulcresswell',
               'Steve Wardle and Ricky Wong', 'Steve wardle', 'UM System Team',
               'stuartwhitehouse', 'glenngreed', 'richardbarnes', 'joemancell',
               'stevewardle', 'roddysharp', 'samcusworth', 'Enrico Olivier',
               'enricoolivier', 'Ricky Olivier', 'Rich Gilham',
               'richardgilham', 'Richard Gilham'
              ],
}


#IGNORE_CODE_AREAS_IN_STATISTICS = ['ext_kgo_update', ]
IGNORE_CODE_AREAS_IN_STATISTICS = []

TECHNICAL_GROUPS= ['rose_stem', 'rose-stem', 'meta-data', 'meta_data', 
                   'fcm_make', 'fcm-make', 'optimisation', 'coupling',
                   'technical', 'utils', 'fieldcalc', 'mule', 'scripts']

SCIENCE_GROUPS= ['radiation', 'convection', 'stochastic_physics', 
                   'lsp_cloud', 'dynamics', 'ukca', 'bl_jules', 'gwd',
                   'idealised']

SYNONYM_WARNINGS = []


class InvalidParsingMethodError(Exception):
    """Exception when --method specified is not valid."""
    def __init__(self, method):
        self.method = method

    def __str__(self):
        return "The trunk log parsing method '%s' does not exist."%(self.method)

    __repr__ = __str__


def _parse_options():
    """Parse command line options."""
    usage = "Usage: %prog [options] <from-revision> <to-revision>"
    parser = OptionParser(usage=usage)
    parser.add_option("--first", dest='first', action="store",
                      help="First trunk revision (argument takes precedence)")
    parser.add_option("--last", dest='last', action="store",
                      help="Last trunk revision  (argument takes precedence)")
    parser.add_option("--method", dest='method', action="store",
                      help="Mapping method of log", default='UM')
    parser.add_option("--trunk", dest='trunk', action="store",
                      default="fcm:um.x_tr",
                      help="URL of trunk to search")
    parser.add_option("--abort-on-fail", dest="abort", action="store_true",
                      help="Abort when parsing fails rather than ignoring" +
                      " changeset")
    parser.add_option("--short", "-s", dest="short", action="store_true",
                      help="Print commits/ticket information only")
    parser.add_option("--wiki", "-W", dest="wiki", action="store_true",
                      help="Print wiki-formatted log")
    (opts, args) = parser.parse_args()
    if len(args) == 2:
        opts.first = args[0]
        opts.last = args[1]

    if not opts.first or not opts.last:
        parser.error("First and last revision are mandatory arguments")

# If we give it a keyword as the first argument, we don't want to inclusively
# include that log message in the parsing, so convert to a revision number and
# add one.
    if re.search(r'^\D', opts.first):
        stdout = _run_command("%s kp %s" % (FCM, opts.trunk))[1]
        for line in stdout:
            if re.search(r':%s\s*\]' % (opts.first), line):
                original = opts.first
                opts.first = int(re.sub(r'.*=\s*', r'', line)) + 1
                if not opts.short:
                    print "[WARN] Changing first revision to" +\
                          " %s+1 (r%s)" % (original, opts.first)

    thismodule = sys.modules[__name__]
    if not hasattr(thismodule, opts.method):
        raise InvalidParsingMethodError(opts.method)
    opts.method = getattr(thismodule, opts.method)

    return opts


def frequency_information(parser):
    """Prints information about the most and least frequent values for
    a variety of attributes."""
    frequent, members = parser.most_frequent('severity', 'commits')
    if len(frequent) == 1:
        print "The most common severity is: %s (%s commits)" % (
            ", ".join(frequent), len(members) / len(frequent))
    else:
        print "The most common severities are: %s (%s commits each)" % (
            ", ".join(frequent), len(members) / len(frequent))

    frequent, members = parser.most_frequent('code_area', 'tickets')
    if len(frequent) == 1:
        print "The most common code area is: %s (%s tickets)" % (
            ", ".join(frequent), len(members) / len(frequent))
    else:
        print "The most common code areas are: %s (%s tickets each)" % (
            ", ".join(frequent), len(members) / len(frequent))

    frequent, members = parser.most_frequent('reviewer', 'commits')
    if len(frequent) == 1:
        print "The most frequent committer to the trunk is" + \
            ": %s (%s commits)" % (", ".join(frequent), len(members) /
            len(frequent))
    else:
        print "The most frequent committers to the trunk" + \
            " are: %s (%s commits each)" % (", ".join(frequent), len(members) /
            len(frequent))

    frequent, members = parser.least_frequent('reviewer', 'commits')
    if len(frequent) == 1:
        print "The least frequent committer to the trunk" + \
            " is: %s (%s commits each)" % (", ".join(frequent), len(members) /
            len(frequent))
    else:
        print "The least frequent committers to the trunk" + \
            " are: %s (%s commits each)" % (", ".join(frequent), len(members) /
            len(frequent))
    print "-"*80


def problem_tickets(parser):
    """Prints information about problem tickets which required more than
    one commit."""
    largest, members = parser.largest('commits', 'tickets')
    print "The largest number of commits per ticket is: %s" % (
          ", ".join(largest))
    print "-"*80
    tickets_by_num_commits = parser.group('commits', 'tickets')
    for i in range(int(largest[0]), 1, -1):
        if str(i) in tickets_by_num_commits:
            tickets = tickets_by_num_commits[str(i)]
            print "%s tickets with %s commits:" % (len(tickets), i)
            for ticket in tickets:
                print "  Ticket %s by %s, " % (ticket,
                      ticket.author) + \
                    "reviewed by %s" % (ticket.reviewer)
                print '    "%s"' % (ticket.initial_log)
            print "-"*80


def trivial_and_non_trivial_tickets(parser):
    """Print information about trivial and non-trivial tickets."""
    members = parser.group('severity', 'tickets', {'severity': 'trivial'})
    count_members = sum(len(v) for v in members.itervalues())
    if count_members > 0:
        print "There were %s trivial tickets" % (count_members)
        auth, tickets = parser.most_frequent('author', 'tickets', {
            'severity': 'trivial'}, invert=False)
        print "The most productive trivial authors were: %s (%s tickets)" % (
            ", ".join(auth), len(tickets))
        print "-"*80

    notmembers = parser.group('severity', 'tickets', {'severity': 'trivial'},
                              invert=True)
    count_notmembers = sum(len(v) for v in notmembers.itervalues())

    if count_notmembers > 0:
        print "There were %s non-trivial tickets" % (count_notmembers)
        auth, tickets = parser.most_frequent('author', 'tickets', {
                                             'severity': 'trivial'},
                                             invert=True)
        print "The most productive non-trivial authors were: %s (%s tickets)"\
              % (", ".join(auth), len(tickets))
        print "-"*80

def weekday_counter(fromdate, todate):
    daygenerator = (fromdate + timedelta(x + 1) for x in xrange((todate - fromdate).days))
    working_days = sum(1 for day in daygenerator if day.weekday() < 5)
    return working_days

def commits_per_ticket_info(parser, short=False):
    """Print statistics about commits and tickets."""
    num_commits = float(len(parser.commits))
    num_tickets = float(len(parser.tickets))

    tickets_by_num_commits = parser.group('commits', 'tickets')
    if '1' in tickets_by_num_commits:
        tickets_with_single_commit = len(tickets_by_num_commits['1'])
    else:
        tickets_with_single_commit = 0

    if '2' in tickets_by_num_commits:
        tickets_with_two_commits = len(tickets_by_num_commits['2'])
    else:
        tickets_with_two_commits = 0

    tickets_with_multiple_commits = num_tickets - tickets_with_single_commit
    tickets_with_more_than_two_commits = tickets_with_multiple_commits - \
                                         tickets_with_two_commits
    percentage_tickets_with_multiple_commits = tickets_with_multiple_commits \
                                               / num_tickets * 100.0
    percentage_tickets_with_more_than_two_commits = \
         tickets_with_more_than_two_commits / num_tickets * 100.0
    num_commits_ignore = 0
    num_tickets_ignore = 0
    for area in IGNORE_CODE_AREAS_IN_STATISTICS:
        area_commits = check_for_synonyms(parser.group('code_area', 'commits'))
        area_tickets = check_for_synonyms(parser.group('code_area', 'tickets'))
        if area in area_commits:
            matching_commits = len(area_commits[area])
            matching_tickets = len(area_tickets[area])
            num_commits_ignore += matching_commits
            num_tickets_ignore += matching_tickets

    commits_per_ticket = float(num_commits - num_commits_ignore) / \
        float(num_tickets - num_tickets_ignore)
    commits_per_ticket_squared = commits_per_ticket / (float(num_tickets - num_tickets_ignore))
    delta_time = parser.commits[-1].datetime - parser.commits[0].datetime
    days = float(delta_time.days + 1)    # Round up

    work_days = weekday_counter(parser.commits[0].datetime, parser.commits[-1].datetime)
    if work_days < 1:
        work_days = 1

    if short:
        # How buggy the release is (using bugs found during the release)
        print "{0:.3f} commits per ticket".format(commits_per_ticket)
        # Attempt to normalise the above for size of release
        print "{0:.3f} milli-commits per square ticket".format(commits_per_ticket_squared*1000)
        # A measure of risk we're accepting on the trunk
        print "{0:.3f}% of tickets had multiple commits".format(
                                     percentage_tickets_with_multiple_commits)
    else:
        print "%s tickets and %s commits total in" % (int(num_tickets),
              int(num_commits)), delta_time
        print "Release lasted %s working days" % (work_days)
        if len(IGNORE_CODE_AREAS_IN_STATISTICS) > 0:
            print "Ignoring code areas: %s" % (','.join(
                IGNORE_CODE_AREAS_IN_STATISTICS))
            print "%s tickets and %s " % (num_commits_ignore,
                  num_commits_ignore) + "commits ignored in the commits per " \
                  + "ticket statistic"
        print "{0:3f} commits per calendar day, {1:.3f} tickets per calendar day".format(
              num_commits/days, num_tickets/days)
        print "{0:3f} commits per working day, {1:.3f} tickets per working day".format(
              num_commits/work_days, num_tickets/work_days)
        print "-"*80
        print "{0:.3f} commits per ticket".format(commits_per_ticket)
        print "{0:.3f} milli-commits per square ticket".format(commits_per_ticket_squared*1000)
        print "{0:.3f}% of tickets had multiple commits".format(
                                percentage_tickets_with_multiple_commits)
        print "{0:.3f}% of tickets had more than two commits".format(
                                percentage_tickets_with_more_than_two_commits)
        print "-"*80


def check_for_synonyms(group):
    """Check for synonyms in group returned by parser, alert user, and
    collapse."""
    to_delete = []
    for i in group:
        j = re.sub(r'_', r'-', i)
        if j in group and i != j:
            group[j] += group[i]
            if i not in SYNONYM_WARNINGS or j not in SYNONYM_WARNINGS:
                SYNONYM_WARNINGS.append(i)
                SYNONYM_WARNINGS.append(j)
                print "[WARN] %s and %s treated as synonymous" % (i, j)
            to_delete.append(i)
    for i in group:
        j = i.lower()
        if j in group and i != j:
            group[j] += group[i]
            if i not in SYNONYM_WARNINGS or j not in SYNONYM_WARNINGS:
                SYNONYM_WARNINGS.append(i)
                SYNONYM_WARNINGS.append(j)
                print "[WARN] %s and %s treated as synonymous" % (i, j)
            to_delete.append(i)
    for i in to_delete:
        del(group[i])
    return group


def code_areas_by_ticket(parser):
    """Print number and membership of tickets in each code area."""
    code_areas = check_for_synonyms(parser.group('code_area', 'tickets'))
    
    science_tickets = []
    technical_tickets = []
    orphaned_tickets = []
    orphaned_groups = {}
    
    for i in code_areas:
        print "%s: %s tickets %s" % (i, len(code_areas[i]), code_areas[i])
        if i in TECHNICAL_GROUPS:
            technical_tickets = technical_tickets + code_areas[i]
        elif i in SCIENCE_GROUPS:
            science_tickets = science_tickets + code_areas[i]
        else:
            orphaned_groups[i] = True
            orphaned_tickets = orphaned_tickets + code_areas[i]

    print "-"*80
    print "%s technical tickets"%(len(technical_tickets))
    print "%s science tickets"%(len(science_tickets))
    print "%s orphaned tickets"%(len(orphaned_tickets)), orphaned_groups.keys()
    print "-"*80
    return len(technical_tickets), len(science_tickets)


def tickets_per_month(parser):
    """Print number of tickets in each month."""
    months = {}
    order = []
    for i in parser.tickets:
        k = CALENDAR[str(i.initial_datetime.month)]
        if k in months:
            months[k] += 1
        else:
            months[k] = 1
    for j in months:
        print "%s: %s tickets" % (j, months[j])

    print "-"*80


def times_of_day(parser):
    """Print information on commits by time of day."""
    times = {}
    for i in parser.commits:
        for name, j in TIMES_OF_DAY.iteritems():
            if i.datetime.hour >= j[0] and i.datetime.hour < j[1]:
                if name in times:
                    times[name].append(i)
                else:
                    times[name] = [i]

    for i in times:
        print "%s: %s commits" % (i, len(times[i]))

    print "-"*80


def tickets_by_type(parser):
    """Print information about tickets by type."""
    ticket_types = check_for_synonyms(parser.group('ticket_type', 'tickets'))
    for i, j in ticket_types.iteritems():
        print "%s: %s tickets" % (i, len(j))
    print "-"*80


def commits_by_week(parser):
    """Print number of commits by week of year number."""
    week_nums = {}
    order = []

    for i in parser.commits:
        k = i.datetime.isocalendar()[1]
        if k in week_nums:
            week_nums[k].append(i)
        else:
            week_nums[k] = [i]
        if k not in order:
            order.append(k)
    for week in order:
        print "Week %s: %s commits" % (week, len(week_nums[week]))
    print "-"*80
#    commits_by_week_excluding_umsys(parser, week_nums)


def commits_by_week_excluding_umsys(parser, total_week_nums=None):
    """Print number of commits by week of year number."""
    week_nums = {}
    commits = parser.group('author', 'commits')
    for author in commits:
        if author in GROUPS['UM sys']:
            pass
        else:
            for i in commits[author]:
                k = i.datetime.isocalendar()[1]
                if k in week_nums:
                    week_nums[k].append(i)
                else:
                    week_nums[k] = [i]
    print "Commits by week excluding UM sys authors:"
    for week in sorted(week_nums.keys()):
        if total_week_nums is None:
            print "Week %s: %s commits" % (week, len(week_nums[week]))
        else:
            print "Week %s: %s commits (%s%%)" % (week, len(week_nums[week]),
                                                 float(len(week_nums[week]))
                                                 / float(len(
                                                 total_week_nums[week]))
                                                 * 100)
    print "-"*80


def break_regression(parser):
    """Print tickets which broke regression."""
    tickets = check_for_synonyms(parser.group('regression', 'tickets',
                                 {'regression': 'regression'}, invert=True))
    if tickets:
        for i in tickets:
            print "%s: %s tickets %s" % (i, len(tickets[i]), tickets[i])
        print "The remaining tickets maintained regression."
    else:
        print "All tickets maintained regression."
    print "-"*80


def count_by_group(parser):
    """Count tickets by group."""
    group = {}
    multi = []
    for i in parser.tickets:
        found = False
        for grp, authors in GROUPS.iteritems():
            if i.author in authors:
                if grp in group:
                    group[grp].append(i)
                else:
                    group[grp] = [i]
                found = True
        if not found:
            multi.append(i)
            
    groups_by_number = {}
    for grp in group:
        num_tickets = len(group[grp])
        groups_by_number[grp] = num_tickets
    groups_in_order = sorted(groups_by_number, key=groups_by_number.get)

    for grp in reversed(groups_in_order):
        print "%s: %s tickets %s" % (grp, len(group[grp]), group[grp])

    text = 'Multiple authors: '
    for j in multi:
        text = text + "%s (%s) " % (j, j.author)
    if len(multi) > 0:
        print text
    print "-"*80
    return len(group['UM sys'])


def count_umsys(parser):
    """Count number of tickets from UM sys."""
    group = {}
    for i in parser.tickets:
        found = False
        for grp, authors in GROUPS.iteritems():
            if i.author in authors:
                if grp in group:
                    group[grp].append(i)
                else:
                    group[grp] = [i]
                found = True
    return len(group['UM sys'])
    

def finalise(parser, tech, sci, umsys):
    """Final output which depends on other routines."""
    
    num_tickets = float(len(parser.tickets))

    print "{0:.3f}% of tickets are science".format(sci/num_tickets*100.0)
    print "{0:.3f}% of tickets are technical".format(tech/num_tickets*100.0)
    print "{0:.3f}% of tickets come from the UM system team ({1} tickets)".format(
        umsys/num_tickets*100.0, umsys)
    enhancements = parser.group('ticket_type', 'tickets')['enhancement']
    science_enhancements = []
    for i in enhancements:
        if i.code_area in SCIENCE_GROUPS and i.severity != 'trivial':
          science_enhancements.append(i)
    print "{0:.3f}% of tickets are non-trivial science enhancements ({1} tickets)".format(
           len(science_enhancements)/num_tickets*100, len(science_enhancements))
    print "-"*80


def unique_authors(parser):
    """Print number of unique authors. This doesn't cater for the difference
    between e.g. Tom and Thomas, but does deal with multiple ticket authors
    linked by "&" or "and"."""
    authors = {}
    numeq1 = 0
    numeq2 = 0
    num2t5 = 0
    numgt5 = 0
    numgt10 = 0
    for i in parser.tickets:
        auth = re.split("\s+and|&|/\s*",i.author)
        for a in auth:
            if a not in authors:
                authors[a] = 1
            else:
                authors[a] += 1
    for author in authors:
        if authors[author] == 1:
            numeq1 += 1
        if authors[author] == 2:
            numeq2 += 1
        if authors[author] > 2 and authors[author] < 6:
            num2t5 += 1
        if authors[author] > 5 and authors[author] < 11:
            numgt5 += 1
        if authors[author] > 10:
            numgt10 += 1
    print "There are {0} unique author names".format(len(authors))
    print "There are {0} unique authors who contributed a single ticket".format(
            numeq1)
    print "There are {0} unique authors who contributed 2 tickets".format(
            numeq2)
    print "There are {0} unique authors who contributed 3 to 5 tickets".format(
            num2t5)
    print "There are {0} unique authors who contributed 6 to 10 tickets".format(
            numgt5)
    print "There are {0} unique authors who contributed > 10 tickets".format(
            numgt10)
    print "-"*80
    return authors


def umsys_authors(parser):
    """Print UM System Team member authorship in descending order."""
    um_sys_authored = {}
    um_sys_reviewed = {}
    um_sys_reviewed_multi_commits = {}
    
    for i in parser.tickets:
        auth = re.split("\s*and|&|/\s*",i.author)
        for a in auth:
            if a in GROUPS["UM sys"]:
                if a in um_sys_authored:
                    um_sys_authored[a] += 1
                else:
                    um_sys_authored[a] = 1
        if i.reviewer in GROUPS["UM sys"]:
            if i.reviewer in um_sys_reviewed:
                um_sys_reviewed[i.reviewer] += 1
            else:
                um_sys_reviewed[i.reviewer] = 1
            if i.reviewer in um_sys_reviewed_multi_commits:
                if i.num_commits > 1:
                    um_sys_reviewed_multi_commits[i.reviewer] += 1
            else:
                um_sys_reviewed_multi_commits[i.reviewer] = 0
    
    for j in reversed(sorted(um_sys_authored, key=um_sys_authored.get)):
        print "%s: %s tickets authored"%(j, um_sys_authored[j])
    print "-"*80
    for j in reversed(sorted(um_sys_reviewed, key=um_sys_reviewed.get)):
        print "%s: %s tickets reviewed (%s with multiple commits)"%(j, 
            um_sys_reviewed[j], um_sys_reviewed_multi_commits[j])
    print "-"*80


def main():
    """Main program block, providing an example of how to use the
    TrunkParser class."""

    # Parse command line options
    opts = _parse_options()

    # Print information about what we're going to parse
    if not opts.short:
        print "-"*80
        print "Parsing %s between %s and %s (inclusive)" % (opts.trunk,
                                                            opts.first,
                                                            opts.last)

    # Parse the trunk and create the parser
    parser = TrunkLogParser([opts.first, opts.last], trunk=opts.trunk,
                            method=opts.method, abort_on_fail=opts.abort)
    if not opts.short:
        print "-"*80

    if len(parser.commits) == 0:
        print "No commits detected in this range."
        return

    # Subroutines which print various things
    if opts.short:
        commits_per_ticket_info(parser, short=True)
    else:
        frequency_information(parser)
        problem_tickets(parser)
        trivial_and_non_trivial_tickets(parser)
        tech, sci = code_areas_by_ticket(parser)
        tickets_per_month(parser)
        times_of_day(parser)
        tickets_by_type(parser)
        commits_by_week(parser)
        break_regression(parser)
        unique_authors(parser)
        if 'um' in opts.trunk:
            umsys = count_umsys(parser)
            umsys_authors(parser)
            finalise(parser, tech, sci, umsys)
        commits_per_ticket_info(parser)


    if opts.wiki:
        print parser.trac_format_log()

if __name__ == '__main__':
    main()
