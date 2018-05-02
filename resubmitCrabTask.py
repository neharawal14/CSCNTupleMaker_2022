#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2014.09.14
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import time
from time import gmtime, strftime
import math

#define function for parsing options
def parseOptions():
    global observalbesTags, modelTags, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-t', '--taskDir', dest='TASKDIR', type='string',default='', help='crab task dir')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quite = 0):
    #    print cmd
    status, output = commands.getstatusoutput(cmd)
    if (status !=0 and not quite):
        print 'Error in processing command:\n   ['+cmd+']'
        print 'Output:\n   ['+output+'] \n'
    return output

def resubmit():

    # parse the arguments and options
    global opt, args
    parseOptions()

    taskdir = opt.TASKDIR

    for subdir in [os.path.join(taskdir, d) for d in os.listdir(taskdir)]:
      if (not os.path.isdir(subdir)): continue
      if (not 'crab_' in subdir): continue
      #if ((('Run2015B' in subdir) or ('50ns' in subdir))): continue
      print subdir
      cmd = 'crab resubmit -d '+str(subdir)
      #cmd = 'crab kill -d '+str(subdir)
      output = processCmd(cmd) 
      #print output

if __name__ == "__main__": 
    resubmit() 
