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
    parser.add_option('-t', '--tag', dest='TAG', type='string',default='', help='tag to be appended to the results, default is an empty string')
    parser.add_option('-d', '--datasets', dest='DATASETS', type='string', default='datasets_Fall15_25ns_MiniAODv1.txt', help='txt file with datasets to run over')
    parser.add_option('-c', '--cfg', dest='CONFIGFILE', type='string', default='/scratch/osg/mhl/CSC_2016/localRecoStudy/CMSSW_6_2_0_SLHC23_patch1/src/UFCSCSoftware/UFCSCRootMaker/test/UFCSCRootMaker_template.py', help='configuration template')

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

def submitAnalyzer():

    # parse the arguments and options
    global opt, args
    parseOptions()

    cfgtemplate = opt.CONFIGFILE

    # save working dir
    currentDir = os.getcwd()

    tag = opt.TAG
    outDir='resultsAna_'+tag 

    if (not os.path.isdir(outDir)):
        cmd = 'mkdir '+outDir
        processCmd(cmd)
        cmd = 'mkdir '+outDir+'/cfg/'
        processCmd(cmd)

    # get the datasets
    print '[Gathering Dataset Information]'
    datasets = []
    cross_section = {}
    nfiles = {}
    nevents = {}
    datasetfiles = {}

    with open(opt.DATASETS, "r") as datasetfile:
        for line in datasetfile:

            if (line.startswith('#')): continue

            #if ('ZZTo4L_13TeV_powheg' in line): continue

            dataset = line.split()[0]
            dataset = dataset.rstrip()
            dataset = dataset.lstrip()

            datasets.append(dataset)
            cross_section[dataset] = float(line.split()[1])

            cmd = './das_client.py --query="file dataset='+dataset+'" --limit=0 | grep ".root"'
#            cmd = './das_client.py --query="file dataset='+dataset+' instance=prod/phys03" --limit=0 | grep ".root"'

            output = processCmd(cmd)
            while ('error' in output):
                time.sleep(1.0);
                output = processCmd(cmd)
            datasetfiles[dataset] =  output.split()
            nfiles[dataset] = len(datasetfiles[dataset])
 
            cmd = './das_client.py --query="dataset dataset='+dataset+' | grep dataset.nevents" --limit=0'
            output = processCmd(cmd)
            while ('error' in output):
                time.sleep(1.0);
                output = processCmd(cmd)
            nevents[dataset] = output

#            print cross_section, nfiles, nevents 
            print dataset,'xs:',cross_section[dataset],'nfiles:',nfiles[dataset],'nevents:',nevents[dataset]


    # submit the jobs
    print '[Submitting jobs]'
    jobCount=0


    for dataset in datasets:

        filename = dataset.split('/')[1]+'_'+dataset.split('/')[2]
        if (len(filename)>99):
          newfilename = filename.split('TuneZ2star')[0]+'PU140'
          filename = newfilename

        cfgfile = dataset.lstrip('/')
        cfgfile = cfgfile.replace('/','_')+'.py'

        cmd = 'cp '+cfgtemplate+' '+outDir+'/cfg/'+cfgfile
        output = processCmd(cmd)

        filelist = ''
        for f in range(0,5):
            if ((f+1)>nfiles[dataset]): continue
            filelist += '"'
            if (os.path.isfile('/cms/data'+datasetfiles[dataset][f])):
                filelist += 'file:/cms/data'
            filelist += datasetfiles[dataset][f]
            filelist += '",'
        filelist = filelist.rstrip(',')

        cmd = "sed -i 's~DUMMYFILELIST~"+filelist+"~g' "+outDir+'/cfg/'+cfgfile
        output = processCmd(cmd)

        cmd  = "sed -i 's~DUMMYFILENAME~"+filename+"~g' "+outDir+'/cfg/'+cfgfile
        output = processCmd(cmd)


        if (('PromptReco-v4' in cfgfile) or ('Run2015C_25ns-05Oct' in cfgfile) or ('16Dec' in cfgfile)):
            cmd = "sed -i 's~\"PAT\"~\"RECO\"~g' "+outDir+'/cfg/'+cfgfile
            output = processCmd(cmd)

        if (('MCRUN2_74' in cfgfile)):
            cmd = "sed -i 's~slimmedAddPileupInfo~addPileupInfo~g' "+outDir+'/cfg/'+cfgfile
            output = processCmd(cmd)

        crabcfgfile = 'crabConfig_'+filename+'.py'             

        cmd = 'cp crabConfig_TEMPLATE.py '+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~JOBTAG~"+tag+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~CFGFILE~"+outDir+"/cfg/"+cfgfile+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~OUTFILENAME~"+filename+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~DATASETNAME~"+dataset+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = 'crab submit -c '+outDir+'/cfg/'+crabcfgfile
        print cmd     

        output = processCmd(cmd) 
        while ('error' in output): 
            time.sleep(1.0); 
            output = processCmd(cmd) 
            if ('error' not in output):
                print 'Submitted after retry'
        print output

# run the submitAnalyzer() as main() 
if __name__ == "__main__": 
    submitAnalyzer() 
