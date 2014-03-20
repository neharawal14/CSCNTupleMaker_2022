#!/usr/bin/python
import sys, os, pwd, commands
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array


def parseOptions():
    
    usage = ('usage: %prog [options] \n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
    parser.add_option('-m','--isMC', dest='isMC', type='int', default=0 ,help='isMC default:0')
    parser.add_option('-f','--file', dest='file', type='string', default='cscRootMaker.root' ,help='file default:blank')
    parser.add_option('-n','--maxEvents', dest='maxEvents', type='int', default=1000 ,help='maxEvents default:1000')
    parser.add_option('-d','--outDir', dest='outDir', type='string', default='/home/msnowball/public_html/CSC/' ,help='out directory default:CSC')
    parser.add_option('-j','--jobName',dest='jobName',type='string', default='cscAna',help='name of job and output files')

    parser.add_option('--isDigi', dest='isDigi', type='int', default=0 ,help='isDigi default:0')
    parser.add_option('--isLocalReco', dest='isLocalReco', type='int', default=1 ,help='isLocalReco default:1')
    parser.add_option('--isFullReco', dest='isFullReco', type='int', default=1 ,help='isFullReco default:1')
        
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()
    


class Analysis():

    def __init__(self):

        self.hists1D = {}
        self.hists2D = {}
        self.totalEvents = 0

        self.defineHistos()
        

    def doAnalysis(self,file):
        global opt, args
  
        tfile = ROOT.TFile(file,"READ")
        if not tfile:
            raise RunTimeError,"No input file specified or root file could not be found!"
        
        if opt.isMC:
            tree = tfile.Get("cscRootMaker/Events")
        else:
            tree = tfile.Get("cscRootMaker/Events")


        #Analysis Loop
        for i in range( tree.GetEntries() ):
            tree.GetEntry(i)

            if i%100 == 0:
                print "Event ",i
            if self.totalEvents > opt.maxEvents:
                break
            self.totalEvents+=1

            if opt.isLocalReco or opt.isFullReco:
            #CSCSegments
                CSC_SegmentCounter = [0] * 600
                CSC_SegmentMap = [''] * 600
                for n in range(tree.cscSegments_nSegments):
                    serialRecord = tree.cscSegments_ID_chamberSerial[n]
                    CSC_SegmentCounter[serialRecord] = CSC_SegmentCounter[serialRecord]+1
                    CSC_SegmentMap[serialRecord]+=str(n)+'.'

                for n in range(0,600):
                    if CSC_SegmentCounter[n] == 1:
                        line = CSC_SegmentMap[n].split('.')
                        Identifier=int(line[0])
                        nLayers=0
                        missingLayers = []
                        for m in range(len(tree.cscSegments_recHitRecord_layer[Identifier])):
                            layer = tree.cscSegments_recHitRecord_layer[Identifier][m]
                            missingLayers.append(layer)
                            nLayers+=1
                        if nLayers > 0:
                            self.hists1D['OneSegmentChambers_nRecHitLayers_Norm'].Fill(nLayers)
                            #if not all 6 layers are reconstructed, plot missing layers
                            if nLayers < 6:
                                for j in range(len(missingLayers)):
                                    self.hists1D['OneSegmentChambers_missingLayer_Norm'].Fill(missingLayers[j])
                                    







    def defineHistos(self):
        
        self.hists1D['OneSegmentChambers_nRecHitLayers_Norm'] = ROOT.TH1F("1SegmentChambers_nLayers","; N Layers; Fraction of Events", 9,-0.5,8.5) 
        self.hists1D['OneSegmentChambers_missingLayer_Norm'] = ROOT.TH1F("1SegmentChambers_missingLayer","; Missing Layer(s); Fraction of Events", 9,-0.5,8.5) 


    def writeHistos(self,Histos1D,Histos2D):
        
        ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
        setTDRStyle(False)
        for key in Histos1D:
            c = ROOT.TCanvas("c","c",700,700)
            c.cd()
            normalized = 'Norm' in key
            if normalized:
                Histos1D[key].Scale(1/Histos1D[key].Integral())
            Histos1D[key].Draw("HIST")
            c.SaveAs(opt.outDir+'/'+str(Histos1D[key].GetName())+'.eps')
            c.SaveAs(opt.outDir+'/'+str(Histos1D[key].GetName())+'.png')
            del c


    def writeHistosToRoot(self,Histos1D,Histos2D):
        
        ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
        setTDRStyle(False)
        outFile = ROOT.TFile(opt.jobName+'.root',"RECREATE")
        
        for key in Histos1D:
            normalized = 'Norm' in key
            if normalized:
                Histos1D[key].Scale(1/Histos1D[key].Integral())
            outFile.cd()
            Histos1D[key].Write()
        outFile.Write()
        outFile.Close()




    def endjob(self,singleFile):

        if self.totalEvents > 0:
            if singleFile:
                self.writeHistos(self.hists1D,self.hists2D)
            else:
                self.writeHistosToRoot(self.hists1D,self.hists2D)
        




#Main  
if __name__ == "__main__":

    global opt, args
    parseOptions()

    myClass = Analysis()
    singleFile = False

    if opt.file.endswith(".root"):
        singleFile = True
    elif opt.file.endswith(".txt"):
        singleFile = False
    else:
        raise RuntimeError, "opt.file: file name does not end with .root or .txt!"


    # Loop for parallel or single file 
    if singleFile:
        myClass.doAnalysis(opt.file)
    else:
        lines = open(opt.file,"r")
        for line in lines:
            f = line.split()
            if not f[0].endswith(".root"): continue
            if len(f) < 1: continue
            print "Opening file",f[0]
            myClass.doAnalysis(f[0])
            

    myClass.endjob(singleFile)
