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
    parser.add_option('-n','--maxEvents', dest='maxEvents', type='int', default='1000' ,help='maxEvents default:1000')
    parser.add_option('-d','--outDir', dest='outDir', type='string', default='/home/msnowball/public_html/CSC/' ,help='out directory default:CSC')
    
    
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()
    


def doAnalysis():
  global opt, args
  
  tfile = ROOT.TFile(opt.file,"READ")
  if not tfile:
      raise RunTimeError,"No input file specified or root file could not be found!"

  if opt.isMC:
      tree = tfile.Get("cscRootMaker/Events")
  else:
      tree = tfile.Get("cscRootMaker/Events")



  hists1D = {}
  hists2D = {}
  defineHistos(hists1D,hists2D)
  
  #Analysis Loop
  for i in range( tree.GetEntries() ):
    tree.GetEntry(i)

    if i%100 == 0:
        print "Event ",i
    if i > opt.maxEvents:
        break


    #CSC Segments
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
            for m in range(len(tree.cscSegments_recHitRecord_layer[Identifier])):
                layer = tree.cscSegments_recHitRecord_layer[Identifier][m]
                nLayers+=1
            if nLayers > 0:
                hists1D['OneSegmentChambers_nRecHitLayers'].Fill(nLayers)
                


  #Write Histograms to .eps/.png
  writeHistos(hists1D,hists2D)



#Declare Hists
def defineHistos(Histos1D,Histos2D):

    Histos1D['OneSegmentChambers_nRecHitLayers'] = ROOT.TH1I("1DSegments_nLayers","; N Layers; N Events", 9,-0.5,8.5) 

#Write Hists to files
def writeHistos(Histos1D,Histos2D):

    ROOT.gROOT.ProcessLine(".L tdrstyle.cc")

    for key in Histos1D:
        c = ROOT.TCanvas("c","c",700,700)
        c.cd()
        Histos1D[key].Draw("HIST")
        c.SaveAs(opt.outDir+'/'+str(Histos1D[key].GetName())+'.eps')
        c.SaveAs(opt.outDir+'/'+str(Histos1D[key].GetName())+'.png')
        del c

        

#Main  
if __name__ == "__main__":
  parseOptions()
  doAnalysis()
  
