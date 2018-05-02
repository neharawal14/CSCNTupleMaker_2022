import os, sys, glob, fnmatch, re
from optparse import OptionParser

def find_files(directory, pattern):
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.join(root, basename)
        yield filename

def main():
   usage = "usage: %prog [options] arg"
   parser = OptionParser(usage)
   parser.add_option("-i", "--indir", dest="indir",
                     help="input directory with .root files")
   parser.add_option("-o", "--outdir",
                     help="output directory", default=".", dest="outdir")
   (options, args) = parser.parse_args()

   indir = str(options.indir)
   outdir = options.outdir

   print 'Input Directory is '+indir
   print 'Output Directory is '+outdir

   samples = []

   for filename in find_files(indir, '*.root'):
     
     #if ('MuonEG_Run2015C_' in filename): continue
     if (not 'DoubleEG_Run2015C_50ns' in filename): continue
     
     fullsample = filename.split('/')
     sampledir = ''
     for i in range(len(fullsample)-1):
       sampledir+=fullsample[i]+'/'
     sample = fullsample[len(fullsample)-1]
     sample = re.sub('_[0-9]*.root','',sample)

     if (not sample in samples): samples.append(sample)
     else: continue

     print sample

     os.system('hadd -f ' + outdir + '/' + sample +'.root ' + sampledir + '/' + sample + '_*' + '.root')
     
if __name__ == "__main__":
   main()
