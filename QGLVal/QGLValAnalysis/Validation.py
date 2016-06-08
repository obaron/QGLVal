

import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob

from os.path import join,exists

print 'Python version', sys.version_info
if sys.version_info < (2, 7):
    raise "Must use python 2.7 or greater. Have you forgotten to do cmsenv?"

workdir = 'work'
fileListDir = join(workdir,'files')

#define samples paths
lPath="/afs/cern.ch/user/o/obaron/QGLVal/CMSSW_7_6_3/src/QGLVal/QGLValAnalysis" #local path
#lPath = "/mnt/t3nfs01/data01/shome/grauco/JetMET/CMSSW_7_6_3_patch2/src/ttDM/ne"
#path = "/scratch/grauco/B2G_25ns/"
#path = "/scratch/decosa/07Oct2015/"
t3Path = '//pnfs/psi.ch/cms/trivcat/store/user/grauco/JetMet_76X_v7/'
#t3Path = '/store/user/grauco/B2G_Fw76X/'
#t3Ls = 'xrdfs  xrootd-cms.infn.it ls -u'
t3Ls = 'xrdfs t3dcachedb03.psi.ch ls -u'

#define samples, one folder for each mass value
samples = []
#samples.append("QCD_pythia")
#samples.append("QCD_herwig")

#samples.append("TEST")
#samples.append("QCD")
#samples.append("ZeroBias") 
samples.append("QCD_Pt50to80")
samples.append("QCD_Pt85to120")
samples.append("QCD_Pt600to800")
samples.append("QCD_Pt800to1000")
samples.append("QCD_Pt15to30")
samples.append("QCD_Pt30to50")
samples.append("QCD_Pt120to170")
samples.append("QCD_Pt170to300")
samples.append("QCD_Pt300to470")
samples.append("QCD_Pt470to600")
samples.append("QCD_Pt1000to1400")
samples.append("QCD_Pt1400to1800")
samples.append("QCD_Pt1800to2400")
samples.append("QCD_Pt2400to3200")

usage = 'usage: %prog -l lumi'
parser = optparse.OptionParser(usage)
parser.add_option('-c', '--channel', dest='channel', type='string', default = 'singleH', help='Channel to analyze: singleH or singleZ')
parser.add_option('-C', '--cat', dest='cat', type='string', default = 'cat2', help='Category to analyze: cat0 or cat1 or cat2')
parser.add_option('-s', '--sys', dest='sys', type='string', default = 'noSys', help='Systematics')
parser.add_option('', '--sync', dest='sync', type='string', default = 'noSync', help='Synchro exercise')
parser.add_option('-d', '--isData', dest='isData', type='string', default = 'MC', help='is Data or MC?')
parser.add_option('-g','--gdb', dest='gdb', action='store_true', default=False)
parser.add_option('-n','--dryrun', dest='dryrun', action='store_true', default=False)
parser.add_option('-m','--mode', dest='mode', default='local', choices=['local','t3se'])
parser.add_option('--t3batch', dest='t3batch', action='store_true', default=False)

isData="MC"
(opt, args) = parser.parse_args()

if opt.sys not in ["noSys", "jesUp", "jesDown", "jerUp", "jerDown", "metUnclUp", "metUnclDown"]:
    parser.error('Please choose an allowed value for sys: "noSys", "jesUp", "jesDown", "jerUp", "jerDown","metUnclUp", "metUnclDown"')

# Create working area if it doesn't exist
if not exists(fileListDir):
    os.makedirs(fileListDir)

for s in samples:
    if (s.startswith("Zero") or s.startswith("ZeroBias") or  s.startswith("Jet")): isData="DATA"
    print s
    print str(opt.sync)
    #print cmd
    #os.system(cmd)
    if opt.mode == 'local':
        print 'Local mode not yet supported'

        sPath = join(lPath,s,'*.root')
        # print ' '.join([lLs,sPath])
        # Get the complete list of files
        # listing = subprocess.check_output(lLs.split()+[sPath])
        files = glob.glob(sPath)
        print 'Sample',s,'Files found',len(files)

    elif opt.mode == 't3se':
        # Build the full path of sample files
        sT3Path = join(t3Path,s)
        print ' '.join([t3Ls,sT3Path])

        # Get the complete list of files
        listing = subprocess.check_output(t3Ls.split()+[sT3Path])
        files = listing.split()
        print 'Sample',s,'Files found',len(files)

    # Save it to a semi-temp file
    sampleFileList = join(fileListDir,s+'.txt')
    print sampleFileList
    with open(sampleFileList,'w') as sl:
        sl.write('\n'.join(files))

    outDirs = ['res_2','trees']

    for d in outDirs:
        if exists(d): continue
        os.makedirs(d)

    cmd = 'Validation '+ s + ' ' + sampleFileList  + ' ' + opt.channel + ' ' + opt.cat + ' ' + opt.sys + ' ' + opt.sync + ' ' + isData
#    cmd = "Validation "+ s + " " + path + s  + " " + opt.channel + " " +opt.cat + " " + opt.sys + " " + opt.sync + " " + isData

    if opt.gdb:
        cmd = 'gdb --args '+cmd
    elif opt.t3batch:
        jid = '%s_%s_%s_%s' % (s,opt.channel,opt.cat,opt.sys)
        cmd = 'qexe.py -w ' + workdir + ' ' + jid+' -- '+cmd
    print cmd
 
    if opt.dryrun:
        print 'Dry Run (command will not be executed)'
        continue



    subprocess.call(cmd,shell=True)

