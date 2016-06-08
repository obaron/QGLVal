#!/bin/env python

import optparse

usage = 'usage: %prog'
parser = optparse.OptionParser(usage)

parser.add_option('-r', '--resdir', dest='resdir', type='string', default = '.', help='res directory')
parser.add_option('', '--store', dest='store', action='store_true', default = True, help='Store')

(opt, args) = parser.parse_args()

import time
import sys
import ROOT
import copy
import commands, os.path
import numpy as n

from toPlot import samples
import QCD_plot

settings = {}
store = []

settings.update(QCD_plot.settings)
store += QCD_plot.store

outhistos = 'output/scaled'
for d in [outhistos]:
    if os.path.exists(d): continue
    print 'Creating',d
    os.makedirs(d)

histoSig = {}
from services import Histo

for var,(title) in settings.iteritems():

    for s in samples.itervalues():
        nEntries = 0

        if(hasattr(s, "components")):
            histos = []
            notFound = []
            for c in s.components:

                filename = opt.resdir+'./res/'+c.label+ "_cat2_singleH.root"

                infile_nEvt = ROOT.TFile.Open(filename)
                nEntries = infile_nEvt.Get("h_cutFlow").GetBinContent(1)
                
                infile = ROOT.TFile.Open(filename)
                hin = infile.Get(var)
                if not isinstance(hin, ROOT.TH1):
                    notFound.append( (c.label,filename) )
                    continue
                
                htmp = hin.Clone()
                
                htmp.Sumw2()
                
                if( nEntries != 0):
                    htmp.Scale((1./nEntries) * c.sigma * 1000.* float(2.6) )

                if(htmp.Integral()!=0):
                    htmp.SetMarkerStyle(20)
                    htmp.SetMarkerSize(1.2)
                    htmp.SetLineWidth(2)
                    htmp.SetLineColor(ROOT.kBlack)
                histos.append(htmp)

            if notFound:
                raise RuntimeError('Failed to retrieve %s' % notFound)

            h1 = histos[0]
            [h1.Add(hi) for hi in histos[1:]]
        
            h = Histo.fromTH1(h1)

        outfilename = "%s/QCD_test.root" % (outhistos)
        outfile = ROOT.TFile(outfilename, "UPDATE")

        outfile.cd()
        if var in store:
            print "--> storing var", var
            h.GetHisto().SetName(var)
            h.GetHisto().Write()
        outfile.Close()
