import ROOT
import collections
import os, commands
#from SingleTopWJets import *                                                                                                                                                                      
#from ttDM import *                                                                                                                                                                                
#from TT import *
#from DYJetsToLL import *                                                                                                                                                                          
#from ZJets import *
#from WJets import *
#from QCD import *
from QCD_HT import *
#from Data import *
#from otherBkg import *                                                                                                                                                              
#from BprimeBToHB1800 import*                                                                                                                                                                      
#from BprimeBToHB1800 import*                                                                                                                                                                      
samples = collections.OrderedDict()

#samples["QCD"] = QCD
samples["QCD_HT"]=QCD_HT   
