from ROOT import TH1D, TRandom1, TCanvas, TH1F, TMinuit
import ROOT as r

import numpy as np
import os
import sys
import re
import glob
#import pydot
import h5py
from scipy.optimize import *
import math
import pylab
from sklearn.metrics import roc_curve, roc_auc_score, auc
from skopt import *
import matplotlib.pyplot as plt
#%matplotlib inline


def setuphisto(qs, minnumber, check=0):
    '''
    This is a function to return a variable-bin histogram so that every bin has the same number of events.
    
    qs = a vector of doubles which you want to use to setup the binning.  I usually use the *background* events to setup the binning (since s/b is more sensitive to stat fluctuations when b is small than when s is small).
    minnumber = number of events in each bin.
    check: if 1, the output histogram is filled with qs and some debugging information is printed.

    The output is a vector of two 1D histograms with identical binning: one for signal and one for background.

    '''

    #out = r.Vector(out)

    qs = np.sort(qs)

    length = qs.size
    qs2 = np.empty((length))

    for i in xrange(0,length):
        qs2[i]=qs[i]


    #Now, determine the bin positions.
    binsnumber=(length-length%minnumber)/minnumber;
    
    #bins = np.empty(binsnumber+1)
    bins = np.array([])
    myindex = 0
    #print qs2[0]
    bins = np.concatenate((bins,[qs2[0]]))
    #bins = np.array(bins)
    #print bins.size
    for i in xrange(minnumber, length, minnumber):

        #std::cout << myindex << " a " << bins[myindex] << " " << qs2[i] << std::endl;

        if bins[myindex] == qs2[i]:
            myindex+=1
            
        else:
            bins = np.concatenate((bins,[qs2[i]]))
            myindex+=1
    
    
    bins = np.concatenate((bins, [qs2[length-1]]))
    #print bins[0]
    signal = TH1F("","",bins.size-1,bins)
    background = TH1F("","",bins.size-1,bins)

    #out.push_back(signal);
    #out.push_back(background);

    return signal, background



#Calculate AUC from likelihood 
#def f_auc_sd(npar, deriv, f, exponents, flag):

#loading files with tau values for signal and background
fnsig="H2bb_sd.txt"
fnbkg="g2bb_sd.txt"
j=0
fsig=pylab.loadtxt(fnsig)
fbkg=pylab.loadtxt(fnbkg)
        
#values of exponents set at function call from minimization routine

exponents = [-1,0,0,-1,1]

a = exponents[0]
b = exponents[1]
c = exponents[2]
d = exponents[3]
e = exponents[4]
        
if a==0 and b==0 and c==0 and d==0 and e==0:
    print a,b,c,d,e, 0.5
    #return 0.5
        
gs=0
gb=0
        
bkg_obs = np.zeros((100000))
        
#loop over 100,000 sig and bkg events, measure product observable on them
for j in xrange(0, 100000):
    T11H = fsig[j][0]  
    T12H = fsig[j][1]
    T13H = fsig[j][2]
    T22H = fsig[j][3] 
    T23H = fsig[j][4]
            
    T11g = fbkg[j][0]  
    T12g = fbkg[j][1]
    T13g = fbkg[j][2]
    T22g = fbkg[j][3] 
    T23g = fbkg[j][4]
                
    ts=(T11H**a)*(T12H**b)*(T13H**c)*(T22H**d)*(T23H**e)
    tb=(T11g**a)*(T12g**b)*(T13g**c)*(T22g**d)*(T23g**e)
    bkg_obs[j] = tb
    if(ts>gs):
        gs=ts
                    
    if(tb>gb):
        gb=tb
                    
#find if the product observable measured on sig or bkg obtained a higher max val 
        
max=0
        
if(gs>gb):
    max=gs
else:
    max=gb
            
max = max+(1.0/10000.0)
        
#set upper limit of range of histograms to max(=max+epsilon)

#ridiculous number of bins
#hs = TH1D("hs",("Convolved"), 1000, 0., max)
#hb = TH1D("hb",("Convolved"), 1000, 0., max)
hs, hb = setuphisto(bkg_obs, 10)
hs.SetLineColor(2)
hb.SetLineColor(9)
        
#(inefficiently) loop over events again to fill histograms for product observable measured on sig and bkg
for j in xrange(0, 100000):
    T11H = fsig[j][0]  
    T12H = fsig[j][1]
    T13H = fsig[j][2]
    T22H = fsig[j][3] 
    T23H = fsig[j][4]
            
    T11g = fbkg[j][0]  
    T12g = fbkg[j][1]
    T13g = fbkg[j][2]
    T22g = fbkg[j][3] 
    T23g = fbkg[j][4]
                
    ts=(T11H**a)*(T12H**b)*(T13H**c)*(T22H**d)*(T23H**e)
    tb=(T11g**a)*(T12g**b)*(T13g**c)*(T22g**d)*(T23g**e)
    hs.Fill(ts)
    hb.Fill(tb)
                

lr = np.empty((1001,2))
xb = 0
xs = 0
        
#calculate likelihood ratio
for k in xrange(0, 1001):
    xb=hb.GetBinContent(k)
    xs=hs.GetBinContent(k)
            
    if(xs==0):
        lr[k][0]=k
        lr[k][1]=0.0
                    
    else:
        lr[k][0]=k
        lr[k][1]=(xb/(float)(xs+xb))
                    
#sort likelihood ratio in increasing order
for p in xrange(0, 1001):
    for q in xrange(0, 1001-1):
        if(lr[q][1]>lr[q+1][1]):
            bint=lr[q+1][0]
            lt=lr[q+1][1]
            lr[q+1][0]=lr[q][0]
            lr[q+1][1]=lr[q][1]
            lr[q][0]=bint
            lr[q][1]=lt
                            
l = TH1D("l",("Likelihood"), 1001, 0., 1)
ls = TH1D("ls",("ls"), 1001, 0., 1)
lb = TH1D("lb",("lb"), 1001, 0., 1)      
        
        
#get signal and background likelihoods
for k in xrange(0, 1001):
    
    l.SetBinContent(k,lr[k][1])
    bin = (int)(lr[k][0])
    ls.SetBinContent(k, hs.GetBinContent(bin))
    lb.SetBinContent(k, hb.GetBinContent(bin))
                
ls.Scale(1.0/ls.Integral())
lb.Scale(1.0/lb.Integral())

sums=0
sumb=0
        
x = np.empty((1001))
y = np.empty((1001))
        
cut=0
        
#calculate fpr and tpr
for k in xrange(0, 1001):
    sums+=ls.GetBinContent(k)
    sumb+=lb.GetBinContent(k)
    x[k]=sumb
    y[k]=sums
            
        #calculate AUC
area = auc(x,y)
        
print a,b,c,d,e,area
        
hs.Draw()
hb.Draw("same")
        
        #hs.Delete()
        #hb.Delete()
        #ls.Delete()
        #l.Delete()
        #lb.Delete()
    
        #return for minimizaton routine
        #return (1.0-area)
