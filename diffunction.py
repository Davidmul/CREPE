from __future__ import division
import numpy as np
import pylab as pyl
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import argparse
import os

def findinterpolpoints2low(QE):
  lowQE = (QE[0] - 2*(QE[1] - QE[0])) #locked off
  return lowQE

def findinterpolpointslow(QE):
  lowQE = (QE[0] - (QE[1] - QE[0])) #locked off
  return lowQE
  
#def findinterpolpointshigh(QE,ne):
#  highQE = (QE[ne] + (QE[ne] - QE[ne-1])) #locked off
#  highQE=highQE.clip(min=0)
#  return highQE

#def findinterpolpoints2high(QE,ne):
#  highQE = (QE[ne] + 2*(QE[ne] - QE[ne-1])) #locked off
#  highQE=highQE.clip(min=0)
#  return highQE

def findinterpolpointslowenergy(QE):
    if  QE[0]==0:
        lowQE=0
    else: 
        lowQE = (QE[0]**2)/QE[1] #checked
    return lowQE

def findinterpolpointshighenergy(QE,ne):
    #print 'finalpoint is', QE[ne]
    #print 'second final point is', QE[ne-1]
    if  QE[ne]==0 or QE[ne-1]==0:
        highQE=0
    else:
        highQE = (QE[ne]**2)/QE[ne-1] #checked
    if highQE<0:
	#highQE=0
	return highQE
    else:
	return highQE

#function finding the nearest energy value to 1.4GHz and 151 MHz.
def find_nearest(targetarray,value,bfeldy,E):
    n=0
    outputarray = []
    outputenergy = []
    #E = E/1.6e-12
    for makefreq in E:
        freqrange = 16*(((E))**2)*(bfeldy*1e6)#calculating frequency range E is ALREADY in Ge
        #freqrange = 16*(((E)/1e9)**2)*(bfeldy[n]*1e6)#calculating frequency range E is in Ge
        #print 'freqrange is', freqrange
        idx = (np.abs(freqrange-value)).argmin()
        outputarray = np.append(outputarray,targetarray[idx,n])
        outputenergy = np.append(outputenergy,E[idx])
        n=n+1
        #print n
    return outputarray, outputenergy



def plotting(a,b,c):
        elapsedtime = c*0.001
        plt.plot(a,b,'k-')
        plt.xlabel('Radius (kpc)', fontsize=20)
        plt.xlim([0.1,10.1])
	plt.ylabel('Spectral Index', fontsize=20)
	plt.ylim([-1.5,-0.48])
        plt.title(str(elapsedtime)+'Myrs elapsed', fontsize=20)
        plt.savefig('spectralevolutioniter'+str(c))
        plt.close()
        os.system('mv spectralevolutioniter'+str(c)+'.png makingplotsfortalk')

#######remaining functions are just dianostic and plotting functions######
def createfinalplot(a,b,n):
	obsr, obsy, error = np.loadtxt('iringspectraloutfileerrorbeamlargererror.txt',unpack=True, usecols=[0,1,2]) 
	plt.errorbar(obsr, obsy, yerr=error, fmt='o')
	plt.plot(a,b, 'k-')
	plt.xlabel('Radius (kpc)', fontsize=20)
	plt.xlim([0.1,15.1])
	plt.ylabel('Spectral Index', fontsize=20)
	plt.ylim([-2.5,-0.45])
        if n<1000:
            plt.savefig('finalspectrum00'+str(n))
            plt.close()
            os.system('mv finalspectrum00'+str(n)+'.png makingplotsfortalk')
        elif 1000<n<10000:
            plt.savefig('finalspectrum0'+str(n))
            plt.close()
            os.system('mv finalspectrum0'+str(n)+'.png makingplotsfortalk')
        else:
            plt.savefig('finalspectrum'+str(n))
            plt.close()
            os.system('mv finalspectrum'+str(n)+'.png makingplotsfortalk')

def createfinalplotHBALBA(a,b,n):
	obsr, obsy, error = np.loadtxt('iringspectraloutfileerrorbeamlargererror.txt',unpack=True, usecols=[0,1,2]) 
	plt.errorbar(obsr, obsy, yerr=error, fmt='o')
	plt.plot(a,b, 'k-')
	plt.xlabel('Radius (kpc)', fontsize=20)
	plt.xlim([0.1,15.1])
	plt.ylabel('Spectral Index', fontsize=20)
	plt.ylim([-2.5,-0.45])
        if n<1000:
            plt.savefig('HBALBAfinalspectrum00'+str(n))
            plt.close()
            os.system('mv HBALBAfinalspectrum00'+str(n)+'.png makingplotsfortalk')
        elif 1000<n<10000:
            plt.savefig('HBALBAfinalspectrum0'+str(n))
            plt.close()
            os.system('mv HBALBAfinalspectrum0'+str(n)+'.png makingplotsfortalk')
        else:
            plt.savefig('HBALBAfinalspectrum'+str(n))
            plt.close()
            os.system('mv HBALBAfinalspectrum'+str(n)+'.png makingplotsfortalk')



def outputfiles(radius,specindexhbavla,specindexlbahba,LOFARLBACRE,LOFARHBACRE,GMRT330ACRE,GMRT610ACRE,VLA1400CRE,VLA3000CRE,n):
	#print 'Creating output files'
	outfile = open('finaloutputhbavla'+str(n)+'.txt', 'w')
	np.savetxt(outfile, np.transpose((radius, specindexhbavla)), fmt='%s %s')
	os.system('mv finaloutputhbavla'+str(n)+'.txt makingplotsfortalk')
        outfile1 = open('finaloutputlbahba'+str(n)+'.txt', 'w')
	np.savetxt(outfile1, np.transpose((radius, specindexlbahba)), fmt='%s %s')
	os.system('mv finaloutputlbahba'+str(n)+'.txt makingplotsfortalk')
	outfile3 = open('LOFARLBACREdistribution'+str(n)+'.txt', 'w')
        np.savetxt(outfile3, np.transpose((radius,LOFARLBACRE )), fmt='%s %s')
	os.system('mv LOFARLBACREdistribution'+str(n)+'.txt makingplotsfortalk')
	outfile4 = open('LOFARHBACREdistribution'+str(n)+'.txt', 'w')
        np.savetxt(outfile4, np.transpose((radius, LOFARHBACRE)), fmt='%s %s')
        os.system('mv LOFARHBACREdistribution'+str(n)+'.txt makingplotsfortalk')
        outfile5 = open('GMRT330CREdistribution'+str(n)+'.txt', 'w')
	np.savetxt(outfile5, np.transpose((radius, GMRT330ACRE)), fmt='%s %s')
	os.system('mv GMRT330CREdistribution'+str(n)+'.txt makingplotsfortalk')
	outfile6 = open('GMRT610CREdistribution'+str(n)+'.txt', 'w')
        np.savetxt(outfile6, np.transpose((radius, GMRT610ACRE)), fmt='%s %s')
	os.system('mv GMRT610CREdistribution'+str(n)+'.txt makingplotsfortalk')
	outfile7 = open('VLA1400CREdistribution'+str(n)+'.txt', 'w')
        np.savetxt(outfile7, np.transpose((radius, VLA1400CRE)), fmt='%s %s')
        os.system('mv VLA1400CREdistribution'+str(n)+'.txt makingplotsfortalk')
        outfile8 = open('VLA3000CREdistribution'+str(n)+'.txt', 'w')
        np.savetxt(outfile8, np.transpose((radius, VLA3000CRE)), fmt='%s %s')
        os.system('mv VLA3000CREdistribution'+str(n)+'.txt makingplotsfortalk')

#for testing purposes
#def outputfiles(radius,LOFARCRE,VLACRE,n):
#	#print 'Creating output files'
#	outfile2 = open('LOFARCREdistribution'+str(n)+'.txt', 'w')
#        np.savetxt(outfile2, np.transpose((radius, LOFARCRE)), fmt='%s %s')
#	os.system('mv LOFARCREdistribution'+str(n)+'.txt makingplotsfortalk')
#	outfile3 = open('VLACREdistribution'+str(n)+'.txt', 'w')
#        np.savetxt(outfile3, np.transpose((radius, VLACRE)), fmt='%s %s')
#        os.system('mv VLACREdistribution'+str(n)+'.txt makingplotsfortalk')




def fitscalelengthforVLA(sourcex,sourcey,n):
	expfunc = lambda p,x: p[0]*np.exp(p[1]*x)
	experr = lambda p,x,y: expfunc(p,x)-y
	p0=[20000,-0.5] #Initial guesses
	fitout=leastsq(experr,p0[:],args=(sourcex,sourcey))
	paramsout=fitout[0]
	outfile = open('VLAfitsresults'+str(n)+'.txt', 'w')
        print >>outfile, 'Fitted Parameters for 1.4GHZ:\na = %.2f , b = %.2f' % (paramsout[0],paramsout[1])
        print >>outfile, 'Scale length at VLA frequency is found to be %.2f kpc' % (-1*(1/paramsout[1]))
	#print('Fitted Parameters for 1.4GHZ:\na = %.2f , b = %.2f' % (paramsout[0],paramsout[1]))
 	#print ('Scale length at VLA frequency is found to be %.2f kpc' % (-1*(1/paramsout[1])))
	xsmooth=np.linspace(sourcex[0],sourcex[-1])
	xsmooth=np.linspace(sourcex[0],sourcex[-1])
	plt.rc('font',family='serif')
	fig1=plt.figure(1)
	frame1=fig1.add_axes((.1,.3,.8,.6)) 
	plt.plot(sourcex,sourcey,'r.')
	plt.plot(xsmooth,expfunc(paramsout,xsmooth),'b-')
	plt.yscale('log')
	frame1.set_xticklabels([]) #We will plot the residuals below, so no x-ticks on this plot
	plt.title('Scale Length Fit for VLA Frequency')
	plt.ylabel('N(R)', fontsize=20)
	plt.xlabel('Radius kpc', fontsize=20)
	plt.xlim(np.min(sourcex),np.max(sourcex))
	plt.grid(True)

	from matplotlib.ticker import MaxNLocator
	plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower')) #Removes lowest ytick label
 
	frame2=fig1.add_axes((.1,.1,.8,.2))
	plt.plot(sourcex,expfunc(paramsout,sourcex)-sourcey,'k')
	plt.ylabel('Residuals')
	plt.grid(True)
	plt.savefig('VLAscalelength'+str(n)+'.png',dpi=100)
	plt.close()
	plt.xlabel('Radius kpc', fontsize = 20)
	plt.xlim(np.min(sourcex),np.max(sourcex))
	os.system('mv VLAscalelength'+str(n)+'.png makingplotsfortalk')
	
	scalelengthoutfile3 = open('VLACREscalelength'+str(n)+'.txt', 'w')
	np.savetxt(scalelengthoutfile3, np.transpose((xsmooth, expfunc(paramsout,xsmooth))), fmt='%s %s')
	os.system('mv VLACREscalelength'+str(n)+'.txt makingplotsfortalk')
	os.system('mv VLAfitsresults'+str(n)+'.txt makingplotsfortalk')


def fitscalelengthforLOFAR(sourcex,sourcey,n):
        expfunc = lambda p,x: p[0]*np.exp(p[1]*x)
        experr = lambda p,x,y: expfunc(p,x)-y
        p0=[20000,-0.5] #Initial guesses
        fitout=leastsq(experr,p0[:],args=(sourcex,sourcey))
        paramsout=fitout[0]
        outfile = open('LOFARfitsresults'+str(n)+'.txt', 'w')
        print >>outfile, 'Fitted Parameters for 151MHZ:\na = %.2f , b = %.2f' % (paramsout[0],paramsout[1])
        print >>outfile, 'Scale length at LOFAR frequency is found to be %.2f kpc' % (-1*(1/paramsout[1]))
        #print('Fitted Parameters for 151MHZ:\na = %.2f , b = %.2f' % (paramsout[0],paramsout[1]))
        #print ('Scale length is found to be %.2f kpc' % (-1*(1/paramsout[1])))
        xsmooth=np.linspace(sourcex[0],sourcex[-1])
        xsmooth=np.linspace(sourcex[0],sourcex[-1])
        plt.rc('font',family='serif')
        fig1=plt.figure(1)
        frame1=fig1.add_axes((.1,.3,.8,.6))
        plt.plot(sourcex,sourcey,'r.')
        plt.plot(xsmooth,expfunc(paramsout,xsmooth),'b-')
        plt.yscale('log')
        frame1.set_xticklabels([]) #We will plot the residuals below, so no x-ticks on this plot
        plt.title('Scale Length Fit at 150 MHz')
        plt.ylabel('N(R)', fontsize=20)
        plt.xlabel('Radius kpc',fontsize=20)
        plt.xlim(np.min(sourcex),np.max(sourcex))
        plt.grid(True)

        from matplotlib.ticker import MaxNLocator
        plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower')) #Removes lowest ytick label

        frame2=fig1.add_axes((.1,.1,.8,.2))
        plt.plot(sourcex,expfunc(paramsout,sourcex)-sourcey,'k')
        plt.ylabel('Residuals')
        plt.grid(True)
        plt.savefig('LOFARscalelength'+str(n)+'.png',dpi=100)
        plt.close()
        plt.xlabel('Radius kpc',fontsize=20)
        plt.xlim(np.min(sourcex),np.max(sourcex))
        os.system('mv LOFARscalelength'+str(n)+'.png makingplotsfortalk')

        scalelengthoutfile3 = open('LOFARCREscalelength'+str(n)+'.txt', 'w')
        np.savetxt(scalelengthoutfile3, np.transpose((xsmooth, expfunc(paramsout,xsmooth))), fmt='%s %s')
        os.system('mv LOFARCREscalelength'+str(n)+'.txt makingplotsfortalk')
        os.system('mv LOFARfitsresults'+str(n)+'.txt makingplotsfortalk')               

def timeplots(time,center,extended,inter,arm):
    plt.rc('font',family='serif')
    fig1=plt.figure(1)
    plt.plot(time,center,'b--')
    plt.plot(time,extended,'b-')
    plt.plot(time,inter,'r--')
    plt.plot(time,arm,'r-')
    plt.title('Evolution of Spectral Index')
    plt.ylabel('Spectral Index', fontsize=20)
    plt.xlabel('Time (Myrs)',fontsize=20)
    plt.xlim(np.min(time),np.max(time))
    plt.savefig('spectralindexevolution.png',dpi=100)
    plt.close()
    os.system('mv spectralindexevolution.png makingplotsfortalk')
    outfile = open('centerplot.txt', 'w')
    np.savetxt(outfile, np.transpose((time, center)), fmt='%s %s')
    os.system('mv centerplot.txt makingplotsfortalk')
    outfile2 = open('extendedplot.txt', 'w')
    np.savetxt(outfile2, np.transpose((time, extended)), fmt='%s %s')
    os.system('mv extendedplot.txt makingplotsfortalk')
    outfile3 = open('interarmplot.txt', 'w')
    np.savetxt(outfile3, np.transpose((time, inter)), fmt='%s %s')
    os.system('mv interarmplot.txt makingplotsfortalk')
    outfile4 = open('extendedplotD.txt', 'w')
    np.savetxt(outfile4, np.transpose((time, arm)), fmt='%s %s')
    os.system('mv extendedplotD.txt makingplotsfortalk')

def plotenergyevo(finale, energy,n):
	plt.plot(finale,energy,'r')
	#plt.yscale('log')
	plt.xscale('log')
        plt.ylabel('N(E)', fontsize=20)
        plt.xlabel('Energy GeV',fontsize=20)
        plt.xlim(np.min(finale),np.max(finale))
	plt.ylim(np.min(energy),np.max(energy))
	plt.savefig('ENERGYPLOT'+str(n)+'.png',dpi=100)
	plt.close()
	outfile = open('ENERGYresults'+str(n)+'.txt', 'w')
	np.savetxt(outfile, np.transpose((finale, energy)), fmt='%s %s')
	os.system('mv ENERGYresults'+str(n)+'.txt makingplotsfortalk')
	os.system('mv ENERGYPLOT'+str(n)+'.png makingplotsfortalk')
