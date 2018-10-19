'''
Created on 16.11.2017

@author: stoll
'''
import numpy
from numpy import NaN as NaN
import detect_peaks

def findFreeSpace(FFTdata, Fdata, Finterval, level):
    
    if Finterval[0] < Fdata[0]:
        Finterval[0] = Fdata[0]
    if Finterval[1] > Fdata[-1]:
        Finterval[1] = Fdata[-1]
    
    start_ind = numpy.where(Fdata <= Finterval[0])[0][-1]
    end_ind = numpy.where(Fdata >= Finterval[1])[0][0]
    
    
    problem_ind = numpy.where(FFTdata[numpy.arange(start_ind,end_ind)]>level)+start_ind
    
    problem_ind = numpy.insert(problem_ind, 0, start_ind)
    problem_ind = numpy.append(problem_ind, end_ind)
    
    maxspace = 0
    maxspacepos = 0
    
    if (len(problem_ind) < 1):
        return Finterval

    for i in range(len(problem_ind)-1):
        
        space = Fdata[problem_ind[i+1]]-Fdata[problem_ind[i]-1]
        if space > maxspace:
            maxspace = space
            maxspacepos = problem_ind[i]
            
    if maxspace == 0:
        freespace = [0, 0]
    else:
        freespace = [Fdata[maxspacepos], Fdata[maxspacepos]+maxspace]
    
    return freespace


def findRealPeaks(A_1, tf, f, f_valid, sens_FindPeaks, sens_TFOrigin, sens_SumUpPeaks, maxamp, problvl):
    
    nt = A_1.shape[0]
    nf = A_1.shape[1]
    
    factor = 1.5
    
    peakamp = [0]*nt
    peakloc = [0]*nt
    collectivepeaks = [[]]*nt
    ispeakfromTF = [[]]*nt
    gap = [[]]*nt
    numrealpeaks = [0]*nt
    peakloc_ind = [False]*nf
    
    
    def TFpeakCheck(peakAmp, TFAmp):
        return ((peakAmp/(factor*(sens_TFOrigin/50)*TFAmp)) < problvl) & (TFAmp > -2*(sens_TFOrigin/50)+4)
    
    for i in range(nt):
        peakloci = detect_peaks.detect_peaks(A_1[i][f_valid],mph=0.15*numpy.square((sens_FindPeaks/50)-2)*maxamp,mpd=3)
        peakampi = A_1[i][f_valid][peakloci]
        peakloc[i] = f[f_valid][peakloci]
        peakamp[i] = peakampi
        collectivepeaks[i] = []
        
        if len(peakloci)==0:
            ispeakfromTF[i] = [NaN]
            gap[i] = []
        elif len(peakloci)==1:
            numrealpeaks[i] = 1;
            collectivepeaks[i] = [[peakloci[0]]];
            ispeakfromTF[i] = [[((peakampi[0]/(factor*(sens_TFOrigin/50)*tf[peakloci[0]])) < problvl) & (tf[peakloci[0]] > -2*(sens_TFOrigin/50)+4)]]
            peakloc_ind[peakloci[0]] = True;
            gap[i] = []
        else:
            numrealpeaks[i] = 1
            collectivepeaks[i] = [[peakloci[0]]]
            ispeakfromTF[i] = [[((peakampi[0]/(factor*(sens_TFOrigin/50)*tf[peakloci[0]])) < problvl) & (tf[peakloci[0]] > -2*(sens_TFOrigin/50)+4)]]
            gap[i] = [];
            for j in range(len(peakloci)-1):
                freesp = findFreeSpace(A_1[i], f, f[[peakloci[j], peakloci[j+1]]], problvl);
                gap[i].append(freesp)
                if (f[peakloci[j+1]] - f[collectivepeaks[i][-1][-1]] > -10*(sens_SumUpPeaks/50+20)) & (gap[i][j][0] != 0):
                    numrealpeaks[i] = numrealpeaks[i] + 1
                    collectivepeaks[i].append([peakloci[j+1]])
                    ispeakfromTF[i].append([((peakampi[0]/(factor*(sens_TFOrigin/50)*tf[peakloci[0]])) < problvl) & (tf[peakloci[0]] > -2*(sens_TFOrigin/50)+4)])
                    
                else:
                    collectivepeaks[i][-1].append(peakloci[j+1])
                    ispeakfromTF[i][-1].append(((peakampi[0]/(factor*(sens_TFOrigin/50)*tf[peakloci[0]])) < problvl) & (tf[peakloci[0]] > -2*(sens_TFOrigin/50)+4))
                    
                peakloc_ind[peakloci[j]] = True
        
    realpeakloc = []
    isrealpeakfromTF = [[]]*nt
    for i in range(nt):
        realpeakloc.append([0]*numrealpeaks[i])
        if numrealpeaks[i]>0:
            for j in range(numrealpeaks[i]):
                realpeakloc[i][j]=numpy.mean(f[f_valid][collectivepeaks[i][j]])
                if len(ispeakfromTF[i])>0:
                    if numpy.mean(ispeakfromTF[i][j])<= .5:
                        isrealpeakfromTF[i].append(0)
                        pass
                    else:
                        isrealpeakfromTF[i].append(1)
                        pass
                    pass
                else:
                    isrealpeakfromTF[i].append(.95);
                    pass
                pass
            pass
        else:
            isrealpeakfromTF[i] = [.95]
            pass
            
        pass
    meannumrealpeaks = numpy.mean(numpy.array(numrealpeaks)[numpy.where(numpy.array(numrealpeaks)>0)])
        
    return peakloc_ind, collectivepeaks, numrealpeaks, meannumrealpeaks, realpeakloc, isrealpeakfromTF

def analyseFDAmplitude(A_1, f, iValidF, realpeakloc, numrealpeaks, multiplepeaks, sens_UnproblematicLFD, sens_UnproblematicHFD, sens_UnproblematicNFD, problvl):
    
    nt = A_1.shape[0];
    
    minfreespaceLFD = -50*(sens_UnproblematicLFD/50)+105; # min frequency range that has to be free so that the low frequency domain is considered unproblematic
    minfreespaceHFD = -50*(sens_UnproblematicHFD/50)+105; # min frequency range that has to be free so that the high frequency domain is considered unproblematic
    minfreespaceNFD = -30*(sens_UnproblematicNFD/50)+31;  # min frequency range (to each side of the respective peak) that has to be free so that the peak neighboring frequency domain is considered unproblematic
    
    LFDsize = 120;  # determines the size of the examined low frequency domain if there is no peak
    HFDsize = 120;  # determines the size of the examined high frequency domain if there is no peak
    NFDsize = 80;   # determines how far the examined peak neighboring frequency domain reaches to each side of the peak in Hz

    firstpeakpos = []
    lastpeakpos = []
    freespace_lfd = []
    freespace_hfd = []    
    islfdunproblematic = []
    ishfdunproblematic = []
    freespace_r = []
    freespace_l = []
    isnfdunproblematic = []
    
    for i in range(nt):
        if numrealpeaks[i]>0:
            firstpeakpos.append(realpeakloc[i][0])
            lastpeakpos.append(realpeakloc[i][-1])
        else:
            firstpeakpos.append(LFDsize)
            lastpeakpos.append(f[-1]-HFDsize)
        
       
        freespace_lfd.append(findFreeSpace(A_1[i][iValidF], f[iValidF], [f[iValidF][0], firstpeakpos[i]], problvl))
        if (freespace_lfd[i][1]-freespace_lfd[i][0]) >= minfreespaceLFD:
            islfdunproblematic.append(True)
        else:
            islfdunproblematic.append(False)
        
        freespace_hfd.append(findFreeSpace(A_1[i][iValidF], f[iValidF], [lastpeakpos[i], f[iValidF][-1]], problvl))
        if (freespace_hfd[i][1]-freespace_hfd[i][0]) >= minfreespaceHFD:
            ishfdunproblematic.append(True)
        else:
            ishfdunproblematic.append(False)
        
        freespace_r.append([])
        freespace_l.append([])
        isnfdunproblematic.append(True)
        if ~multiplepeaks:
            if numrealpeaks[i]==1:
                freespace_r[i] = findFreeSpace(A_1[i][iValidF], f[iValidF], [realpeakloc[i][0], realpeakloc[i][0]+NFDsize], problvl)
                freespace_l[i] = findFreeSpace(A_1[i][iValidF], f[iValidF], [realpeakloc[i][0]-NFDsize, realpeakloc[i][0]], problvl);
                if (freespace_r[i][1]-freespace_r[i][0] < minfreespaceNFD) | (freespace_l[i][1]-freespace_l[i][0] < minfreespaceNFD):
                    isnfdunproblematic[i] = False
            elif numrealpeaks[i] > 1:
                isnfdunproblematic[i] = False
    return islfdunproblematic, ishfdunproblematic, isnfdunproblematic


def analyseTimeVariance(A_1, f_step, collectivepeaks, sens_TimevariantBehavior, peakloc_ind, problvl):
    
    nt = A_1.shape[0];
    nf = A_1.shape[1];
    
    
    lowampot_duration = round(-20*(sens_TimevariantBehavior/50)+40)
    ##
    lowampmap = A_1<problvl
    
    lowampot_all = [[False]*nf]*(nt-lowampot_duration)
    
    
    for i in range(nt-lowampot_duration):
        lowampot_all[i][:]=numpy.all(lowampmap[i:i+lowampot_duration][:], axis=0)
        
    lowampot = numpy.any(lowampot_all, axis=0)
    
    timevariantpeaks = peakloc_ind & lowampot
    #timeinvariantpeaks = peakloc_ind & ~lowampot
        
    ispeaktimevariant = [[]]*nt;
    isrealpeaktimevariant= [[]]*nt;
    
    
    for i in range(nt):
        npeaks = len(collectivepeaks[i])
        if npeaks > 0:
            ispeaktimevariant[i] = [[]]*npeaks
            isrealpeaktimevariant[i] = [0]*npeaks
            for j in range(npeaks):
                ispeaktimevariant[i][j] = timevariantpeaks[collectivepeaks[i][j]]
                isrealpeaktimevariant[i][j] = numpy.round(numpy.mean(ispeaktimevariant[i][j]))
        else:
            ispeaktimevariant[i] = []
            isrealpeaktimevariant[i] = []
        
    return isrealpeaktimevariant
    pass



