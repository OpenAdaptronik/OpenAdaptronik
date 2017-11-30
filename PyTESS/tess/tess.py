'''
Created on 10.11.2017

@author: stoll
'''

import csv

from numpy import NaN, argwhere, transpose
import numpy
from scipy import disp
import scipy.io
import waterfalls
import detect_peaks
import helpers

import matplotlib.pyplot as plt


if __name__ == '__main__':
    doplots = False # for debugging
    measurementdataFILENAME = 'C:\\Users\\stoll\\Desktop\\openadaptronik\\TESS\\Data\\measurement_data.mat'
    tfFILENAME = 'TF.mat'
    analysisweightsFILENAME = 'C:\\Users\\stoll\\Desktop\\openadaptronik\\TESS\\Data\\analysisweights.csv'
    DesiredAmpLevel = 8
    MinFrequency = 3;               # general minimum frequency taken into account in Hz
    TFMode = 0
    #analysis sensitivity parameters
    sens_FindPeaks = 50;            # sensitivity of finding peaks in measurement data in %
                                    # | 0% -> peaks are virtually not found 
                                    # | 100% -> peaks are found very often
    
    sens_SumUpPeaks = 50;           # sensitivity of summing up peaks to real peaks in %
                                    # | 0% -> peaks are virtually never summed up to a real peak 
                                    # | 100% -> peaks are very likely to be summed up to a real peak 
    
    sens_TFOrigin = 50;             # sensitivity of finding that a real peak's origin is the TF (and not the excitation) in %
                                    # | 0% -> finding that the origin of a real peak is the excitation is very likely 
                                    # | 100% -> finding that the origin of a real peak is the TF is very likely
    
    sens_TimevariantBehavior = 50;  # sensitivity of considering the vibration behavior to be timevariant
                                    # | 0% -> behavior is virtually never considered timevariant 
                                    # | 100% -> behavior is very easily considered timevariant
                                    
    sens_UnproblematicLFD = 50;     # sensitivity of considering the amplitude level in the low frequency domain (LFD) of A_1 to be unproblematic
                                    # | 0% -> the amp level in the LFD is virtually never considered unproblematic
                                    # | 100% -> the amp level in the LFD is very easily considered unproblematic
                                    
    sens_UnproblematicNFD = 50;     # sensitivity of considering the amplitude level in the (peak-) neighboring frequency domain (NFD) of A_1 to be unproblematic
                                    # | 0% -> the amp level in the NFD is virtually never considered unproblematic
                                    # | 100% -> the amp level in the NFD is very easily considered unproblematic
                                    
    sens_UnproblematicHFD = 50;     # sensitivity of considering the amplitude level in the high frequency domain (HFD) of A_0 to be unproblematic
                                    # | 0% -> the amp level in the HFD is virtually never considered unproblematic
                                    # | 100% -> the amp level in the HFD is very easily considered unproblematic
    SV = numpy.zeros(17)
    
    # reading data from file ... put database request here:
    measurementdata = scipy.io.loadmat(measurementdataFILENAME);
    
    # solution identifiers
    SID = []
    SID.append('System Verstimmen, passiv')
    SID.append('System Verstimmen, schaltbar')
    SID.append('System Verstimmen, semi-aktiv')
    SID.append('Daempfung erhoehen, passiv')
    SID.append('Daempfung erhoehen, schaltbar')
    SID.append('Daempfung erhoehen, semi-aktiv')
    SID.append('Aktorik in Struktur einbringen')
    SID.append('Tilger, passiv')
    SID.append('Tilger, schaltbar')
    SID.append('Neutralisator, passiv')
    SID.append('Neutralisatior, schaltbar')
    SID.append('Adaptiver Tilger/Neutralisator')
    SID.append('Elastische Lagerung, passiv')
    SID.append('Elastische Lagerung, schaltbar')
    SID.append('Elastische Lagerung, semi-aktiv')
    SID.append('Elastische Lagerung, aktiv')
    SID.append('Inertialmassenaktor IMA')
    
    
    # CSV-Datei mit den Bewertungsgewichten/Parametern einlesen
    analysisweights = []
    with open(analysisweightsFILENAME) as csvfile:
        rd = csv.reader(csvfile)
        for row in rd:
            analysisweights.append(row)
    del analysisweights [0]
    analysisweights = transpose([[float(j) for j in i] for i in analysisweights])
    
    # hier noch checken, ob gleich formatiert, abgetastet, ...
    t = measurementdata['t'][0]
    a_0 = measurementdata['a_0'][0]
    a_1 = measurementdata['a_1'][0]
    A_0, f, tout = waterfalls.waterfall(t,a_0,512)
    A_1, f, tout = waterfalls.waterfall(t,a_1,512)
    tfmat = numpy.divide(A_1,A_0)
    tf = numpy.average(tfmat, 0, A_0)
    
    iValidF = numpy.array(numpy.where(f >= MinFrequency))
    
    ifmin = argwhere(f>=MinFrequency)[0];
    
    #f_fft_rec_max_mod = f[iValidF]
    fft_rec_max_mod = A_1[:, iValidF].flatten()
    
    top10 = [0]*10
    for i in range(10):
        imax = numpy.nanargmax(fft_rec_max_mod, axis=0)
        mmax = numpy.nanmax(fft_rec_max_mod, axis=0)
        top10[i] = fft_rec_max_mod[imax]
        fft_rec_max_mod[imax]=NaN
    maxamp = numpy.median(top10)
    problvl = (DesiredAmpLevel/100)*maxamp
    
    peakloc_ind = [False]*len(f)
    factor = 1.5
    realpeakloc = []
    numrealpeaks = numpy.array([0]*len(tout))
    isrealpeakfromTF = []
    peakamp = [0]*len(tout)
    peakloc = [0]*len(tout)
    gap = [[]]*len(tout)
    collectivepeaks = [[]]*len(tout)
    ispeakfromTF = [[]]*len(tout)
    
    for i in range(len(tout)):
        peakloci = detect_peaks.detect_peaks(A_1[i],mph=0.15*numpy.square((sens_FindPeaks/50)-2)*maxamp,mpd=3)
        peakampi = A_1[i][peakloci]
        peakloc[i] = f[peakloci]
        peakamp[i] = peakampi
        collectivepeaks[i] = []
        
        if len(peakloci)==0:
            #collectivepeaks[i].append(0)
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
            for j in range(len(peakloci)-1):
                freesp = helpers.findfreespace(A_1[i], f, f[[peakloci[j], peakloci[j+1]]], problvl);
                gap[i].append(freesp)
                if (f[peakloci[j+1]] - collectivepeaks[i][-1] > -10*(sens_SumUpPeaks/50+20)) & (gap[i][j][0] != 0):
                    numrealpeaks[i] = numrealpeaks[i] + 1
                    collectivepeaks[i].append([peakloci[j+1]])
                    ispeakfromTF[i].append([((peakampi[0]/(factor*(sens_TFOrigin/50)*tf[peakloci[0]])) < problvl) & (tf[peakloci[0]] > -2*(sens_TFOrigin/50)+4)])
                    pass
                else:
                    collectivepeaks[i][-1].append(peakloci[j+1])
                    ispeakfromTF[i][-1].append(((peakampi[0]/(factor*(sens_TFOrigin/50)*tf[peakloci[0]])) < problvl) & (tf[peakloci[0]] > -2*(sens_TFOrigin/50)+4))
                    
                    pass
                peakloc_ind[peakloci[j]] = True;
            pass
        

    isrealpeakfromTF = [[]]*len(tout)
    for i in range(len(tout)):
        realpeakloc.append([0]*numrealpeaks[i])
        if numrealpeaks[i]>0:
            for j in range(numrealpeaks[i]):
                realpeakloc[i][j]=numpy.mean(f[collectivepeaks[i][j]])
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
    
    meannumrealpeaks = numpy.mean(numrealpeaks[numpy.where(numrealpeaks>0)])

    if meannumrealpeaks > 1.1:
        multiplepeaks = True
        SV = SV+analysisweights[0]
    else:
        multiplepeaks = False
        SV = SV+analysisweights[1]
        pass
    meanpeaksfromtf = numpy.zeros(len(tout))
    
    
    for i in range(len(tout)):
        meanpeaksfromtf[i] = numpy.nanmean(isrealpeakfromTF[i])
        pass
    
    if numpy.nanmean(meanpeaksfromtf) >= 0.95:
        allpeaksfromtf = True
        SV = SV + analysisweights[2]
    else:
        allpeaksfromtf = False
        SV = SV + analysisweights[3]
        
        
    lowampot = [False]*len(f) 
    
    lowampot_duration = round(-20*(sens_TimevariantBehavior/50)+40)
    ##
    lowampmap = A_1<problvl
    for pl in range(len(f)):
        tvpeak = False
        low_duration = 0
        for tl in range(len(tout)):
            if lowampmap[tl][pl]:
                low_duration += 1
                if low_duration*(tout[1]-tout[0])>lowampot_duration:
                    tvpeak = True
                    pass
                pass
            pass
        lowampot[pl] = tvpeak
        pass
    
    timevariantpeaks = numpy.array(peakloc_ind) & numpy.array(lowampot)
    timeinvariantpeaks = numpy.array(peakloc_ind) & ~numpy.array(lowampot)
    
    ispeaktimevariant=[[[]]]
    isrealpeaktimevariant=[[]]
    
    for i in range(len(numrealpeaks)):
        ispeaktimevariant.append([])
        isrealpeaktimevariant.append([])
        if numrealpeaks[i]>0:
            for j in range(numrealpeaks[i]):
                ispeaktimevariant[i].append(collectivepeaks[i][j])
                isrealpeaktimevariant[i].append(numpy.round(numpy.mean(collectivepeaks[i][j])))
    
    meanpeakstimevariant=[]           
    for i in range(len(isrealpeaktimevariant)):
        meanpeakstimevariant.append(numpy.nanmean(isrealpeaktimevariant[i]))
        
    if numpy.nanmean(meanpeakstimevariant) >= 0.05:
        timevariantbehavior = True
        SV = SV + analysisweights[4]
    else:
        timevariantbehavior = False
        SV = SV + analysisweights[5]
        
        
    # Analysis of Amplitude Level

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
    
    for i in range(len(tout)):
        if numrealpeaks[i]>0:
            firstpeakpos.append(realpeakloc[i][0])
            lastpeakpos.append(realpeakloc[i][-1])
        else:
            firstpeakpos.append(LFDsize)
            lastpeakpos.append(f[-1]-HFDsize)
        
       
        freespace_lfd.append(helpers.findfreespace(A_1[i][iValidF], f[iValidF], [MinFrequency, firstpeakpos[i]], problvl))
        if (freespace_lfd[i][1]-freespace_lfd[i][0]) >= minfreespaceLFD:
            islfdunproblematic.append(True)
        else:
            islfdunproblematic.append(False)
        
        freespace_hfd.append(helpers.findfreespace(A_1[i][iValidF], f[iValidF], [lastpeakpos[i], f[-1]], problvl))
        if (freespace_hfd[i][1]-freespace_hfd[i][0]) >= minfreespaceHFD:
            ishfdunproblematic.append(True)
        else:
            ishfdunproblematic.append(False)
        
        freespace_r.append([])
        freespace_l.append([])
        isnfdunproblematic.append(True)
        if ~multiplepeaks:
            if numrealpeaks[i]==1:
                freespace_r[i] = helpers.findfreespace(A_1[i][iValidF], f[iValidF], [realpeakloc[i][0], realpeakloc[i][0]+NFDsize], problvl)
                freespace_l[i] = helpers.findfreespace(A_1[i][iValidF], f[iValidF], [realpeakloc[i][0]-NFDsize, realpeakloc[i][0]], problvl);
                if (freespace_r[i][1]-freespace_r[i][0] < minfreespaceNFD) | (freespace_l[i][1]-freespace_l[i][0] < minfreespaceNFD):
                    isnfdunproblematic[i] = False
            elif numrealpeaks[i] > 1:
                isnfdunproblematic[i] = False
    
    if numpy.mean(islfdunproblematic) >= .95:
        lfdunproblematic = True
        SV = SV + analysisweights[6]
    else:
        lfdunproblematic = False
        SV = SV + analysisweights[7]
        
    if ~multiplepeaks:
        if numpy.mean(isnfdunproblematic)>=.95:
            nfdunproblematic = True
            SV = SV + analysisweights[8]
        else:
            nfdunproblematic = False
            SV = SV + analysisweights[9]
    
    if numpy.mean(ishfdunproblematic) >= .95:
        hfdunproblematic = True
        SV = SV + analysisweights[10]
    else:
        hfdunproblematic = False
        SV = SV + analysisweights[11]
    
    srt = numpy.argsort(SV)
    lst = numpy.sort(SV)
    
            
    disp(' ')
    disp(' ')
    disp('     ----------------- Ergebnisse der Analyse -----------------')
    disp(['            beste Strategie: ', SID[srt[-1]]])
    disp(['       zweitbeste Strategie: ', SID[srt[-2]]])
    disp(['       drittbeste Strategie: ', SID[srt[-3]]])
    disp(['       viertbeste Strategie: ', SID[srt[-4]]])
    disp(['      fuenftbeste Strategie: ', SID[srt[-5]]])
    disp('     ----------------------------------------------------------')
    disp(' ')
    disp(' ')
        
    if doplots:
        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(t,a_1)
        plt.title('a_1')
        plt.xlabel('Time (s)')
        
        plt.subplot(2,1,2)
        plt.plot(t,a_0)
        plt.title('a_0')
        plt.xlabel('Time (s)')
        
        plt.figure()
        plt.plot(f,tf)
        plt.title('TF')
        plt.xlabel('Frequency (Hz)')
        
        plt.figure()
        plt.pcolor(f, tout, A_1)
        plt.title('A_1')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Time (s)')
        
        plt.figure()
        plt.pcolor(f, tout, A_0)
        plt.title('A_0')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Time (s)')

        plt.show(block = True)
    
    pass
    #disp(timevariantpeaks)
    #disp(timeinvariantpeaks)
    #disp(isrealpeakfromTF)    
        
    #disp(peakloc_ind)
    #disp(numrealpeaks)
    #disp(peakamp)
    #disp(peakloc)
    #disp(gap)
    #disp(collectivepeaks)
    #disp(ispeakfromTF)
    #disp(realpeakloc)
    #disp(SV)
    pass    