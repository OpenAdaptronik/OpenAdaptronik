import numpy
import matplotlib.pyplot as plt
def waterfall(t,signal,nfft,overlap = .75):
    t0 = t[0]
    ts = t[1]-t0
    fs = 1/ts 
    n = len(signal)
    
    k = numpy.arange(nfft)
    T = nfft/fs
    frq = k/T
    frq = frq[range(int(nfft/2))]
    n0 = 0
    n1 = n0+nfft
    y = signal[n0:n1]
    Y = numpy.fft.fft(y)
    Y = numpy.multiply(Y,[2.0]*len(Y))
    Y = numpy.divide(Y, [nfft]*len(Y))
    Y = Y[range(int((nfft)/2))]
    wf = [numpy.absolute(Y)]
    tout = [t0]
    n0 = int(n0+(nfft*(1-overlap)))
    n1 = n0+nfft
    
    
    while n1 < n:
        y = signal[n0:n1]
        Y = numpy.fft.fft(y)
        Y = numpy.divide(Y, [nfft]*len(Y))
        Y = Y[range(int((nfft)/2))]
        wf = numpy.concatenate([wf, [numpy.absolute(Y)]],0)
        tout = tout + [t[n0]]
        n0 = int(n0+(nfft*(1-overlap)))
        n1 = n0+nfft
    return wf, frq, tout

    
def plotwaterfall(t,signal,nfft):
    wf, frq, tout = waterfall(t, signal, nfft)
    plt.pcolor(wf)
    plt.show()
    return wf, frq, tout
