class gene(object):
    """docstring for gene"""
    def __init__(self, arg):
        self.mirna=[]
        for i in arg:
            self.mirna.append(int(i[1]))

        self.degradome=[]
        for i in arg:
            self.degradome.append(int(i[0]))





################################################################



def peakdetect(y_axis, x_axis = None, lookahead = 30, delta=0):

    from numpy import NaN, Inf, arange, isscalar, array, asarray
    import numpy as np

    """
    this function is a slightly modified version of 
    'https://gist.github.com/sixtenbe/1178136' gist.
    Converted from/based on a MATLAB script at: 
    http://billauer.co.il/peakdet.html
    
    function for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively
    
    keyword arguments:
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- (optional) A x-axis whose values correspond to the y_axis list
        and is used in the return to specify the postion of the peaks. If
        omitted an index of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 200) 
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the function from picking up false peaks towards to end of
        the signal. To work well delta should be set to delta >= RMSnoise * 5.
        (default: 0)
            delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the function
    
    return -- two lists [max_peaks, min_peaks] containing the positive and
        negative peaks respectively. Each cell of the lists contains a tupple
        of: (position, peak_value) 
        to get the average peak value do: np.mean(max_peaks, 0)[1] on the
        results to unpack one of the lists into x, y coordinates do: 
        x, y = zip(*tab)
    """
    max_peaks = []
    min_peaks = []
    dump = []   #Used to pop the first hit which almost always is false

    avg=int(sum(y_axis)/float(len(y_axis)))

    # store data length for later use
    length = len(y_axis)
    x_axis=range(len(y_axis))
    
    #perform some checks
    if lookahead < 1:
        raise ValueError, "Lookahead must be '1' or above in value"
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError, "delta must be a positive number"
    
    #maxima and minima candidates are temporarily stored in
    #mx and mn respectively
    mn, mx = np.Inf, -np.Inf
    
    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x
        
        ####look for max####
        if y < mx-delta and mx != np.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if max(y_axis[index:index+lookahead]) < mx and mx>avg:
                max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue
            #else:  #slows shit down this does
            #    mx = ahead
            #    mxpos = x_axis[np.where(y_axis[index:index+lookahead]==mx)]

        ####look for min####
        if y > mn+delta and mn != -np.Inf:
            #Minima peak candidate found 
            #look ahead in signal to ensure that this is a peak and not jitter
            if min(y_axis[index:index+lookahead]) > mn:
                min_peaks.append([mnpos, mn])
                dump.append(False)
                #set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
            #else:  #slows shit down this does
            #    mn = ahead
            #    mnpos = x_axis[np.where(y_axis[index:index+lookahead]==mn)]

    try:
        if dump[0]:
            max_peaks.pop(0)
        else:
            min_peaks.pop(0)
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass

    return y_axis, max_peaks, avg


###################################################################################

def peak_width(y_axis,max_peaks,avg):
    peaks=[]
    for peak in max_peaks:
        start, stop = peak[0],peak[0]
        try:
            while y_axis[start-1]<= y_axis[start] and y_axis[start-1]>avg:
                start -= 1
        except:
            pass
        try:
            while y_axis[stop+1]<=y_axis[stop] and y_axis[stop+1]>avg:
                stop += 1
        except:
            pass
        if (stop-start)>19:
            peaks.append((start,stop))
    return peaks

###################################################################################


if __name__ == "__main__":
    dic=dict()
    peak=dict()
    fh=open('bti_germ.cov')

    for line in fh:
        tabs=line.split()
        try:
            # tabs[2] is the miRNA coverage, the tabs[3] is degradome coverage
            dic[tabs[0]].append((tabs[2],tabs[3])) 
        except:
            dic[tabs[0]]=[]
            dic[tabs[0]].append((tabs[2],tabs[3]))

    
    for key in dic.keys():
        c=gene(dic[key])
        a,b,c = peakdetect(c.mirna)
        peaks_bed = peak_width(a,b,c)
        if len(peaks_bed)>0:
            for i in peaks_bed:
                print "%s\t%i\t%i" % (key,i[0],i[1])



#     from math import pi
#     import pylab
    
#     i = 10000
#     x = np.linspace(0,3.7*pi,i)
#     y = (0.3*np.sin(x) + np.sin(1.3 * x) + 0.9 * np.sin(4.2 * x) + 0.06 * 
#     np.random.randn(i))
#     y *= -1
    
#     _max, _min = peakdetect(y, x, 750, 0.30)
#     xm = [p[0] for p in _max]
#     ym = [p[1] for p in _max]
#     xn = [p[0] for p in _min]
#     yn = [p[1] for p in _min]
    
#     plot = pylab.plot(x, y)
#     pylab.hold(True)
#     pylab.plot(xm, ym, 'r+')
#     pylab.plot(xn, yn, 'g+')
    
    
#     pylab.show()