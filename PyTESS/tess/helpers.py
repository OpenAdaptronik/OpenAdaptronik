'''
Created on 16.11.2017

@author: stoll
'''
import numpy

def findfreespace(FFTdata, Fdata, Finterval, level):
    start_ind = numpy.where(Fdata >= Finterval[0])[0][0]
    end_ind = numpy.where(Fdata <= Finterval[1])[0][-1]
    problem_ind = numpy.where(FFTdata[numpy.arange(start_ind,end_ind)]>level)+start_ind
    
    maxspace = 0
    maxspacepos = 0
    
    if (len(problem_ind[0]) < 1):
        return Finterval

    for i in range(len(problem_ind[0])-1):
        
        space = problem_ind[0][i+1]-problem_ind[0][i]-2
        if space > maxspace:
            maxspace = space;
            maxspacepos = problem_ind[0][i];
            
    if maxspace == 0:
        freespace = [0, 0]
    else:
        freespace = [Fdata[maxspacepos], Fdata[maxspacepos]+maxspace]
    
    return freespace