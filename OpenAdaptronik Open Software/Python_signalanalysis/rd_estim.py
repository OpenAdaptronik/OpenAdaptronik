import numpy as np

def RD_frame(data,trigger_level,mirror,trigger_condition,RD_order,counter_in,RD_sequence):
    
    counter = counter_in
    D = RD_sequence
    trigger_event_time = []
    trigger_event_amp = []
    
    for N in range(len(data)-RD_order-3):
        
        """take sequence, walk through the frame"""
        u = data[N:N+RD_order+1]
        
        trigger = 0
        
        if trigger_condition == 1:
            if (u[1]<trigger_level and u[2]>trigger_level):
                trigger = 1.0
            elif mirror == 1 and u[1]<trigger_level and u[2]>trigger_level:
                trigger= -1.0
                
        if trigger_condition == 2:
            if u[1]>trigger_level:
                trigger = 1.0
            elif mirror == 1 and u[1]<-trigger_level:
                trigger = -1.0
                
        if trigger_condition == 3:
            if u[2]>trigger_level:
                if (u[3]<u[2] and u[2]>u[1]) or \
                (u[3]>u[2] and u[2]<u[1]):
                    trigger = 1.0
            elif mirror == 1 and u[2]<-trigger_level:
                if (u[3]<u[2] and u[2]>u[1]) or \
                (u[3]>u[2] and u[2]<u[1]):
                    trigger = -1.0
        
        if trigger != 0:
            
            #D = (1.0-1.0/counter)*D + 1.0/counter*trigger*u[0:len(D)] 
            D = D + trigger*u[0:len(D)]
            counter = counter+1.0
            trigger_event_time = np.append(trigger_event_time,float(N))
            trigger_event_amp = np.append(trigger_event_amp,trigger*trigger_level)
            
    D = D / counter
            
            
    return D, counter, trigger_event_time, trigger_event_amp
     
            
                
             
        
    
    
  