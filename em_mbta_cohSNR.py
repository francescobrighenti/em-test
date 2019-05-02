""" 
 Written by Francesco Brighenti (francesco.brighenti.phys@gmail.com) and Giuseppe Greco, 
 based on Harry-Fairhurst arXiv:1012.4939. Produces a file txt with right ascension, 
 declination, coherent SNR and GPS time of the maximum for galaxies from the GLADE cathalog 
 within 90% C.L. of Bayestar skymap for a given Binary Neutron Star candidate trigger on GraceDB. 
 The list is in descending order of cohSNR, so the higher a galaxy is in the list the more likely
 it is that she is the host of the BNS event.

"""

import .virgotools as vt
import .em_mbta_v05 as em

delta_t=1./4096 # time resolution

galFile="BAYESgalaxies.txt"
""" makes list of galaxies objects from file with their coordinates """
galaxies=get_galaxies(galFile)

with open("GCNheader.txt") as GCNheader:
    data=GCNheader.read()
    TriggerTime=

tStart=TriggerTime-0.5 #### per esempio...
duration=

""" prepares the lists with GW Events for each ifo for the specified
    start-time and  duration""" 
eventsH=get_gw_events(tStart,duration,"H")
eventsL=get_gw_events(tStart,duration,"L")
eventsV=get_gw_events(tStart,duration,"V")


""" re-organizes the events to make lists of triple detection for the same template"""
triple_coincidences=[[evH,evL,evV] for evH in eventsH for evL in eventsL for evV in eventsV 
                            if (evH.parameters.tmplt_index==evL.parameters.tmplt_index)
                            &  (evH.parameters.tmplt_index==evV.parameters.tmplt_index)]



""" Computes the ranking """
galaxies_ranking=rank_galaxies(galaxies,triple_coincidences,TriggerTime)


""" Prints to file galaxies_ranking.txt """
np.savetxt('galaxies_ranking_'{time}'.txt'.format(time=TriggerTime),galaxies_ranking,delimiter=',')


############################################################################################
#                                                                                          #
#                                  Functions definitions                                   #
#                                                                                          #
############################################################################################


def rank_galaxies(galaxies,triple_coincidences,TriggerTime):
    """ Returns a 2D array with ra, dec, cohSNR, gps_time_max_cohSNR ordered for descending cohSNR"""
    ranking=np.array([])
    for galaxy in galaxies:
        cohSNR=np.array([compute_cohSNR(triple,galaxy,TriggerTime) for triple in triple_coincidences])
        max_cohSNR=cohSNR[np.where(cohSNR == max(cohSNR[:,0]))[0][0]] #extract the maximum value among all 
                                                                  #the triple coincidences
        galaxy_SNR=np.append((np.array([galaxy.ra,galaxy.dec,]),max_cohSNR))
        ranking=np.vstack((ranking,galaxy_SNR))
    final_rank=ranking[ranking[:,2].argsort()[::-1]] # sorts the array in descending order of cohSNR 
    return final_rank



def compute_proj_matrix(triple_ev,galaxy,TriggerTime):
    """ Computes the projection matrix M^{mu,nu} (see Harry-Fairhurst).
        Input:  list of events in triple coincidence (sublist of the triple_coinc list);
                interest galaxy and trigger time (is sufficient as reference time for the scope)"""
    galaxy.get_antenna_patterns(TriggerTime)
    wp=np.array([
        [np.sqrt(triple_ev[0].parameters.sigma_sq) * galaxy.antenna_patterns_H[0]],
        [np.sqrt(triple_ev[1].parameters.sigma_sq) * galaxy.antenna_patterns_L[0]],
        [np.sqrt(triple_ev[2].parameters.sigma_sq) * galaxy.antenna_patterns_V[0]]])
    wc=np.array([
        [np.sqrt(triple_ev[0].parameters.sigma_sq) * galaxy.antenna_patterns_H[1]],
        [np.sqrt(triple_ev[1].parameters.sigma_sq) * galaxy.antenna_patterns_L[1]],
        [np.sqrt(triple_ev[2].parameters.sigma_sq) * galaxy.antenna_patterns_V[1]]])
    
    a=np.dot(wp,wp)
    b=np.dot(wc,wc)
    c=np.dot(wp.wc)
    det=a*b-c**2
    zeros=np.zeros((2,2))
    block=np.array([[b, -c],
                    [-c, a] ])
    
    M=(1./det)*np.vstack((np.hstack((block,zeros)),np.hstack((zeros,block))))
    
    return M

   
   
def make_MFO_HLV_array(triple_ev,galaxy,TriggerTime):
    """ Prepares the MFO array that is necessary for cohSNR. Time series are arranged
        according to the time delays between detectors for the given galaxy and trigger time """
    reference_ifo=None
    max_SNR=0
    t_ref=0

    ### extracts the MFO time series af the events and identifies the reference ifo

    for event in triple_ev:
        if event.parameters.rwSNR > max_SNR : # selects the reference ifo according to the value of re-weightedSNR
            max_SNR=event.parameters.rwSNR
            reference_ifo=event.ifo
            t_ref=event.parameters.GTime

    # computes the time of travel between the reference ifo and the 3 ifos
    galaxy.get_relative_time_delays(reference_ifo,TriggerTime)
    [delay_Href,delay_Lref,delay_Vref]=galaxy.time_delays
    max_delay=max(galaxy.time_delays[galaxy.time_delays!=0])
    min_delay=min(galaxy.time_delays[galaxy.time_delays!=0])

    Hrw_timeseries=np.array([
        [galaxy.antenna_patterns_H[0]*triple_ev[0].mfo_data[:,0]],  # Fplus_H  h_phase
        [galaxy.antenna_patterns_H[1]*triple_ev[0].mfo_data[:,0]],  # Fcross_H h_phase
        [galaxy.antenna_patterns_H[0]*triple_ev[0].mfo_data[:,1]],  # Fplus_H  h_quadr
        [galaxy.antenna_patterns_H[1]*triple_ev[0].mfo_data[:,1]]]) # Fcross_H h_quadr 
    Lrw_timeseries=np.array([
        [galaxy.antenna_patterns_L[0]*triple_ev[1].mfo_data[:,0]],
        [galaxy.antenna_patterns_L[1]*triple_ev[1].mfo_data[:,0]],
        [galaxy.antenna_patterns_L[0]*triple_ev[1].mfo_data[:,1]],
        [galaxy.antenna_patterns_L[1]*triple_ev[1].mfo_data[:,1]]])
    Vrw_timeseries=np.array([
        [galaxy.antenna_patterns_V[0]*triple_ev[2].mfo_data[:,0]],
        [galaxy.antenna_patterns_V[1]*triple_ev[2].mfo_data[:,0]],
        [galaxy.antenna_patterns_V[0]*triple_ev[2].mfo_data[:,1]],
        [galaxy.antenna_patterns_V[1]*triple_ev[2].mfo_data[:,1]]])

    #### appends zeros to time series according to the time delay 
    #### in order to prepare them to be summed together     
    
    shift_max=abs(int(np.around(max_delay/delta_t)))
    shift_min=abs(int(np.around(min_delay/delta_t)))

    if max_delay>0 and min_delay>0:
        if delay_Href==0:
            Hrw_timeseries=np.hstack((Hrw_timeseries,np.zeros((4,shift_max))))
        elif delay_Href==max_delay:
            Hrw_timeseries=np.hstack((np.zeros((4,shift_max)),Hrw_timeseries))
        elif delay_Href==min_delay:
            Hrw_timeseries=np.hstack((np.zeros((4,shift_min)),Hrw_timeseries,
                                    np.zeros((4,shift_max-shift_min))))
        if delay_Lref==0:
            Lrw_timeseries=np.hstack((Lrw_timeseries,np.zeros((4,shift_max))))
        elif delay_Lref==max_delay:
            Lrw_timeseries=np.hstack((np.zeros((4,shift_max)),Lrw_timeseries))
        elif delay_Lref==min_delay:
            Lrw_timeseries=np.hstack((np.zeros((4,shift_min)),Lrw_timeseries,
                                    np.zeros((4,shift_max-shift_min))))
        if delay_Vref==0:
            Vrw_timeseries=np.hstack((Vrw_timeseries,np.zeros((4,shift_max))))
        elif delay_Vref==max_delay:
            Vrw_timeseries=np.hstack((np.zeros((4,shift_max)),Vrw_timeseries))
        elif delay_Vref==min_delay:
            Vrw_timeseries=np.hstack((np.zeros((4,shift_min)),Vrw_timeseries,
                                    np.zeros((4,shift_max-shift_min))))
        
    elif max_delay>0 and min_delay<0:
        if delay_Href==0:
            Hrw_timeseries=np.hstack((np.zeros((4,shift_min)),Hrw_timeseries,np.zeros((4,shift_max))))
        elif delay_Href==max_delay:
            Hrw_timeseries=np.hstack((np.zeros((4,shift_max+shift_min)),Hrw_timeseries))
        elif delay_Href==min_delay:
            Hrw_timeseries=np.hstack((Hrw_timeseries,np.zeros((4,shift_max+shift_min))))
        if delay_Lref==0:
            Lrw_timeseries=np.hstack((np.zeros((4,shift_min)),Lrw_timeseries,np.zeros((4,shift_max))))
        elif delay_Lref==max_delay:
            Lrw_timeseries=np.hstack((np.zeros((4,shift_max+shift_min)),Lrw_timeseries))
        elif delay_Vref==min_delay:
            Lrw_timeseries=np.hstack((Lrw_timeseries,np.zeros((4,shift_max+shift_min))))
        if delay_Vref==0:
            Vrw_timeseries=np.hstack((np.zeros((4,shift_min)),Vrw_timeseries,np.zeros((4,shift_max))))
        elif delay_Vref==max_delay:
            Vrw_timeseries=np.hstack((np.zeros((4,shift_max+shift_min)),Vrw_timeseries))
        elif delay_Vref==min_delay:
            Vrw_timeseries=np.hstack((Vrw_timeseries,np.zeros((4,shift_max+shift_min))))

    elif max_delay<0 and min_delay<0:
        if delay_Href==0:
            Hrw_timeseries=np.hstack((np.zeros((4,shift_min)),Hrw_timeseries))
        elif delay_Href==max_delay:
            Hrw_timeseries=np.hstack((np.zeros((4,shift_min-shift_max)),Hrw_timeseries,
                                    np.zeros((4,shift_max))))
        elif delay_Href==min_delay:
            Hrw_timeseries=np.hstack((Hrw_timeseries,np.zeros((4,shift_min))))
        if delay_Lref==0:
            Lrw_timeseries=np.hstack((np.zeros((4,shift_min)),Lrw_timeseries))        
        elif delay_Lref==max_delay:
            Lrw_timeseries=np.hstack((np.zeros((4,shift_min-shift_max)),Lrw_timeseries,
                                    np.zeros((4,shift_max))))
        elif delay_Lref==min_delay:
            Lrw_timeseries=np.hstack((Lrw_timeseries,np.zeros((4,shift_min))))
        if delay_Vref==0:
            Vrw_timeseries=np.hstack((np.zeros((4,shift_min)),Vrw_timeseries))
        elif delay_Vref==max_delay:
            Vrw_timeseries=np.hstack((np.zeros((4,shift_min-shift_max)),Vrw_timeseries,
                                    np.zeros((4,shift_max))))
        elif delay_Vref==min_delay:
            Vrw_timeseries=np.hstack((Vrw_timeseries,np.zeros((4,shift_min))))
    else:# raise error
         pass

    MFO_HLV=Hrw_timeseries+Lrw_timeseries+Vrw_timeseries
   
    return MFO_HLV,time_origin


   
def compute_cohSNR(triple_ev,galaxy,TriggerTime):
    """ computes the coherent snr time series and returns 
        the value and the (approximate) gps time of its maximum """
    MFO_HLV,time_origin = make_MFO_HLV_array(triple_ev,galaxy,TriggerTime)
    M=compute_proj_matrix(triple_ev,galaxy,TriggerTime)

    cohSNR=np.einsum('it,ij,jt->t',MFO_HLV,M,MFO_HLV)  # as Eq.(2.26) in Harry-Fairhurst. t can be viewed as time index

    max_cohSNR=max(cohSNR)
    time_steps=np.where(cohSNR == max_cohSNR)[0][0]
    time_of_max_coh_SNR=time_origin+delta_t*time_steps # time_origin is the gps time associated to the first element of the MFO_HLV array
        

    return max_cohSNR,time_of_max_coh_SNR

   
   
def get_gw_events(tStart,duration,ifo):
    """ Returns a list of GW events, instancies of the GWEvent class, for given input parameters"""
    events=[]
    ffl=vt.FrameFile("Mbta%c1_RX_Clstr.ffl"%(ifo))
    inputFile=ffl.get_frame(TriggerTime) ### dove entra l'informazione dell'ifo??
    iFile=PyFd.FrFileINew(inputFile) ### come trovare l'input file?
    event=PyFd.FrEventReadTF(iFile,"Mbta%c_00-Chi2OK"%(ifo), tStart, duration, 0, 0)
    while(event):
        ev=em.GWEvent(event,inputFile,ifo)
        ev.get_template_params()
        if ev.parameters.mass1 <= 3 and ev.parameters.mass2 <= 3 : # selects only BNS candidates
            ev.get_mfo()
            events.append(ev)
            event=event[0].next
        else:
            event=event[0].next
    return events

   
   
def get_galaxies(galFile):
    """ Input: comma-separated file with ra and dec coordinates of galaxies withing skymap;     
        Output: list of galaxies, instancies of the Galaxy class defined in em_mbta"""
    infile=open(galFile,'r')
    coordinates=np.array([[float(coords.split(',')[0]),float(coords.split(',')[1])]
                        for coords in infile])
    infile.close()
    galaxies=[]
    for ra, dec in coordinates:
        gal=em.Galaxy(galFile,ra,dec)
        galaxies.append(gal)
    return galaxies
