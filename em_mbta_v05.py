import numpy as np
import astropy
import collections
import PyFd
import pycbc.detector


H1=pycbc.detector.Detector('H1')
L1=pycbc.detector.Detector('L1')
V1=pycbc.detector.Detector('V1')


class Galaxy(object):

    def __init__(self,galFile,ra,dec):
        self.source=str(galFile)
        self.ra=ra
        self.dec=dec
        self.time_delays=None
        self.antenna_patterns_H=None
        self.antenna_patterns_L=None
        self.antenna_patterns_V=None
        

    def get_relative_time_delays(self, reference_ifo, gps_time):
        """ Input: GPS Trigger time in seconds and str among (H,L,V) to indicate the reference ifo 
                   for computing time shifts, i.e. that with greatest coincSNR.
            Output: numpy array with
                     time-delay-from-ref-to-H, time-delay-from-ref-to-L,
                     time-delay-from-ref-to-V.
            The time delays are the delays from the reference detector for a signal with the given 
            sky location; i.e. return t1-t2 where t1 is the arrival time at the "current detector"
            (H for dtH, V for dtV etc) and t2 is the arrival time in the reference ifo.
        """
        ref_ifo=pycbc.detector('{ifo}1'.format(ifo=reference_ifo))

        def dtH(ra,dec):
            return H1.time_delay_from_other_detector(other_detector=ref_ifo, 
                                right_ascension=ra,declination=dec, t_gps=gps_time)
        def dtL(ra,dec):
            return L1.time_delay_from_other_detector(other_detector=ref_ifo, right_ascension=ra,
                                                    declination=dec, t_gps=gps_time)
        def dtV(ra,dec):
            return V1.time_delay_from_other_detector(other_detector=ref_ifo, right_ascension=ra,
                                                    declination=dec, t_gps=gps_time)
            
        self.time_delays=np.array([dtH(self.ra,self.dec),dtL(self.ra,self.dec),
                                    dtV(self.ra,self.dec)])

    def get_antenna_patterns(self, gps_time):
        """ Return the detectors response functions in an array [Fplus,Fcross] for
            each detector at the given gps time for the galaxies in the skymap"""

        def HFpc(ra,dec):
            """ return an array with the H1's fplus and fcross polarization factors 
                for this sky location """
            Fplus,Fcross=H1.antenna_pattern(right_ascension=ra, declination=dec, polarization=0,
                                            t_gps=gps_time)
            return np.array([Fplus,Fcross])

        def LFpc(ra,dec):
            Fplus,Fcross=L1.antenna_pattern(right_ascension=ra, declination=dec, polarization=0,
                                            t_gps=gps_time)
            return np.array([Fplus,Fcross])

        def VFpc(ra,dec):
            Fplus,Fcross=V1.antenna_pattern(right_ascension=ra, declination=dec, polarization=0,
                                            tgps=gps_time)
            return np.array([Fplus,Fcross])


        self.antenna_patterns_H=HFpc(self.ra,self.dec)
        self.antenna_patterns_L=LFpc(self.ra,self.dec)
        self.antenna_patterns_V=VFpc(self.ra,self.dec)


class GWEvent(object):

    Params=collections.namedtuple('Parameters',['tmplt_index','mass1','mass2',
                                    'MChirp','spin1z','spin2z','rwSNR','sigma_sq'])    
                                    #GTime is the vector time origin GPS time [s]

    def __init__(self,event,inputFile,ifo):
        self.event=event
        self.source=str(inputFile)
        self.ifo=str(ifo)
        self.i_file=PyFd.FrFileINew(self.source)
        self.vect=None
        self.mfo_data=None
        self.parameters=None

    def get_mfo(self):
        """
            Output: 2D numpy array with the matched filtering output in phase and quadrature
        """
        PyFd.FrEventReadData(self.i_file,self.event);
        self.vect=PyFd.FrEventGetVectF(self.event,"{ifo}1:MFO_measured_P".format(ifo=self.ifo))
        Quad_vect=PyFd.FrEventGetVectF(self.event,"{ifo}1:MFO_measured_Q".format(ifo=self.ifo))
        Phase_vect=self.vect

        A=Phase_vect[0].nData
        B=Quad_vect[0].nData

        if A == B:
            self.mfo_data=np.array([[Phase_vect[0].dataF[i],Quad_vect[0].dataF[i]]
                                for i in np.arange(A)])
        else:
            print("Lenghts of MFO time series are not the same for the two polarizations; The longest one will be cut to match the length of the shortest")
            self.mfo_data=np.array([[Phase_vect[0].dataF[i],Quad_vect[0].dataF[j]]
                                for i,j in zip(np.arange(A),np.arange(B))])

    def get_template_params(self):
        """ Output: named tuple containing the values of the following parameters
            for the template of the event in consideration: template index, 
            mass1, mass2, MChirp, spin1z, spin2z,re-weighted SNR and 
            the sensibility of the detector (sigma_sq).
        """
        self.parameters=self.Params(self.event[0].parameters[0],
                                self.event[0].parameters[2],self.event[0].parameters[3],
                                self.event[0].parameters[4],self.event[0].parameters[5],
                                self.event[0].parameters[6],self.event[0].parameters[19],
                                self.event[0].parameters[18])
