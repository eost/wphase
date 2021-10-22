'''
A simple class that deals with sac files

Written by Z. Duputel, December 2013
        available on GitHub: https://github.com/eost/sacpy.git
'''

############################################################################
#
#                  W phase source inversion package                 
#                               -------------
#
#        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#                      
# (c) California Institute of Technology and Universite de Strasbourg / CNRS 
#                                  April 2013
#
#    Neither the name of the California Institute of Technology (Caltech) 
#    nor the names of its contributors may be used to endorse or promote 
#    products derived from this software without specific prior written 
#    permission
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################


import os,sys
import numpy  as np
import shutil as sh
import scipy.signal as signal
from copy     import deepcopy
from datetime import datetime, timedelta

NVHDR = 6
ITIME = 1

def unpack_c(chararray,strip=True):
    S = ''
    for c in chararray: S+= c.decode('utf-8')
    if strip: S = S.strip()
    return S

def pack_c(char,size):
    S = deepcopy(char)
    c_size = len(char)
    for i in range(size-c_size):
        S = S+' '
    return np.array(S,dtype='c')


class SacError(Exception):
    """
    Raised if the SAC file is corrupted
    """
    pass


class Sac(object):
    '''
    A simple sac class
    '''
    
    def __init__(self,filename=None,datflag=True):
        '''
        Constructor
        Args:
            * filename: read sac filename (optional)
        '''        
        self.delta  =  -12345.
        self.depmin =  -12345.
        self.depmax =  -12345.
        self.scale  =  -12345.
        self.odelta =  -12345.
        self.b      =  -12345.
        self.e      =  -12345.
        self.o      =  -12345.
        self.a      =  -12345.
        self.internal1 = -12345.
        self.t      = np.ones((10,),dtype='float32')*-12345. 
        self.f      =  -12345.
        self.resp   = np.ones((10,),dtype='float32')*-12345. 
        self.stla   =  -12345.
        self.stlo   =  -12345.
        self.stel   =  -12345. 
        self.stdp   =  -12345. 
        self.evla   =  -12345.
        self.evlo   =  -12345.
        self.evel   =  -12345.
        self.evdp   =  -12345.
        self.mag    =  -12345.
        self.user   = np.ones((10,),dtype='float32')*-12345.
        self.dist   =  -12345.
        self.az     =  -12345.
        self.baz    =  -12345.
        self.gcarc  =  -12345.
        self.internal2 = -12345.
        self.internal3 = -12345.
        self.depmen =  -12345.
        self.cmpaz  =  -12345.
        self.cmpinc =  -12345.
        self.xminimum = -12345.
        self.xmaximum = -12345.
        self.yminimum = -12345.
        self.ymaximum = -12345.
        self.nzyear =  -12345
        self.nzjday =  -12345
        self.nzhour =  -12345
        self.nzmin  =  -12345
        self.nzsec  =  -12345
        self.nzmsec =  -12345
        self.nvhdr  =  NVHDR
        self.norid  =  -12345
        self.nevid  =  -12345
        self.npts   =  -12345
        self.internal4 = -12345
        self.nwfid  = -12345
        self.nxsize = -12345
        self.nysize = -12345
        self.iftype =  ITIME
        self.idep   = -12345
        self.iztype = -12345
        self.iinst  = -12345
        self.istreg = -12345
        self.ievreg = -12345
        self.ievtyp = -12345
        self.iqual  = -12345
        self.isynth = -12345
        self.imagtyp = -12345
        self.imagsrc = -12345
        self.leven  = -12345
        self.lpspol = -12345
        self.lovrok = -12345
        self.lcalda = -12345
        self.kstnm  = '-12345'
        self.kevnm  = '-12345'
        self.khole  = '-12345'
        self.ko     = '-12345'     
        self.ka     = '-12345' 
        self.kt = []
        for i in range(10):
            self.kt.append('-12345')
            self.kf     = '-12345'
        self.kuser  = []
        for i in range(3):
            self.kuser.append('-12345')
        self.kcmpnm = '-12345'
        self.knetwk = '-12345'
        self.kdatrd = '-12345'
        self.kinst  = '-12345'
        self.id     = self.knetwk+'_'+self.kstnm+'_'+self.khole+'_'+self.kcmpnm
        self.depvar =  np.array([])

        # Read sac file if filename is specified
        if filename is not None:
            assert os.path.exists(filename), filename+' not found'
            self.read(filename,datflag)

        # Spectrum flag
        self.spec = False

        # All done
        

    def read(self,FILE,npts=None,datflag=True):
        '''
        Read sac file
        Args:
           * FILE: input sac file name
           * npts: number of data points to be read
           * datflag: True: read data, False: read header only
        '''
        # Open file
        fid     = open(FILE,'rb')
        
        # Check endianness
        fid.seek(316,0)
        npts = np.fromfile(fid,'<i4',1)[0]
        fid.seek(0,2)
        fsize = fid.tell()
        if fsize==632+4*npts:
            ftype='<f4'
            itype='<i4'
        elif fsize==632+4*npts.byteswap():
            ftype='>f4'
            itype='>i4'
        else:
            raise SacError("Number of points in header and length of trace inconsistent !")
        
        # Read header
        fid.seek(0,0)
        self.delta     = np.fromfile(fid,ftype,1)[0]
        self.depmin    = np.fromfile(fid,ftype,1)[0]
        self.depmax    = np.fromfile(fid,ftype,1)[0]
        self.scale     = np.fromfile(fid,ftype,1)[0]
        self.odelta    = np.fromfile(fid,ftype,1)[0]
        self.b         = np.fromfile(fid,ftype,1)[0]
        self.e         = np.fromfile(fid,ftype,1)[0]
        self.o         = np.fromfile(fid,ftype,1)[0]
        self.a         = np.fromfile(fid,ftype,1)[0]
        self.internal1 = np.fromfile(fid,ftype, 1)[0]
        self.t         = np.fromfile(fid,ftype,10)
        self.f         = np.fromfile(fid,ftype,1)[0]
        self.resp      = np.fromfile(fid,ftype,10)
        self.stla      = np.fromfile(fid,ftype,1)[0]
        self.stlo      = np.fromfile(fid,ftype,1)[0]
        self.stel      = np.fromfile(fid,ftype,1)[0]
        self.stdp      = np.fromfile(fid,ftype,1)[0]
        self.evla      = np.fromfile(fid,ftype,1)[0]
        self.evlo      = np.fromfile(fid,ftype,1)[0]
        self.evel      = np.fromfile(fid,ftype,1)[0]
        self.evdp      = np.fromfile(fid,ftype,1)[0]
        self.mag       = np.fromfile(fid,ftype,1)[0]
        self.user      = np.fromfile(fid,ftype,  10)
        self.dist      = np.fromfile(fid,ftype,1)[0]
        self.az        = np.fromfile(fid,ftype,1)[0]
        self.baz       = np.fromfile(fid,ftype,1)[0]
        self.gcarc     = np.fromfile(fid,ftype,1)[0]
        self.internal2 = np.fromfile(fid,ftype,1)[0]
        self.internal3 = np.fromfile(fid,ftype,1)[0]
        self.depmen    = np.fromfile(fid,ftype,1)[0]
        self.cmpaz     = np.fromfile(fid,ftype,1)[0]
        self.cmpinc    = np.fromfile(fid,ftype,1)[0]
        self.xminimum  = np.fromfile(fid,ftype,1)[0]
        self.xmaximum  = np.fromfile(fid,ftype,1)[0]
        self.yminimum  = np.fromfile(fid,ftype,1)[0]
        self.ymaximum  = np.fromfile(fid,ftype,1)[0]
        fid.seek(7*4,1)
        self.nzyear    = np.fromfile(fid,itype,1)[0]
        self.nzjday    = np.fromfile(fid,itype,1)[0]
        self.nzhour    = np.fromfile(fid,itype,1)[0]
        self.nzmin     = np.fromfile(fid,itype,1)[0]
        self.nzsec     = np.fromfile(fid,itype,1)[0]
        self.nzmsec    = np.fromfile(fid,itype,1)[0]
        self.nvhdr     = np.fromfile(fid,itype,1)[0]
        self.norid     = np.fromfile(fid,itype,1)[0]
        self.nevid     = np.fromfile(fid,itype,1)[0]
        self.npts      = np.fromfile(fid,itype,1)[0]
        self.internal4 = np.fromfile(fid,itype,1)[0]
        self.nwfid     = np.fromfile(fid,itype,1)[0]
        self.nxsize    = np.fromfile(fid,itype,1)[0]
        self.nysize    = np.fromfile(fid,itype,1)[0]
        fid.seek(4,1);
        self.iftype    = np.fromfile(fid,itype,1)[0]
        self.idep      = np.fromfile(fid,itype,1)[0]
        self.iztype    = np.fromfile(fid,itype,1)[0]
        fid.seek(4,1);
        self.iinst     = np.fromfile(fid,itype,1)[0]
        self.istreg    = np.fromfile(fid,itype,1)[0]
        self.ievreg    = np.fromfile(fid,itype,1)[0]
        self.ievtyp    = np.fromfile(fid,itype,1)[0]
        self.iqual     = np.fromfile(fid,itype,1)[0]
        self.isynth    = np.fromfile(fid,itype,1)[0]
        self.imagtyp   = np.fromfile(fid,itype,1)[0]
        self.imagsrc   = np.fromfile(fid,itype,1)[0]
        fid.seek(8*4,1);
        self.leven     = np.fromfile(fid,itype,1)[0]
        self.lpspol    = np.fromfile(fid,itype,1)[0]
        self.lovrok    = np.fromfile(fid,itype,1)[0]
        self.lcalda    = np.fromfile(fid,itype,1)[0]
        fid.seek(4,1);
        self.kstnm     = unpack_c(np.fromfile(fid,'c',8))
        self.kevnm     = unpack_c(np.fromfile(fid,'c',16),False)
        self.khole     = unpack_c(np.fromfile(fid,'c',8))
        self.ko        = unpack_c(np.fromfile(fid,'c',8))
        self.ka        = unpack_c(np.fromfile(fid,'c',8))
        for i in range(10):
            self.kt[i] = unpack_c(np.fromfile(fid,'c',8))
        self.kf = unpack_c(np.fromfile(fid,'c',8))
        for i in range(3):
            self.kuser[i] = unpack_c(np.fromfile(fid,'c',8))
        self.kcmpnm = unpack_c(np.fromfile(fid,'c',8))
        self.knetwk = unpack_c(np.fromfile(fid,'c',8))
        self.kdatrd = unpack_c(np.fromfile(fid,'c',8))
        self.kinst  = unpack_c(np.fromfile(fid,'c',8))
        self.e = self.b + float(self.npts-1) * self.delta
        if self.khole=='' or self.khole=='-12345':
            self.khole = '--'
        self.id = self.knetwk+'_'+self.kstnm+'_'+self.khole+'_'\
                     +self.kcmpnm                        

        # Don't read waveform
        if not datflag: 
            fid.close()
            # All done
            return

        # Read waveform
        fid.seek(632,0);
        if npts is None or npts < 0 or npts > self.npts:
            npts = self.npts
        else:
            self.npts = int(npts)            
        if self.npts > 0:
            self.depvar = np.fromfile(fid,ftype,self.npts)
        fid.close()

        # Re-assign min/max amplitudes and end time
        self.depmin  = self.depvar.min()
        self.depmax  = self.depvar.max()
        self.e       = self.b + float(self.npts - 1) * self.delta
        
        # All done


    def write(self,FILE):
        '''
        Write sac file
        Args:
           * FILE: output sac file name
        '''

        # Check that we are in the time domain
        assert not self.spec, "Can only save seismograms in the time-domain"
        
        # convert to list
        if type(self.depvar)==list:
            self.depvar = np.array(self.depvar)
        
        # Dummy variables
        dumi = np.array(-12345    ,dtype='int32')
        dumf = np.array(-12345.0  ,dtype='float32')
        dumc = np.array('-12345  ',dtype='c')
        
        # Re-assign min/max amplitudes and end time
        self.depmin  = self.depvar.min()
        self.depmax  = self.depvar.max()
        self.e = self.b + float(self.npts - 1) * self.delta

        # Write file
        fid = open(FILE,'wb')
        np.array(self.delta,dtype='float32').tofile(fid)
        np.array(self.depmin,dtype='float32').tofile(fid)
        np.array(self.depmax,dtype='float32').tofile(fid)
        np.array(self.scale,dtype='float32').tofile(fid)
        np.array(self.odelta,dtype='float32').tofile(fid)
        np.array(self.b,dtype='float32').tofile(fid)
        np.array(self.e,dtype='float32').tofile(fid)
        np.array(self.o,dtype='float32').tofile(fid)
        np.array(self.a,dtype='float32').tofile(fid)
        np.array(self.internal1,dtype='float32').tofile(fid)
        np.array(self.t,dtype='float32').tofile(fid)
        np.array(self.f,dtype='float32').tofile(fid)
        np.array(self.resp,dtype='float32').tofile(fid)
        np.array(self.stla,dtype='float32').tofile(fid)
        np.array(self.stlo,dtype='float32').tofile(fid)
        np.array(self.stel,dtype='float32').tofile(fid)
        np.array(self.stdp,dtype='float32').tofile(fid)
        np.array(self.evla,dtype='float32').tofile(fid)
        np.array(self.evlo,dtype='float32').tofile(fid)
        np.array(self.evel,dtype='float32').tofile(fid)
        np.array(self.evdp,dtype='float32').tofile(fid)
        np.array(self.mag,dtype='float32').tofile(fid)
        np.array(self.user,dtype='float32').tofile(fid)
        np.array(self.dist,dtype='float32').tofile(fid)
        np.array(self.az,dtype='float32').tofile(fid)
        np.array(self.baz,dtype='float32').tofile(fid)
        np.array(self.gcarc,dtype='float32').tofile(fid)
        np.array(self.internal2,dtype='float32').tofile(fid)
        np.array(self.internal3,dtype='float32').tofile(fid)
        np.array(self.depmen,dtype='float32').tofile(fid)
        np.array(self.cmpaz,dtype='float32').tofile(fid)
        np.array(self.cmpinc,dtype='float32').tofile(fid)
        np.array(self.xminimum,dtype='float32').tofile(fid)
        np.array(self.xmaximum,dtype='float32').tofile(fid)
        np.array(self.yminimum,dtype='float32').tofile(fid)
        np.array(self.ymaximum,dtype='float32').tofile(fid)  
        for i in range(7):
            dumf.tofile(fid)
        np.array(self.nzyear,dtype='int32').tofile(fid)  
        np.array(self.nzjday,dtype='int32').tofile(fid)  
        np.array(self.nzhour,dtype='int32').tofile(fid)
        np.array(self.nzmin,dtype='int32').tofile(fid)  
        np.array(self.nzsec,dtype='int32').tofile(fid)  
        np.array(self.nzmsec,dtype='int32').tofile(fid)  
        np.array(self.nvhdr,dtype='int32').tofile(fid)  
        np.array(self.norid,dtype='int32').tofile(fid)  
        np.array(self.nevid,dtype='int32').tofile(fid)
        np.array(self.npts,dtype='int32').tofile(fid) 
        np.array(self.internal4,dtype='int32').tofile(fid)  
        np.array(self.nwfid,dtype='int32').tofile(fid)
        np.array(self.nxsize,dtype='int32').tofile(fid)
        np.array(self.nysize,dtype='int32').tofile(fid)
        dumi.tofile(fid)
        np.array(self.iftype,dtype='int32').tofile(fid)
        np.array(self.idep,dtype='int32').tofile(fid)
        np.array(self.iztype,dtype='int32').tofile(fid)
        dumi.tofile(fid)
        np.array(self.iinst,dtype='int32').tofile(fid)
        np.array(self.istreg,dtype='int32').tofile(fid)
        np.array(self.ievreg,dtype='int32').tofile(fid)
        np.array(self.ievtyp,dtype='int32').tofile(fid)
        np.array(self.iqual,dtype='int32').tofile(fid)
        np.array(self.isynth,dtype='int32').tofile(fid)
        np.array(self.imagtyp,dtype='int32').tofile(fid)
        np.array(self.imagsrc,dtype='int32').tofile(fid)
        for i in range(8):
            dumi.tofile(fid)
        np.array(self.leven,dtype='int32').tofile(fid)
        np.array(self.lpspol,dtype='int32').tofile(fid)
        np.array(self.lovrok,dtype='int32').tofile(fid)
        np.array(self.lcalda,dtype='int32').tofile(fid)
        dumi.tofile(fid)
        pack_c(self.kstnm,8).tofile(fid)  
        pack_c(self.kevnm,16).tofile(fid)  
        pack_c(self.khole,8).tofile(fid)
        pack_c(self.ko,8).tofile(fid)
        pack_c(self.ka,8).tofile(fid)
        for i in range(10):
            pack_c(self.kt[i],8).tofile(fid)
        pack_c(self.kf,8).tofile(fid)
        for i in range(3):
            pack_c(self.kuser[i],8).tofile(fid)
        pack_c(self.kcmpnm,8).tofile(fid)
        pack_c(self.knetwk,8).tofile(fid)
        pack_c(self.kdatrd,8).tofile(fid)
        pack_c(self.kinst,8).tofile(fid)

        # Write data
        np.array(self.depvar,dtype='float32').tofile(fid)
        fid.close()
                
        # All done

        
    def rsac(self,FILE,npts=None,datflag=True):
        '''
        Clone of self.read()
        '''
        sys.stderr.write('FutureWarning: sac.rsac will be replaced by sac.read in the future\n')
        self.read(FILE,npts=None,datflag=True)

        # All done

        
    def wsac(self,FILE):
        '''
        Clone of self.write
        '''
        sys.stderr.write('FutureWarning: sac.wsac will be replaced by sac.write in the future\n')
        self.write(FILE)

        # All done

        
    def getnzdatetime(self):
        '''
        Get the reference datetime
        '''
        # Reference datetime
        nzjday = timedelta(int(self.nzjday)-1)
        nzusec = int(self.nzmsec*1e3)
        nztime = datetime(self.nzyear,1,1,self.nzhour,self.nzmin,self.nzsec,nzusec)+nzjday

        # All done
        return nztime

    def getdatetime(self,time_param):
        '''
        Get a datetime object from a sac time parameter
        Args:
            * time_param: sac time parameter defined as a time 
                difference (in seconds) with respect to nztime
                (e.g., o, b, e, t[0], ...)
        '''
        # Check that time_param is set
        assert int(time_param) != -12345, 'Parameter not set'
        # Get nzdatetime
        nztime = self.getnzdatetime()
        # Define dt_o
        dt_o = timedelta(seconds=int(time_param))
        # Get otime
        otime = nztime + dt_o
        # All done 
        return otime

    def getodatetime(self):
        '''
        Get a origin datetime object
        '''
        # All done
        return self.getdatetime(self.o)

    def getbdatetime(self):
        '''
        Get a begin datetime object
        '''
        # All done
        return self.getdatetime(self.b)

    def getedatetime(self):
        '''
        Get an end datetime object
        '''
        # Recompute end time (just to be sure everything is all right)
        self.e = self.b + float(self.npts-1) * self.delta
        # All done
        return self.getdatetime(self.e)

    def getarrivaldatetimes(self):
        '''
        Get a dictionary of arrival time objects
        '''
        # Loop on travel times
        phase_datetimes = {}
        for i,t in enumerate(self.t):
            if int(t)!=-12345:
                phase_datetimes[self.kt[i]] = self.getdatetime(t)
        # All done
        return phase_datetimes

    def setotime(self,otime):
        '''
        Set o 
        Arg:
            * otime: datetime instance
        '''
        # Get nzdatetime
        nztime = self.getnzdatetime()
                
        # Time difference is o
        self.o  = np.float32((otime-nztime).total_seconds())

        # All done
        
    def setarrivaltimes(self,phase_dict):
        '''
        Set t and kt 
        Arg:
          * phase_dict: phase pick dictionary {name: arrival_datetime)                
        '''

        # Get nzdatetime
        nztime = self.getnzdatetime()

        # Loop over phase names and arrival times 
        i = 0
        for pname, ptime in phase_dict.items():

            # Phase name
            assert len(pname)<8, 'phase name %s is too long'%(pname)
            self.kt[i] = pname

            # Arrival time                        
            self.t[i]  = np.float32((ptime-nztime).total_seconds())
            i += 1
                        
        # All done

    
    def integrate(self):
        '''
        Integration using the traperoidal rule
        '''

        # Integration
        w  = self.depvar.copy()
        wi = self.depvar.cumsum()
        wi = (2*wi[1:]-(w[0]+w[1:]))*self.delta/2.
        self.depvar = wi.copy()

        # Re-assign b, e, npts, min/max amplitudes
        self.b += self.delta/2.
        self.npts -= 1
        self.e = self.b + float(self.npts - 1) * self.delta
        self.depmin  = self.depvar.min()
        self.depmax  = self.depvar.max()
        
        # All done
        return

    def derivate(self):
        '''
        Derivate using a two point difference operator
        '''
        # Two-point derivative
        self.depvar = (self.depvar[1:]-self.depvar[:-1])/self.delta

        # Adjusting the number of points and begin/end times
        self.npts -= 1
        self.b += 0.5 * self.delta
        self.e = self.b + float(self.npts - 1) * self.delta

        # All done


    def isempty(self):
        '''
        Check if important attributes are there
        '''    
        if (self.npts < 0) or (self.delta < 0) or (not self.depvar.size):
            return True

        # All done
        return False

    def resetdepmindepmax(self):
        '''
        Reset depmin/depmax attributes
        '''
        self.depmin = self.depvar.min()
        self.depmax = self.depvar.max() 
        
    def filter(self, freq, order=4, btype='lowpass'):
        '''
        Bandpass filter the data using a butterworth filter
        Args:
            * freq: A scalar or length-2 sequence giving the critical frequencies (in Hz)
            * order:  Order of the filter.
            * btype: {'lowpass', 'highpass', 'bandpass', 'bandstop'}, optional
              (default is 'lowpass')
        '''
        
        # Check that headers are correct
        assert not self.isempty(), 'Some sac attributes are missing (e.g., npts, delta, depvar)'

        # Filter design
        if type(freq) is list:
            freq = np.array(freq)
        Wn = freq * 2. * self.delta # Normalizing frequencies
        sos = signal.butter(order, Wn, btype, output='sos')
        
        # Filter waveform
        depvar = signal.sosfilt(sos, self.depvar)
        self.depvar = depvar.astype('float32')

        # Reset depmin/depmax
        self.resetdepmindepmax()

        # All done
        return

    def pad(self,tmin = None, tmax = None):
        '''
        Padding data with zeros
        if tmin < self.b (beginning), adding zeros at the beginning
        if tmax > self.e (end), adding zeros at the end
        '''
        # Get trace beginning and end
        self.e = self.b + float(self.npts - 1) * self.delta
        tb = self.b
        te = self.e

        # Set the pad width
        nbeg = 0
        nend = 0
        if tmin is not None and tmin < tb:
            nbeg = int(np.ceil((tb-tmin)/self.delta))
        if tmax is not None and tmax > te:
            nend = int(np.ceil((tmax-te)/self.delta))

        # Zero padding
        gout = np.pad(self.depvar,((nbeg,nend),),mode="constant",constant_values=0.)
        self.npts = len(gout)
        self.b = self.b - nbeg * self.delta
        self.e = self.e + nend * self.delta
        self.depvar = gout.copy()

        # Reset depmin/depmax
        self.resetdepmindepmax()

        # All done
        return


    def cut(self,beg,end):
        '''
        Cut data given beg/end datetime
        '''

        # Get time vector
        nztime = self.getnzdatetime()
        b  = self.b
        dt = self.delta
        time = np.arange(self.npts)*dt + b

        # Extract data between beg and end
        tbeg = (beg-nztime).total_seconds()
        tend = (end-nztime).total_seconds()
        i = np.where((time>=tbeg)*(time<=tend))[0]

        # Only
        self.b    = time[i[0]]
        self.e    = time[i[-1]]
        self.npts = len(i)
        self.depvar = self.depvar[i]

        # Reset depmin/depmax
        self.resetdepmindepmax()

        # All done
        return

    def fft(self,n=None):
        '''
        Compute fourier transform and return the seismogram spectrum
        Output: Seismogram spectrum in the frequency domain (type: seismogram)        
        Args:
            n: Number of points for the fft (default is npts)
        '''
        spectrum = self.copy()
        spectrum.spec = True
        spectrum.depvar = np.fft.rfft(self.depvar,n=n) 
        if n is not None:
            spectrum.npts = n
        
        # All done
        return spectrum

    def ifft(self):
        '''
        Compute the inverse fourrier transform and returns the seismogram spectrum
        Output: Seismogram in the time domain (type: seismogram)
        '''
        seis = self.copy()
        seis.spec = False
        seis.depvar=None
        seis.depvar = np.fft.irfft(self.depvar) 
        
        # All done
        return seis

    def freq(self):
        '''
        Returns the frequency vector of the current data
        '''
        freq = np.fft.rfftfreq(self.npts,d=self.delta)
        # All done        
        return freq

    def evalresp(self,PZ):
        '''
        Return frequency response
        Args:
            * PZ: dictionary including 'poles', 'zeros' and 'Const'
        '''
        # Check that the PZ dictionary is complete
        assert 'zeros' in PZ, 'zeros key must be specified in the PZ dictionary'
        assert 'poles' in PZ, 'poles key must be specified in the PZ dictionary'
        assert 'Const' in PZ, 'Const key must be specified in the PZ dictionary'

        # Evaluate the response in the frequency domain
        s = 2.j*np.pi*self.freq()
        resp = np.ones(s.shape,dtype=np.complex128)*PZ['Const']
        for z in PZ['zeros']: resp *= s-z
        for p in PZ['poles']: resp /= s-p

        # All done
        return resp

    def convresp(self,PZ):
        '''
        Convolve with instrument response
        Args:
            * PZ: dictionary including 'poles', 'zeros' and 'Const'        
        '''
        npts = self.npts
        # Trivial dtrend
        self.depvar -= self.depvar[0]+np.arange(npts)*(self.depvar[-1]-self.depvar[0])/(npts-1)
        # Zero padding
        self.pad(tmax=self.b+2.*self.npts*self.delta)
        # Evaluate the instrument response from Poles and Zeros
        resp = self.evalresp(PZ)
        # Convolve with the instrument response
        self.depvar = np.fft.irfft(resp*np.fft.rfft(self.depvar))[:npts]
        self.npts = npts
        self.e    = self.b + float(self.npts)*self.delta
        self.depmin = self.depvar.min()
        self.depmax = self.depvar.max()        
        # Reset depmin/depmax
        self.resetdepmindepmax()

        # All done
        return

    def deconvresp(self,PZ,filtfreq=None):
        '''
        Deconvolve with instrument response
        Args:
            * PZ: dictionary including 'poles', 'zeros' and 'Const'
            * filtfreq: list (or array) of 4 frequencies f1<f2<f3<f4 to
                cope with zero response at zero frequency. It will
                apply a cosine high-pass cosine taper between f1 and f2 
                and a low-pass cosine taper between f3 and f4
        ''' 

        # Zero padding
        npts = self.npts
        npts2 = int(2**np.ceil(np.log2(abs(npts))))
        self.pad(tmax=self.b+npts2*self.delta)

        # Design cosine-tapered filter
        freq = self.freq()
        filt = np.ones(freq.shape)
        if filtfreq is not None:
            # Check filtfreq input
            filtfreq=np.array(filtfreq)
            assert filtfreq.size==4, 'filtfreq must be a list or array including 4 frequencies'
            assert (np.sort(filtfreq)==filtfreq).all(), 'f1<f2<f3<f4 in filtfreq not verified'
            f1,f2,f3,f4 = filtfreq

            # Taper bounds
            i1 = np.where(freq<f1)[0]
            i2 = np.where((freq>=f1)*(freq<=f2))[0]
            i3 = np.where((freq>=f3)*(freq<=f4))[0]
            i4 = np.where(freq>f4)[0]
            df1 = f2-f1
            df2 = f4-f3

            # filter
            filt[i1] *= 0.
            filt[i4] *= 0.
            filt[i2] *= 0.5*(1+np.cos(np.pi/df1*(freq[i2]-f1-df1)))
            filt[i3] *= 0.5*(1+np.cos(np.pi/df2*(freq[i3]-f3)))

        # Evaluate instrument response from Poles and Zeros
        resp = self.evalresp(PZ)

        # Remove instrument and filter
        if freq[0] == 0.:
            resp[0] = 1.
            filt[0] = 0.
        F = np.fft.rfft(self.depvar)*(filt/resp)
        self.depvar = np.fft.irfft(F)[:npts]
        self.npts = npts

        # Reset depmin/depmax
        self.resetdepmindepmax()
    
        # All done
        return

    def rmean(self):
        '''
        Remove means
        '''
        self.depvar -= self.depvar.mean()

        # All done
        return

    def detrend(self,Npoints=None,poly_order=1):
        '''
        Remove linear trend along axis from data
        Args:
            * Npoints: (optional) remove linear trend from linear fit only on the first Npoints
            * poly_order: Polynomial order
        '''
        # Linear fit
        x = self.time()-self.b
        y = self.depvar
        if Npoints:
            x = x[:Npoints]
            y = self.depvar[:Npoints]
        p = np.poly1d(np.polyfit(x,y,poly_order))

        # Remove trend
        self.depvar -= p(self.time()-self.b)

        # Reset depmin/depmax
        self.resetdepmindepmax()
        
        # All done
        return

    def taper(self,alpha=0.1):
        '''
        Applies a symmetric cosine taper to each end of data
        Args:
            * alpha: taper width (value between 0. and 1.)
        '''
        H = signal.tukey(self.npts,alpha=alpha,sym=True)
        self.depvar *= H
        # All done
        return

    def time(self,datetime_format=False):
        '''
        Returns the time vector of the current data relative to nztime
        '''
        time = np.arange(self.npts)*self.delta + self.b
        if datetime_format:
            r = self.getnzdatetime()
            time = np.array([r+timedelta(seconds=t) for t in time])

        # All done
        return time

    def plot(self,ptype=None,xlog=False,ylog=False,**kwargs):
        '''
        Plot the seismogram or spectrum
        Args: All arguments are optional
            - ptype: plot type can be None, 'amp' for absolute amplitude, 'pha' for the phase, 
                     'real' for the real part or 'imag' for the imaginary part.
            - xlog: if True use a log scale on the x axis
            - ylog: if True use a log scale on the y axis
            - *kwargs* can be used to set line properties in pyplot commands (see help of plt.plot)
        examples:
                s.plot(color='r') or s.plot(color='red') will plot the seismogram with a red line
        Use plt.show() to show the corresponding figure
        '''

        # Import the matplotlib module
        import matplotlib.pyplot as plt

        # Check attributes        
        assert not self.isempty(),'Some attributes are missing (e.g., npts, delta, depvar)'

        # Time or frequency vector
        if self.spec is False: # Time vector
            x = self.time()
            xlabel = 'Time, sec'
        else: # Frequency vector
            x = self.freq()
            xlabel = 'Freq., Hz'  
        # What do we want to plot?
        ylabel = 'Amplitude'        
        if ptype is None and not self.spec: # Standard seismogram plot
            y = self.depvar
        elif (ptype is None and self.spec) or ptype == 'amp':  # Amplitude
            y = np.abs(self.depvar)
        elif ptype == 'pha':     # Phase
            y = np.angle(self.depvar)
            ylabel = 'Phase'            
        elif ptype == 'real': # Real part
            y = np.real(self.depvar)
            ylabel = 'Real part amplitude'
        elif ptype == 'imag':
            y = np.imag(self.depvar)
            ylabel = 'Imag. part amplitude'            
        else:
            print('Error: ptype should be None, amp, pha, real or imag')
            return 1        

        # Do we use log scale?
        plotf = plt.plot  # Default: no log scale
        if xlog and ylog: # loglog scale
            plotf=plt.loglog
        elif xlog:        # x log scale
            plotf=plt.semilogx
        elif ylog:        # y log scale
            plotf=plt.semilogy
        
        # Plot seismogram
        lines = plotf(x,y,**kwargs)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        # All done
        return lines    

    def __add__(self, other):
        '''
        Addition operation.         
        other can be:
          - sacpy.sac object
          - list or ndarray
          - real number (float or int)
        '''        

        # Check if the operation can be done
        accepted=(self.__class__,int,float,list,np.ndarray,np.float32,np.float64)
        assert isinstance(other,accepted), 'Unsuported type'
        
        # Copy current object
        res  = self.copy()
        flag = False

        # Adding two sac files
        if isinstance(other,self.__class__):
            assert self.npts  == other.npts,  'Header field mismatch: npts'
            assert self.delta == other.delta, 'Header field mismatch: delta'
            assert self.b     == other.b,     'Header field mismatch: b'
            assert self.e     == other.e,     'Header field mismatch: e'          
            res.depvar += other.depvar
            flag = True

        # Adding array or list
        if isinstance(other,(list,np.ndarray)):
            assert len(other)==self.npts, 'Header field mismatch: npts'
            res.depvar += other
            flag = True

        # Adding real number
        if isinstance(other,(int,float,np.float32,np.float64)):
            res.depvar += other
            flag = True

        # Re-assign min and max amplitudes
        res.depmin  = res.depvar.min()
        res.depmax  = res.depvar.max()
        
        # Check that operation was done
        assert flag, 'Operation could not be completed'
        
        # All done
        return res

    def __sub__(self, other):
        '''
        Substraction operation.         
        other can be:
          - sacpy.sac object
          - list or ndarray
          - real number (float or int)
        '''        

        # Check if the operation can be done
        accepted=(self.__class__,int,float,list,np.ndarray,np.float32,np.float64)
        assert isinstance(other,accepted), 'Unsuported type'
        
        # Copy current object
        res  = self.copy()
        flag = False

        # Adding two sac files
        if isinstance(other,self.__class__):
            assert self.npts  == other.npts,  'Header field mismatch: npts'
            assert self.delta == other.delta, 'Header field mismatch: delta'
            assert self.b     == other.b,     'Header field mismatch: b'
            assert self.e     == other.e,     'Header field mismatch: e'          
            res.depvar -= other.depvar
            flag = True

        # Adding array or list
        if isinstance(other,(list,np.ndarray)):
            assert len(other)==self.npts, 'Header field mismatch: npts'
            res.depvar -= other
            flag = True

        # Adding real number
        if isinstance(other,(int,float,np.float32,np.float64)):
            res.depvar -= other

        # Re-assign min and max amplitudes
        res.depmin  = res.depvar.min()
        res.depmax  = res.depvar.max()
        
        # Check that operation was done
        assert flag, 'Operation could not be completed'
        
        # All done
        return res    

    def __mul__(self, other):
        '''
        Multiplication operation.         
        other can be:
          - sacpy.sac object
          - list or ndarray
          - real number (float or int)
        '''        

        # Check if the operation can be done
        accepted=(self.__class__,int,float,list,np.ndarray,np.float32,np.float64)
        assert isinstance(other,accepted), 'Unsuported type'
        
        # Copy current object
        res  = self.copy()
        flag = False

        # Multiplying two sac files
        if isinstance(other,self.__class__):
            assert self.npts  == other.npts,  'Header field mismatch: npts'
            assert self.delta == other.delta, 'Header field mismatch: delta'
            assert self.b     == other.b,     'Header field mismatch: b'
            assert self.e     == other.e,     'Header field mismatch: e'          
            res.depvar *= other.depvar
            flag = True

        # Multiplying by an array or a list
        if isinstance(other,(list,np.ndarray)):
            assert len(other)==self.npts, 'Header field mismatch: npts'
            res.depvar *= other
            flag = True

        # Multiplying by a real number
        if isinstance(other,(int,float,np.float32,np.float64)):
            res.depvar *= other
            flag = True

        # Re-assign min and max amplitudes
        res.depmin  = res.depvar.min()
        res.depmax  = res.depvar.max()
        
        # Check that operation was done
        assert flag, 'Operation could not be completed'
        
        # All done
        return res

    def __div__(self, other):
        '''
        Multiplication operation.         
        other can be:
          - sacpy.sac object
          - list or ndarray
          - real number (float or int)
        ''' 

        # Check if the operation can be done
        accepted=(self.__class__,int,float,list,np.ndarray,np.float32,np.float64)
        assert isinstance(other,accepted), 'Unsuported type'
        
        # Copy current object
        res  = self.copy()
        flag = False

        # Dividing by a sac file
        if isinstance(other,self.__class__):
            assert self.npts  == other.npts,  'Header field mismatch: npts'
            assert self.delta == other.delta, 'Header field mismatch: delta'
            assert self.b     == other.b,     'Header field mismatch: b'
            assert self.e     == other.e,     'Header field mismatch: e'          
            res.depvar /= other.depvar
            flag = True

        # Dividing by an array or a list
        if isinstance(other,(list,np.ndarray)):
            assert len(other)==self.npts, 'Header field mismatch: npts'
            res.depvar /= other
            flag = True

        # Dividing by a real number
        if isinstance(other,(int,float,np.float32,np.float64)):
            res.depvar /= other
            flag = True

        # Re-assign min and max amplitudes
        res.depmin  = res.depvar.min()
        res.depmax  = res.depvar.max()
        
        # Check that operation was done
        assert flag, 'Operation could not be completed'
        
        # All done
        return res

    def copy(self):
        '''
        Returns a copy of the sac object
        '''
        # All done
        return deepcopy(self)                

def zero_pad_start(t,sac,t0):
    tmin = t[0]
    dt   = sac.delta
    tpad = np.arange(tmin,t0-dt,-dt)
    if tpad[-1]>t0:
        tpad = np.append(tpad,tpad[-1]-dt)
    tpad = tpad[1:]
    tpad = tpad[::-1]
    tout = np.append(tpad,t)
    gout = np.append(0.0*tpad,sac.depvar)
    # all done
    return tout,gout

