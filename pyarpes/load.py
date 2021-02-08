import numpy as np
import pandas as pd
import h5py as h5
from scipy.constants import pi
#from scipy.constants import m_e
#from scipy.constants import hbar
#from scipy.constants import e

def infoThemis(fname):
    """
    Exracting metadata from a data file from the laser-based ARPES setup at KTH.
    [M. H. Berntsen,Berntsen, O. Götberg, O. Tjernberg, Rev. Sci. Instrum. 82, 095113 (2011)]

    Input:
            fname = filename of file, including path
    """

    with h5.File(fname,'r') as f:
        columns = ['Name','Lens mode','Kinetic energy','Pass energy','Type','Comment']
        name = list(f.keys())
        lensmode = [] 
        Ekin = [] 
        Ep = [] 
        comment = [] 
        data_type = []

        for key in f.keys():
            lensmode.append(f[key].attrs['lensmode'])
            Ekin.append(f[key].attrs['Kinetic Energy'])
            Ep.append(f[key].attrs['Pass Energy'])
            comment.append(f[key].attrs['TITLE'])
            for dset in f[key]:
                data_type.append(dset)


        
        
        df = pd.DataFrame(np.array([name,lensmode,Ekin,Ep,data_type,comment]).transpose(),columns=columns)
        df['Comment'] = df['Comment'].str.replace('\n','; ')
        pd.set_option('display.max_colwidth',None)

    return df



def loadThemis(fname,dataset,**kwargs):
    """
    Loading data from the laser-based ARPES setup at KTH.
    [M. H. Berntsen,Berntsen, O. Götberg, O. Tjernberg, Rev. Sci. Instrum. 82, 095113 (2011)]


    Input:  
            fname = file name (including path) of file to load.
                    If no argument is given, a dialoge is opened where 
                    the file can be selected.
            dataset = name of data set to load. If no dataset is given
                    a window is presented where the dataset can be selected.
            bins =  [optional] 3-tuple with number of elements to combine into one bin,
                    first number along x, second along y and third along z (time/energy)

    Output: 
            (Single Events):
                xscale = DLD pixel x-position of electron event
                yscale = DLD pixel y-position of electron event
                tEscale = detection time of electron event (flight time relative t_0)

            (Converted):
                xscale  = x-coordinate for electron event on detector in degrees (anglex)
                yscale  = y-coordinate for electron event on detector in degrees (angley)
                tEscale = kinetic energy of electron event (in eV)

    Usage: xscale,yscale,Ek (or t) = loadThemis()

    """
    # Check if the selected dataset is SingleEvents or Converted


    with h5.File(fname,'r') as f:
        
        # Determine the type of the data: single events or converted data
        for dset in f[dataset]:
            data_type = dset
        
        #data = np.array(f.get('{}/{}'.format(dataset,data_type)))
        #if data.size <= 1:
        #    data_type = 'conversion'
        #    data = np.array(f.get('{}/{}'.format(dataset,data_type)))
            

        if data_type == 'events':
            # Use standard binning if bins is not defined
            try:
                bins = kwargs['bins']
            except NameError:
                bins = (1024,1024,500)
    
            tfactor = f[dataset][data_type].attrs['FIELD_2_FACTOR']
            toffset = f[dataset][data_type].attrs['FIELD_2_OFFSET']
            data = np.array(f.get('{}/{}'.format(dataset,data_type)))
            df_data = pd.DataFrame(data)
            time = np.array(df_data['time'])*tfactor-toffset
            xscale = np.array(df_data.x)
            yscale = np.array(df_data.y)
            df_data_short = df_data.sample(20000)  # trying to bin all data crashes the kernel... out of memory? Therefore, load only part of data...

            # Get bin edges for binning in time
            #tmp,binedges = np.histogram(time,bins[2])
            # Get bin edges for binning in 
            H, edges = np.histogramdd(np.array(df_data_short),bins)

            #xscale = np.array(df_data.x)
            #yscale = np.array(df_data.y)
            
            #data,tE = np.histogram(time,bins[2])
            #tE = tE[:-1] + np.diff(tE)/2
        else: # the data is converted
            data = np.array(f.get('{}/{}'.format(dataset,data_type)))
            # Get attributes for selected dataset
            min_energy=f[dataset][data_type].attrs['minimumenergy']
            max_energy=f[dataset][data_type].attrs['maximumenergy']
            max_angle=f[dataset][data_type].attrs['maximumangle']
            # Create axis vectors
            xscale=np.linspace(-(max_angle/pi*180),(max_angle/pi*180),data.shape[0])
            yscale=np.linspace(-(max_angle/pi*180),(max_angle/pi*180),data.shape[1])
            tE=np.linspace(min_energy,max_energy,data.shape[2])
            
   
    return H, edges #xscale,yscale,tE,data