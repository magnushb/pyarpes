import numpy as np
from scipy.constants import pi
from scipy.constants import m_e
from scipy.constants import hbar
from scipy.constants import e
import math as math

def angle2k(angle,Ekin):
    """
    Calculates k-value from angle and kinetic energy.

    Input:  
            angle = single value, vector or matrix (np.array) with angles in degrees        
            Ekin = kinetic energy in eV (np.array)
            
    Output: 
            k = kvalues (np.array) in inverse Ångstrom

    Usage:  k = angletok(angle,Ekin)
    """

    k = np.multiply(np.sqrt(2*m_e*Ekin*e),np.sin(angle*pi/180)/hbar*1e-10)

    return k

def cart2pol(x,y):
        """
        Converts cartesian coordinates (x,y) to polar (rho, theta).
        If the coordinates cover the origo of the cartesian coordinate system, 
        the angle of this point in polar coordinates is set to zero.

        Input: 
                x = numpy array with x coordinates
                y = numpy array with y coordinates

        Output:
                theta = polar angles (in degrees) for the (x,y) pair of the input
                rho = radial distance to (x,y) point from origo

        Usage:  theta,rho = cart2pol(x,y)
        """

        theta = np.arctan2(y,x)*180/pi
        rho = np.hypot(x,y)

        return theta, rho


def pol2cart(theta,rho):
        """
        Converts polar coordinates (theta,rho) to cartesian (x, y).

        Input: 
                theta = numpy array with the polar angle (in degrees) to the point
                rho = numpy array with the radial distance to the point

        Output:
                x = x coordinate for the point
                y = y coordinate for the point

        Usage:  x,y = pol2cart(theta,rho)
        """

        x = np.multiply(rho,np.cos(theta/180*pi))
        y = np.multiply(rho,np.sin(theta/180*pi))

        return x, y


def angleangle2kk(X,Y,Ekin):
        """
        Converts an 2D angle meshgrid to k-k grid. Useful for plotting constant energy slices for 2D systems.
        Note that if you have a kz dispersion (3D system), the measured angle-angle-Energy slice does not correspond to a
        constant kz value.

        Input:
                x_angle = 2D numpy array with x-angle values (meshgrid-style)
                y_angle = 2D numpy array with y-angle values (meshgrid-style)
                Ekin = kinetic energy in eV for the angle slice (single number or meshgrid-style with one value)

        Output:
                kx = k_x values (meshgrid-style)
                ky = k_y values (meshgrid-style)

        Usage:  kx,ky = anglegrid2kk(x_angle,y_angle,Ekin)
        """

        th,rho = cart2pol(X,Y)
        k = angle2k(rho,Ekin)
        kx,ky = pol2cart(th,k)

        return kx,ky


def angleE2kE(angle,Ek):
        """
        Converts an angle vector and a corresponding energy vector into a k-matrix and a corresponding kinetic energy matrix.

        This uses the angle and energy vectors that you get from a measured slice using a hemispherical analyzer
        and return a 2D k-matrix and Ek-matrix that can be used to plot the intensity in the slice as a k-Ek slice.

        Input:
                angle = numpy array with angle values (1D vector) in degrees
                Ek = numpy array with kinetic energies (1D vector) in eV

        Output: 
                k = 2D array with k-values, each row corresponding to the angle vector converted using one single Ek value
                E = 2D array with energies in eV, where each column is identical

        Usage:  k,energy = angle2kgrid(angle,Ek)
        """
        k,energy = np.meshgrid(angle,Ek)
        k = k.astype('float64')
        
        for u in range(Ek.size):
                k[u,:]=angle2k(angle,Ek[u])

        return k,energy


def Eph2Eftof(Eph,workfunction):
        """
        Function for calculating the Fermi-level electron time-of-flight in DriftMode
        of the Themis1000 analyzer based on the photon energy and the workfunction of
        the sample.

        The sample-to-detector distance is fixed at 880.5 mm. 
        Flight time is in nano seconds (same unit as Themis acquisition software).
        
        Input:
                Eph = photon energy in eV
                workfunction = the workfunction of the sample, in eV

        Output:
                t = time-of-flight for the electrons at the Fermi level, in nanoseconds
    
        Usage: t = Eph2Eftof(10.8,4.5)
    
        """
    
        Ekin = Eph - workfunction
        v = np.sqrt(2*Ekin*e/m_e)
        t = 880.5e-3/v/1e-9
        
        return t


def tof2Ek(tof):
        """
        Function for calculating the kinetic energy of electrons 
        based on the measured time-of-flight through the Themis1000 analyzer
        in DriftMode.

        Input: 
                tof = time-of-flight, in nanoseconds
        
        Output:
                Ek = kinetic energy of electron, in eV

        Usage: Ek = tof2Ek(tof)
        """

        v = 880.5e-3/tof/1e-9
        Ek = 0.5*m_e*v**2/e
        
        return Ek


def Ek2tof(Ek):
        """
        Function for calculating the time-of-flight through the Themis1000 
        analyzer in DriftMode for electrons with a given kinetic energy. 

        Input: 
               Ek = kinetic energy of electron, in eV 
        
        Output:
                tof = time-of-flight, in nanoseconds

        Usage: Ek = Ek2tof(Ek)
        """

        v = np.sqrt(2*Ek*e/m_e)
        t = 880.5e-3/v/1e-9

        return t