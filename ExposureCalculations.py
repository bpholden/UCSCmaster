import numpy as np
from x_gaussslit import *
    
def getI2_K(unc):


    A = 4.47
    B = -1.58
    return 10 ** ( A + B*np.log10(unc))
    
def getI2_M(unc):


    A = 4.14
    B = -1.73
    return 10 ** ( A + B*np.log10(unc))
    
def getEXPMeter(i2, bv):
    delta = 4.52
    epsilon = -0.196
    squiggle = 0.262
    x = delta + epsilon*bv + squiggle*bv**2
    x += 0.05
    return i2 * 10**x

def getEXPMeter_Rate(v, bv, el, seeing, decker="W"):
    alpha = -0.908
    beta = 0.0852
    Const = -21.8
    if seeing == 0:
        apflog( "Warning: AVG_FWHM seems to be 0. Using 15 instead.",level="Warn")

        
        seeing = 15
    # seeing  = 13.99
    light = x_gaussslit(slit_size[decker][0]/seeing, slit_size[decker][1]/seeing, 0, 0)
    # light = 0.442272
    
    VC = v - 2.5*np.log10(light)
    x = (-1/2.5) * (VC + alpha*bv + beta*(1/np.cos(np.radians(90-el))) + Const)
    return (10 ** x)


def getEXPTime(cnts, v, bv, el, seeing, decker="W"):
    alpha = -0.0311
    beta = 0.158
    Const = -11.2
    if seeing <= 0:
        seeing = np.array(15)
    # seeing  = 13.99
    light = x_gaussslit(slit_size[decker][0]/seeing, slit_size[decker][1]/seeing, 0, 0)
    # light = 0.442272
    
    #if el < 15.0:
     #   el = 15.0
        # bogus exposure time but the APF does not work this low anyway
        
    if len(bv) != len(el): print "Error: getEXPTime arrays don't match"
    VC = v - 2.5*np.log10(light)
    x = (-1/2.5) * (VC + alpha*bv + beta*(1/np.cos(np.radians(90-el))) + Const)
    cnt_rate = 10**x
    time = 0
    time = cnts/cnt_rate
    return time
