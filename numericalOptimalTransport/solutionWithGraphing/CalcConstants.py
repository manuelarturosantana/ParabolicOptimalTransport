import numpy as np
from math import sqrt

def calcDt(uVal, dx, del1):
    if (.5 * del1) < uVal:
        denom = .5 * del1
    else:
        denom = uVal
    big = (1/ (dx**2 * denom))
    return 1 / (2 * big)

def calcDx(max3rd, max4th, del1, maxdx):
    third = (3 * del1) / (2 * max3rd)
    fourth = sqrt((6 * del1) / max4th)
    if maxdx < third and maxdx < fourth:
        return maxdx
    elif third < fourth:
        return third
    else:
        return fourth

def deltacal(startX, EndX, C, D, fun, gun):
    fDomain = np.linspace(startX, EndX, 10000)
    u_xx = 1
    fVals = [fun(x) for x in fDomain]
    fMin, fMax = min(fVals), max(fVals)
    gDomain = np.linspace(C, D, 10000)
    gVals = [gun(x) for x in gDomain]
    gMin, gMax = min(gVals), max(gVals)
    delta1 = u_xx * (fMin/fMax)*(gMin/gMax)
    delta2 = u_xx *(fMax/fMin) *(gMax/gMin)
    print(delta1)
    print(delta2)
    gDerv = np.gradient(gVals)
    fDerv = np.gradient(fVals)
    gDmax, gDmin = max(abs(gDerv)), min(abs(gDerv))
    fDmax, fDmin = max(abs(fDerv)), min(abs(fDerv))
    gSec = np.gradient(gDerv)
    gSmax = max(abs(gSec))
    uTripM = 0 # This is only for our hard coded initial function.
    eVt = (u_xx * gMax) /fMin
    Fmax = fDmax / fMin
    Gmax = gDmax / gMin
    GDmax = (Gmax * gSmax + gDmax**2) * (1/ gMin**2)
    cont1 = ( Fmax + 2* delta2* Gmax + sqrt(Fmax + (2* delta2*Gmax)**2 + 4 * delta2**2 *GDmax))
    cont2 = (uTripM / u_xx) + Fmax + Gmax*u_xx
    if cont1 > cont2:
        phi_x = delta2 *( cont1 * eVt + Fmax + Gmax * delta2)
    else:
        phi_x = delta2 * (cont2 * eVt + Fmax + Gmax * delta2)
    return delta1, delta2, phi_x