import numpy as np
from scipy.integrate import quad
from math import log, sqrt
def u_0(x): #Our initial function.
    return .5 * ((x) ** 2)

def f(x):
    #This function defines f(x) from -1 to 1
    return (1 / (exp(1) - exp(-1))) * exp(x)

def g(x):
    #This function define g(x) from 2 to 4
    return (1 / 2) * (x - 3) ** 2 + (1 / 3)

#Here we have a sample dictionary for the left and right bounds for f and g.
boundDict = {'A' : -1, 'B' : 1, 'C' : 2, 'D' : 4}

def discSol(uList, index, x, deltaX, deltaT, fOfX, gOfX):
    '''
    Takes in a list from the solution of the previous row, and returns a
    value for the point x_j at the next time step.
    '''
    u_xx = (uList[index + 1] - 2 * uList[index] + uList[index - 1]) / (deltaX ** 2)
    u_x = (uList[index + 1] - uList[index - 1]) / (2 * deltaX)
    return (log(u_xx) - log(fOfX(x) / gOfX(u_x))) * deltaT + uList[index], u_xx


def leftBoundary(uList, index, deltaX, deltaT, x, boundDict, fOfX, gOfX):
    '''
    Takes in a list from the previous time step and using our approximation at the left
    boundary returns a value for the left boundary condition.
    '''
    u_xx = (2 * (uList[index + 1] - uList[index] - boundDict['C'] * deltaX)) / (deltaX ** 2)
    # u'(x) = C at this point so we can replace that in g(x)
    return (log(u_xx) - log(fOfX(x) / gOfX(boundDict['C']))) * deltaT + uList[index], u_xx


def rightBoundary(uList, index, deltaX, deltaT, x, D, fOfX, gOfX):
    '''
        Takes in a list from the previous time step and using our approximation at the right
        boundary returns a value for the right boundary condition.
    '''
    u_xx = (2*(uList[index - 1] - uList[index] + D * deltaX)) / deltaX ** 2
    # u'(x) = D at this point so we can replace that in g(x)
    return (log(u_xx) - log(fOfX(x) / gOfX(D))) * deltaT + uList[index], u_xx


def firstRow(xVals, boundDict, dx):
    '''
    Uses the initial condition to calculate the first row of the finite difference scheme.
    Additionally the minimun value of \Delta U^n_j is calculated so it be used to calculate \Delta t
    '''
    uSol = []
    for val in xVals: #create the first row based of the initial condition.
        uSol.append(u_0(val))
    for i in range(len(uSol)):
        if i == 0: #left boundary condition estimate for the second derivative
            u_xx = 2 * (uSol[i + 1] - uSol[i] - boundDict['A'] * dx)
            secU = u_xx
        elif i == len(uSol) - 1: #right boundary condition estimate for the second derivative
            u_xx = 2 * (uSol[i - 1] - uSol[i] + boundDict['B'] * dx)
            if secU > u_xx:
                secU = u_xx
        else: #interitor bounds estimate for the second derivative.
            u_xx = uSol[i - 1] + uSol[i + 1] - 2 * uSol[i]
            if secU > u_xx:
                secU = u_xx
    return uSol, (secU / dx**2) #Returns the first row and the minimun second derivative approximation in python.


def calcOtherRow(xVals, rowBefore, deltaX, deltaT, boundDict, foFX, goFX,):
    '''
    Uses the finite difference scheme and the row before to calculate the next row.
    Additionally it returns the smallest second derivative approximation to calculate the next time step.
    '''
    uSol = [] #List to contain the approximation of the next row.
    minUVal = 0
    for index in range(len(xVals)):
        if index == 0:#left boundary condition estimate for the second derivative
            uVal, minUVal = leftBoundary(rowBefore, index, deltaX, deltaT, xVals[0], boundDict['C'],foFX, goFX)
            uSol.append(uVal)
        elif index == (len(xVals) - 1): #right boundary condition estimate for the second derivative
            uVal, u_xxVal = rightBoundary(rowBefore, index, deltaX, deltaT, xVals[-1], boundDict['D'],foFX, goFX)
            uSol.append(uVal)
            if u_xxVal < minUVal:
                minUVal = u_xxVal
        else: #interior bounds estimate for the second derivative 
            uVal, u_xxVal =  discSol(rowBefore, index, xVals[index], deltaX, deltaT, foFX, goFX)
            uSol.append(uVal)
            if u_xxVal < minUVal:
                minUVal = u_xxVal
    return uSol, minUVal


def calcFInt(xValues, boundDict, f):
    '''
    Calculates the integral of f at each boundary point to be used in the tolerance checker.
    '''
    fList = []
    for i in range(len(xValues)):
        fVal, err = quad(f, boundDict['A'], xValues[i])
        fList.append(fVal)
    return fList

def intTolChecker(xValues, uSol, tolerance, boundDict, gFunc, fList):
    '''
    Compares the integrals F(x) to G(u'(x)) using a subset of the grid points to see if the
    difference of the two integrals is within tolerance.
    '''
    estUx = [boundDict['C']] # The boundary condition u'(A) = C
    dx = xValues[1] - xValues[0]
    for i in range(len(uSol) - 1): #This loop estimates u'(x) using our finite difference approximation.
        if i != 0:
            u_x = (uSol[i + 1] - uSol[i - 1]) / (2 * dx)
            estUx.append(u_x)
    estUx.append(boundDict['D']) #The boundary condition u'(B) = D
    gList = [] # This collects the CDF of g
    for i in range(len(xValues)):
        if i % 10 == 0: #This uses the subset of the grid points equally spaced to check the tolerance.
            gVal, err = quad(gFunc, -1, estUx[i])
            gList.append(gVal)
    maxError = 0
    for value in range(len(gList)):
        error = abs(fList[value * 10] - gList[value]) #The indexing on fList is to account for the spacing of G
        if error > maxError:
            maxError = error
    return maxError > tolerance


def intTolCheckerBetter(xValues, uSol, boundDict, tolerance, gFunc, fList):
    '''
    Compares the integrals F(x) to G(u'(x)) using all grid points to see if the
    difference of the two integrals is within tolerance.
    '''
    estUx = [boundDict['C']] # The boundary condition u'(A) = C
    dx = xValues[1] - xValues[0]
    for i in range(1, len(uSol) - 1): #This loop estimate u'(x) using our finite difference approximation.
        if i != 0:
            u_x = (uSol[i + 1] - uSol[i - 1]) / (2 * dx)
            estUx.append(u_x)
    estUx.append(boundDict['D'])  #The boundary condition u'(B) = D
    gList = []
    for i in range(len(xValues)): #Calculates the G(u'(x)) on every grid point
        gVal, err = quad(gFunc, -1, estUx[i])
        gList.append(gVal)
    maxError = 0
    for value in range(len(fList)): # Checks to see if the error is within tolerance
        error = abs(fList[value] - gList[value])
        if error > maxError:
            maxError = error
    return maxError > tolerance

def calcDt(uVal, dx, del1):
    '''
    from the stability analysis calculate Delta t given the smallest second spacial derivative
    approximation, size of dx, and the lower bound on the second derivative.
    '''
    if (.5 * del1) < uVal: #compares the lower bound on the second derivative to the approxmation
        minVal = .5 * del1
    else:
        minVal = uVal
    return (1 / 2) * dx ** 2 * minVal #returns delta T

def calcDx(max3rd, max4th, del1, maxdx):
    '''
    Uses the maximum of the 3rd and 4th spacial derivative, and the lower bound on the second spacial derivative
    to calculate the largest value dx can be. If this value is too large to produce good accuracy a limiting value
    is passed in to assure dx stays low enough.
    '''
    third = (3 * del1) / (2 * max3rd)
    fourth = sqrt((6 * del1) / max4th)
    if maxdx < third and maxdx < fourth: #Here it returns the smallest dx for the best resolution on the solution
        return maxdx
    elif third < fourth:
        return third
    else:
        return fourth

def deltaPhiCal(boundDict, u_02Der, u_03Der, fun, gun):
    '''
    Calculates the upper and lower bound on second spacial derivative, (delta1 and delta2) and the lower bound
    on the third spacial derivative (psi). Assumes the second and third derivative of the initial condition can be worked
    out analytically and passed into the function.
    '''
    fDomain = np.linspace(boundDict['A'], boundDict['B'], 10000)
    fVals = [fun(x) for x in fDomain] #Finds the minimum and maximum values of f
    fMin, fMax = min(fVals), max(fVals)
    gDomain = np.linspace(boundDict['C'], boundDict['D'], 10000)
    gVals = [gun(x) for x in gDomain] #Finds the minimum and maximum values of g
    gMin, gMax = min(gVals), max(gVals)
    delta1 = u_02Der * (fMin/fMax) * (gMin/gMax) #calculates delta1, and delta2 from the stability analysis
    delta2 = u_02Der *(fMax/fMin) * (gMax/gMin)
    gDerv = np.gradient(gVals) #calculates the derivative of f and g to be used in calculating the third derivative.
    fDerv = np.gradient(fVals)
    gDmax, gDmin = max(abs(gDerv)), min(abs(gDerv))
    fDmax, fDmin = max(abs(fDerv)), min(abs(fDerv))
    gSec = np.gradient(gDerv) # Calculates the maximum second derivate of g.
    gSmax = max(abs(gSec))
    eVt = (u_02Der * gMax) /fMin # Calculates e raised to the first time derivative of the parabolic problem.
    Fmax = fDmax / fMin # Here F and G are logarithmic deriviates used to calculate the third spacial derivative.
    Gmax = gDmax / gMin
    GDmax = (Gmax * gSmax + gDmax ** 2) * (1/ gMin ** 2) # A calculation for G'
    #These two constants come from the stability analysis to be compared when calculating the third spcial
    #derivative.
    cont1 = ( Fmax + 2 * delta2 * Gmax + sqrt(Fmax + (2 * delta2 * Gmax) ** 2 + 4 * delta2 ** 2 * GDmax))
    cont2 = (u_03Der / u_02Der) + Fmax + Gmax * u_02Der
    if cont1 > cont2:
        psi = delta2 *( cont1 * eVt + Fmax + Gmax * delta2) #After the constants are compared, psi is calculated.
    else:
        psi  = delta2 * (cont2 * eVt + Fmax + Gmax * delta2)
    return delta1, delta2, psi

def numsolv(maxdx, tolerance, boundDict, fKey, gKey, functionDict):
    '''
    This function is used to iterate through the finite difference approximation up to a tolerance. It takes in
    the tolerance, an upper bound for dx, a dictionary containing the four boundary points. It returns the final approximations which can then be used to
    estimate u_x and u_xx
    '''
    f = functionDict[fKey] #initializing the functions for f and g
    g = functionDict[gKey]
    Gamma = 10 # Here we abitrarily choose a value for the upper bound on the fourth derivative.
    del1, del2, Psi = deltaPhiCal(boundDict, f, g)
    dx = calcDx(Psi, Gamma, del1, maxdx)
    numStep = int((boundDict['A'] - boundDict['B']) / dx) #Here we make sure dx gives us an integer number of grid points.
    xVals = np.linspace(boundDict['A'], boundDict['B'], numStep)
    dx = xVals[1]- xVals[0] #The new delta x based on the finite number of grid points
    rowbefore, minU_xx = firstRow(xVals, boundDict, dx)
    intFList = calcFInt(xVals, boundDict, f)
    #Checking for tolerence based on a subset of the grid points.
    while intTolChecker(xVals, currentRow, tolerance, g, intFList):
        deltaT = calcDt(minU_xx, dx, del1)
        rowbefore = currentRow
        currentRow, minU_xx = calcOtherRow(xVals, rowbefore, dx, deltaT, boundDict, f, g)
    #Checking for tolerance based on all the grid points.
    while intTolCheckerBetter(xVals, currentRow,tolerance, g, intFList):
        deltaT = calcDt(minU_xx, dx, del1)
        rowbefore = currentRow
        currentRow, minU_xx = calcOtherRow(xVals, rowbefore, dx, deltaT, boundDict, f,g)
    return currentRow

