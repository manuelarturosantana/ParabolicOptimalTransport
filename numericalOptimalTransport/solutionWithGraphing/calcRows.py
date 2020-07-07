from math import log
from functions import u_0


def firstRow(xVals, C, D, dx):
    uSol = []
    uxxs = []
    for val in xVals:
        uSol.append(u_0(val))
    for i in range(len(uSol)):
        if i == 0:
            u_xx = 2 * (uSol[i + 1] - uSol[i] - C * dx)
            uxxs.append(u_xx/(dx**2))
            secU = u_xx
        elif i == len(uSol) - 1:
            u_xx = 2 * (uSol[i - 1] - uSol[i] + D * dx)
            uxxs.append(u_xx/(dx**2))
            if secU > u_xx:
                secU = u_xx
        else:
            u_xx = uSol[i - 1] + uSol[i + 1] - 2 * uSol[i]
            uxxs.append(u_xx/(dx**2))
            if secU > u_xx:
                secU = u_xx
    return uSol, (secU / dx**2), uxxs


def calcOtherRow(xVals, rowBefore, deltaX, deltaT, Cleft, Dright, foFX, goFX,):
    uSol = []
    uxxs = []
    minUVal = 0
    for index in range(len(xVals)):
        if index == 0:
            uVal, minUVal = (leftBoundary(rowBefore, index, deltaX, deltaT, xVals[0], Cleft,foFX, goFX))
            uSol.append(uVal)
            uxxs.append(minUVal)
        elif index == (len(xVals) - 1):
            uVal, u_xxVal =(rightBoundary(rowBefore, index, deltaX, deltaT, xVals[-1], Dright,foFX, goFX))
            uSol.append(uVal)
            uxxs.append(u_xxVal)
            if u_xxVal < minUVal:
                minUVal = u_xxVal
        else:
            uVal, u_xxVal = (discSol(rowBefore, index, xVals[index], deltaX, deltaT, foFX, goFX))
            uSol.append(uVal)
            uxxs.append(u_xxVal)
            if u_xxVal < minUVal:
                minUVal = u_xxVal
    return uSol, minUVal, uxxs

def discSol(uList, index, x, deltaX, deltaT, fOfX, gOfX):
    u_xx = (uList[index + 1] - 2 * uList[index] + uList[index - 1]) / (deltaX ** 2)
    u_x = (uList[index + 1] - uList[index - 1]) / (2 * deltaX)
    return (log(u_xx) - log(fOfX(x) / gOfX(u_x))) * deltaT + uList[index], u_xx


def leftBoundary(uList, index, deltaX, deltaT, x, C, fOfX, gOfX):  # A --> C
    u_xx = (2 * (uList[index + 1] - uList[index] - C * deltaX)) / (deltaX ** 2)
    return (log(u_xx) - log(fOfX(x) / gOfX(C))) * deltaT + uList[index], u_xx  # u'(x) = C at this point so we can replace that in g(x)


def rightBoundary(uList, index, deltaX, deltaT, x, D, fOfX, gOfX):  # B --> D
    u_xx = (2*(uList[index - 1] - uList[index] + D * deltaX)) / deltaX ** 2
    return (log(u_xx) - log(fOfX(x) / gOfX(D))) * deltaT + uList[index], u_xx