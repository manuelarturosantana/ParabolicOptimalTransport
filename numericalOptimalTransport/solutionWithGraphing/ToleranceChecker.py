from scipy.integrate import quad

def calcFInt(xValues, f):
    fList = []
    for i in range(len(xValues)):
        fVal, err = quad(f,-1, xValues[i])
        fList.append(fVal)
    return fList

def intTolChecker(xValues, uSol, tolerance, gFunc, fList):
    estUx = [-1]
    dx = xValues[1] - xValues[0]
    for i in range(len(uSol) - 1):
        if i != 0:
            u_x = (uSol[i + 1] - uSol[i - 1]) / (2 * dx)
            estUx.append(u_x)
    estUx.append(1)
    gList = []
    for i in range(len(xValues)):
        if i % 10 == 0:
            gVal, err = quad(gFunc, -1, estUx[i])
            gList.append(gVal)
    maxError = 0
    #Currently the F list is the full thing will the g List is just the small thing
    for value in range(len(gList)):
        error = abs(fList[value * 10] - gList[value])
        if error > maxError:
            maxError = error
    return maxError > tolerance


def intTolCheckerBetter(xValues, uSol, tolerance, gFunc, fList):
    estUx = [-1]
    dx = xValues[1] - xValues[0]
    for i in range(1, len(uSol) - 1):
        if i != 0:
            u_x = (uSol[i + 1] - uSol[i - 1]) / (2 * dx)
            estUx.append(u_x)
    estUx.append(1)
    gList = []
    for i in range(len(xValues)):
        gVal, err = quad(gFunc, -1, estUx[i])
        gList.append(gVal)
    maxError = 0
    for value in range(len(fList)):
        error = abs(fList[value] - gList[value])
        if error > maxError:
            maxError = error
    return maxError > tolerance