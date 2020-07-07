import numpy as np
from CalcConstants import calcDx, deltacal, calcDt
from calcRows import firstRow, calcOtherRow,
from ToleranceChecker import calcFInt, intTolChecker, intTolCheckerBetter

def numsolv(startX, EndX, maxdx, tolerance, C, D, fKey, gKey, functionDict):
    '''
    Takes in all the parts of the equation, the tolerance, and a maxdx if the stability theory allows it to be
    too large. Returns a dictionary containing the following listed by key
    uSolution: A list of the last approximation after the program ran to tolerance, and then to the next
    iterations % 1000 == 0
    xVals: A np.array of the x points used, good for graphing
    dtVals: The Values of dt as they changed over time.
    numOfStep: The number of iterations used to get within tolerance
    graphSol: A multidimensional list of the calculated solution every thousand timesteps.
    u_xx: A list of the approximated second derivative
    minu_xx: A list of the minimum approximated second derivatives at each time step.
    '''
    f = functionDict[fKey]
    g = functionDict[gKey]
    Gamma = 10
    Kappa = 10
    del1, del2, Psi = deltacal(startX, EndX, C, D, f, g)
    dx = calcDx(Psi, Gamma, del1, maxdx)
    numStep = int((EndX - startX) / dx)
    xVals = np.linspace(startX, EndX, numStep)
    dx = xVals[1]- xVals[0]
    graphSol = []
    dtVals = []
    uxxVals = []
    minUxxVals =[]
    totaltime = 0
    timestep = [0]
    rowbefore, minU_xx, uxxs = firstRow(xVals, C, D, dx)
    uxxVals.append(uxxs)
    graphSol.append(rowbefore)
    deltaT = calcDt(minU_xx, dx, del1)
    print(deltaT)
    currentRow, minU_xx, uxxs = calcOtherRow(xVals, rowbefore, dx, deltaT, C, D, f,g)  # This one is because of the way that SolT is written we need at least two rows to check it.
    iterations = 1
    totaltime = totaltime + deltaT
    INT_F_LIST = calcFInt(xVals, f)
    while intTolChecker(xVals, currentRow, tolerance, g, INT_F_LIST):
        deltaT = calcDt(minU_xx, dx, del1)
        rowbefore = currentRow
        currentRow, minU_xx, uxxs = calcOtherRow(xVals, rowbefore, dx, deltaT, C, D, f, g)
        totaltime += deltaT
        iterations += 1
        minUxxVals.append(minU_xx)
        if iterations % 1000 == 0:
            graphSol.append(currentRow)
            timestep.append(totaltime)
            dtVals.append(deltaT)
            uxxVals.append(uxxs)
    #Here we use a more accurate tolerance checker after its passed the "easier" one
    while intTolCheckerBetter(xVals, currentRow,tolerance, g, INT_F_LIST):
        deltaT = calcDt(minU_xx, dx, del1)
        rowbefore = currentRow
        currentRow, minU_xx, uxxs = calcOtherRow(xVals, rowbefore, dx, deltaT, C, D,f,g)
        totaltime += deltaT
        iterations += 1
        minUxxVals.append(minU_xx)
        if iterations % 1000 == 0:
            graphSol.append(currentRow)
            timestep.append(totaltime)
            dtVals.append(deltaT)
            uxxVals.append(uxxs)
    while iterations % 1000 != 0: #after within tolerance this pushes the code to run until we can add the solution to our graph
        deltaT = calcDt(minU_xx, dx, del1)
        rowbefore = currentRow
        currentRow, minU_xx, uxxs = calcOtherRow(xVals, rowbefore, dx, deltaT, C, D,f,g)
        totaltime += deltaT
        iterations += 1
        minUxxVals.append(minU_xx)
        if iterations % 1000 == 0:
            graphSol.append(currentRow)
            timestep.append(totaltime)
            dtVals.append(deltaT)
            uxxVals.append(uxxs)
    graphSol = np.asarray(graphSol)  # This is because the plotting only takes a numpy array
    solDict = {'uSolution' : currentRow, 'xVals': xVals, 'dtVals': dtVals, 'numOfStep': iterations, 'graphSol': graphSol,
               'timestep': timestep, 'u_xx': uxxVals, 'minu_xx': minUxxVals}
    return solDict