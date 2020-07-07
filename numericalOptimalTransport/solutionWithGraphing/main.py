
from time import time








start = time()
solution, xVals, deltaTVals, rounds, graphing, tVals, secondD, ratios = numsolv(startX=-1, EndX=1, maxdx=.01, tolerance =0.01, C= -1, D=1,
                                                   fKey='littleStep', gKey='funPiece', functionDict=myFunctionDict)

end = time()
print(f"The total time was {end - start} \n\n")
print(solution)
print(f"The value of t is {tVals[len(graphing)-1]}")
print(rounds)