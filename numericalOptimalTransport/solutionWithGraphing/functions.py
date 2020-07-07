from math import exp, log, sin, cos
def exponential(x):
    return (1 / (exp(1) - exp(-1))) * exp(x) #\frac{1}{e^1 - e^{-1}}e^x

def oneHalfXSquared(x):
    return (1 / 2) * x ** 2 + (1/3)

def oneFourthX(x): #\frac{1}{4} (x+2)
    return (1/4)* (x+2)

def cosine(x):
    return (1/(2*sin(1)))*cos(x)

def littleStep(x):
    if x <= 0:
        return 0.75
    else:
        return 0.25

def bigStep(x):
    if x <= -0.5:
        return 0.3
    elif x <= 0:
        return 0.6
    elif x <= 0.5:
        return 0.2
    else:
        return 0.9

def constant(x): #\frac{1}{2}
    return .5

def funPiece(x):
    if x <= 0:
        return (1/(2 - 2 * exp(-1))) * exp(x)
    else:
        return (3/14) * (x - 2) ** 2
def logarithmic(x):
    return (log(x+2))/(3*log(3)+2)+2/(3*log(3)+2) #\frac{\log(x+2)}{3\log(3)+ 2} + \frac{2}{3\log(3) + 2}

def intFunPiece(x):
    if x <= 0:
        return (1/(2 - 2 * exp(-1))) * exp(x)
    else:
        return (1/14) * (x-2) ** 3 + intFunPiece(0) - intFunPiece(-1)

def didntwork(x):
    return (9/20) * x + (1/2)

def aLittleChange(x):
    return (1/32)*x + (1/2)

def hiFreq(x):
    alpha = 50 / (sin(100) + 200) #\frac{50}{\sin(100)+ 200}\cos(100x)+2
    return alpha * (cos(100 * x) + 2)

def testa(x):
    return 1/(2 * (x +((exp(2)+1)/ (exp(2)-1))))

def testb(x):
    return 1/(-2 *(x - ((exp(2)+1)/(exp(2) - 1))))


myFunctionDict = {
    'exponential' : exponential,
    'oneHalfXSquared' : oneHalfXSquared,
    'oneFourthX' : oneFourthX,
    'cos': cosine,
    'littleStep': littleStep,
    'bigStep' : bigStep,
    'constant': constant,
    'log' : logarithmic,
    'funPiece' : funPiece,
    'hiFreq': hiFreq,
    'a': testa,
    'b': testb,
    'nine20th': didntwork,
    'stepEnd': aLittleChange
}

def u_0(x):
    return .5 * ((x) ** 2)