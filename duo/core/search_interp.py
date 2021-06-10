import numpy as np
import duo.core.duo_exception as de

#------------------------------------------------------------
#------------------------------------------------------------
def LinLin(x0, y0, x1, y1, m):
    return (y1 - y0) / (x1 - x0) * (m - x0) + y0

#------------------------------------------------------------
#------------------------------------------------------------
def LogLog(x0, y0, x1, y1, m):
    m = np.log(m)

    x0 = np.log(x0)
    y0 = np.log(y0)

    x1 = np.log(x1)
    y1 = np.log(y1)

    result = (y1 - y0) / (x1 - x0) * (m - x0) + y0
    result = np.exp(result)

    return result

#------------------------------------------------------------
#------------------------------------------------------------
def InterpXSLinearSearch(xsList, energy):
    if energy < xsList[0].energy:
        return 0.0
    elif energy > xsList[-1].energy:
        return 0.0
    else:
        for i in range(len(xsList)):
            if energy == xsList[i].energy:
                return xsList[i].microXS

            if energy < xsList[i].energy:
                if xsList[i - 1].microXS == 0.0 or xsList[i].microXS == 0.0:
                    # use lin-lin interpolation
                    return LinLin(xsList[i - 1].energy, xsList[i - 1].microXS, xsList[i].energy, xsList[i].microXS, energy)
                else:
                    # still use lin-lin interpolation
                    return LinLin(xsList[i - 1].energy, xsList[i - 1].microXS, xsList[i].energy, xsList[i].microXS, energy)

#------------------------------------------------------------
#------------------------------------------------------------
def InterpXSBinarySearch(xsList, energy):
    if energy < xsList[0].energy:
        return 0.0
    elif energy > xsList[-1].energy:
        return 0.0
    else:
        low = 0
        high = len(xsList) - 1

        while high - low > 1:
            mid = (low + high) // 2
            if energy < xsList[mid].energy:
                high = mid
            elif energy > xsList[mid].energy:
                low = mid
            else: # energy == xsList[mid].energy:
                return xsList[mid].microXS

        # use lin-lin interpolation
        return LinLin(xsList[low].energy, xsList[low].microXS, xsList[high].energy, xsList[high].microXS, energy)



