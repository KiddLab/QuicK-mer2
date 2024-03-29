import numpy
from numpy import median

def lowess(x, y, f=2./3., iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest

Lowess smoother: Robust locally weighted regression.
The lowess function fits a nonparametric regression curve to a scatterplot.
The arrays x and y contain an equal number of elements; each pair
(x[i], y[i]) defines a data point in the scatterplot. The function returns
the estimated (smooth) values of y.

The smoothing span is given by f. A larger value for f will result in a
smoother curve. The number of robustifying iterations is given by iter. The
function will run faster with a smaller number of iterations."""
    n = len(x)
    r = int(numpy.ceil(f*n))
    h = [numpy.sort(numpy.abs(x-x[i]))[r] for i in range(n)]
    w = numpy.clip(numpy.abs(([x]-numpy.transpose([x]))/h),0.0,1.0)
    w = 1-w*w*w
    w = w*w*w
    yest = numpy.zeros(n)
    delta = numpy.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:,i]
            b = numpy.array([sum(weights*y), sum(weights*y*x)])
            A = numpy.array([[sum(weights), sum(weights*x)],
                             [sum(weights*x), sum(weights*x*x)]])

#            beta = numpy.linalg.solve(A,b)            
            # changed to deal with rare cases with a singular matrix error
            beta = numpy.linalg.lstsq(A,b,rcond=-1)[0]
               
            
            yest[i] = beta[0] + beta[1]*x[i]
        residuals = y-yest
        s = numpy.median(abs(residuals))
        delta = numpy.clip(residuals/(6*s),-1,1)
        delta = 1-delta*delta
        delta = delta*delta
    return yest
