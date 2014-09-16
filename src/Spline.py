import numpy as np
import sys
import matplotlib.pyplot as plt
class Spline(object):
    """ 
    :param array uk: The knot vector containing all knot points :math:`\{u_i\}`. It must be an array of dimension (n,) which contains floats.
    :param array d: The control point vector containing all control points :math:`\{\\textbf{d_i}\}` . It must be an array of dimension (2,n-2) which contains floats.
    :raises TypeError: if **uk** is not an array of floats
    :raises TypeError: if **d** is not an array of floats
    :raises ValueError: if relation between **uk** and **d** is wrong
    :raises ValueError: if there are not enough elements (at least 6) in **uk**
    """

    def __init__(self,uk,d):

        if not isinstance(uk,np.ndarray):  
            raise TypeError("uk must be a numpy array")
        if not issubclass(uk.dtype.type,float):
            raise TypeError("uk must be an np.array of floats")
        if not issubclass(d.dtype.type,float):
            raise TypeError("d must be an np.array of floats")       
        if len(uk) != d.shape[1] + 2:
            raise ValueError("relation between uk and d is wrong")
        if uk.size < 7:
            raise ValueError("uk must at least contain 6 elements")
        if not (uk == sorted(uk)).all():
            raise ValueError("uk must be a sorted np.array")
    
        self.d = np.array(d,dtype='float')
        self.uk = np.array(uk,dtype='float')
    
    def __call__(self,u):
        """
        :param np.array.float u: np.array with all point where one wants to evaluate the spline
        :raises TypeError: if u is not an np.array of floats
        :raises ValueError: if any value in u is outside boundaries of uk
        :returns np.array.dim(2xn).floats Val: x and y values for the spline. x values in
                                    the first row and y values in the second.
        """
        if not isinstance(u,np.ndarray):  
            raise TypeError("uk must be a numpy array")
        if not issubclass(u.dtype.type,float):
            raise TypeError("uk must contain floats")
        if max(u) > self.uk[-2] or min(u) < self.uk[2]:
            raise ValueError('u is out of bounds.')
         
        n = len(u)
        Val = np.zeros((2, n))
        for i in range(n):
            tmp = self._eval(u[i])
            Val[0, i] = tmp[0]
            Val[1, i] = tmp[1]
        return Val
 
    def _eval(self,u):
        """
        :param float u: point where one wants to evaluate the spline
        :raises ValueError: if u is outside boundaries of uk
        :raises TypeError: if u is not a float
        :returns np.array.dim(2x1).floats: x and y values for the spline. x values in
                                        the first row and y values in the second.
        """
        if not issubclass(type(u),float):
            raise TypeError('u must be of type float')
        if u > self.uk[-2] or u < self.uk[2]:
            raise ValueError('u is out of bounds.')

        
        index=self._findHotInterval(u)
        #alpha = (uk[index] - u) / (uk[index] - uk[index+1])
        #value = alpha * self._findD(u,index,index) + (1-alpha) * self._findD(u,index+1,index+1)
        xValue = self._findD(u,index+1,index,0)
        yValue = self._findD(u,index+1,index,1)        

        return np.array([[xValue], [yValue]])
        
    def _findD(self,u,leftMost,rightMost,coord):
        if rightMost - leftMost == 2:
            return self.d[coord,leftMost]
        alpha = (self.uk[rightMost+1] - u)/(self.uk[rightMost+1] - self.uk[leftMost-1])
        return alpha * self._findD(u,leftMost-1,rightMost,coord) + (1 - alpha) * self._findD(u,leftMost,rightMost+1,coord)
            
    def _findHotInterval(self,u):
        """
        :param float u: point where one wants to find the hot intervall
        :raises ValueError: if u is outside boundaries of uk
        :raises TypeError: if u is not a float
        :returns int j: index for the left side of the intervall
        """

        if not isinstance(u,float):
            raise TypeError("u must be a float")            
        if (u < self.uk[0]).any() or (u >= self.uk[-1]).any():
            raise ValueError("u is outside boundaries of uk")
        j = np.searchsorted(self.uk[1:],u)
            
        return int(j)

    def plot(self,uStart=None,uStop=None,numPoints=1000):
        """
        :param float uStart: Point where plotting starts. If none are chosen it starts in the first point where the spline is defined.
        :param float uStop: Point where plotting ends. If none are chosen it stops in the last point where the spline is defined.
        :param int numPoints: number of point where one wants to evaluate the spline
        :raises ValueError: if any value in u is outside boundaries of uk
        :raises TypeError: if uStart or uStop is not a float
        :returns np.array.float u: np.array with all point where the spline has been evaluated
        :returns np.array.dim(2xn).floats val: x and y values for the spline. x values in the first row and y values in the second.
        """

        if (uStart == None): uStart = float(self.uk[2])
        if (uStop == None): uStop = float(self.uk[-3]-sys.float_info.epsilon)

        if (type(uStart) != float): raise TypeError('uStart is not a float')
        if (type(uStop) != float): raise TypeError('uStop is not a float')

        u_values = np.linspace(uStart,uStop,numPoints)
        try:
            coords = self.__call__(u_values)
        except ValueError, e:
            raise e

        plt.plot(coords[0,:],coords[1,:])
        plt.show()

        return coords

    def getBaseFunc(self,j):

        """
        :param int j: index for base function, negative values are accepted
        :returns function N: Base function
        :raises TypeError: if j is not a integer
        :raises ValueError: if j is outside uk
        """
        uk=self.uk

        if type(j) != int:
            raise TypeError('j must be an integer value')

        if j > len(uk)-2:
            raise ValueError('j is out of bounds')

        def N(u):
            try:
                index = self._findHotInterval(u)
            except ValueError or TypeError:
                raise

            if (index -j > 4) or (index < j):
                nVal=0
            elif (index - j == 3):
                try :
                    nVal=(uk[j+4]-u)**3/((uk[j+4]-uk[j+1])*(uk[j+4]-uk[j+2])*(uk[j+4]-uk[j+3]))
                except ZeroDivisionError:
                    nVal = 0
            elif (index - j == 2):
                try:
                    nVal = (u-uk[j])*(uk[j+3]-u)**2/((uk[j+3]-uk[j])*(uk[j+3]-uk[j+1])*(uk[j+3]-uk[j+2]))
                    nVal += (uk[j+4] - u)*(u-uk[j+1])*(uk[j+3]-u)/((uk[j+4]-uk[j+1])*(uk[j+3]-uk[j+1])*(uk[j+3]-uk[j+2]))
                    nVal += (uk[j+4]-u)**2*(u-uk[j+2])/((uk[j+4]-uk[j+1])*(uk[j+4]-uk[j+2])*(uk[j+3]-uk[j+2]))
                except ZeroDivisionError:
                    nVal = 0
            elif (index - j == 1):
                try:
                    nVal = (u-uk[j])**2*(uk[j+2]-u)/((uk[j+3]-uk[j])*(uk[j+2]-uk[j])*(uk[j+2]-uk[j+1]))
                    nVal += (u-uk[j])*(uk[j+3]-u)*(u-uk[j+1])/((uk[j+3]-uk[j])*(uk[j+3]-uk[j+1])*(uk[j+2]-uk[j+1]))
                    nVal += (uk[j+4]-u)*(u-uk[j+1])**2/((uk[j+4]-uk[j+1])*(uk[j+3]-uk[j+1])*(uk[j+2]-uk[j+1]))
                except ZeroDivisionError:
                    nVal = 0
            elif (index == j):
                try:
                    nVal=(u-uk[j])**3/((uk[j+3]-uk[j])*(uk[j+2]-uk[j])*(uk[j+1]-uk[j]))
                except ZeroDivisionError:
                    nVal = 0
            return nVal

        return N
