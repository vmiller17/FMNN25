class Spline

def __init__(self,uk,d) #Laubinot
""" 
    :param array.float uk: 
    :param array.float d:
    :raises TypeError: if uk is not an array of floats
    :raises TypeError: if d is not an array of floats
    :raises ValueError: if relation between uk and d is wrong

"""
    
def __call__(self,u) #Eli  
""" 
    :param array.float u: array with all point where one wants to evaluate the spline
    :raises TypeError: if u is not an array of floats
    :raises ValueError: if any value in u is outside boundaries of uk
    :returns array.dim(2xn).floats Val: x and y values for the spline. x values in
                                    the first row and y values in the second.
                                     
"""

def _eval(self,u) #Eli
"""
    :param float u: point where one wants to evaluate the spline
    :raises ValueError: if u is outside boundaries of uk
    :raises TypeError: if u is not a float
    :returns array.dim(2x1).floats: x and y values for the spline. x values in
                                    the first row and y values in the second.
    
"""

def _findHotIntervall(self,u) #Laubinot
"""
    :param float u: point where one wants to find the hot intervall
    :raises ValueError: if u is outside boundaries of uk
    :raises TypeError: if u is not a float
    :returns int j: index for the left side of the intervall
    
"""

def plot(self,uStart=uk[2],uStop=uk[-3]-kinds.default_float_kind.MIN*100,numPoints=1000) #FeelFree
"""
    :param float uStart: point where ploting starts, default uk[2] 
    :param float uStop: point where ploting ends, default uk[2]
    :param int numPoints: number of point where one wants to evaluate the spline
    :raises ValueError: if any value in u is outside boundaries of uk
    :raises TypeError: if uStart or uStop is not a float
    :returns array.float u: array with all point where the spline has been evaluated
    :returns array.dim(2xn).floats val: x and y values for the spline. x values in
                                    the first row and y values in the second.
    
"""

def getBaseFunc(self,j) #Hanna

"""
    :param int j: index for base function, negative values are accepted
    :returns function N: Base function
    :raises TypeError: if j is not a integer
    :raises ValueError: if j is outside uk
"""