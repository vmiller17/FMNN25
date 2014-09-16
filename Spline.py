class Spline

    def __init__(self,uk,d) #Labinot
        """ 
        :param array.float uk: 
        :param array.float d:
        :raises TypeError: if uk is not an array of floats
        :raises TypeError: if d is not an array of floats
        :raises ValueError: if relation between uk and d is wrong
        """

        if not issubclass(uk.dtype.type,float):
            raise TypeError("uk must be an array of floats")
        if not issubclass(d.dtype.type,float):
            raise TypeError("d must be an array of floats")       
        if len(uk) != d.shape[1]:
            raise ValueError("relation between uk and d is wrong")
        if uk.size < 3:
            raise ValueError("uk must at least contain 3 elements")
        if not (uk == sorted(uk)).all():
            raise ValueError("uk must be a sorted array")
    
        self.d = array(d,dtype='float')
        self.uk = array(uk,dtype='float')
    
    def __call__(self,u): #Eli  
        """ 
        :param array.float u: array with all point where one wants to evaluate the spline
        :raises TypeError: if u is not an array of floats
        :raises ValueError: if any value in u is outside boundaries of uk
        :returns array.dim(2xn).floats Val: x and y values for the spline. x values in the first row and y values in the second.

        """

    def _eval(self,u): #Eli
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

        if not isinstance(u,float):
            raise TypeError("u must be a float")            
        if (u < self.uk[0]).any() or (u >= self.uk[-1]).any():
            raise ValueError("u is outside boundaries of uk")
        sortUk=sort(self.uk)            
        j = searchsorted(sortUk[1:],u)
            
        return int(j)


    def plot(self,uStart=uk[2],uStop=uk[-3]-kinds.default_float_kind.MIN*100,numPoints=1000): #FeelFree
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

    def getBaseFunc(self,j): #Hanna

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
                index = _findHotIntervall(u)
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
                    nVal = (u-uk[j])*(uk[j+3]-u)**2/((uk[j+3]-uk[j])*(uk[j+3]-uk[j+1])(uk[j+3-uk[j+2]]))
                    nVal += (uk[j+4] - u)*(u-uk[j+1])*(uk[j+3]-u)/((uk[j+4]-uk[j+1])*(uk[j+3]-uk[j+2])*(uk[j+3]-uk[j+2]))
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
                    nVal=(u-uk(j))**3/((uk[j+3]-uk[j])*(uk[j+2]-uk[j])*(uk[j+1]-uk[j]))
                except ZeroDivisionError:
                    nVal = 0
            return nVal

        return N
