import sys
sys.path = sys.path + ['../']
from nose.tools import raises,with_setup
import numpy as np
import Spline

class TestFindHotInterval:
    def setUp(self):
        self.testSpline = Spline.Spline(np.linspace(0,2),np.vstack((np.linspace(0,1),np.linspace(0,1))))

    def tearDown(self):
        del self.testSpline

    #@with_setup(setupFunc,tearDownFunc)
    @raises(ValueError)
    def testNotInInterval(self):
        self.testSpline._findHotInterval(3.)

    #@with_setup(setupFunc,tearDownFunc)
    @raises(TypeError)
    def testNotFloatIn(self):
        self.testSpline._findHotInterval(1)

    #@with_setup(setupFunc,tearDownFunc)
    def testFindLeftmostInterval(self):
        j = self.testSpline._findHotInterval(0.)
        assert j == 0

    #@with_setup(setupFunc,tearDownFunc)
    @raises(ValueError)
    def testRightEdgeNotPartOfInterval(self):
        j = self.testSpline._findHotInterval(2.)

    #@with_setup(setupFunc,tearDownFunc)
    def testReturnIsInt(self):
        j = self.testSpline._findHotInterval(1.)
        assert type(j) == int

    #@with_setup(setupFunc,tearDownFunc)
    def testReturnsAllIntervals(self):
        testIntervals = np.linspace(0,2.)
        testIntervals = testIntervals + 0.02
        j = 0
        for i in testIntervals[:-1]:
            Interval = self.testSpline._findHotInterval(i)
            assert j == Interval
            j += 1

class TestCallMethod:
    def setUp(self):
        self.testSpline = Spline.Spline(np.array([0,1,2,3,4,5,6,7],dtype='float64'),np.array([[-6,-4,-2,2,4,6],[-6,-4,-2,2,4,6]],dtype='float64'))
    def tearDown(self):
        del self.testSpline

    @raises(TypeError)
    def testNotArray(self):
        self.testSpline(0.0)

    @raises(TypeError)
    def testNotFloats(self):
        self.testSpline(np.array([1],dtype='int64'))

    @raises(ValueError)
    def testNotInInterval(self):
        self.testSpline(np.array([-1],dtype='float64'))

    def testReturnsArray(self):
        value = self.testSpline(np.array([3.5],dtype='float64'))
        assert type(value) == np.ndarray

    def testReturnsFloats(self):
        value = self.testSpline(np.array([3.5],dtype='float64'))
        assert value.dtype == 'int64'

    def testReturnsCorrectSize(self):
        values = self.testSpline(np.array([1.75,3.5,5.25],dtype='float64'))
        assert values.shape == (2,3)

    def testEvaluateOnePoint(self):
        value = self.testSpline(np.array([3.5],dtype='float64'))
        assert value[0,0] == 0
        assert value[1,0] == 0

    def testEvaluateManyPoints(self):
        values = self.testSpline(np.array([1.75,3.5,5.25],dtype='float64'))
        assert values[0,0] == -3
        assert values[1,0] == -3
        assert values[0,1] == 0
        assert values[1,1] == 0
        assert values[0,2] == 3
        assert values[1,2] == 3
        
class TestEval:

    def setUp(self):
        uk = np.array([0,1,2,3,4,5,6,7])
        d = np.array([[-6,-4,-2,2,4,6],[-6,-4,-2,2,4,6]])
        self.testSpline = Spline.Spline(uk,d)

    def tearDown(self):
        del self.testSpline
        
    @raises(ValueError)
    def testValueError(self):
        self.testSpline._eval(3.)
        
    @raises(TypeError)
    def testTypeError(self):
        self.testSpline._eval(1)
        
    def testOutputDimensions(self):
        out = self.testSpline._eval(1.)
        assert (2,1) == out.shape
        
    def testOutputType(self):
        out = self.testSpline._eval(1.)
        assert type(out) == np.ndarray        
        
        
    def testEvalOnePoint(self):
        out = self.testSpline._eval(3.5)
        assert out[0,0] == 0
        assert out[1,0] == 0
