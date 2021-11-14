from delaunay import Delaunay2d
import random
import numpy

if __name__ == '__main__':
    xyPoints = [numpy.array([random.random(), random.random()]) for i in range(10)]
    delaunay = Delaunay2d(xyPoints)
    delaunay.show()
