import numpy
import math
import copy
import matplotlib.pyplot as plt


class Delaunay2d:
    epsilone = 1.23456789e-14

    def __init__(self, points):
        self.points = points[:]
        self.triangles = []
        self.edge2Triangles = {}
        self.boundaryEdges = set()
        self.appliedBoundaryEdges = None
        self.holes = None

        cg = numpy.zeros((2,), numpy.float64)
        for pt in points:
            cg += pt
        cg /= len(points)

        def distanceSquare(pt):
            d = pt - cg
            return numpy.dot(d, d)

        self.points.sort(key=distanceSquare)

        # создание первого треугольника
        # если получили нулевую область, отбросить точки
        area = 0.0
        index = 0
        stop = False
        while not stop and index + 2 < len(points):
            area = self.getArea(index, index + 1, index + 2)
            if abs(area) < self.epsilone:
                del self.points[index]
            else:
                stop = True
        if index <= len(self.points) - 3:
            tri = [index, index + 1, index + 2]
            self.makeCounterClockwise(tri)
            self.triangles.append(tri)

            e01 = (tri[0], tri[1])
            self.boundaryEdges.add(e01)
            e12 = (tri[1], tri[2])
            self.boundaryEdges.add(e12)
            e20 = (tri[2], tri[0])
            self.boundaryEdges.add(e20)

            e01 = self.makeKey(e01[0], e01[1])
            self.edge2Triangles[e01] = [0, ]

            e12 = self.makeKey(e12[0], e12[1])
            self.edge2Triangles[e12] = [0, ]

            e20 = self.makeKey(e20[0], e20[1])
            self.edge2Triangles[e20] = [0, ]

        else:
            return

        for i in range(3, len(self.points)):
            self.addPoint(i)


    def show(self, width=500, height=500):
        xmin = min([p[0] for p in self.points])
        ymin = min([p[1] for p in self.points])
        xmax = max([p[0] for p in self.points])
        ymax = max([p[1] for p in self.points])
        padding = 5
        w = width - 2 * padding
        h = height - 2 * padding

        for e in self.edge2Triangles:
            i1, i2 = e
            xp1 = padding + int(w * (self.points[i1][0] - xmin) / (xmax - xmin))
            yp1 = padding + int(h * (ymax - self.points[i1][1]) / (ymax - ymin))
            xp2 = padding + int(w * (self.points[i2][0] - xmin) / (xmax - xmin))
            yp2 = padding + int(h * (ymax - self.points[i2][1]) / (ymax - ymin))
            plt.plot(
                [xp1, xp2],
                [yp1, yp2],
                color="black"
            )

        plt.savefig('demo.png', bbox_inches='tight')

    def getTriangles(self):
        return self.triangles

    def getEdges(self):
        return self.edge2Triangles.keys()

    def getArea(self, ip0, ip1, ip2):
        d1 = self.points[ip1] - self.points[ip0]
        d2 = self.points[ip2] - self.points[ip0]
        return (d1[0] * d2[1] - d1[1] * d2[0])

    def isEdgeVisible(self, ip, edge):
        area = self.getArea(ip, edge[0], edge[1])
        if area < self.epsilone:
            return True
        return False

    def makeCounterClockwise(self, ips): # Изменение порядка узлов
        area = self.getArea(ips[0], ips[1], ips[2])
        if area < -self.epsilone:
            ip1, ip2 = ips[1], ips[2]
            ips[1], ips[2] = ip2, ip1

    def flipOneEdge(self, edge):
        res = set()

        tris = self.edge2Triangles.get(edge, [])
        if len(tris) < 2:
            return res

        iTri1, iTri2 = tris
        tri1 = self.triangles[iTri1]
        tri2 = self.triangles[iTri2]

        iOpposite1 = -1
        iOpposite2 = -1
        for i in range(3):
            if not tri1[i] in edge:
                iOpposite1 = tri1[i]
            if not tri2[i] in edge:
                iOpposite2 = tri2[i]

        da1 = self.points[edge[0]] - self.points[iOpposite1]
        db1 = self.points[edge[1]] - self.points[iOpposite1]
        da2 = self.points[edge[0]] - self.points[iOpposite2]
        db2 = self.points[edge[1]] - self.points[iOpposite2]
        crossProd1 = self.getArea(iOpposite1, edge[0], edge[1])
        crossProd2 = self.getArea(iOpposite2, edge[1], edge[0])
        dotProd1 = numpy.dot(da1, db1)
        dotProd2 = numpy.dot(da2, db2)
        angle1 = abs(math.atan2(crossProd1, dotProd1))
        angle2 = abs(math.atan2(crossProd2, dotProd2))

        if angle1 + angle2 > math.pi * (1.0 + self.epsilone):
            newTri1 = [iOpposite1, edge[0], iOpposite2]  # triangle a
            newTri2 = [iOpposite1, iOpposite2, edge[1]]  # triangle b

            self.triangles[iTri1] = newTri1
            self.triangles[iTri2] = newTri2

            del self.edge2Triangles[edge]

            e = self.makeKey(iOpposite1, iOpposite2)
            self.edge2Triangles[e] = [iTri1, iTri2]

            e = self.makeKey(iOpposite1, edge[1])
            v = self.edge2Triangles[e]
            for i in range(len(v)):
                if v[i] == iTri1:
                    v[i] = iTri2
            res.add(e)

            e = self.makeKey(iOpposite2, edge[0])
            v = self.edge2Triangles[e]
            for i in range(len(v)):
                if v[i] == iTri2:
                    v[i] = iTri1
            res.add(e)

            res.add(self.makeKey(iOpposite1, edge[0]))
            res.add(self.makeKey(iOpposite2, edge[1]))

        return res

    def flipEdges(self):
        edgeSet = set(self.edge2Triangles.keys())

        continueFlipping = True

        while continueFlipping:
            newEdgeSet = set()
            for edge in edgeSet:
                newEdgeSet |= self.flipOneEdge(edge)

            edgeSet = copy.copy(newEdgeSet)
            continueFlipping = (len(edgeSet) > 0)

    def addPoint(self, ip):
        boundaryEdgesToRemove = set()
        boundaryEdgesToAdd = set()

        for edge in self.boundaryEdges:

            if self.isEdgeVisible(ip, edge):
                newTri = [edge[0], edge[1], ip]
                newTri.sort()
                self.makeCounterClockwise(newTri)
                self.triangles.append(newTri)

                e = list(edge[:])
                e.sort()
                iTri = len(self.triangles) - 1
                self.edge2Triangles[tuple(e)].append(iTri)

                e1 = [ip, edge[0]]
                e1.sort()
                e1 = tuple(e1)
                e2 = [edge[1], ip]
                e2.sort()
                e2 = tuple(e2)
                v1 = self.edge2Triangles.get(e1, [])
                v1.append(iTri)
                v2 = self.edge2Triangles.get(e2, [])
                v2.append(iTri)
                self.edge2Triangles[e1] = v1
                self.edge2Triangles[e2] = v2

                boundaryEdgesToRemove.add(edge)
                boundaryEdgesToAdd.add((edge[0], ip))
                boundaryEdgesToAdd.add((ip, edge[1]))

        for bedge in boundaryEdgesToRemove:
            self.boundaryEdges.remove(bedge)
        for bedge in boundaryEdgesToAdd:
            bEdgeSorted = list(bedge)
            bEdgeSorted.sort()
            bEdgeSorted = tuple(bEdgeSorted)
            if len(self.edge2Triangles[bEdgeSorted]) == 1:
                self.boundaryEdges.add(bedge)

        flipped = True
        while flipped:
            flipped = self.flipEdges()

    def makeKey(self, i1, i2):
        if i1 < i2:
            return (i1, i2)
        return (i2, i1)