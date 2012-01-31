from __future__ import division
import math
import random
import heapq
import time

import numpy as np
import visual as vis
import networkx as nx
import scipy.spatial as sps
#import matplotlib.pyplot as plt


def scaleA(a, n):
    """scales a given A configuration using n
    this should be used before the a is given to forwardMap"""

    r = 1.0 / (n - 1)
    aScaled = a * r
    return aScaled

def unScaleA(a, n):
    """scales a given A configuration using n
    this should be used before the a is given to forwardMap"""

    r = 1.0 / (n - 1)
    aUnScaled = a / r
    return aUnScaled

def convertFromCubeCoordinates(coordinates, n):
    newCoordinates = np.dot(np.diag(np.array([200, 200, 20])), 
                            coordinates - (.5 * np.ones_like(coordinates)))
    return scaleA(newCoordinates, n)

def convertToCubeCoordinates(a, n):
    aUnScaled = unScaleA(a, n)
    newCoordinates = np.dot(np.diag(np.array([1.0 / 200, 1.0 / 200, 1.0 / 20])), 
                            aUnScaled + (np.array([100, 100, 10])))
    return newCoordinates

def forwardMap(aScaled, n):
    """Computes the forward kinematic map for 
    a given a(configuration) and n for a discrete kinematic elastic chain."""

    #useful constants
    r = 1.0 / (n - 1)

    #check if input a is nonsingular
    rem1 = aScaled[1] % math.pi * (n - 1)
    rem2 = aScaled[2] % math.pi
    if (np.allclose(rem1, 0.0) or np.allclose((n - 1) * math.pi - rem1, 0.0)) \
            and (np.allclose(rem2, 0.0) or np.allclose(math.pi - rem2, 0.0)):
        print "Invalid Costate tried"
        return None, None, None

    p = np.zeros((n + 1, 3))
    x = np.zeros((n + 1, 3))
    u = np.zeros(n)
    JInvT = np.eye(3)
    p[0] = aScaled

    #first step has r = 0
    #recursion step
    u[0] = -p[0][2]
    x[1] = np.array([0, 0, u[0]]) + x[0]
    p[1] = p[0]

    for i in range(n)[1:]:
        #constants used in the calculations
        rcosxi3 = r * math.cos(x[i][2])
        rsinxi3 = r * math.sin(x[i][2])
        JInvT[2][0] = rsinxi3
        JInvT[2][1] = -rcosxi3

        #recursion step
        u[i] = -np.dot(np.array([-rsinxi3, rcosxi3, 1]), p[i])
        x[i + 1] = np.array([rcosxi3, rsinxi3, u[i]]) + x[i]
        p[i + 1] = np.dot(JInvT, p[i])
    return u, x, p

def extractXCurve(point):
    """after forwardMap has been called, call this 
    on the x parameter to get a visual curve of the elastica"""
    vcurve = vis.curve()
    for part in point.x:
        vcurve.append([part[0], part[1]])
    return vcurve

def planner(qInit, qGoal, rho, iterations):
    """Single Query Bi-Directional Lazy (SBL) planner by Sanchez, Latombe.
    Plans a path between qInit and qGoal."""
    try:
        tInit = Tree((10, 10), qInit)
    except InputError:
        print 'Initial point is not stable.'
        raise
    try:
        tGoal = Tree((10, 10), qGoal)
    except InputError:
        print 'Goal point is not stable.'
        raise
    curves = []
    for i in range(iterations)[1:]:
        if (i % 50) == 0:
            tInit.reconfigureArrays()
            tGoal.reconfigureArrays()
        choice = random.random()
        if choice < 0.5:
            milestoneTree, notMilestoneTree = tInit, tGoal
        else:
            milestoneTree, notMilestoneTree = tGoal, tInit
        milestone = expandTree(milestoneTree, rho)
        path = connectTrees(milestone, milestoneTree, notMilestoneTree, rho)
        if (i % 5) == 0:
            for curve in curves:
                curve.visible = False
            curves = []
            curves.extend(displayGraph(tInit.graph.edges(), vis.color.green))
            curves.extend(displayGraph(tGoal.graph.edges(), vis.color.blue))
            #time.sleep(.2) #for animating
        if path is not None:

            #print path[0]
            for curve in curves:
                curve.visible = False
            displayGraph(tInit.graph.edges(), vis.color.green)
            displayGraph(tGoal.graph.edges(), vis.color.blue)
            #tInit.graph.add_edges_from(tGoal.graph.reverse().edges())
            if choice < .5:
                completePath = conglomeratePath(path[0], tInit, path[1], tGoal)
            else:
                completePath = conglomeratePath(path[1], tInit, path[0], tGoal)
            #for start, end in tGoal.graph.edges():
                #intermediateConfigurations = tGoal.graph[start][end]['intermediateConfigurations']
                #intermediateConfigurations.reverse()
                #tInit.addEdge(end, start, weight=tGoal.graph[start][end]['weight'],
                                #intermediateConfigurations=intermediateConfigurations)
            #completePath = nx.shortest_path(tInit.graph, qInit, qGoal)
            #if choice < .5:
                #completePath = joinPaths(tInit, path[1], tGoal, path[0])
            #else:
                #completePath = joinPaths(tInit, path[0], tGoal, path[1])
            #pathCurve = vis.curve(color=vis.color.orange)
            #for node in completePath:
        if choice < 0.5:
            milestoneTree, notMilestoneTree = tInit, tGoal
        else:
            milestoneTree, notMilestoneTree = tGoal, tInit
                #pathCurve.append(node.coordinates)
            return completePath
    return False

def plannerB(pointInit, bEnd, rho, iterations):
    """RRT planner to get a startpoint to the bEnd """
    try:
        tInit = Tree((10, 10), pointInit)
    except InputError:
        print 'Initial point is not stable.'
        raise
    curves = []
    for i in range(iterations)[1:]:
        if (i % 50) == 0:
            tInit.reconfigureArrays()
        milestone = expandTreeB(tInit, rho)
        path = connectTreesB(milestone, tInit, bEnd, rho)
        if (i % 5) == 0:
            for curve in curves:
                curve.visible = False
            curves = []
            curves.extend(displayGraph(tInit.graph.edges(), vis.color.green))
            #time.sleep(.2) #for animating
        if path is not None:

            #print path[0]
            for curve in curves:
                curve.visible = False
            displayGraphB(tInit.graph.edges(), vis.color.green)
            #tInit.graph.add_edges_from(tGoal.graph.reverse().edges())
            completePath = conglomeratePathB(path, tInit)
            #for start, end in tGoal.graph.edges():
                #intermediateConfigurations = tGoal.graph[start][end]['intermediateConfigurations']
                #intermediateConfigurations.reverse()
                #tInit.addEdge(end, start, weight=tGoal.graph[start][end]['weight'],
                                #intermediateConfigurations=intermediateConfigurations)
            #completePath = nx.shortest_path(tInit.graph, qInit, qGoal)
            #if choice < .5:
                #completePath = joinPaths(tInit, path[1], tGoal, path[0])
            #else:
                #completePath = joinPaths(tInit, path[0], tGoal, path[1])
            #pathCurve = vis.curve(color=vis.color.orange)
            #for node in completePath:
                #pathCurve.append(node.coordinates)
            return completePath
    return False

def joinPaths(firstTree, firstHalfPath, secondTree, secondHalfPath):
    firstPath = firstTree.graph.subgraph(firstHalfPath)
    #print firstPath
    secondPath = secondTree.graph.subgraph(secondHalfPath).reverse()
    #print secondPath
    firstPath.add_edges_from(secondPath.edges())
    return firstPath

def connectTrees(milestone, milestoneTree, notMilestoneTree, rho):
    """Attempts to connect the most recently added-to tree to 
    the other tree.  Tries two configurations, tests for distance, and 
    then connects."""
    newMilestone = notMilestoneTree.findCloseNewMilestone(milestone)
    if newMilestone is None:
        newMilestone = notMilestoneTree.getUniformSample()
    pathSuccess = None
    #try twice, once nearby, once random
    for i in range(2):
        distance = l2Distance(milestone, newMilestone)
        if distance < rho:
            milestoneTree.addEdge(milestone, newMilestone, weight=distance)

            #get the path from the root to the other root to test
            node = newMilestone
            firstHalfPath = nx.shortest_path(milestoneTree.graph, source=milestoneTree.root, target=node)
            #print firstHalfPath
            secondHalfPath = nx.shortest_path(notMilestoneTree.graph, source=notMilestoneTree.root, target=node)
            pathSuccess = testPath(milestoneTree, firstHalfPath, notMilestoneTree, secondHalfPath)

            #check if the path is clear.  if it is, it will be not None,
            #otherwise it will be None
            if pathSuccess is not None:
                return firstHalfPath, secondHalfPath
        if pathSuccess is None or i is 1:
            return None
        newMilestone = notMilestoneTree.getUniformSample()

def connectTreesB(milestone, beginTree, bEnd, rho, epsilon=.001):
    """Attempts to connect the most recently added-to tree to 
    the other tree.  Tries two configurations, tests for distance, and 
    then connects."""
    distance = l2DistanceBp2Undefined(milestone, bEnd)
    if distance < rho:
        pathResults = edgeBetweenB(milestone, bEnd, epsilon)
        if pathResults:
            beginTree.addEdgeB(milestone, pathResults[1], intermediateConfigurations=pathResults[0])
            #get the path from the root to the other root to test
            node = pathResults[1]
            path = nx.shortest_path(beginTree.graph, source=beginTree.root, target=node)
            return path
        else:
            #no path exist between this node and the endpoint
            return None

def expandTree(tree, rho, iterations=100):
    possibleSample = None
    while True:
        milestone = tree.getRandomSample()
        for i in range(iterations)[1:]:
            possibleSample = giveSampleFromCube(milestone, rho / i)
            if possibleSample.isStable():
                tree.addEdge(milestone, possibleSample)
                #print len(tree.graph.nodes())
                return possibleSample

def expandTreeB(tree, rho, iterations=100, epsilon=.001):
    while True:
        milestone = tree.getRandomSample()
        for i in range(iterations)[1:]:
            possibleB = giveSampleFromCubeB(milestone, rho / i)
            pathResults = edgeBetweenB(milestone, possibleB, epsilon)
            if pathResults:
                tree.addEdgeB(milestone, pathResults[1], intermediateConfigurations=pathResults[0])
                #print len(tree.graph.nodes())
                return pathResults[1]

def giveSampleFromCube(point, radius):
    """Samples an (approximately) uniformly random sample from within a cube 
    around point, while staying within the bounds of [0, 1]^n, 
    where n is the dimension of point's coordinate"""
    coordinates = []
    for i in range(len(point.coordinates)):
        newCoordinate = 2.0 * radius * (random.random() - 0.5)
        newCoordinate += point.coordinates[i]
        if newCoordinate > 1.0:
            newCoordinate = 1.0
        if newCoordinate < 0.0:
            newCoordinate = 0.0
        coordinates.append(newCoordinate)
    return Point(coordinates)

def giveSampleFromCubeB(point, radius):
    """Samples an (approximately) uniformly random sample from within a cube 
    around point, while staying within the bounds of [-1, 1]^(n - 1), 
    where n is the dimension of point's coordinate.  The last 
    coordinate is untouched."""
    coordinates = []
    while not coordinates:
        for i in range(len(point.coordinates) - 1):
            newCoordinate = 2.0 * radius * (random.random() - 0.5)
            newCoordinate += point.b[i]
            if newCoordinate > 1.0:
                newCoordinate = 1.0
            if newCoordinate < -1.0:
                newCoordinate = -1.0
            coordinates.append(newCoordinate)
        newCoordinate = 2.0 * radius * (random.random() - 0.5)
        newCoordinate += point.coordinates[i]
        coordinates.append(newCoordinate)
        if np.linalg.norm(coordinates[:-1]) > 1.0:
            coordinates = []
    return coordinates

def testPath(milestoneTree, firstHalfPath, notMilestoneTree, secondHalfPath, epsilon=0.003):
    """Checks a potential path for collision.  Maintains a priority queue 
    of the segments to be checked so that the ones most likely to collide are 
    checked first.  Collision free segments take the longest time to check."""
    #initialize heap for edge testing order.
    #weight is negative edge weight since the python heapq implementation is a minheap 
    #but we want a maxheap
    segmentsToCheck = []
    #print firstHalfPath
    #try:
        #firstHalfPath.reverse()
    #except AttributeError:
        #plt.figure()
        #nx.draw(milestoneTree.graph)
        #plt.figure()
        #nx.draw(notMilestoneTree.graph)
        #plt.show()

    #print 'firstHalfPath: ', firstHalfPath
    #print 'secondHalfPath: ', secondHalfPath
    print 'testing Path....'
    for prevNode, node in zip(firstHalfPath[:-1], firstHalfPath[1:]):
        #print milestoneTree.graph.edges()
        #print milestoneTree.graph[node]
        #print prevNode, node
        heapq.heappush(segmentsToCheck, (-milestoneTree.graph[prevNode][node]['weight'], prevNode, node, milestoneTree, notMilestoneTree))
    #keep these edges reversed so that the edges that they reference appear in the 
    #appropriate tree.  at the end we will flip these edges accordingly
    #print secondHalfPath
    #try:
        #secondHalfPath.reverse()
    #except AttributeError:
        #plt.figure()
        #nx.draw(milestoneTree.graph)
        #plt.figure()
        #nx.draw(notMilestoneTree.graph)
        #plt.show()

    for prevNode, node in zip(secondHalfPath[:-1], secondHalfPath[1:]):
        heapq.heappush(segmentsToCheck, (-notMilestoneTree.graph[prevNode][node]['weight'], prevNode, node, notMilestoneTree, milestoneTree))

    i = 0
    while segmentsToCheck:
        #print i
        i += 1
        currentEdge = heapq.heappop(segmentsToCheck)
        #print currentEdge
        testResults = testSegment(currentEdge)
        #this means that collision has been detected
        if testResults < 0.0:
            currentEdge[3].removeEdge(currentEdge[1], currentEdge[2], currentEdge[4], firstHalfPath[-1])
            print 'collision resolved. returning to sampling'
            return None
        #means it hasn't been shown to be safe yet
        if testResults > epsilon:
            heapq.heappush(segmentsToCheck, (-testResults, currentEdge[1], currentEdge[2], currentEdge[3], currentEdge[4]))
    return True

def testSegment(edgeData):
    """Bisects each segment of the edge. Essentially doubles the number of 
    checked collisions for this edge.  Returns negative if a collision has been 
    detected; otherwise positive."""
    #print edgeData
    #print edgeData[3].graph.edges()
    edgeDict = edgeData[3].graph[edgeData[1]][edgeData[2]]

    #start bisecting between start node
    prevNode = edgeData[1]
    configurationsChecked = []
    #print '# of intermediateConfigurations: ', len(edgeDict['intermediateConfigurations'])
    for node in edgeDict['intermediateConfigurations']:
        #print 'intermediateConfiguration: ', node
        newNode = bisectionConfiguration(prevNode, node)
        if newNode.isStable():
            configurationsChecked.append(newNode)
        else:
            return -1.0
        configurationsChecked.append(node)
        prevNode = node

    #end by bisecting between end node
    newNode = bisectionConfiguration(prevNode, edgeData[2])
    #print 'new node: ', newNode
    if newNode.isStable():
        configurationsChecked.append(newNode)
    else:
        print 'is not stable: ', newNode
        print 'of this edge: ', edgeData[1], edgeData[2]
        return -1.0
    edgeDict['weight'] = edgeDict['weight'] / 2.0
    edgeDict['intermediateConfigurations'] = configurationsChecked
    edgeData[3].graph[edgeData[1]][edgeData[2]]['intermediateConfigurations'] = configurationsChecked
    edgeData[3].graph[edgeData[1]][edgeData[2]]['weight'] = edgeDict['weight']
    #print 'weight: ', edgeDict['weight']
    return edgeDict['weight']

def bisectionConfiguration(node1, node2):
    """Returns a Point in which the coordinates bisect 
    the coordinates of node1 and node2"""
    return Point((node1.coordinates + node2.coordinates) / 2.0)

class Tree(object):
    def __init__(self, arrayShape, root):
        self.graph = nx.DiGraph()
        self.densityArray = np.zeros(arrayShape)

        #3d array h to be filled with points in the corresponding bins
        #try:
            #self.pointArray = [[[[] for k in range(arrayShape[2])] \
                #for j in range(arrayShape[1])] for i in range(arrayShape[0])]
        #except IndexError:
        self.pointArray = [[[] for j in range(arrayShape[1])] \
                for i in range(arrayShape[0])]

        self.graph.add_node(root)
        self.numberInArray = 0
        self.arrayDimensions = [0, 1]
        self.maxArrayDimension = len(root.coordinates)
        self.sortedListOfCells = []
        self.setInArray = set()
        self.addNodeToArray(root)
        self.root = root
        stable = self.root.isStable()
        if not stable:
            print 'root of tree: ', root, ' is not stable.'
            raise InputError('Tree root is not stable')

    def addEdge(self, startPoint, endPoint, weight=None, intermediateConfigurations=None):
        """Add an edge to the tree.  Adds in edge attributes weight 
        and stores the collision checking 
        data and beginning/intermediate/end configurations."""
        if weight is None:
            weight = l2Distance(startPoint, endPoint)
        self.graph.add_edge(startPoint, endPoint)
        #print "edges: ", self.graph.edges() 
        #print "graph type: ", type(self.graph)
        #print "edges from startPoint: ", self.graph[startPoint]
        #print "edge type: ", type(self.graph[startPoint])
        #print "edge dict: ", self.graph[startPoint][endPoint]
        self.graph[startPoint][endPoint]['weight'] = weight
        if intermediateConfigurations is None:
            self.graph[startPoint][endPoint]['intermediateConfigurations'] = []
        else:
            self.graph[startPoint][endPoint]['intermediateConfigurations'] = intermediateConfigurations
        self.addNodeToArray(startPoint)
        self.addNodeToArray(endPoint)

    def addEdgeB(self, startPoint, endPoint, weight=None, intermediateConfigurations=None):
        """Add an edge to the tree.  Adds in edge attributes weight 
        and stores the collision checking 
        data and beginning/intermediate/end configurations."""
        if weight is None:
            weight = l2DistanceB(startPoint, endPoint)
        #print "edges: ", self.graph.edges() 
        #print "graph type: ", type(self.graph)
        #print "edges from startPoint: ", self.graph[startPoint]
        #print "edge type: ", type(self.graph[startPoint])
        #print "edge dict: ", self.graph[startPoint][endPoint]
        self.graph.add_edge(startPoint, endPoint)
        self.graph[startPoint][endPoint]['weight'] = weight
        if intermediateConfigurations is None:
            self.graph[startPoint][endPoint]['intermediateConfigurations'] = []
        else:
            self.graph[startPoint][endPoint]['intermediateConfigurations'] = intermediateConfigurations
        self.addNodeToArray(startPoint)
        self.addNodeToArray(endPoint)


    def removeEdge(self, startPoint, endPoint, otherTree, linkingNode):
        """Removes the edge in this tree corresponding to startPoint 
        and endPoint, and then transfers any nodes below the edge on 
        this tree to the other tree, at insertionPoint."""
        print 'collision detected. transferring edges......'
        graphBelow = nx.traversal.dfs_tree(self.graph, endPoint)
        pathToConnection = nx.shortest_path(graphBelow, source=endPoint, target=linkingNode)
        #graphToReverse = nx.traversal.dfs_tree(self.graph, pathToConnection[1])
        for start, end in zip(pathToConnection[:-1], pathToConnection[1:]):
            intermediateConfigurations = self.graph[start][end]['intermediateConfigurations']
            intermediateConfigurations.reverse()
            #print 'edge transfer: ', start, end, intermediateConfigurations
            otherTree.addEdge(end, start, weight=self.graph[start][end]['weight'],
                            intermediateConfigurations=intermediateConfigurations)
            #otherTree.addNodeToArray(start)
            #otherTree.addNodeToArray(end)
            #self.removeNodeFromArray(end)
            self.graph.remove_edge(start, end)
            graphBelow.remove_edge(start, end)
        #curves = displayGraph(graphBelow.edges(), vis.color.red)
        #time.sleep(.1)
        #for curve in curves:
            #curve.visible = False
        for start, end in graphBelow.edges():
            #if (start, end) in edgesToReverse:
                #print 'edge found in edgeToReverse: ', start, end
                #intermediateConfigurations = self.graph[start][end]['intermediateConfigurations']
                #intermediateConfigurations.reverse()
                #print 'edge transfer: ', start, end, intermediateConfigurations
                #otherTree.addEdge(end, start, weight=self.graph[start][end]['weight'],
                                #intermediateConfigurations=intermediateConfigurations)
                #otherTree.addNodeToArray(start)
                #time.sleep(3)
            #else:
            intermediateConfigurations = self.graph[start][end]['intermediateConfigurations']
            #print 'edge transfer: ', start, end, intermediateConfigurations
            #print 'not reversed'
            otherTree.addEdge(start, end, weight=self.graph[start][end]['weight'],
                            intermediateConfigurations=intermediateConfigurations)
            #otherTree.addNodeToArray(start)
            self.graph.remove_edge(start, end)
            graphBelow.remove_edge(start, end)

        for node in graphBelow.nodes():
            otherTree.addNodeToArray(node)
            self.removeNodeFromArray(node)
            self.graph.remove_node(node)
        #self.graph.remove_edges_from(graphBelow.edges())
        #self.reconfigureArrays() #trying to fix the issues
        #print 'edges Below:', graphBelow
        #print 'edges edges:', graphBelow.edges()
        #print 'edges remaining in tree: ', self.graph.edges()
        #print 'nodes remaining in tree: ', self.graph.nodes()

    def addNodeToArray(self, node):
        """Add a node to both densityArray and pointArray."""
        if node in self.setInArray:
            return None
        else:
            indices = self.getArrayIndices(node)
            self.densityArray[indices] += 1
            if self.densityArray[indices] == 1:
                self.sortedListOfCells.append(indices)
            self.sortedListOfCells.sort(key=lambda indicies: self.densityArray[indicies])
            self.numberInArray += 1
            self.pointArray[indices[0]][indices[1]].append(node)
            self.setInArray.add(node)

    def removeNodeFromArray(self, node):
        """Remove a node from both densityArray and pointArray.
        Raises a KeyError if the node is not in pointArray."""
        if node not in self.setInArray:
            return None
        else:
            indices = self.getArrayIndices(node)
            self.densityArray[indices] -= 1
            if self.densityArray[indices] == 0:
                self.sortedListOfCells.remove(indices)
            self.sortedListOfCells.sort(key=lambda indicies: self.densityArray[indicies])
            self.numberInArray -= 1
            self.pointArray[indices[0]][indices[1]].remove(node)
            self.setInArray.remove(node)

    def getArrayIndices(self, node):
        """Gets the array indices for the current array configuration 
        in the tree, corresponding to the input node. Returned 
        as a tuple."""
        firstDimensionTick = 1.0 / self.densityArray.shape[0]
        secondDimensionTick = 1.0 / self.densityArray.shape[1]
        #print self.arrayDimensions
        #print self.densityArray.shape
        #print node
        #print node.coordinates
        firstIndex = int(self.densityArray.shape[0] * (node.coordinates[self.arrayDimensions[0]] % firstDimensionTick))
        secondIndex = int(self.densityArray.shape[1] * (node.coordinates[self.arrayDimensions[1]] % secondDimensionTick))
        return (firstIndex, secondIndex)

    def reconfigureArrays(self):
        """Rotates the dimensions that the arrays project onto.
        This ensures that the points diffuse in all dimensions. 
        Then clears the current arrays, and inputs all of the nodes 
        with respect to the new array structure."""
        #print self.maxArrayDimension
        self.arrayDimensions[0] = (self.arrayDimensions[0] + len(self.arrayDimensions)) % self.maxArrayDimension
        self.arrayDimensions[1] = (self.arrayDimensions[1] + len(self.arrayDimensions)) % self.maxArrayDimension
        self.densityArray = np.zeros_like(self.densityArray)
        self.pointArray = [[[] for j in range(self.densityArray.shape[1])] \
                for i in range(self.densityArray.shape[0])]
        self.setInArray = set()
        self.sortedListOfCells = []
        self.numberInArray = 0
        map(self.addNodeToArray, self.graph.nodes())

    def getRandomSample(self):
        """Picks a node according to the probability density associated 
        with the densityArray."""
        #currently using a linear distribution on the sorted list of cells
        choice = random.random()
        n = len(self.sortedListOfCells)
        #print n
        cellChoice = int(math.floor(n * (1 - math.sqrt(1 - choice))))
        #print cellChoice
        indicies = self.sortedListOfCells[cellChoice]
        return random.choice(self.pointArray[indicies[0]][indicies[1]])

    def getUniformSample(self):
        """Picks a node according to a uniform probability 
        distribution on the nodes."""
        #print self.graph.nodes()
        return random.choice(self.graph.nodes())

    def findCloseNewMilestone(self, milestone):
        """Gives the closest node to the milestone that is in the same array cell"""
        indices = self.getArrayIndices(milestone)
        if self.pointArray[indices[0]][indices[1]]:
            closeMilestone = min(self.pointArray[indices[0]][indices[1]], key=lambda point : lInfDistance(point, milestone))
        else:
            closeMilestone = None
        #print closeMilestone
        return closeMilestone

def expandTowardsB(p, q, epsilon):
    Z = p.Zn()
    deltab = q - p.b
    deltaa = np.linalg.solve(Z, deltab)
    deltaamag = np.linalg.norm(deltaa)
    if deltaamag < epsilon:
        newPoint = Point(None, a=deltaa + p.a)
        needToContinueExpanding = False
    else:
        newPoint = Point(None, a=deltaa / deltaamag * epsilon + p.a)
        needToContinueExpanding = True
    if newPoint.isStable():
        return newPoint, needToContinueExpanding, True
    else:
        return newPoint, needToContinueExpanding, False

def edgeBetweenB(p, q, epsilon):
    intermediateConfigurations = []
    currentp = p
    needToContinueExpanding = False
    while not needToContinueExpanding:
        currentp, needToContinueExpanding, stable = expandTowardsB(currentp, q, epsilon)
        if not stable:
            return False
        else:
            intermediateConfigurations.append(currentp)
            #add in last point as a node
            if not needToContinueExpanding:
                currentp, needToContinueExpanding, stable = expandTowardsB(currentp, q, epsilon)
                if not stable:
                    return False
                else:
                    if needToContinueExpanding:
                        intermediateConfigurations.append(currentp)
    return intermediateConfigurations, currentp

class Point(object):
    """Provides a class that is the basic object(node of manipulation) 
    in the planning class."""
    #TODO: add in conversion from script A coordinates to [0,1]^n coordinates and 
    #back.  Just store each copy separately.
    
    def __init__(self, coordinates, a=None):
        if a is not None:
            self.a = a
            self.coordinates = convertToCubeCoordinates(a, len(a))
        else:
            if coordinates is None:
                coordinates = np.array([])
            self.coordinates = np.array(coordinates)
            self.a = None
        self.Z = None
        self.stable = None

    #override container functions to provide functionality
    def __str__(self):
        return str(self.coordinates)

    def __repr__(self):
        return str(self.coordinates)

    def __getitem__(self, key):
        return Point(self.coordinates[key])

    def __len__(self):
        return len(self.coordinates)

    def Zn(self):
        if self.Z is None:
            Zi = np.zeros((self.n + 1, 3, 3))
            duda1 = -self.x[0][1]
            duda2 = self.x[0][0]
            duda3 = -1.0
            Zi[1, 2, :] = np.array([duda1, duda2, duda3])
            for i in range(self.n + 1)[1:]:
                rcosxi3 = self.r * math.cos(self.x[i][2])
                rsinxi3 = self.r * math.sin(self.x[i][2])
                duda1 = -self.x[i][1] - rsinxi3 - self.a[0] * (Zi[i][1][0] + \
                                rcosxi3 * Zi[i][2][0]) + self.a[1] * (Zi[i][0][0] - 
                                rsinxi3 * Zi[i][2][0])
                duda2 = self.x[i][0] + rcosxi3 - self.a[0] * (Zi[i][1][1] + \
                                rcosxi3 * Zi[i][2][1]) + self.a[1] * (Zi[i][0][1] - 
                                rsinxi3 * Zi[i][2][1])
                duda3 = - self.a[0] * (Zi[i][1][2] + rcosxi3 * Zi[i][2][2]) \
                    + self.a[1] * (Zi[i][0][2] - rsinxi3 * Zi[i][2][2]) - 1.0
                Zi[i, :, :] = Zi[i - 1, :, :] + np.vstack((rsinxi3 * Zi[i - 1, 2, :], 
                            rcosxi3 * Zi[i - 1, 2, :], np.array([duda1, duda2, duda3])))
            self.Z = Zi[self.n, :, :]
            return Zi[self.n, :, :]
        else:
            return self.Z

    def isStable(self):
        """Collision checking for a configuration.
        Discrete Riccatti recursion through the state to 
        make sure that the configuration is locally minimizing.
        Check to see that the given trajectory is stable.
        Returns True if it is stable, False otherwise.
        """
        if self.stable is not None:
            return self.stable

        #useful constants
        n = 40
        self.n = n
        r = 1.0 / (n - 1)
        self.r = r
        #calculate trajectories first
        if self.a is None:
            self.a = convertFromCubeCoordinates(self.coordinates, n)
        self.u, self.x, self.p = forwardMap(self.a, n)
        self.b = self.x[-1]

        
        #set up A matrix quickly
        A = np.zeros((12, 13))
        for i in range(4): #e3's
            A[3 * i + 2][i] = 1.0
        for i in range(9): #-I's and I's within J's
            A[i][4 + i] = -1.0
            A[i + 3][4 + i] = 1.0
        for i in range(3): #last column of J's
            A[3 + 3 * i][6 + 3 * i] = -r * (n - 3 + i) * math.sin(self.x[n - 3 + i][2])
            A[4 + 3 * i][6 + 3 * i] = r * (n - 3 + i) * math.cos(self.x[n - 3 + i][2])
        #debugging purposes
        #print "A is: "
        #print A

        #set up B matrix quickly
        B = np.zeros((12, 3))
        for i in range(3):
            B[i][i] = -1.0
        B[2][0] = r * (n - 4) * math.sin(self.x[n - 4][2])
        B[2][1] = -r * (n - 4) * math.cos(self.x[n - 4][2])
        #debugging purposes
        #print "B is: "
        #print B

        #set up M matrix
        M = np.zeros((13, 13))
        for i in range(4):
            M[i][i] = 1.0
        for i in range(3):
            M[5 + 3 * i][5 + 3 * i] = self.Q22(n - 3 + i, r)

        #debugging purposes
        #print "M is: "
        #print M

        u, s, vh = np.linalg.svd(A)
        #print "u is: "
        #print u
        #print "s is: "
        #print s
        #print "vh is: "
        #print vh
        #print s.shape
        #print np.diag(s).shape
        #print np.rank(A)


        #N is the null space of A
        count = 0
        for i in range(13):
            if np.allclose(np.zeros(12), np.dot(A, vh[i])):
                count += 1
        if count!=1:
            print "Nullspace of A is larger than expected!!!"
        N = vh[12]

        #test for Positive definiteness of M on the null space of A
        L = np.dot(N, np.dot(M, N))
        positiveDefinite = True
        if not L.shape:
            if L <= 0:
                positiveDefinite = False
        else:
            print "ERROR! L is not a scalar!"
            return False
            for eig in np.linalg.eigvals(L):
                if eig <= 0:
                    positiveDefinite = False

        if not positiveDefinite:
            return False
        else:
            APseudo = np.linalg.pinv(A)
            if not L.shape: #don't handle other cases currently
                K = 1.0 / L * np.dot(N, M)
                A_p_B = np.dot(APseudo, B)
                I_NK = np.add(np.eye(13), -np.outer(N, K))
                P_i_1 = np.dot(A_p_B.T, np.dot(I_NK.T, np.dot(M, np.dot(I_NK, A_p_B))))

                iterates = range(n - 4)
                iterates.reverse()
                for i in iterates:
                    s_i_1 = 1 + P_i_1[2][2]
                    if s_i_1 <= 0:
                        return False
                    else:
                        P_i = np.add(self.Corner22(self.Q22(i, r)), np.dot(self.J(i, r).T, np.dot(np.add(P_i_1, \
                                -np.dot(P_i_1, np.dot(self.Corner22(1.0 / s_i_1), P_i_1))), self.J(i, r))))
                        P_i_1 = P_i.copy()
                        

                self.stable = True
                edgeCurves = pairUpXCurve(self.x)
                tree = BoundingTree(edgeCurves[1:], n, 1, n-1)
                collision = checkForCollision(tree)
                return not collision

    def Corner22(self,input):
        """ Q matrix """
        ret = np.zeros((3,3))
        ret[2][2] = input
        return ret


    def Q22(self, i, r):
        """ Q matrix """
        return -r*(self.p[0][0]*math.cos(self.x[i][2]) + self.p[0][1]*math.sin(self.x[i][2]))

    def J(self, i, r):
        """ J matrix """
        ret_J = np.eye(3)
        ret_J[0][2] = -r * math.sin(self.x[i][2])
        ret_J[1][2] = r * math.cos(self.x[i][2])
        return ret_J

def lInfDistance(point1, point2, memo={}):
    """memoized l2Distance function"""
    try:
        return memo[point1, point2]
    except KeyError:
        try:
            return memo[point2, point1]
        except KeyError:
            distance = sps.distance.chebyshev(point1.coordinates, point2.coordinates)
            memo[point1, point2] = distance
            memo[point2, point1] = distance
            return distance

def l2Distance(point1, point2, memo={}):
    """memoized l2Distance function"""
    try:
        return memo[point1, point2]
    except KeyError:
        try:
            return memo[point2, point1]
        except KeyError:
            distance = sps.distance.euclidean(point1.coordinates, point2.coordinates)
            memo[point1, point2] = distance
            memo[point2, point1] = distance
            return distance

def l2DistanceB(point1, point2, memo={}):
    """memoized l2Distance function"""
    try:
        return memo[point1, point2]
    except KeyError:
        try:
            return memo[point2, point1]
        except KeyError:
            distance = sps.distance.euclidean(point1.b, point2.b)
            memo[point1, point2] = distance
            memo[point2, point1] = distance
            return distance

def l2DistanceBp2Undefined(point1, point2):
    """memoized l2Distance function"""
    distance = sps.distance.euclidean(point1.b, point2)
    return distance

def pairUpXCurve(point):
    edgesInCurve = []
    for firstPoint, secondPoint in zip(point[:-1], point[1:]):
        edgesInCurve.append((firstPoint[:-1], secondPoint[:-1]))
    return edgesInCurve

class BoundingTree(object):
    def __init__(self, edgesInCurve, n, indexStart, indexEnd):
        if len(edgesInCurve) is 1: 
            self.center = (edgesInCurve[0][0] + edgesInCurve[0][1]) / 2.0
            self.radius = 1.0 / (2 * (n - 1))
            self.curveSegments = edgesInCurve
            self.left = None
            self.right = None
            #self.sphere = vis.sphere(pos=self.center, color=vis.color.red, radius=self.radius)
            self.sphere = None
            self.isSphere = False
            #self.sphere = vis.sphere(pos=self.center, color=vis.color.green, radius=self.radius) #, opacity=(1 - (1.0 / math.log((level + 1))))
            #self.sphere.visible = False
            self.index = indexStart
        else:
            i = int(math.ceil(len(edgesInCurve) / 2.0))
            self.left = BoundingTree(edgesInCurve[:i], n, indexStart, indexStart + i - 1)
            self.right = BoundingTree(edgesInCurve[i:], n, indexStart + i, indexEnd)
            self.center = (self.left.center + self.right.center) / 2.0
            self.radius = max(self.left.radius, self.right.radius) + np.linalg.norm(self.center - self.left.center)
            self.curveSegments = edgesInCurve
            #self.sphere = vis.sphere(pos=self.center, color=vis.color.green, radius=self.radius) #, opacity=(1 - (1.0 / math.log((level + 1))))
            #self.sphere.visible = False
            self.isSphere = True
            self.indexStart, self.indexEnd = indexStart, indexEnd

    def unDisplaySpheres(self):
        def oneLevelUndisplay(root):
            if root is not None:
                try:
                    root.sphere.visible = False
                except AttributeError:
                    pass
                oneLevelUndisplay(root.left)
                oneLevelUndisplay(root.right)
        oneLevelUndisplay(self)

def checkForCollision(root):
    listOfCollisionsToCheck = [root]
    tupleType = type((1, 2))
    while listOfCollisionsToCheck:
        currentNode = listOfCollisionsToCheck.pop()
        if type(currentNode) is tupleType:
            #print 'type checked!'
            node1, node2 = currentNode
            fun = checkCollisionSingleQuery(node1, node2)
            #print fun
            check1, check2 = fun
            #both subdivide
            if check1 == 'subdivide' and check2 == 'subdivide':
                listOfCollisionsToCheck.append((node1.left, node2.left))
                listOfCollisionsToCheck.append((node1.right, node2.left))
                listOfCollisionsToCheck.append((node1.left, node2.right))
                listOfCollisionsToCheck.append((node1.right, node2.right))
            #one is at bottom, subdivide other
            elif check1 == 'subdivide' and check2 == 'ignore':
                listOfCollisionsToCheck.append((node1.left, node2))
                listOfCollisionsToCheck.append((node1.right, node2))
            #one is at bottom, subdivide other
            elif check1 == 'ignore' and check2 == 'subdivide':
                listOfCollisionsToCheck.append((node1, node2.left))
                listOfCollisionsToCheck.append((node1, node2.right))
            #collision detected between two lines
            elif check1 == 'collision' and check2 == 'collision':
                return True
        else:
            #print currentNode
            if currentNode.isSphere:
                node1, node2 = currentNode.left, currentNode.right
                #check collisions between children
                check1, check2 = checkCollisionSingleQuery(node1, node2)
                if check1 == 'collision' and check2 == 'collision':
                    return True
                elif check1 == 'subdivide' and check2 == 'subdivide':
                    listOfCollisionsToCheck.append((node1, node2))
                listOfCollisionsToCheck.append(node1)
                listOfCollisionsToCheck.append(node2)
            #else, do nothing as we have a single line
    return False

def checkCollisionSingleQuery(node1, node2):
    if node1.isSphere:
        if node2.isSphere:
            if (node1.radius + node2.radius) > np.linalg.norm(node1.center - node2.center):
                return 'subdivide', 'subdivide'
            else:
                return None, None
        else:
            return 'subdivide', 'ignore'
    else:
        if node2.isSphere:
            return 'ignore', 'subdivide'
        else:
            #check to see if they are contiguous edges
            #if (np.allclose(node1.curveSegments[0][1], node2.curveSegments[0][0])) \
                  #or (np.allclose(node1.curveSegments[0][0], node2.curveSegments[0][1])):
            if math.fabs(node1.index - node2.index) <= 1:
                return None, None
            else:
                A1 = node1.curveSegments[0][1] - node1.curveSegments[0][0]
                A2 = -node2.curveSegments[0][1] + node2.curveSegments[0][0]
                A = np.array([A1, A2]).T
                b = -node1.curveSegments[0][0] + node2.curveSegments[0][0] 
                t, s = np.dot(np.linalg.pinv(A), b)
                if (t < 0.0) or (t > 1.0) or (s < 0.0) or (s > 1.0):
                    return None, None
                else:
                    #curves = []
                    #curve = vis.curve(color=vis.color.red)
                    #curve.append(node1.curveSegments[0][0])
                    #curve.append(node1.curveSegments[0][1])
                    #curves.append(curve)
                    #curve = vis.curve(color=vis.color.red)
                    #curve.append(node2.curveSegments[0][0])
                    #curve.append(node2.curveSegments[0][1])
                    #curves.append(curve)
                    #print 't, s: ', t, s
                    #print 'curve1: ', node1.curveSegments[0][0], node1.curveSegments[0][1]
                    #print 'curve2: ', node2.curveSegments[0][0], node2.curveSegments[0][1]
                    #print 'collision found!'
                    return 'collision', 'collision'

def displayGraph(edges, color):
    arrows = []
    #print edges
    for p1, p2 in edges:
        arrow = vis.arrow(pos=p1.coordinates, axis=(p2.coordinates - p1.coordinates), color=color, shaftwidth=.0025)
        #curve = vis.curve(color=color)
        #print 'visual object coordinates: ', p1.coordinates
        #time.sleep(4)
        #curve.append([float(p1[0].coordinates), float(p1[1].coordinates), float(p1[2].coordinates)])
        #curve.append([float(p2[0].coordinates), float(p2[1].coordinates), float(p2[2].coordinates)])
        arrows.append(arrow)
    return arrows

def displayGraphB(edges, color):
    arrows = []
    #print edges
    for p1, p2 in edges:
        arrow = vis.arrow(pos=p1.b, axis=(p2.b - p1.b), color=color, shaftwidth=.0025)
        #curve = vis.curve(color=color)
        #print 'visual object coordinates: ', p1.coordinates
        #time.sleep(4)
        #curve.append([float(p1[0].coordinates), float(p1[1].coordinates), float(p1[2].coordinates)])
        #curve.append([float(p2[0].coordinates), float(p2[1].coordinates), float(p2[2].coordinates)])
        arrows.append(arrow)
    return arrows

def conglomeratePath(firstHalfPath, firstTree, secondHalfPath, secondTree):
    fullPath = []
    if firstHalfPath:
        for prevNode, node in zip(firstHalfPath[:-1], firstHalfPath[1:]):
            fullPath.append(prevNode)
            for intNode in firstTree.graph[prevNode][node]['intermediateConfigurations']:
                #print intNode
                fullPath.append(intNode)
        fullPath.append(firstHalfPath[-1])
    secondPath = []
    if secondHalfPath:
        for prevNode, node in zip(secondHalfPath[:-1], secondHalfPath[1:]):
            secondPath.append(prevNode)
            for intNode in secondTree.graph[prevNode][node]['intermediateConfigurations']:
                #print intNode
                secondPath.append(intNode)
        secondPath.reverse()
        fullPath.extend(secondPath)
    return fullPath

def conglomeratePathB(path, tree):
    fullPath = []
    if path:
        for prevNode, node in zip(path[:-1], path[1:]):
            fullPath.append(prevNode)
            for intNode in tree.graph[prevNode][node]['intermediateConfigurations']:
                #print intNode
                fullPath.append(intNode)
        fullPath.append(path[-1])
    return fullPath

def animatePath(path, graphScene, elasticaScene):
    elasticaScene.select()
    #n = len(path[0].x)
    #previousCollision = True
    for intermediateConfiguration in path:
        #elasticaScene.select()
        #tree = BoundingTree(edgeCurves[1:], n, 1, n-1)
        #collision = checkForCollision(tree)
        #print 'collision checking!: ', collision
        print 'A coordinates: ', intermediateConfiguration.a
        vcurve = extractXCurve(intermediateConfiguration)
        vcurve.visible = True
        graphScene.select()
        esphere = vis.sphere(pos=intermediateConfiguration.coordinates, color=vis.color.orange,
                radius=.01)
        time.sleep(.065)
        #if previousCollision and not collision:
            #time.sleep(1)
        #elif not previousCollision and collision:
            #time.sleep(1)
        #previousCollision = collision
        #tree.unDisplaySpheres()
        #if c_i != self.Full_Path[-1]:
        esphere.visible = False
        vcurve.visible = False
        elasticaScene.select()

class InputError(Exception):
    """Exceptions raised about input errors."""
    def __init__(self, msg):
        self.msg = msg

if __name__ == '__main__':
    scene1 = vis.display(center=[.5, .5, .5])
    scene1.autoscale = 0
    scene1.range = (.7,.7,.7)
    spherewait = vis.sphere()
    spherewait.visible = False
    scene2 = vis.display(center=[0.0, 0.0, 0.0], autocenter=False)
    scene2.autoscale = 0
    scene2.range = (1,1,1)
    spherewait = vis.sphere()
    spherewait.visible = False
    scene1.select()
    time.sleep(1)
    #a_1 = raw_input("Please input a_1_f ")

    while True:
        p1 = Point(np.random.rand(3))
        p2 = Point(np.random.rand(3))
        print 'start: ', p1
        print 'end: ', p2
        scene1.select()
        #time.sleep(.1)
        try:
            while not p1.isStable():
                p1 = Point(np.random.rand(3))
            while not p2.isStable():
                p2 = Point(np.random.rand(3))
            vis.sphere(pos=p1.b, color=vis.color.green, radius=.03)
            vis.sphere(pos=p2.b, color=vis.color.blue, radius=.03)
            fullPath = plannerB(p1, p2.b, .35, 1000)
            print 'path found'
            animatePath(fullPath, scene1, scene2)
            time.sleep(2)
            for obj in scene1.objects:
                obj.visible = False
            for obj in scene2.objects:
                obj.visible = False
        except InputError:
            for obj in scene1.objects:
                obj.visible = False
            for obj in scene2.objects:
                obj.visible = False
            continue

