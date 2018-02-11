#-------------------------------------------------------------------------------
# Name:        Point in Polygon coursework
# Purpose:
#
# Author:      Woravich Kumthonkittikul
#
# Created:     07/12/2017
# Copyright:   (c) Y520 2017
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#import library
import math
import numpy
import csv
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

#define needed general functions, which do not have to be involved with object. So every segment of object is supported.
def pointOnLine (px,py,x1,y1,x2,y2):
    """Check whether point p(px,py) does lie on the line between (x1,y1) and (x2,y2) or not"""
    if [px,py] == [x1,y1] or [px,py] == [x2,y2]: #check if it is on the end points or not
        return True

    if abs((x2-x1)*(py-y1)-(px-x1)*(y2-y1)) < 10**-7: #check if the slope of p to an end point is equal to slope of given line or not
#using 10**-7 incase of facing with the absolute zero issue.
        if (px-x1)*(px-x2) <= 0: #check if px is between x1 or x2 or not
            return True
    return False

def lineCrossing(x1,y1,x2,y2,x3,y3,x4,y4):
    """Check wheter two lines are crossing to one another or not.
      Return 'Crossing' if they are crossing
      Return 'Touching lines' if they are touching
      Return 'Same extended lines' if they are parallel and pass through same points when they are extended
      Otherwise return 'Not Crossing'"""

    #line crossing equations
    a = (x2-x1)*(y3-y2)-(y2-y1)*(x3-x2)
    b = (x2-x1)*(y4-y2)-(y2-y1)*(x4-x2)
    c = (x4-x3)*(y1-y4)-(y4-y3)*(x1-x4)
    d = (x4-x3)*(y2-y4)-(y4-y3)*(x2-x4)
    if a*b < 0 and c*d < 0: #crossing condition
        return "Crossing"
    #10**-7
    elif abs(a) < 10**-7 and abs(b) < 10**-7 and abs(c) < 10**-7 and abs(d)<10**-7: #if a=0 & b=0 & c= 0 & d=0, the lines are the same extended lines
        return "Same extended lines"
    elif abs(a) < 10**-7 or abs(b) < 10**-7 or abs(c) < 10**-7 or abs(d)<10**-7: #if a=0 or b=0 or c= 0 or d=0, the lines are likely to be touching lines"
        #This is for solving spacial cases.
        if pointOnLine(x1,y1,x3,y3,x4,y4) or pointOnLine(x2,y2,x3,y3,x4,y4) or pointOnLine(x3,y3,x1,y1,x2,y2) or pointOnLine(x4,y4,x1,y1,x2,y2):
            #if there is any end points lie on the another lines, they will be touching lines
            return "Touching lines"
    return "Not Crossing"

class Geom(object):
    """Parrent Geometry class"""
    def getStartPoint(self):
        """Get start point"""
        return self.coords[0]
    def getEndPoint(self):
        """Get end point"""
        return self.coords[-1]
    def getNumPoints(self):
        """Get the total number of vertices(points) in the object"""
        return len(self.coords)
    def addPoint(self,point):
        """Add new end point"""
        if type(point) == numpy.ndarray: #in case the input is array, eg - Point.coords
            self.coords = [self.coords,point]
        elif isinstance(point,Point): #in case the input is a Point object, eg - Point
            self.coords = [self.coords,point.coords]
        else: #report an error if they are not array or point object
            raise NotImplementedError("An Added point should be an Point object or an array")
    def mbrINT(self,geo_object):
        """Check whether the Minimum Bounding Rectangle of this object intersects with the Minimum Bounding Rectangle of the input object or not."""
        if not(isinstance(geo_object,Geom)): #if the input is not a Geom object, report an error
            raise NotImplementedError("The input has to be Geometry Object.")
        else:
            (xMin1,xMax1,yMin1,yMax1) = self.bound #object MBR ----- the attribute 'bound' is defined below
            (xMin2,xMax2,yMin2,yMax2) = geo_object.bound #the input MBR
            if((xMax1 >= xMin2 and xMax1<=xMax2) or (xMin1 >= xMin2 and xMin1<=xMax2)) and ((yMax1 >= yMin2 and yMax1<=yMax2) or (yMin1 >= yMin2 and yMin1<=yMax2)):
                #return true if they intersect one another
                return True
            else:
                return False
    def mbrInside(self,geo_object):
        """Check whether the Minimum Bounding Rectangle of this object is inside with the Minimum Bounding Rectangle of the input object or not."""
        if not(isinstance(geo_object,Geom)):
            raise NotImplementedError("The input has to be Geometry Object.") #if the input is not a Geom object, report an error
        else:
            (xMin1,xMax1,yMin1,yMax1) = self.bound #object MBR
            (xMin2,xMax2,yMin2,yMax2) = geo_object.bound #the input MBR
            if xMax1<= xMax2 and xMin1>= xMin2 and yMax1 <= yMax2 and yMin1>= yMin2:
                # return true if object MBR is completely inside the input MBR
                return True
            else:
                return False


    #additional feature
    def inPolygon(self,polygon):
        """Check whether the object is completely inside the input polygon or not (touching boundary is considered as inside).
        Work with all Geom object (Point,Line,Polygon)"""
        if not isinstance(polygon,Polygon): #if the input is not a Polygon object, report an error
            raise NotImplementedError("The input has to be Polygon Object.")
        else:
            if not(self.mbrInside(polygon)): #if object MBR is completely inside the polygon MBR
                return False
            for ver in self.coords: #checking all vertices have to be inside the polygon
                if not polygon.containPoint(ver[0],ver[1])[0]: # --- funtion 'containPoint' is defined below
                    return False
            # chceking all segments of line and polygon are crossing the input polygon or not. if cross, the object is not completely inside the input polygon.
            if isinstance(self,Point):
                return True
            elif isinstance(self,Line):
                nSeg = len(self.coords)-1
            else:
                nSeg = len(self.coords)
            for i in range(nSeg):
                st = i%len(self.coords)
                en = (i+1)%len(self.coords)
                for j in range(len(polygon.coords)):
                    st_poly = j%len(polygon.coords)
                    en_poly = (j+1)%len(polygon.coords)
                    cross=lineCrossing(self.coords[st,0],self.coords[st,1],self.coords[en,0],self.coords[en,1],polygon.coords[st_poly,0],polygon.coords[st_poly,1],polygon.coords[en_poly,0],polygon.coords[en_poly,1])
                    if cross == "Crossing":
                        return False
            return True




    @property
    def maxX(self):
        """Maximum X."""
        return max(self.coords[:,0])
    @property
    def maxY(self):
        """Maximum Y."""
        return max(self.coords[:,1])
    @property
    def minX(self):
        """Minimum X."""
        return min(self.coords[:,0])
    @property
    def minY(self):
        """Minimum Y."""
        return min(self.coords[:,1])
    @property
    def bound(self):
        """The Minimum Bounding Rectangle of the object, formatted as (min X, max X, min Y, max Y)"""
        return (self.minX,self.maxX,self.minY,self.maxY)

class Point(Geom):
    """ A simple point class"""
    def __init__(self, x=0.0, y =0.0, z=float('nan')):
        self.__coords=numpy.array([x,y,z],dtype=float)
        self.__coords.shape = (1,3)

    @property
    def coords(self):
        """ The array coordinate of the point"""
        return self.__coords
    @property
    def x(self):
        """x attribute"""
        return self.__coords[0,0]
    @property
    def y(self):
        """y attribute"""
        return self.__coords[0,1]
    @property
    def z(self):
        """z attribute"""
        return self.__coords[0,2]

    def addPoint(self,point): #overwrite the parent method
        """Not work for point object"""
        return "Can't add a point to a point"

    def getX(self):
        """Return X"""
        return self.__coords[0,0]
    def getY(self):
        """Return Y"""
        return self.__coords[0,1]
    def getZ(self):
        """Return Z"""
        return self.__coords[0,2]

    def setX(self,x):
        """set x"""
        self.__coords[0,0] = x
    def setY(self,y):
        """set y"""
        self.__coords[0,1] = y
    def setZ(self,z):
        """set z"""
        self.__coords[0,2] = z

    # This is for PIP problem could identified points at boudary
    def inPolygon2(self,polygon):
        """Check wheter the point is inside the input polygon or not.
      Return True,'Inside' if it is inside
      Return True,'Boundary' if it is on the boudary of the input polygon
      Return False, 'Outside' if it is outside
      """
        if not isinstance(polygon,Polygon): #if the input is not a Polygon object, report an error
            raise NotImplementedError("The input has to be Polygon Object.")
        return polygon.containPoint(self.x,self.y) #using polygon method, which is defined below




class Line(Geom):
    """ A simple line class"""
    def __init__(self,points = []):
        if (type(points) == list): #the input of this class can be a list of array eg. [Point.coords,Point2.coords,Point3.coords], a list of point objects eg.[Point,Point2,Point3], a list of list eg. [[1,1],[2,3],[5,10]]
            if all(type(pt) == numpy.ndarray and (pt.shape[-1]==3 or pt.shape[-1]==2) for pt in points):
                coords = numpy.array(numpy.vstack(points),dtype=float)
            elif all(type(pt) == list and (len(pt)==3 or len(pt)==2) for pt in points):
                coords = numpy.array(numpy.vstack(points),dtype=float)
            elif all(isinstance(pt,Point) for pt in points):
                coords = numpy.array(numpy.vstack([pt.coords for pt in points]),dtype=float)
            else :
                raise NotImplementedError("Line class needs a list of arrays or a list of Point objects")
        else:
            raise NotImplementedError("Line class needs a list of arrays or a list of Point objects")
        if len(coords) < 2: #the input should have at least 2 points
            raise NotImplementedError("Line class needs at least two points")
        if coords.shape[-1]<3: #set z to nan if no z value
            temp = coords
            coords = numpy.zeros((len(coords),3))
            coords[:,:-1] = temp
            coords[:,-1] = [float("nan")]*len(coords)
        self.__coords = coords
        self.__coords.shape = (len(self.__coords),3)

    @property
    def coords(self):
        """coords array of line"""
        return self.__coords


    @coords.setter
    def coords(self,points): #same as __init__
        """set coords array of line"""
        if (type(points) == list):
            if all(type(pt) == numpy.ndarray and (pt.shape[-1]==3 or pt.shape[-1]==2) for pt in points):
                coords = numpy.array(numpy.vstack(points),dtype=float)
            elif all(type(pt) == list and (len(pt)==3 or len(pt)==2) for pt in points):
                coords = numpy.array(numpy.vstack(points),dtype=float)
            elif all(isinstance(pt,Point) for pt in points):
                coords = numpy.array(numpy.vstack([pt.coords for pt in points]),dtype=float)
            else :
                raise NotImplementedError("Line class needs a list of arrays or a list of Point objects")
        else:
            raise NotImplementedError("Line class needs a list of arrays or a list of Point objects")
        if len(coords) < 2:
            raise NotImplementedError("Line class needs at least two points")
        if coords.shape[-1]<3:
            temp = coords
            coords = numpy.zeros((len(coords),3))
            coords[:,:-1] = temp
            coords[:,-1] = [float("nan")]*len(coords)
        self.__coords = coords
        self.__coords.shape = (len(self.__coords),3)

class Polygon(Geom):
    """ A simple polygon class"""
    def __init__(self,points = []): #similar to line
        if (type(points) == list):
            if all(type(pt) == numpy.ndarray and (pt.shape[-1]==3 or pt.shape[-1]==2) for pt in points):
                coords = numpy.array(numpy.vstack(points),dtype=float)
            elif all(type(pt) == list and (len(pt)==3 or len(pt)==2) for pt in points):
                coords = numpy.array(numpy.vstack(points),dtype=float)
            elif all(isinstance(pt,Point) for pt in points):
                coords = numpy.array(numpy.vstack([pt.coords for pt in points]),dtype=float)
            else :
                raise NotImplementedError("Polygon class needs a list of arrays or a list of Point objects")
        else:
            raise NotImplementedError("Polygon class needs a list of arrays or a list of Point objects")
        if len(coords) < 3:
            raise NotImplementedError("Polygon class needs at least three points")

        if coords.shape[-1]<3: #the input should have at least 3 points
            temp = coords
            coords = numpy.zeros((len(coords),3))
            coords[:,:-1] = temp
            coords[:,-1] = [float("nan")]*len(coords)
        self.__coords = coords
        self.__coords.shape = (len(self.__coords),3)

    @property
    def coords(self):
        """coords array of polygon"""
        return self.__coords

    @coords.setter
    def coords(self,points): #same as __init__
        """set coords array of polygon"""
        if (type(points) == list):
            if all(type(pt) == numpy.ndarray and (pt.shape[-1]==3 or pt.shape[-1]==2) for pt in points):
                coords = numpy.array(numpy.vstack(points),dtype=float)
            elif all(type(pt) == list and (len(pt)==3 or len(pt)==2) for pt in points):
                coords = numpy.array(numpy.vstack(points),dtype=float)
            elif all(isinstance(pt,Point) for pt in points):
                coords = numpy.array(numpy.vstack([pt.coords for pt in points]),dtype=float)
            else :
                raise NotImplementedError("Polygon class needs a list of arrays or a list of Point objects")
        else:
            raise NotImplementedError("Polygon class needs a list of arrays or a list of Point objects")
        if len(coords) < 3:
            raise NotImplementedError("Polygon class needs at least three points")
        if coords.shape[-1]<3:
            temp = coords
            coords = numpy.zeros((len(coords),3))
            coords[:,:-1] = temp
            coords[:,-1] = [float("nan")]*len(coords)
        self.__coords.shape = coords
        self.__coords.shape = (len(self.__coords),3)

    def getEndPoint(self):
        """Return end point(equal to start point)"""
        return self.coords[0]
    #for PIP prpblem work with a pair of coordinates
    def containPoint(self,px,py):
        """Check wheter the input point p(px,py) is inside the polygon or not.
      Return True,'Inside' if it is inside
      Return True,'Boundary' if it is on the boudary of the polygon
      Return False, 'Outside' if it is outside
      """
        if not (self.minX <= px <= self.maxX and self.minY <= py <= self.maxX): #if outside MBR it's outside
            return (False, "Outside")
        else:
            #Ray casting algorithm
            countInt = 0
            endRayX, endRayY = px, self.minY-1
            endpointList = numpy.vstack([self.coords[1:],self.coords[0]])
            for i,start,end in zip(range(self.getNumPoints()), self.coords, endpointList):
                cross = lineCrossing(px,py,endRayX,endRayY,start[0],start[1],end[0],end[1])
                if cross == "Not Crossing":
                    continue
                elif cross == "Crossing":
                    countInt += 1
                elif cross == "Touching lines": #spacial case
                    #check the Ray if it cross boudary at the vertex or just pass through vertex.
                    if px == end[0]:
                        if py != end[1]:
                            if i != self.getNumPoints()-1:
                                nextX = endpointList[i+1,0]
                            else:
                                nextX = endpointList[0,0]
                            if (px-start[0])*(px-nextX) < 0:
                                countInt += 1
                        else:
                            return (True,"Boundary")
                    elif px != start[0]:
                        return (True,"Boundary")
                else:
                    if (py-start[1])*(py-end[1]) <= 0:
                        return (True,"Boundary")
            if countInt%2 == 0 :
                return (False, "Outside")
            else:
                return (True, "Inside")


def plotPIPbyCSV(csvPoint,csvPolygon):
    """Plot the point and polygon from csv files regarding to the PIP problem. And return a list of Point objects
        The fuction work for CSV files with or without header.
        x and y coordinates need to be stored at first and second column respectively
    """
    poly_pts = []
    with open(csvPolygon) as csvfile: #load CVS
        datareader = csv.reader(csvfile)
        for i,row in enumerate(datareader):
            try: # use only number so can cut the header
                float(row[0])
                float(row[1])
            except ValueError:
                print "row %d:"%(i) ,row[0],row[1],"is not included into a polygon"
                continue
            poly_pts.append(Point(row[0],row[1])) #create a list of vertices(Point objects)

    poly1 = Polygon(poly_pts) # create polygon
    pts = []
    pts_X = []
    pts_Y = []
    ans = []
    point_outside = []
    point_inside = []
    point_bound = []
    #add point to list and classify point
    with open(csvPoint) as csvfile:
        datareader = csv.reader(csvfile)
        for i,row in enumerate(datareader):
            try: # use only number so can cut the header
                float(row[0])
                float(row[1])
            except ValueError:
                print "row %d:"%(i) ,row[0],row[1],"is not included into a set of points"
                continue
            pts.append(Point(row[0],row[1]))
            pts_X.append(float(row[0]))
            pts_Y.append(float(row[1]))
            pt = Point(row[0],row[1])
            if pt.inPolygon2(poly1)[1] == "Outside":
                point_outside.append(pt.coords)
            elif pt.inPolygon2(poly1)[1] == "Inside":
                point_inside.append(pt.coords)
            else:
                point_bound.append(pt.coords)

    point_outside = numpy.array(numpy.vstack(point_outside))
    point_inside = numpy.array(numpy.vstack(point_inside))
    point_bound = numpy.array(numpy.vstack(point_bound))

    #setting plotting environment
    polygon_mbr = [[poly1.maxX,poly1.maxY],[poly1.maxX,poly1.minY],[poly1.minX,poly1.minY],[poly1.minX,poly1.maxY]]

    fig1= plt.figure()
    ax1 = fig1.add_subplot(111,aspect="equal")
    polygon = plt.Polygon(poly1.coords[:,0:2], fc="#9595bc")
    ax1.add_patch(polygon)
    ax1.add_patch(plt.Polygon(polygon_mbr,fill = False))
    plt.plot(point_inside[:,0],point_inside[:,1],'bo', label = "INSIDE")
    plt.plot(point_outside[:,0],point_outside[:,1], "ro", label = "OUTSIDE")
    plt.plot(point_bound[:,0], point_bound[:,1], "yo", label = "BOUNDARY")
    plt.title("Point in Polygon")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    axes = plt.gca()
    width = axes.get_xlim()[1]-axes.get_xlim()[0]
    height = axes.get_ylim()[1]-axes.get_ylim()[0]
    plt.xlim((axes.get_xlim()[0]-width*0.15,axes.get_xlim()[1]+width*0.45 ))
    plt.ylim((axes.get_ylim()[0]-height*0.15,axes.get_ylim()[1]+height*0.45 ))
    plt.legend(numpoints = 1)
    plt.grid(alpha=0.7, linestyle="dashed")
    plt.show()
    return pts



#Testing
csvPoint = r'E:\Work\GIS\T1_GIS Principles and Technology\Python Coursework\testData\testPoints.csv' #change path
csvPolygon = r'E:\Work\GIS\T1_GIS Principles and Technology\Python Coursework\testData\testPoly.csv' #change path
plotPIPbyCSV(csvPoint,csvPolygon)





#Additional Features

class PointBuilder:
    """Interactive plotting class"""
    def __init__(self, inside,outside,boundary,poly):
        if isinstance(inside,patches.Polygon):
            self.inside = inside
            self.newPoly = inside
            self.xy = self.newPoly.get_xy()
            self.pt = outside
            self.x = list(outside.get_xdata())
            self.y = list(outside.get_ydata())
        else:
            self.inside = inside
            self.x_in = list(inside.get_xdata())
            self.y_in = list(inside.get_ydata())
            self.outside = outside
            self.x_out = list(outside.get_xdata())
            self.y_out = list(outside.get_ydata())
            self.boundary = boundary
            self.x_bound = list(boundary.get_xdata())
            self.y_bound = list(boundary.get_ydata())
        self.pt_list = []
        self.poly = poly

        self.cid = inside.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        #print('click', event)

        print event.xdata, event.ydata, " has been added to the list"
        self.pt_list.append(Point(event.xdata,event.ydata))
        if isinstance(self.inside,patches.Polygon):
            if event.inaxes!=self.newPoly.axes: return
            self.x.append(event.xdata)
            self.y.append(event.ydata)
            self.pt.set_data(self.x,self.y)
            self.xy[len(self.pt_list)-1,0] = event.xdata
            self.xy[len(self.pt_list)-1,1] = event.ydata
            self.newPoly.set_xy(self.xy)
            self.xy = self.newPoly.get_xy()
            if len(self.pt_list) <= 3:
                self.newPoly.figure.canvas.draw()
            self.pt.figure.canvas.draw()

        else:
            if event.inaxes!=self.inside.axes: return
            if self.pt_list[-1].inPolygon2(self.poly)[1] == "Inside":
                self.x_in = list(self.inside.get_xdata())
                self.y_in = list(self.inside.get_ydata())
                self.x_in.append(event.xdata)
                self.y_in.append(event.ydata)
                self.inside.set_data(self.x_in, self.y_in)
                self.inside.figure.canvas.draw()
            elif self.pt_list[-1].inPolygon2(self.poly)[1] == "Outside":
                self.x_out = list(self.outside.get_xdata())
                self.y_out = list(self.outside.get_ydata())
                self.x_out.append(event.xdata)
                self.y_out.append(event.ydata)
                self.outside.set_data(self.x_out, self.y_out)
                self.outside.figure.canvas.draw()
            else:
                self.x_bound = list(self.boundary.get_xdata())
                self.y_bound = list(self.boundary.get_ydata())
                self.x_bound.append(event.xdata)
                self.y_bound.append(event.ydata)
                self.boundary.set_data(self.x_bound, self.y_bound)
                self.boundary.figure.canvas.draw()

def createPointsByMouse(extent):
    """Plot a set of points created by clicking to the graph and input polygon, according to PIP problem. And return a list of Point objects from a mouse cursor.
     Input a extent as a polygon object or a list of the MBR (in formatt [min X, max X, min Y, max Y]).
"""
    if isinstance(extent,Polygon):
        (minX,maxX,minY,maxY) = extent.bound
        polygon = extent
    elif len(extent) == 4 and all(type(arg) in [float,int] for arg in extent) :
        (minX,maxX,minY,maxY) = extent
        polygon = Polygon([[maxX,maxY],[maxX,minY],[minX,minY],[minX,maxY]])
    else:
        raise NotImplementedError("The extent input must be a polygon object or a list of minimum bounding box(in format [x min, x max, y min, y max]).")

    polygon_mbr = [[polygon.maxX,polygon.maxY],[polygon.maxX,polygon.minY],[polygon.minX,polygon.minY],[polygon.minX,polygon.maxY]]
    fig = plt.figure()
    ax = fig.add_subplot(111,aspect="equal")
    ax.set_title('click to on the map to create point')
    poly = plt.Polygon(polygon.coords[:,0:2], fc="#9595bc")
    ax.add_patch(poly)
    #ax.add_patch(plt.Polygon(polygon_mbr,fill = False))
    inside, = ax.plot([], [],'bo', label = "INSIDE")
    outside, = ax.plot([], [],'ro', label = "OUTSIDE")  # empty line
    boundary, = ax.plot([], [],'yo', label = "BOUNDARY")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.legend(numpoints = 1)
    plt.grid(alpha=0.7, linestyle="dashed")
    axes = plt.gca()
    width = axes.get_xlim()[1]-axes.get_xlim()[0]
    height = axes.get_ylim()[1]-axes.get_ylim()[0]
    plt.xlim((axes.get_xlim()[0]-width*0.15,axes.get_xlim()[1]+width*0.45 ))
    plt.ylim((axes.get_ylim()[0]-height*0.15,axes.get_ylim()[1]+height*0.45 ))
    pointbuilder = PointBuilder(inside,outside,boundary,polygon)
    plt.show()
    return pointbuilder.pt_list

def createLineByMouse(extent):
    """Plot a line created by clicking to the graph. And return a Line object.
     Input a extent as a polygon object or a list of the MBR coordinates(in formatt [min X, max X, min Y, max Y]).
"""
    if isinstance(extent,Polygon):
        (minX,maxX,minY,maxY) = extent.bound
    elif len(extent) == 4 and all(type(arg) in [float,int] for arg in extent) :
        (minX,maxX,minY,maxY) = extent
    else:
        raise NotImplementedError("The extent input must be a polygon object or a list of minimum bounding box(in format [x min, x max, y min, y max]).")

    fig = plt.figure()
    ax = fig.add_subplot(111,aspect="equal")
    ax.set_title('click to add your line, close the window when it is finished')
    line, = ax.plot([], [],'-bo', label = "Line")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.legend(numpoints = 1)
    plt.grid(alpha=0.7, linestyle="dashed")
    plt.xlim((minX,maxX))
    plt.ylim((minY,maxY))
    polygon = Polygon([[maxX,maxY],[maxX,minY],[minX,minY],[minX,maxY]])
    pointbuilder = PointBuilder(line,line,line,polygon)
    plt.show()
    if len(pointbuilder.pt_list) < 2:
        return False
    else:
        return Line(pointbuilder.pt_list)

def createPolygonByMouse(extent):
    """Plot a polygon created by clicking to the graph. And return a Line object.
     Input a extent as a polygon object or a list of the MBR coordinates (in formatt [min X, max X, min Y, max Y]).
"""
    if isinstance(extent,Polygon):
        (minX,maxX,minY,maxY) = extent.bound
    elif len(extent) == 4 and all(type(arg) in [float,int] for arg in extent) :
        (minX,maxX,minY,maxY) = extent
    else:
        raise NotImplementedError("The extent input must be a polygon object or a list of minimum bounding box(in format [x min, x max, y min, y max]).")

    fig = plt.figure()
    ax = fig.add_subplot(111,aspect="equal")
    ax.set_title('click to add your polygon, close the window when it is finished')
    p = plt.Polygon([[0,0],[0,0]], fc="#9595bc", label = "Polygon")
    poly = ax.add_patch(p)
    pt, = ax.plot([], [],'-bo')
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.legend(numpoints = 1)
    plt.grid(alpha=0.7, linestyle="dashed")
    plt.xlim((minX,maxX))
    plt.ylim((minY,maxY))
    polygon = Polygon([[maxX,maxY],[maxX,minY],[minX,minY],[minX,maxY]])
    pointbuilder = PointBuilder(poly,pt,poly,polygon)
    plt.show()
    if len(pointbuilder.pt_list) < 3:
        return False
    else:
        return Polygon(pointbuilder.pt_list)


#examples use of additional funtions
poly = createPolygonByMouse([0,100,0,100])
line = createLineByMouse([0,100,0,100])

try:
    pt_list = createPointsByMouse(poly)
except ValueError:
    pt_list = createPointsByMouse([1,100,1,100])

#testing inPolygon attribute
try:
    print "Your just created line is inside your just created polygon: ", line.inPolygon(poly)
except ValueError:
    pass
