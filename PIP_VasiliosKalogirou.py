##This program uses the ray-casting algorithm the check if a point is inside, outside, or on the boundaries of the polygon.
##The polygon is formed by a set of coordinates that are imported as a csv file.
##The points to be tested are imported from a csv file and then plotted using the matplotlib library.
##The program can take a random point as an input point to test, but it does not plot it.
##The program contains helpful comments to explain how it was built and handles almost all the mathematical errors apart from the case when the line segment 
##    between two consecutive vertices is vertical and its slope is indefinite. A warning message is showed in this case by the program.
##The user needs to place the two files in a folder and set the directory and the names of the files manually.
##This program does not handle possible errors inside the files so the user must make sure all the coordinates are in the correct type and format.

#Step 1: Import the necessary libraries
import numpy as np #import numpy library 
import math # import math library
import matplotlib.pyplot as plt #import matplotlib library
import os # import os library
import csv #import csv library to read and write csv files

#Step 2: Create the necessary classes that will be used in this program.
class Geom(object):    # Geometry class - inherits from the object class. Points, lines and polygons all inherit from this class.
    def getStartPoint(self): 
        return self.coords[0]
    def getEndPoint(self):
        return self.coords[-1]
    def getNumPoints(self):
        return self.coords.shape[0]
    def addPoint(self, point):
        self.coords = np.vstack([self.coords, point])
    @property
    def minX(self):#This property returns the minimum X coordinate of every instance that inherits from the class Geom. Applies to points, lines and polygons.
        return min(self.coords[:,0])
    @property
    def maxX(self):#This property returns the maximum X coordinate of every instance that inherits from the class Geom. Applies to points, lines and polygons.
        return max(self.coords[:,0])
    @property
    def minY(self):#This property returns the minimum Y coordinate of every instance that inherits from the class Geom. Applies to points, lines and polygons.
        return min(self.coords[:,1])
    @property
    def maxY(self):#This property returns the maximum Y coordinate of every instance that inherits from the class Geom. Applies to points, lines and polygons.
        return max(self.coords[:,1])

class Point(Geom):    # Point class
    """A simple point class"""
    def __init__(self, x=0, y=0, z=float('nan')): #This is the constructor of the point class.
        self.__coords = np.array([x,y,z], dtype=float)
        self.__coords.shape = (1,3)
    @property
    def x(self):
        return self.__coords[0,0]
    @property
    def y(self):
        return self.__coords[0,1]
    @property
    def z(self):
        return self.__coords[0,2]
    @x.setter
    def x(self, x):
        self.__coords[0,0] = x
    @y.setter   
    def y(self, y):
        self.__coords[0,1] = y
    @z.setter
    def z(self, z):
        self.__coords[0,2] = z
    @property
    def coords(self):
        return self.__coords
    def addPoint(self, point):
        return "Can't add a point to a point"

class Line(Geom):   # Line class
    def __init__(self, points = []): #This is the constructor of the line class
        self.__coords = np.vstack(points)
    @property
    def coords(self):
        return self.__coords
    @coords.setter
    def coords(self, points):
        self.__coords = np.vstack(points)

class Polygon(Line, Geom): # Polygon class
    def getEndPoint(self): #This is the constructor of the line class
        return self.getStartPoint()
    def testMBR(self, testPoint): #Function that tests if the point is inside the MBR
        if (testPoint.x < self.minX) or (testPoint.x > self.maxX) or (testPoint.y<self.minY) or (testPoint.y>self.maxY):
            return 0
        else:
            #print "The Point is inside the Minimum Bounding Rectangle"
            return 1
    def testPIP(self, testPoint):
        specialCase1=False #Checks if the testing point lies on a vertex of the polygon.
        specialCase2=False #Checks if the testing point lies on an edge of the polygon.
        test = False #test is a boolean variables that switches from false to true every time the two testing lines intersect. 
        xA = testPoint.x #assign the x value of the testing point to the x coordinate of point A
        yA = testPoint.y #assign the y value of the testing point to the y coordinate of point A
        xB = self.maxX + 10 #assign a value x which is outside of the MBR to point B
        yB = self.maxY + 10 #assign a value y which is outside of the MBR to point B

        #Check the intersections with the line segments that form the polygon. -1 is used not to take into account the last point which is the first one as well.
        for i in range(0,self.getNumPoints()-1): 
            if (testPoint.x==self.coords[i][0]) and (testPoint.y==self.coords[i][1]):
                specialCase1=True
                #print "Special Case 1: The testing point lies on a vertex"
                return "boundary"
                break
            else:
                xC = self.coords[i][0] #The line is formed by two vertices. C is the first vertex of the line. Assign its X coordinate.
                yC = self.coords[i][1] #The line is formed by two vertices. C is the first vertex of the line. Assign its Y coordinate.
                xD = self.coords[i+1][0] #The line is formed by two vertices. D is the second vertex of the line. Assign its X coordinate.
                yD = self.coords[i+1][1] #The line is formed by two vertices. D is the second vertex of the line. Assign its Y coordinate.
                lp1 = Point(xC,yC) #Create point C - essentially this is the first vertex of a line segment in the polygon
                lp2 = Point(xD,yD) #Create point D - essentially this is the second vertex of a line segment in the polygon
                ls = Line([lp1.coords,lp2.coords]) #Create the Line Segment which is formed by two consecutive points, C and D
                #We first check if the line is vertical or horizontal so that we don't calculate a slope value of zero or infity that may cause errors.
                #Check if the line segment is vertical: slope = infinity AND the point lies on the line segment
                if (xA==xC and xA==xD and yA>ls.minY and yA<ls.maxY):
                    specialCase2=True
                    #print "Special Case 2: The point lies on an edge. This line segment is a vertical line."
                    return "boundary"
                    break
                #Check if the line segment is horizontal: slope = infinity AND the point lies on the line segment
                elif (yA==yC and yA==yD and xA>ls.minX and xA<ls.maxX):
                    specialCase2=True
                    #print "Special Case 2: The point lies on an edge. This line segment is a horizontal line."
                    return "boundary"
                    break
                else:
                    m=(yD-yC)/(xD-xC) #Calculate the slope of the line which is formed by two consecutive vertices
                    beta=yD-m*xD #Calculate the beta parameter of the line which is formed by two consecutive vertices
                    if ((yA-(m*xA)-beta)==0 and xA>ls.minX and xA<ls.maxX and yA>ls.minY and yA<ls.maxY): #Check if the point lies on an edge.
                        specialCase2=True
                        #print "Special Case 2: The point lies on an edge"
                        return "boundary"
                        break
                    else: #The following parameters are used to implement the ray-casting algorithm.
                        param1 = (xB-xA)*(yC-yB)-(yB-yA)*(xC-xB)
                        param2 = (xB-xA)*(yD-yB)-(yB-yA)*(xD-xB)
                        param3 = (xD-xC)*(yA-yD)-(yD-yC)*(xA-xD)
                        param4 = (xD-xC)*(yB-yD)-(yD-yC)*(xB-xD)
                        if (param1*param2 <0) and (param3*param4<0):
                            test=not(test) #Switch the status of the test boolean.
        if specialCase1==False and specialCase2==False: #If the point is not on a vertex or on a boundary then check if it is inside or outside
            if test==True:
                return "inside"
            else:
                return "outside"

#Step 3: Set the directory for the csv files. This is a variable that the user must change.
os.chdir(r'C:\Users\Vasilis\Dropbox\5_UCL\1st_Semester\2.GIS Principles and Technology\PointInPolygon\testData')

#Step 4: Read the coordinates of the polygon
with open ('testPoly2.csv') as csvfile: #Open the csv file that contains the coordinates of the polygon.
    dataReader = csv.reader(csvfile) #dataReader is a 'csv.reader' variable used to read what is inside the csv file.
    rList= dataReader.next() #Create a random list to read the coordinates of the first point. rList is a temporary variable
    startPoint=Point(x=rList[0], y =rList[1]) #Create the first point which will be used to create the polygon. startPoint is a temporary variable
    myPolygon=Polygon(startPoint.coords) #create the polygon using the first point
    for line in dataReader:
        myPoint = Point(x=line[0],y=line[1])
        myPolygon.addPoint(myPoint.coords) #Add the point to the polygon using the addPoint method.
    myPolygon.addPoint(myPolygon.getEndPoint()) #Close the polygon.
    del rList, startPoint #delete the temporary variables, dataReader
    
#Step 5: Create these 6 lists to append the coordinates of the points in them.
    #These will be used later to plot the points differently according to their status.
outX=[]
outY=[]
inX=[]
inY=[]
boundX=[]
boundY=[]
#Step 6: Read the points from the csv file given and test if they are inside, outside, or on the boundaries of the polygon.
with open ('testPoints.csv') as testingPoints: #Open the csv file that contains the coordinates of the polygon.
    s=1 #Index used to show which point is being tested.
    dataReader = csv.reader(testingPoints) #dataReader is a 'csv.reader' variable used to read what is inside the csv file.
    r = dataReader.next() #r is a temporal variable used to ignore the header of the file.
    for line in dataReader:
        tpoint=Point(x=line[0],y=line[1])
        #First check if the point is inside the Minimum Bounding Rectangle
        if myPolygon.testMBR(tpoint) == 0:
            print "Point" + str(s) + ": " + "is outside of the polygon (outside of the MBR)"
            outX.append(tpoint.x) #Append the x coordinate of the point in the outX list.
            outY.append(tpoint.y) #Append the y coordinate of the point in the outY list.
        else:
                if myPolygon.testPIP(tpoint) == "boundary":
                    print "Point" + str(s) + ": " + "lies on the boundary of the polygon"
                    boundX.append(tpoint.x) #Append the x coordinate of the point in the boundX list.
                    boundY.append(tpoint.y) #Append the y coordinate of the point in the boundY list.
                elif myPolygon.testPIP(tpoint) == "inside":
                    print "Point" + str(s) + ": " + "is inside the polygon"
                    inX.append(tpoint.x) #Append the x coordinate of the point in the inX list.
                    inY.append(tpoint.y) #Append the y coordinate of the point in the inY list.
                else:
                    print "Point" + str(s) + ": " + "is outside of the polygon"
                    outX.append(tpoint.x) #Append the x coordinate of the point in the outX list.
                    outY.append(tpoint.y) #Append the y coordinate of the point in the outX list.
        s=s+1
    del r, dataReader, s, line #Delete the temporary variables
    
#Step 6: Test your own point if you want. Uncomment the next two lines of code, write the coordinates of the point you want to test and run the program.
p1=Point(x=6, y=9)
if myPolygon.testMBR(p1) == 0:
    print "Point entered: " + "is outside of the polygon (outside of the MBR)"
else:
    if myPolygon.testPIP(p1) == "boundary":
        print "Point entered: " + "lies on the boundary of the polygon"
    elif myPolygon.testPIP(p1) == "inside":
        print "Point entered: " + "is inside the polygon"
    else:
        print "Point entered: " + "is outside of the polygon"

#Step 7: PlotPlot 
polygonPlot,=plt.plot(myPolygon.coords[:,0],myPolygon.coords[:,1],'b-', linewidth=1) #Plot the polygon
polygonPlot,=plt.plot([myPolygon.minX, myPolygon.minX, myPolygon.maxX, myPolygon.maxX,myPolygon.minX],
                      [myPolygon.minY, myPolygon.maxY,myPolygon.maxY,myPolygon.minY,myPolygon.minY],'r-',linewidth=2) #Plot the Minimum Boundary Rectangle
polygonVertices=plt.plot(myPolygon.coords[:,0],myPolygon.coords[:,1],'bo', ms=4) #Plot the Vertices of the Polygon - optional
pointsInside, =plt.plot(inX,inY,'g^', ms=8) #Plot the points that are inside the polygon as green triangles
pointsOutside, =plt.plot(outX, outY, 'ro', ms=8) #Plot the points that are outside the polygon as red circles
boundaryPoints, = plt.plot(boundX, boundY, 'ms', ms = 8) #Plot the points that lie on the boundaries as magenta squares
plt.axis([myPolygon.minX-5,myPolygon.maxX+5,myPolygon.minY-5,myPolygon.maxY+5]) #Enlarge the extent of the plot
plt.xlabel("X") #Label the x axis
plt.ylabel("Y") #Label the y axis
plt.legend([pointsInside, pointsOutside, boundaryPoints], ["Inside", "Outside", "Boundary"], loc=4) #Put a legend
plt.show() #Show the plot
