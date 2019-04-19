# -----------------------------------
# Author:
# Elliot Steissberg
# ERDC-EL
# Version 2.4
# March 18, 2019
#
# Revised by Todd Steissberg, HEC
# April 15, 2019
# -----------------------
"""
The velocity path is in fort64.h5 or equivalent, the elevation path is in fort63.h5 or equivalent, 
and all others are found in the main {project name}.p04.hdf project hdf file or equivalent
"""

import sys
import os
import h5py
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy import interpolate

velocitiesPath = "Datasets/Depth-averaged Velocity (64)/Values"
timesPath = "Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Time Date Stamp"
waterSurfaceElevationsPath = "Datasets/Water Surface Elevation (63)/Values"

"""The values dictionary will contain the following keys:
Times:                  Time stamp of each output time step
X:                      East-West (X) coordinate of ground surface at the nodes
Y:                      North-South (Y) coordinate of ground surface at the nodes
Z:                      Elevation of the ground surface at the nodes
U:                      East-West velocity component at the nodes
V:                      North-South velocity component at the nodes
WaterSurfaceElevations: Water surface elevation at the nodes
Depths:                 Depths at the nodes
"""
valuesDict = {}

def quiverPlot(X, Y, U, V, timeStamp, scale):
    figure, axis = plt.subplots()
    axis.set_title(timeStamp)
    plt.axis('equal')
    M = np.hypot(U, V)
    Q = axis.quiver(X, Y, U, V, M, units='xy', scale=scale)  # pivot='tip', width=1, scale=1)
    quiverKey = axis.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
                               coordinates='figure')


def plotVelocities(projectFilename, outputPath, scale=1.0, fileExtension="png"):
    X = valuesDict["X"]
    Y = valuesDict["Y"]
    for timeStamp in valuesDict["Times"]:
        outputFilename = createOutputFilename(outputPath, projectFilename, timeStamp,
                                              fileExtension, customLabel="velocities")
        timeIndex = valuesDict["Times"].index(timeStamp)
        U = np.array(valuesDict["U"][timeIndex])
        V = np.array(valuesDict["V"][timeIndex])
        print("Saving plot to %s" % outputFilename)
        quiverPlot(X, Y, U, V, timeStamp, scale)
        plt.savefig(outputFilename)
        plt.close()

def plotDepths(projectFilename, outputPath, fileExtension="png", xmin=0.0, xmax=10000.0,
               dx=4.0, ymin=0.0, ymax=10000.0, dy=4.0, zmin=0.0, zmax=100.0, dz=1.0):
    X = valuesDict["X"]
    Y = valuesDict["Y"]

    # Prepare uniform interpolation grid
    xi = np.linspace(xmin, xmax, int(dx))
    yi = np.linspace(ymin, ymax, int(dy))
    Xi, Yi = np.meshgrid(xi, yi)

    for timeStamp in valuesDict["Times"]:
        outputFilename = createOutputFilename(outputPath, projectFilename, timeStamp,
                                              fileExtension, customLabel="depths")
        timeIndex = valuesDict["Times"].index(timeStamp)
        Z = np.array(valuesDict["Depths"][timeIndex])

        # Interpolate to grid
        interpToGrid = interpolate.interp2d(X, Y, Z, kind='cubic')
        Zi = interpToGrid(xi, yi)

        # Plot depths
        print("Saving plot to %s" % outputFilename)
        figure, axis = plt.subplots()
        axis.set_title(timeStamp)
        plt.axis('equal')
        axis.pcolor(Xi, Yi, Zi)
        plt.savefig(outputFilename)
        plt.close()


def getValues(projectFilename, fort64Filename, fort63Filename, tinFilename):
    """Get and compute values from the four specified input files"""

    # Open input files
    projectFile = h5py.File(projectFilename, "r")
    fort63File = h5py.File(fort63Filename, "r")
    fort64File = h5py.File(fort64Filename, "r")
    tinFile = open(tinFilename, "r")

    # Velocities
    velocities = fort64File[velocitiesPath]
    valuesDict["U"] = [list(np.transpose(velocities[i])[0]) for i in range(len(velocities))]
    valuesDict["V"] = [list(np.transpose(velocities[i])[1]) for i in range(len(velocities))]

    # Time
    timesBinary = list(projectFile[timesPath])
    timesAscii = list(map(lambda t: t.decode("ASCII"), timesBinary))
    valuesDict["Times"] = timesAscii

    # Elevation of ground surface
    valuesDict["WaterSurfaceElevations"] = list(fort63File[waterSurfaceElevationsPath])

    waterSurfaceElevations = valuesDict["WaterSurfaceElevations"]
    numElevationPoints = len(waterSurfaceElevations[0])

    # Assemble a list of x, y, z coordinates and depths at the nodes
    i = 1
    tinLines = tinFile.readlines()
    numPoints = int(tinLines[1].strip().split("\t")[1])
    dataLines = tinLines[2:(numPoints + 2)]
    xList = []
    yList = []
    zList = []
    nodeList = []
    for line in dataLines:
        data = line.strip().split("\t")
        i, x, y = data[0].split(" ")
        z = data[1]
        x, y, z = list(map(float, [x, y, z]))
        nodeList.append(int(i))
        xList.append(x)
        yList.append(y)
        zList.append(z)

    valuesDict["Nodes"] = np.array(nodeList)
    valuesDict["X"] = np.array(xList)
    valuesDict["Y"] = np.array(yList)
    valuesDict["Z"] = np.array(zList)

    # Create the calculated depth array for each time step
    waterDepthsList = []
    # Iterate over time steps
    for i in range(len(waterSurfaceElevations)):
        groundElevation = valuesDict["Z"][i]
        waterDepths = []
        for waterSurfaceElevation in waterSurfaceElevations[i]:
            if waterSurfaceElevation > -999.0 and groundElevation > -999.0:
                depth = waterSurfaceElevation - groundElevation
                waterDepths.append(depth)
            else:
                waterDepths.append(0.0)
        waterDepthsList.append(waterDepths)
    valuesDict["Depths"] = waterDepthsList

def createOutputFilename(outputPath, projectFilename, timeStamp, fileExtension, customLabel=""):
    pathParts = os.path.splitext(projectFilename)[0].split("/")
    baseFilename = pathParts[-1]
    timeLabel = parseDate(timeStamp).replace("/", "_").replace(" ", "_").replace(":", "")
    if customLabel != "":
        outputFilename = "%s/%s_%s_%s.%s" % (outputPath, baseFilename, timeLabel, customLabel, fileExtension)
    else:
        outputFilename = "%s/%s_%s.%s" % (outputPath, baseFilename, timeLabel, fileExtension)
    return outputFilename


def writeToAscii(projectFilename, header, dataFormatString, outputPath, fileExtension):
    """Creates an output file for each time step in TecPlot format
    from the dictionary created by getValues"""
    for timeStamp in valuesDict["Times"]:
        outputFilename = createOutputFilename(outputPath, projectFilename, timeStamp, fileExtension)
        timeIndex = valuesDict["Times"].index(timeStamp)
        outputFile = open(outputFilename, "w")
        lines = [header]
        for i in range(len(valuesDict["X"])):
            x = valuesDict["X"][i]
            y = valuesDict["Y"][i]
            z = valuesDict["Z"][i]
            u = valuesDict["U"][timeIndex][i]
            v = valuesDict["V"][timeIndex][i]
            depth = valuesDict["Depths"][timeIndex][i]
            lines.append(dataFormatString % (x, y, z, u, v, depth))
        print("Writing to %s" % outputFilename)
        outputFile.write("\n".join(lines))
        outputFile.close()


def writeToCsv(projectFilename, outputPath):
    """Creates an output file for each time step in comma-delimited format
    from the dictionary created by getValues"""
    header = "X,Y,Z,U,V,Depth"
    dataFormatString = "%f,%f,%f,%f,%f,%f"
    writeToAscii(projectFilename, header, dataFormatString, outputPath, "csv")


def writeToTecPlot(projectFilename, outputPath):
    """Creates an output file for each time step in TecPlot format
    from the dictionary created by getValues"""
    header = 'VARIABLES = "X" "Y" "Z" "U" "V" "DEPTH"'
    dataFormatString = "%f %f %f %f %f %f"
    writeToAscii(projectFilename, header, dataFormatString, outputPath, "dat")


def parseDate(timeString):
    """(Hopefully) Works on the times from main output, not from optional files.
     Note: the naming discrepancy in the datasets needs to be fixed by the RAS team."""
    months = {"JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
              "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12}

    monthStr = "undefined"
    dayStr = "undefined"
    yearStr = "undefined"

    dateStr, timeStr = str(timeString).split(" ")

    for i in months.keys():
        if i in dateStr:
            monthStr = months[i]
            dayStr = dateStr.split(i)[0]
            yearStr = dateStr.split(i)[1]

    # Create output date-time string
    # Only keep hour and minute from time string
    hourMinuteStr = ":".join(timeStr.split(":")[0:2])
    dateTime = "%s/%s/%s %s" % (dayStr, monthStr, yearStr, hourMinuteStr)
    return dateTime


if __name__ == "__main__":
    # This section gets filenames with either command-line args (semicolon-separated)
    # or with input from running it directly.

    # try:
    #     projectFilename, fort64Filename, fort63Filename, tinFilename = sys.argv[1].split(";")
    # except:
    #     projectFilename = input("Main project HDF5 file (p*.hdf): ")
    #     fort64Filename = input("ADCIRC HDF5 velocity file (fort64.h5): ")
    #     fort63Filename = input("ADCIRC HDF5 water surface elevation file (fort63.h5): ")
    #     tinFilename = input("TIN grid file: ")

    # Enter full paths to project HDF5 output file and ADCIRC files generated by RasMapper
    # projectFilename = "C:/Users/q0hectes/Documents/IdeaProjects/HDF5utilityForDaveSmith/src/Muncie.p04.hdf"
    # fort63Filename = "C:/Users/q0hectes/Documents/IdeaProjects/HDF5utilityForDaveSmith/src/fort63.h5"
    # fort64Filename = "C:/Users/q0hectes/Documents/IdeaProjects/HDF5utilityForDaveSmith/src/fort64.h5"
    # tinFilename = "C:/Users/q0hectes/Documents/IdeaProjects/HDF5utilityForDaveSmith/src/2D interior Area_TIN.grd"
    # outputPath = "C:/Users/q0hectes/Documents/IdeaProjects/HDF5utilityForDaveSmith/src"
    # vectorScale = 0.01
    projectFilename = "F:/HDF_utility_for_Dave_Smith/ModelForTesting/HEC-RAS/SF_Clearwater.p01.hdf"
    fort63Filename = "F:/HDF_utility_for_Dave_Smith/ModelForTesting/HEC-RAS/P1/fort63.h5"
    fort64Filename = "F:/HDF_utility_for_Dave_Smith/ModelForTesting/HEC-RAS/P1/fort64.h5"
    tinFilename = "F:/HDF_utility_for_Dave_Smith/ModelForTesting/HEC-RAS/P1/2D domain_TIN.grd"
    outputPath = "F:/HDF_utility_for_Dave_Smith/ModelForTesting/HEC-RAS/P1"
    vectorScale = 0.1

    # First, create the dictionary of relevant values from the read files,
    # and then calculate the depth values from the tin or grid file.
    getValues(projectFilename, fort64Filename, fort63Filename, tinFilename)

    meanDepth = np.mean(np.transpose(valuesDict["Depths"])[0])
    print(meanDepth)

    writeToTecPlot(projectFilename, outputPath)
    plotVelocities(projectFilename, outputPath, scale=vectorScale)
    # plotDepths(projectFilename, outputPath, fileExtension="png", xmin=404000.0, xmax=413000.0,
    #            dx=100.0, ymin=1800000.0, ymax=1806000.0, dy=100.0, zmin=-10.0, zmax=20.0, dz=1.0)

