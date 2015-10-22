"""
Bond Angle Calculation and Bond Distances
by Irenaeus Chan
"""

import sys
import string
import math

backboneAtoms = ("N ", "CA", "C ")

def dotProduct(vector1, vector2):
	return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]

def bondAngle(dotproduct, vector1, vector2):
	return math.degrees(math.acos(dotproduct/(vectorMagnitude(vector1) * vectorMagnitude(vector2))))

def vectorCalculation(coord1, coord2):
	return [coord1[0]-coord2[0], coord1[1]-coord2[1], coord1[2]-coord2[2]]

def vectorMagnitude(vector):
	return ((vector[0])**2 + (vector[1])**2 + (vector[2])**2)**0.5

def readFile(file_name):
	info = []
	with open(file_name, "r") as stream:
		for line in stream:
			if line[0:4] == "ATOM" and line[13:15] in backboneAtoms:
				info.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])
	return info

def calculateNCC(atom):
	i = 0
	with open("NCC.txt", "w") as output:
		while (i < len(atom)):
			oldVector = vectorCalculation(atom[i], atom[i+1])
			newVector = vectorCalculation(atom[i+1], atom[i+2])
			output.write(str(bondAngle(dotProduct(oldVector, newVector), oldVector, newVector))+'\n')
			i+=3

def calculateCCN(atom):
	i = 1
	with open("CCN.txt", "w") as output:
		while (i < len(atom)):
			oldVector = vectorCalculation(atom[i], atom[i+1])
			newVector = vectorCalculation(atom[i+1], atom[i+2])
			output.write(str(bondAngle(dotProduct(oldVector, newVector), oldVector, newVector))+'\n')
			i+=3

def calculateCNC(atom):
	i = 2
	with open("CNC.txt", "w") as output:
		while (i < len(atom)):
			oldVector = vectorCalculation(atom[i], atom[i+1])
			newVector = vectorCalculation(atom[i+1], atom[i+2])
			output.write(str(bondAngle(dotProduct(oldVector, newVector), oldVector, newVector))+'\n')
			i+=3

def calculateNC(atom):
	i = 0
	with open("N-C.txt", "w") as output:
		while (i < len(atom)):
			vector = vectorCalculation(atom[i], atom[i+1])
			output.write(str(vectorMagnitude(vector))+ '\n')
			i+=3

def calculateCC(atom):
	i = 1
	with open("C-C.txt", "w") as output:
		while (i < len(atom)):
			vector = vectorCalculation(atom[i], atom[i+1])
			output.write(str(vectorMagnitude(vector))+ '\n')
			i+=3

def calculateCN(atom):
	i = 2
	with open("C-N.txt", "w") as output:
		while (i < len(atom)):
			vector = vectorCalculation(atom[i], atom[i+1])
			output.write(str(vectorMagnitude(vector))+ '\n')
			i+=3


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print "ERROR: No file was provided"
		sys.exit(1)
	elif (sys.argv[1].endswith(".pdb")):
		print "Succesful..."
	else:
		print "ERROR: File type is incorrect."
		sys.exit(1)

	atom = readFile(sys.argv[1])
	
	calculateNC(atom)
	calculateCC(atom)
	calculateCN(atom)
	calculateNCC(atom)
	calculateCCN(atom)
	calculateCNC(atom)
	