#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		if self.banded:
			score, alignment1, alignment2 = self.alignBanded(seq1,seq2)
		else:
			score, alignment1, alignment2 = self.alignUnbanded(seq1,seq2)

		return {'align_cost':score, 'seqi_first100':alignment1[:101], 'seqj_first100':alignment2[:101]}

	# unrestricted algorithm
	def alignUnbanded(self, seq1, seq2):
		# get lengths
		xLength = len(seq1) + 1
		if xLength > self.MaxCharactersToAlign + 1:
			xLength = self.MaxCharactersToAlign + 1
		yLength = len(seq2) + 1
		if yLength > self.MaxCharactersToAlign + 1:
			yLength = self.MaxCharactersToAlign + 1
		# initialize matrix, fill first row and column
		matrix = [[self.Tuple(math.inf, 0) for x in range(yLength)] for y in range(xLength)]
		for i in range(xLength-1):
			matrix[i][0].value = i*5
			matrix[i][0].prev = 1
		for j in range(yLength-1):
			matrix[0][j].value = j*5
			matrix[0][j].prev = 3
		# loop through remaining and update
		for i in range(1, xLength):
			for j in range(1, yLength):
				# store potential values
				leftValue = matrix[i-1][j].value + 5
				# if match, subtract 3 is option
				if seq1[i - 1] == seq2[j - 1]:
					diagonalValue = matrix[i - 1][j - 1].value - 3
				else:
					diagonalValue = matrix[i - 1][j - 1].value + 1
				upValue = matrix[i][j - 1].value + 5
				# update with minimum possible, priority: up, diagonal, left (opposite to project description because seq1 and seq2 are switched)
				if upValue < matrix[i][j].value:
					matrix[i][j].value = upValue
					matrix[i][j].prev = 3
				if diagonalValue < matrix[i][j].value:
					matrix[i][j].value = diagonalValue
					matrix[i][j].prev = 2
				if leftValue < matrix[i][j].value:
					matrix[i][j].value = leftValue
					matrix[i][j].prev = 1
		# get alignments
		alignment1 = ""
		alignment2 = ""
		i = xLength-1
		j = yLength-1
		while i > 0 and j > 0:
			tempTuplePrev = matrix[i][j].prev
			# if left pointer
			if tempTuplePrev == 1:
				alignment1 = seq1[i-1] + alignment1
				alignment2 = "-" + alignment2
				i -= 1
			# if diagonal pointer
			if tempTuplePrev == 2:
				alignment1 = seq1[i-1] + alignment1
				alignment2 = seq2[j-1] + alignment2
				i -= 1
				j -= 1
			# if up pointer
			if tempTuplePrev == 3:
				alignment1 = "-" + alignment1
				alignment2 = seq2[j-1] + alignment2
				j -= 1

		return matrix[xLength-1][yLength-1].value, alignment1, alignment2

	# banded algorithm
	def alignBanded(self, seq1, seq2):
		# get lengths
		xLength = len(seq1) + 1
		if xLength > self.MaxCharactersToAlign + 1:
			xLength = self.MaxCharactersToAlign + 1
		yLength = len(seq2) + 1
		if yLength > self.MaxCharactersToAlign + 1:
			yLength = self.MaxCharactersToAlign + 1
		# maximum is the shorter of the two, going diagonal
		if yLength > xLength + 3:
			yLength = xLength + 3
		# initialize matrix, fill first row and column
		matrix = [[self.Tuple(math.inf, 0) for x in range(7)] for y in range(yLength)]
		# fill top row
		for j in range(6,2, -1):
			matrix[0][j].value = j*5
			matrix[0][j].prev = 1
		# fill left column
		j = 0
		for i in range(3, -1, -1):
			matrix[j][i].value = j*5
			matrix[j][i].prev = 3
			j += 1
		# set default final index
		finalXIndex = 6
		# loop
		for i in range(1, yLength):
			start = i-3
			for a in range(7):
				j = start + a  # j is index in seq, a is index in array
				# check if out of bounds
				if j < 0:
					continue
				if j > xLength-1:
					if i == yLength - 1:
						finalXIndex = a-1
						break
					else:
						continue
				# not out of bounds
				leftValue = math.inf
				if a > 0:
					leftValue = matrix[i][a-1].value + 5
				# check for match
				if seq1[j-1] == seq2[i-1]:
					diagonalValue = matrix[i-1][a].value - 3
				else:
					diagonalValue = matrix[i-1][a].value + 1
				upValue = math.inf
				if a < 6:
					upValue = matrix[i-1][a+1].value + 5

				if leftValue < matrix[i][a].value:
					matrix[i][a].value = leftValue
					matrix[i][a].prev = 1
				if diagonalValue < matrix[i][a].value:
					matrix[i][a].value = diagonalValue
					matrix[i][a].prev = 2
				if upValue < matrix[i][a].value:
					matrix[i][a].value = upValue
					matrix[i][a].prev = 3

		# get alignments
		alignment1 = ""
		alignment2 = ""
		i = yLength-1
		a = finalXIndex
		j = xLength-1
		while i > 0 and j > 0:
			tempTuplePrev = matrix[i][a].prev
			# if left pointer
			if tempTuplePrev == 1:
				alignment1 = seq1[j-1] + alignment1
				alignment2 = "-" + alignment2
				a -= 1
				j -= 1
			# if diagonal pointer
			if tempTuplePrev == 2:
				alignment1 = seq1[j-1] + alignment1
				alignment2 = seq2[i-1] + alignment2
				i -= 1
				j -= 1
			# if up pointer
			if tempTuplePrev == 3:
				alignment1 = "-" + alignment1
				alignment2 = seq2[i-1] + alignment2
				i -= 1
				a += 1

		return matrix[yLength - 1][finalXIndex].value, alignment1, alignment2


# tuple to store value and back pointer
	class Tuple:
		def __init__(self, value, prev):
			self.value = value
			self.prev = prev # 1 is left, 2 is diagonal, 3 is up, 0 is undefined
