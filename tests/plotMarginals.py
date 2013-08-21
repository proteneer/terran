import pylab as pl
import numpy as np

def loadData(filename):
    for line in csv.reader(open(filename), delimiter=' '):
        if line:
            yield line


# plot the marginals
for i in range(7):

	# import the data
	X = []

	filename = 'X' + str(i)

	for data in loadData(filename):
		X.append([float(x) for x in data])

	curves = []

	for d in range(2):
		filename = 'curve'+str(d)
		curve = []
		for data in loadData(filename):
			curve.append(float(data))
		curves.append(curve)


