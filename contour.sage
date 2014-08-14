import scipy.interpolate as sip
import numpy as np

indexstr = []
indexstr = raw_input('> ')
filename = "output_psi_" + indexstr + ".txt"
print(filename)
f = open(filename, "r")
A = []
tmp = []

for i in range(100):
	tmp = f.readline().split()
	tmp = map(float, tmp)
	A.append(tmp)

f.close()

#p = matrix_plot(A)

#import matplotlib.pyplot as plt
#p = plt.contour(A)

#from pprint import pprint
#pprint(p)

#b = []
#for i in range(100):
#	b.append((i/100).N())
#
#g = sip.interp2d(b,b,A,kind='linear')
#
#def fg(x,y):
#    return g(x,y).max()

b = np.linspace(0,1,100)
A = np.array(A)

g = sip.RectBivariateSpline(b,b,A)

def fg(x,y):
    return g(x,y).tolist()[0][0]

p = contour_plot(fg, (0,1), (0,1),
        contours=20, fill=False, cmap='hsv')
p.show()
