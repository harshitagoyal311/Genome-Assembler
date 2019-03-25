import numpy as np
from matplotlib import pyplot as plt

# file1 = open("test0", "r")
# file2 = open("test1", "r")
# file3 = open("test2", "r")
# file4 = open("test3", "r")

# n1 = file1.readline()
# n2 = file2.readline()
# n3 = file3.readline()
# n4 = file4.readline()
# d1 = file1.readline()
# d2 = file2.readline()
# d3 = file3.readline()
# d4 = file4.readline()

data1 = np.loadtxt("out0", skiprows=2)
data2 = np.loadtxt("out1", skiprows=2)
data3 = np.loadtxt("out2", skiprows=2)
data4 = np.loadtxt("out3", skiprows=2)
data5 = np.loadtxt("out4", skiprows=2)
data6 = np.loadtxt("out5", skiprows=2)
data7 = np.loadtxt("out6", skiprows=2)
data8 = np.loadtxt("out7", skiprows=2)

x1 = data1.T[1]
y1 = data1.T[2]

x2 = data2.T[1]
y2 = data2.T[2]

x3 = data3.T[1]
y3 = data3.T[2]

x4 = data4.T[1]
y4 = data4.T[2]

x5 = data5.T[1]
y5 = data5.T[2]

x6 = data6.T[1]
y6 = data6.T[2]

x7 = data7.T[1]
y7 = data7.T[2]

x8 = data8.T[1]
y8 = data8.T[2]


plt.scatter(x1, y1, color = "b", s=10)
plt.scatter(x2, y2, color = "g", s=10)
plt.scatter(x3, y3, color = "r", s=10)
plt.scatter(x4, y4, color = "c", s=10)
plt.scatter(x5, y5, color = "m", s=10)
plt.scatter(x6, y6, color = "y", s=10)
plt.scatter(x7, y7, color = "k", s=10)
plt.scatter(x8, y8, color = "0.75", s=10)

plt.show()
print("Plotting done")