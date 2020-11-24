import numpy as np
import matplotlib.pyplot as plt

xedges = [0, 1, 3, 5]
yedges = [0, 2, 3, 4, 6, 9]

fig, ax = plt.subplots()
x = np.random.normal(2, 1, 100)
y = np.random.normal(1, 1, 100)
H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))
H = H.T  # Let each row list bins with common y range.
plt.imshow(H, interpolation='nearest', origin='low',
         extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.savefig("test.png")

print H

#print xedges

print (xedges[1:] + xedges[:-1]) / 2

resp = []
rms = []

for i in range(len(H[0])):
    #print(H[:,i])
    #print(np.mean(H[:,i]))
    #print(np.sqrt(np.mean(H[:,i]**2)))
    resp.append(np.mean(H[:,i]))
    rms.append(np.sqrt(np.mean(H[:,i]**2)))
