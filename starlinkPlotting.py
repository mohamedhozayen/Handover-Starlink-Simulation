# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 03:24:22 2022

@author: Mohamed Hozayen


"""

#import sys
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

# Read CSV
#csvFileName = sys.argv[1]
df = pd.read_csv('Starlink_Sat_Pos.csv')
df = pd.DataFrame(df, columns=['Sat','X','Y','Z'])

X, Y, Z = df['X'], df['Y'], df['Z']


st=0
end=len(X)

# Plot X,Y,Z
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X[st:end], Y[st:end], Z[st:end], color='white', edgecolors='grey', alpha=0.5)
ax.scatter(X[st:end], Y[st:end], Z[st:end], c='red')
plt.show()


n1=0
n2=1
earth_radius=6378
for i in range(0, len(X)-1):
    dxs=(X[n1]-X[n2])**2
    dys=(Y[n1]-Y[n2])**2
    hs=(Z[n1]-Z[n2])**2
    print((hs+dxs+dys)**(1/2))
    n1+=1
    n2+=1

# 6928.139 distance from origin
# 556.139 distance from earth
# 604.4024 distance between two adjacent satellite

f=1/2
x1=X[3]*f
y1=Y[3]*f
z1=Z[3]*f
dxs=(x1-0)**2
dys=(y1-0)**2
hs=(z1-0)**2
print((hs+dxs+dys)**(1/2)-earth_radius)

h=550
dxs=(6341.906185-(6341.906185))**2
dys=(-2393.439206-(-2393.439206))**2
hs=(1.432061e+03-(1.432061e+03-h))**2
print((hs+dxs+dys)**(1/2))

v=earth_radius*.707107
dxs=(v-0)**2
dys=(v-0)**2
hs=(0-0)**2
print((hs+dxs+dys)**(1/2)-earth_radius)


x=1337.4490938501156
y= 4161.6663053555785
z=4644.4211140063035
dxs=(x-0)**2
dys=(y-0)**2
hs=(z-0)**2
print((hs+dxs+dys)**(1/2))