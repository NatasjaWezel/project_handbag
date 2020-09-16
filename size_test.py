import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import matplotlib

sizes = [6.863417982155114, 1.0295126973232671, 8.579272477693893, 2.40219629375429, 1.0295126973232671]
x = [1, 3, 4, 2, 1]
y = [2, 3, 4, 1, 2]
z = [2, 1, 4, 2, 1]

bigger_sizes = [100 * i for i in sizes]
print(bigger_sizes)

norm = plt.Normalize(min(sizes), max(sizes))
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","blue","red"])


##########################################################################################################
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(5):
    p = ax.scatter(x[i], y[i], z[i], s=bigger_sizes[i], c=sizes[i], cmap=cmap, norm=norm)

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

ax.set_xlim(0, 5)
ax.set_ylim(0, 5)
ax.set_zlim(0, 5)

fig.colorbar(p)
plt.show()


###########################################################################################################
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

p = ax.scatter(x, y, z, s=bigger_sizes, c=sizes, cmap=cmap, norm=norm)

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

ax.set_xlim(0, 5)
ax.set_ylim(0, 5)
ax.set_zlim(0, 5)

fig.colorbar(p)
plt.show()