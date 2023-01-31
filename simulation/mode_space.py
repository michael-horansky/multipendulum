
from class_multipendulum import *
import matplotlib.cm as cm


def scalar_field(x, y):
    return(x ** 2 + y ** 2)






n = 256
x = np.linspace(-3., 3., n)
y = np.linspace(-3., 3., n)

xm, ym = np.meshgrid(x, y)
zm = scalar_field(xm, ym)

plt.pcolormesh(xm, ym, zm)
plt.colorbar()
plt.show()
