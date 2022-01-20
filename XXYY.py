'''
Check if xx+yy=2(pm+mp) (yes!)
'''
import numpy as np

sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]])
sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]])
sigma_y = 1.0j * np.array([[0.0, -1.0], [1.0, 0.0]])
sigma_p = 0.5*(sigma_x+1.0j*sigma_y)
sigma_m = 0.5*(sigma_x-1.0j*sigma_y)

print(sigma_p)
print(sigma_m)

xx = np.kron(sigma_x, sigma_x)
yy = np.kron(sigma_y, sigma_y)

pm = np.kron(sigma_p, sigma_m)
mp = np.kron(sigma_m, sigma_p)

print(xx+yy-2*(pm+mp))