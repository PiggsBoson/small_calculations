#Finding the eigenvectors of the system part of Hamiltonian
import numpy as np
from numpy import linalg as LA

sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]])
sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]])
sigma_y = 1.0j * np.array([[0.0, - 1.0], [1.0, 0.0]])
rt17 = np.sqrt(17)
rrt17p1 = np.sqrt(rt17 + 1)
rrt17m1 = np.sqrt(rt17 - 1)
coupling = np.kron(sigma_x,sigma_x) + np.kron(sigma_y,sigma_y) + np.kron(sigma_z,sigma_z) #isotropic coupling
dipole = np.kron(sigma_x,sigma_x) + np.kron(sigma_y,sigma_y)
# print(coupling)

H_As = -1/2*sigma_z+2*sigma_x #control with positive value
v0 = np.array([1/4*(rt17+1), -1]) #eigenvector 1
v0n = np.sqrt((17+rt17)/8)
v1 = np.array([1/4*(rt17-1), 1]) #eigenvector 2
v1n = np.sqrt((17-rt17)/8)
psi_a0 = v0/v0n #normalized eigenvector 1
psi_a1 = v1/v1n #normalized eigenvector 2
# print(psi_a0@H_As@psi_a0/(rt17/2) )
# print(H_As@psi_a0/(rt17/2) + psi_a0)
# print(LA.norm(psi_a0))
# print(psi_a1@H_As@psi_a1 )
# print(H_As@psi_a1/(rt17/2) - psi_a1)
# print(LA.norm(psi_a1))

U_A = np.sqrt(1/(2*rt17))*np.array([[rrt17p1, rrt17m1], [-4.0/rrt17p1, 4.0/rrt17m1]])
# U_A2 = np.kron(U_A,U_A)
print(U_A.conj().T @H_As@U_A/(-rt17/2) )
print(U_A.conj().T @psi_a0 )
print(U_A.conj().T @psi_a1 )
# resultA = U_A2.T@coupling@U_A2 - coupling
# filteredA = [[x if np.absolute(x) > 1e-6 else 0 for x in v] for v in resultA]
# # print(filteredA)

# H_Bs = -1/2*sigma_z-2*sigma_x #control with negative value
# U_B = np.sqrt(1/(2*rt17))*np.array([[rrt17p1, -rrt17m1], [4.0/rrt17p1, 4.0/rrt17m1]])
# U_B2 = np.kron(U_B,U_B)
# # print(U_B.conj().T @H_Bs@U_B/(rt17/2) )
# resultB = U_B2.T@coupling@U_B2 - coupling
# filteredB = [[x if np.absolute(x) > 1e-6 else 0 for x in v] for v in resultB]
# # print(filteredB)

# resultAdp = U_A2.T@dipole@U_A2 - dipole
# filteredAdp = [[x if np.absolute(x) > 1e-6 else 0 for x in v] for v in resultAdp]
# # print(filteredAdp)