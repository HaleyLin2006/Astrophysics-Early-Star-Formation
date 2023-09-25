import numpy as np
from scipy.special import legendre
from numpy.polynomial import Polynomial

# Define the degree of the Legendre polynomial
degree = 3

# Compute the Legendre polynomial of the specified degree
legendre_poly = legendre(degree)

# Convert the Legendre polynomial coefficients to a NumPy Polynomial
legendre_poly_np = Polynomial(legendre_poly)

# Calculate the derivative of the Legendre polynomial
derivative_poly_np = legendre_poly_np.deriv()(np.pi/2)

# Print the original polynomial and its derivative
print(f"Legendre Polynomial P_{degree}(x): {legendre_poly_np}")
print(f"Derivative of P_{degree}(x): {derivative_poly_np}")
