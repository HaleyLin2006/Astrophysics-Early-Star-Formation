import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

# Generate sample data in polar coordinates
theta = np.linspace(0, 2 * np.pi, 100)  # Angle values
r = np.linspace(0, 100, 100)  # Random radial distances
data = np.random.random(100)  # Random data values

# Convert polar coordinates to Cartesian coordinates
x = r * np.cos(theta)
y = r * np.sin(theta)

# Define a grid in (x, y) coordinates
x_grid = np.linspace(min(x), max(x), 300)
y_grid = np.linspace(min(y), max(y), 300)
x_grid, y_grid = np.meshgrid(x_grid, y_grid)

# Interpolate the data onto the grid
data_interp = griddata((x, y), data, (x_grid, y_grid), method='linear')

# Plot the interpolated data using imshow()
plt.imshow(data_interp, extent=(0, max(x), 0, max(y)),
           cmap='viridis', origin='lower', aspect='auto')
plt.colorbar(label='Data')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Interpolated Data in Cartesian Coordinates')
plt.show()
