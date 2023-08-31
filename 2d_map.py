import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

# Example data in (r, theta) format
theta_values = np.linspace(0, 2 * np.pi, 100)
r_values = np.linspace(0, 5, 100)
data_values = []  # List of lists, each sublist represents data at a theta value

# Fill data_values with example values (you should replace this with your actual data)
for theta in range(len(theta_values)):
    data = 100*np.sin(3 * theta_values[theta]) / r_values
    print(data)  # Example data at a specific theta
    data_values.append(data)

# Convert (r, theta) to (x, y) coordinates
x_values = []
y_values = []
for i, theta in enumerate(theta_values):
    x_values.extend(r_values * np.cos(theta))
    y_values.extend(r_values * np.sin(theta))

# Create a flattened (x, y) array and corresponding flattened data array
x_flat = np.array(x_values)
y_flat = np.array(y_values)
data_flat = np.array(data_values).flatten()

# Create a grid in (x, y) coordinates
x_grid = np.linspace(min(x_flat), max(x_flat), 300)
y_grid = np.linspace(min(y_flat), max(y_flat), 300)
xx, yy = np.meshgrid(x_grid, y_grid)

# Interpolate data onto the (x, y) grid
data_interp = griddata((x_flat, y_flat), data_flat, (xx, yy), method='linear')

# Plot the interpolated data using imshow()
plt.imshow(data_interp, extent=[min(x_flat), max(x_flat), min(y_flat), max(y_flat)],
           cmap='viridis', origin='lower', aspect='auto')

# Add labels, title, and color bar
plt.colorbar(label='Value')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Data Dependent on (r, theta) Coordinates in (x, y) Plane')

plt.show()
