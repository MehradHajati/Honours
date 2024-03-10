import numpy as np
import matplotlib.pyplot as plt

# Define the Three-Hump Camel function
def three_hump_camel(x, y):
    return 2*x**2 - 1.05*x**4 + (x**6) / 6 + x*y + y**2

# Generate grid of x and y values for the contour
x = np.linspace(-5, 5, 400)
y = np.linspace(-5, 5, 400)
X, Y = np.meshgrid(x, y)

# Compute Z values for the contour
Z = three_hump_camel(X, Y)

# Sample lists of x and y values
list_x = []
list_y = []

with open("PS_many.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        list_x.append(float(x))
        list_y.append(float(y))



# Create the contour plot
plt.figure(figsize=(8, 6))
contour = plt.contour(X, Y, Z, levels=np.logspace(-1, 3, 50), cmap='viridis')
plt.xlabel('x')
plt.ylabel('y')

# Plot the points from the lists
plt.scatter(list_x, list_y, color='black')

# Connect the points with lines
#for i in range(len(list_x)-1):
#    plt.plot([list_x[i], list_x[i+1]], [list_y[i], list_y[i+1]], 'k-')

plt.show()