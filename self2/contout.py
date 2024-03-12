import numpy as np
import matplotlib.pyplot as plt

# out of bounds counter

out = 0

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
list_x_1 = []
list_y_1 = []
list_x_2 = []
list_y_2 = []
list_x_3 = []
list_y_3 = []
list_x_4 = []
list_y_4 = []
list_x_5 = []
list_y_5 = []

''' with open("PS_one_particle1.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        list_x_1.append(float(x))
        list_y_1.append(float(y))
        
with open("PS_one_particle2.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        list_x_2.append(float(x))
        list_y_2.append(float(y))
        
with open("PS_one_particle3.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        list_x_3.append(float(x))
        list_y_3.append(float(y))

with open("PS_one_particle4.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        list_x_4.append(float(x))
        list_y_4.append(float(y))
    
with open("PS_one_particle5.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        list_x_5.append(float(x))
        list_y_5.append(float(y))
 '''

''' with open("PS_group_best.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        list_x_1.append(float(x))
        list_y_1.append(float(y)) '''
        
'''with open("PS_output.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        x = float(x)
        y = float(y)
        if(x > 0.25 or x < -0.25 or y > 0.25 or y < -0.25):
            out = out + 1 
        list_x_1.append(x)
        list_y_1.append(y)'''
        
with open("SA_one_best.txt", 'r') as file:
    # Iterate over each line in the file
    for line in file:
        # Split the line into components based on whitespace or another delimiter
        x, y = line.split()
        list_x_1.append(float(x))
        list_y_1.append(float(y))

# Create the contour plot
plt.figure(figsize=(8, 6))
contour = plt.contour(X, Y, Z, levels=np.logspace(-1, 3, 50), cmap='viridis')
plt.xlabel('x')
plt.ylabel('y')

# Plot the points from the lists
plt.scatter(list_x_1, list_y_1, color='black')
''' plt.scatter(list_x_2, list_y_2, color='blue')
plt.scatter(list_x_3, list_y_3, color='green')
plt.scatter(list_x_4, list_y_4, color='red')
plt.scatter(list_x_5, list_y_5, color='yellow') '''

# Connect the points with lines
for i in range(len(list_x_1)-1):
   plt.plot([list_x_1[i], list_x_1[i+1]], [list_y_1[i], list_y_1[i+1]], 'k-')
   '''plt.plot([list_x_2[i], list_x_2[i+1]], [list_y_2[i], list_y_2[i+1]], 'b-')
   plt.plot([list_x_3[i], list_x_3[i+1]], [list_y_3[i], list_y_3[i+1]], 'g-')
   plt.plot([list_x_4[i], list_x_4[i+1]], [list_y_4[i], list_y_4[i+1]], 'r-')
   plt.plot([list_x_5[i], list_x_5[i+1]], [list_y_5[i], list_y_5[i+1]], 'y-') '''

print(out)
plt.show()
