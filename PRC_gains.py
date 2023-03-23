import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Set Seaborn theme
matplotlib.style.use('seaborn')

# Open the file for reading
with open("PRC_gains_data.txt", "r") as f:
    # Read in the data from the file
    data = []
    for line in f:
        row = [float(val) for val in line.split()]
        data.append(row)
data = np.array(data).T.tolist()

print(data)
# Extract the data columns for plotting
x = [row[0] for row in data]
y1 = [row[1] for row in data]
y2 = [row[2] for row in data]
y3 = [row[4] for row in data]
y4 = [row[5] for row in data]

# Create the first plot with second and third columns as a function of the first one
plt.figure()
plt.plot(data[0], data[1], label="LSB (6MHz)")
plt.plot(data[0], data[2], label="USB (6MHz)")
plt.xlabel("PRM RoC (m)")
plt.ylabel("PRC optical gain")
plt.title("PRM RoC detuning impact on the 6MHz SB gain")
plt.legend()

# Create the second plot with fifth and sixth columns as a function of the first one
plt.figure()
plt.plot(data[0], data[3], label="LSB (56MHz)")
plt.plot(data[0], data[4], label="USB (56MHz)")
plt.xlabel("PRM RoC (m)")
plt.ylabel("PRC optical gain")
plt.title("PRM RoC detuning impact on the 56MHz SB gain")
plt.legend()
plt.show()
