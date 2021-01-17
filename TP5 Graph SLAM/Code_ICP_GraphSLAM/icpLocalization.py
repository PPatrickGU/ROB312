"""
 Simple ICP localisation demo
 Compute position of each scan using ICP
 with respect to the previous one
 author: David Filliat
"""

import readDatasets as datasets
import matplotlib.pyplot as plt
import icp
import numpy as np
import copy

# Reading data
#scanList = datasets.read_fr079(0)
scanList = datasets.read_u2is(0)

# Copy for reference display
odomScanList = copy.deepcopy(scanList)

# Parameters for scan processing
minScan = 1
step = 5
maxScan = len(scanList)-step

# Init displays
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(14, 7))

c = np.random.rand(3,)
ax1.scatter(odomScanList[minScan]["x"], odomScanList[minScan]["y"], color=c, s=1)
ax1.scatter(odomScanList[minScan]["pose"][0], odomScanList[minScan]["pose"][1], color=c, s=3)
ax1.axis([-5.5, 12.5, -12.5, 6.5])
ax1.set_title('Pose from raw odometry')
ax2.scatter(scanList[minScan]["x"], scanList[minScan]["y"], color=c, s=1)
ax2.scatter(scanList[minScan]["pose"][0], scanList[minScan]["pose"][1], color=c, s=3)
ax2.axis([-5.5, 12.5, -12.5, 6.5])
ax2.set_title('Pose after ICP correction')
plt.pause(0.1)


for a in range(minScan, maxScan, step):
    s1 = scanList[a]
    s2 = scanList[a+step]
    # perform ICP
    R, t, error = icp.icp(s1, s2, 200, 1e-7)

    # correct future scans
    for b in range((a+step), maxScan, step):
        scanList[b] = datasets.transform_scan(scanList[b], R, t)

    # Display
    c = np.random.rand(3,)
    ax1.scatter(odomScanList[a+step]["x"], odomScanList[a+step]["y"], color=c, s=1)
    ax1.scatter(odomScanList[a+step]["pose"][0], odomScanList[a+step]["pose"][1], color=c, s=3)
    ax2.scatter(scanList[a+step]["x"], scanList[a+step]["y"], color=c, s=1)
    ax2.scatter(scanList[a+step]["pose"][0], scanList[a+step]["pose"][1], color=c, s=3)
    plt.pause(0.1)

plt.savefig('ICPLocalization.png')
print("Press Q in figure to finish...")
plt.show()
