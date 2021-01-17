"""
 Incremental ICP SLAM  - Basic implementation for teaching purpose only...
 Computes position of each scan with respect to closest one in the
 current map and add scan to the map if it is far enough from all the
 existing ones
 author: David Filliat
"""

import readDatasets as datasets
import matplotlib.pyplot as plt
import icp
import numpy as np
import copy


# Reading data
scanList = datasets.read_fr079(0)
# scanList = datasets.read_u2is(0)

# Parameters for scan processing
minScan = 0
step = 3
maxScan = len(scanList) - step

# Parameters for map building
distThresholdAdd = 0.01

# Copy for reference display and map init
odomScanList = copy.deepcopy(scanList)
map = [scanList[minScan]]

# Init displays
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(14, 7))
c = np.random.rand(3,)
ax1.scatter(odomScanList[minScan]["x"],
            odomScanList[minScan]["y"], color=c, s=1)
ax1.axis([-5.5, 12.5, -12.5, 6.5])
ax1.set_title('Pose from raw odometry')
ax2.scatter(map[0]["x"], map[0]["y"], color=c, s=1)
ax2.axis([-5.5, 12.5, -12.5, 6.5])
ax2.set_title('Map after incremental ICP SLAM')
plt.pause(0.1)

# Perform incremental SLAM
for i in range(minScan + step, maxScan, step):

    # get list of map scan sorted by distance
    sorteddist, sortedId = datasets.find_closest_scan(map, scanList[i])
    refScanId = sortedId[0]
    print('Matching new scan to reference scan ' + str(refScanId))

    # perform ICP with closest scan
    R, t, error = icp.icp(map[refScanId], scanList[i], 200, 1e-7)

    # Correct all future scans odometry pose
    for j in range(i, maxScan, step):
        scanList[j] = datasets.transform_scan(scanList[j], R, t)

    # Add scan to map if it is far enough
    if np.linalg.norm(scanList[i]["pose"][0:2] -
                      map[refScanId]["pose"][0:2]) > distThresholdAdd:
        map.append(scanList[i])
        print('Added to map, new size : ' + str(len(map)))

        # Map display
        ax2.cla()
        for scan in map:
            c = np.random.rand(3,)
            ax2.scatter(scan["x"], scan["y"], color=c, s=1)
            ax2.scatter(scan["pose"][0],
                        scan["pose"][1], color=c, s=3)
        ax2.axis([-5.5, 12.5, -12.5, 6.5])
        ax2.set_title('Map after incremental ICP SLAM')

    # Scan isplay
    c = np.random.rand(3,)
    ax1.scatter(odomScanList[i]["x"],
                odomScanList[i]["y"], color=c, s=1)
    ax1.scatter(odomScanList[i]["pose"][0],
                odomScanList[i]["pose"][1], color=c, s=3)
    plt.pause(0.1)

plt.savefig('ICPIncrementalSLAM.png')
print("Press Q in figure to finish...")
plt.show()
