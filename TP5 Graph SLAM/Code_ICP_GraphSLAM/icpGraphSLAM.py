"""
 Graph SLAM with ICP - Basic implementation for teaching purpose only...
 Computes position of each scan with respect to several scans in the
 current map and optimizes the positions to reduce errors, the add scan
 to the map if it is far enough of all the existing ones
 author: David Filliat
"""

import readDatasets as datasets
import matplotlib.pyplot as plt
import icp
import numpy as np
import copy


# Reading data
# scanList = datasets.read_fr079(0)
scanList = datasets.read_u2is(0)

# Parameters for scan processing
minScan = 0
step = 5
maxScan = len(scanList) - step

# Parameters for map building
distThresholdAdd = 0.25
distThresholdMatch = 1
maxICPError = 0.4

# Copy for reference display and map init
odomScanList = copy.deepcopy(scanList)
map = [scanList[minScan]]

# Initialize graph of relative positions
# simple representation using weight matrixes
maxsize = int(np.around((maxScan - minScan)/step))
Graphtx = np.zeros((maxsize, maxsize))
Graphty = np.zeros((maxsize, maxsize))
Graphtheta = np.zeros((maxsize, maxsize))

# Init displays
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(14, 7))
c = np.random.rand(3,)
ax1.scatter(odomScanList[minScan]["x"],
            odomScanList[minScan]["y"], color=c, s=1)
ax1.scatter(odomScanList[minScan]["pose"][0],
            odomScanList[minScan]["pose"][1], color=c, s=3)
ax1.axis([-5.5, 12.5, -12.5, 6.5])
ax1.set_title('Poses from raw odometry')
ax2.scatter(map[0]["x"], map[0]["y"], color=c, s=1)
ax2.scatter(map[0]["pose"][0],
            map[0]["pose"][1], color=c, s=3)
ax2.axis([-5.5, 12.5, -12.5, 6.5])
ax2.set_title('Map after ICP Graph SLAM')
plt.pause(0.1)


# Process scans
for i in range(minScan + step, maxScan, step):

    print('Processing new scan')

    # get list of map scan sorted by distance
    sorteddist, sortedId = datasets.find_closest_scan(map, scanList[i])

    # Keep only the ones below the distance threshold, or the closest one
    closeScans = sortedId[sorteddist < distThresholdMatch]
    if len(closeScans) == 0:
        closeScans = [sortedId[0]]

    # perform ICP with closest scan to correct future odometry
    R, t, error = icp.icp(map[closeScans[0]], scanList[i], 200, 1e-7)

    # Correct all future scans odometry pose
    for j in range(i, maxScan, step):
        scanList[j] = datasets.transform_scan(scanList[j], R, t)

    # --- Add scan to map and update graph if needed
    if np.linalg.norm(scanList[i]["pose"][0:2] -
                      map[closeScans[0]]["pose"][0:2]) > distThresholdAdd:

        map.append(scanList[i])
        print('Adding new scan with links to : ' + str(closeScans))

        # Get ref to last scan in map (i.e. new scan)
        id2 = len(map) - 1
        s2 = map[-1]

        # ---- Build graph
        edgeNB=0
        for idi in closeScans:
            print("Adding edge")
            # take the reference scan among the closest map scan
            s1 = map[idi]

            # compute position of new scan wrt the ref scan
            Ri, ti, error = icp.icp(s1, s2, 200, 1e-7)

            if error < maxICPError or edgeNB == 0:
                edgeNB += 1
                # compute absolute position of new scan
                si = datasets.transform_scan(s2, Ri, ti)

                # compute relative pose with ref scan
                deltatheta = icp.angle_wrap(si["pose"][2] - s1["pose"][2])
                deltat = si["pose"][0:2] - s1["pose"][0:2]

                # Add relative position in the graph
                Graphtx[idi, id2] = deltat[0]
                Graphty[idi, id2] = deltat[1]
                Graphtheta[idi, id2] = deltatheta

                # Make graph symetric
                Graphtx[id2, idi] = - Graphtx[idi, id2]
                Graphty[id2, idi] = - Graphty[idi, id2]
                Graphtheta[id2, idi] = - Graphtheta[idi, id2]

        # --- Optimize graph until updates fall below threshold
        updateMax = 1
        updateNB = 0
        while updateMax > 1e-5:
            updateMax = 0
            updateNB += 1

            # Recompute each scan pose from its neighbors
            for k in range(len(map)):
                # Create a list of scan pose computed through neighbor pose and
                # relative position
                newPoseList = []
                for l in range(len(map)):
                    if Graphtx[l, k] != 0:
                        newAngle = icp.angle_wrap(map[l]["pose"][2] + Graphtheta[l, k])
                        newX = map[l]["pose"][0] + Graphtx[l, k]
                        newY = map[l]["pose"][1] + Graphty[l, k]
                        newPoseList.append(newX)
                        newPoseList.append(newY)
                        newPoseList.append(newAngle)

                # Compute new pose as mean of positions from neighbors and update map
                newPose = np.array([np.mean(newPoseList[0::3]),
                                    np.mean(newPoseList[1::3]),
                                    icp.mean_angle(newPoseList[2::3])])
                update = abs(newPose.reshape(-1) - map[k]["pose"].reshape(-1))
                map[k] = datasets.update_scan_pose(map[k], newPose)
                updateMax = max(max(update), updateMax)

        print('New map size : ' + str(len(map)) + ', Graph Updates : ' + str(updateNB))

    # Display
    c = np.random.rand(3,)
    ax1.scatter(odomScanList[i]["x"],
                odomScanList[i]["y"], color=c, s=1)
    ax1.scatter(odomScanList[i]["pose"][0],
                odomScanList[i]["pose"][1], color=c, s=3)
    ax2.cla()
    for scan in map:
        c = np.random.rand(3,)
        ax2.scatter(scan["x"], scan["y"], color=c, s=1)
        ax2.scatter(scan["pose"][0],
                    scan["pose"][1], color=c, s=3)
    ax2.axis([-5.5, 12.5, -12.5, 6.5])
    ax2.set_title('Map after ICP Graph SLAM')
    plt.pause(0.1)

plt.savefig('ICPGraphSLAM.png')
print("Press Q in figure to finish...")
plt.show()
