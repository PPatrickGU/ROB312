"""
Reading a laser scan dataset
author: David Filliat
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import icp


def create_scan(ranges, angles, pose):
    """
    create_scan creates a scan dict from its components
    - input
    ranges : np.array of scan distances
    angles : np.array of angles of each reading
    pose : absolute pose as np.array [x,y,theta]
    x,y: absolute position of each scan point, convenient for icp and plot
    - output
    scan : a dict with ranges, angles, pose, x & y coordinates of scan points
    """

    scan = {'ranges': np.array(ranges),
            'angles': np.array(angles),
            'pose': pose.reshape(-1),
            'x': pose[0] + np.multiply(ranges, np.cos(angles + pose[2])),
            'y': pose[1] + np.multiply(ranges, np.sin(angles + pose[2])),
            }

    return scan


def transform_scan(scan, R, t):
    """
    Change the pose of a scan with rotation matrix and translation vector
    input
    - scan : scan structure from 'create_scan' function
    - R : 2x2 rotation matrix
    - t : 2x1 translation vector
    """
    pose = scan["pose"]
    newXY = np.matmul(R, pose[0:2].reshape(2,-1)) + t
    newTheta = pose[2] + math.atan2(R[1, 0], R[0, 0])
    
    newscan = dict()
    newscan["pose"] = np.array([newXY[0,0], newXY[1,0], newTheta]).reshape(-1)
    newscan["ranges"] = scan["ranges"]
    newscan["angles"] = scan["angles"]

    newscan["x"] = newscan["pose"][0] + \
        np.multiply(newscan["ranges"], np.cos(newscan["angles"] +
                                              newscan["pose"][2]))
    newscan["y"] = newscan["pose"][1] + \
        np.multiply(newscan["ranges"], np.sin(newscan["angles"] +
                                              newscan["pose"][2]))

    return newscan


def update_scan_pose(scan, newPose):
    """
    Update the pose of a scan to a new pose
    - scan : scan structure from 'create_scan' function
    - newPose : new pose ...
    """

    newscan = dict()
    newscan["pose"] = newPose.reshape(-1)
    newscan["ranges"] = scan["ranges"]
    newscan["angles"] = scan["angles"]

    newscan["x"] = newscan["pose"][0] + \
        np.multiply(newscan["ranges"], np.cos(newscan["angles"] +
                                              newscan["pose"][2]))
    newscan["y"] = newscan["pose"][1] + \
        np.multiply(newscan["ranges"], np.sin(newscan["angles"] +
                                              newscan["pose"][2]))

    return newscan
    

def find_closest_scan(map, scan):
    """
    Return map scan ids sorted according to distance to scan
    """

    def distance(scan1, scan2):
        """
        Computes distance between two scan pose
        """
        dist = np.linalg.norm(scan1["pose"][0:2] - scan2["pose"][0:2]) + abs(icp.angle_wrap(scan1["pose"][2] - scan2["pose"][2]))/15
        # The following is to prevent matching scans with too large rotation
        # differences
        if abs(icp.angle_wrap(scan1["pose"][2]-scan2["pose"][2])) > np.pi/3:
            dist = dist * 2
            
        return dist

    distances = [distance(scan, previous_scan) for previous_scan in map]
    distances = np.array(distances).reshape(-1)
    sorted_id = np.argsort(distances)

    return distances[sorted_id], sorted_id


def read_u2is_odom_entry(file):
    """
    read the next odometry entry
    format specific to 'u2is' dataset
    """

    secs = file.readline()
    secs = int(secs[10:])
    nsecs = file.readline()
    nsecs = int(nsecs[11:])
    odomTime = [secs, nsecs]
    odomData = []
    for i in range(13):
        odomData.append(float(file.readline()))

    return odomTime, odomData


def read_u2is_laser_entry(file):
    """
    read the next  laser scan entry
    format specific to 'u2is' dataset
    """

    secs = file.readline()
    secs = int(secs[10:])
    nsecs = file.readline()
    nsecs = int(nsecs[11:])
    laserTime = [secs, nsecs]
    line = file.readline()
    laserData = line[9:-2].split(',')
    laserData = np.array([float(i) for i in laserData])

    # print("Reading laser : " + str(laserTime))

    return laserTime, laserData


def read_u2is(number):
    """
    Reading and formating u2is dataset
    - input : number of scans to read
    - output : list of dict with scans
    """

    if number == 0 or number > 855:
        number = 845

    print('Reading u2is dataset')

    fileLaser = open('dataset/u2is/laser_filt.txt', 'r')
    fileOdom = open('dataset/u2is/odom_filt.txt', 'r')

    scanList = []

    # reading first odom
    odomTime, odomData = read_u2is_odom_entry(fileOdom)

    angles = np.arange(-2.35619449615, 2.35619449615, 0.00436332309619)

    for i in range(number):

        # Reading raw data
        laserTime, laserData = read_u2is_laser_entry(fileLaser)

        # time synchronization
        while odomTime[0] > laserTime[0]:
            laserTime, laserData = read_u2is_laser_entry(fileLaser)

        if odomTime[0] == laserTime[0]:
            while odomTime[1] > laserTime[1]:
                laserTime, laserData = read_u2is_laser_entry(fileLaser)

        # Replace max range readings and remove borders
        laserData[laserData > 20.0] = np.inf
        laserData[0:80] = np.inf
        laserData[-80:] = np.inf

        # Converting quaternion to yaw angle
        siny_cosp = 2.0 * (odomData[6] * odomData[5] +
                           odomData[3] * odomData[4])
        cosy_cosp = 1.0 - 2.0 * (odomData[4] * odomData[4] +
                                 odomData[5] * odomData[5])
        yaw = math.atan2(siny_cosp, cosy_cosp)

        # Compute position of laser
        pose = np.array([odomData[0] + 0.1*np.cos(yaw),
                odomData[1] + 0.1*np.sin(yaw), yaw])

        # Create scanlist
        scanList.append(create_scan(laserData[0::2], angles[0::2], pose))

        # reading next odom
        odomTime, odomData = read_u2is_odom_entry(fileOdom)

    fileLaser.close()
    fileOdom.close()
    print('Finished reading ' + str(len(scanList)) + ' scans')

    return scanList


def read_fr079_odom_entry(file):
    """
    read the next odometry entry
    format specific to 'fr079' dataset
    """

    secs = file.readline()
    secs = int(secs[10:])
    nsecs = file.readline()
    nsecs = int(nsecs[11:])
    odomTime = [secs, nsecs]
    odomData = []
    for i in range(13):
        odomData.append(float(file.readline()))

    return odomTime, odomData


def read_fr079_laser_entry(file):
    """
    read the next  laser scan entry
    format specific to 'fr079' dataset
    """

    secs = file.readline()
    secs = int(secs[10:])
    nsecs = file.readline()
    nsecs = int(nsecs[11:])
    laserTime = [secs, nsecs]
    line = file.readline()
    laserData = line[9:-2].split(',')
    laserData = np.array([float(i) for i in laserData])

    # print("Reading laser : " + str(laserTime))

    return laserTime, laserData


def read_fr079(number):
    """
    Reading and formating u2is dataset
    - input : number of scans to read
    - output : list of dict with scans
    """

    if number == 0 or number > 4919:
        number = 4919

    print('Reading FR079 dataset')

    file = open('dataset/fr079/laserData.txt', 'r')

    scanList = []

    # discard first line
    line = file.readline()

    angles = np.arange(-math.pi / 2, math.pi / 2, math.pi / 359)

    for i in range(number):

        # Reading raw data
        line = file.readline()
        rawData = line[11:-20].split(' ')
        rawData = np.array([float(i) for i in rawData])

        # Create scanlist
        scanList.append(create_scan(rawData[0:360], angles, rawData[360:363]))

    file.close()
    print('Finished reading ' + str(len(scanList)) + ' scans')

    return scanList

if __name__ == '__main__':
    read_fr079(10)
