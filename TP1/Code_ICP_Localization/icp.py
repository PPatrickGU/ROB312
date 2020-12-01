"""
 Simple 2D ICP implementation
 author: David Filliat
"""

import numpy as np
from scipy.spatial import KDTree
import math
import heapq


# A few helper function

def angle_wrap(a):
    """
    Keep angle between -pi and pi
    """
    return np.fmod(a + np.pi, 2*np.pi ) - np.pi


def mean_angle(angleList):
    """
    Compute the mean of a list of angles
    """

    mcos = np.mean(np.cos(angleList))
    msin = np.mean(np.sin(angleList))

    return math.atan2(msin, mcos)

def get_index1(lst=None, item=''):
    return [index for (index,value) in enumerate(lst) if value == item]

def icp(model, data, maxIter, thres):
    """
    ICP (iterative closest point) algorithm
    Simple ICP implementation for teaching purpose
    - input
    model : scan taken as the reference position
    data : scan to align on the model
    maxIter : maximum number of ICP iterations
    thres : threshold to stop ICP when correction is smaller
    - output
    R : rotation matrix
    t : translation vector
    meandist : mean point distance after convergence
    """

    print('Running ICP, ', end='')

    # Various inits
    olddist = float("inf")  # residual error
    maxRange = 10  # limit on the distance of points used for ICP

    # Create array of x and y coordinates of valid readings for reference scan
    valid = model["ranges"] < maxRange
    ref = np.array([model["x"], model["y"]])
    ref = ref[:, valid] #size = [2, 450+]

    # Create array of x and y coordinates of valid readings for processed scan
    valid = data["ranges"] < maxRange
    dat = np.array([data["x"], data["y"]])
    dat = dat[:, valid]

    # ----------------------- TODO ------------------------
    # Filter data points too close to each other
    # Put the result in dat_filt

    delete_list = [] #(2,461)
    # dat_filt = dat
    for i in range(dat.shape[1]-1):
        if np.linalg.norm((dat.T)[i+1] - (dat.T)[i]) < 0.02:
            delete_list.append(i)
    dat_filt = np.delete(dat, delete_list, 1)
    print(dat_filt.shape)

    # Initialize transformation to identity
    R = np.eye(2)
    t = np.zeros((2, 1))

    # Main ICP loop
    for iter in range(maxIter):

        # ----- Find nearest Neighbors for each point, using kd-trees for speed
        tree = KDTree(ref.T)
        distance, index = tree.query(dat_filt.T)
        meandist = np.mean(distance)
        #
        # # ----------------------- TODO ------------------------
        # # filter points matchings, keeping only the closest ones
        # # you have to modify :
        # # - 'dat_matched' with the points
        # # - 'index' with the corresponding point index in ref

        dat_matched = dat_filt.copy()

        XX = 0.9
        # #ata_index = np.sort(np.argsort(distance)[:int(len(distance) * XX)])
        data_index = np.argsort(distance)[:int(len(distance) * XX)]
        index = index[data_index]
        dat_matched = dat_filt.copy()[:, data_index]

        # # normal shooting matchings, keeping only the best ones
        # index = np.zeros(dat_filt.shape[1], dtype = int)
        # cos = np.zeros(dat_filt.shape[1])
        # for i in range(dat_filt.shape[1] - 1):
        #     vector_tan = dat_filt.T[i + 1] - dat_filt.T[i]
        #     for j in range(ref.shape[1]):
        #         vector_nor = dat_filt.T[i] - ref.T[j]
        #         min = np.inf
        #         cos_tmp = abs(np.dot(vector_tan, vector_nor.T) / (
        #                 np.linalg.norm(vector_tan) * np.linalg.norm(vector_nor)))
        #         if cos_tmp < min:
        #             index[i] = j
        #             cos [i] = cos_tmp
        # XX = 0.6
        # data_index = np.sort(np.argsort(cos)[:int(len(cos) * XX)])
        # index = index[data_index]
        # dat_matched = dat_filt.copy()[:, data_index]



        # ----- Compute transform

        # Compute point mean
        mdat = np.mean(dat_matched, 1)
        mref = np.mean(ref[:, index], 1)

        # Use SVD for transform computation
        C = np.transpose(dat_matched.T-mdat) @ (ref[:, index].T - mref)
        u, s, vh = np.linalg.svd(C)
        Ri = vh.T @ u.T
        Ti = mref - Ri @ mdat

        # Apply transformation to points
        dat_filt = Ri @ dat_filt
        dat_filt = np.transpose(dat_filt.T + Ti)

        # Update global transformation
        R = Ri @ R
        t = Ri @ t + Ti.reshape(2, 1)

        # Stop when no more progress
        if abs(olddist-meandist) < thres:
            break

        # store mean residual error to check progress
        olddist = meandist

    print("finished with mean point corresp. error {:f}".format(meandist))

    return R, t, meandist
