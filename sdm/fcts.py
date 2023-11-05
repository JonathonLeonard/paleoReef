import os
import warnings
warnings.filterwarnings("ignore")

import struct
import subprocess
import pygplates
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import pyplot as plt
from scipy.spatial import cKDTree

from scipy.spatial.transform import Rotation as scpRot

def cartesianToPolarCoords(XYZ, useLonLat=True):
    
    X, Y, Z = XYZ[:, 0], XYZ[:, 1], XYZ[:, 2]
    R = (X ** 2 + Y ** 2 + Z ** 2) ** 0.5
    theta = np.arctan2(Y, X)
    phi = np.arccos(Z / R)

    if useLonLat is True:
        theta, phi = np.degrees(theta), np.degrees(phi)
        lon, lat = theta - 180, 90 - phi
        lon[lon < -180] = lon[lon < -180] + 360
        return lon, lat
    else:
        return R, theta, phi
    
def getCellArea(ds, res):

    dlat = res*np.pi/180.
    dlon = res*np.pi/180.
    a = 6378.137
    e = 0.08181919

    lat = ds.latitude.data
    lon = ds.longitude.data

    mlon, mlat = np.meshgrid(lon, lat)
    mlonlat = np.dstack([mlon.flatten(), mlat.flatten()])[0]

    rlat = mlat.flatten()*np.pi/180.
    cellArea = a**2 * np.cos(rlat) * (1. - e**2) * dlat * dlon / (1. - e**2 * np.sin(rlat)**2)**2

    return cellArea


def polarToCartesian(radius, theta, phi, useLonLat=True):
    if useLonLat:
        theta, phi = np.radians(theta + 180.0), np.radians(90.0 - phi)
    X = radius * np.cos(theta) * np.sin(phi)
    Y = radius * np.sin(theta) * np.sin(phi)
    Z = radius * np.cos(phi)

    # Return data either as a list of XYZ coordinates or as a single XYZ coordinate
    if type(X) == np.ndarray:
        return np.stack((X, Y, Z), axis=1)
    else:
        return np.array([X, Y, Z])

def getPlateIDs(vtime, lonlat, velshare, tree=None):

    # Read plate IDs from gPlates exports
    velfile = velshare+str(int(vtime))+"Ma.xy"
    data = pd.read_csv(
        velfile,
        sep=r"\s+",
        engine="c",
        header=None,
        na_filter=False,
        dtype=float,
        low_memory=False,
    )
    data = data.drop_duplicates().reset_index(drop=True)
    gplateID = data.iloc[:, -1].to_numpy().astype(int)
    if tree is None:
        llvel = data.iloc[:, 0:2].to_numpy()
        tree = cKDTree(llvel)
        dist, ids = tree.query(lonlat, k=1)
        return gplateID[ids], tree
    else:
        dist, ids = tree.query(lonlat, k=1)
        return gplateID[ids]

def quaternion(axis, angle):
        return [
            np.sin(angle / 2) * axis[0],
            np.sin(angle / 2) * axis[1],
            np.sin(angle / 2) * axis[2],
            np.cos(angle / 2),
        ]

def getRotations(time, deltaTime, plateIds, rotation, plate):

    rotationModel = pygplates.RotationModel(rotation)
    topoFeature = pygplates.FeatureCollection(plate)

    rotations = {}
    for plateId in np.unique(plateIds):
        stageRotation = rotationModel.get_rotation(
            int(time - deltaTime), int(plateId), int(time)
        )
        stageRotation = stageRotation.get_euler_pole_and_angle()
        axisLatLon = stageRotation[0].to_lat_lon()
        axis = polarToCartesian(1, axisLatLon[1], axisLatLon[0])
        angle = stageRotation[1]
        rotations[plateId] = scpRot.from_quat(quaternion(axis, angle))
    return rotations

def movePlates(sXYZ, plateIds, rotations):
    newXYZ = np.copy(sXYZ)
    for idx in np.unique(plateIds):
        rot = rotations[idx]
        newXYZ[plateIds == idx] = rot.apply(newXYZ[plateIds == idx])
    return newXYZ

def clusterED(
        time, mvxyz, erodep, output=False, cwd=".",
    ):

    nproc = 4
    clustngbh = 6
    clustdist = 10.0e3

    if output:
        print("\ndbscan MPI")
    dims = [len(mvxyz), 3]
    linepts = mvxyz.ravel()
    lgth = len(linepts)
    fbin = "nodes" + str(time) + ".bin"
    with open(fbin, mode="wb") as f:
        f.write(struct.pack("i" * 2, *[int(i) for i in dims]))
        f.write(struct.pack("f" * (lgth), *[float(i) for i in linepts]))
    fnc = "clusters" + str(time) + ".nc"
    mpi_args = [
        "mpirun",
        "-np",
        str(nproc),
        "dbscan",
        "-i",
        fbin,
        "-b",
        "-m",
        "2",
        "-e",
        str(clustdist),
        "-o",
        fnc,
    ]
    runSubProcess(mpi_args, output, cwd)
    if output:
        print("\nGet global ID of clustered vertices")
    cluster = xr.open_dataset(fnc)
    isClust = cluster.cluster_id.values > 0
    clustPtsX = cluster.position_col_X0.values[isClust]
    clustPtsY = cluster.position_col_X1.values[isClust]
    clustPtsZ = cluster.position_col_X2.values[isClust]
    clustPts = np.vstack((clustPtsX, clustPtsY))
    clustPts = np.vstack((clustPts, clustPtsZ)).T
    ptree = cKDTree(mvxyz)
    dist, ids = ptree.query(clustPts, k=1)
    isCluster = np.zeros(len(mvxyz), dtype=int)
    isCluster[ids] = 1
    idCluster = isCluster > 0
    ptsCluster = mvxyz[idCluster]
    ctree = cKDTree(ptsCluster)
    _, clustNgbhs = ctree.query(ptsCluster, k=clustngbh)
    clustNgbhs = clustNgbhs[:, 1:]
    args = [
        "rm",
        fbin,
        fnc,
    ]
    runSubProcess(args, output, cwd)

    # Get erodep of nearest neighbours
    edInCluster = erodep[idCluster]
    neighbourED = edInCluster[clustNgbhs]

    # For points in cluster, set new heights to the maximum height of
    # nearest neighbours
    clustED = erodep.copy()
    neighbourED.partition(1, axis=1)
    clustED[idCluster] = np.mean(
        neighbourED[:, -int(clustngbh / 2) :], axis=1
    )

    return clustED

def runSubProcess(args, output=True, cwd="."):
    # Launch a subprocess
    p = subprocess.Popen(
        args,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )

    # Capture and re-print OpenMC output in real-time
    lines = []
    while True:
        # If OpenMC is finished, break loop
        line = p.stdout.readline()
        if not line and p.poll() is not None:
            break

        lines.append(line)
        if output:
            # If user requested output, print to screen
            print(line, end="")

    # Raise an exception if return status is non-zero
    if p.returncode != 0:
        # Get error message from output and simplify whitespace
        output = "".join(lines)
        if "ERROR: " in output:
            _, _, error_msg = output.partition("ERROR: ")
        elif "what()" in output:
            _, _, error_msg = output.partition("what(): ")
        else:
            error_msg = "dbscan aborted unexpectedly."
        error_msg = " ".join(error_msg.split())

        raise RuntimeError(error_msg)

def interpData(data, sXYZ, mvxyz, ngbh=1):
    # Build the kdtree
    ptree = cKDTree(mvxyz)
    distNbghs, idNbghs = ptree.query(sXYZ, k=ngbh)
    if ngbh == 1:
        return data[idNbghs]

    # Inverse weighting distance...
    weights = np.divide(
        1.0, distNbghs, out=np.zeros_like(distNbghs), where=distNbghs != 0,
    )
    onIDs = np.where(distNbghs[:, 0] == 0)[0]
    temp = np.sum(weights, axis=1)
    tmp = np.sum(weights * data[idNbghs], axis=1)
    # Data
    interpd = np.divide(tmp, temp, out=np.zeros_like(temp), where=temp != 0)
    if len(onIDs) > 0:
        interpd[onIDs] = data[idNbghs[onIDs, 0]]
    return interpd
