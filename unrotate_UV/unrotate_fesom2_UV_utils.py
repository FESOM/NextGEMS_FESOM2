import math as mt
import numpy as np

def compute_face_coords(mesh):
    """ Compute coordinates of elements (triangles)
    Parameters:
    -----------
    mesh: mesh object
        fesom mesh object
    Returns:
    --------
    face_x: numpy array
        x coordinates
    face_y: numpy array
        y coordinates
    """
    first_mean = mesh.x2[mesh.elem].mean(axis=1)
    j = np.where(np.abs(mesh.x2[mesh.elem][:, 0] - first_mean) > 100)[0]
    cyclic_elems = mesh.x2[mesh.elem][j].copy()
    new_means = np.where(cyclic_elems > 0, cyclic_elems, cyclic_elems + 360).mean(
        axis=1
    )
    new_means[new_means > 180] = new_means[new_means > 180] - 360
    first_mean[j] = new_means
    face_x = first_mean
    face_y = mesh.y2[mesh.elem].mean(axis=1)
    return face_x, face_y

def scalar_r2g(al, be, ga, rlon, rlat):
    """
    Converts rotated coordinates to geographical coordinates.
    Parameters
    ----------
    al : float
        alpha Euler angle
    be : float
        beta Euler angle
    ga : float
        gamma Euler angle
    rlon : array
        1d array of longitudes in rotated coordinates
    rlat : array
        1d araay of latitudes in rotated coordinates
    Returns
    -------
    lon : array
        1d array of longitudes in geographical coordinates
    lat : array
        1d array of latitudes in geographical coordinates
    """

    rad = mt.pi / 180
    al = al * rad
    be = be * rad
    ga = ga * rad
    rotate_matrix = np.zeros(shape=(3, 3))
    rotate_matrix[0, 0] = np.cos(ga) * np.cos(al) - np.sin(ga) * np.cos(be) * np.sin(al)
    rotate_matrix[0, 1] = np.cos(ga) * np.sin(al) + np.sin(ga) * np.cos(be) * np.cos(al)
    rotate_matrix[0, 2] = np.sin(ga) * np.sin(be)
    rotate_matrix[1, 0] = -np.sin(ga) * np.cos(al) - np.cos(ga) * np.cos(be) * np.sin(
        al
    )
    rotate_matrix[1, 1] = -np.sin(ga) * np.sin(al) + np.cos(ga) * np.cos(be) * np.cos(
        al
    )
    rotate_matrix[1, 2] = np.cos(ga) * np.sin(be)
    rotate_matrix[2, 0] = np.sin(be) * np.sin(al)
    rotate_matrix[2, 1] = -np.sin(be) * np.cos(al)
    rotate_matrix[2, 2] = np.cos(be)

    rotate_matrix = np.linalg.pinv(rotate_matrix)

    rlat = rlat * rad
    rlon = rlon * rad

    # Rotated Cartesian coordinates:
    xr = np.cos(rlat) * np.cos(rlon)
    yr = np.cos(rlat) * np.sin(rlon)
    zr = np.sin(rlat)

    # Geographical Cartesian coordinates:
    xg = rotate_matrix[0, 0] * xr + rotate_matrix[0, 1] * yr + rotate_matrix[0, 2] * zr
    yg = rotate_matrix[1, 0] * xr + rotate_matrix[1, 1] * yr + rotate_matrix[1, 2] * zr
    zg = (
        rotate_matrix[2, 0] * xr + rotate_matrix[2, 1] * yr + rotate_matrix[2, 2] * zr
    )

    # Geographical coordinates:
    lat = np.arcsin(zg)
    lon = np.arctan2(yg, xg)

    a = np.where((np.abs(xg) + np.abs(yg)) == 0)
    if a:
        lon[a] = 0

    lat = lat / rad
    lon = lon / rad

    return (lon, lat)

def vec_rotate_r2g(al, be, ga, lon, lat, urot, vrot, flag):
    """
    Rotate vectors from rotated coordinates to geographical coordinates.
    Parameters
    ----------
    al : float
        alpha Euler angle
    be : float
        beta Euler angle
    ga : float
        gamma Euler angle
    lon : array
        1d array of longitudes in rotated or geographical coordinates (see flag parameter)
    lat : array
        1d array of latitudes in rotated or geographical coordinates (see flag parameter)
    urot : array
        1d array of u component of the vector in rotated coordinates
    vrot : array
        1d array of v component of the vector in rotated coordinates
    flag : 1 or 0
        flag=1  - lon,lat are in geographical coordinate
        flag=0  - lon,lat are in rotated coordinate
    Returns
    -------
    u : array
        1d array of u component of the vector in geographical coordinates
    v : array
        1d array of v component of the vector in geographical coordinates
    """

    #   first get another coordinate
    if flag == 1:
        (rlon, rlat) = scalar_g2r(al, be, ga, lon, lat)
    else:
        rlon = lon
        rlat = lat
        (lon, lat) = scalar_r2g(al, be, ga, rlon, rlat)

    #   then proceed...
    rad = mt.pi / 180
    al = al * rad
    be = be * rad
    ga = ga * rad

    rotate_matrix = np.zeros(shape=(3, 3))
    rotate_matrix[0, 0] = np.cos(ga) * np.cos(al) - np.sin(ga) * np.cos(be) * np.sin(al)
    rotate_matrix[0, 1] = np.cos(ga) * np.sin(al) + np.sin(ga) * np.cos(be) * np.cos(al)
    rotate_matrix[0, 2] = np.sin(ga) * np.sin(be)
    rotate_matrix[1, 0] = -np.sin(ga) * np.cos(al) - np.cos(ga) * np.cos(be) * np.sin(
        al
    )
    rotate_matrix[1, 1] = -np.sin(ga) * np.sin(al) + np.cos(ga) * np.cos(be) * np.cos(
        al
    )
    rotate_matrix[1, 2] = np.cos(ga) * np.sin(be)
    rotate_matrix[2, 0] = np.sin(be) * np.sin(al)
    rotate_matrix[2, 1] = -np.sin(be) * np.cos(al)
    rotate_matrix[2, 2] = np.cos(be)

    rotate_matrix = np.linalg.pinv(rotate_matrix)
    rlat = rlat * rad
    rlon = rlon * rad
    lat = lat * rad
    lon = lon * rad

    #   vector in rotated Cartesian
    txg = -vrot * np.sin(rlat) * np.cos(rlon) - urot * np.sin(rlon)
    tyg = -vrot * np.sin(rlat) * np.sin(rlon) + urot * np.cos(rlon)
    tzg = vrot * np.cos(rlat)

    #   vector in geo Cartesian
    txr = (
        rotate_matrix[0, 0] * txg
        + rotate_matrix[0, 1] * tyg
        + rotate_matrix[0, 2] * tzg
    )
    tyr = (
        rotate_matrix[1, 0] * txg
        + rotate_matrix[1, 1] * tyg
        + rotate_matrix[1, 2] * tzg
    )
    tzr = (
        rotate_matrix[2, 0] * txg
        + rotate_matrix[2, 1] * tyg
        + rotate_matrix[2, 2] * tzg
    )

    #   vector in geo coordinate
    v = (
        -np.sin(lat) * np.cos(lon) * txr
        - np.sin(lat) * np.sin(lon) * tyr
        + np.cos(lat) * tzr
    )
    u = -np.sin(lon) * txr + np.cos(lon) * tyr

    u = np.array(u)
    v = np.array(v)

    return (u, v)

def scalar_g2r(al, be, ga, lon, lat):
    """
    Converts geographical coordinates to rotated coordinates.
    Parameters
    ----------
    al : float
        alpha Euler angle
    be : float
        beta Euler angle
    ga : float
        gamma Euler angle
    lon : array
        1d array of longitudes in geographical coordinates
    lat : array
        1d array of latitudes in geographical coordinates
    Returns
    -------
    rlon : array
        1d array of longitudes in rotated coordinates
    rlat : array
        1d araay of latitudes in rotated coordinates
    """

    rad = mt.pi / 180
    al = al * rad
    be = be * rad
    ga = ga * rad

    rotate_matrix = np.zeros(shape=(3, 3))

    rotate_matrix[0, 0] = np.cos(ga) * np.cos(al) - np.sin(ga) * np.cos(be) * np.sin(al)
    rotate_matrix[0, 1] = np.cos(ga) * np.sin(al) + np.sin(ga) * np.cos(be) * np.cos(al)
    rotate_matrix[0, 2] = np.sin(ga) * np.sin(be)
    rotate_matrix[1, 0] = -np.sin(ga) * np.cos(al) - np.cos(ga) * np.cos(be) * np.sin(
        al
    )
    rotate_matrix[1, 1] = -np.sin(ga) * np.sin(al) + np.cos(ga) * np.cos(be) * np.cos(
        al
    )
    rotate_matrix[1, 2] = np.cos(ga) * np.sin(be)
    rotate_matrix[2, 0] = np.sin(be) * np.sin(al)
    rotate_matrix[2, 1] = -np.sin(be) * np.cos(al)
    rotate_matrix[2, 2] = np.cos(be)

    lat = lat * rad
    lon = lon * rad

    # geographical Cartesian coordinates:
    xr = np.cos(lat) * np.cos(lon)
    yr = np.cos(lat) * np.sin(lon)
    zr = np.sin(lat)

    # rotated Cartesian coordinates:
    xg = rotate_matrix[0, 0] * xr + rotate_matrix[0, 1] * yr + rotate_matrix[0, 2] * zr
    yg = rotate_matrix[1, 0] * xr + rotate_matrix[1, 1] * yr + rotate_matrix[1, 2] * zr
    zg = rotate_matrix[2, 0] * xr + rotate_matrix[2, 1] * yr + rotate_matrix[2, 2] * zr

    # rotated coordinates:
    rlat = np.arcsin(zg)
    rlon = np.arctan2(yg, xg)

    a = np.where((np.abs(xg) + np.abs(yg)) == 0)
    if a:
        lon[a] = 0

    rlat = rlat / rad
    rlon = rlon / rad

    return (rlon, rlat)