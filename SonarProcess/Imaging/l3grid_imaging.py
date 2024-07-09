import numpy as np
import pandas as pd




def sub2ind(shape, i, j, k):
    return (i * shape[1] + j) * shape[2] + k


def l3grid_imaging(grd_in, grd_out):
    x = grd_in['x'].reshape(-1).astype('float32')
    y = grd_in['y'].reshape(-1).astype('float32')
    z = grd_in['z'].reshape(-1).astype('float32')
    xg1 = grd_out['x']
    yg1 = grd_out['y']
    zg1 = grd_out['z']

    out_shape = xg1.shape

    dx = np.abs(xg1[0][1][0] - xg1[0][0][0])
    dy = np.abs(yg1[1][0][1] - yg1[0][0][0])
    dz = np.abs(zg1[0][0][1] - zg1[0][0][0])

    # The following three lines are added to ensure correct indexing in the
    # calculation of interpolation weight
    xg = np.squeeze(xg1[0, :, 0])
    yg = np.squeeze(yg1[:, 0, 0])
    zg = np.squeeze(zg1[0, 0, :])

    xmin = np.min(xg)
    ymin = np.min(yg)
    zmin = np.min(zg)
    xmax = np.max(xg)
    ymax = np.max(yg)
    zmax = np.max(zg)

    Ia = grd_in['Ia'].reshape(-1).astype('float32')
    Id = grd_in['Id'].reshape(-1).astype('float32')
    Ia_filt = grd_in['Ia_filt'].reshape(-1).astype('float32')
    Id_filt = grd_in['Id_filt'].reshape(-1).astype('float32')
    Kp = grd_in['Kp'].reshape(-1).astype('float32')

    # use only finite points within the grid
    ii = np.where(
        np.isfinite(Id) & (Id > 1e-9) & (Ia > 1e-9) & np.isfinite(Ia) & (x >= xmin) & (x <= xmax) & (y >= ymin) & (
                    y <= ymax) & (z >= zmin) & (z <= zmax))[0]

    if ii.shape[0] > 0:
        x = x[ii]
        y = y[ii]
        z = z[ii]
        Ia = Ia[ii]
        Id = Id[ii]
        Ia_filt = Ia_filt[ii]
        Id_filt = Id_filt[ii]
        Kp = Kp[ii]

        # grd_in = {
        #         'x' : x,        #.reshape(-1).astype('float32'),
        #         'y' : y,        #.reshape(-1).astype('float32'),
        #         'z' : z,        #.reshape(-1).astype('float32'),
        #         'Ia' : Ia,
        #         'Id' : Id,        #.reshape(-1).astype('float32'),
        #         # 'Ia_filt' : Ia_filt,
        #         # 'Id_filt' : Id_filt,
        #         'Kp' : Kp
        #     }

        i = np.array([np.floor((x - xmin) / dx),
                      np.ceil((x - xmin) / dx)]).astype(np.int32)
        j = np.array([np.floor((y - ymin) / dy),
                      np.ceil((y - ymin) / dy)]).astype(np.int32)
        k = np.array([np.floor((z - zmin) / dz),
                      np.ceil((z - zmin) / dz)]).astype(np.int32)

        grd_out['Ia'] = grd_out['Ia'].reshape(-1)
        grd_out['Id'] = grd_out['Id'].reshape(-1)
        grd_out['Ia_filt'] = grd_out['Ia_filt'].reshape(-1)
        grd_out['Id_filt'] = grd_out['Id_filt'].reshape(-1)
        grd_out['Kp'] = grd_out['Kp'].reshape(-1)
        grd_out['w'] = grd_out['w'].reshape(-1)

        for n in range(2):
            wx = 1 - np.abs(xg[i[n, :]] - x) / dx
            for m in range(2):
                wy = 1 - np.abs(yg[j[m, :]] - y) / dy
                for l in range(2):
                    wz = 1 - np.abs(zg[k[l, :]] - z) / dz
                    w = np.sqrt(np.power(wx, 2) + np.power(wy, 2) + np.power(wz, 2))
                    p = sub2ind(out_shape, j[m, :], i[n, :], k[l, :])
                    grd_out['Ia'][p] = grd_out['Ia'][p] + w * Ia  # averaged backscatter intensity
                    grd_out['Id'][p] = grd_out['Id'][p] + w * Id  # differenced backscatter intensity
                    grd_out['Ia_filt'][p] = grd_out['Ia_filt'][
                                                p] + w * Ia_filt  # filtered averaged backscatter intensity
                    grd_out['Id_filt'][p] = grd_out['Id_filt'][
                                                p] + w * Id_filt  # filtered differenced backscatter intensity
                    grd_out['Kp'][p] = grd_out['Kp'][p] + w * Kp
                    grd_out['w'][p] = grd_out['w'][p] + w  # weight function

        grd_out['Ia'] = grd_out['Ia'].reshape(out_shape)
        grd_out['Id'] = grd_out['Id'].reshape(out_shape)
        grd_out['Ia_filt'] = grd_out['Ia_filt'].reshape(out_shape)
        grd_out['Id_filt'] = grd_out['Id_filt'].reshape(out_shape)
        grd_out['Kp'] = grd_out['Kp'].reshape(out_shape)
        grd_out['w'] = grd_out['w'].reshape(out_shape)

    return grd_out
