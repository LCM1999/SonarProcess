import numpy as np

def covis_coords_darrell(origin,rng,azim,yaw,roll,pitch,central_head,declination):
    # COVIS central pointing direction re true north
    ang_rot = np.pi/180*(central_head + declination - 360)
    
    # rotation matrix for central heading
    ROT = np.array([[np.cos(ang_rot), np.sin(ang_rot), 0],
                    [-np.sin(ang_rot), np.cos(ang_rot), 0 ],
                    [0, 0, 1]])
    
    # Matrices giving coordinate transformations
    R = np.array([[np.cos(roll), 0, -np.sin(roll)],
                  [0, 1, 0],
                  [np.sin(roll), 0, np.cos(roll)]])
    Y = np.array([[np.cos(yaw), -np.sin(yaw), 0],
                  [np.sin(yaw), np.cos(yaw), 0],
                  [0, 0, 1]])
    P = np.array([[1, 0, 0],
                  [0, np.cos(pitch), np.sin(pitch)],
                  [0, -np.sin(pitch), np.cos(pitch)]])
    
    # Receiver coordinates
    xr = np.array([rng]).T @ np.sin([azim])
    yr = np.array([rng]).T @ np.cos([azim])
    zr = np.zeros(xr.shape)
    rr = np.array([
        xr.T.ravel(),
        yr.T.ravel(),
        zr.T.ravel()
    ])
    
    # World coordinates
    rw = np.linalg.solve((R@P@Y), rr)
    rw = ROT@rw
    #xw = rw[0,:].reshape(len(rng),len(azim))
    xw = rw[0,:].reshape(len(azim),len(rng)).T
    yw = rw[1,:].reshape(len(azim),len(rng)).T
    zw = rw[2,:].reshape(len(azim),len(rng)).T

    x0 = origin[0]
    y0 = origin[1]
    z0 = origin[2]

    xw = xw + x0
    yw = yw + y0
    zw = zw + z0

    return xw, yw, zw