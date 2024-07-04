import numpy as np

def sub2ind(shape, i, j):
    return i * shape[1] + j

def l2grid(x,y,v,xg,yg,vg,wg):
    sp = xg.shape
    dx = np.abs(xg[0,1]-xg[0,0])
    dy = np.abs(yg[1,0]-yg[0,0])

    # The following two lines are added to ensure correct indexing in the
    # calculation of interpolation weight
    xg = xg[0,:] 
    yg = yg[:,0]

    xmin = np.min(xg)
    ymin = np.min(yg)
    xmax = np.max(xg)
    ymax = np.max(yg)

    x = x.reshape(-1)
    y = y.reshape(-1)
    v = v.reshape(-1)

    vg = vg.reshape(-1)
    wg = wg.reshape(-1)

    # use only points within the grid
    ii = np.where((x>=xmin) & (x<=xmax) & (y>=ymin) & (y<=ymax) & (np.isfinite(v)) & (np.isfinite(x)))[0]
    if ii.shape[0] == 0:
        return
    x=x[ii]
    y=y[ii]
    v=v[ii]

    i = np.array([
        np.floor((x-xmin)/dx),
        np.ceil((x-xmin)/dx)
    ]).astype(np.int32)
    j = np.array([
        np.floor((y-ymin)/dy),
        np.ceil((y-ymin)/dy)
    ]).astype(np.int32)

    for n in range(2):
       wx = 1-np.abs(xg[i[n,:]]-x)/dx
       for m in range(2):
          wy = 1-np.abs(yg[j[m,:]]-y)/dy
          w = np.sqrt(np.power(wx,2)+np.power(wy,2))
          k = sub2ind(sp,j[n,:],i[m,:])
          vg[k] = vg[k] + w*v
          wg[k] = wg[k] + w

    vg = vg.reshape(sp)
    wg = wg.reshape(sp)

    return vg, wg
