import scipy as sp
import numpy as np
from scipy import random
PI = np.pi


def ComputeXYZdeg(ra,dec,R,ra0,dec0) :
    return ComputeXYZ(np.radians(ra),np.radians(dec),R,np.radians(ra0),np.radians(dec0))

def ComputeXYZ(ra,dec,R,ra0,dec0) :
    '''
    XYZ of a point P (ra,dec,R) in a frame with
    observer at O, Z along OP, X along ra0, Y along dec0

    angles in radians
    tested that ra,dec, R = box.ComputeRaDecR(R0,ra0,dec0,X,Y,Z)
    x,y,z = box.ComputeXYZ(ra[0],dec[0],R,ra0,dec0)
    print(x-X,y-Y,z-R0-Z)#        prints ~1E-13  for random inputs
    '''

    theta0 = PI/2 - dec0 # polar angle or colatitude
                #......         unit vector in X, Y and Z directions
    cosra0 = np.cos(ra0)
    sinra0 = np.sin(ra0)
    costheta0 = np.cos(theta0)
    sintheta0 = np.sin(theta0)
    xVec = np.array([-sinra0, cosra0 ,0])
    yVec = np.array([-costheta0 * cosra0, -costheta0 * sinra0 , sintheta0])
    zVec = np.array([sintheta0 * cosra0, sintheta0 * sinra0, costheta0])
                # coordinates of OP
    theta = PI/2 - dec
    cosra = np.cos(ra)
    sinra = np.sin(ra)
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    OP = R * np.array([sintheta*cosra,sintheta*sinra,costheta])

    #print('op',OP)
    #print(xVec,yVec,zVec)
    X = np.dot(OP,xVec)
    Y = np.dot(OP,yVec)
    Z = np.dot(OP,zVec)
    return X,Y,Z

def ComputeRaDecR(R0,ra0,dec0,X,Y,Z) :
    '''
    ra, dec, R of P(X,Y,Z) in a box with center located in C(R0,ra0,dec0)

    angles in radian -pi < ra < pi   -pi/2 < dec < pi/2
    R0, ra0, dec0 define the position of the box in a (U,V,W) frame
        centered at the oberver position O
    box oriented such that Z along OC, X along ra, Y along dec
    R0,ra0,dec0 are scalar   X, Y, Z can be scalars, 1D or 2D arrays
    '''
    theta0 = np.pi/2 - dec0 # polar angle or colatitude
                #......         unit vector in X, Y and Z directions
    cosra0 = np.cos(ra0)
    sinra0 = np.sin(ra0)
    costheta0 = np.cos(theta0)
    sintheta0 = np.sin(theta0)
    # unit vector along X Y and Z in the (U,V,W) frame
    xVec = np.array([-sinra0, cosra0 ,0])[:, None,None] # (3,1,1)
    yVec = np.array([-costheta0 * cosra0, -costheta0 * sinra0 , sintheta0])[:, None,None]
    zVec = np.array([sintheta0 * cosra0, sintheta0 * sinra0, costheta0])[:, None,None]
    #print(xVec.shape, yVec.shape, zVec.shape)
    #print(X.shape, Y.shape, Z.shape)
                #......         vector OP
    OPVec = X * xVec + Y * yVec + (R0+Z) * zVec # broadcasting due to [:, None,None] above
                    # allows to run with X and Y that are 1D or 2D arrays
    U = OPVec[0]    # (1,1) (1,N) or (N1,N2) depending on X,Y,Z shape
    V = OPVec[1]
    W = OPVec[2]
                #......         ra dec of OP
                # arctan2(y,x): angle relative to x axis in [-pi,pi]
    ra = np.arctan2(V,U)    # angle raltive to U
    dec = np.arctan2(W, np.sqrt(U*U + V*V))  # angle relative to (u,V) plane
                                        # in [-pi2/,pi/2] since sqrt(U*U+V*V) > 0
    R = np.sqrt(U*U+V*V+W*W)
    #print(ra.shape[0],ra.shape[1])
    if (np.isscalar(X+Y+Z)):
        return ra[0,0],dec[0,0],R[0,0]
    #if (np.isscalar((X+Y+Z)[0])) : # X,Y,Z 1D arrays
    #   we may get X+Y+Z a zero length array, in which case, this fails
    if (ra.shape[0]==1) : # X,Y,Z 1D arrays
        return ra[0],dec[0],R[0]
    return ra, dec, R

def Compute_cos_min(LX,LY,Rmax) :
    #
    #  angle for a point at X=Xmax and Rmax from observer
    sinx_max = (LX/2) / Rmax
    cosx_min = np.sqrt(1 - sinx_max**2)
    tanx_max = sinx_max / cosx_min
    #  angle for a point at Y=Ymax and Rmax from observer
    siny_max = (LY/2) / Rmax
    cosy_min = np.sqrt(1 - siny_max**2)
    tany_max = siny_max / cosy_min
    #
    sin_max = (np.sqrt(LX*LX+LY*LY)/2) / Rmax
    cos_min = np.sqrt(1 - sin_max**2)
    return tanx_max,tany_max,cos_min

def box_center(LX,LY,Rmax,Rmin,margin) :
    # we want to include a cone Rmin < R < Rmax in a LX x LY box minus margin
    # returns the resulting tan_max, the required LZ and the box_center distance
    tanx_max,tany_max,cos_min = Compute_cos_min(LX-2*margin,LY-2*margin,Rmax)
    LZ = Rmax - (cos_min * Rmin) + 2 * margin
    R0 = Rmax -LZ/2 + margin
    return R0,LZ,tanx_max,tany_max


def box_limitOld(LX,LY,LZ,R0,margin) :
    # compute R_min, R_max, tanx_max and tany_max
    # to be inside the box with a margin
    # for a box located at R0 from the observer
    #    we want that in any direction R_min < R_QSO < R_max
    Rmax = R0 + LZ/2  - margin
    #  angle for a point at X=Xmax and Rmax from observer
    sinx_max = (LX/2 -margin) / Rmax
    cosx_min = np.sqrt(1 - sinx_max**2)
    tanx_max = sinx_max / cosx_min
#    print(sinx_max, tanx_max)    # prov
    siny_max = (LY/2 -margin) / Rmax
    cosy_min = np.sqrt(1 - siny_max**2)
    tany_max = siny_max / cosy_min
    cos_min = min(cosx_min,cosy_min)
#    sin_max = max(sinx_max,siny_max)
#    cos_min = np.sqrt(1 - sin_max**2)
    Rmin = (R0 - LZ/2 + margin) / cos_min

#    print(R0-LZ/2, Rmin, Rmax)
    #print(0, Rmin-R0+LZ/2, Rmax-R0+LZ/2)
    #print(sin_alpha)
    return Rmin,Rmax,tanx_max,tany_max

def box_limit(LX,LY,LZ,R0,margin) :
    # compute R_min, R_max, tanx_max and tany_max
    # to be inside the box with a margin
    # for a box located at R0 from the observer
    #    we want that in any direction R_min < R_QSO < R_max
    Rmax = R0 + LZ/2  - margin
    #  angle for a point at X=Xmax and Rmax from observer
    if (LX != LY) :
        print("case LX =",LX,"!= LY=",LY,"not implemented yet")
        exit(0)
    sinx_max = (LX/2 -margin) / Rmax
    cosx_min = np.sqrt(1 - sinx_max**2)
    tanx_max = sinx_max / cosx_min
#    print(sinx_max, tanx_max)    # prov
    siny_max = (LY/2 -margin) / Rmax
    cosy_min = np.sqrt(1 - siny_max**2)
    tany_max = siny_max / cosy_min

    sin_max = (np.sqrt(LX*LX+LY*LY)/2 -margin) / Rmax
    cos_min = np.sqrt(1 - sin_max**2)
    Rmin = (R0 - LZ/2 + margin) / cos_min

#    print(R0-LZ/2, Rmin, Rmax)
    #print(0, Rmin-R0+LZ/2, Rmax-R0+LZ/2)
    #print(sin_alpha)
    return Rmin,Rmax,tanx_max,tany_max



def sample_box(xbox,ybox,zbox,rho,threshold=2.5):
#                                                               obsolete
#
    dmax = 3    # we compute rho over +- dmax cells  should be a parameter <==

#.......................    select pixels with rho > threshold*rms
    rms = rho.std()
    iqso = sp.where(rho > threshold*rms)
#    print("threshold*rms=", threshold*rms)
    iqso = sp.array(iqso)        # (3,N_QSO) 3 raws N col
#.......................   remove QSO at less than dmax cells from the box edges
    imax = xbox.size - dmax
    jmax = ybox.size - dmax
    kmax = zbox.size - dmax
    jqso = iqso.T
    jqso = jqso[jqso[:,0]>dmax]
    jqso = jqso[jqso[:,0]<imax]
    jqso = jqso[jqso[:,1]>dmax]
    jqso = jqso[jqso[:,1]<jmax]
    jqso = jqso[jqso[:,2]>dmax]
    jqso = jqso[jqso[:,2]<kmax]
    iqso = jqso.T
#    print(iqso)
#    print(sp.where(iqso))

#   needs to implement z distribution
#   taking into account that z depends not only on Z but also on X and Y
#   must also select in a cone, not in the full box, this requires to know z0
#   must randomize the position inside the cell

    print("found: ",len(iqso[0])," quasars")

#    return [x+y+z for x,y,z in zip(xbox[iqso[0]],ybox[iqso[1]],zbox[iqso[2]])]
    return [qso(x,y,z) for x,y,z in zip(xbox[iqso[0]],ybox[iqso[1]],zbox[iqso[2]])]


def ComputeRaDecR2(x, y, z, ra0, dec0):
    # all angles should be in radians
    # x, y, z are arrays of same dimension, and ra0 dec0 floats
    # Compute (ra,dec,R0) in base R of a point P(x,y,z) where (x,y,z)
    # are the coordinates in a base R' rotated by ra0 and dec0 with
    # respect to base R :
    # P_in_R' = Rot(dec0).Rot(ra0).P_in_R
    numra = (np.cos(ra0)*x - np.sin(dec0)*np.sin(ra0)*y
             + np.cos(dec0)*np.sin(ra0)*z)
    denomra = (-np.sin(ra0)*x - np.sin(dec0)*np.cos(ra0)*y
               + np.cos(dec0)*np.cos(ra0)*z)
    numdec = np.cos(dec0)*y + np.sin(dec0)*z
    R = np.sqrt(x**2 + y**2 + z**2)

    ra = np.zeros(x.shape)
    msk = np.where((numra > 0) & (denomra > 0))
    ra[msk] = np.arctan(numra[msk] / denomra[msk])
    msk = np.where((numra > 0) & (denomra < 0))
    ra[msk] = np.arctan(numra[msk] / denomra[msk]) + np.pi
    msk = np.where((numra < 0) & (denomra < 0))
    ra[msk] = np.arctan(numra[msk] / denomra[msk]) + np.pi
    msk = np.where((numra < 0) & (denomra > 0))
    ra[msk] = np.arctan(numra[msk] / denomra[msk]) + 2*np.pi

    msk = np.where((numra == 0) & (denomra > 0))
    ra[msk] = 0
    msk = np.where((numra > 0) & (denomra == 0))
    ra[msk] = np.pi / 2
    msk = np.where((numra == 0) & (denomra < 0))
    ra[msk] = np.pi
    msk = np.where((numra < 0) & (denomra == 0))
    ra[msk] = 3*np.pi / 2
    msk = np.where((numra == 0) & (denomra == 0))
    ra[msk] = 0

    dec = np.arcsin(numdec/R)
    return ra, dec, R


def ComputeXYZ2(ra, dec, R, ra0, dec0):
    # all angles should be in radians
    x = R * (np.cos(ra0)*np.cos(dec)*np.sin(ra)
             - np.sin(ra0)*np.cos(dec)*np.cos(ra))
    y = R * (-np.sin(ra0)*np.sin(dec0)*np.cos(dec)*np.sin(ra)
             + np.cos(dec0)*np.sin(dec)
             - np.cos(ra0)*np.sin(dec0)*np.cos(dec)*np.cos(ra))
    z = R * (np.cos(dec0)*np.sin(ra0)*np.cos(dec)*np.sin(ra)
             + np.sin(dec0)*np.sin(dec)
             + np.cos(ra0)*np.cos(dec0)*np.cos(dec)*np.cos(ra))
    return x, y, z
