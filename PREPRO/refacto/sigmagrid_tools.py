from progressbar import progressbar
import numpy as np
import toolsf

def scoord2z(point_type, zeta, topo,theta_s, theta_b,N,hc,scoord='new2008'):
    '''
    scoord2z finds z at either rho or w points (positive up, zero at rest surface)

    Inputs:
      point_type        'r' or 'w'
      zeta               sea surface height
      topo              array of depths (e.g., from grd file)
      theta_s           surface focusing parameter
      theta_b           bottom focusing parameter
      N                 number of vertical rho-points
      hc                critical depth
      scoord            'new2008' :new scoord 2008  or 'old1994' for Song scoord
    
    Outputs:
      z                 sigma coordinates
      Cs                Cs parameter
    '''
    def CSF(sc,theta_s,theta_b):
        '''
        Allows use of theta_b > 0 (July 2009)
        '''
        one64 = np.float64(1)
        if theta_s > 0.:
            csrf = ((one64 - np.cosh(theta_s * sc))
                       / (np.cosh(theta_s) - one64))
        else:
            csrf = -sc ** 2
        sc1 = csrf + one64
        if theta_b > 0.:
            Cs = ((np.exp(theta_b * sc1) - one64)
                / (np.exp(theta_b) - one64) - one64)
        else:
            Cs = csrf
        return Cs

    N = np.float64(N)
    cff1 = 1. / np.sinh(theta_s)
    cff2 = 0.5 / np.tanh(0.5 * theta_s)
    sc_w = (np.arange(N + 1, dtype=np.float64) - N) / N
    sc_r = ((np.arange(1, N + 1, dtype=np.float64)) - N - 0.5) / N

    if 'w' in point_type:
        sc = sc_w
        N += 1. # add a level
    else:
        sc = sc_r
        
    if len(np.array(zeta).shape)>2: # case zeta is 3-D (in time)
        z  = np.empty((int(zeta.shape[0]),) + (int(N),) + topo.shape, dtype=np.float64)
    else:
        z  = np.empty((int(N),) + topo.shape, dtype=np.float64)

    if scoord in 'new2008':
        Cs = CSF(sc,theta_s,theta_b)
    elif scoord in 'old1994':
        Cs = (1. - theta_b) * cff1 * np.sinh(theta_s * sc) + \
           theta_b * (cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5)

    if scoord in 'new2008':
        hinv = 1. / (topo + hc)
        cff = hc * sc
        cff1 = Cs
            
        if len(np.array(zeta).shape)>2:
            for t in range(zeta.shape[0]):
                if 'w' in point_type:
                    z[t,0]=-topo
                    start=1
                else:
                    start=0
                for k in np.arange(start,N, dtype=int):
                    z[t,k] = zeta[t] + (zeta[t] + topo) * (cff[k] + cff1[k] * topo) * hinv
        else:
            for k in np.arange(N, dtype=int):
                z[k] = zeta + (zeta + topo) * (cff[k] + cff1[k] * topo) * hinv

    elif scoord in 'old1994':
        hinv = 1. / topo
        cff  = hc * (sc - Cs)
        cff1 = Cs
        cff2 = sc + 1

        if len(np.array(zeta).shape)>2:
            for t in range(zeta.shape[0]):
                for k in np.arange(N,dtype=int) + 1:
                    z0      = cff[k-1] + cff1[k-1] * topo
                    z[t,k-1, :] = z0 + zeta[t,:] * (1. + z0 * hinv)
        else:
            for k in np.arange(N,dtype=int) + 1:
                z0      = cff[k-1] + cff1[k-1] * topo
                z[k-1, :] = z0 + zeta * (1. + z0 * hinv)
    else:
        raise Exception("Unknown scoord, should be 'new2008' or 'old1994'")

    if sc_r is None:
        sc_r = sc_r
    return z.squeeze(), np.float32(Cs)

####################

def ztosigma(vin,Z,zcroco):
    '''
    This fonction perform the z to sigma transformation for
    3D (Z,Y,X) or 4D (T,Z,Y,X) variables

    Input:
      Vin       Input variables to put on sigma grid (3 or 4D)
      Z         Input depth values (1D)
      zcroco    CROCO vertical level (3D)

    output:
       vout     Vin projection on zcroco
    '''
# Do a vertical interpolation from z levels to sigma CROCO levels
    if len(zcroco.shape)>3:
        [T,N,M,L]=np.shape(zcroco)
        four_dim=True
    else:
        [N,M,L]=np.shape(zcroco)
        four_dim=False

    [Nz]=np.shape(Z)
#
# Find the grid position of the nearest vertical levels
#
    i1=np.arange(0,L)
    j1=np.arange(0,M)
    if four_dim:
        t1=np.arange(0,T)
        [jmat,tmat,imat]=np.meshgrid(j1,t1,i1)
        VAR=np.reshape(vin,T*Nz*M*L)
        vout=np.zeros((T,N,M,L))
    else:
        [imat,jmat]=np.meshgrid(i1,j1)
        VAR=np.reshape(vin,Nz*M*L)
        vout=np.zeros((N,M,L))

    for ks in progressbar(range(N),' Sigma layer : ', 40):
        if four_dim:
            sigmalev=zcroco[:,ks,:,:]
            thezlevs=np.zeros((T,M,L),dtype=int)-1
        else:
            sigmalev=zcroco[ks,:,:]
            thezlevs=np.zeros((M,L),dtype=int)-1

        for kz in range(Nz):
            thezlevs[np.where(sigmalev>Z[kz])]=thezlevs[np.where(sigmalev>Z[kz])]+1

        if four_dim:
            pos= L*M*Nz*tmat+ L*M*thezlevs + L*jmat + imat
        else:
            pos= L*M*thezlevs + L*jmat + imat

        z1=Z[thezlevs]
        z2=Z[thezlevs+1]
        v1=VAR[pos]
        v2=VAR[pos+L*M]

        if four_dim:
            vout[:,ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
        else:
            vout[ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
    return vout

####################
