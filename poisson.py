import ampere_pole,grid,matrices
from numpy import linalg as la
import numpy as np
import sys
import time
import point,sphtrig
from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid,zeros,arange,ones_like,array,roll,hstack
from matplotlib.pyplot import figure,subplot,pcolormesh,colorbar,contour,show,contourf,rgrids,ylim
from scipy import interpolate
import time

if __name__ == "__main__":
    Sigma_const = 8. # scalar for the time being
    interp = True

    # AMPERE stuff
    amp = ampere_pole.ampere_pole('data/AMPERE_north_2015-08-15T20-50-00Z.txt')

    grd = grid.generate_from_data(amp.get_grid())
    jr = amp.get_data()
    if interp:
        ntheta = 100
        nphi = 720

        # complete full circle with AMPERE data so that interpolated grid doesn't have a huge hole
        # note, however, that we still perform the solve on a cell-centered grid, i.e., in the case of interpolation
        # here we make sure that the resultant grid doesn't wrap around
        p_tmp = hstack((grd.p,2.*pi+grd.p[:,[0]]))
        t_tmp = hstack((grd.t,grd.t[:,[0]]))
        jr_tmp = hstack((jr,jr[:,[0]]))

        f = interpolate.RectBivariateSpline(grd.t[:,0],grd.p[0,:],jr,kx=1,ky=1)  # note, assume uniform ampere grid here
        # this here defines cell vortices (and wraps around)
        new_t = linspace(t_tmp.min(),t_tmp.max(),ntheta)
        new_p = linspace(p_tmp.min(),p_tmp[0,:].max(),nphi+1)
        # now redefine as cell centers
        new_p = 0.5*(new_p[:-1]+new_p[1:])
        grd = grid.grid()
        grd.t,grd.p = meshgrid(new_t,new_p,indexing='ij')
        jr = f(grd.t[:,0],grd.p[0,:])

    jr[abs(jr)<=0.2]=0.

    dt = grd.t[1:,:]-grd.t[:-1,:]
    dp = grd.p[:,1:]-grd.p[:,:-1]

    (Nt,Np) = grd.t.shape
    ############

    # define conductance as uniform Pedersen for now
    SigmaP = ones_like(jr)*Sigma_const
    cosd     = -2.*cos(grd.t)/sqrt(1.+3.*cos(grd.t)**2)
    F11          = sin(grd.t)*SigmaP/cosd**2
    F22          = SigmaP


    K = lambda i,j: j*Nt+i
    A=zeros((Nt*Np,Nt*Np))    # we don't need this now that we went to the sparse matrix case, but keep it here for consistency for now
    from scipy.sparse import coo_matrix

    Jr = zeros(Nt*Np)

    nnz = Np*(Nt-2)*5 + Np + Np*(Np+1) # this is the number of non-zeros in the matrix. Note, this depends on the stencil and boundary conditions
    data=zeros(nnz)
    II = zeros(nnz)
    JJ = zeros(nnz)
    count = 0
    t0=time.clock()
    # inner block
    for j in arange(0,Np):
        for i in arange(1,Nt-1):

            jm1 = (j-1)%Np
            jp1  = (j+1)%Np
            j     = j%Np  # not needed (J<=Np-1) but do it for symmetry

            # note, we use the above definitions of j's except in reference to dp array, which has Np-1 points
            # in those cases, we take the modulo in place

            ft = 1./(dt[i,j]+dt[i-1,j])/sin(grd.t[i,j])
            fp = 1./(dp[i,j%(Np-1)]+dp[i, (j-1)%(Np-1)])/sin(grd.t[i,j])**2

            A[K(i,j),K(i,j)] = -ft*( (F11[i,j]+F11[i+1,j])/dt[i,j]+(F11[i,j]+F11[i-1,j])/dt[i-1,j] ) - \
                                 fp*( (F22[i,j]+F22[i,jp1])/dp[i,j%(Np-1)]+(F22[i,j]+F22[i,jm1])/dp[i, (j-1)%(Np-1)] )
            data[count] = A[K(i,j),K(i,j)]
            II[count] = K(i,j)
            JJ[count] = K(i,j)
            count+=1

            A[K(i,j),K(i+1,j)] = ft*(F11[i,j]+F11[i+1,j])/dt[i,j]
            data[count] = A[K(i,j),K(i+1,j)]
            II[count]  = K(i,j)
            JJ[count] = K(i+1,j)
            count+=1

            A[K(i,j),K(i-1,j)] = ft*(F11[i,j]+F11[i-1,j])/dt[i-1,j]
            data[count]  = A[K(i,j),K(i-1,j)]
            II [count] = K(i,j)
            JJ[count]  = K(i-1,j)
            count+=1

            A[K(i,j),K(i,jp1)] = fp*(F22[i,j]+F22[i,jp1])/dp[i, j%(Np-1)]
            data[count] = A[K(i,j),K(i,jp1)]
            II[count] = K(i,j)
            JJ[count] = K(i,jp1)
            count+=1

            A[K(i,j),K(i,jm1)] = fp*(F22[i,j]+F22[i,jm1])/dp[i, (j-1)%(Np-1)]
            data[count] = A[K(i,j),K(i,jm1)]
            II[count] = K(i,j)
            JJ[count] = K(i,jm1)
            count+=1

            Jr[K(i,j)] = jr[i,j]

    # low lat boundary
    for j in arange(0,Np):
        A[K(Nt-1,j),K(Nt-1,j)] = 1.
        data[count] = A[K(Nt-1,j),K(Nt-1,j)]
        II[count] = K(Nt-1,j)
        JJ[count] = K(Nt-1,j)
        count+=1

        Jr[K(Nt-1,j)] = 0. # FIXME: make arbitrary BC

    # pole boundary
    for j in arange(0,Np):
        A[K(0,j),K(0,j)] = 1.
        data[count] = A[K(0,j),K(0,j)]
        II[count] = K(0,j)
        JJ[count] =  K(0,j)
        count+=1

        for jj in arange(0,Np): 
            A[K(0,j),K(1,jj)] = -dp[1,jj%(Np-1)]/(2.*pi)
            data[count] = A[K(0,j),K(1,jj)]
            II[count] = K(0,j)
            JJ[count] = K(1,jj)
            count+=1
            # note, not setting J because it's initializaed to zero anyway
            
    print('Time spend constructing matrix (s):',time.clock()-t0)
    t0 = time.clock()

    # Now, do the solve
    import scipy
    from scipy.sparse import linalg
    M = coo_matrix( (array(data),(array(II),array(JJ))),shape=(Np*Nt,Np*Nt) )
    #    pot,status=scipy.sparse.linalg.lgmres(A,Jr)    # this becomes very slow for large matrices
    #    solve=scipy.sparse.linalg.factorized(A)   # this is faster but still slow. need to use sparse matrix as implemented now
    #    pot=solve(Jr)

    pot=scipy.sparse.linalg.spsolve(M.tocsc(),Jr)

    t0=time.clock()
    print('Time spend on solve (s):',time.clock()-t0)

    pot*=6.5**2*1.e3

#    print "Convergence status: ", status
    print "Potential min/max (kV): ",pot.min(),pot.max()

    # Plotting

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,50.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,pot.reshape(Np,Nt).T,51)
    colorbar().set_label('Potential, kV')

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,50.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,jr,45)
    colorbar().set_label('Current, microA/m^2')
    contour(grd.p+pi+pi/2,grd.t,pot.reshape(Np,Nt).T,13,colors='black')
    

