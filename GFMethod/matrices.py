from numpy import cos, sin,zeros,arange,pi
import point,sphtrig,util

class matrices():

    def __init__(self,grid_s,grid_d):
        """ Initiate the grids. Then compute the M matrices, since they only depend on the grids.
        
        Arguments
        grid_s -- the source grid.
        grid_d -- the destination grid.
        """

        self.grid_s = grid_s
        self.grid_d = grid_d

        (self.M1,self.M2,self.dM1,self.dM2) = self.M()

    def M(self):
        M1 = zeros((2*self.grid_d.nt*self.grid_d.np,self.grid_s.nt*self.grid_s.np),order='F')
        M2 = zeros((2*self.grid_d.nt*self.grid_d.np,self.grid_s.nt*self.grid_s.np),order='F')
        dM1 = zeros((2*self.grid_d.nt*self.grid_d.np,self.grid_s.nt*self.grid_s.np),order='F')
        dM2 = zeros((2*self.grid_d.nt*self.grid_d.np,self.grid_s.nt*self.grid_s.np),order='F')

        for iv in arange(self.grid_d.nt):
            for jv in arange(self.grid_d.np):
                v = iv*self.grid_d.np+jv

                for iu in arange(self.grid_s.nt):
                    for ju in arange(self.grid_s.np):
                        u = iu*self.grid_s.np+ju

                        pd = point.point(self.grid_d.t[iv,jv],self.grid_d.p[iv,jv])
                        ps = point.point(self.grid_s.t[iu,ju],self.grid_s.p[iu,ju])
                     
                        tp = sphtrig.theta_p(ps,pd)
                        dA = self.grid_s.A[iu,ju]
                        M1[2*v,u] = dA/4./pi*(cos(ps.t)-cos(pd.t)*cos(tp))/sin(pd.t)/(1.-cos(tp))
                        M1[2*v+1,u]  = -dA/4./pi*sin(ps.t)*sin(ps.p-pd.p)/(1.-cos(tp))

                        dM1[2*v,u] = -dA/4./pi* sin(ps.t)*sin(ps.p-pd.p)/(1-cos(tp)) *  ( 
                            (cos(ps.t)-cos(pd.t)*cos(tp))/(1-cos(tp)) + cos(pd.t)  
                        )

                        dM1[2*v+1,u] = -dA/4./pi*sin(ps.t)/(1-cos(tp)) * (
                            -sin(ps.p-pd.p)**2/(1-cos(tp))*sin(pd.t)*sin(ps.t) - cos(ps.p-pd.p)
                        )
                        
                        M2[2*v,u] = -M1[2*v+1,u]
                        M2[2*v+1,u]    =  M1[2*v,u]
                        dM2[2*v,u] = -dM1[2*v+1,u]
                        dM2[2*v+1,u]    =  dM1[2*v,u]

        return(M1,M2,dM1,dM2)

    def L(self,Sigma,jr):
        L1 = zeros((self.grid_d.nt*self.grid_d.np,self.grid_s.nt*self.grid_s.np))
        L2 = zeros((self.grid_d.nt*self.grid_d.np,self.grid_s.nt*self.grid_s.np))
        divJ = zeros(self.grid_s.nt*self.grid_s.np)
#################
        # Move inside the loop once conductances are defined
        cosd = 1.  # fix this
        Sigma_P = 8. # fix 
        Sigma_H = 0. # fix 
        Rm = -Sigma_P/(Sigma_P**2+Sigma_H**2)*cosd
#################

#        dM1dphi = util.dPhi(Rm*self.M1)   # DIVIDE by dPhi eventually
#        dM2dphi = util.dPhi(Rm*self.M2)   # DIVIDE by dPhi eventually
#        dA = 0.0019350268407194113
        #(self.grid_s.t[1]-self.grid_s.t[0])*(self.grid_s.p[1]-self.grid_s.p[0])

        for iv in arange(self.grid_d.nt):
            for jv in arange(self.grid_d.np):
                v = iv*self.grid_d.np+jv

                for iu in arange(self.grid_s.nt):
                    for ju in arange(self.grid_s.np):
                        u = iu*self.grid_s.np+ju

                        if iu==self.grid_s.nt-1:
                            dt = self.grid_s.t[iu,ju]-self.grid_s.t[iu-1,ju]
                        else:
                            dt = self.grid_s.t[iu+1,ju]-self.grid_s.t[iu,ju]

                        if ju==self.grid_s.np-1:
                            dp = self.grid_s.p[iu,ju]-self.grid_s.p[iu,ju-1]
                        else:
                            dp = self.grid_s.p[iu,ju+1]-self.grid_s.p[iu,ju]


                        divJ[u] = jr[iu,ju]
                        pd = point.point(self.grid_d.t[iv,jv],self.grid_d.p[iv,jv])
                        ps = point.point(self.grid_s.t[iu,ju],self.grid_s.p[iu,ju])

                        L1[v,u] = -Rm*self.dM1[2*v,u]/sin(pd.t)
                        L2[v,u] = -Rm*self.dM2[2*v,u]/sin(pd.t) + Sigma_P*(
                            ( 
                                (abs( self.grid_s.t[iu,ju] -  pd.t ) <=dt/2.) & (abs( self.grid_s.p[iu,ju] -  pd.p) <=dp/2.)
                            )*1. - sin(pd.t)*self.grid_s.A[iu,ju]/4./pi)
        return(L1,L2,divJ)




