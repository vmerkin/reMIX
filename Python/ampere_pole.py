from numpy import loadtxt,pi,sin,concatenate,zeros,hstack,ones

class ampere_pole():
    def __init__(self,file_name):
        self.fileName = file_name
        self.data=loadtxt(self.fileName)

    def get_grid(self):
        p = self.data[:24,:]*pi/12.
        t = self.data[24:48,:]*pi/180.

        t = hstack( ( zeros((t.shape[0],1)),t ) )
        p = hstack( ( p[:,[0]],p ) )
        
        # FIXME: hardwiring uniform spherical grid for now
        # need to generalized for arbirary grids
        dt = t[0,1]-t[0,0]
        dp = p[1,0]-p[0,0]
        A = sin(t)*dt*dp

        # note, always have theta first, phi second
        # return area as well
        return(t.T,p.T,A.T)  


    def get_data(self):
        tmp = self.data[48:,:]
        tmp = hstack( ( ones((tmp.shape[0],1))*tmp[:,0].mean(),tmp ) )
        return tmp.T 

