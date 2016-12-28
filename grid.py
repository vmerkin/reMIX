from numpy import pi,linspace,meshgrid

class grid():
    def __init__(self):
        self.t = None
        self.p = None
        self.nt = None
        self.np = None
        self.A = None    # area elements

class uniform(grid):
    def __init__(self,nt,np):
        # note, these are cell centers
        self.t,self.p = meshgrid(linspace(0,pi/4,nt+2)[1:nt+1],linspace(0,2*pi,np+2)[1:np+1],indexing='ij')
        self.nt = nt
        self.np = np

class generate_from_data(grid):
    def __init__(self,data):
        self.t = data[0]
        self.p = data[1]
        self.A = data[2]
        self.nt = data[1].shape[0]
        self.np = data[1].shape[1]

