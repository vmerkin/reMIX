from numpy import loadtxt,pi

class ampere():
    def __init__(self,file_name):
        self.fileName = file_name
        self.data=loadtxt(self.fileName)

    def get_grid(self):
        p = self.data[:24,:]*pi/12.
        t = self.data[24:48,:]*pi/180.
        return(t.T,p.T)  # note, always have theta first, phi second

    def get_data(self):
        return self.data[48:,:].T
