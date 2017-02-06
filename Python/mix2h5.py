import mix
import h5py
import os.path
from numpy import array

filein = '/Users/merkivg1/work/Events/August2010/data/Aug2010_mix_2010-08-04T00-00-00Z.hdf'
fileout = filein[:-3]+'h5'

if __name__ == "__main__":
    if os.path.exists(fileout): os.remove(fileout)

    ########################
    # read MIX data in
    ########################
    (simtime,
     x,y,
     psi_n,psi_s,
     fac_n,fac_s,
     sigmap_n,sigmap_s,
     sigmah_n,sigmah_s,
     energy_n,energy_s,
     flux_n,flux_s) = mix.get_data(filein)

    f=h5py.File(fileout,'w')

    # note, cutting out periodic boundary
    f.attrs.create('SimTime',array(simtime),dtype=int)
    f.create_dataset('Grid X',data=x[:-1,:].T,dtype='float32')
    f.create_dataset('Grid Y',data=y[:-1,:].T,dtype='float32')
    f.create_dataset('Potential North [kV]',data=psi_n[:-1,:].T,dtype='float32')
    f.create_dataset('Potential South [kV]',data=psi_s[:-1,:].T,dtype='float32')
    f.create_dataset('FAC North [microAm2]',data=-fac_n[:-1,:].T,dtype='float32')
    f.create_dataset('FAC South [microAm2]',data=fac_s[:-1,:].T,dtype='float32')
    f.create_dataset('Pedersen conductance North [S]',data=sigmap_n[:-1,:].T,dtype='float32')
    f.create_dataset('Pedersen conductance South [S]',data=sigmap_s[:-1,:].T,dtype='float32')
    f.create_dataset('Hall conductance North [S]',data=sigmah_n[:-1,:].T,dtype='float32')
    f.create_dataset('Hall conductance South [S]',data=sigmah_s[:-1,:].T,dtype='float32')
    f.create_dataset('Average energy North [keV]',data=energy_n[:-1,:].T,dtype='float32')
    f.create_dataset('Average energy South [keV]',data=energy_s[:-1,:].T,dtype='float32')
    f.create_dataset('Number flux North [1cm2 s]',data=flux_n[:-1,:].T,dtype='float32')
    f.create_dataset('Number flux South [1cm2 s]',data=flux_s[:-1,:].T,dtype='float32')

    f.close()
