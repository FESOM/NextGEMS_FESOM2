import xarray as xr
import pandas as pd
import os
import numpy as np
from collections import namedtuple
import math as mt
from joblib import Parallel, delayed
import sys
from unrotate_fesom2_UV_utils import *

try:
    mesh_path = sys.argv[1]
    ufile = sys.argv[2]
    vfile = sys.argv[3]
    ofolder = sys.argv[4]
    n_jobs = sys.argv[5]
except:
    raise ValueError("provide proper input parameters")
    
# mesh_path = '/work/ab0995/a270046/fesom2-meshes/NEMO_ecmwf/NEMO/'
# ufile = '/work/bk1040/DYAMOND/.input_winter_data/IFS-FESOM2-4km/u.fesom.2020.nc'
# vfile = '/work/bk1040/DYAMOND/.input_winter_data/IFS-FESOM2-4km/v.fesom.2020.nc'
# ofolder = '/mnt/lustre01/work/ab0995/a270088/DYAMOND/ROTATE_UV/IFS-FESOM2-4km/'


nodes = pd.read_csv(os.path.join(mesh_path,'nod2d.out'), skiprows=1, delim_whitespace=True, names=['n','x','y','f'])
elem = np.loadtxt(os.path.join(mesh_path,'elem2d.out'), skiprows=1)
elem = (elem-1).astype('int')

fmesh = namedtuple("mesh", "x2 y2 elem")
mesh = fmesh(x2=nodes.x.values.astype('float32'), y2=nodes.y.values.astype('float32'), elem=elem)

face_x, face_y = compute_face_coords(mesh)

u = xr.open_dataset(ufile)
v = xr.open_dataset(vfile)

def fesoro(time):
    
    u3d = np.zeros((1,u.u.shape[1], u.u.shape[2])).astype('float32')
    v3d = np.zeros((1,u.u.shape[1], u.u.shape[2])).astype('float32')
    
    u0 = u.u[time,:,:].data
    v0 = v.v[time,:,:].data
    
    for i in range(u.u.shape[2]):
        ut = u0[:,i]
        vt = v0[:,i]
        urot, vrot = vec_rotate_r2g(50, 15, -90, face_x, face_y, ut, vt, flag=1)
        u3d[0,:,i] = urot
        v3d[0,:,i] = vrot
#         print(i)
    out_u = xr.Dataset({'u':(['time','elem','nz1'], u3d.astype('float32'), u.u.attrs)},
          coords={'time':[u.time[time].data],
                 }, attrs=u.attrs)
    out_v = xr.Dataset({'v':(['time','elem','nz1'], v3d.astype('float32'), v.v.attrs)},
          coords={'time':[v.time[time].data],
                 }, attrs=v.attrs)
    
    out_u.to_netcdf(f'{ofolder}/u_{str(time).zfill(10)}.nc', encoding={'u': {'dtype': 'float32', '_FillValue': None}})
    out_v.to_netcdf(f'{ofolder}/v_{str(time).zfill(10)}.nc', encoding={'v': {'dtype': 'float32', '_FillValue': None}})
    print(f'time processed {str(time).zfill(10)}')
    

Parallel(n_jobs=int(n_jobs))(delayed(fesoro)(i) for i in range(u.u.shape[0]))

uu = xr.open_mfdataset(f'{ofolder}/u_*.nc', chunks={'time':30})
uu.to_netcdf(f'{ofolder}/fesom_unrotated_uout.nc')
uu.close()

vv = xr.open_mfdataset(f'{ofolder}/v_*.nc', chunks={'time':30})
vv.to_netcdf(f'{ofolder}/fesom_unrotated_vout.nc')
vv.close()



    
    
    
    
    
    