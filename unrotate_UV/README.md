Unrotate FESOM2 data
====================

FESOM2 vector data outputed in rotated coordinates, that is internally used for computation. We have to unrotate them back to geographical ones to be able to plot.

This example does not use [pyfesom2](https://github.com/FESOM/pyfesom2/tree/master/pyfesom2), although with it code would be shorter.

Files
-----

- **unrotate_fesom2_UV.py** - main script.
- **unrotate_fesom2_UV_utils.py** - a small module with necessary functions. If you copy the main script, this additional file should come with it.
- **Rotate_uv_testing.ipynb** - notebook where the script was developed, can use for interactive work

Usage
-----
on mistral you need

```
module load python3/unstable
```

On other machines dependencies listed in the beggining of the file should be installed (the only "non standard" one is joblib).

```bash
python unrotate_fesom2_UV.py path_to_fesom2_mesh_folder path_to_U_file path_to_V_file path_to_output_dir number_of_processes
```
Example:

```
python unrotate_fesom2_UV.py /work/ab0995/a270046/fesom2-meshes/NEMO_ecmwf/NEMO/ /work/bk1040/DYAMOND/.input_winter_data/IFS-FESOM2-4km/u.fesom.2020.nc /work/bk1040/DYAMOND/.input_winter_data/IFS-FESOM2-4km/v.fesom.2020.nc /mnt/lustre01/work/ab0995/a270088/DYAMOND/ROTATE_UV/IFS-FESOM2-4km/test/ 5
```

You will get in the output folder 2 files **fesom_unrotated_uout.nc** and **fesom_unrotated_vout.nc**, and a bunch of files with individual time steps. You should take of removing them by yourself.

The more resources (cores/memory) you have, the more processes you can use. 

Performance
-----------

There is a lot of room for improvement, but on the 320 timesteps of ORCA25 mesh with 20 processes on instance with 24 cores and 65 Gb of memory conversion takes 16 minutes, plus 5 minutes for saving the output. Currently this is OK. One should look into vectorisation of rotation functions if we want to make it effective, this developent is planned as part of pyfesom2.
