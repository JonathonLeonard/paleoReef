import os
import subprocess
import numpy as np
import pandas as pd

from mpi4py import MPI

times = pd.read_csv('data/bristol_sim_list.csv')['time'].values

MPIcomm = MPI.COMM_WORLD
MPIrank = MPIcomm.Get_rank()
MPIsize = MPIcomm.Get_size()

splits = np.array_split(times, MPIsize)
steps = splits[MPIrank]

for k in range(0,len(steps)):
    folder = 'photozoan.res2/proj_'+str(steps[k])+'Ma'
    if not os.path.isdir(folder):
        subprocess.call(['Rscript', 'projectModel.r', str(steps[k]), 'res2'])