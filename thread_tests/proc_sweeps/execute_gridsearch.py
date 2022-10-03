from re import sub
import subprocess

for num_proc in range(1, 7):
    subprocess.run(['mpirun', 
                    '-n', 
                    str(num_proc), 
                    'python', 
                    '../Scripts/mpiMapSpecAngleScan.py', 
                    str(12)])