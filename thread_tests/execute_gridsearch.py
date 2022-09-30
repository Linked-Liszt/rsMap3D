from re import sub
import subprocess

thread_search_axis = [48, 24, 12, 6, 2, 1]
proc_search_axis = [6, 3, 2, 1]

for num_proc in proc_search_axis:
    for num_thread in thread_search_axis:
        subprocess.run(['mpirun', 
                        '-n', 
                        str(num_proc), 
                        'python', 
                        '../Scripts/mpiMapSpecAngleScan.py', 
                        str(num_thread)])