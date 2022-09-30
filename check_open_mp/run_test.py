import os
import subprocess
os.environ["OMP_NUM_THREADS"] = "5"

subprocess.run('./test_openmp')
