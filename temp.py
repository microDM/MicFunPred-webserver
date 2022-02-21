import subprocess

cmd = "McFunPred_run_pipeline.py -h"
subBlast = subprocess.Popen(cmd,stdin=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
try:
    output,error = subBlast.communicate()
except subprocess.SubprocessError:
    print('dsds')