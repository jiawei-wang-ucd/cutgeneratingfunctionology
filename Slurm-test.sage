#! /home/wangjw/sage/sage-8.3/sage
#SBATCH --array=0-1 --time=01:00:00
#
# --array: Specify the range of the array tasks.
# --time: Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds" etc.
import os
task_id = os.getenv("SLURM_ARRAY_TASK_ID")
if task_id:
    task_id = int(task_id)

print("Sage hello from task {}".format(task_id))

import sys
sys.path = [''] + sys.path

import igp
from igp import *

factory=1

for i in range (1000):
    for j in range(1000):
        for k in range(1000):
            for l in range(1000):
                factory*=1
print("results : {}".format(factory))


