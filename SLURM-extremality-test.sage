#! /home/wangjw/sage/sage-8.3/sage
#SBATCH --array=0-9 --time=00:10:00
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

functions = [
chen_4_slope,
dg_2_step_mir,
dg_2_step_mir_limit,
drlm_2_slope_limit,
drlm_3_slope_limit,
drlm_backward_3_slope,
gj_2_slope,
gj_2_slope_repeat,
gj_forward_3_slope,
gmic]


factory = functions[task_id]
print("Extremality test for {}".format(factory))
h = factory()
print(extremality_test(h))
for epsilon in [1/10,1/100]:
    new_h=two_slope_fill_in_extreme(h,epsilon)
    print(extremality_test(new_h))

