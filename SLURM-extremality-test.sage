#! /home/mkoeppe/sage/sage
#SBATCH --array=1-14 --time 10
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

functions = [bcdsp_arbitrary_slope,
bhk_gmi_irrational,
bhk_irrational,
bhk_slant_irrational,
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
extremality_test(h)

