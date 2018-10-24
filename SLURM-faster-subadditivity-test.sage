#! /home/wangjw/sage/sage-8.3/sage
#SBATCH --array=0-11 --time 72:00:00
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

logging.disable(logging.INFO)

functions = [bcdsp_arbitrary_slope,
chen_4_slope,
dg_2_step_mir,
drlm_backward_3_slope,
gj_forward_3_slope,
kzh_7_slope_1,
kzh_7_slope_2,
kzh_7_slope_3,
kzh_7_slope_4,
kzh_10_slope_1,
kzh_28_slope_1,
kzh_28_slope_2]



factory = functions[task_id]
h = factory()
previous_f=None

with open('Coin_performance_%s.csv' % task_id, mode='w') as file:
    performance_table = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    performance_table.writerow(['base_function', 'two_slope_fill_in_epsilon', 'perturbation_epsilon', '# breakpoints','# slopes','# vertices','# additive vertices','add_v/v ratio','delta pi min','is subadditive','.min(constant/BFS)', '.min(lp<=20/BFS)','.min(lp<=50/BFS)','.min(lp<=100/BFS)','.min(lp<=500/BFS)','.min(lp<=10000/BFS)','.min(constant/DFS)', '.min(lp<=20/DFS)','.min(lp<=50/DFS)','.min(lp<=100/DFS)','.min(lp<=500/DFS)','.min(lp<=10000/DFS)','naive subadditivity test (s)'])
    for two_slope_fill_in_epsilon in [0,QQ(1)/10,QQ(1)/20,QQ(1)/50,QQ(1)/100, QQ(1)/1000]:
        if two_slope_fill_in_epsilon==0:
            new_f=h
        else:
            new_f=symmetric_2_slope_fill_in(h,two_slope_fill_in_epsilon)
        if not previous_f:
            previous_f=h
        else:
            if new_f==previous_f:
                continue
            previous_f=new_f
        for perturbation_epsilon in [0,QQ(1)/1000,QQ(1)/500,QQ(1)/100,QQ(1)/50]:
            global fn
            fn=function_random_perturbation(new_f,perturbation_epsilon)
            v=number_of_vertices(fn)
            add_v=number_of_additive_vertices(fn)
            delta_pi_min=minimum_of_delta_pi(fn)
            if delta_pi_min>=0:
                is_sub=True
            else:
                is_sub=False
            res_0_D=measure_T_min(fn,0,'DFS',solver='Coin',number=1,repeat=1)
            res_0_D.append(float(res_0_D[1]/v))
            res_20_D=measure_T_min(fn,20,'DFS',solver='Coin',number=1,repeat=1)
            res_20_D.append(float(res_20_D[1]/v))
            res_50_D=measure_T_min(fn,50,'DFS',solver='Coin',number=1,repeat=1)
            res_50_D.append(float(res_50_D[1]/v))
            res_100_D=measure_T_min(fn,100,'DFS',solver='Coin',number=1,repeat=1)
            res_100_D.append(float(res_100_D[1]/v))
            res_500_D=measure_T_min(fn,500,'DFS',solver='Coin',number=1,repeat=1)
            res_500_D.append(float(res_500_D[1]/v))
            res_10000_D=measure_T_min(fn,10000,'DFS',solver='Coin',number=1,repeat=1)
            res_10000_D.append(float(res_10000_D[1]/v))
            res_0_B=measure_T_min(fn,0,'BFS',solver='Coin',number=1,repeat=1)
            res_0_B.append(float(res_0_B[1]/v))
            res_20_B=measure_T_min(fn,20,'BFS',solver='Coin',number=1,repeat=1)
            res_20_B.append(float(res_20_B[1]/v))
            res_50_B=measure_T_min(fn,50,'BFS',solver='Coin',number=1,repeat=1)
            res_50_B.append(float(res_50_B[1]/v))
            res_100_B=measure_T_min(fn,100,'BFS',solver='Coin',number=1,repeat=1)
            res_100_B.append(float(res_100_B[1]/v))
            res_500_B=measure_T_min(fn,500,'BFS',solver='Coin',number=1,repeat=1)
            res_500_B.append(float(res_500_B[1]/v))
            res_10000_B=measure_T_min(fn,10000,'BFS',solver='Coin',number=1,repeat=1)
            res_10000_B.append(float(res_10000_B[1]/v))
            naive_time=sage_timeit('minimum_of_delta_pi(fn)',globals(),seconds=True,number=1,repeat=1)
            performance_table.writerow([factory,two_slope_fill_in_epsilon,perturbation_epsilon,len(fn.end_points()),number_of_slopes(fn),v,add_v,float(add_v/v),delta_pi_min,is_sub,res_0_B,res_20_B,res_50_B,res_100_B,res_500_B,res_10000_B,res_0_D,res_20_D,res_50_D,res_100_D,res_500_D,res_10000_D,naive_time])



