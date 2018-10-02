# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

import random
import csv
import time
from sage.misc.sage_timeit import sage_timeit

default_function_name_list=['gj_forward_3_slope','drlm_backward_3_slope','chen_4_slope','kzh_7_slope_1','kzh_10_slope_1']
default_two_slope_fill_in_epsilon_list=[1/(i*10) for i in range(1,11)]
default_perturbation_epsilon_list=[i/100 for i in range(3)]

def write_performance_table(function_name_list,two_slope_fill_in_epsilon_list,perturbation_epsilon_list):
    with open('performance.csv', mode='w') as file:
        performance_table = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        performance_table.writerow(['base_function', 'two_slope_fill_in_epsilon', 'perturbation_epsilon', '# breakpoints','# slopes','# vertices','# additive vertices','add_v/v ratio','delta pi min','is subadditive','# .min nodes (constant)', '# .min nodes (affine)', '.min n/v ratio (constant)','.min n/v ratio (affine)', '# .is_subadditive nodes (constant)', '# .is_subadditive nodes (affine)', '.is_subadditive n/v ratio (constant)', '.is_subadditive n/v ratio (affine)', 'time .min (constant) (s)','time .min (affine) (s)', 'time .is_subadditive (constant) (s)','time .is_subadditive (affine) (s)', 'time subadditivity_test (s)'])
        for name in function_name_list:
            global base_fn
            base_fn=eval(name)()
            v=number_of_vertices(base_fn)
            add_v=number_of_additive_vertices(base_fn)
            m,n_m_c,t_m_c=measure_T_min(base_fn,bound="constant")
            m,n_m_a,t_m_a=measure_T_min(base_fn,bound="affine")
            is_sub,n_i_c,t_i_c=measure_T_is_subadditive(base_fn,bound="constant")
            is_sub,n_i_a,t_i_a=measure_T_is_subadditive(base_fn,bound="affine")
            performance_table.writerow([name,None,None,len(base_fn.end_points()),number_of_slopes(base_fn),v,add_v,float(add_v/v),m,is_sub,n_m_c,n_m_a,float(n_m_c/v),float(n_m_a/v),n_i_c,n_i_a,float(n_i_c/v),float(n_i_a/v),t_m_c,t_m_a,t_i_c,t_i_a,sage_timeit('subadditivity_test(base_fn)',globals(),seconds=True)])
            for fill_in_epsilon in two_slope_fill_in_epsilon_list:
                for perturb_epsilon in perturbation_epsilon_list:
                    global fn
                    fn=test_function_from_two_slope_fill_in_extreme_functions(base_fn,fill_in_epsilon,perturb_epsilon)
                    add_v=number_of_additive_vertices(fn)
                    v=number_of_vertices(fn)
                    m,n_m_c,t_m_c=measure_T_min(fn,bound="constant")
                    m,n_m_a,t_m_a=measure_T_min(fn,bound="affine")
                    is_sub,n_i_c,t_i_c=measure_T_is_subadditive(fn,bound="constant")
                    is_sub,n_i_a,t_i_a=measure_T_is_subadditive(fn,bound="affine")
                    performance_table.writerow([name,float(fill_in_epsilon),float(perturb_epsilon),len(fn.end_points()),number_of_slopes(fn),v,add_v,float(add_v/v),m,is_sub,n_m_c,n_m_a,float(n_m_c/v),float(n_m_a/v),n_i_c,n_i_a,float(n_i_c/v),float(n_i_a/v),t_m_c,t_m_a,t_i_c,t_i_a,sage_timeit('subadditivity_test(fn)',globals(),seconds=True)])

def measure_T_min(fn,bound="constant"):
    global T
    T=SubadditivityTestTree(fn)
    if bound=="constant":
        t=sage_timeit('T.minimum(max_number_of_bkpts=0)',globals(),seconds=True)
        return T.min,T.number_of_nodes(),t
    else:
        t=sage_timeit('T.minimum(max_number_of_bkpts=100000)',globals(),seconds=True)
        return T.min,T.number_of_nodes(),t

def measure_T_is_subadditive(fn,bound="constant"):
    global T
    T=SubadditivityTestTree(fn)
    if bound=="constant":
        t=sage_timeit('T.is_subadditive(max_number_of_bkpts=0)',globals(),seconds=True)
        return T._is_subadditive,T.number_of_nodes(),t
    else:
        t=sage_timeit('T.is_subadditive(max_number_of_bkpts=100000)',globals(),seconds=True)
        return T._is_subadditive,T.number_of_nodes(),t

def test_function_from_two_slope_fill_in_extreme_functions(fn,fill_in_epsilon=1/10,perturb_epsilon=1/10,seed=1,N=10,stay_in_unit_interval=False):
    """
    Return a perturbation function of the 2 slope approximation function of the given function fn.
    """
    f1=two_slope_fill_in_extreme(fn,fill_in_epsilon)
    return function_random_perturbation(f1,epsilon=perturb_epsilon,seed=seed,N=N,stay_in_unit_interval=stay_in_unit_interval)

def function_random_perturbation(fn,epsilon=1/10,seed=1,N=10,stay_in_unit_interval=False):
    """
    Return a random perturbation of the given function fn.
    """
    random.seed(seed)
    bkpts=[]
    values=[]
    if epsilon==0:
        return fn
    for bkpt in fn.end_points():
        bkpts.append(bkpt)
        pert=random.randint(-N,N)
        value=fn(bkpt)+pert*epsilon/N
        if stay_in_unit_interval:
            if value>1:
                value=1
            if value<0:
                value=0
        values.append(value)
    return piecewise_function_from_breakpoints_and_values(bkpts, values)
    
def histogram_delta_pi(fn,sampling='vertices',q=5,epsilon=1/10000):
    """
    The histogram of the values of delta pi over given points in the complex.
    """
    values=[]
    if sampling=='vertices':
        bkpts=fn.end_points()
        bkpts2=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
        for x in bkpts:
            for y in bkpts:
                values.append(delta_pi(fn,x,y))
        for z in bkpts2[1:-1]:
            for x in bkpts[1:-1]:
                y=z-x
                if 0<y<1 and y not in bkpts:
                    val=delta_pi(fn,x,y)
                    #symmetry
                    values=values+[val,val]
    elif sampling=='grid':
        for i in range(q+1):
            for j in range(q+1):
                x=i/q
                y=j/q
                values.append(delta_pi(fn,x,y))
    else:
        raise ValueError, "Can't recognize sampling."
    return np.histogram(values,bins=[0,epsilon,1/2,1,3/2,2], density=False)

def number_of_additive_vertices(fn):
    counter=0
    bkpts=fn.end_points()
    bkpts2=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    for x in bkpts:
        for y in bkpts:
            if delta_pi(fn,x,y)==0:
                counter+=1
    for z in bkpts2[1:-1]:
        for x in bkpts[1:-1]:
            y=z-x
            if 0<y<1 and y not in bkpts and delta_pi(fn,x,y)==0:
                counter+=2
    return counter

def number_of_vertices(fn):
    """
    Return the number of vertices of the complex delta_pi.
    """
    bkpts=fn.end_points()
    bkpts2=fn.end_points()[1:-1]+[1+bkpt for bkpt in fn.end_points()[:-1]]
    counter=len(bkpts)^2
    for z in bkpts2:
        for x in bkpts:
            y=z-x
            if 0<y<1 and y not in bkpts:
                #symmetry
                counter+=2
    return counter

def additive_vertices_ratio(fn):
    return number_of_additive_vertices(fn)/number_of_vertices(fn)


    
