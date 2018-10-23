# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

import random
import csv
import time
from sage.misc.sage_timeit import sage_timeit

random.seed(500)

default_function_name_list=['gj_forward_3_slope','drlm_backward_3_slope','chen_4_slope','kzh_7_slope_1','kzh_10_slope_1']
default_two_slope_fill_in_epsilon_list=[1/(i*10) for i in range(1,11)]
default_perturbation_epsilon_list=[i/100 for i in range(3)]
default_max_number_of_bkpts=[0,10,20,40,100,400,1000,10000,100000]


def write_performance_table(function_name_list,two_slope_fill_in_epsilon_list,perturbation_epsilon_list):
    with open('performance.csv', mode='w') as file:
        performance_table = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        performance_table.writerow(['base_function', 'two_slope_fill_in_epsilon', 'perturbation_epsilon', '# breakpoints','# slopes','# vertices','# additive vertices','add_v/v ratio','delta pi min','is subadditive','# .min nodes (constant)', '# .min nodes (affine)', '# .min nodes (mixed)','.min n/v ratio (constant)','.min n/v ratio (affine)','.min n/v ratio (mixed)', '# .is_subadditive nodes (constant)', '# .is_subadditive nodes (affine)', '# .is_subadditive nodes (mixed)', '.is_subadditive n/v ratio (constant)', '.is_subadditive n/v ratio (affine)', '.is_subadditive n/v ratio (mixed)','time .min (constant) (s)','time .min (affine) (s)','time .min (mixed) (s)', 'time .is_subadditive (constant) (s)','time .is_subadditive (affine) (s)', 'time .is_subadditive (mixed) (s)','time subadditivity_test (s)'])
        for name in function_name_list:
            global base_fn
            base_fn=eval(name)()
            v=number_of_vertices(base_fn)
            add_v=number_of_additive_vertices(base_fn)
            m,n_m_c,t_m_c=measure_T_min(base_fn,bound="constant",number=2,repeat=2)
            m,n_m_a,t_m_a=measure_T_min(base_fn,bound="affine",number=2,repeat=2)
            m,n_m_m,t_m_m=measure_T_min(base_fn,bound="mixed",number=2,repeat=2)
            is_sub,n_i_c,t_i_c=measure_T_is_subadditive(base_fn,bound="constant",number=2,repeat=2)
            is_sub,n_i_a,t_i_a=measure_T_is_subadditive(base_fn,bound="affine",number=2,repeat=2)
            is_sub,n_i_m,t_i_m=measure_T_is_subadditive(base_fn,bound="mixed",number=2,repeat=2)
            performance_table.writerow([name,None,None,len(base_fn.end_points()),number_of_slopes(base_fn),v,add_v,float(add_v/v),m,is_sub,n_m_c,n_m_a,n_m_m, float(n_m_c/v),float(n_m_a/v),float(n_m_m/v),n_i_c,n_i_a,n_i_m,float(n_i_c/v),float(n_i_a/v),float(n_i_m/v),t_m_c,t_m_a,t_m_m,t_i_c,t_i_a,t_i_m,sage_timeit('subadditivity_test(base_fn,stop_if_fail=True)',globals(),number=2,repeat=2,seconds=True)])
            for fill_in_epsilon in two_slope_fill_in_epsilon_list:
                for perturb_epsilon in perturbation_epsilon_list:
                    global fn
                    fn=test_function_from_two_slope_fill_in_extreme_functions(base_fn,fill_in_epsilon,perturb_epsilon)
                    add_v=number_of_additive_vertices(fn)
                    v=number_of_vertices(fn)
                    m,n_m_c,t_m_c=measure_T_min(fn,bound="constant",number=1,repeat=1)
                    m,n_m_a,t_m_a=measure_T_min(fn,bound="affine",number=1,repeat=1)
                    is_sub,n_i_c,t_i_c=measure_T_is_subadditive(fn,bound="constant",number=1,repeat=1)
                    is_sub,n_i_a,t_i_a=measure_T_is_subadditive(fn,bound="affine",number=1,repeat=1)
                    performance_table.writerow([name,None,None,len(fn.end_points()),number_of_slopes(fn),v,add_v,float(add_v/v),m,is_sub,n_m_c,n_m_a,n_m_m, float(n_m_c/v),float(n_m_a/v),float(n_m_m/v),n_i_c,n_i_a,n_i_m,float(n_i_c/v),float(n_i_a/v),float(n_i_m/v),t_m_c,t_m_a,t_m_m,t_i_c,t_i_a,t_i_m,sage_timeit('subadditivity_test(fn,stop_if_fail=True)',globals(),number=1,repeat=1,seconds=True)])

def measure_T_min(fn,max_number_of_bkpts,search_method,**kwds):
    global f
    f=fn
    t2=sage_timeit('T=SubadditivityTestTree(f)',globals(),seconds=True)
    def time_min(max_number_of_bkpts=max_number_of_bkpts,search_method=search_method,**kwds):
        global T
        T=SubadditivityTestTree(f)
        T.minimum(max_number_of_bkpts=max_number_of_bkpts,search_method=search_method,**kwds)
    global proc
    proc = time_min
    t1=sage_timeit('proc()',globals(),seconds=True,**kwds)
    return [T.min,T.number_of_nodes(),t1-t2]

def measure_T_is_subadditive(fn,max_number_of_bkpts,search_method,**kwds):
    global f
    f=fn
    t2=sage_timeit('T=SubadditivityTestTree(f)',globals(),seconds=True)
    def time_limit(max_number_of_bkpts=max_number_of_bkpts,search_method=search_method,**kwds):
        global T
        T=SubadditivityTestTree(f)
        T.is_subadditive(stop_if_fail=True, max_number_of_bkpts=max_number_of_bkpts,search_method=search_method,**kwds)
    global proc
    proc = time_limit
    t1=sage_timeit('proc()',globals(),seconds=True,**kwds)
    return [T.number_of_nodes(),t1-t2]

def function_random_perturbation(fn,epsilon,number_of_bkpts_ratio=10):
    """
    Return a random perturbation of the given function fn. Randomly perturb function values at randomly chosen 1/10 breakpoints.
    """
    values=fn.values_at_end_points()
    n=len(fn.end_points())
    if epsilon==0:
        return fn
    bkpts_number=n//number_of_bkpts_ratio
    pert_bkpts=random.sample(range(1, n-1), bkpts_number)
    for i in pert_bkpts:
        if random.randint(0,1)==0:
            values[i]+=epsilon
        else:
            values[i]-=epsilon
    return piecewise_function_from_breakpoints_and_values(fn.end_points(), values)
    
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

def minimum_of_delta_pi(fn):
    """
    Return the min of delta_pi of fn. (Quatratic complexity)
    """
    global_min=10000
    for x in fn.end_points():
        for y in fn.end_points():
            delta=delta_pi(fn,x,y)
            if delta<global_min:
                global_min=delta
    for z in fn.end_points():
        for x in fn.end_points():
            y=z-x
            delta=delta_pi(fn,x,y)
            if delta<global_min:
                global_min=delta
    for z in fn.end_points():
        for x in fn.end_points():
            z=1+z
            y=z-x
            delta=delta_pi(fn,x,y)
            if delta<global_min:
                global_min=delta
    return global_min

def is_goal_reached(fn,goal=0,stop_if_fail=True,keep_exact_solutions=True):
    """
    Return if delta_pi of fn can reach goal-epsilon. (Quatratic complexity)
    """
    exact_solutions=set()
    superior_solutions=set()
    for x in fn.end_points():
        for y in fn.end_points():
            delta=delta_pi(fn,x,y)
            if keep_exact_solutions and delta==goal:
                exact_solutions.add((x,y))
            if delta<goal:
                superior_solutions.add((x,y))
                if stop_if_fail:
                    return True,superior_solutions
    for z in fn.end_points():
        for x in fn.end_points():
            y=z-x
            delta=delta_pi(fn,x,y)
            if keep_exact_solutions and delta==goal:
                exact_solutions.add((x,y))
            if delta<goal:
                superior_solutions.add((x,y))
                if stop_if_fail:
                    return True,superior_solutions
    for z in fn.end_points():
        for x in fn.end_points():
            z=1+z
            y=z-x
            delta=delta_pi(fn,x,y)
            if keep_exact_solutions and delta==goal:
                exact_solutions.add((x,y))
            if delta<goal:
                superior_solutions.add((x,y))
                if stop_if_fail:
                    return True,superior_solutions
    if len(superior_solutions)>0:
        return True,superior_solutions
    else:
        return False,superior_solutions


    
