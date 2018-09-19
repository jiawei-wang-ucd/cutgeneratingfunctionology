# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

import random 

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
    
def histogram_delta_pi(fn,sampling='vertices',q=1,epsilon=1/10000):
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


    