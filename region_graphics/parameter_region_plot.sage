from sage.plot.contour_plot import equify
from sage.plot.contour_plot import ContourPlot
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options, suboptions
from sage.plot.colors import rgbcolor, get_cmap
from sage.misc.misc import xsrange
import operator

@options(plot_points=100, incol='blue', outcol=None, bordercol=None, borderstyle=None, borderwidth=None,frame=False,axes=True, legend_label=None, aspect_ratio=1, alpha=1)
def region_plot_patch(f, xrange, yrange, plot_points, incol, outcol, bordercol, borderstyle, borderwidth, alpha, **options):
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    import numpy

    if not isinstance(f, (list, tuple)):
        f = [f]

    f = [equify(g) for g in f]

    g, ranges = setup_for_eval_on_grid(f, [xrange, yrange], plot_points)
    xrange,yrange=[r[:2] for r in ranges]

    xy_data_arrays = numpy.asarray([[[func(x, y) for x in xsrange(*ranges[0], include_endpoint=True)]
                                     for y in xsrange(*ranges[1], include_endpoint=True)]
                                    for func in g],dtype=float)
    xy_data_array=numpy.abs(xy_data_arrays.prod(axis=0))
    # Now we need to set entries to negative iff all
    # functions were negative at that point.
    neg_indices = (xy_data_arrays<0).all(axis=0)
    xy_data_array[neg_indices]=-xy_data_array[neg_indices]

    from matplotlib.colors import ListedColormap
    incol = rgbcolor(incol)
    if outcol:
        outcol = rgbcolor(outcol)
        cmap = ListedColormap([incol, outcol])
        cmap.set_over(outcol, alpha=alpha)
    else:
        outcol = rgbcolor('white')
        cmap = ListedColormap([incol, outcol])
        cmap.set_over(outcol, alpha=0)
    cmap.set_under(incol, alpha=alpha)

    g = Graphics()

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        options['aspect_ratio'] = 'automatic'

    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore=['xmin', 'xmax']))
    g.add_primitive(ContourPlot(xy_data_array, xrange,yrange,
                                dict(contours=[-1e-20, 0, 1e-20], cmap=cmap, fill=True, **options)))

    if bordercol or borderstyle or borderwidth:
        cmap = [rgbcolor(bordercol)] if bordercol else ['black']
        linestyles = [borderstyle] if borderstyle else None
        linewidths = [borderwidth] if borderwidth else None
        g.add_primitive(ContourPlot(xy_data_array, xrange, yrange,
                                    dict(linestyles=linestyles, linewidths=linewidths,
                                         contours=[0], cmap=[bordercol], fill=False, **options)))

    return g

def linesplot(K, fill=False, color="blue", linewidth=1, legend_label=None, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1):
    x,y = var('x,y')
    leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
    if leq:
        print "WARNING: equation list is not empty!"
    if fill:
        g = region_plot([ lhs(x, y) < 0 for lhs in lin ], (x, xmin, xmax), (y, ymin, ymax), incol=color, plot_points=1000)
        return g
    g = Graphics()
    for l in lin:
        g += implicit_plot(l(x, y) == 0, (x, xmin, xmax), (y, ymin, ymax), color=color, linewidth=linewidth)
    g += line([(0,0),(0,1)], color = color, legend_label=legend_label, zorder=-10)
    return g

def regionplot(K, color="blue", alpha=0.5, legend_label=None, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, plot_points=1000):
    x,y = var('x,y')
    leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
    if leq:
        print "WARNING: equation list is not empty!"
    g = region_plot_patch([ lhs(x, y) < 0 for lhs in lin ], (x, xmin, xmax), (y, ymin, ymax), incol=color, alpha=alpha, plot_points=plot_points, bordercol=color)
    g += line([(0,0),(0,1)], color = color, legend_label=legend_label, alpha = alpha, zorder=-10)
    return g

def regionplot_two_parameters(fun_name="drlm_backward_3_slope", var_name=('f','b'), var_value=[1/12-1/30, 2/12], \
                              xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, plot_points=1000):
    """
    sage: logging.disable(logging.INFO)
    sage: g, fname, fparam = regionplot_two_parameters("drlm_backward_3_slope", ('f','b'), [1/12-1/30, 2/12])
    sage: g.save(fname+".pdf", title=fname+fparam, legend_loc=7)
    sage: g, fname, fparam = regionplot_two_parameters("drlm_backward_3_slope", ('f','b'), [1/12+1/30, 2/12])
    sage: g, fname, fparam = regionplot_two_parameters("gj_2_slope", ('f','lam'), [3/5, 1/6], plot_points=100)
    """
    K = SymbolicRealNumberField(var_value, var_name)
    if len(var_name) == 2:
        h = eval(fun_name)(K.gens()[0], K.gens()[1], field=K, conditioncheck=False)
    else:
        raise NotImplementedError, "Not 2 parameters. Not implemented."
    reg_cf = regionplot(K, color="orange", alpha=0.5, legend_label="construction", \
                        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
    minimality_test(h)
    reg_cfm = regionplot(K, color="green", alpha=0.5, legend_label="min_test", \
                        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
    extremality_test(h)
    reg_cfe = regionplot(K, color="blue", alpha=0.5, legend_label="ext_test", \
                        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
    K = SymbolicRealNumberField(var_value, var_name)
    if len(var_name) == 2:
        h = eval(fun_name)(K.gens()[0], K.gens()[1], field=K, conditioncheck=True)
    else:
        raise NotImplementedError, "Not 2 parameters. Not implemented."
    reg_ct  = regionplot(K, color="red", alpha=0.5, legend_label="conditioncheck", \
                        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
    p = point([K._values],color = "white", size = 10, zorder=5)
    g = reg_cf+reg_ct+reg_cfm+reg_cfe+p
    fun_param = "(%s=%s, %s=%s)" % (var_name[0], var_value[0], var_name[1], var_value[1]) 
    g.show(show_legend=False, title=fun_name+fun_param) #, axes_labels=var_name)
    return g, fun_name, fun_param

