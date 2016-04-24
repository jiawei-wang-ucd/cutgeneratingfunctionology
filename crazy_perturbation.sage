# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

# Reminder: need coerce all input to common RNF.

class CrazyPiece:
    def __init__(self, interval, generators=None, cosets=[]):
        # assume interval is open, represented as (a, b)
        self.interval = interval
        # assume that generators QQ linearly independent
        self.generators = generators
        if not is_QQ_linearly_independent(*generators):
            logging.warn("Generators are not linearly independent over Q.")
        if len(generators) < 2:
            logging.warn("The group is not dense.")
        self.hermite_form_generators = find_hermite_form_generators(generators)
        # cosets is a list of (coset_repr, shift),
        # assume coset_repr represent distinct cosets;
        # assume shifts are non-zero
        self.cosets = cosets
        for i in range(len(cosets)):
            (r, s) = cosets[i]
            if s == 0:
                logging.warn("Have shift = 0.")
            for j in range(i):
                rr = cosets[j][0]
                if is_in_ZZ_span(r-rr, generators):
                    logging.warn("Not unique coset representative.")

    def __call__(self, x):
        x = fractional(x)
        ###assert self.interval[0] < x < self.interval[1]
        for (r,s) in self.cosets:
            if is_in_ZZ_span(x-r, self.generators):
                return s
        return x.parent().zero()

    def __neg__(self):
        new_cosets = [(r, -s) for (r, s) in self.cosets]
        return CrazyPiece(self.interval, self.generators, new_cosets)

    def __mul__(self, other):
        """
        Multiply `self` by a scalar.
        """
        # assume scalar multiplication
        new_cosets = [(r, s*other) for (r, s) in self.cosets]
        return CrazyPiece(self.interval, self.generators, new_cosets)

    __rmul__ = __mul__

    def __div__(self, other):
        return self * (1 / other)


class PiecewiseCrazyFunction:
    # assume that all inputs are elements from a same RNF.
    def __init__(self, pwl, crazy_pieces):
        self.pwl = pwl
        # assume that crazy pieces' intervals are disjoint.
        self.crazy_pieces = crazy_pieces

    def find_crazy_piece(self, x, xeps=0):
        x = fractional(x)
        if x == 1 and xeps == 1:
            x = x.parent().zero()
        elif x == 0 and xeps == -1:
            x = x.parent().one()
        for cp in self.crazy_pieces:
            if (cp.interval[0] < x < cp.interval[1]) or \
               (xeps == 1) and (x == cp.interval[0]) or \
               (xeps == -1) and (x == cp.interval[1]):
                return cp
        return None

    def __add__(self,other):
        if isinstance(other, PiecewiseCrazyFunction):
            # assume that intervals of crazy pieces from self and from other are disjoint.
            return PiecewiseCrazyFunction(self.pwl + other.pwl, self.crazy_pieces + other.crazy_pieces)
        else:
            # assume other is FastPiecewise
            return PiecewiseCrazyFunction(self.pwl + other, self.crazy_pieces)

    def __neg__(self):
        return PiecewiseCrazyFunction(-self.pwl, [-cp  for cp in self.crazy_pieces])

    def __mul__(self, other):
        """
        Multiply `self` by a scalar.
        """
        # assume scalar multiplication
        return PiecewiseCrazyFunction(self.pwl * other, [cp * other for cp in self.crazy_pieces])

    __rmul__ = __mul__

    def __div__(self, other):
        return self * (1 / other)

    def __sub__(self, other):
        return self + (-other)

    @cached_method
    def __call__(self, x):
        crazy_piece = self.find_crazy_piece(x, 0)
        if crazy_piece is None:
            return self.pwl(x)
        else:
            return self.pwl(x) + crazy_piece(x)

    def plot(self, rgbcolor='magenta'):
        pwl = self.pwl
        g = pwl.plot(rgbcolor=rgbcolor)
        for crazy_piece in self.crazy_pieces:
            a, b = crazy_piece.interval[0], crazy_piece.interval[1]
            pwla = pwl.limit(a, 1)
            pwlb = pwl.limit(b, -1)
            shifts = [s for (r, s) in crazy_piece.cosets]
            max_shift = max(shifts)
            min_shift = min(shifts)
            g += polygon([(a, pwla+min_shift), (a, pwla+max_shift), (b, pwlb+max_shift), (b, pwlb+min_shift)], color=rgbcolor, alpha=0.1)
            for s in shifts:
                g += line([(a, pwla+s), (b, pwlb+s)], color=rgbcolor)
        return g

    def range(self):
        pwl = self.pwl
        lb = min(flatten(pwl.limits_at_end_points()))
        ub = max(flatten(pwl.limits_at_end_points()))
        for crazy_piece in self.crazy_pieces:
            pwla = pwl.limit(crazy_piece.interval[0],1)
            pwlb = pwl.limit(crazy_piece.interval[1],-1)
            max_shift = max([s for (r, s) in crazy_piece.cosets])
            min_shift = min([s for (r, s) in crazy_piece.cosets])
            if min_shift + min(pwla, pwlb) < lb:
                lb = min_shift + min(pwla, pwlb)
            if max_shift + max(pwla, pwlb) > ub:
                ub = max_shift + max(pwla, pwlb)
        return (lb, ub)

    def end_points(self):
        bkpts = copy(self.pwl.end_points())
        for crazy_piece in self.crazy_pieces:
            bkpts += [crazy_piece.interval[0], crazy_piece.interval[1]]
        return uniq(bkpts)

    def limit(self,  x0, epsilon):
        if epsilon != 0:
            raise NotImplementedError()
        else:
            return self(x0)

def is_in_ZZ_span(x, generators):
    # assume that all inputs are elements from a same RNF.
    # generators are linearly independent over Q
    lgens = [g.list() for g in generators]
    lx = x.list()
    #if rank(matrix(QQ,lgens+[lx])) != len(generators):
    #    return False
    try:
        s = matrix(QQ,lgens).solve_left(matrix(QQ,lx))
        return all((si in ZZ) for si in s[0])
    except ValueError:
        return False

def find_hermite_form_generators(generators):
    lgens = [g.list() for g in generators]
    return (matrix(QQ,lgens).hermite_form())

def find_epsilon_for_crazy_perturbation(fn, cp, show_plots=False):
    """
    EXAMPLE::

        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: bkpts = h.end_points()
        sage: t1 = bkpts[10]-bkpts[6]
        sage: t2 = bkpts[13]-bkpts[6]
        sage: f = bkpts[37]
        sage: ucl = bkpts[17]
        sage: ucr = bkpts[18]
        sage: generators = [t1, t2]
        sage: pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
        sage: crazy_piece_1 = CrazyPiece((ucl, ucr), generators, [(ucl, 1), (ucr, -1)])
        sage: crazy_piece_2 = CrazyPiece((f-ucr, f-ucl), generators, [(f-ucr, 1), (f-ucl, -1)])
        sage: cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])
        sage: find_epsilon_for_crazy_perturbation(h, cp)
        0.0002639108814623441?
    """
    # assume fn is a subadditive pwl function, cp (crazy perturbation) is a non_zero PiecewiseCrazyFunction with cp(0)=cp(f)=0.
    bkpt = uniq(copy(fn.end_points())+cp.end_points())
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    type_1_vertices = [(x, y, x+y) for x in bkpt for y in bkpt if x <= y]
    type_2_vertices = [(x, z-x, z) for x in bkpt for z in bkpt2 if x < z < 1+x]
    vertices = set(type_1_vertices + type_2_vertices)
    g = plot_2d_complex(fn)
    m = 3
    for (x, y, z) in vertices:
        for (xeps, yeps, zeps) in [(0,0,0)]+list(nonzero_eps):
            deltafn = delta_pi_general(fn, x, y, (xeps, yeps, zeps))
            if deltafn > 0:
                possible =  True
                if deltafn < m:
                    m = deltafn
            elif (xeps, yeps, zeps) == (0, 0, 0): # 0-d face
                possible = (delta_pi(cp, x, y) == 0)
            elif xeps!= 0 and yeps!=0 and zeps!=0: # 2-d face, 6 cases
                possible = ((delta_pi_general(cp.pwl, x, y, (xeps, yeps, zeps))==0)\
                            and (cp.find_crazy_piece(x, xeps) is None) \
                            and (cp.find_crazy_piece(y, yeps) is None) \
                            and (cp.find_crazy_piece(z, zeps) is None) )
            elif xeps == 0: # vertical 1-d face, 2 cases
                possible = (cp(x) + cp.pwl.limit(y, yeps) == cp.pwl.limit(fractional(z), zeps)) \
                           and check_move_on_crazy_pieces((1, x), cp.find_crazy_piece(y, yeps), cp.find_crazy_piece(z, zeps))
            elif yeps == 0: # horizontal 1-d face, 2 cases
                possible = (cp(y) + cp.pwl.limit(x, xeps) == cp.pwl.limit(fractional(z), zeps)) \
                           and check_move_on_crazy_pieces((1, y), cp.find_crazy_piece(x, xeps), cp.find_crazy_piece(z, zeps))
            elif zeps == 0: # diagonal 1-d face, 2 cases
                possible = (cp.pwl.limit(x, xeps) + cp.pwl.limit(y, yeps) == cp(fractional(z))) \
                           and check_move_on_crazy_pieces((-1, z), cp.find_crazy_piece(x, xeps), cp.find_crazy_piece(y, yeps))
            else:
                raise ValueError, "Something Wrong in enumerating limits"
            if not possible:
                #return (x, y, z, xeps, yeps, zeps)
                if not show_plots:
                    return 0
                else:
                    m = 0
                    g += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)), r=0.01)
                    g += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)), r=0.01)
    range_cp = cp.range()
    M = max([abs(range_cp[0]), abs(range_cp[1])]) * 3
    C = max([abs(fi._slope) for fi in cp.pwl._functions])
    epsilon = m / max([M, C, 1/1000])
    if epsilon == 0 and show_plots:
        g.show(figsize=40)
    return epsilon

def check_move_on_crazy_pieces((move_sign, move_dist), cp1, cp2):
    if (cp1 is None) and (cp2 is not None) or (cp1 is not None) and (cp2 is None):
        return False
    elif (cp1 is not None) and (cp2 is not None):
        #compare if the groups on cp1 and cp2 are the same
        # TODO: set up these subgroups in a high-level way in Sage, and compare.
        if not cp1.hermite_form_generators == cp2.hermite_form_generators:
            logging.warn("Different groups. Open question.")
            return False
        if move_sign == 1:
            return (all((s == cp2(r + move_dist)) for (r, s) in cp1.cosets) \
                    and all((s == cp1(r - move_dist)) for (r, s) in cp2.cosets))
        else: # move_sign == -1:
            return (all((s == - cp2(move_dist - r)) for (r, s) in cp1.cosets) \
                    and all((s == - cp1(move_dist - r)) for (r, s) in cp2.cosets))
    else: # (cp1 is None) and (cp2 is None)
        return True


def random_test_number(fn):
    if randint(0, 5)==0:
        # Pick f
        try:
            return find_f(fn.pwl)
        except AttributeError:
            return find_f(fn)
    breakpoints = fn.end_points()
    if randint(0, 5) == 0:
        # Pick a breakpoint
        return breakpoints[randint(0, len(breakpoints)-1)]
    # Pick a point from the interior of some interval
    crazy_pieces = []
    intervals = []
    try:
        crazy_pieces = fn.crazy_pieces
        intervals = fn.pwl.intervals()
    except AttributeError:
        intervals = fn.intervals()
    if crazy_pieces and randint(0, 1) == 0:
        # Pick from crazy piece
        crazy_piece = crazy_pieces[randint(0, len(crazy_pieces)-1)]
        if randint(0, 0) == 0: # Always...
            # Pick from support of microperiodic
            cosets = crazy_piece.cosets
            coset = cosets[randint(0, len(cosets)-1)][0]
        else:
            coset = ZZ(randint(0, 12345678)) / ZZ(randint(1, 1234567))
        generators = crazy_piece.generators
        x = coset + sum(randint(0, 12345678) * gen for gen in generators)
        assert generators[0] < 1
        x = x - floor(x / generators[0] * generators[0])
        i = floor((1 - x) / generators[0])
        x = x + randint(0, i) * generators[0]
        return x
    interval = intervals[randint(0, len(intervals)-1)]
    denom = 12345678
    x = interval[0] + ZZ(randint(0, denom)) / denom * (interval[1] - interval[0])
    return x

def random_6_tuple(fn):
    # FIXME: should do limits!
    if randint(0, 1) == 0:
        x = random_test_number(fn)
        y = random_test_number(fn)
        z = x + y
        return x, y, z, 0, 0, 0
    else:
        x = random_test_number(fn)
        z = randint(0, 1) + random_test_number(fn)
        y = fractional(z - x)
        return x, y, z, 0, 0, 0

def minimality_test_randomized(fn, orig_function=None, max_iterations=None):
    """
    EXAMPLE::

        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: bkpts = h.end_points()
        sage: t1 = bkpts[10]-bkpts[6]
        sage: t2 = bkpts[13]-bkpts[6]
        sage: f = bkpts[37]
        sage: ucl = bkpts[17]
        sage: ucr = bkpts[18]
        sage: generators = [t1, t2]
        sage: pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
        sage: crazy_piece_1 = CrazyPiece((ucl, ucr), generators, [(ucl, 1), (ucr, -1)])
        sage: crazy_piece_2 = CrazyPiece((f-ucr, f-ucl), generators, [(f-ucr, 1), (f-ucl, -1)])
        sage: cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])
        sage: eps = find_epsilon_for_crazy_perturbation(h, cp)
        sage: hcp = h + eps * cp
        sage: minimality_test_randomized(hcp, h, max_iterations=10)
        True
    """
    smallest_delta = 10
    num_it = 0
    while max_iterations is None or num_it < max_iterations:
        num_it = num_it + 1
        x, y, z, xeps, yeps, zeps = random_6_tuple(fn)
        delta = delta_pi_general(fn, x, y, (xeps, yeps, zeps))
        if delta < 0:
            logging.warning("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
            return False
        if 0 < delta and orig_function is not None:
            if delta_pi_general(orig_function, x, y, (xeps, yeps, zeps)) == 0:
                logging.warning("Lost additivity: pi(%s%s) + pi(%s%s) - pi(%s%s) > 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
                return False
        if 0 == delta and orig_function is not None:
            if delta_pi_general(orig_function, x, y, (xeps, yeps, zeps)) != 0:
                logging.info("New additivity: pi(%s%s) + pi(%s%s) - pi(%s%s) = 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        if 0 < delta < smallest_delta:
            smallest_delta = delta
            logging.info("After {} tries, smallest Delta pi now: {}".format(num_it, delta))
    return True
