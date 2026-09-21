#########################
# PD quadratization     #
#########################

PD_QUADRATIZATION_VERSION = "2026-08-02-apelrot-qbee-style"


def _pd_quad_as_poly_system(system, parameters=None, inputs=None, time_var=None,
                            allow_laurent=False):
    if isinstance(system, PolynomializationSystem):
        return system
    return pd_polynomialization_system(system, parameters=parameters, inputs=inputs,
                                       time_var=time_var,
                                       allow_laurent=allow_laurent)


def _pd_quad_zero_tuple(n):
    return tuple([0] * n)


def _pd_quad_unit_tuple(n, i):
    return tuple([1 if i == j else 0 for j in range(n)])


def _pd_quad_tuple_add(left, right):
    return tuple([left[i] + right[i] for i in range(len(left))])


def _pd_quad_tuple_sub(left, right):
    return tuple([left[i] - right[i] for i in range(len(left))])


def _pd_quad_tuple_degree(monomial):
    return sum(monomial)


def _pd_quad_tuple_abs_degree(monomial):
    return sum([abs(item) for item in monomial])


def _pd_quad_tuple_to_expr(monomial, variables):
    ans = SR(1)
    for var_item, power in zip(variables, monomial):
        if power != 0:
            ans = ans * (var_item ** power)
    return ans


def _pd_quad_expr_to_terms(expr, variables, allow_laurent=False):
    expr = _pd_expr(expr).expand()
    n = len(variables)
    variable_index = {}
    for i, x in enumerate(variables):
        variable_index[x] = i

    if expr == 0:
        return {}

    if _pd_is_add(expr):
        terms = _pd_operands(expr)
    else:
        terms = [expr]

    ans = {}
    for term in terms:
        coefficient, monomial = _pd_quad_parse_term(term, variables,
                                                    variable_index,
                                                    allow_laurent)
        if coefficient != 0:
            ans[monomial] = ans.get(monomial, SR(0)) + coefficient

    cleaned = {}
    for monomial, coefficient in ans.items():
        if coefficient != 0:
            cleaned[monomial] = coefficient
    return cleaned


def _pd_quad_parse_term(term, variables, variable_index, allow_laurent):
    coefficient = SR(1)
    monomial = [0] * len(variables)

    if _pd_is_mul(term):
        factors = _pd_operands(term)
    else:
        factors = [term]

    for factor in factors:
        if factor in variable_index:
            monomial[variable_index[factor]] += 1
            continue

        operands = _pd_operands(factor)
        if _pd_is_pow(factor) and len(operands) == 2 and operands[0] in variable_index:
            power = operands[1]
            if not _pd_is_integer(power):
                raise ValueError("Non-integer power in a polynomial term: {}".format(factor))
            power_value = _pd_integer_value(power)
            if power_value < 0 and not allow_laurent:
                raise ValueError("Laurent monomial is not allowed here: {}".format(factor))
            monomial[variable_index[operands[0]]] += power_value
            continue

        if _pd_depends_on(factor, variables):
            raise ValueError("Expression is not a monomial in the selected variables: {}".format(factor))
        coefficient = coefficient * factor

    if (not allow_laurent) and any([power < 0 for power in monomial]):
        raise ValueError("Laurent monomial is not allowed here: {}".format(term))
    return coefficient, tuple(monomial)


def _pd_quad_rhs_dicts(system, allow_laurent=False):
    variables = system.state_variables()
    rhs = system.rhs(True)
    return [_pd_quad_expr_to_terms(fx, variables, allow_laurent=allow_laurent)
            for fx in rhs]


def _pd_quad_max_total_degree(rhs_dicts):
    degree = 0
    for rhs in rhs_dicts:
        for monomial in rhs:
            degree = max(degree, _pd_quad_tuple_degree(monomial))
    return degree


def _pd_quad_generate_nonnegative_tuples(n, total):
    if n == 0:
        if total == 0:
            return [tuple()]
        return []
    if n == 1:
        return [(total,)]
    ans = []
    for first in range(total + 1):
        for rest in _pd_quad_generate_nonnegative_tuples(n - 1, total - first):
            ans.append(tuple([first] + list(rest)))
    return ans


def _pd_quad_apelrot_monomials(n, max_degree):
    ans = []
    for degree in range(2, max_degree):
        ans.extend(_pd_quad_generate_nonnegative_tuples(n, degree))
    return ans


def _pd_quad_build_square_map(original_dim, introduced_monomials):
    generalized_monomials = []
    for i in range(original_dim):
        generalized_monomials.append(_pd_quad_unit_tuple(original_dim, i))
    generalized_monomials.extend(list(introduced_monomials))

    squares = {}
    squares[_pd_quad_zero_tuple(original_dim)] = _pd_quad_zero_tuple(len(generalized_monomials))

    for i, monomial in enumerate(generalized_monomials):
        squares[monomial] = _pd_quad_unit_tuple(len(generalized_monomials), i)

    for i, left in enumerate(generalized_monomials):
        for j, right in enumerate(generalized_monomials):
            product = _pd_quad_tuple_add(left, right)
            rewritten = [0] * len(generalized_monomials)
            rewritten[i] += 1
            rewritten[j] += 1
            rewritten = tuple(rewritten)
            if product not in squares:
                squares[product] = rewritten
            elif _pd_quad_tuple_abs_degree(rewritten) < _pd_quad_tuple_abs_degree(squares[product]):
                squares[product] = rewritten

    return squares


def _pd_quad_lie_derivative_dict(monomial, vector_field):
    ans = {}
    for i, power in enumerate(monomial):
        if power != 0:
            for rhs_monomial, coefficient in vector_field[i].items():
                new_monomial = list(_pd_quad_tuple_add(monomial, rhs_monomial))
                new_monomial[i] -= 1
                new_monomial = tuple(new_monomial)
                ans[new_monomial] = ans.get(new_monomial, SR(0)) + power * coefficient
    cleaned = {}
    for monomial, coefficient in ans.items():
        if coefficient != 0:
            cleaned[monomial] = coefficient
    return cleaned


def _pd_quad_rewrite_dict(rhs_dict, square_map, generalized_variables):
    ans = SR(0)
    missing = []
    for monomial, coefficient in rhs_dict.items():
        if monomial not in square_map:
            missing.append(monomial)
        else:
            ans = ans + coefficient * _pd_quad_tuple_to_expr(square_map[monomial],
                                                            generalized_variables)
    if len(missing) > 0:
        raise ValueError("Quadratization does not cover monomials: {}".format(missing))
    return ans.expand()


def _pd_quad_apply(system, introduced_monomials, method_name,
                   new_vars_name="z_", start_new_vars_with=0):
    variables = system.state_variables()
    n = len(variables)
    rhs_dicts = _pd_quad_rhs_dicts(system, allow_laurent=True)

    quad_variables = []
    quad_substitutions = {}
    for idx, monomial in enumerate(introduced_monomials):
        z = var(new_vars_name + str(start_new_vars_with + idx))
        quad_variables.append(z)
        quad_substitutions[z] = _pd_quad_tuple_to_expr(monomial, variables)

    generalized_variables = list(variables) + list(quad_variables)
    square_map = _pd_quad_build_square_map(n, introduced_monomials)

    all_rhs_dicts = list(rhs_dicts)
    for monomial in introduced_monomials:
        all_rhs_dicts.append(_pd_quad_lie_derivative_dict(monomial, rhs_dicts))

    quad_rhs = [_pd_quad_rewrite_dict(rhs_dict, square_map, generalized_variables)
                for rhs_dict in all_rhs_dicts]
    return QuadratizationSystem(system, generalized_variables, quad_rhs,
                                quad_variables, quad_substitutions,
                                method_name=method_name)


class QuadratizationSystem:
    def __init__(self, source_system, variables, rhs, quad_variables,
                 quad_substitutions, method_name="qbee"):
        self.source_system = source_system
        self.variables = list(variables)
        self.equations = {}
        for x, fx in zip(variables, rhs):
            self.equations[x] = _pd_expr(fx).expand()
        self.quad_variables = list(quad_variables)
        self.quad_substitutions = dict(quad_substitutions)
        self.method_name = method_name
        self.x0 = source_system.extended_initial_conditions()
        self.T = source_system.T

    def state_variables(self):
        return list(self.variables)

    def generated_variables(self):
        return list(self.quad_variables)

    def introduced_variables(self):
        ans = []
        if hasattr(self.source_system, "introduced_variables"):
            ans.extend(self.source_system.introduced_variables())
        for z in self.quad_variables:
            ans.append((z, self.quad_substitutions[z]))
        return ans

    def rhs(self):
        return [self.equations[x] for x in self.variables]

    def extended_initial_conditions(self):
        if self.source_system.x0 is None:
            return None
        values = list(self.source_system.extended_initial_conditions())
        subs = {}
        for x, value in zip(self.source_system.state_variables(), values):
            subs[x] = value
        for z in self.quad_variables:
            value = self.quad_substitutions[z]
            for key in subs:
                value = _pd_subs(value, key, subs[key])
            values.append(value)
            subs[z] = value
        return values

    def initial_problem(self):
        if "Initial_problem" not in globals():
            raise NameError("Initial_problem is not loaded. Load fdm.sage before calling initial_problem().")
        return Initial_problem(self.state_variables(), self.rhs(),
                               self.extended_initial_conditions(), self.T)

    def is_quadratic(self):
        variables = self.state_variables()
        for fx in self.rhs():
            terms = _pd_quad_expr_to_terms(fx, variables, allow_laurent=False)
            for monomial in terms:
                if _pd_quad_tuple_degree(monomial) > 2:
                    return False
        return True

    def print_substitutions(self):
        print("Introduced variables:")
        for w, expr in self.introduced_variables():
            print(str(w) + " = " + str(expr))

    def print_equations(self):
        for x, fx in zip(self.variables, self.rhs()):
            print(str(x) + "' = " + str(fx))

    def print(self, with_substitutions=True):
        if with_substitutions:
            self.print_substitutions()
            print("")
        self.print_equations()

    def __len__(self):
        return len(self.variables)


class _PDQBeeSearchSystem:
    def __init__(self, rhs_dicts):
        self.dim = len(rhs_dicts)
        self.rhs_diff = {}
        for i, rhs in enumerate(rhs_dicts):
            current = set()
            unit = _pd_quad_unit_tuple(self.dim, i)
            for monomial in rhs:
                current.add(_pd_quad_tuple_sub(monomial, unit))
            self.rhs_diff[i] = current

        self.vars = []
        self.var_set = set()
        self.squares = set()
        self.nonsquares = set()
        self.add_var(_pd_quad_zero_tuple(self.dim))
        for i in range(self.dim):
            self.add_var(_pd_quad_unit_tuple(self.dim, i))

    def copy(self):
        ans = _PDQBeeSearchSystem.__new__(_PDQBeeSearchSystem)
        ans.dim = self.dim
        ans.rhs_diff = {}
        for key, value in self.rhs_diff.items():
            ans.rhs_diff[key] = set(value)
        ans.vars = list(self.vars)
        ans.var_set = set(self.var_set)
        ans.squares = set(self.squares)
        ans.nonsquares = set(self.nonsquares)
        return ans

    def add_var(self, monomial):
        monomial = tuple(monomial)
        if monomial in self.var_set:
            return

        for i, power in enumerate(monomial):
            if power != 0:
                for rhs_monomial in self.rhs_diff[i]:
                    self.nonsquares.add(_pd_quad_tuple_add(monomial, rhs_monomial))

        for known in self.vars:
            self.squares.add(_pd_quad_tuple_add(known, monomial))
        self.squares.add(_pd_quad_tuple_add(monomial, monomial))
        self.vars.append(monomial)
        self.var_set.add(monomial)
        self.nonsquares = set([item for item in self.nonsquares if item not in self.squares])

    def introduced_vars(self):
        ans = []
        for monomial in self.vars:
            if _pd_quad_tuple_abs_degree(monomial) >= 2 or sum(monomial) < 0:
                ans.append(monomial)
        return ans

    def new_vars_count(self):
        return len(self.introduced_vars())

    def is_quadratized(self):
        return len(self.nonsquares) == 0

    def smallest_nonsquare(self):
        return min([(prod([abs(power) + 1 for power in monomial]), monomial)
                    for monomial in self.nonsquares])[1]


def _pd_quad_get_decompositions(monomial):
    if len(monomial) == 0:
        return set([(tuple(), tuple())])
    previous = _pd_quad_get_decompositions(tuple(monomial[:-1]))
    ans = set()
    last = monomial[-1]
    sign = 1
    if last < 0:
        sign = -1
    for left, right in previous:
        for amount in range(abs(last) + 1):
            a = tuple(list(left) + [sign * amount])
            b = tuple(list(right) + [last - sign * amount])
            if str(a) <= str(b):
                ans.add((a, b))
            else:
                ans.add((b, a))
    return ans


def _pd_quad_search_score(system):
    total = sum([_pd_quad_tuple_abs_degree(monomial) for monomial in system.nonsquares])
    return total + system.dim * len(system.vars) + 10 * system.new_vars_count()


def _pd_quad_next_generation(system):
    if system.is_quadratized():
        return []
    nonsquare = system.smallest_nonsquare()
    ans = []
    for decomposition in _pd_quad_get_decompositions(nonsquare):
        candidate = system.copy()
        for item in decomposition:
            candidate.add_var(item)
        if candidate.new_vars_count() > system.new_vars_count():
            ans.append(candidate)
    return sorted(ans, key=_pd_quad_search_score)


def _pd_quad_qbee_beam(rhs_dicts, max_new_vars=30, max_nodes=5000, beam_width=40):
    start = _PDQBeeSearchSystem(rhs_dicts)
    if start.is_quadratized():
        return []

    best = None
    frontier = [start]
    nodes = 0

    while len(frontier) > 0 and nodes < max_nodes:
        next_frontier = []
        for system in frontier:
            nodes += 1
            if nodes > max_nodes:
                break
            if system.is_quadratized():
                if best is None or system.new_vars_count() < best.new_vars_count():
                    best = system
                continue
            if system.new_vars_count() >= max_new_vars:
                continue
            if best is not None and system.new_vars_count() >= best.new_vars_count():
                continue
            next_frontier.extend(_pd_quad_next_generation(system))

        if best is not None:
            next_frontier = [item for item in next_frontier
                             if item.new_vars_count() < best.new_vars_count()]
        frontier = sorted(next_frontier, key=_pd_quad_search_score)[:beam_width]

    if best is None:
        raise ValueError("QBee-style quadratization was not found within max_nodes/max_new_vars.")
    return best.introduced_vars()


def _pd_quad_qbee_greedy(rhs_dicts, max_new_vars=100, max_steps=1000):
    system = _PDQBeeSearchSystem(rhs_dicts)
    steps = 0
    while (not system.is_quadratized()) and steps < max_steps:
        if system.new_vars_count() >= max_new_vars:
            break
        generation = _pd_quad_next_generation(system)
        if len(generation) == 0:
            break
        system = generation[0]
        steps += 1
    if not system.is_quadratized():
        raise ValueError("Greedy QBee-style quadratization was not found within the limits.")
    return system.introduced_vars()


def _pd_quad_qbee_bnb(rhs_dicts, max_new_vars=30, max_nodes=5000):
    start = _PDQBeeSearchSystem(rhs_dicts)
    if start.is_quadratized():
        return []

    best = [None]
    nodes = [0]

    try:
        greedy_vars = _pd_quad_qbee_greedy(rhs_dicts, max_new_vars=max_new_vars,
                                           max_steps=max_nodes)
        greedy_system = _PDQBeeSearchSystem(rhs_dicts)
        for monomial in greedy_vars:
            greedy_system.add_var(monomial)
        if greedy_system.is_quadratized():
            best[0] = greedy_system
    except Exception:
        pass

    def visit(system):
        nodes[0] += 1
        if nodes[0] > max_nodes:
            return
        if system.is_quadratized():
            if best[0] is None or system.new_vars_count() < best[0].new_vars_count():
                best[0] = system
            return
        if system.new_vars_count() >= max_new_vars:
            return
        if best[0] is not None and system.new_vars_count() >= best[0].new_vars_count():
            return
        for candidate in _pd_quad_next_generation(system):
            visit(candidate)

    visit(start)
    if best[0] is None:
        raise ValueError("QBee branch-and-bound quadratization was not found within max_nodes/max_new_vars.")
    return best[0].introduced_vars()


def pd_quadratize(system, method="qbee", polynomialize_first=False,
                  parameters=None, inputs=None, time_var=None,
                  allow_laurent=False, polynomialization_upper_bound=20,
                  new_vars_name="z_", start_new_vars_with=0,
                  max_new_vars=30, max_nodes=5000, beam_width=40,
                  print_result=False):
    poly_system = _pd_quad_as_poly_system(system, parameters=parameters,
                                          inputs=inputs, time_var=time_var,
                                          allow_laurent=allow_laurent)
    if (not poly_system.is_polynomial()) and polynomialize_first:
        poly_system = polynomialize(system, upper_bound=polynomialization_upper_bound,
                                    parameters=parameters, inputs=inputs,
                                    time_var=time_var,
                                    allow_laurent=allow_laurent)
    if not poly_system.is_polynomial():
        raise ValueError("Quadratization requires a polynomialized system.")

    method_key = str(method).lower()
    is_apelrot = method_key in ["apelrot", "abh1", "appelroth"]

    rhs_dicts = _pd_quad_rhs_dicts(poly_system, allow_laurent=(not is_apelrot))

    if is_apelrot:
        max_degree = _pd_quad_max_total_degree(rhs_dicts)
        if any([any([power < 0 for power in monomial])
                for rhs in rhs_dicts for monomial in rhs]):
            raise ValueError("The Apelrot method requires ordinary polynomial RHS, not Laurent monomials.")
        introduced = _pd_quad_apelrot_monomials(len(poly_system.state_variables()),
                                                max_degree)
        result = _pd_quad_apply(poly_system, introduced, "apelrot",
                                new_vars_name=new_vars_name,
                                start_new_vars_with=start_new_vars_with)
    elif method_key == "qbee":
        introduced = _pd_quad_qbee_beam(rhs_dicts, max_new_vars=max_new_vars,
                                        max_nodes=max_nodes,
                                        beam_width=beam_width)
        result = _pd_quad_apply(poly_system, introduced, "qbee",
                                new_vars_name=new_vars_name,
                                start_new_vars_with=start_new_vars_with)
    elif method_key == "qbee_greedy":
        introduced = _pd_quad_qbee_greedy(rhs_dicts, max_new_vars=max_new_vars,
                                          max_steps=max_nodes)
        result = _pd_quad_apply(poly_system, introduced, "qbee_greedy",
                                new_vars_name=new_vars_name,
                                start_new_vars_with=start_new_vars_with)
    elif method_key == "qbee_bnb":
        introduced = _pd_quad_qbee_bnb(rhs_dicts, max_new_vars=max_new_vars,
                                       max_nodes=max_nodes)
        result = _pd_quad_apply(poly_system, introduced, "qbee_bnb",
                                new_vars_name=new_vars_name,
                                start_new_vars_with=start_new_vars_with)
    else:
        raise ValueError("Unknown quadratization method: {}".format(method))

    if print_result:
        result.print()
    return result


def pd_polynomialize_and_quadratize(system, method="qbee",
                                    polynomialization_upper_bound=20,
                                    parameters=None, inputs=None, time_var=None,
                                    allow_laurent=False,
                                    new_poly_vars_name="w_",
                                    new_quad_vars_name="z_",
                                    max_new_vars=30, max_nodes=5000,
                                    beam_width=40, print_result=False):
    poly_system = polynomialize(system,
                                upper_bound=polynomialization_upper_bound,
                                parameters=parameters, inputs=inputs,
                                time_var=time_var,
                                allow_laurent=allow_laurent,
                                new_vars_name=new_poly_vars_name)
    return pd_quadratize(poly_system, method=method,
                         allow_laurent=allow_laurent,
                         new_vars_name=new_quad_vars_name,
                         max_new_vars=max_new_vars,
                         max_nodes=max_nodes,
                         beam_width=beam_width,
                         print_result=print_result)


def pd_quadratize_problem(problem, method="qbee",
                          polynomialize_first=False,
                          polynomialization_upper_bound=20,
                          parameters=None, inputs=None, time_var=None,
                          allow_laurent=False,
                          new_vars_name="z_",
                          start_new_vars_with=0,
                          max_new_vars=30, max_nodes=5000,
                          beam_width=40, print_result=False):
    quad_system = pd_quadratize(problem, method=method,
                                polynomialize_first=polynomialize_first,
                                polynomialization_upper_bound=polynomialization_upper_bound,
                                parameters=parameters, inputs=inputs,
                                time_var=time_var,
                                allow_laurent=allow_laurent,
                                new_vars_name=new_vars_name,
                                start_new_vars_with=start_new_vars_with,
                                max_new_vars=max_new_vars,
                                max_nodes=max_nodes,
                                beam_width=beam_width,
                                print_result=print_result)
    return [quad_system.initial_problem(), quad_system]


quadratize = pd_quadratize
polynomialize_and_quadratize = pd_polynomialize_and_quadratize
