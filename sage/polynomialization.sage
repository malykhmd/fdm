#########################
# Polynomialization for FDM
#########################

PD_POLYNOMIALIZATION_VERSION = "2026-08-02-worklist-exp-products"


def _pd_as_list(value):
    if type(value) == type([]) or type(value) == type(()):
        return list(value)
    return [value]


def _pd_expr(expr):
    return SR(expr)


def _pd_op_name(expr):
    try:
        op = expr.operator()
    except Exception:
        return ""
    if op is None:
        return ""
    return getattr(op, "__name__", str(op))


def _pd_is_add(expr):
    return "add" in _pd_op_name(expr)


def _pd_is_mul(expr):
    return "mul" in _pd_op_name(expr)


def _pd_is_pow(expr):
    return "pow" in _pd_op_name(expr)


def _pd_is_exp(expr):
    return _pd_op_name(expr) == "exp"


def _pd_function_argument(expr, name):
    operands = _pd_operands(expr)
    if _pd_op_name(expr) == name and len(operands) == 1:
        return operands[0]
    return None


def _pd_companion_substitutions(expr):
    expr = _pd_expr(expr)
    ans = set()

    arg = _pd_function_argument(expr, "sin")
    if arg is not None:
        ans.add(cos(arg))

    arg = _pd_function_argument(expr, "cos")
    if arg is not None:
        ans.add(sin(arg))

    arg = _pd_function_argument(expr, "sinh")
    if arg is not None:
        ans.add(cosh(arg))

    arg = _pd_function_argument(expr, "cosh")
    if arg is not None:
        ans.add(sinh(arg))

    return ans


def _pd_operands(expr):
    try:
        return list(expr.operands())
    except Exception:
        return []


def _pd_variables(expr):
    try:
        return set(expr.variables())
    except Exception:
        return set()


def _pd_depends_on(expr, variables):
    return len(_pd_variables(expr).intersection(set(variables))) > 0


def _pd_is_integer(value):
    try:
        return bool(value in ZZ)
    except Exception:
        pass
    try:
        return bool(SR(value).is_integer())
    except Exception:
        return False


def _pd_integer_value(value):
    if not _pd_is_integer(value):
        return None
    try:
        return ZZ(value)
    except Exception:
        return value


def _pd_is_nonnegative_integer(value):
    ivalue = _pd_integer_value(value)
    return ivalue is not None and ivalue >= 0


def _pd_is_negative_integer(value):
    ivalue = _pd_integer_value(value)
    return ivalue is not None and ivalue < 0


def _pd_reconstruct(expr, operands):
    if len(operands) == 0:
        return expr
    if _pd_is_add(expr):
        return sum(operands)
    if _pd_is_mul(expr):
        ans = SR(1)
        for item in operands:
            ans = ans * item
        return ans
    if _pd_is_pow(expr) and len(operands) == 2:
        return operands[0] ** operands[1]
    try:
        return expr.operator()(*operands)
    except Exception:
        return expr


def _pd_subs(expr, old, new):
    try:
        return expr.subs({old: new})
    except Exception:
        return expr.subs(old == new)


def _pd_expr_key(expr):
    try:
        return (len(str(expr)), str(expr))
    except Exception:
        return (0, "")


def _pd_same_expr(left, right):
    try:
        if str(left) == str(right):
            return True
    except Exception:
        pass
    try:
        return bool(SR(left) == SR(right))
    except Exception:
        return False


def _pd_simplify(expr):
    return SR(expr)


def _pd_key_from_value(dictionary, value):
    value_string = str(value)
    for key in dictionary:
        try:
            if str(dictionary[key]) == value_string:
                return key
        except Exception:
            pass
    return None


def _pd_find_power_value(dictionary, base, exp):
    for key in dictionary:
        value = dictionary[key]
        operands = _pd_operands(value)
        if _pd_is_pow(value) and len(operands) == 2:
            if _pd_same_expr(operands[0], base):
                ratio = _pd_integer_ratio(exp, operands[1])
                if ratio is not None:
                    return key, ratio
    return None, None


def _pd_integer_ratio(numerator, denominator):
    numerator = SR(numerator)
    denominator = SR(denominator)
    for n in range(-20, 21):
        if n != 0 and _pd_same_expr(numerator, n * denominator):
            return ZZ(n)
    if _pd_same_expr(numerator, denominator):
        return ZZ(1)
    return None


def _pd_exp_argument(expr):
    operands = _pd_operands(expr)
    if _pd_is_exp(expr) and len(operands) == 1:
        return operands[0]
    return None


def _pd_known_expression(system, expr):
    expr_string = str(expr)
    for x in system.variables.state:
        try:
            if str(x) == expr_string:
                return True
        except Exception:
            pass
    for value in system.substitutions.values():
        try:
            if str(value) == expr_string:
                return True
        except Exception:
            pass
    return False


def _pd_exp_argument_buildable_from_values(values, argument):
    pseudo = {}
    for i, value in enumerate(values):
        pseudo[i] = value

    key, ratio = _pd_find_exp_value(pseudo, argument)
    if key is not None:
        return True

    argument_operands = _pd_operands(argument)
    if _pd_is_add(argument) and len(argument_operands) > 1:
        for item in argument_operands:
            key, ratio = _pd_find_exp_value(pseudo, item)
            if key is None:
                return False
        return True

    return False


def _pd_new_substitution_candidates(system, expressions):
    normalized = set()
    for item in expressions:
        candidate = system._normalize_candidate(item)
        if not _pd_known_expression(system, candidate):
            normalized.add(candidate)
        for companion in _pd_companion_substitutions(candidate):
            companion = system._normalize_candidate(companion)
            if not _pd_known_expression(system, companion):
                normalized.add(companion)

    result = []
    for item in sorted(list(normalized), key=_pd_expr_key):
        if _pd_known_expression(system, item):
            continue
        exp_argument = _pd_exp_argument(item)
        if exp_argument is not None:
            known_values = list(system.substitutions.values()) + result
            if _pd_exp_argument_buildable_from_values(known_values, exp_argument):
                continue
        if _pd_is_pow(item):
            found, ratio = _pd_find_power_value(
                system.substitutions, _pd_operands(item)[0], _pd_operands(item)[1])
            if found is not None:
                continue
        result.append(item)
    return result


def _pd_find_exp_value(dictionary, argument):
    for key in dictionary:
        value = dictionary[key]
        stored_argument = _pd_exp_argument(value)
        if stored_argument is None:
            continue
        ratio = _pd_integer_ratio(argument, stored_argument)
        if ratio is not None:
            return key, ratio
    return None, None


def _pd_rewrite_exp_product(dictionary, argument):
    operands = _pd_operands(argument)
    if not (_pd_is_add(argument) and len(operands) > 1):
        return None

    result = SR(1)
    for item in operands:
        exp_var, ratio = _pd_find_exp_value(dictionary, item)
        if exp_var is None:
            return None
        result = result * (exp_var ** ratio)
    return result


def _pd_is_laurent_monom(expr, laurent_variables):
    return _pd_is_laurent_base(expr, laurent_variables)


def _pd_is_laurent_base(expr, laurent_variables):
    expr = _pd_expr(expr)
    if expr in laurent_variables:
        return True
    operands = _pd_operands(expr)
    if _pd_is_mul(expr):
        return all([_pd_is_laurent_base(item, laurent_variables) for item in operands])
    if _pd_is_pow(expr) and len(operands) == 2:
        base, exp = operands
        return _pd_is_integer(exp) and _pd_is_laurent_base(base, laurent_variables)
    return False


def pd_is_polynomial_expr(expr, variables, laurent_variables=None):
    variables = _pd_as_list(variables)
    laurent_variables = set([] if laurent_variables is None else laurent_variables)
    expr = _pd_expr(expr)

    if not _pd_depends_on(expr, variables):
        return True

    operands = _pd_operands(expr)
    if len(operands) == 0:
        return True
    if _pd_is_add(expr) or _pd_is_mul(expr):
        return all([pd_is_polynomial_expr(item, variables, laurent_variables)
                    for item in operands])
    if _pd_is_pow(expr) and len(operands) == 2:
        base, exp = operands
        if _pd_is_nonnegative_integer(exp):
            return pd_is_polynomial_expr(base, variables, laurent_variables)
        if _pd_is_negative_integer(exp) and _pd_is_laurent_base(base, laurent_variables):
            return True
        return False
    return False


def pd_find_nonpolynomial_terms(expr, variables, laurent_variables=None):
    expr = _pd_expr(expr)
    variables = _pd_as_list(variables)
    laurent_variables = set([] if laurent_variables is None else laurent_variables)

    if pd_is_polynomial_expr(expr, variables, laurent_variables):
        return set()
    if not _pd_depends_on(expr, variables):
        return set()

    operands = _pd_operands(expr)
    if _pd_is_add(expr) or _pd_is_mul(expr):
        ans = set()
        for item in operands:
            ans = ans.union(pd_find_nonpolynomial_terms(item, variables, laurent_variables))
        return ans

    ans = set([expr])
    for item in operands:
        ans = ans.union(pd_find_nonpolynomial_terms(item, variables, laurent_variables))
    return ans


class PDVariables:
    def __init__(self, variables, parameters=None, inputs=None,
                 new_var_base_name="w_", start_new_vars_with=0,
                 allow_laurent=True, laurent_variables=None):
        self.state = list(variables)
        self.parameter = set([] if parameters is None else parameters)
        self.input = set([] if inputs is None else inputs)
        self.generated = []
        self.base_name = new_var_base_name
        self.start_id = start_new_vars_with
        if laurent_variables is None:
            self.laurent = set(self.state).union(self.input) if allow_laurent else set()
        else:
            self.laurent = set(laurent_variables)
        self.allow_laurent = allow_laurent

    def create(self):
        index = self.start_id + len(self.generated)
        name = self.base_name + str(index)
        new_var = var(name)
        self.state.append(new_var)
        self.generated.append(new_var)
        if self.allow_laurent:
            self.laurent.add(new_var)
        return new_var

    def copy(self):
        variables = PDVariables([], [])
        variables.state = list(self.state)
        variables.parameter = set(self.parameter)
        variables.input = set(self.input)
        variables.generated = list(self.generated)
        variables.base_name = self.base_name
        variables.start_id = self.start_id
        variables.laurent = set(self.laurent)
        variables.allow_laurent = self.allow_laurent
        return variables


class PolynomializationSystem:
    def __init__(self, variables, rhs, parameters=None, inputs=None,
                 time_var=None, x0=None, T=None, allow_laurent=True,
                 laurent_variables=None):
        variables = _pd_as_list(variables)
        rhs = _pd_as_list(rhs)
        if len(variables) != len(rhs):
            raise ValueError("The number of variables and equations must be equal.")

        self.variables = PDVariables(variables, parameters, inputs,
                                     allow_laurent=allow_laurent,
                                     laurent_variables=laurent_variables)
        self.equations = {}
        for x, fx in zip(variables, rhs):
            self.equations[x] = _pd_expr(fx).expand()
        self.substitutions = {}
        self.time_var = time_var
        self.x0 = None if x0 is None else _pd_as_list(x0)
        self.T = T

    def copy(self):
        system = PolynomializationSystem([], [])
        system.variables = self.variables.copy()
        system.equations = dict(self.equations)
        system.substitutions = dict(self.substitutions)
        system.time_var = self.time_var
        system.x0 = None if self.x0 is None else list(self.x0)
        system.T = self.T
        return system

    def state_variables(self):
        return list(self.variables.state)

    def generated_variables(self):
        return list(self.variables.generated)

    def introduced_variables(self):
        return [(w, self.substitutions[w]) for w in self.variables.generated]

    def rhs(self, use_polynomial=True):
        if use_polynomial:
            return [self.polynomial_rhs(x) for x in self.variables.state]
        return [self.equations[x] for x in self.variables.state]

    def polynomial_rhs(self, x):
        fx = self._try_convert_to_polynomial(self.equations[x])
        if fx is None:
            return self.equations[x]
        return fx.expand()

    def is_polynomial(self):
        return all([self._try_convert_to_polynomial(self.equations[x]) is not None
                    for x in self.variables.state])

    def add_new_var(self, substitution, new_var=None):
        substitution = _pd_expr(substitution)
        existing = _pd_key_from_value(self.substitutions, substitution)
        if existing is not None:
            return existing
        if new_var is None:
            new_var = self.variables.create()
        elif new_var not in self.variables.state:
            self.variables.state.append(new_var)
            self.variables.generated.append(new_var)
            if self.variables.allow_laurent:
                self.variables.laurent.add(new_var)

        self.substitutions[new_var] = substitution
        self.equations[new_var] = self.lie_derivative(substitution)
        return new_var

    def lie_derivative(self, expr):
        expr = _pd_expr(expr)
        ans = SR(0)
        operands = _pd_operands(expr)

        if _pd_is_pow(expr) and len(operands) == 2:
            base, exp = operands
            if not _pd_is_nonnegative_integer(exp):
                base_derivative = self.lie_derivative(base)
                diff_variables = list(self.variables.state) + list(self.variables.input)
                exp_derivative = self.lie_derivative(exp) if _pd_depends_on(exp, diff_variables) else SR(0)
                ans = expr * exp * base_derivative / base
                if exp_derivative != 0:
                    ans = ans + expr * log(base) * exp_derivative
                return ans

        for x in self.variables.state:
            if x in _pd_variables(expr):
                ans = ans + diff(expr, x) * self.equations[x]

        for u in list(self.variables.input):
            if u in _pd_variables(expr):
                du = var(str(u) + "_dot")
                self.variables.input.add(du)
                if self.variables.allow_laurent:
                    self.variables.laurent.add(du)
                ans = ans + diff(expr, u) * du

        if self.time_var is not None:
            try:
                ans = ans + diff(expr, self.time_var)
            except Exception:
                pass

        return ans

    def available_substitutions(self):
        variables = set(self.variables.state).union(self.variables.input)
        ans = set()
        for x in self.variables.state:
            if self._try_convert_to_polynomial(self.equations[x]) is None:
                ans = ans.union(pd_find_nonpolynomial_terms(
                    self.equations[x], list(variables), self.variables.laurent))

        normalized = set()
        for item in ans:
            normalized.add(item)
        return _pd_new_substitution_candidates(self, normalized)

    def _normalize_candidate(self, expr):
        expr = _pd_expr(expr)
        operands = _pd_operands(expr)
        if _pd_is_pow(expr) and len(operands) == 2:
            base, exp = operands
            if _pd_is_negative_integer(exp) and not (base in self.variables.laurent):
                return SR(1) / base
        return expr

    def _try_convert_to_polynomial(self, expr):
        replaced = self._rewrite_with_substitutions(_pd_expr(expr))
        variables = list(self.variables.state) + list(self.variables.input)
        if pd_is_polynomial_expr(replaced, variables, self.variables.laurent):
            return replaced
        return None

    def _rewrite_with_substitutions(self, expr):
        expr = _pd_expr(expr)
        replacement = _pd_key_from_value(self.substitutions, expr)
        if replacement is not None:
            return replacement

        exp_argument = _pd_exp_argument(expr)
        if exp_argument is not None:
            exp_var, ratio = _pd_find_exp_value(self.substitutions, exp_argument)
            if exp_var is not None:
                return exp_var ** ratio
            exp_product = _pd_rewrite_exp_product(self.substitutions, exp_argument)
            if exp_product is not None:
                return exp_product

        operands = _pd_operands(expr)
        if _pd_is_pow(expr) and len(operands) == 2:
            base = self._rewrite_with_substitutions(operands[0])
            exp = operands[1]

            pow_var, ratio = _pd_find_power_value(self.substitutions, operands[0], exp)
            if pow_var is not None:
                return pow_var ** ratio

            inv_var = _pd_key_from_value(self.substitutions, SR(1) / operands[0])
            if inv_var is not None and _pd_is_negative_integer(exp):
                return inv_var ** (-_pd_integer_value(exp))
            return base ** exp

        if len(operands) == 0:
            return expr

        new_operands = [self._rewrite_with_substitutions(item) for item in operands]
        replaced = _pd_reconstruct(expr, new_operands)

        replacement = _pd_key_from_value(self.substitutions, replaced)
        if replacement is not None:
            return replacement

        exp_argument = _pd_exp_argument(replaced)
        if exp_argument is not None:
            exp_var, ratio = _pd_find_exp_value(self.substitutions, exp_argument)
            if exp_var is not None:
                return exp_var ** ratio
            exp_product = _pd_rewrite_exp_product(self.substitutions, exp_argument)
            if exp_product is not None:
                return exp_product

        return replaced

    def extended_initial_conditions(self):
        if self.x0 is None:
            return None
        values = list(self.x0)
        subs = {}
        for x, value in zip(self.variables.state[:len(self.x0)], self.x0):
            subs[x] = value
        for w in self.variables.generated:
            value = self.substitutions[w]
            for key in subs:
                value = _pd_subs(value, key, subs[key])
            values.append(value)
            subs[w] = value
        return values

    def initial_problem(self):
        if "Initial_problem" not in globals():
            raise NameError("Initial_problem is not loaded. Load fdm.sage before calling initial_problem().")
        return Initial_problem(self.state_variables(), self.rhs(True),
                               self.extended_initial_conditions(), self.T)

    def print_substitutions(self):
        print("Introduced variables:")
        for w in self.variables.generated:
            print(str(w) + " = " + str(self.substitutions[w]))

    def print_equations(self, use_polynomial=True):
        for x, fx in zip(self.variables.state, self.rhs(use_polynomial)):
            print(str(x) + "' = " + str(fx))

    def print(self, use_polynomial=True, with_substitutions=True):
        if with_substitutions:
            self.print_substitutions()
            print("")
        self.print_equations(use_polynomial)

    def __len__(self):
        return len(self.variables.state)


def pd_system_from_problem(problem, parameters=None, inputs=None, time_var=None,
                           allow_laurent=True, laurent_variables=None):
    data = problem.list()
    f = data[0]
    x = data[1]
    x0 = data[2]
    T = data[3]
    return PolynomializationSystem(x, f, parameters=parameters, inputs=inputs,
                                   time_var=time_var, x0=x0, T=T,
                                   allow_laurent=allow_laurent,
                                   laurent_variables=laurent_variables)


def pd_system_from_pairs(system, parameters=None, inputs=None, time_var=None,
                         allow_laurent=True, laurent_variables=None):
    x = [item[0] for item in system]
    f = [item[1] for item in system]
    return PolynomializationSystem(x, f, parameters=parameters, inputs=inputs,
                                   time_var=time_var,
                                   allow_laurent=allow_laurent,
                                   laurent_variables=laurent_variables)


def pd_polynomialization_system(system, parameters=None, inputs=None, time_var=None,
                                allow_laurent=True, laurent_variables=None):
    if isinstance(system, PolynomializationSystem):
        return system.copy()
    if hasattr(system, "list"):
        return pd_system_from_problem(system, parameters=parameters, inputs=inputs,
                                      time_var=time_var,
                                      allow_laurent=allow_laurent,
                                      laurent_variables=laurent_variables)
    return pd_system_from_pairs(system, parameters=parameters, inputs=inputs,
                                time_var=time_var,
                                allow_laurent=allow_laurent,
                                laurent_variables=laurent_variables)


def pd_apply_substitution(system, substitution):
    new_system = system.copy()
    substitution = system._normalize_candidate(substitution)
    new_system.add_new_var(substitution)
    return new_system


def _pd_poly_algo_step(system, upper_bound, best_nvars):
    if system.is_polynomial():
        return (len(system.variables.generated), system)

    if len(system.variables.generated) >= upper_bound:
        return (+Infinity, None)
    if len(system.variables.generated) >= best_nvars - 1:
        return (+Infinity, None)

    substitutions = system.available_substitutions()
    if len(substitutions) == 0:
        return (+Infinity, None)

    best_system = None
    min_nvars = best_nvars
    for substitution in substitutions:
        next_system = pd_apply_substitution(system, substitution)
        nvars, candidate = _pd_poly_algo_step(next_system, upper_bound, min_nvars)
        if nvars < min_nvars:
            min_nvars = nvars
            best_system = candidate

    return (min_nvars, best_system)


def _pd_poly_greedy(system, upper_bound, max_layers=100, trace=False):
    current = system.copy()
    for layers in range(0, max_layers + 1):
        if trace:
            print("pd greedy: layer {}, generated {}".format(
                layers, len(current.variables.generated)))
        if current.is_polynomial():
            return current
        if layers > max_layers:
            return None
        if len(current.variables.generated) >= upper_bound:
            return None
        substitutions = current.available_substitutions()
        if trace:
            print("pd greedy: candidates {}".format([str(s) for s in substitutions]))
        if len(substitutions) == 0:
            return None
        for substitution in substitutions:
            if len(current.variables.generated) >= upper_bound:
                return None
            current = pd_apply_substitution(current, substitution)
    return None


def _pd_poly_worklist(system, upper_bound, max_steps=10000, trace=False):
    current = system.copy()
    queue = list(current.variables.state)
    head = 0
    steps = 0

    while head < len(queue):
        steps = steps + 1
        if steps > max_steps:
            raise RuntimeError("Polynomialization stopped: max_steps was exceeded.")
        if len(current.variables.generated) >= upper_bound:
            return None

        x = queue[head]
        head = head + 1
        rhs = current._rewrite_with_substitutions(current.equations[x])
        try:
            rhs = rhs.expand()
        except Exception:
            pass
        current.equations[x] = rhs

        variables = list(current.variables.state) + list(current.variables.input)
        terms = pd_find_nonpolynomial_terms(rhs, variables, current.variables.laurent)
        candidates = _pd_new_substitution_candidates(current, terms)

        if trace:
            print("pd worklist: equation {}, generated {}, candidates {}".format(
                x, len(current.variables.generated), [str(c) for c in candidates]))

        if len(candidates) == 0:
            continue

        for candidate in candidates:
            if len(current.variables.generated) >= upper_bound:
                return None
            new_var = current.add_new_var(candidate)
            if new_var is not None:
                queue.append(new_var)
        queue.append(x)

    for x in list(current.variables.state):
        current.equations[x] = current._rewrite_with_substitutions(current.equations[x])

    return current if current.is_polynomial() else None


def polynomialize(system, upper_bound=10, new_vars_name="w_", start_new_vars_with=0,
                  parameters=None, inputs=None, time_var=None,
                  allow_laurent=True, laurent_variables=None, search="worklist",
                  fallback_upper_bound=50, max_steps=10000, trace=False):
    start_system = pd_polynomialization_system(system, parameters=parameters,
                                               inputs=inputs, time_var=time_var,
                                               allow_laurent=allow_laurent,
                                               laurent_variables=laurent_variables)
    start_system.variables.base_name = new_vars_name
    start_system.variables.start_id = start_new_vars_with

    if search == "worklist":
        result = _pd_poly_worklist(start_system, upper_bound, max_steps=max_steps, trace=trace)
        if result is None and fallback_upper_bound is not None and fallback_upper_bound > upper_bound:
            result = _pd_poly_worklist(start_system, fallback_upper_bound,
                                       max_steps=max_steps, trace=trace)
    elif search == "branch":
        nvars, result = _pd_poly_algo_step(start_system, upper_bound, +Infinity)
    else:
        result = _pd_poly_greedy(start_system, upper_bound, trace=trace)
        if result is None and fallback_upper_bound is not None and fallback_upper_bound > upper_bound:
            result = _pd_poly_greedy(start_system, fallback_upper_bound, trace=trace)
    if result is None:
        raise ValueError("Polynomialization was not found within the given upper_bound.")
    return result


def pd_polynomialize_problem(problem, upper_bound=10, new_vars_name="w_",
                             start_new_vars_with=0, parameters=None, inputs=None,
                             time_var=None, allow_laurent=True,
                             laurent_variables=None, search="worklist",
                             fallback_upper_bound=50, max_steps=10000,
                             trace=False, print_result=False):
    poly_system = polynomialize(problem, upper_bound=upper_bound,
                                new_vars_name=new_vars_name,
                                start_new_vars_with=start_new_vars_with,
                                parameters=parameters, inputs=inputs,
                                time_var=time_var,
                                allow_laurent=allow_laurent,
                                laurent_variables=laurent_variables,
                                search=search,
                                fallback_upper_bound=fallback_upper_bound,
                                max_steps=max_steps, trace=trace)
    poly_problem = poly_system.initial_problem()
    if print_result:
        poly_system.print()
    return [poly_problem, poly_system]


def pd_problem(x, f, x0, T):
    if "Initial_problem" not in globals():
        raise NameError("Initial_problem is not loaded. Load fdm.sage before calling pd_problem().")
    return Initial_problem(x, f, x0, T)


def pd_print_problem(problem):
    data = problem.list()
    f = data[0]
    x = data[1]
    x0 = data[2]
    T = data[3]
    print("Variables:")
    print("  {}".format(x))
    print("Equations:")
    for xi, fi in zip(x, f):
        print("  {}' = {}".format(xi, fi))
    print("Initial values:")
    for xi, xi0 in zip(x, x0):
        print("  {}(0) = {}".format(xi, xi0))
    print("T = {}".format(T))


EquationSystem = PolynomializationSystem


def eq_list_to_eq_system(system, parameters=None, inputs=None, time_var=None,
                         allow_laurent=True, laurent_variables=None):
    return pd_system_from_pairs(system, parameters=parameters, inputs=inputs,
                                time_var=time_var, allow_laurent=allow_laurent,
                                laurent_variables=laurent_variables)


def available_substitutions(system):
    return set(system.available_substitutions())


def apply_substitution(system, substitution):
    return pd_apply_substitution(system, substitution)


def find_nonpolynomial_terms(expr, variables, laurent_variables=None):
    return pd_find_nonpolynomial_terms(expr, variables, laurent_variables)


def is_pow_with_negative_integer_exp(expr):
    operands = _pd_operands(expr)
    return _pd_is_pow(expr) and len(operands) == 2 and _pd_is_negative_integer(operands[1])


def is_nonpolynomial_function(expr, variables):
    return not pd_is_polynomial_expr(expr, variables)
