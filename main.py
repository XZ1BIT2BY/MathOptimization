import re
import math

EPS = 1e-9

# ---------- Разбор входных данных ----------
def parse_equation(eq_line):
    s = eq_line.replace(' ', '')
    left, right = re.split('<=|>=|=', s)
    sign = re.findall('<=|>=|=', s)[0]
    coeffs = {}
    for num, var in re.findall(r'([+-]?\d*)(x\d+)', left):
        if num in ('', '+'):
            val = 1
        elif num == '-':
            val = -1
        else:
            val = int(num)   # при желании можно заменить на float(num)
        coeffs[var] = val
    return coeffs, sign, float(right)

def read_lp(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    goal = lines[0].lower()
    obj_str = lines[1]

    obj = {}
    for num, var in re.findall(r'([+-]?\d*)(x\d+)', obj_str):
        if num in ('', '+'):
            v = 1
        elif num == '-':
            v = -1
        else:
            v = int(num)
        obj[var] = v

    constraints = []
    free_vars = set()

    for ln in lines[2:]:
        low = ln.lower()
        if low.startswith('free'):
            for var in re.findall(r'(x\d+)', ln):
                free_vars.add(var)
            continue
        constraints.append(parse_equation(ln))

    free_vars = sorted(free_vars, key=lambda x: int(re.findall(r'\d+', x)[0]))
    return goal, obj, constraints, free_vars

# ---------- Обработка свободных переменных ----------
def transform_free_variables(obj, constraints, free_vars):
    """
    Каждую свободную xk заменяем на xk_p - xk_m, обе >= 0.
    """
    free_map = {}  # 'x2' -> ('x2_p', 'x2_m')
    for v in free_vars:
        free_map[v] = (f"{v}_p", f"{v}_m")

    def add_term(dic, var, coef):
        dic[var] = dic.get(var, 0) + coef

    # Целевая функция
    new_obj = {}
    for var, coef in obj.items():
        if var in free_map:
            vp, vm = free_map[var]
            add_term(new_obj, vp, coef)
            add_term(new_obj, vm, -coef)
        else:
            add_term(new_obj, var, coef)

    # Ограничения
    new_constraints = []
    for coeffs, sign, rhs in constraints:
        new_c = {}
        for var, coef in coeffs.items():
            if var in free_map:
                vp, vm = free_map[var]
                add_term(new_c, vp, coef)
                add_term(new_c, vm, -coef)
            else:
                add_term(new_c, var, coef)
        new_constraints.append((new_c, sign, rhs))

    return new_obj, new_constraints, free_map

# ---------- Приведение задачи к канонической форме ----------
def build_canonical(goal, obj, constraints):
    # сначала соберём все переменные из уже преобразованных ограничений/цели
    orig_vars = sorted(
        {v for c, _, _ in constraints for v in c} | set(obj.keys()),
        key=lambda x: int(re.findall(r'\d+', x)[0])
    )

    A_rows = []
    b = []
    row_types = []
    slack_vars = []
    art_vars = []
    added_per_row = []

    for i, (coeffs, sign, rhs) in enumerate(constraints):
        # НОРМАЛИЗАЦИЯ ПО ПРАВОЙ ЧАСТИ:
        # если rhs < 0, умножаем всё на -1 и переворачиваем знак
        if rhs < -EPS:
            rhs = -rhs
            new_coeffs = {var: -coef for var, coef in coeffs.items()}
            coeffs = new_coeffs
            if sign == "<=":
                sign = ">="
            elif sign == ">=":
                sign = "<="
            # '=' не меняем

        row = [coeffs.get(v, 0.0) for v in orig_vars]
        added = []

        if sign == '<=':
            s = f"s{i+1}"
            added.append(s)
            slack_vars.append(s)
            row += [1.0]
        elif sign == '>=':
            s = f"s{i+1}"
            a = f"a{i+1}"
            added.extend([s, a])
            slack_vars.append(s)
            art_vars.append(a)
            row += [-1.0, 1.0]
        elif sign == '=':
            a = f"a{i+1}"
            added.append(a)
            art_vars.append(a)
            row += [1.0]

        A_rows.append(row)
        b.append(rhs)
        row_types.append(sign)
        added_per_row.append(added)

    added_all = []
    for added in added_per_row:
        for var in added:
            if var not in added_all:
                added_all.append(var)
    all_vars = orig_vars + added_all

    m = len(A_rows)
    n_orig = len(orig_vars)
    for i in range(m):
        current_len = len(A_rows[i])
        target_len = n_orig + len(added_all)
        if current_len < target_len:
            this_added = added_per_row[i]
            added_vals = {}
            for j, var in enumerate(this_added):
                added_vals[var] = A_rows[i][n_orig + j]
            full_added_part = [added_vals.get(var, 0.0) for var in added_all]
            A_rows[i] = A_rows[i][:n_orig] + full_added_part

    c = [float(obj.get(v, 0.0)) for v in orig_vars] + [0.0] * len(added_all)

    return all_vars, A_rows, b, c, slack_vars, art_vars, row_types

# ---------- Отображение таблицы ----------
def print_tableau(tableau, basic_vars, all_vars, phase, step):
    m = len(tableau) - 1
    header = ["Базис"] + all_vars + ["Свободный член"]
    print(f"\n=== {phase}: итерация {step} ===")
    print("# Текущее состояние симплекс-таблицы")
    print(" | ".join(f"{h:>10}" for h in header))
    print("-" * (11 * len(header)))
    for i in range(m):
        row = tableau[i]
        print(f"{basic_vars[i]:>8} | " + " | ".join(f"{row[j]:10.3f}" for j in range(len(row))))
    print("-" * (11 * len(header)))
    last = tableau[-1]
    print(f"{'F':>8} | " + " | ".join(f"{last[j]:10.3f}" for j in range(len(last))))
    print("=" * (11 * len(header)))

# ---------- Базовые операции симплекс-метода ----------
def pivot(tableau, row_idx, col_idx):
    piv = tableau[row_idx][col_idx]
    if abs(piv) < EPS:
        raise ValueError("Опорный элемент близок к нулю.")
    tableau[row_idx] = [v / piv for v in tableau[row_idx]]
    for i in range(len(tableau)):
        if i == row_idx:
            continue
        factor = tableau[i][col_idx]
        tableau[i] = [
            tableau[i][j] - factor * tableau[row_idx][j]
            for j in range(len(tableau[0]))
        ]

def find_entering_col(obj_row):
    candidates = [(j, val) for j, val in enumerate(obj_row[:-1])]
    min_val = min(val for j, val in candidates)
    if min_val >= -EPS:
        return None
    for j, val in candidates:
        if val == min_val:
            return j
    return None

def find_leaving_row(tableau, col):
    ratios = []
    for i, row in enumerate(tableau[:-1]):
        coeff = row[col]
        if coeff > EPS:
            ratios.append((i, row[-1] / coeff))
        else:
            ratios.append((i, float('inf')))
    min_ratio = min(r for i, r in ratios)
    if math.isinf(min_ratio):
        return None
    for i, r in ratios:
        if abs(r - min_ratio) < 1e-12:
            return i
    return None

def simplex_iterations(tableau, basic_vars, all_vars, phase_name):
    step = 0
    while True:
        print_tableau(tableau, basic_vars, all_vars, phase_name, step)
        enter = find_entering_col(tableau[-1])
        if enter is None:
            print("# Оптимальное решение достигнуто на данной фазе.")
            break
        leave = find_leaving_row(tableau, enter)
        if leave is None:
            raise ValueError(f"{phase_name}: неограниченная область допустимых решений.")
        print(f"# Переход: переменная {basic_vars[leave]} "
              f"заменяется на {all_vars[enter]} (строка {leave}, столбец {enter})")
        pivot(tableau, leave, enter)
        basic_vars[leave] = all_vars[enter]
        step += 1
    return tableau, basic_vars

# ---------- Фаза I ----------
def phase_one(all_vars, A, b, c, art_vars, row_types):
    print("\n--- ПОИСК ДОПУСТИМОГО БАЗИСА ---")
    m = len(A)
    n = len(all_vars)
    basic_vars = []

    for i, rt in enumerate(row_types):
        found = None
        for j, v in enumerate(all_vars):
            if v.startswith('s') and abs(A[i][j]) > EPS:
                found = v
                break
        if found:
            basic_vars.append(found)
        else:
            found_a = None
            for v in art_vars:
                j = all_vars.index(v)
                if abs(A[i][j]) > EPS:
                    found_a = v
                    break
            if found_a is None:
                for j, v in enumerate(all_vars):
                    col = [A[r][j] for r in range(m)]
                    if abs(A[i][j] - 1) < EPS and sum(
                        abs(x) for k, x in enumerate(col) if k != i
                    ) < EPS:
                        found_a = v
                        break
            if found_a is None:
                raise RuntimeError("Не удалось сформировать начальный базис.")
            basic_vars.append(found_a)

    tableau = []
    for i in range(m):
        tableau.append([A[i][j] for j in range(n)] + [b[i]])

    obj_aux = [1.0 if v in art_vars else 0.0 for v in all_vars] + [0.0]
    tableau.append(obj_aux)

    for i, bv in enumerate(basic_vars):
        if bv in art_vars:
            tableau[-1] = [
                tableau[-1][k] - tableau[i][k] for k in range(len(tableau[0]))
            ]

    print("\n# Запуск симплекс-итераций для Поиска допустимого базиса")
    tableau, basic_vars = simplex_iterations(
        tableau, basic_vars, all_vars, "Поиск допустимого базиса"
    )

    W = tableau[-1][-1]
    if abs(W) > 1e-6:
        raise ValueError(
            f"Поиск допустимого базиса: допустимого решения не существует (W* = {W})"
        )

    keep_idx = [j for j, v in enumerate(all_vars) if v not in art_vars]
    new_all_vars = [all_vars[j] for j in keep_idx]
    new_tableau = []
    for i in range(len(tableau)):
        new_row = [tableau[i][j] for j in keep_idx] + [tableau[i][-1]]
        new_tableau.append(new_row)

    new_basic = []
    for bv in basic_vars:
        new_basic.append(bv if bv in new_all_vars else new_all_vars[0])

    return new_all_vars, new_tableau, new_basic

# ---------- Фаза II ----------
def phase_two(all_vars, tableau, basic_vars, orig_c, goal):
    print("\n--- ПОИСК ОПТИМАЛЬНОГО РЕШЕНИЯ ---")
    n = len(all_vars)
    obj_row = [-orig_c[j] for j in range(n)] + [0.0]
    tableau[-1] = obj_row[:]

    cB = []
    for bv in basic_vars:
        if bv in all_vars:
            cB.append(orig_c[all_vars.index(bv)])
        else:
            cB.append(0.0)

    reduced = tableau[-1][:]
    for i, cb in enumerate(cB):
        if abs(cb) > EPS:
            reduced = [
                reduced[j] - cb * tableau[i][j] for j in range(len(reduced))
            ]
    tableau[-1] = reduced

    print("\n# Запуск симплекс-итераций для поиска оптимального решения")
    tableau, basic_vars = simplex_iterations(
        tableau, basic_vars, all_vars, "Поиск оптимального решения"
    )
    return tableau, basic_vars

# ---------- Извлечение результата ----------
def extract_solution(tableau, basic_vars, all_vars, free_map=None):
    if free_map is None:
        free_map = {}

    sol = {v: 0.0 for v in all_vars}
    m = len(tableau) - 1
    for i in range(m):
        sol[basic_vars[i]] = tableau[i][-1]
    Z = tableau[-1][-1]

    # Восстанавливаем исходные свободные переменные xk = xk_p - xk_m
    for orig, (vp, vm) in free_map.items():
        sol[orig] = sol.get(vp, 0.0) - sol.get(vm, 0.0)

    return sol, Z

# ---------- Основная программа ----------
def main():
    goal, obj, constraints, free_vars = read_lp("input.txt")

    print("=== ВХОДНЫЕ ДАННЫЕ ===")
    print("# Тип задачи (max/min):")
    print("Тип задачи:", goal)
    print("\n# Целевая функция:")
    print("Целевая функция:", obj)
    print("\n# Ограничения:")
    for c in constraints:
        print(" ", c)
    print("\n# Свободные переменные:", free_vars)

    # --- Обработка свободных переменных ---
    if free_vars:
        obj, constraints, free_map = transform_free_variables(obj, constraints, free_vars)
        print("\n# После разложения свободных переменных:")
        print("Новая целевая функция:", obj)
        print("Новые ограничения:")
        for c in constraints:
            print(" ", c)
    else:
        free_map = {}

    all_vars, A, b, c_vec, slacks, arts, row_types = build_canonical(goal, obj, constraints)
    print("\n=== ПОДГОТОВКА КАНОНИЧЕСКОЙ ФОРМЫ ===")
    print("# Список всех переменных:")
    print("Переменные:", all_vars)
    print("Искусственные переменные:", arts)

    if arts:
        all_vars_p, tableau_p, basic_p = phase_one(all_vars, A, b, c_vec, arts, row_types)
    else:
        print("\n# Искусственные переменные отсутствуют, базис формируется по slack-переменным.")
        basic_p = []
        for i in range(len(A)):
            found = None
            for j, v in enumerate(all_vars):
                if v.startswith('s') and abs(A[i][j]) > EPS:
                    found = v
                    break
            if found is None:
                for j, v in enumerate(all_vars):
                    col = [A[r][j] for r in range(len(A))]
                    if abs(A[i][j] - 1) < EPS and sum(
                        abs(col[k]) for k in range(len(col)) if k != i
                    ) < EPS:
                        found = v
                        break
            if found is None:
                found = f"b{i}"
            basic_p.append(found)
        tableau_p = [row[:] + [b_i] for row, b_i in zip(A, b)]
        tableau_p.append([0.0] * (len(all_vars)) + [0.0])
        all_vars_p = all_vars[:]

    orig_c = []
    for v in all_vars_p:
        if v in all_vars:
            orig_c.append(c_vec[all_vars.index(v)])
        else:
            orig_c.append(0.0)

    try:
        tableau_final, basic_final = phase_two(all_vars_p, tableau_p, basic_p, orig_c, goal)
    except ValueError as e:
        print("\n--- РЕЗУЛЬТАТ РЕШЕНИЯ ---")
        print(str(e))
        return

    sol_all, Z = extract_solution(tableau_final, basic_final, all_vars_p, free_map)

    print("\n--- РЕЗУЛЬТАТ РЕШЕНИЯ ---")
    print("# Оптимальные значения переменных:")

    # Выводим только исходные x1, x2, ... (без _p, _m)
    x_vars = sorted(
        [v for v in sol_all.keys()
         if v.startswith('x')
         and not v.endswith('_p')
         and not v.endswith('_m')],
        key=lambda name: int(re.findall(r'\d+', name)[0])
    )
    for v in x_vars:
        print(f"{v:>4} = {sol_all.get(v, 0.0):8.4f}")

    s_vars = [v for v in sol_all.keys() if v.startswith('s')]
    for v in sorted(s_vars, key=lambda name: int(re.findall(r'\d+', name)[0])):
        print(f"{v:>4} = {sol_all.get(v, 0.0):8.4f}")

    print(f"\n# Оптимальное значение целевой функции:")
    print(f"F* = {Z:8.4f}")

if __name__ == "__main__":
    main()
