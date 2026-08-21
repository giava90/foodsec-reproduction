"""
functions_supply.py
-----------------
Extensions on top of the original `functions.py`, which is imported
unchanged and NOT modified.

Purpose
=======
Add a "nutrient supplied off-farm" (e.g. fortification / supplementation)
counterfactual to the land-use optimization.

The original "without nutrient" analysis removed a nutrient by DROPPING its
column, which (a) changed the MAR denominator from 17 to 16 nutrients, so the
resulting curve was not point-for-point comparable to the full model, and
(b) carried the optimizer's infeasibility flags at low land.

Here, a "supplied" nutrient is handled differently:
  * Its adequacy is guaranteed (its NAR is set to 1) WITHOUT requiring any land
    to produce it -- i.e. the land-linked production constraint for that
    nutrient is removed, so the optimizer no longer allocates land to satisfy it.
  * It STAYS in the MAR denominator (still 17 nutrients), so the "supplied"
    curve and the full-model curve measure the same thing and can be overlaid.

This is the correct way to draw the "B12 supplied exogenously" curve: the black
(full) curve's bend should relax, rather than vanishing into an infeasible gap.

Everything not related to the supplied nutrient is left exactly as the original author wrote it.
"""

import numpy as np
from scipy.optimize import linprog

# Reuse the originals untouched.
from functions import from_land_to_intake, mar_f  # noqa: F401


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _resolve_supplied_idx(product_composition, supplied_nutrients):
    """Turn a list of nutrient names into column indices in product_composition.
    Names not present are ignored (with a printed note)."""
    idx = []
    cols = list(product_composition.columns)
    for name in supplied_nutrients:
        if name in cols:
            idx.append(cols.index(name))
        else:
            print(f"[functions_supply] supplied nutrient '{name}' not in columns; ignored.")
    return set(idx)


def mar_f_with_supply(x, y, supplied_idx=()):
    """Mean adequacy ratio where supplied nutrients count as fully adequate
    (ratio = 1) but REMAIN in the denominator, so the metric matches the
    full model (same number of nutrients)."""
    assert len(x) == len(y), "The two vectors must have the same length"
    ratio = x / y
    ratio[ratio > 1] = 1
    for i in supplied_idx:
        ratio[i] = 1.0
    return ratio.sum() / len(x)


# ---------------------------------------------------------------------------
# MAR vs land, with a nutrient supplied off-farm
# (mirror of optimize_land_usage_max_mar_fixing_caloric_intakes)
# ---------------------------------------------------------------------------
def optimize_land_usage_max_mar_fixing_caloric_intakes_with_supply(
        land, product_composition, dri,
        max_kcal="default",
        min_kcal="default",
        calorie_intake_per_macro={"prot": {"min": 0.15, "max": 0.35},
                                  "carb": {"min": 0.50, "max": 0.65}},
        calories_from_carb=4,
        calories_from_prot=4,
        yearly_ul_in_mg={},
        land_ul_for_animals=False,
        supplied_nutrients=()):
    """Same LP as the original optimize_land_usage_max_mar_fixing_caloric_intakes,
    except that for each nutrient in `supplied_nutrients` the land-linked
    adequacy constraint (a_i <= sum_j hat_c_ij x_j) is removed. Its a_i is then
    only bounded by a_i <= 1 and, because the objective maximizes sum(a_i),
    it reaches 1 for free -- no land is spent to satisfy it."""
    n_products = len(product_composition)
    n_nutrients = len(dri)
    supplied_idx = _resolve_supplied_idx(product_composition, supplied_nutrients)

    if type(max_kcal) is str and max_kcal == "default":
        max_kcal = dri[0] * 1.1
    if type(min_kcal) is str and min_kcal == "default":
        min_kcal = dri[0] * 0.9

    one_over_r = 1 / dri
    hat_c = product_composition.values * one_over_r
    hat_c = hat_c.astype(float)

    c = np.concatenate((np.zeros(n_products), -np.ones(n_nutrients)))
    A = []
    b = []
    for i in range(n_nutrients):
        # a_i <= 1  (kept for every nutrient, including supplied)
        constraint = np.zeros(n_nutrients)
        constraint[i] = 1
        constraint = np.concatenate((np.zeros(n_products), constraint))
        A.append(constraint)
        b.append(1)
        # a_i <= sum_j hat_c_ij x_j  (SKIPPED for supplied nutrients)
        if i in supplied_idx:
            continue
        constraint = np.zeros(n_nutrients)
        constraint[i] = 1
        coeff_i = hat_c[:, i]
        constraint = np.concatenate((-coeff_i, constraint))
        A.append(constraint)
        b.append(0)

    # Total land constraint
    A.append(np.concatenate((np.ones(n_products), np.zeros(n_nutrients))))
    b.append(land)

    # Caloric constraints
    if max_kcal:
        constraint = np.concatenate((product_composition["kcal"].values, np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(max_kcal)
    if min_kcal:
        constraint = np.concatenate((product_composition["kcal"].values, np.zeros(n_nutrients)))
        A.append(-constraint)
        b.append(-min_kcal)

    # Macro-nutrient caloric shares
    if "carb" in product_composition.columns and "carb" in calorie_intake_per_macro:
        constraint = np.concatenate(([calories_from_carb * g for g in product_composition["carb"].values],
                                     np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(calorie_intake_per_macro["carb"]["max"] * max_kcal)
        A.append(-constraint)
        b.append(-calorie_intake_per_macro["carb"]["min"] * min_kcal)
    if "prot" in product_composition.columns and "prot" in calorie_intake_per_macro:
        constraint = np.concatenate(([calories_from_prot * g for g in product_composition["prot"].values],
                                     np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(calorie_intake_per_macro["prot"]["max"] * max_kcal)
        A.append(-constraint)
        b.append(-calorie_intake_per_macro["prot"]["min"] * min_kcal)

    # Upper limits (act on production only; harmless to keep for supplied nutrients)
    for nutrient, value in yearly_ul_in_mg.items():
        if value is None:
            continue
        constraint = np.concatenate((product_composition[nutrient].values, np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(value)

    x_bounds = [(0, land) for _ in range(n_products)] + [(0, None) for _ in range(n_nutrients)]
    if land_ul_for_animals:
        for animal, ul in land_ul_for_animals.items():
            if animal in product_composition.index:
                x_bounds[product_composition.index.get_loc(animal)] = (0, ul)

    res = linprog(c, A_ub=A, b_ub=b, bounds=x_bounds, method="highs")
    return res


def compute_max_mar_fixing_caloric_intakes_with_supply(
        land_values, product_composition, annual_dietary_requirements,
        max_kcal="default", min_kcal="default",
        calorie_intake_per_macro={"prot": {"min": 0.15, "max": 0.35},
                                  "carb": {"min": 0.50, "max": 0.65}},
        calories_from_carb=4,
        calories_from_prot=4,
        yearly_ul_in_mg={},
        land_ul_for_animals=False,
        failure=-1,
        supplied_nutrients=()):
    """Mirror of the original compute_max_mar_fixing_caloric_intakes, but with a
    nutrient supplied off-farm. Returns the same [mar_value, res, control]
    structure per land value, so it plugs straight into the existing plotting
    code. MAR keeps the supplied nutrient in the denominator (set to 1)."""
    if failure not in [-1, "remove_constraints"]:
        raise ValueError("failure must be either -1 or 'remove_constraints'")

    supplied_idx = _resolve_supplied_idx(product_composition, supplied_nutrients)
    mar_values = []
    n_products = len(product_composition)

    for land in land_values:
        res = optimize_land_usage_max_mar_fixing_caloric_intakes_with_supply(
            land, product_composition, annual_dietary_requirements,
            max_kcal=max_kcal, min_kcal=min_kcal,
            calorie_intake_per_macro=calorie_intake_per_macro,
            calories_from_carb=calories_from_carb,
            calories_from_prot=calories_from_prot,
            yearly_ul_in_mg=yearly_ul_in_mg,
            land_ul_for_animals=land_ul_for_animals,
            supplied_nutrients=supplied_nutrients)

        if not res.success and failure == -1:
            mar_values.append([-0.1, res, -1])
            continue
        elif not res.success and failure == "remove_constraints":
            if type(min_kcal) == str and min_kcal == "default":
                min_kcal_local = annual_dietary_requirements[0] * 0.9
            else:
                min_kcal_local = min_kcal
            step = 0.01
            i = 0
            new_min_kcal = min_kcal_local
            while new_min_kcal is not False and new_min_kcal > 0 and not res.success:
                i += 1
                new_min_kcal = min_kcal_local - i * step * min_kcal_local
                res = optimize_land_usage_max_mar_fixing_caloric_intakes_with_supply(
                    land, product_composition, annual_dietary_requirements,
                    max_kcal="default", min_kcal=new_min_kcal,
                    calorie_intake_per_macro=calorie_intake_per_macro,
                    calories_from_carb=calories_from_carb,
                    calories_from_prot=calories_from_prot,
                    yearly_ul_in_mg=yearly_ul_in_mg,
                    land_ul_for_animals=land_ul_for_animals,
                    supplied_nutrients=supplied_nutrients)
            if not res.success:
                mar_values.append([-0.1, res, -1])
                continue
            else:
                control = 0
        else:
            control = 1

        intake = from_land_to_intake(res.x[:n_products], y=False,
                                     composition=product_composition, edible_portion=False)
        mar_value = mar_f_with_supply(intake, annual_dietary_requirements, supplied_idx)
        mar_values.append([mar_value, res, control])
    return mar_values


# ---------------------------------------------------------------------------
# Minimum land to full adequacy, with a nutrient supplied off-farm
# (mirror of min_land_usage + max_mar_min_land_fixing_caloric_intakes)
# ---------------------------------------------------------------------------
def min_land_usage_with_supply(
        product_composition, dri,
        max_kcal="default", min_kcal="default",
        calorie_intake_per_macro={"prot": {"min": 0.15, "max": 0.35},
                                  "carb": {"min": 0.50, "max": 0.65}},
        calories_from_carb=4,
        calories_from_prot=4,
        yearly_ul_in_mg={},
        land_ul_for_animals=False,
        land=np.inf,
        supplied_nutrients=()):
    """Same as the original min_land_usage, but the >=1 adequacy constraint is
    dropped for each supplied nutrient -- so the minimum land returned is the
    land needed to satisfy every OTHER nutrient, with the supplied one assumed
    covered off-farm."""
    n_products = len(product_composition)
    n_nutrients = len(dri)
    supplied_idx = _resolve_supplied_idx(product_composition, supplied_nutrients)

    if type(max_kcal) is str and max_kcal == "default":
        max_kcal = dri[0] * 1.1
    if type(min_kcal) is str and min_kcal == "default":
        min_kcal = dri[0] * 0.9

    c = np.ones(n_products)
    one_over_r = 1 / dri
    hat_c = product_composition.values * one_over_r
    hat_c = hat_c.astype(float)

    A = []
    b = []
    for i in range(n_nutrients):
        if i in supplied_idx:
            continue  # supplied off-farm -> no land required for this nutrient
        coeff_i = hat_c[:, i]
        A.append(-coeff_i)
        b.append(-1)

    if max_kcal:
        A.append(product_composition["kcal"].values)
        b.append(max_kcal)
    if min_kcal:
        A.append(-product_composition["kcal"].values)
        b.append(-min_kcal)

    if "carb" in product_composition.columns and "carb" in calorie_intake_per_macro:
        constraint = np.array([calories_from_carb * g for g in product_composition["carb"].values])
        A.append(constraint)
        b.append(calorie_intake_per_macro["carb"]["max"] * max_kcal)
        A.append(-constraint)
        b.append(-calorie_intake_per_macro["carb"]["min"] * min_kcal)
    if "prot" in product_composition.columns and "prot" in calorie_intake_per_macro:
        constraint = np.array([calories_from_prot * g for g in product_composition["prot"].values])
        A.append(constraint)
        b.append(calorie_intake_per_macro["prot"]["max"] * max_kcal)
        A.append(-constraint)
        b.append(-calorie_intake_per_macro["prot"]["min"] * min_kcal)

    for nutrient, value in yearly_ul_in_mg.items():
        if value is None:
            continue
        A.append(np.array(product_composition[nutrient].values))
        b.append(value)

    if land != np.inf:
        A.append(np.ones(n_products))
        b.append(land)

    x_bounds = [(0, land) for _ in range(n_products)]
    if land_ul_for_animals:
        for animal, ul in land_ul_for_animals.items():
            if animal in product_composition.index:
                x_bounds[product_composition.index.get_loc(animal)] = (0, ul)

    res = linprog(c, A_ub=A, b_ub=b, bounds=x_bounds, method="highs")
    return res


def max_mar_min_land_fixing_caloric_intakes_with_supply(
        product_composition, annual_dietary_requirements,
        max_kcal="default", min_kcal="default",
        calorie_intake_per_macro={"prot": {"min": 0.15, "max": 0.35},
                                  "carb": {"min": 0.50, "max": 0.65}},
        calories_from_carb=4,
        calories_from_prot=4,
        yearly_ul_in_mg={},
        land_ul_for_animals=False,
        failure=-1,
        land=np.inf,
        supplied_nutrients=()):
    """Mirror of the original max_mar_min_land_fixing_caloric_intakes, with a
    nutrient supplied off-farm. Returns [mar_value, res, control]; sum(res.x)
    is the minimum land to reach full adequacy of all OTHER nutrients."""
    if failure not in [-1, "remove_constraints"]:
        raise ValueError("failure must be either -1 or 'remove_constraints'")

    supplied_idx = _resolve_supplied_idx(product_composition, supplied_nutrients)
    n_products = len(product_composition)

    res = min_land_usage_with_supply(
        product_composition, annual_dietary_requirements,
        max_kcal=max_kcal, min_kcal=min_kcal,
        calorie_intake_per_macro=calorie_intake_per_macro,
        calories_from_carb=calories_from_carb,
        calories_from_prot=calories_from_prot,
        yearly_ul_in_mg=yearly_ul_in_mg,
        land_ul_for_animals=land_ul_for_animals,
        land=land, supplied_nutrients=supplied_nutrients)

    if not res.success and failure == -1:
        return [-0.1, res, -1]
    elif not res.success and failure == "remove_constraints":
        if type(min_kcal) == str and min_kcal == "default":
            min_kcal_local = annual_dietary_requirements[0] * 0.9
        else:
            min_kcal_local = min_kcal
        step = 0.01
        i = 0
        while min_kcal_local > 0 and not res.success:
            i += 1
            new_min_kcal = min_kcal_local - i * step * min_kcal_local
            res = min_land_usage_with_supply(
                product_composition, annual_dietary_requirements,
                max_kcal="default", min_kcal=new_min_kcal,
                calorie_intake_per_macro=calorie_intake_per_macro,
                calories_from_carb=calories_from_carb,
                calories_from_prot=calories_from_prot,
                yearly_ul_in_mg=yearly_ul_in_mg,
                land_ul_for_animals=land_ul_for_animals,
                land=land, supplied_nutrients=supplied_nutrients)
        if not res.success:
            return [-0.1, res, -1]
        control = 0
    else:
        control = 1

    intake = from_land_to_intake(res.x[:n_products], y=False,
                                 composition=product_composition, edible_portion=False)
    mar_value = mar_f_with_supply(intake, annual_dietary_requirements, supplied_idx)
    return [mar_value, res, control]
