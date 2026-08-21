import pdb

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from itertools import product
import random



def mar_f(x,y):
    """ 
    Function to calculate the mean adequacy ratio
    Before computing the index, we assert that the two vectors have the same length.
    Then, we compute the ratio of the two vectors, cap the ratio at 1, and return the mean of the capped ratios.
    x: a vector of intakes
    y: a vector of dietary reference intake 
    """
    assert len(x) == len(y), "The two vectors must have the same length"
    ratio = x/y
    ratio[ratio > 1] = 1
    return sum(ratio)/len(x)

def from_land_to_intake(land, y, composition, edible_portion=False):
    """
    Function to convert land usage to intake
    land: a vector of land usage
    y: a vector of yield
    edible_portion: a vector of edible fraction
    composition: a matrix of composition
    """
    if edible_portion is not False and y is not False:
        y = y*edible_portion
    else:
        y = np.ones(len(composition.index))
    return np.dot(composition.values.T, land*y)


def from_product_to_intake(product, quantity, composition):
    """
    Function to convert land usage to intake
    product: the product to consider in the diet
    quantity: the quantity of the product
    composition: a matrix of composition
    """
    nutrition = composition.loc[product]*quantity
    nutrition = nutrition.values
    return nutrition

def optimize_land_usage(land, product_composition, dri, y,e, max_kcal="default", min_kcal="default", penalty_lambda=1.0):
    """
    Function to optimize the land usage
    land: the total land available
    product_composition: a matrix of composition
    dri: a vector of dietary reference intake (the order of the nutrients in the matrix of composition must be the same)
    y: a vector of yield
    e: a vector of edible portion
    max_kcal: the maximum number of kcal. If "no_limit", the limit is set to infinity. If "default", the limit is set to 120% of the DRI
    min_kcal: the minimum number of kcal. If "no_limit", the limit is set to 0. If "default", the limit is set to 90% of the DRI
    penalty_lambda: the penalty for exceeding the DRI
    """
    n_products = len(product_composition)
    n_nutrients = len(dri)
    n = product_composition["kcal"].values
    if type(max_kcal) is str:
        if max_kcal == "no_limit":
            max_kcal = dri[0]*10
        else:
            max_kcal = dri[0]*1.2
    if type(min_kcal) is str:   
        if min_kcal == "no_limit":
            min_kcal = 0
        else:
            min_kcal =  dri[0]*0.9
    # Define the coefficients of the objective function
    one_over_r = 1/dri
    c_over_r = product_composition.values * one_over_r
    sum_c_over_r = c_over_r.sum(axis=1)
    c_basic = - (y*e)*sum_c_over_r
    # concatenate the values of c with the penalty values n_products times
    c = np.concatenate((c_basic, penalty_lambda*np.ones(n_nutrients)))
    
    # Define the coefficients of the inequality constraints
    A = [np.concatenate((y*e*n, np.zeros(n_nutrients))), 
        -np.concatenate((y*e*n, np.zeros(n_nutrients)))]
    b = [max_kcal, 
        -min_kcal]
    for i in range(n_nutrients):
        a_temp = np.zeros(n_nutrients)
        a_temp[i] = -1
        a_temp = np.concatenate(((y*e)*c_over_r[:,i], a_temp))
        A.append(a_temp)
        b.append(+1)


    # Define the coefficients of the equality constraints
    A_eq = [1 for i in range(n_products)]
    A_eq = np.concatenate((A_eq, np.zeros(n_nutrients)))
    A_eq = np.array([A_eq])
    b_eq = np.array([land])

    # Define the bounds of the variables
    x_bounds = [(0, land) for i in range(n_products)]
    x_bounds_penalty = [(0, None) for i in range(n_nutrients)]
    # Solve the problem
    res = linprog(c, A_ub=A, b_ub=b, A_eq=A_eq, b_eq=b_eq, bounds=x_bounds+x_bounds_penalty, method='highs')
    return res



def optimize_land_usage_max_mar(land, product_composition, dri, penalty_lambda=1.0, max_kcal="no_limit", min_kcal="no_limit", 
                                product_categories=None, recommended_caloric_intake_per_category=None, tolerance=0.05):
    """
    Function to optimize the land usage by maximizing the mean adequacy ratio while keeping the difference between the DRI small for each nutrient
    land: the total land available
    product_composition: a matrix of composition 
    dri: a vector of dietary reference intake (the order of the nutrients in the matrix of composition must be the same)
    penalty_lambda: the penalty for not having a similar ratio of each nutrient
    max_kcal: the maximum number of kcal. If "no_limit", the limit is set to infinity. If "default", the limit is set to 120% of the DRI
    min_kcal: the minimum number of kcal. If "no_limit", the limit is set to 0. If "default", the limit is set to 90% of the DRI
    product_categories: a dictionary {product:category}
    recommended_caloric_intake_per_category: a dictionary {category:caloric_intake_percentage}
    tolerance: the tolerance on the caloric intake recommended per category
    """
    n_products = len(product_composition)
    n_nutrients = len(dri)
    # Define the coefficients of the objective function
    one_over_r = 1/dri
    hat_c = product_composition.values * one_over_r
    hat_c = hat_c.astype(float)
    tilde_c = hat_c.sum(axis=1)

    # Coefficients for the objective function
    c = np.concatenate((-tilde_c, penalty_lambda*np.ones(n_nutrients)))

    # Define the coefficients of the inequality constraints
    A = []
    b = [] 
    #Then the constraints to have similar ratio of each nutrient
    av_hat_c = tilde_c / n_nutrients # average coefficient for each nutrient
    for i in range(n_nutrients):
        # first constraint
        coeff_i =  hat_c[:,i] - av_hat_c
        constraint = np.zeros(n_nutrients)
        constraint[i] = -1
        constraint = np.concatenate((coeff_i, constraint))
        A.append(constraint)
        b.append(0)
        # second constraint
        constraint = np.zeros(n_nutrients)
        constraint[i] = -1
        constraint = np.concatenate((-coeff_i, constraint))
        A.append(constraint)
        b.append(0)

    # Total land constraint
    A.append(np.concatenate((np.ones(n_products), np.zeros(n_nutrients))))
    b.append(land)

    # Caloric constraint
    if max_kcal != "no_limit":
        if max_kcal == "default":
            max_kcal = dri[0]*1.2
        constraint = np.concatenate((product_composition["kcal"].values, np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(max_kcal)
    if min_kcal != "no_limit":
        if min_kcal == "default":
            min_kcal =  dri[0]*0.9
        constraint = np.concatenate((product_composition["kcal"].values, np.zeros(n_nutrients)))
        A.append(-constraint)
        b.append(-min_kcal)

    #Lancet caloric recommendation breakdown by category
    # we add a tolerance of 5% on the caloric intake
    if product_categories is not None and recommended_caloric_intake_per_category is not None:
        for category, caloric_intake in recommended_caloric_intake_per_category.items():
            # get the products in the category
            products_in_category = [product for product in product_categories if product_categories[product] == category]
            # sort the products kcal in the category by product_composition.index
            coeff_i = [values["kcal"] if product in products_in_category else 0 for product, values in product_composition.iterrows()]
            if sum(coeff_i) == 0:
                print("No product in the category", category)
                print("Removing the constraint for this category")
                continue
            if len([product  for product, values in product_composition.iterrows() if product in products_in_category]) <=1:
                print("Only one product in the category", category)
                print("Removing the constraint for this category")
                continue
            print("Category", category, "caloric intake", caloric_intake)
            print([product  for product, values in product_composition.iterrows() if product in products_in_category])
            print("----------------")
            constraint = coeff_i - product_composition["kcal"].values * caloric_intake * (1 + tolerance)
            constraint = np.concatenate((constraint, np.zeros(n_nutrients)))
            A.append(constraint)
            b.append(0)
            constraint = coeff_i - product_composition["kcal"].values * caloric_intake * (1 - tolerance) 
            constraint = np.concatenate((constraint, np.zeros(n_nutrients)))
            A.append( - constraint)
            b.append(0)
    
    # Define the bounds of the variables
    x_bounds = [(0, land) for i in range(n_products)] + [(0, None) for i in range(n_nutrients)]
    # Solve the problem
    res = linprog(c, A_ub=A, b_ub=b, bounds=x_bounds, method='highs')

    return res


def minimize_objective(c, a, L):
    n = len(c)  # number of x variables
    m = len(a[0])  # number of k terms in absolute value sums
    
    # Coefficients for the objective function: -c_i for x_i and 1 for u_k
    c_obj = [-ci for ci in c] + [1] * m
    
    # Constraints matrix and vector
    A = []
    b = []
    
    # Constraints for u_k >= sum(x_i * a_ik) and u_k >= -sum(x_i * a_ik)
    for k in range(m):
        constraint1 = [a[i][k] for i in range(n)] + [-1 if j == k else 0 for j in range(m)]
        constraint2 = [-a[i][k] for i in range(n)] + [-1 if j == k else 0 for j in range(m)]
        A.append(constraint1)
        A.append(constraint2)
        b.append(0)
        b.append(0)
    
    # Constraint for sum(x_i) <= L
    A.append([1] * n + [0] * m)
    b.append(L)
    
    # Bounds for variables x_i (no explicit bounds) and u_k (non-negative)
    bounds = [(None, None)] * n + [(0, None)] * m
    
    # Use linprog to solve the problem
    result = linprog(c_obj, A_ub=A, b_ub=b, bounds=bounds, method='highs')
    
    # Extract the optimized values
    x_opt = result.x[:n]
    u_opt = result.x[n:]
    
    return x_opt, u_opt, result.fun


def compute_max_mar(land_values, product_composition, annual_dietary_requirements, penalty_lambda=1.0):
    """
    This function computes the maximum mean adequacy ratio for a list of land values

    Parameters
    ----------
    land_values : list
        The list of land values
    product_composition : pandas.DataFrame
        The composition of the products
    annual_dietary_requirements : list 
        annual_dietary_requirements
    penalty_lambda : float
    Returns
    -------
    list
        The maximum mean adequacy ratio for each land value
    """
    mar_values = []
    n_products = len(product_composition)
    for land in land_values:
        res = optimize_land_usage_max_mar(land, product_composition, annual_dietary_requirements, penalty_lambda=penalty_lambda)
        if not res.success:
            mar_values.append(-0.1)
            continue
        annual_dietary_intake_from_solution = from_land_to_intake(res.x[:n_products], y=False, composition=product_composition, edible_portion= False)
        mar_value = mar_f(annual_dietary_intake_from_solution, annual_dietary_requirements)
        mar_values.append(mar_value)
    return mar_values


def compute_max_mar_basic(land_values, product_composition, annual_dietary_requirements, max_kcal="default", min_kcal="default"):
    """
    This function computes the maximum mean adequacy ratio for a list of land values

    Parameters
    ----------
    land_values : list
        The list of land values
    product_composition : pandas.DataFrame
        The composition of the products
    annual_dietary_requirements : list 
        annual_dietary_requirements
    Returns
    -------
    list
        The maximum mean adequacy ratio for each land value
    """
    mar_values = []
    n_products = len(product_composition)
    for land in land_values:
        res = optimize_land_usage_max_mar_basic(land, product_composition, annual_dietary_requirements, max_kcal, min_kcal)
        if not res.success:
            mar_values.append(-0.1)
            continue
        annual_dietary_intake_from_solution = from_land_to_intake(res.x[:n_products], y=False, composition=product_composition, edible_portion= False)
        mar_value = mar_f(annual_dietary_intake_from_solution, annual_dietary_requirements)
        mar_values.append(mar_value)
    return mar_values

def optimize_land_usage_max_mar_basic(land, product_composition, dri, max_kcal="default", min_kcal="default"):
    """
    Function to optimize the land usage by maximizing the mean adequacy ratio while keeping the difference between the DRI small for each nutrient
    land: the total land available
    product_composition: a matrix of composition 
    dri: a vector of dietary reference intake (the order of the nutrients in the matrix of composition must be the same)
    max_kcal: the maximum number of kcal. If "no_limit", the limit is set to infinity. If "default", the limit is set to 120% of the DRI
    min_kcal: the minimum number of kcal. If "no_limit", the limit is set to 0. If "default", the limit is set to 90% of the DRI
    Returns
    -------
    res: the result of the optimization
    """
    n_products = len(product_composition)
    n_nutrients = len(dri)
    if type(max_kcal) is str:
        if max_kcal == "default":
            max_kcal = dri[0]*1.1
    if type(min_kcal) is str:   
        if min_kcal == "default":
            min_kcal =  dri[0]*0.9
    # Define the coefficients of the objective function
    one_over_r = 1/dri
    hat_c = product_composition.values * one_over_r
    hat_c = hat_c.astype(float)

    c = np.concatenate((np.zeros(n_products), -np.ones(n_nutrients)))
    # Define the coefficients of the inequality constraints
    A = []
    b = []
    for i in range(n_nutrients):
        # first constraint a_i <= 1
        constraint = np.zeros(n_nutrients)
        constraint[i] = 1
        constraint = np.concatenate((np.zeros(n_products), constraint))
        A.append(constraint)
        b.append(1)
        # second constraint a_i <= sum_j hat_c_ij
        constraint = np.zeros(n_nutrients)
        constraint[i] = 1
        coeff_i = hat_c[:,i]
        constraint = np.concatenate((-coeff_i, constraint))
        A.append(constraint)
        b.append(0)
    # Total land constraint
    A.append(np.concatenate((np.ones(n_products), np.zeros(n_nutrients))))
    b.append(land)
    # Caloric constraint
    if max_kcal:
        constraint = np.concatenate((product_composition["kcal"].values, np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(max_kcal)
    # Define the bounds of the variables
    x_bounds = [(0, land) for i in range(n_products)] + [(0, None) for i in range(n_nutrients)]
    # Solve the problem
    res = linprog(c, A_ub=A, b_ub=b, bounds=x_bounds, method='highs')

    return res



def compute_max_mar_fixing_caloric_intakes(land_values, product_composition, annual_dietary_requirements, 
                                           max_kcal="default", min_kcal="default",
                                           calorie_intake_per_macro={"prot": {"min": 0.15, "max" : 0.35}, # we add a 5% to the lower bound to not account for fats
                                                                            "carb": {"min": 0.50, "max" : 0.65}},# we add a 5% to the lower bound to not account for fats
                                           calories_from_carb = 4,
                                           calories_from_prot = 4,
                                           yearly_ul_in_mg = {},
                                           land_ul_for_animals = False,
                                           failure = -1):
    """
    This function computes the maximum mean adequacy ratio for a list of land values

    Parameters
    ----------
    land_values : list
        The list of land values
    product_composition : pandas.DataFrame
        The composition of the products
    annual_dietary_requirements : list 
        annual_dietary_requirements
    max_kcal: float
        the maximum number of kcal. If "no_limit", the limit is set to infinity. If "default", the limit is set to 120% of the DRI
    min_kcal: float
        the minimum number of kcal. If "no_limit", the limit is set to 0. If "default", the limit is set to 90% of the DRI
    calorie_intake_per_macro: dict
        a dictionary with the minimum and maximum caloric intake per macro-nutrient
    calories_from_carb: float
        the number of calories per gram of carbohydrates
    calories_from_prot: float
        the number of calories per gram of proteins
    yearly_ul_in_mg: dict
        a dictionary with the yearly upper limits for each nutrient
    land_ul_for_animals: dict
        the upper limit for the land usage for each animal product, if any.
        If False, no upper limit is set.
    failure: int or str
        if -1, the function will return -0.1 if the optimization fails
        if "remove_constraints", the function will try to remove the constraints and optimize again if the optimization fails
    Returns
    -------
    list
        The maximum mean adequacy ratio for each land value
    """
    # check that failure is either -1 or "remove_constraints"
    if failure not in [-1, "remove_constraints"]:
        raise ValueError("failure must be either -1 or 'remove_constraints'")
    mar_values = []
    n_products = len(product_composition)
    for land in land_values:
        res = optimize_land_usage_max_mar_fixing_caloric_intakes(land, product_composition, annual_dietary_requirements, 
                                                                 max_kcal=max_kcal, 
                                                                 min_kcal=min_kcal,
                                                                 calorie_intake_per_macro=calorie_intake_per_macro,
                                                                 calories_from_carb = calories_from_carb,
                                                                 calories_from_prot = calories_from_prot,
                                                                 yearly_ul_in_mg = yearly_ul_in_mg,
                                                                 land_ul_for_animals = land_ul_for_animals)
        if not res.success and failure == -1:
            control = -1
            mar_values.append([-0.1, res,control])
            continue
        elif not res.success and failure == "remove_constraints":
            print(res.message)
            print("We might have too many constraints")
            if type(min_kcal) == str:
                if min_kcal == "default":
                    min_kcal =  annual_dietary_requirements[0]*0.9
            step = 0.01
            i=0
            new_min_kcal = min_kcal
            while new_min_kcal > 0 and not res.success:
                i+=1
                new_min_kcal = min_kcal - i*step*min_kcal
                res = optimize_land_usage_max_mar_fixing_caloric_intakes(land, 
                                                                    product_composition, 
                                                                    annual_dietary_requirements, 
                                                                    max_kcal="default",
                                                                    min_kcal=new_min_kcal,
                                                                    calorie_intake_per_macro=calorie_intake_per_macro,
                                                                    calories_from_carb = calories_from_carb,
                                                                    calories_from_prot = calories_from_prot,
                                                                    yearly_ul_in_mg = yearly_ul_in_mg,
                                                                    land_ul_for_animals = land_ul_for_animals
                                                                    )
                
            if not res.success:
                print("... STILL not successful")
                control = -1
                mar_values.append([-0.1, res, control])
                continue
            else:
                control = 0
        else:
            control = 1
        annual_dietary_intake_from_solution = from_land_to_intake(res.x[:n_products], y=False, composition=product_composition, edible_portion= False)
        mar_value = mar_f(annual_dietary_intake_from_solution, annual_dietary_requirements)
        mar_values.append([mar_value, res, control])
    return mar_values



def max_mar_min_land_fixing_caloric_intakes(product_composition, annual_dietary_requirements, 
                                           max_kcal="default", min_kcal="default",
                                           calorie_intake_per_macro={"prot": {"min": 0.15, "max" : 0.35}, # we add a 5% to the lower bound to not account for fats
                                                                            "carb": {"min": 0.50, "max" : 0.65}},# we add a 5% to the lower bound to not account for fats
                                           calories_from_carb = 4,
                                           calories_from_prot = 4,
                                           yearly_ul_in_mg = {},
                                           land_ul_for_animals = False,
                                           failure = -1,
                                           land = np.inf):
    """
    This function computes the minimum land needed to get the MAR =1.

    Parameters
    ----------
    product_composition : pandas.DataFrame
        The composition of the products
    annual_dietary_requirements : list 
        annual_dietary_requirements
    max_kcal: float
        the maximum number of kcal. If "no_limit", the limit is set to infinity. If "default", the limit is set to 120% of the DRI
    min_kcal: float
        the minimum number of kcal. If "no_limit", the limit is set to 0. If "default", the limit is set to 90% of the DRI
    calorie_intake_per_macro: dict
        a dictionary with the minimum and maximum caloric intake per macro-nutrient
    calories_from_carb: float
        the number of calories per gram of carbohydrates
    calories_from_prot: float
        the number of calories per gram of proteins
    yearly_ul_in_mg: dict
        a dictionary with the yearly upper limits for each nutrient
    land_ul_for_animals: dict
        the upper limit for the land usage for each animal product, if any.
        If False, no upper limit is set.
    failure: int or str
        if -1, the function will return -0.1 if the optimization fails
        if "remove_constraints", the function will try to remove the constraints and optimize again if the optimization fails
    land: float
        the maximum land that can be used
    Returns
    -------
    list
        The maximum mean adequacy ratio for each land value
    """
    # check that failure is either -1 or "remove_constraints"
    if failure not in [-1, "remove_constraints"]:
        raise ValueError("failure must be either -1 or 'remove_constraints'")
    n_products = len(product_composition)
    res = min_land_usage(product_composition, annual_dietary_requirements, 
                        max_kcal=max_kcal, 
                        min_kcal=min_kcal,
                        calorie_intake_per_macro=calorie_intake_per_macro,
                        calories_from_carb = calories_from_carb,
                        calories_from_prot = calories_from_prot,
                        yearly_ul_in_mg = yearly_ul_in_mg,
                        land_ul_for_animals = land_ul_for_animals,
                        land = land)
    if not res.success and failure == -1:
        control = -1
        print(res.message)
        return [-0.1, res,control]
    elif not res.success and failure == "remove_constraints":
        print(res.message)
        print("We might have too many constraints")
        if type(min_kcal) == str:
            if min_kcal == "default":
                min_kcal =  annual_dietary_requirements[0]*0.9
        step = 0.01
        i=0
        while min_kcal > 0 and not res.success:
            i+=1
            new_min_kcal = min_kcal - i*step*min_kcal
            res = min_land_usage(product_composition, 
                                annual_dietary_requirements, 
                                max_kcal="default",
                                min_kcal=new_min_kcal,
                                calorie_intake_per_macro=calorie_intake_per_macro,
                                calories_from_carb = calories_from_carb,
                                calories_from_prot = calories_from_prot,
                                yearly_ul_in_mg = yearly_ul_in_mg,
                                land_ul_for_animals = land_ul_for_animals,
                                land = land
                                )
        
        if not res.success:
            print("... STILL not successful")
            control = -1
            return [-0.1, res, control]
        else:
            control = 0
    else:
        control = 1
    annual_dietary_intake_from_solution = from_land_to_intake(res.x[:n_products], y=False, composition=product_composition, edible_portion= False)
    mar_value = mar_f(annual_dietary_intake_from_solution, annual_dietary_requirements)
    return [mar_value, res, control]

def optimize_land_usage_max_mar_fixing_caloric_intakes(land, product_composition, dri, 
                                                       max_kcal="default", # a maximum caloric intake of 110% of the DRI
                                                       min_kcal="default", 
                                                       calorie_intake_per_macro={"prot": {"min": 0.15, "max" : 0.35}, # we add a 5% to the lower bound to not account for fats
                                                                                 "carb": {"min": 0.50, "max" : 0.65}},# we add a 5% to the lower bound to not account for fats
                                                       calories_from_carb = 4,
                                                       calories_from_prot = 4,
                                                       yearly_ul_in_mg = {},
                                                       land_ul_for_animals = False): 
    """
    Function to optimize the land usage by maximizing the mean adequacy ratio while keeping the difference between the DRI small for each nutrient
    land: the total land available
    product_composition: a matrix of composition 
    dri: a vector of dietary reference intake (the order of the nutrients in the matrix of composition must be the same)
    max_kcal: the maximum number of kcal. If False, the limit is set to infinity. If "default", the limit is set to 120% of the DRI
    min_kcal: the minimum number of kcal. If False, the limit is set to 0. If "default", the limit is set to 90% of the DRI
    calorie_intake_per_macro: a dictionary with the minimum and maximum caloric intake per macro-nutrient
    calories_from_carb: the number of calories per gram of carbohydrates
    calories_from_prot: the number of calories per gram of proteins
    yearly_ul_in_mg: a dictionary with the yearly upper limits for each nutrient
    land_ul_for_animals: the upper limit for the land usage for each animal product, if any.
    If None, no upper limit is set.
    Returns
    -------
    res: the result of the optimization
    """
    n_products = len(product_composition)
    n_nutrients = len(dri)
    if type(max_kcal) is str:
        if max_kcal == "default":
            max_kcal = dri[0]*1.1
    if type(min_kcal) is str:   
        if min_kcal == "default":
            min_kcal =  dri[0]*0.9
    # Define the coefficients of the objective function
    one_over_r = 1/dri
    hat_c = product_composition.values * one_over_r
    hat_c = hat_c.astype(float)

    c = np.concatenate((np.zeros(n_products), -np.ones(n_nutrients)))
    # Define the coefficients of the inequality constraints
    A = []
    b = []
    for i in range(n_nutrients):
        # first constraint a_i <= 1
        constraint = np.zeros(n_nutrients)
        constraint[i] = 1
        constraint = np.concatenate((np.zeros(n_products), constraint))
        A.append(constraint)
        b.append(1)
        # second constraint a_i <= sum_j hat_c_ij
        constraint = np.zeros(n_nutrients)
        constraint[i] = 1
        coeff_i = hat_c[:,i]
        constraint = np.concatenate((-coeff_i, constraint))
        A.append(constraint)
        b.append(0)
    # Total land constraint
    A.append(np.concatenate((np.ones(n_products), np.zeros(n_nutrients))))
    b.append(land)
    # Caloric constraint
    # The maximum caloric intake should not exceed max_kcal
    if max_kcal:
        constraint = np.concatenate((product_composition["kcal"].values, np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(max_kcal)
    # The minimum caloric intake should be at least min_kcal 
    if min_kcal:
        constraint = np.concatenate((product_composition["kcal"].values, np.zeros(n_nutrients)))
        A.append(-constraint)  
        b.append(-min_kcal)
    #The maximum and minimum caloric intake per macro-nutrient should be respected
    if "carb" in product_composition.columns and "carb" in calorie_intake_per_macro:
        constraint = np.concatenate(([calories_from_carb*grams_of_carbs for grams_of_carbs in product_composition["carb"].values], np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(calorie_intake_per_macro["carb"]["max"]*max_kcal)
        A.append(-constraint)
        b.append(-calorie_intake_per_macro["carb"]["min"]*min_kcal)
    if "prot" in product_composition.columns and "prot" in calorie_intake_per_macro:
        constraint = np.concatenate(([calories_from_prot*grams_of_prot for grams_of_prot in product_composition["prot"].values], np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(calorie_intake_per_macro["prot"]["max"]*max_kcal)
        A.append(-constraint)
        b.append(-calorie_intake_per_macro["prot"]["min"]*min_kcal)
    # Define the maxim intake for each nutrient
    for nutrient, value in yearly_ul_in_mg.items():
        if value == None: # no upper limit
            continue
        constraint = np.concatenate((product_composition[nutrient].values, np.zeros(n_nutrients)))
        A.append(constraint)
        b.append(value)
    #Define the bounds of the variables
    x_bounds = [(0, land) for i in range(n_products)] + [(0, None) for i in range(n_nutrients)]
    # If there is an upper limit for the land usage for eggs, we set the upper limit for the land usage for eggs
    if land_ul_for_animals:
        for animal, ul in land_ul_for_animals.items():
            if animal in product_composition.index:
                x_bounds[product_composition.index.get_loc(animal)] = (0, ul)                

    # Solve the problem
    res = linprog(c, A_ub=A, b_ub=b, bounds=x_bounds, method='highs')

    return res



def min_land_usage(product_composition, dri, 
                    max_kcal="default", # a maximum caloric intake of 110% of the DRI
                    min_kcal="default", 
                    calorie_intake_per_macro={"prot": {"min": 0.15, "max" : 0.35}, # we add a 5% to the lower bound to not account for fats
                                                "carb": {"min": 0.50, "max" : 0.65}},# we add a 5% to the lower bound to not account for fats
                    calories_from_carb = 4,
                    calories_from_prot = 4,
                    yearly_ul_in_mg = {},
                    land_ul_for_animals = False,
                    land = np.inf): 
    """
    Function to optimize the land usage by maximizing the mean adequacy ratio while keeping the difference between the DRI small for each nutrient
    product_composition: a matrix of composition 
    dri: a vector of dietary reference intake (the order of the nutrients in the matrix of composition must be the same)
    max_kcal: the maximum number of kcal. If False, the limit is set to infinity. If "default", the limit is set to 120% of the DRI
    min_kcal: the minimum number of kcal. If False, the limit is set to 0. If "default", the limit is set to 90% of the DRI
    calorie_intake_per_macro: a dictionary with the minimum and maximum caloric intake per macro-nutrient
    calories_from_carb: the number of calories per gram of carbohydrates
    calories_from_prot: the number of calories per gram of proteins
    yearly_ul_in_mg: a dictionary with the yearly upper limits for each nutrient
    land_ul_for_animals: the upper limit for the land usage for each animal product, if any.
    If None, no upper limit is set.
    land: the total land available

    Returns
    -------
    res: the result of the optimization
    """
    n_products = len(product_composition)
    n_nutrients = len(dri)
    if type(max_kcal) is str:
        if max_kcal == "default":
            max_kcal = dri[0]*1.1
    if type(min_kcal) is str:   
        if min_kcal == "default":
            min_kcal =  dri[0]*0.9
    # Coefficients for the objective function
    c = np.ones(n_products)
    # Coefficients of the inequality constraints
    one_over_r = 1/dri
    hat_c = product_composition.values * one_over_r
    hat_c = hat_c.astype(float)
    A = []
    b = []
    for i in range(n_nutrients):
        # constraint sum_j hat_c_ij x_j >=1
        constraint = np.ones(n_products)
        coeff_i = hat_c[:,i]
        A.append(-coeff_i)
        b.append(-1)
    # Caloric constraint
    # The maximum caloric intake should not exceed max_kcal
    if max_kcal:
        constraint = product_composition["kcal"].values
        A.append(constraint)
        b.append(max_kcal)
    # The minimum caloric intake should be at least min_kcal 
    if min_kcal:
        constraint = product_composition["kcal"].values
        A.append(-constraint)  
        b.append(-min_kcal)
    #The maximum and minimum caloric intake per macro-nutrient should be respected
    if "carb" in product_composition.columns and "carb" in calorie_intake_per_macro:
        constraint = np.array([calories_from_carb*grams_of_carbs for grams_of_carbs in product_composition["carb"].values])
        A.append(constraint)
        b.append(calorie_intake_per_macro["carb"]["max"]*max_kcal)
        A.append(-constraint)
        b.append(-calorie_intake_per_macro["carb"]["min"]*min_kcal)
    if "prot" in product_composition.columns and "prot" in calorie_intake_per_macro:
        constraint = np.array([calories_from_prot*grams_of_prot for grams_of_prot in product_composition["prot"].values])
        A.append(constraint)
        b.append(calorie_intake_per_macro["prot"]["max"]*max_kcal)
        A.append(-constraint)
        b.append(-calorie_intake_per_macro["prot"]["min"]*min_kcal)
    # Define the maxim intake for each nutrient
    for nutrient, value in yearly_ul_in_mg.items():
        if value == None: # no upper limit
            continue
        constraint = np.array(product_composition[nutrient].values)
        A.append(constraint)
        b.append(value)
    # Total land constraint
    if land != np.inf:
        A.append(np.ones(n_products))
        b.append(land)
    #Define the bounds of the variables
    x_bounds = [(0, land) for i in range(n_products)]
    # If there is an upper limit for the land usage for eggs, we set the upper limit for the land usage for eggs
    if land_ul_for_animals:
        for animal, ul in land_ul_for_animals.items():
            if animal in product_composition.index:
                x_bounds[product_composition.index.get_loc(animal)] = (0, ul)                

    # Solve the problem
    res = linprog(c, A_ub=A, b_ub=b, bounds=x_bounds, method='highs')

    return res

def creating_product_composition_table(country, df_crop_yield, 
                                       df_composition, 
                                       df_yield_nutrients_per_animal_per_ha, 
                                       country_to_macro_area, 
                                       macro_area_to_lsu_per_ha, 
                                       baseline_lsu_per_ha=1.0, 
                                       verbose=False):
    """
    This function creates the product composition table for a specific country.
    It takes the country name, the crop yield dataframe, the composition dataframe, the animal composition dataframe,
    the country to macro area dataframe and the macro area to lsu per ha dataframe as input.
    It returns the product composition table for the country. 
    First it creates the crop yield dataframe.
    Then, it create the animal composition dataframe.
    Finally, it merges the two dataframes to create the product composition table.
    The function also checks if the products are available in the country and if not, it drops them from the dataframe.
    All animals are considered to be available in the country.

    Parameters:
    country (str): The name of the country.
    df_crop_yield (pd.DataFrame): The crop yield dataframe.
    df_composition (pd.DataFrame): The composition dataframe.
    df_yield_nutrients_per_animal_per_ha (pd.DataFrame): The animal composition dataframe.
    country_to_macro_area (dict): The country to macro area dict.
    macro_area_to_lsu_per_ha (dict): The macro area to lsu per ha dict.
    baseline_lsu_per_ha (float): The baseline lsu per ha assumed when the df_yield_nutrients_per_animal_per_ha was created. Default is 1.0.
    verbose (bool): If True, the function will print the products that are not available in the country. Default is False.

    Returns:
    df_products_sg (pd.DataFrame): The product composition table for the country. If not enough products are available, it returns None.
    """
    # consider the yield for a specific country
    df_crop_yield_sg = df_crop_yield[df_crop_yield["country"]==country]
    # drop the rows that have NaN in the "mean" column
    df_crop_yield_sg = df_crop_yield_sg.dropna(subset=["mean"])
    # consider the composition of the products that are available in the country
    av_products = df_crop_yield_sg["crop"].values
    # get the intersection of the available products and the products in the composition table
    products = set(av_products).intersection(set(df_composition.index))
    av_products = df_crop_yield["crop"].values
    products = list(products)
    if len(products) <= 1:
        print("Not enough products to optimize the land usage")
        return None
    if len(products) != len(av_products):
        # drop the products that are not available in the country
        df_crop_yield_sg = df_crop_yield_sg[df_crop_yield_sg["crop"].isin(products)]
        if verbose:
            print(country, "not all products are available")
            print("The products that are not available are", set(av_products).difference(set(products)))        
    # order the dataframe of composition and yield by the products
    df_composition_sg = df_composition.loc[products]
    df_crop_yield_sg = df_crop_yield_sg.set_index("crop")
    df_crop_yield_sg = df_crop_yield_sg.loc[products]
    # extract the composition of the edible portion
    e = df_composition_sg["edible.prop"].values #edible_portion
    # drop the fat and edible.prop columns as this will not be used in the optimization
    df_composition_sg = df_composition_sg.drop(columns=["fat", "edible.prop"])
    # dataframe with edible composition per ha per 1kg of crop
    df_effective_composition_per_ha_sg = df_composition_sg * (df_crop_yield_sg["mean"].values * e)[:, np.newaxis]

    # create the animal composition dataframe
    # get the macro area for the country
    lsu_per_ha = country_to_macro_area[country]
    # get the lsu per ha for the macro area
    lsu_per_ha = macro_area_to_lsu_per_ha[lsu_per_ha]
    # update the yield of the animals according to the macro area
    df_yield_nutrients_per_animal_per_ha = df_yield_nutrients_per_animal_per_ha.copy()
    df_yield_nutrients_per_animal_per_ha = df_yield_nutrients_per_animal_per_ha * (lsu_per_ha / baseline_lsu_per_ha)
    # append the df_composition_sg to the df_yield_nutrients_per_animal_per_ha dataframe
    df_products_sg = pd.concat([df_yield_nutrients_per_animal_per_ha, df_effective_composition_per_ha_sg])
    # preparing input data for the optimization
    return df_products_sg

def test_fake_yield_df(df_crop_yield, countries, fname, save=False):
    # create a fake table with the columns product, year, country, regione, yield_per_ha
    # use the product for the product, year from 2015 to 2025, country from countries, region use country and a dashed line with _1 to _4, yield_per_ha random values between 0.5 and 5
    data = []
    years = range(2015, 2026)
    product_set = set()
    for product in list(df_crop_yield["crop"].unique()):
        if product in product_set:
            continue
        if product == None:
            continue
        product_set.add(product)
        p_values = df_crop_yield[df_crop_yield["crop"]==product]["mean"].values
        if len(p_values) == 0:
            continue
        mu= p_values.mean()
        if len(p_values) ==1:
            sigma = mu
        else:
            sigma = p_values.std()
        for country in countries:
            for year in years:
                for region_suffix in range(1, 5):
                    region = f"{country}_{region_suffix}"
                    yield_per_ha = round(random.normalvariate(mu, sigma), 2)
                    yield_per_ha = max(0, yield_per_ha)
                    data.append({
                        "crop": product,
                        "year": year,
                        "country": country,
                        "region": region,
                        "mean": yield_per_ha
                    })
    df_fake_yield = pd.DataFrame(data)
    # save the dataframe to a csv file
    if save:
        df_fake_yield.to_csv("../data/processed/"+fname+".csv", index=False)
    return df_fake_yield


def preprocess_region_yield(df_yield_region):
    # Group by all keys except season_name to make an average of the yield values across season_name across each crop, year, country and region
    df_clean = (
        df_yield_region
        .groupby(["crop", "country", "region", "year"], as_index=False)
        .agg({"yield": "mean"})
    )
    # note that the yield from hv are in different units compared to FAO: we need to multiply by 1000
    df_clean["yield"] *=1000
    return df_clean


def create_full_dataframe(crops, years, countries, country_to_region):
    rows = []

    for country in countries:
        regions = country_to_region.get(country, [country])
        for crop, year, region in product(crops, years, regions):
            rows.append((crop, country, region, year))

    df = pd.DataFrame(rows, columns=["crop", "country", "region", "year"])
    return df

def fill_yield(df_full, df_yield_region, df_yield_country):
    # Merge region-level yields
    df = df_full.merge(
        df_yield_region,
        on=["crop", "country", "region", "year"],
        how="left",
        suffixes=("", "_region")
    )

    # Merge country-level yields
    df = df.merge(
        df_yield_country,
        on=["crop", "country", "year"],
        how="left",
        suffixes=("", "_country")
    )

    # Fill missing region yield with country yield
    df["yield"] = df["yield"].fillna(df["yield_per_ha"])

    # Clean up
    df = df.drop(columns=["yield_per_ha"])

    return df

def min_land_to_reach_nar_given_solution(product_composition, dri, nutrient,solution):
    one_over_r = 1/dri
    hat_c = product_composition.values * one_over_r
    hat_c = hat_c.astype(float)
    return 1/np.dot(hat_c[:,nutrient],solution)