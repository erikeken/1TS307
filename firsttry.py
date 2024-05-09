from pyomo.environ import *
import sys

model = ConcreteModel()

factories = [1, 2]
shops = [1, 2, 3, 4]
distribution_centers = [1, 2]

""" Parameters """
cost_factory_to_shop = {
    (1, 1): 1, (1, 2): 0.5, (1, 3): 1, (1, 4): 0.75,
    (2, 1): 1.5, (2, 2): 1, (2, 3): 1, (2, 4): 0.25
}
cost_factory_to_dc = {
    (1, 1): 0.5, (1, 2): 0.25,
    (2, 1): 0.15, (2, 2): 0.75
}
cost_dc_to_shop = {
    (1, 1): 1.25, (1, 2): 0.75, (1, 3): 0.75, (1, 4): 0.25,
    (2, 1): 0.5, (2, 2): 0.1, (2, 3): 0.6, (2, 4): 0.4
}

factory_capacity = {1: 30000, 2: 50000}
shop_demand = {1: 2000, 2: 5000, 3: 15000, 4: 20000}

""" Decision Variables """
model.x = Var(factories, shops, within=NonNegativeReals) # Factory -> shop
model.y = Var(factories, distribution_centers, within=NonNegativeReals) # Factory -> distribution center
model.z = Var(distribution_centers, shops, within=NonNegativeReals) # Distribution center -> shop

""" Objective """
def objective(model):
    func = sum(cost_factory_to_shop[i, k] * model.x[i, k] for i in factories for k in shops) + \
           sum(cost_factory_to_dc[i, j] * model.y[i, j] for i in factories for j in distribution_centers) + \
           sum(cost_dc_to_shop[j, k] * model.z[j, k] for j in distribution_centers for k in shops)
    return func

model.objective = Objective(rule=objective, sense=minimize)

""" Constraints """
# Factory capacity constraint
def factory_capacity_(model, i):
    return sum(model.x[i, k] for k in shops) + sum(model.y[i, j] for j in distribution_centers) <= factory_capacity[i]

model.factory_capacity_con = Constraint(factories, rule=factory_capacity_)

# Shop demand constraint
def shop_demand_(model, k):
    return sum(model.x[i, k] for i in factories) + sum(model.z[j, k] for j in distribution_centers) == shop_demand[k]

model.shop_demand_con = Constraint(shops, rule=shop_demand_)

# Distribution center balance
def distribution_center_(model, j):
    return sum(model.y[i, j] for i in factories) == sum(model.z[j, k] for k in shops)

model.distribution_center_con = Constraint(distribution_centers, rule=distribution_center_)

# Solver
solver = SolverFactory('glpk')
result = solver.solve(model, tee=True)

# Print results
print("Status:", result.solver.status)
print("Total Cost:", model.objective())

def printout(file=None):
    for i in factories:
        for j in distribution_centers:
            print(f"Factory {i} to DC {j}: {model.y[i, j].value}", file=file)
    for j in distribution_centers:
        for k in shops:
            print(f"DC {j} to Shop {k}: {model.z[j, k].value}", file=file)
    for i in factories:
        for k in shops:
            print(f"Factory {i} to Shop {k}: {model.x[i, k].value}", file=file)

with open("output.txt", 'w') as f:
    printout(f)

# Display decision variable amounts
printout()