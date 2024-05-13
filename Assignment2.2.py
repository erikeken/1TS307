import random
import numpy as np

assembly_times = [0, 29, 3, 5, 22, 8, 6, 14, 2, 9, 5, 23, 30, 23, 2, 19, 16, 2, 19, 29, 2, 19, 10, 16, 23, 5, 5, 2, 40, 29, 30, 31, 1, 2, 2]

precedence_constraints = {
    2: [1], 3: [1], 4: [2, 3], 5: [1], 6: [5], 7: [5], 8: [4], 9: [4], 10: [1],
    11: [4], 12: [10], 13: [8], 14: [7], 15: [14, 13], 16: [15], 17: [21], 18: [15],
    19: [18], 20: [17], 21: [19], 22: [21], 23: [22], 24: [23], 25: [24], 26: [25],
    27: [26, 34], 28: [27, 13], 29: [28, 2], 30: [29, 12, 6], 31: [30], 32: [31],
    33: [28], 34: [32, 27], 35: [33]
}

max_cycle_time = 61

# parameters
population_size = 1000
crossover_rate = 0.8
mutation_rate = 0.1
generations = 1000

def initialize_population(size, num_parts):
    population = []
    for _ in range(size):
        individual = list(np.random.permutation(num_parts) + 1)
        population.append(individual)
    return population

def is_valid(individual):
    for part, predecessors in precedence_constraints.items():
        for pre in predecessors:
            if individual.index(pre) > individual.index(part):
                return False
    return True

def evaluate_fitness(individual, assembly_times, max_cycle_time):
    workload = [0]
    station = 0
    for part in individual:
        if workload[station] + assembly_times[part] > max_cycle_time:
            station += 1
            workload.append(0)
        workload[station] += assembly_times[part]
    num_stations = len(workload)
    workload_deviation = np.std(workload)
    return num_stations, workload, workload_deviation

def selection(population, fitness, num_parents):
    fitness_sorted = sorted(zip(population, fitness), key=lambda x: x[1][0])
    selected_parents = [individual for individual, _ in fitness_sorted[:num_parents]]
    return selected_parents

def crossover(parents, crossover_rate):
    offspring = []
    for _ in range(len(parents) // 2):
        if random.random() < crossover_rate:
            parent1, parent2 = random.sample(parents, 2)
            point = random.randint(1, len(parent1) - 1)
            child1 = parent1[:point] + [p for p in parent2 if p not in parent1[:point]]
            child2 = parent2[:point] + [p for p in parent1 if p not in parent2[:point]]
            offspring.extend([child1, child2])
        else:
            offspring.extend(random.sample(parents, 2))
    return [child for child in offspring if is_valid(child)]

def mutation(offspring, mutation_rate):
    for child in offspring:
        if random.random() < mutation_rate:
            idx1, idx2 = random.sample(range(len(child)), 2)
            child[idx1], child[idx2] = child[idx2], child[idx1]
    return [child for child in offspring if is_valid(child)]

def replace_population(population, offspring, fitness, fitness_offspring):
    combined = list(zip(population + offspring, fitness + fitness_offspring))
    combined_sorted = sorted(combined, key=lambda x: x[1][0])
    new_population = [individual for individual, _ in combined_sorted[:len(population)]]
    return new_population

def find_best_solution(population, fitness):
    best_index = np.argmin([fit[0] for fit in fitness])
    return population[best_index], fitness[best_index]

# population
population = initialize_population(population_size, len(assembly_times) - 1)

for gen in range(generations):
    fitness = [evaluate_fitness(ind, assembly_times, max_cycle_time) for ind in population]
    
    # selection
    num_parents = population_size // 2
    parents = selection(population, fitness, num_parents)
    
    # crossover
    offspring = crossover(parents, crossover_rate)
    
    # mutation
    offspring = mutation(offspring, mutation_rate)
    
    # evaluate fitness of offspring
    fitness_offspring = [evaluate_fitness(ind, assembly_times, max_cycle_time) for ind in offspring]
    
    # replacement
    population = replace_population(population, offspring, fitness, fitness_offspring)

# best solution
best_solution, best_fitness = find_best_solution(population, [evaluate_fitness(ind, assembly_times, max_cycle_time) for ind in population])

# Display results
num_stations, workload, workload_deviation = best_fitness
print(f"Minimum number of stations: {num_stations}")
print(f"Parts assembled at each station: {best_solution}")
print(f"Total workload at each station: {workload}")
print(f"Workload deviation: {workload_deviation}")
