import random
import numpy as np

# Task times in seconds
task_times = {
    1: 29, 2: 3, 3: 3, 4: 22, 5: 6, 6: 14, 7: 2, 8: 8, 9: 5, 10: 10,
    11: 23, 12: 30, 13: 23, 14: 7, 15: 19, 16: 19, 17: 2, 18: 19, 19: 16,
    20: 29, 21: 5, 22: 10, 23: 16, 24: 23, 25: 5, 26: 5, 27: 2, 28: 40,
    29: 2, 30: 5, 31: 1, 32: 1, 33: 2, 34: 2, 35: 5
}

# Precedence constraints (not used in this simple example, but can be integrated)
precedence_constraints = {
    2: [1], 3: [1], 4: [1], 5: [2], 6: [5], 7: [6, 14], 8: [4],
    9: [8], 10: [1], 11: [4], 12: [10], 13: [9], 14: [1], 15: [13],
    16: [15], 17: [1], 18: [7], 19: [14, 18], 20: [19], 21: [19],
    22: [21], 23: [22], 24: [23], 25: [21], 26: [25], 27: [26],
    28: [13], 29: [28], 30: [29], 31: [30], 32: [31], 33: [24],
    34: [33], 35: [34]
}

# Parameter exploration settings
params = [
    {'population_size': 50, 'generations': 1000, 'crossover_rate': 0.8, 'mutation_rate': 0.2},
    {'population_size': 100, 'generations': 1500, 'crossover_rate': 0.85, 'mutation_rate': 0.15},
    {'population_size': 75, 'generations': 2000, 'crossover_rate': 0.9, 'mutation_rate': 0.1},
    {'population_size': 200, 'generations': 5000, 'crossover_rate': 0.7, 'mutation_rate': 0.1},
    {'population_size': 500, 'generations': 2000, 'crossover_rate': 0.7, 'mutation_rate': 0.1}
]

cycle_time = 61

def initialize_population(population_size):
    population = []
    for _ in range(population_size):
        chromosome = list(task_times.keys())
        random.shuffle(chromosome)
        population.append(chromosome)
    return population

def fitness(chromosome):
    stations = []
    current_station = []
    current_time = 0

    for task in chromosome:
        task_time = task_times[task]
        if current_time + task_time > cycle_time:
            stations.append(current_station)
            current_station = [task]
            current_time = task_time
        else:
            current_station.append(task)
            current_time += task_time
    stations.append(current_station)

    num_stations = len(stations)
    workload_deviation = sum([abs(cycle_time - sum(task_times[t] for t in station)) for station in stations]) / num_stations
    
    return num_stations, workload_deviation

def selection(population, fitnesses):
    total_fitness = sum(1/f[1] for f in fitnesses)
    selection_probs = [1/f[1]/total_fitness for f in fitnesses]
    selected = random.choices(population, weights=selection_probs, k=len(population))
    if len(selected) % 2 != 0:
        selected.append(random.choice(selected))  # Ensure even number of parents
    return selected

def crossover(parent1, parent2, crossover_rate):
    if random.random() < crossover_rate:
        cut = random.randint(1, len(parent1)-1)
        child1 = parent1[:cut] + [task for task in parent2 if task not in parent1[:cut]]
        child2 = parent2[:cut] + [task for task in parent1 if task not in parent2[:cut]]
        return child1, child2
    else:
        return parent1, parent2

def mutate(chromosome, mutation_rate):
    if random.random() < mutation_rate:
        i, j = random.sample(range(len(chromosome)), 2)
        chromosome[i], chromosome[j] = chromosome[j], chromosome[i]
    return chromosome

def genetic_algorithm(population_size, generations, crossover_rate, mutation_rate):
    population = initialize_population(population_size)
    best_solution = None
    best_fitness = (float('inf'), float('inf'))  # (num_stations, workload_deviation)
    
    for generation in range(generations):
        fitnesses = [fitness(chromosome) for chromosome in population]
        num_stations = [f[0] for f in fitnesses]
        deviations = [f[1] for f in fitnesses]
        
        best_idx = np.argmin(deviations)
        current_best = (num_stations[best_idx], deviations[best_idx])
        if current_best < best_fitness:
            best_solution = population[best_idx]
            best_fitness = current_best
        
        selected = selection(population, fitnesses)
        next_population = []
        
        for i in range(0, len(selected), 2):
            parent1, parent2 = selected[i], selected[i+1]
            child1, child2 = crossover(parent1, parent2, crossover_rate)
            next_population.extend([mutate(child1, mutation_rate), mutate(child2, mutation_rate)])
        
        population = next_population
    
    return best_solution, best_fitness

best_overall_solution = None
best_overall_fitness = (float('inf'), float('inf'))  # (num_stations, workload_deviation)
best_params = None

for param in params:
    population_size = param['population_size']
    generations = param['generations']
    crossover_rate = param['crossover_rate']
    mutation_rate = param['mutation_rate']
    
    best_solution, best_fitness = genetic_algorithm(population_size, generations, crossover_rate, mutation_rate)
    
    if best_fitness < best_overall_fitness:
        best_overall_solution = best_solution
        best_overall_fitness = best_fitness
        best_params = param

    # Report results for each parameter setting
    stations = []
    current_station = []
    current_time = 0
    for task in best_solution:
        task_time = task_times[task]
        if current_time + task_time > cycle_time:
            stations.append(current_station)
            current_station = [task]
            current_time = task_time
        else:
            current_station.append(task)
            current_time += task_time
    stations.append(current_station)

    workload_per_station = [sum(task_times[t] for t in station) for station in stations]
    workload_deviation = sum([abs(cycle_time - w) for w in workload_per_station]) / len(stations)

    print(f"Parameter setting: {param}")
    print(f"Minimum number of stations: {len(stations)}")
    for i, station in enumerate(stations):
        print(f"Station {i+1}: Tasks {station}, Workload {workload_per_station[i]}")
    print(f"Workload deviation: {workload_deviation}\n")

# Report the best overall solution
print("Best solution:")
print(f"Parameters: {best_params}")
stations = []
current_station = []
current_time = 0
for task in best_overall_solution:
    task_time = task_times[task]
    if current_time + task_time > cycle_time:
        stations.append(current_station)
        current_station = [task]
        current_time = task_time
    else:
        current_station.append(task)
        current_time += task_time
stations.append(current_station)

workload_per_station = [sum(task_times[t] for t in station) for station in stations]
workload_deviation = sum([abs(cycle_time - w) for w in workload_per_station]) / len(stations)

print(f"Minimum number of stations: {len(stations)}")
for i, station in enumerate(stations):
    print(f"Station {i+1}: Tasks {station}, Workload {workload_per_station[i]}")
print(f"Workload deviation: {workload_deviation}")
