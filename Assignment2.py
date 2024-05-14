import random
import numpy as np
# parameters
population_size = 500
generations = 2000
crossover_rate = 0.7
mutation_rate = 0.1
cycle_time = 61

# task time s
task_times = {
    1: 29, 2: 3, 3: 3, 4: 22, 5: 6, 6: 14, 7: 2, 8: 8, 9: 5, 10: 10,
    11: 23, 12: 30, 13: 23, 14: 7, 15: 19, 16: 19, 17: 2, 18: 19, 19: 16,
    20: 29, 21: 5, 22: 10, 23: 16, 24: 23, 25: 5, 26: 5, 27: 2, 28: 40,
    29: 2, 30: 5, 31: 1, 32: 1, 33: 2, 34: 2, 35: 5
}


precedence_constraints = {
    2: [1], 3: [1], 4: [1], 5: [2], 6: [5], 7: [6, 14], 8: [4],
    9: [8], 10: [1], 11: [4], 12: [10], 13: [9], 14: [1], 15: [13],
    16: [15], 17: [1], 18: [7], 19: [14, 18], 20: [19], 21: [19],
    22: [21], 23: [22], 24: [23], 25: [21], 26: [25], 27: [26],
    28: [13], 29: [28], 30: [29], 31: [30], 32: [31], 33: [24],
    34: [33], 35: [34]
}


def initialize_population():
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
    workload_deviation = sum([(cycle_time - sum(task_times[t] for t in station)) for station in stations]) / num_stations
    
    return num_stations, workload_deviation

def selection(population, fitnesses):
    total_fitness = sum(1/f for f in fitnesses)
    selection_probs = [1/f/total_fitness for f in fitnesses]
    selected = random.choices(population, weights=selection_probs, k=population_size)
    return selected

def crossover(parent1, parent2):
    if random.random() < crossover_rate:
        cut = random.randint(1, len(parent1)-1)
        child1 = parent1[:cut] + [task for task in parent2 if task not in parent1[:cut]]
        child2 = parent2[:cut] + [task for task in parent1 if task not in parent2[:cut]]
        return child1, child2
    else:
        return parent1, parent2

def mutate(chromosome):
    if random.random() < mutation_rate:
        i, j = random.sample(range(len(chromosome)), 2)
        chromosome[i], chromosome[j] = chromosome[j], chromosome[i]
    return chromosome

def genetic_algorithm():
    population = initialize_population()
    best_solution = None
    best_fitness = (float('inf'), float('inf'))  # (num_stations, workload_deviation)
    
    for generation in range(generations):
        fitnesses = [fitness(chromosome) for chromosome in population]
        num_stations = [f[0] for f in fitnesses]
        deviations = [f[1] for f in fitnesses]
        
        best_idx = np.argmin(num_stations)
        current_best = (num_stations[best_idx], deviations[best_idx])
        if current_best < best_fitness:
            best_solution = population[best_idx]
            best_fitness = current_best
        
        selected = selection(population, deviations)
        next_population = []
        
        for i in range(0, population_size, 2):
            parent1, parent2 = selected[i], selected[i+1]
            child1, child2 = crossover(parent1, parent2)
            next_population.extend([mutate(child1), mutate(child2)])
        
        population = next_population
    
    return best_solution, best_fitness

best_solution, best_fitness = genetic_algorithm()

#  results
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
workload_deviation = sum([(cycle_time - w) for w in workload_per_station]) / len(stations)

print(f"Minimum number of stations: {len(stations)}")
for i, station in enumerate(stations):
    print(f"Station {i+1}: Tasks {station}, Workload {workload_per_station[i]}")
print(f"Workload deviation: {workload_deviation}")