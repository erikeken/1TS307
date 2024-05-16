import numpy as np

# Define problem inputs
CT = 61  # Cycle time
n = 35  # Number of tasks

def convert_string_to_array(data_string):
    lines = data_string.strip().split('\n')
    data = [list(map(int, line.split('\t'))) for line in lines]
    array = np.array(data)
    return array

data_string = """0	1	0	0	1	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	1	0	1	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	1	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0"""
PR = convert_string_to_array(data_string)

t = np.array([
    29, 3, 5, 22, 6, 14, 2, 5, 22, 30, 23, 30, 23, 2, 19, 29,
    2, 2, 19, 29, 6, 10, 16, 23, 5, 5, 5, 40, 2, 5, 5, 1, 40, 2, 2
])

# GA parameters
popsize = 10
r_cr = 0.7  # Crossover rate
n_cr = int(popsize * r_cr)  # Number of crossovers
r_mu = 0.1  # Mutation rate
n_mu = int(popsize * r_mu)  # Number of mutations
r_el = 0.1  # Elitism rate
n_el = int(popsize * r_el)  # Number of elite individuals
G = 0
Max_G = 10  # Maximum generations

# Initialize random population
TP = np.random.rand(popsize, n)

# Encoding and Decoding function
def encode_and_decode(TP):
    TS = np.zeros((popsize, n), dtype=int)
    AsTaSt = np.zeros((popsize, n), dtype=int)
    
    for i in range(popsize):
        PR2 = PR.copy()
        tac = 0
        while tac < n:
            prsum = np.sum(PR2, axis=0)
            clt = np.where(prsum == 0)[0]
            if len(clt) == 0:
                break
            maxup = -1
            slt = -1
            for j in clt:
                if TP[i, j] > maxup:
                    maxup = TP[i, j]
                    slt = j
            TS[i, tac] = slt
            tac += 1
            PR2[slt, :] = 0
            PR2[:, slt] = 1
        
        # Decoding
        stc = 1
        tac = 0
        tt = 0
        while tac < n:
            if tt + t[TS[i, tac]] <= CT:
                AsTaSt[i, tac] = stc
                tt += t[TS[i, tac]]
                tac += 1
            else:
                stc += 1
                tt = 0
    return TS, AsTaSt

# Initial encoding and decoding
TS, AsTaSt = encode_and_decode(TP)
fitness = np.max(AsTaSt, axis=1)

# Main body of GA
while G < Max_G:
    G += 1
    n_TP = np.zeros_like(TP)
    
    # Crossover
    for i in range(0, n_cr, 2):
        # Select parents
        par1, par2 = np.random.choice(popsize, 2, replace=False)
        # Select crossover points
        pos1, pos2 = sorted(np.random.choice(n, 2, replace=False))
        n_TP[i] = TP[par1]
        n_TP[i + 1] = TP[par2]
        n_TP[i, pos1:pos2] = TP[par2, pos1:pos2]
        n_TP[i + 1, pos1:pos2] = TP[par1, pos1:pos2]

    # Mutation
    for i in range(n_cr, n_cr + n_mu):
        idx = np.random.randint(0, popsize)
        pos1, pos2 = np.random.choice(n, 2, replace=False)
        n_TP[i] = TP[idx]
        n_TP[i, [pos1, pos2]] = n_TP[i, [pos2, pos1]]

    # Elitism
    elite_indices = np.argsort(fitness)[:n_el]
    for i, idx in enumerate(elite_indices, start=n_cr + n_mu):
        n_TP[i] = TP[idx]

    # Replace old population
    TP = n_TP

    # Recalculate task sequences and station assignments
    TS, AsTaSt = encode_and_decode(TP)
    fitness = np.max(AsTaSt, axis=1)  # Update fitness

# Find best solution
best_index = np.argmin(fitness)
best_sequence = TS[best_index]
best_stations = AsTaSt[best_index]

# Calculate results
unique_stations = np.unique(best_stations)
num_stations = len(unique_stations)
parts_at_stations = {station: np.where(best_stations == station)[0] + 1 for station in unique_stations}
workloads = {station: t[parts_at_stations[station] - 1].sum() for station in unique_stations}
avg_workload = np.mean(list(workloads.values()))
workload_deviation = sum(abs(CT - workload) for workload in workloads.values()) / num_stations

# Prepare a formatted string for output
output_lines = ["Parts assembled at each station:"]
for station, parts in parts_at_stations.items():
    # Convert the array of parts into a comma-separated string
    parts_list = ', '.join(map(str, parts))
    output_lines.append(f"Station {station}: [ {parts_list} ]")
output = '\n'.join(output_lines)


print("Minimum number of stations:", num_stations)
print("Parts assembled at each station:", output)
print("Total workload at each station:", workloads)
print("Workload deviation:", workload_deviation)