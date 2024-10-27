import numpy as np
import matplotlib.pyplot as plt

# Helper function to compute the Hamming weight
def hamming_weight(x):
    return bin(x).count('1')

# Function to create templates for a given value of a
def create_templates(a):
    templates = {}
    for b in range(16):
        templates[b] = hamming_weight(a ^ b)
    return templates

# Function to simulate the power trace with some noise
def simulate_power_trace(a, b, noise_level=0.1):
    hw = hamming_weight(a ^ b)
    noise = np.random.normal(0, noise_level)
    return hw + noise

# Function to perform the side-channel attack
def side_channel_attack(a, observed_traces):
    templates = create_templates(a)
    results = []
    for trace in observed_traces:
        differences = {b: abs(trace - templates[b]) for b in templates}
        guessed_b = min(differences, key=differences.get)
        results.append(guessed_b)
    return results

# Function to generate a random 4x4 matrix of rank 4 in GF(16)
def generate_random_rank_4_matrix():
    while True:
        matrix = np.random.randint(0, 16, (4, 4))
        if np.linalg.matrix_rank(matrix) == 4:
            return matrix

# Parameters
num_trials = 1000
num_traces_per_trial = 10  # Number of traces to use in each trial
noise_level = 0.1

correct_guesses_total = 0
total_guesses = 0

for _ in range(num_trials):
    matrix = generate_random_rank_4_matrix()
    a = matrix[0, 0]  # We'll attack the element at position (0,0)
    
    true_b_values = np.random.randint(0, 16, num_traces_per_trial)
    observed_traces = [simulate_power_trace(a, b, noise_level) for b in true_b_values]
    
    guessed_b_values = side_channel_attack(a, observed_traces)
    
    correct_guesses = sum(1 for true_b, guessed_b in zip(true_b_values, guessed_b_values) if true_b == guessed_b)
    correct_guesses_total += correct_guesses
    total_guesses += num_traces_per_trial

accuracy = correct_guesses_total / total_guesses
print(f"Accuracy over {num_trials} trials: {accuracy:.2f}")

# Optionally, you can visualize the results
# plt.hist([simulate_power_trace(a, b) for b in range(16)], bins=16, alpha=0.75, label='Simulated Traces')
# plt.xlabel('Power Trace Value')
# plt.ylabel('Frequency')
# plt.title('Histogram of Simulated Power Traces')
# plt.legend()
# plt.show()
