import numpy as np
import matplotlib.pyplot as plt

def hamming_weight(x):
    return bin(x).count('1')

def create_templates(a):
    templates = {}
    for b in range(16):
        templates[b] = hamming_weight(a ^ b)
    return templates

def simulate_power_trace(a, b, noise_level=0.1):
    hw = hamming_weight(a ^ b)
    noise = np.random.normal(0, noise_level)
    return hw + noise

def side_channel_attack(a, observed_traces):
    templates = create_templates(a)
    results = []
    for trace in observed_traces:
        differences = {b: abs(trace - templates[b]) for b in templates}
        guessed_b = min(differences, key=differences.get)
        results.append(guessed_b)
    return results

# Parameters
num_trials = 1000
a = 5  # Known value of a, can be any integer from 0 to 15
num_traces_per_trial = 10  # Number of traces to use in each trial

correct_guesses_total = 0
total_guesses = 0

for _ in range(num_trials):
    true_b_values = np.random.randint(0, 16, num_traces_per_trial)
    observed_traces = [simulate_power_trace(a, b) for b in true_b_values]
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
import numpy as np
import matplotlib.pyplot as plt

def hamming_weight(x):
    return bin(x).count('1')

def create_templates(a):
    templates = {}
    for b in range(16):
        templates[b] = hamming_weight(a ^ b)
    return templates

def simulate_power_trace(a, b, noise_level=0.1):
    hw = hamming_weight(a ^ b)
    noise = np.random.normal(0, noise_level)
    return hw + noise

def side_channel_attack(a, observed_traces):
    templates = create_templates(a)
    results = []
    for trace in observed_traces:
        differences = {b: abs(trace - templates[b]) for b in templates}
        guessed_b = min(differences, key=differences.get)
        results.append(guessed_b)
    return results

# Parameters
num_trials = 1000
a = 5  # Known value of a, can be any integer from 0 to 15
num_traces_per_trial = 10  # Number of traces to use in each trial

correct_guesses_total = 0
total_guesses = 0

for _ in range(num_trials):
    true_b_values = np.random.randint(0, 16, num_traces_per_trial)
    observed_traces = [simulate_power_trace(a, b) for b in true_b_values]
    guessed_b_values = side_channel_attack(a, observed_traces)
    
    correct_guesses = sum(1 for true_b, guessed_b in zip(true_b_values, guessed_b_values) if true_b == guessed_b)
    correct_guesses_total += correct_guesses
    total_guesses += num_traces_per_trial

accuracy = correct_guesses_total / total_guesses
print(f"Accuracy over {num_trials} trials: {accuracy:.2f}")

# Optionally, you can visualize the results
plt.hist([simulate_power_trace(a, b) for b in range(16)], bins=16, alpha=0.75, label='Simulated Traces')
plt.xlabel('Power Trace Value')
plt.ylabel('Frequency')
plt.title('Histogram of Simulated Power Traces')
plt.legend()
plt.show()