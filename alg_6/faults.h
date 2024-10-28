/*
 * faults.h
 *
 *  Created on: Oct 28, 2024
 *      Author: gustavo
 */

#ifndef FAULTS_H_
#define FAULTS_H_

float epsilon = 0.95; // Lower bound for P_ijk
// Track counters for success and failure
int step_2_failures = 0;
// Parameters for fault probability and the signing matrix
float P[v_SNOVA][l_SNOVA][l_SNOVA]; // Probability matrix

int successes_r1 = 0;
int successes_r2 = 0;
#define Q 16

// Gamma thresholds
int Gamma_epsilon_1_sup, Gamma_epsilon_1_inf;
int Gamma_epsilon_2_sup, Gamma_epsilon_2_inf;

void print_P_matrix(float m_p[v_SNOVA][l_SNOVA][l_SNOVA]) {
	for (int i = 0; i < v_SNOVA; i++) {
		printf("i = %d:\n", i);
		for (int j = 0; j < l_SNOVA; j++) {
			printf("  j = %d: ", j);
			for (int k = 0; k < l_SNOVA; k++) {
				printf("%.2f ", m_p[i][j][k]);
			}
			printf("\n"); // Newline after each row in the k dimension
		}
		printf("\n"); // Extra newline after each slice in the j dimension
	}
}

#define Q 16
#define NUM_RUNS 100
#define L l_SNOVA
#define V_SIZE (v_SNOVA * l_SNOVA * l_SNOVA)

// Function to generate a random field element in F_q, excluding fixed_value
uint8_t random_F_q_excluding(uint8_t fixed_value) {
	uint8_t random_value;
	do {
		random_value = rand() % 4;
	} while (random_value == fixed_value);
	return random_value;
}

// Function to generate a probability betwen epsilon and 1.
float random_between_ep_and_1() {
	float random_value;
	do {
		random_value = (float) rand() / RAND_MAX;
	} while (random_value < epsilon);
	return random_value;
}

// Fault injection with Bernoulli distribution
void inject_faults(uint8_t *V, uint8_t *vinegar_byte) {
	uint8_t fixed_value = random_F_q_excluding(0); // Example fixed fault value in F_q

	for (int i = 0; i < v_SNOVA; i++) {
		for (int j = 0; j < l_SNOVA; j++) {
			for (int k = 0; k < l_SNOVA; k++) {
				float res = ((float) rand() / RAND_MAX);
				if (res <= P[i][j][k]) {
					V[(i * l_SNOVA * l_SNOVA) + j * l_SNOVA + k] = fixed_value;
				} else {
					V[(i * l_SNOVA * l_SNOVA) + j * l_SNOVA + k] =
							random_F_q_excluding(fixed_value);
				}
			}
		}
	}
}

/*void compute_stats(uint8_t *V, uint8_t *c_beta) {
 // Loop through each column (beta) of V
 for (int beta = 0; beta < l_SNOVA; beta++) {
 // Loop through each possible value in F_q (Q = 16)
 for (uint8_t x = 0; x < Q; x++) {
 c_beta[beta * Q + x] = 0; // Initialize count for each value x in column beta
 // Loop through each row of V
 for (int i = 0; i < v_SNOVA; i++) {
 // Check if element in V[i][beta] (equivalently V[i * l_SNOVA + beta]) matches x
 if (V[i * l_SNOVA + beta] == x) {
 c_beta[beta * Q + x]++; // Increment count for value x in column beta
 }
 }
 }
 }
 }*/

void compute_stats(uint8_t *V, uint8_t *c_beta) {
	for (int beta = 0; beta < l_SNOVA; beta++) {
		for (uint8_t x = 0; x < Q; x++) {
			c_beta[beta * Q + x] = 0;
			for (int i = 0; i < v_SNOVA; i++) {
				for (int j = 0; j < v_SNOVA; j++) {
					if (V[i * l_SNOVA * l_SNOVA + j * l_SNOVA + beta] == x) {
						c_beta[beta * Q + x]++;
					}
				}
			}
		}
	}
}

// Initialize P matrix with epsilon value
void initialize_P_matrix(int inf, int sup) {
	// By default, each entry can be chosen uniformly randomly.
	for (int i = 0; i < v_SNOVA; ++i) {
		for (int j = 0; j < l_SNOVA; ++j) {
			for (int k = 0; k < l_SNOVA; ++k) {
				P[i][j][k] = (float) 1 / Q;
			}
		}
	}
	// Only columns indexed by inf to sup-1 of V will be affected by a fault.
	for (int i = 0; i < v_SNOVA; ++i) {
		for (int j = 0; j < l_SNOVA; ++j) {
			for (int k = inf; k < sup; ++k) {
				P[i][j][k] = random_between_ep_and_1(epsilon); // Vijk holds a fixed value with probability P[i][j][k]
			}
		}
	}
}

// Calculate Gamma thresholds based on binomial distribution
void calculate_Gamma_thresholds(int lv, float epsilon) {
	float mu = lv * epsilon;
	float sigma = sqrt(lv * epsilon * (1 - epsilon));

	Gamma_epsilon_1_sup = (int) floor(mu + sigma);
	Gamma_epsilon_1_inf = (int) floor(mu - sigma);

	Gamma_epsilon_2_sup = (int) floor(mu + 2 * sigma);
	Gamma_epsilon_2_inf = (int) floor(mu - 2 * sigma);
}

#endif /* FAULTS_H_ */
