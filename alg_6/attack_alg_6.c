#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#include "gf16_matrix_inline.h"
#include "snova.h"
#include "ct_functions.h"
#include "aes/snova_aes.h"

#include "snova_kernel.h"
#include "util/util.h"
#define Q 16
#define NUM_RUNS 100
#define L 4 // Example for l; replace as appropriate for SNOVA
#define V_SIZE (v_SNOVA * l_SNOVA * l_SNOVA)

// Function to generate a random field element in F_q, excluding fixed_value
uint8_t random_F_q_excluding(uint8_t fixed_value) {
	uint8_t random_value;
	do {
		random_value = rand() % Q;
	} while (random_value == fixed_value);
	return random_value;
}

// Function to generate a random field element in F_q, excluding fixed_value
float random_between_ep_and_1(float epsilon) {
	float random_value;
	do {
		random_value = (float) rand() / RAND_MAX;
	} while (random_value <epsilon);
	return random_value;
}

// Parameters for fault probability and the signing matrix
float epsilon = 0.95; // Lower bound for P_ijk
float P[v_SNOVA][l_SNOVA][l_SNOVA]; // Probability matrix

// Track counters for success and failure
int step_2_failures = 0;
int successes_r1 = 0;
int successes_r2 = 0;

// Gamma thresholds
int Gamma_epsilon_1_sup, Gamma_epsilon_1_inf;
int Gamma_epsilon_2_sup, Gamma_epsilon_2_inf;

// Initialize P matrix with epsilon value
void initialize_P_matrix(float epsilon, int inf, int sup) {
	for (int i = 0; i < v_SNOVA; ++i) {
		for (int j = 0; j < l_SNOVA; ++j) {
			for (int k = 0; k < l_SNOVA; ++k) {
				P[i][j][k] = (float)1/Q;
			}
		}
	}
    
    for (int i = 0; i < v_SNOVA; ++i) {
		for (int j = 0; j < l_SNOVA; ++j) {
			for (int k = inf; k < sup; ++k) {
				P[i][j][k] = random_between_ep_and_1(epsilon);
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

// Fault injection with Bernoulli distribution
void inject_faults(uint8_t *V) {
	uint8_t fixed_value = 0x0F; // Example fixed fault value in F_q

	for (int i = 0; i < v_SNOVA; i++) {
		for (int j = 0; j < l_SNOVA; j++) {
			for (int k = 0; k < l_SNOVA; k++) {
				if ((float) rand() / RAND_MAX <= P[i][j][k]) {
					V[(i * l_SNOVA * l_SNOVA) + j * l_SNOVA + k] = fixed_value;
				} else {
					V[(i * l_SNOVA * l_SNOVA) + j * l_SNOVA + k] =
							random_F_q_excluding(fixed_value);
				}
			}
		}
	}
}

// Compute occurrences of each element in F_q for each column beta in V
void compute_stats(uint8_t *V, uint8_t *c_beta) {
	for (int beta = 0; beta < l_SNOVA; beta++) {
		for (uint8_t x = 0; x < Q; x++) {
			c_beta[beta * Q + x] = 0;
			for (int i = 0; i < v_SNOVA; i++) {
				if (V[i * l_SNOVA + beta] == x) {
					c_beta[beta * Q + x]++;
				}
			}
		}
	}
}

bool generate_signature(uint8_t *pt_signature, const uint8_t *digest,
		uint64_t bytes_digest, uint8_t *array_salt, Aalpha_t Aalpha,
		Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, T12_t T12,
		F11_t F11, F12_t F12, F21_t F21,
		const uint8_t pt_public_key_seed[seed_length_public],
		const uint8_t pt_private_key_seed[seed_length_private],
		int8_t *vinegar_in_byte) {

	gf16_t Gauss[m_SNOVA * lsq_SNOVA][m_SNOVA * lsq_SNOVA + 1];
	gf16_t Temp[lsq_SNOVA][lsq_SNOVA];
	gf16_t t_GF16, solution[m_SNOVA * lsq_SNOVA];

	gf16m_t Left_X_tmp, Right_X_tmp;
	gf16_t *Left_X, *Right_X;

	gf16m_t Left[lsq_SNOVA][v_SNOVA], Right[lsq_SNOVA][v_SNOVA];
	gf16m_t X_in_GF16Matrix[n_SNOVA] = { 0 };
	//gf16m_t X_in_GF16Matrix_faulty[n_SNOVA] = {0};
	gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
	gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
	gf16m_t signature_in_GF16Matrix[n_SNOVA] = { 0 };

	// temp
	gf16m_t gf16m_temp0;
	gf16m_t gf16m_temp1;
	gf16m_t gf16m_secret_temp0;

	int attemps = 0;

	Left_X = Left_X_tmp;
	Right_X = Right_X_tmp;
	int flag_redo = 1;
	
	int counter;

	// printf("hash value in GF16: \n");
	// print_gf16(key->hash_in_GF16, m_SNOVA * lsq_SNOVA);
	// printf("===================\n");
	//--------------------------------------
	do {
		
		flag_redo = 0;
		// put hash value in the last column of Gauss matrix
		for (int index = 0; index < (m_SNOVA * lsq_SNOVA); index++) {
			Gauss[index][m_SNOVA * lsq_SNOVA] = hash_in_GF16[index];
		}

		// Proceed with assigning vinegar values to the matrix elements
		counter = 0;
		for (int index = 0; index < v_SNOVA; index++) {
			for (int i = 0; i < rank; ++i) {
				for (int j = 0; j < rank; ++j) {
					set_gf16m(X_in_GF16Matrix[index], i, j,
							((counter & 1) ?
									(vinegar_in_byte[counter >> 1] >> 4) :
									(vinegar_in_byte[counter >> 1] & 0xF)));
					counter++;
				}
			}
		}

		// evaluate the vinegar part of central map
		evaluation_part(Aalpha, Balpha, Qalpha1, Qalpha2, F11, Left, Right,
				X_in_GF16Matrix, Fvv_in_GF16Matrix, gf16m_temp0, gf16m_temp1);

		// add to the last column of Gauss matrix
		add_last_column(Gauss, Fvv_in_GF16Matrix);

		// compute the coefficients of Xo and put into Gauss matrix and compute
		// the coefficients of Xo^t and add into Gauss matrix
		for (int i = 0; i < m_SNOVA; ++i) {
			for (int index = 0; index < o_SNOVA; ++index) {
				clean_temp(Temp);

				for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
					for (int j = 0; j < v_SNOVA; ++j) {
						gf16m_mul(Left[alpha][j], F12[i][j][index],
								gf16m_temp0);
						gf16m_mul(gf16m_temp0, Qalpha2[alpha], Left_X_tmp);
						Left_X = Left_X_tmp;
						Right_X = Balpha[alpha];
						/*
						 Left_X = Left[alpha][j] * F12[i][j][k] * Qalpha2[alpha];
						 Right_X = Balpha[alpha];
						 */
						for (int ti = 0; ti < lsq_SNOVA; ++ti) {
							for (int tj = 0; tj < lsq_SNOVA; ++tj) {
								gf16_t temp3 = 0;
								temp3 =
										gf16_get_mul(
												get_gf16m(Left_X, ti / rank, tj / rank),
												get_gf16m(Right_X, tj % rank, ti % rank));
								Temp[ti][tj] = gf16_get_add(Temp[ti][tj],
										temp3);
							}
						}
					}
				}

				for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
					for (int j = 0; j < v_SNOVA; ++j) {
						Left_X = Aalpha[alpha];
						gf16m_mul(Qalpha1[alpha], F21[i][index][j],
								gf16m_temp0);
						gf16m_mul(gf16m_temp0, Right[alpha][j], Right_X_tmp);
						Right_X = Right_X_tmp;
						/*
						 Left_X = Aalpha[alpha];
						 Right_X = Qalpha1[alpha] * F21[i][k][j] *
						 Right[alpha][j];
						 */
						for (int ti = 0; ti < lsq_SNOVA; ++ti) {
							for (int tj = 0; tj < lsq_SNOVA; ++tj) {
								gf16_t temp2 = 0;
								temp2 =
										gf16_get_mul(
												get_gf16m(Left_X, ti / rank, tj % rank),
												get_gf16m(Right_X, tj / rank, ti % rank));
								Temp[ti][tj] = gf16_get_add(Temp[ti][tj],
										temp2);
							}
						}
					}
				}
				for (int ti = 0; ti < lsq_SNOVA; ++ti) {
					for (int tj = 0; tj < lsq_SNOVA; ++tj) {
						Gauss[i * lsq_SNOVA + ti][index * lsq_SNOVA + tj] =
								Temp[ti][tj];
					}
				}
			}
		}

		// Gauss elimination
		for (int i = 0; i < m_SNOVA * lsq_SNOVA; ++i) {
			if (Gauss[i][i] == 0) {
				for (int j = i + 1; j < m_SNOVA * lsq_SNOVA; ++j) {
					if (Gauss[j][i] != 0) {
						for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
							t_GF16 = Gauss[i][k];
							Gauss[i][k] = Gauss[j][k];
							Gauss[j][k] = t_GF16;
						}
						break;
					}
				}
			}
			if (Gauss[i][i] == 0) {
				flag_redo = 1;
				break;
			}

			t_GF16 = inv(Gauss[i][i]);
			for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
				Gauss[i][k] = gf16_get_mul(Gauss[i][k], t_GF16);
			}

			for (int j = i + 1; j < m_SNOVA * lsq_SNOVA; ++j) {
				if (Gauss[j][i] != 0) {
					t_GF16 = Gauss[j][i];
					for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
						Gauss[j][k] = gf16_get_add(Gauss[j][k],
								gf16_get_mul(Gauss[i][k], t_GF16));
					}
				}
			}
		}

		if (!flag_redo) {
			for (int i = m_SNOVA * lsq_SNOVA - 1; i >= 0; --i) {
				t_GF16 = 0;
				for (int k = i + 1; k < m_SNOVA * lsq_SNOVA; ++k) {
					t_GF16 = gf16_get_add(t_GF16,
							gf16_get_mul(Gauss[i][k], solution[k]));
				}
				solution[i] = gf16_get_add(Gauss[i][m_SNOVA * lsq_SNOVA],
						t_GF16);
			}
		}

		attemps++;

	} while (attemps <= 0);
	if (flag_redo) {
		return 0;
	}

	

	for (int index = 0; index < o_SNOVA; ++index) {
		for (int i = 0; i < rank; ++i) {
			for (int j = 0; j < rank; ++j) {
				set_gf16m(X_in_GF16Matrix[index + v_SNOVA], i, j,
						solution[index * lsq_SNOVA + i * rank + j]);
			}
		}
	}

	//gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};

	for (int i = 0; i < n_SNOVA; i++) {
		// Set the first 'rank' elements (first row) to 1
		for (int j = 0; j < rank; j++) {
			X_in_GF16Matrix[i][j] = 1;
		}
	}

	printf("X_in_GF16Matrix:\n");
	print_all_matrices(X_in_GF16Matrix, n_SNOVA);

	for (int index = 0; index < v_SNOVA; ++index) {
		gf16m_clone(signature_in_GF16Matrix[index], X_in_GF16Matrix[index]);
		for (int i = 0; i < o_SNOVA; ++i) {
			gf16m_mul(T12[index][i], X_in_GF16Matrix[v_SNOVA + i],
					gf16m_secret_temp0);
			gf16m_add(signature_in_GF16Matrix[index], gf16m_secret_temp0,
					signature_in_GF16Matrix[index]);
			/*
			 signature_in_GF16Matrix[index] = signature_in_GF16Matrix[index] +
			 T12[index][i] * X_in_GF16Matrix[v_SNOVA + i];
			 */
		}
	}
	for (int index = 0; index < o_SNOVA; ++index) {
		gf16m_clone(signature_in_GF16Matrix[v_SNOVA + index],
				X_in_GF16Matrix[v_SNOVA + index]);
	}
	printf("signature_in_GF16Matrix:\n");
	print_all_matrices(signature_in_GF16Matrix, n_SNOVA);

	// output signature
	for (int index = 0; index < n_SNOVA * lsq_SNOVA; ++index) {
		((gf16_t*) signature_in_GF16Matrix)[index] = get_gf16m(
				signature_in_GF16Matrix[index / lsq_SNOVA],
				(index % lsq_SNOVA) / l_SNOVA, (index % lsq_SNOVA) % l_SNOVA);
	}
	convert_GF16s_to_bytes(pt_signature, (gf16_t*) signature_in_GF16Matrix,
	n_SNOVA * lsq_SNOVA);
	for (int i = 0; i < bytes_salt; ++i) {
		pt_signature[bytes_signature + i] = array_salt[i];
	}

	SNOVA_CLEAR(gf16m_secret_temp0);
	return 1;

}

// Modified signing algorithm with fault injection and statistics collection
bool sign_with_fault_injection(uint8_t *pt_signature, const uint8_t *digest,
		uint64_t bytes_digest, uint8_t *array_salt, Aalpha_t Aalpha,
		Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, T12_t T12,
		F11_t F11, F12_t F12, F21_t F21,
		const uint8_t pt_public_key_seed[seed_length_public],
		const uint8_t pt_private_key_seed[seed_length_private]) {
	gf16_t Gauss[m_SNOVA * lsq_SNOVA][m_SNOVA * lsq_SNOVA + 1];
	gf16_t Temp[lsq_SNOVA][lsq_SNOVA];
	gf16_t t_GF16, solution[m_SNOVA * lsq_SNOVA];

	gf16m_t Left_X_tmp, Right_X_tmp;
	gf16_t *Left_X, *Right_X;

	gf16m_t Left[lsq_SNOVA][v_SNOVA], Right[lsq_SNOVA][v_SNOVA];
	gf16m_t X_in_GF16Matrix[n_SNOVA] = { 0 };
	//gf16m_t X_in_GF16Matrix_faulty[n_SNOVA] = {0};
	gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
	gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
	gf16m_t signature_in_GF16Matrix[n_SNOVA] = { 0 };

	uint8_t signed_hash[bytes_hash];
	uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1];

	// temp
	gf16m_t gf16m_temp0;
	gf16m_t gf16m_temp1;
	gf16m_t gf16m_secret_temp0;

	Left_X = Left_X_tmp;
	Right_X = Right_X_tmp;
	int flag_redo = 1;
	
	int counter;

	bool signature_success = false;
	createSignedHash(digest, bytes_digest, pt_public_key_seed, array_salt,
			signed_hash);
	convert_bytes_to_GF16s(signed_hash, hash_in_GF16, GF16s_hash);

	uint8_t V[(v_SNOVA * lsq_SNOVA + 1) >> 1];

	do {
		
		signature_success = false;

		inject_faults(V); // Generate values for V with fault injection

		// Proceed with signature generation (assumes a function `generate_signature`)
		signature_success = generate_signature(pt_signature, digest,
				bytes_digest, array_salt, Aalpha, Balpha, Qalpha1, Qalpha2, T12,
				F11, F12, F21, pt_public_key_seed, pt_private_key_seed, V);
		if (!signature_success) {
			step_2_failures++;
			return false;
		}
	} while (!signature_success);

	return signature_success;
}

void read_from_file(const char *filename, uint8_t *pk, size_t nr_bytes_pk,
		uint8_t *sk, size_t nr_bytes_sk, uint8_t *seed_, size_t nr_seed_bytes) {
	FILE *file = fopen(filename, "rb");
	if (file == NULL) {
		perror("Error opening file for reading");
		exit(EXIT_FAILURE);
	}

	// Read pk
	fread(pk, sizeof(uint8_t), nr_bytes_pk, file);
	// Read sk
	fread(sk, sizeof(uint8_t), nr_bytes_sk, file);
	// Read array_salt
	fread(seed_, sizeof(uint8_t), nr_seed_bytes, file);

	fclose(file);
}

// Run experiments to collect statistics
void run_experiments() {
	int bytes_digest = 64;
	srand(time(0)); // Initialize random number generator
	uint8_t pt_signature[bytes_signature + bytes_salt];
	uint8_t V[V_SIZE]; // Vinegar variable storage
	uint8_t digest[bytes_digest];
	uint8_t c_beta[l_SNOVA * Q]; // Counter for occurrences of elements in V columns


	uint8_t array_digest[64];
	uint8_t array_signature1[bytes_signature + bytes_salt];
	uint8_t array_signature2[bytes_signature + bytes_salt];

	uint8_t seed[seed_length];
	uint8_t *pt_private_key_seed;
	uint8_t *pt_public_key_seed;
	uint8_t pk[bytes_pk], sk[bytes_sk];
	uint8_t array_salt[bytes_salt];

	read_from_file("keys.txt", pk, bytes_pk, sk, bytes_sk, seed, seed_length);

	sk_gf16 sk_upk;
	sk_unpack(&sk_upk, sk);

	// Clear Secret!
	SNOVA_CLEAR_BYTE(&sk_upk, sizeof(sk_upk));
	for (int run = 0; run < NUM_RUNS; run++) {
		bool result = sign_with_fault_injection(pt_signature, digest,
				bytes_digest, array_salt, sk_upk.Aalpha, sk_upk.Balpha,
				sk_upk.Qalpha1, sk_upk.Qalpha2, sk_upk.T12, sk_upk.F11,
				sk_upk.F12, sk_upk.F21, sk_upk.pt_public_key_seed,
				sk_upk.pt_private_key_seed);

		if (!result)
			continue; // If Step 2 fails, skip further steps

		compute_stats(V, c_beta);

		for (int beta = 0; beta < l_SNOVA; beta++) {
			for (uint8_t x = 0; x < Q; x++) {
				if (Gamma_epsilon_1_sup >= c_beta[beta * Q + x]
						&& c_beta[beta * Q + x] >= Gamma_epsilon_1_inf) {
					successes_r1++;
					successes_r2++;
					goto next_run;
				}
				if (Gamma_epsilon_2_sup >= c_beta[beta * Q + x]
						&& c_beta[beta * Q + x] >= Gamma_epsilon_2_inf) {
					successes_r2++;
					goto next_run;
				}
			}
		}
		next_run: continue;
	}
}

// Compute and print final results
void compute_results() {
	float FR = (float) step_2_failures / NUM_RUNS;
	float AlgoSR_1 = (float) successes_r1 / (NUM_RUNS - step_2_failures);
	float AlgoSR_2 = (float) successes_r2 / (NUM_RUNS - step_2_failures);
	float AttackSR_1 = (float) successes_r1 / NUM_RUNS;
	float AttackSR_2 = (float) successes_r2 / NUM_RUNS;

	printf("Failure Rate (FR): %.2f\n", FR);
	printf("Algorithm Success Rate 1 (AlgoSR_1): %.2f\n", AlgoSR_1);
	printf("Algorithm Success Rate 2 (AlgoSR_2): %.2f\n", AlgoSR_2);
	printf("Attack Success Rate 1 (AttackSR_1): %.2f\n", AttackSR_1);
	printf("Attack Success Rate 2 (AttackSR_2): %.2f\n", AttackSR_2);
}

int main() {
	snova_init();
	int lv = v_SNOVA * l_SNOVA;
	initialize_P_matrix(epsilon, 0, l);
	calculate_Gamma_thresholds(lv, epsilon); // Calculate thresholds for success conditions
	run_experiments();
	compute_results();
	return 0;
}
