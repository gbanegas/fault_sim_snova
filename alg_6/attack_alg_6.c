#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#include "gf16_matrix_inline.h"

#include "ct_functions.h"
#include "aes/snova_aes.h"

#include "snova_kernel_fault.h"
#include "util/util.h"

int successes_r1 = 0;
int successes_r2 = 0;

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
	uint8_t digest[64];
	uint8_t c_beta[l_SNOVA * Q]; // Counter for occurrences of elements in V columns

	uint8_t seed[seed_length];
	uint8_t pk[bytes_pk], sk[bytes_sk];
	uint8_t array_salt[bytes_salt];

	read_from_file(
			"/home/gustavo/eclipse_code_attacks/snova_attack/src/keys.txt", pk,
			bytes_pk, sk, bytes_sk, seed, seed_length);

	sk_gf16 sk_upk;
	sk_unpack_f(&sk_upk, sk);

	// Clear Secret!
	SNOVA_CLEAR_BYTE(&sk_upk, sizeof(sk_upk));

	uint8_t V[(v_SNOVA * lsq_SNOVA + 1)] = { 0 };
	uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1)] = { 0 };

	// generate the vinegar value
	Keccak_HashInstance hashInstance;
	Keccak_HashInitialize_SHAKE256(&hashInstance);
	Keccak_HashUpdate(&hashInstance, seed,
			8 * seed_length_private);
	Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
	Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
	Keccak_HashFinal(&hashInstance, NULL);
	Keccak_HashSqueeze(&hashInstance, vinegar_in_byte,
			8 * ((v_SNOVA * lsq_SNOVA + 1) >> 1));

	inject_faults(V, vinegar_in_byte); // Generate values for V with fault injection


	printf("vinegar_in_byte:\n");
	for(int i = 0; i < v_SNOVA * l_SNOVA * l_SNOVA; i++){
		printf("%02x,", V[i]);
	}
	printf("\n\n");
	for (int run = 0; run < NUM_RUNS; run++) {
		bool result = sign_with_fault_injection(pt_signature, digest,
				bytes_digest, array_salt, sk_upk.Aalpha, sk_upk.Balpha,
				sk_upk.Qalpha1, sk_upk.Qalpha2, sk_upk.T12, sk_upk.F11,
				sk_upk.F12, sk_upk.F21, sk_upk.pt_public_key_seed,
				sk_upk.pt_private_key_seed, V);

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
	printf("step_2_failures = %d\n", step_2_failures);
	printf("successes_r1 = %d\n", successes_r1);
	printf("successes_r2 = %d\n", successes_r2);
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

void snova_init(void) {
	static int first_time = 1;
	if (first_time) {
		first_time = 0;
		init_gf16_tables();
		gen_S_array_f();

#if OPTIMISATION != 0
        snova_plasma_init();
#endif
	}
}

int main() {
	srand(0);
	snova_init();
	int lv = v_SNOVA * l_SNOVA;
	int l = l_SNOVA;
	printf("Starting matrix P\n");
	initialize_P_matrix(epsilon, 0, l);
	printf("matrix P:\n");
	print_P_matrix(P);
	printf("calculate_Gamma_thresholds\n");
	calculate_Gamma_thresholds(lv, epsilon); // Calculate thresholds for success conditions
	printf("starting experiments\n");
	run_experiments();
	printf("done\n");
	compute_results();
	return 0;
}
