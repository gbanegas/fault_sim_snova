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
void print_c_beta(uint8_t *c_beta) {
	for (int beta = 0; beta < l_SNOVA; beta++) {
		printf("Column beta = %d:\n", beta);
		for (int x = 0; x < Q; x++) {
			printf("  Value %d: %d\n", x, c_beta[beta * Q + x]);
		}
		printf("\n");
	}

}

/*// Compute and print final results
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
 }*/

void snova_init(void) {
	static int first_time = 1;
	if (first_time) {
		first_time = 0;
		init_gf16_tables();
		gen_S_array();

#if OPTIMISATION != 0
        snova_plasma_init();
#endif
	}
}

void sign_digest_esk(uint8_t *pt_signature, const uint8_t *digest,
		uint64_t bytes_digest, uint8_t *array_salt, const uint8_t *esk) {
	snova_init();
	sk_gf16 sk_upk;
	sk_unpack(&sk_upk, esk);

	uint8_t c_beta[l_SNOVA * Q ] = { 0 }; // Counter for occurrences of elements in V columns
	for (int i = 0; i < 100; i++) {
		uint8_t V[(v_SNOVA * lsq_SNOVA + 1) >> 1] = { 0 };
		int success = sign_with_fault_injection(pt_signature, digest,
				bytes_digest, array_salt, sk_upk.Aalpha, sk_upk.Balpha,
				sk_upk.Qalpha1, sk_upk.Qalpha2, sk_upk.T12, sk_upk.F11,
				sk_upk.F12, sk_upk.F21, sk_upk.pt_public_key_seed,
				sk_upk.pt_private_key_seed, V);

		if (!success) {
			compute_stats(V, c_beta);
			printf("Gamma_epsilon_1_sup = %d\n", Gamma_epsilon_1_sup);
			printf("Gamma_epsilon_1_inf = %d\n", Gamma_epsilon_1_inf);
			printf("Gamma_epsilon_2_sup = %d\n", Gamma_epsilon_2_sup);
			printf("Gamma_epsilon_2_inf = %d\n", Gamma_epsilon_2_inf);

			for (int beta = 0; beta < l_SNOVA; beta++) {
				for (uint8_t x = 0; x < Q; x++) {
					// Check for success condition r=1
					//printf("c_beta[%d * Q + %d] = %d\n", beta, x, c_beta[beta * Q + x]);
					if (c_beta[beta * Q + x] <= Gamma_epsilon_1_sup
							&& c_beta[beta * Q + x] >= Gamma_epsilon_1_inf) {
						successes_r1++;
						successes_r2++;  // Count for both r=1 and r=2
						goto next_run;
						// Exit current loop iteration if success condition met
					}
					// Check for success condition r=2
					if (c_beta[beta * Q + x] <= Gamma_epsilon_2_sup
							&& c_beta[beta * Q + x] >= Gamma_epsilon_2_inf) {
						successes_r2++;  // Count for r=2 only
						goto next_run;
						// Exit current loop iteration if success condition met
					}
				}
			}
			next_run: continue;

		}
	}

	// Clear Secret!
	SNOVA_CLEAR_BYTE(&sk_upk, sizeof(sk_upk));
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

int main() {
	srand(0);
	snova_init();

	printf("Starting matrix P\n");
	initialize_P_matrix(0, 2);
	printf("matrix P:\n");
	print_P_matrix(P);
	printf("calculate_Gamma_thresholds\n");
	int lv = l_SNOVA * v_SNOVA;
	calculate_Gamma_thresholds(lv, epsilon); // Calculate thresholds for success conditions
	printf("starting experiments\n");

	uint8_t array_digest[64];
	uint8_t array_signature1[bytes_signature + bytes_salt];
	uint8_t array_signature2[bytes_signature + bytes_salt];

	uint8_t seed[seed_length];
	uint8_t *pt_private_key_seed;
	uint8_t *pt_public_key_seed;
	uint8_t pk[bytes_pk], sk[bytes_sk];
	uint8_t array_salt[bytes_salt];

	uint8_t entropy_input[48];
	for (int i = 0; i < 48; i++) {
		entropy_input[i] = i;
	}
	randombytes_init(entropy_input, NULL, 256);
	/*
	 randombytes(seed, seed_length);


	 randombytes(array_salt, bytes_salt);
	 */

	/**private key seed (32 bytes):
	 767f2e24cc2bc479d09d86dc9abcfde7056a8c266f9ef97ed08541dbd2e1ffa1
	 public key seed (16 bytes):
	 061550234d158c5ec95595fe04ef7a25
	 **/
	/*printf("generate_keys_pack\n");
	 generate_keys_esk(pk, sk, pt_public_key_seed, pt_private_key_seed);*/

	read_from_file(
			"/home/gustavo/eclipse_code_attacks/snova_attack/src/keys.txt", pk,
			bytes_pk, sk, bytes_sk, seed, seed_length);

	pt_public_key_seed = seed;
	pt_private_key_seed = seed + seed_length_public;

	printf("private key seed (%d bytes): \n", seed_length_private);
	print_byte(pt_private_key_seed, seed_length_private);
	printf("public key seed (%d bytes): \n", seed_length_public);
	print_byte(pt_public_key_seed, seed_length_public);

	/*printf("private key size: (%d bytes): \n", bytes_sk);
	 print_byte(sk, bytes_sk);
	 printf("public key size: (%d bytes): \n", bytes_pk);
	 print_byte(pk, bytes_pk);*/

	/*	printf("Message: \n");
	 randombytes(array_digest, 64);
	 print_byte(array_digest, 64);
	 printf("=======================\n");*/
	sign_digest_esk(array_signature1, array_digest, 64, array_salt, sk);

	/*//save_signatures("signatures.bin", array_signature1, bytes_signature + bytes_salt, array_digest, 64);
	 printf("signature (%d byte): \n", bytes_signature + bytes_salt);
	 print_byte(array_signature1, bytes_signature + bytes_salt);
	 printf("=======================\n");*/
	/*run_exp();*/
	//run_experiments();
	compute_results();

	printf("done\n");
	return 0;
}
