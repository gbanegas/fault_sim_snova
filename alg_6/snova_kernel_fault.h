/*
 * snova_kernel_fault.h
 *
 *  Created on: Oct 27, 2024
 *      Author: gustavo
 */

#ifndef SNOVA_KERNEL_FAULT_H_
#define SNOVA_KERNEL_FAULT_H_

#ifndef SNOVA_KNL_H
#define SNOVA_KNL_H

#include "gf16_matrix_inline.h"
#include "snova.h"
#include "ct_functions.h"
#include "aes/snova_aes.h"

#include "snova_utils.h"

#include "shake/snova_shake.h"

#include "deriv_params.h"
#include "gf16_matrix.h"
#include "rng.h"


float epsilon = 0.95; // Lower bound for P_ijk
// Track counters for success and failure
int step_2_failures = 0;
// Parameters for fault probability and the signing matrix
float P[v_SNOVA][l_SNOVA][l_SNOVA]; // Probability matrix


static gf16m_t S[l_SNOVA] = {0};
static uint32_t xS[l_SNOVA][lsq_SNOVA] = {0};
static int S_is_init = 0;
// Gamma thresholds
int Gamma_epsilon_1_sup, Gamma_epsilon_1_inf;
int Gamma_epsilon_2_sup, Gamma_epsilon_2_inf;
// GF[x]/(x^4+x+1) reduction
static inline uint32_t gf16_reduce(uint32_t idx)
{
    uint32_t res, upper;

    res = idx & 0x49249249; // Octal 0o11111111111
    upper = idx >> 12;
    res = res ^ upper ^ (upper << 3);
    upper = res >> 12;
    res = res ^ upper ^ (upper << 3);
    upper = res >> 12;
    res = res ^ upper ^ (upper << 3);

    return res & 0x249;
}

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

// Conversion 4 bit -> 32 bit representation
static inline uint32_t gf16_from_nibble_f(uint8_t idx)
{
    uint32_t middle = idx | idx << 4;
    return (middle & 0x41) | ((middle << 2) & 0x208);
}

// Conversion 32 bit -> 4 bit representation
static inline uint8_t gf16_to_nibble_f(uint32_t idx)
{
    uint32_t res = gf16_reduce(idx);
    res = res | (res >> 4);
    return (res & 0x5) | ((res >> 2) & 0xa);
}

// Conversion 32 bit -> 4 bit representation
static inline uint8_t xgf16_to_nibble_f(uint32_t res)
{
    res = res | (res >> 4);
    return (res & 0x5) | ((res >> 2) & 0xa);
}

// Constant time GF16 inverse
// x^16 - x = 0 implies x^14 = x^-1
static inline uint32_t gf16_inv_f(uint32_t val)
{
    val = gf16_reduce(val);
    uint32_t res2 = gf16_reduce(val * val);
    uint32_t res4 = gf16_reduce(res2 * res2);
    uint32_t res8 = gf16_reduce(res4 * res4);

    return gf16_reduce(res2 * ((res4 * res8) & 0x49249249));
}


/**
 * Generate elements of F16[S]
 */
void gen_S_array_f(void) {
    if (S_is_init) {
        return;
    }

    S_is_init = 1;
    be_aI(S[0], 1);
    be_the_S(S[1]);
    for (int index = 2; index < l_SNOVA; index++) {
        gf16m_mul(S[index - 1], S[1], S[index]);
    }

    for (int index = 0; index < l_SNOVA; index++)
        for (int ij = 0; ij < lsq_SNOVA; ij++)
            xS[index][ij] = gf16_from_nibble_f(S[index][ij]);
}

/**
 * shake256_f
 * @param pt_seed_array - Pointer to the hash input.
 * @param input_bytes - hash lenth.
 * @param pt_output_array - Pointer to the hash output.
 * @param output_bytes - hash input.
 */
void shake256_f(const uint8_t* pt_seed_array, int input_bytes, uint8_t* pt_output_array, int output_bytes) {
    Keccak_HashInstance hashInstance;
    Keccak_HashInitialize_SHAKE256(&hashInstance);
    Keccak_HashUpdate(&hashInstance, pt_seed_array, 8 * input_bytes);
    Keccak_HashFinal(&hashInstance, NULL);
    Keccak_HashSqueeze(&hashInstance, pt_output_array, 8 * output_bytes);
}

/**
 * Convert one byte of data to GF16 representation (using only half of the
 * byte). Example: <bytes 12 34 56 78 9a bc> -> <bytes 02 01 04 03 05 ..... 0c
 * 0b>
 * @param byte_array - input (bytes)
 * @param gf16_array - output (GF16)
 * @param num_of_GF16s - GF16 amount
 */
// void convert_bytes_to_GF16s(const uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s) {
//     uint16_t* GF16_array_2B = (uint16_t*)gf16_array;
//     for (int index = 0; index < num_of_GF16s >> 1; index++) {
//         *(GF16_array_2B++) = (byte_array[index] & 0xF) ^ ((byte_array[index] & 0xF0) << 4);
//     }
//     if (num_of_GF16s % 2 == 1) {
//         gf16_array[num_of_GF16s - 1] = byte_array[num_of_GF16s >> 1] & 0xF;
//     }
// }

void convert_bytes_to_GF16s_f(const uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s) {
    int i;
    int pairs = num_of_GF16s >> 1;

    // Convert each byte into two GF16 values
    for (i = 0; i < pairs; ++i) {
        gf16_array[i * 2] = byte_array[i] & 0x0F;
        gf16_array[i * 2 + 1] = (byte_array[i] >> 4) & 0x0F;
    }

    // Handle the last GF16 value if num_of_GF16s is odd
    if (num_of_GF16s % 2 == 1) {
        gf16_array[num_of_GF16s - 1] = byte_array[pairs] & 0x0F;
    }
}


/**
 * Convert two GF16 values to one byte.
 * Example:
 *  <bytes 02 01 04 03 05 ..... 0c 0b> -> <bytes 12 34 56 78 9a bc>
 * @param byte_array - output (bytes)
 * @param gf16_array - input (GF16)
 * @param num_of_GF16s - GF16 amount
 */
// void convert_GF16s_to_bytes(uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s) {
//     for (int index = 0; index < num_of_GF16s / 2; index++) {
//         byte_array[index] = gf16_array[index << 1] | (gf16_array[(index << 1) + 1] << 4);
//     }
//     if (num_of_GF16s % 2 == 1) byte_array[num_of_GF16s / 2] = gf16_array[num_of_GF16s - 1];
// }
void convert_GF16s_to_bytes_f(uint8_t* byte_array, const gf16_t* gf16_array, int num_of_GF16s) {
    int i;
    int pairs = num_of_GF16s >> 1;

    // Convert pairs of GF16 values into one byte
    for (i = 0; i < pairs; ++i) {
        byte_array[i] = gf16_array[i * 2] | (gf16_array[i * 2 + 1] << 4);
    }

    // Handle the last GF16 value if num_of_GF16s is odd
    if (num_of_GF16s % 2 == 1) {
        byte_array[pairs] = gf16_array[num_of_GF16s - 1];
    }
}

/**
 * pk expand from seed
 *
 * Using AES-CTR encryption as a hash function
 * AES ciphertext padded with zeros.
 * The iv is also padded with zeros.
 * Using input value as the AES key.
 * The ciphertext obtained from AES encryption serves as the output of the hash
 * function.
 * @param pt_public_key_seed - Pointer to the hash input. (Fixed length of 16)
 * @param out_pk - Pointer to the hash output. (Fixed length of
 * bytes_prng_public)
 */
void pk_expand_f(const uint8_t* pt_public_key_seed, uint8_t* out_pk) {
#if PK_EXPAND_SHAKE
    uint64_t pk_bytes[(bytes_prng_public + 7) / 8] __attribute__((aligned(32)));
    snova_shake(pt_public_key_seed, 16, pk_bytes, 8 * ((bytes_prng_public + 7) / 8));
    memcpy(out_pk, pk_bytes, bytes_prng_public);
#else
    AES_128_CTR(out_pk, bytes_prng_public, pt_public_key_seed, 16);
#endif
}

/**
 * Convert one byte of data to GF16 representation (using only half of the
 * byte). cut_in_half Example: <bytes 12 34 56 78 9a bc> -> <bytes 02 04 06 08
 * 0a 0c 01 03 05 07 09 0b>
 * @param byte_array - input (bytes)
 * @param gf16_array - output (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_bytes_to_GF16s_cut_in_half_f(const uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s) {
    int half_GF16s = (num_of_GF16s + 1) >> 1;
    int i;

    // Extract the lower 4 bits of each byte to the first half of gf16_array
    for (i = 0; i < half_GF16s; ++i) {
        gf16_array[i] = byte_array[i] & 0x0F;
    }

    // Extract the upper 4 bits of each byte to the second half of gf16_array
    for (i = 0; i < (num_of_GF16s >> 1); ++i) {
        gf16_array[i + half_GF16s] = byte_array[i] >> 4;
    }
}

/**
 * Convert two GF16 values to one byte.
 * Example:
 *  <bytes 02 04 06 08 0a 0c 01 03 05 07 09 0b> -> <bytes 12 34 56 78 9a bc>
 * @param byte_array - output (bytes)
 * @param gf16_array - input (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_GF16s_to_bytes_merger_in_half_f(uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s) {
    int half_GF16s = (num_of_GF16s + 1) >> 1;
    int i;

    // Combine pairs of GF16 values into one byte
    for (i = 0; i < (num_of_GF16s >> 1); ++i) {
        byte_array[i] = gf16_array[i] | (gf16_array[i + half_GF16s] << 4);
    }

    // If num_of_GF16s is odd, handle the last GF16 value separately
    if (num_of_GF16s & 1) {
        byte_array[i] = gf16_array[i];
    }
}

/**
 * @param c - output
 * @param pt_matrix - input
 */
void gen_a_FqS_f(gf16_t* c, gf16m_t pt_matrix) {
    gf16m_t temp;
    be_aI(pt_matrix, c[0]);
    for (int i = 1; i < rank - 1; ++i) {
        gf16m_scale(S[i], c[i], temp);
        gf16m_add(pt_matrix, temp, pt_matrix);
    }
    gf16m_scale(S[rank - 1], (c[rank - 1] != 0) ? c[rank - 1] : 16 - (c[0] + (c[0] == 0)), temp);
    gf16m_add(pt_matrix, temp, pt_matrix);
    SNOVA_CLEAR(temp);
}

// Constant time version of gen_a_FqS_f
void gen_a_FqS_ct_f(gf16_t* c, gf16m_t pt_matrix) {
    uint32_t xTemp[lsq_SNOVA] = {0};
    uint32_t cX = gf16_from_nibble_f(c[0]);

    for (int ij = 0; ij < l_SNOVA; ij++)
        xTemp[ij * l_SNOVA + ij] = cX;

    for (int i1 = 1; i1 < l_SNOVA - 1; i1++) {
        cX = gf16_from_nibble_f(c[i1]);
        for (int ij = 0; ij < lsq_SNOVA; ij++)
            xTemp[ij] ^= cX * xS[i1][ij];
    }

    uint8_t zero = ct_gf16_is_not_zero(c[rank - 1]);
    uint8_t val = zero * c[rank - 1] + (1 - zero) * (15 + ct_gf16_is_not_zero(c[0]) - c[0]);

    cX = gf16_from_nibble_f(val);
    for (int ij = 0; ij < lsq_SNOVA; ij++)
        xTemp[ij] ^= cX * xS[l_SNOVA - 1][ij];

    for (int ij = 0; ij < lsq_SNOVA; ij++)
        pt_matrix[ij] = gf16_to_nibble_f(xTemp[ij]);

    SNOVA_CLEAR(xTemp);
}

/**
 * Generate the linear map T12
 * @param T12 - output
 * @param seed - input
 */
void gen_seeds_and_T12_f(T12_t T12, const uint8_t* seed) {
    gf16_t* pt_array;
    uint8_t prng_output_private[bytes_prng_private];
    gf16_t GF16_prng_output_private[GF16s_prng_private];

    shake256_f(seed, seed_length_private, prng_output_private, bytes_prng_private);
    convert_bytes_to_GF16s_f(prng_output_private, GF16_prng_output_private, GF16s_prng_private);

    // cout << "T12======" << endl;
    pt_array = GF16_prng_output_private;
    for (int j = 0; j < v_SNOVA; ++j) {
        for (int k = 0; k < o_SNOVA; ++k) {
            gen_a_FqS_ct_f(pt_array, T12[j][k]);
            // --- T12[j][k].printout();
            // printf("T12[%d][%d] = \n", j, k);
            // print_gf16m(key->map2.T12[j][k]);
            pt_array += l_SNOVA;
        }
    }

    // Clear Secret!
    SNOVA_CLEAR(prng_output_private);
    SNOVA_CLEAR(GF16_prng_output_private);
}

/**
 * Generate the random part of public key
 * @param map - P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param pt_public_key_seed - input
 */
void gen_A_B_Q_P_f(map_group1* map, const uint8_t* pt_public_key_seed) {
    uint8_t prng_output_public[((GF16s_prng_public + 1) >> 1)];
    uint8_t temp[lsq_SNOVA * l_SNOVA];
    // ----- pt temp -----
    pk_expand_f(pt_public_key_seed, prng_output_public);
    convert_bytes_to_GF16s_f(prng_output_public, (uint8_t*)map, GF16s_prng_public);
    memcpy(temp, ((uint8_t*)map->Qalpha1), lsq_SNOVA * l_SNOVA);

    // printf("Aalpha: \n");
    for (int index = 0; index < lsq_SNOVA; ++index) {
        be_invertible_by_add_aS(map->Aalpha[index]);
        // print_byte(map->Aalpha[index], sq_rank);
        // Aalpha[index].printout();
    }

    // cout << "B======" << endl;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        be_invertible_by_add_aS(map->Balpha[index]);
        // Balpha[index].printout();
    }
    gf16_t* pt_array = (uint8_t*)(map->Qalpha1) + l_SNOVA * lsq_SNOVA;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        gen_a_FqS_f(pt_array, map->Qalpha2[index]);
        // be_invertible_by_add_aS(map->Qalpha2[index]);
        pt_array += l_SNOVA;
        // Qalpha2[index].printout();
    }

    pt_array = temp;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        gen_a_FqS_f(pt_array, map->Qalpha1[index]);
        // be_invertible_by_add_aS(map->Qalpha1[index]);
        // Qalpha1[index].printout();
        pt_array += l_SNOVA;
    }
}

/**
 * Generate private key (F part)
 * @param map2 - output: F11 F12 F21
 * @param map1 - input: P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param T12 - input
 */
void gen_F_ref_f(map_group2* map2, map_group1* map1, T12_t T12) {
    gf16m_t temp;
    // for (int i = 0; i < m_SNOVA; ++i) {
    //     for (int j = 0; j < v_SNOVA; ++j) {
    //         for (int k = 0; k < v_SNOVA; ++k) {
    //             gf16m_clone(map2->F11[i][j][k], map1->P11[i][j][k]);
    //         }
    //     }
    // }
    memcpy(map2->F11, map1->P11, m_SNOVA * v_SNOVA * v_SNOVA * sq_rank);
    memcpy(map2->F12, map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * sq_rank);
    memcpy(map2->F21, map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * sq_rank);

    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < v_SNOVA; ++j) {
            for (int k = 0; k < o_SNOVA; ++k) {
                // gf16m_clone(map2->F12[i][j][k], map1->P12[i][j][k]);
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(map1->P11[i][j][index], T12[index][k], temp);
                    gf16m_add(map2->F12[i][j][k], temp, map2->F12[i][j][k]);
                    // F12[i][j][k].printout();
                }
            }
        }
    }

    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < o_SNOVA; ++j) {
            for (int k = 0; k < v_SNOVA; ++k) {
                // gf16m_clone(map2->F21[i][j][k], map1->P21[i][j][k]);
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(T12[index][j], map1->P11[i][index][k], temp);
                    gf16m_add(map2->F21[i][j][k], temp, map2->F21[i][j][k]);
                    // F21[i][j][k].printout();
                }
            }
        }
    }

    // Clear Secret!
    SNOVA_CLEAR(temp);
}

/**
 * Generate public key (P22 part)
 * @param outP22 - output
 * @param T12 - input
 * @param P21 - input
 * @param F12 - input
 */
void gen_P22_ref_f(P22_byte_t outP22, T12_t T12, P21_t P21, F12_t F12) {
    gf16m_t temp1, temp2;
    P22_t P22 = {0};
    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < o_SNOVA; ++j) {
            for (int k = 0; k < o_SNOVA; ++k) {
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(T12[index][j], F12[i][index][k], temp1);
                    gf16m_mul(P21[i][j][index], T12[index][k], temp2);
                    gf16m_add(temp1, temp2, temp1);
                    gf16m_add(P22[i][j][k], temp1, P22[i][j][k]);
                    // P22[i][j][k].printout();
                }
            }
        }
    }

    convert_GF16s_to_bytes_f(outP22, (uint8_t*)P22, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);

    // Clear Secret!
    SNOVA_CLEAR(temp1);
    SNOVA_CLEAR(temp2);
}

/**
 * P22 byte to GF16
 * @param P22_gf16s - output
 * @param P22_bytes - input
 */
void input_P22_f(uint8_t* P22_gf16s, const uint8_t* P22_bytes) {
    convert_bytes_to_GF16s_f(P22_bytes, P22_gf16s, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);
}

/**
 * Pack expanded private key. esk = (key_elems, pt_private_key_seed).
 * @param esk - pointer to output expanded private key.
 * @param key_elems - pointer to input snova key elements.
 * @param pt_private_key_seed - pointer to input private key seed.
 */
void sk_pack_f(uint8_t* esk, snova_key_elems* key_elems, const uint8_t* pt_private_key_seed) {
    uint8_t* sk_gf16_ptr = (uint8_t*)(key_elems->map1.Aalpha);
    convert_GF16s_to_bytes_merger_in_half_f(esk, sk_gf16_ptr, (bytes_sk - (seed_length_public + seed_length_private)) * 2);
    memcpy(esk + (bytes_sk - (seed_length_public + seed_length_private)), key_elems->pk.pt_public_key_seed, seed_length_public);
    memcpy(esk + (bytes_sk - seed_length_private), pt_private_key_seed, seed_length_private);
}

/**
 * Unpack expanded secret key. skupk = (esk).
 * @param skupk - pointer to output private key (unpack).
 * @param esk - pointer to input expanded private key.
 */
void sk_unpack_f(sk_gf16* skupk, const uint8_t* esk) {
    convert_bytes_to_GF16s_cut_in_half_f(esk, (uint8_t*)skupk, (bytes_sk - (seed_length_public + seed_length_private)) * 2);
    memcpy(skupk->pt_public_key_seed, esk + (bytes_sk - (seed_length_public + seed_length_private)),
           seed_length_public + seed_length_private);
}

/**
 * Pack public key. pk = (key_elems).
 */
void pk_pack_f(uint8_t* pk, snova_key_elems* key_elems) {
    memcpy(pk, &key_elems->pk, bytes_pk);
}

/**
 * createHashOut
 */
void createSignedHash_f(const uint8_t* digest, uint64_t bytes_digest, const uint8_t* pt_public_key_seed,
                      const uint8_t* array_salt, uint8_t* signed_hash_out) {
    Keccak_HashInstance hashInstance;
    Keccak_HashInitialize_SHAKE256(&hashInstance);
    Keccak_HashUpdate(&hashInstance, pt_public_key_seed, 8 * seed_length_public);
    Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
    Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
    Keccak_HashFinal(&hashInstance, NULL);
    Keccak_HashSqueeze(&hashInstance, signed_hash_out, 8 * bytes_hash);
}

/**
 * Verifies signature.
 * @param pt_digest - pointer to input digest.
 * @param pt_signature - pointer to output signature.
 * @param pk - pointer to output public key.
 * @returns - 0 if signature could be verified correctly and -1 otherwise
 */
int verify_signture_ref_f(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature, const uint8_t* pk) {
    uint8_t hash_in_bytes[bytes_hash];
    uint8_t signed_hash[bytes_hash];
    const uint8_t* pt_salt = pt_signature + bytes_signature;

    gf16m_t Left[lsq_SNOVA][n_SNOVA], Right[lsq_SNOVA][n_SNOVA];
    gf16m_t hash_in_GF16Matrix[m_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA];
    gf16m_t P22[m_SNOVA][o_SNOVA][o_SNOVA];

    map_group1 map1;
    gf16m_t temp1, temp2;

    public_key* pk_stru = (public_key*)pk;

    Keccak_HashInstance hashInstance;
    Keccak_HashInitialize_SHAKE256(&hashInstance);
    Keccak_HashUpdate(&hashInstance, pk_stru->pt_public_key_seed, 8 * seed_length_public);
    Keccak_HashUpdate(&hashInstance, pt_digest, 8 * bytes_digest);
    Keccak_HashUpdate(&hashInstance, pt_salt, 8 * bytes_salt);
    Keccak_HashFinal(&hashInstance, NULL);
    Keccak_HashSqueeze(&hashInstance, signed_hash, 8 * bytes_hash);

#if (o_SNOVA * l_SNOVA) & 0x1 == 1
    signed_hash[bytes_hash - 1] &= 0x0f;
#endif

    convert_bytes_to_GF16s_f(pt_signature, (gf16_t*)signature_in_GF16Matrix, GF16s_signature);
    // generate PRNG part of public key
    gen_A_B_Q_P_f(&map1, pk_stru->pt_public_key_seed);
    // read  P22
    input_P22_f((uint8_t*)P22, pk_stru->P22);
    // evaluate signature GF16Matrix array
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
        for (int index = 0; index < n_SNOVA; ++index) {
            // Left[alpha][index]
            gf16m_transpose(signature_in_GF16Matrix[index], temp1);
            gf16m_mul(temp1, map1.Qalpha1[alpha], temp2);
            gf16m_mul(map1.Aalpha[alpha], temp2, Left[alpha][index]);
            // Right[alpha][index]
            gf16m_mul(map1.Qalpha2[alpha], signature_in_GF16Matrix[index], temp2);
            gf16m_mul(temp2, map1.Balpha[alpha], Right[alpha][index]);

            /*
            Left[alpha][index] = Aalpha[alpha] *
            (signature_in_GF16Matrix[index].transpose()) * Qalpha1[alpha];
            Right[alpha][index] = Qalpha2[alpha] *
            signature_in_GF16Matrix[index] * Balpha[alpha];
            */
        }
    }

    for (int i = 0; i < m_SNOVA; ++i) {
        // gf16m_add(hash_in_GF16Matrix[i], hash_in_GF16Matrix[i],
        // hash_in_GF16Matrix[i]);
        gf16m_set_zero(hash_in_GF16Matrix[i]);
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
            for (int dj = 0; dj < v_SNOVA; ++dj) {
                for (int dk = 0; dk < v_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj], map1.P11[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2, hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj] * P11[i][dj][dk] * Right[alpha][dk];
                }
            }
            for (int dj = 0; dj < v_SNOVA; ++dj) {
                for (int dk = 0; dk < o_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj], map1.P12[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk + v_SNOVA], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2, hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj] * P12[i][dj][dk] *
                    // Right[alpha][dk+v_SNOVA];
                }
            }
            for (int dj = 0; dj < o_SNOVA; ++dj) {
                for (int dk = 0; dk < v_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj + v_SNOVA], map1.P21[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2, hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj + v_SNOVA] * P21[i][dj][dk] *
                    // Right[alpha][dk];
                }
            }
            for (int dj = 0; dj < o_SNOVA; ++dj) {
                for (int dk = 0; dk < o_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj + v_SNOVA], P22[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk + v_SNOVA], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2, hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj + v_SNOVA] * P22[i][dj][dk] *
                    // Right[alpha][dk + v_SNOVA];
                }
            }
        }
    }

    // GF16Matrix array ==>  GF16 array ==> byte array
    for (int index = 0; index < m_SNOVA * lsq_SNOVA; ++index) {
        ((gf16_t*)signature_in_GF16Matrix)[index] =
            get_gf16m(hash_in_GF16Matrix[index / lsq_SNOVA], (index % lsq_SNOVA) / l_SNOVA, (index % lsq_SNOVA) % l_SNOVA);
    }

    convert_GF16s_to_bytes_f(hash_in_bytes, (gf16_t*)signature_in_GF16Matrix, m_SNOVA * lsq_SNOVA);
    int result = 0;
    for (int i = 0; i < bytes_hash; ++i) {
        if (hash_in_bytes[i] != signed_hash[i]) {
            result = -1;
            break;
        }
    }

    return result;
}

void evaluation_part_f(gf16m_t *Aalpha, gf16m_t *Balpha, gf16m_t *Qalpha1, gf16m_t *Qalpha2, gf16m_t(*F11)[24][24], gf16m_t Left[16][24], gf16m_t Right[16][24], gf16m_t X_in_GF16Matrix[29], gf16m_t Fvv_in_GF16Matrix[5], gf16m_t gf16m_temp0, gf16m_t gf16m_temp1) {
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
        for (int index = 0; index < v_SNOVA; ++index) {
            gf16m_transpose(X_in_GF16Matrix[index], gf16m_temp0);
            gf16m_mul(gf16m_temp0, Qalpha1[alpha], gf16m_temp1);
            gf16m_mul(Aalpha[alpha], gf16m_temp1, Left[alpha][index]);
            gf16m_mul(Qalpha2[alpha], X_in_GF16Matrix[index], gf16m_temp1);
            gf16m_mul(gf16m_temp1, Balpha[alpha], Right[alpha][index]);
            /*
                Left[alpha][index] = Aalpha[alpha] *
                (X_in_GF16Matrix[index].transpose()) * Qalpha1[alpha];
                Right[alpha][index] = Qalpha2[alpha] * X_in_GF16Matrix[index] *
                Balpha[alpha];
                */
        }
    }
    for (int i = 0; i < m_SNOVA; ++i) {
        gf16m_set_zero(Fvv_in_GF16Matrix[i]);
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
            for (int j = 0; j < v_SNOVA; ++j) {
                for (int k = 0; k < v_SNOVA; ++k) {
                    gf16m_mul(Left[alpha][j], F11[i][j][k], gf16m_temp0);
                    gf16m_mul(gf16m_temp0, Right[alpha][k], gf16m_temp1);
                    gf16m_add(Fvv_in_GF16Matrix[i], gf16m_temp1, Fvv_in_GF16Matrix[i]);
                    /*
                        Fvv_in_GF16Matrix[i] = Fvv_in_GF16Matrix[i] +
                        Left[alpha][j] * F11[i][j][k] * Right[alpha][k];
                        */
                }
            }
        }
    }
}

void write_signatures_to_csv_f(const char *filename, gf16m_t X_in_GF16Matrix_fault[n_SNOVA][rank]) {
    FILE *fp;
    fp = fopen(filename, "w");

    if (fp == NULL) {
        printf("Error opening the file for writing!\n");
        return;
    }

    // Loop through the rows and columns of the matrix and write to the file
    for (int i = 0; i < n_SNOVA; i++) {
        for (int j = 0; j < rank; j++) {
            fprintf(fp, "%d", X_in_GF16Matrix_fault[i][j]); // Write the element

            // Add a comma after each element except the last one in a row
            if (j < rank - 1) {
                fprintf(fp, ",");
            }
        }
        // New line after each row
        fprintf(fp, "\n");
    }

    // Close the file
    fclose(fp);
    printf("Data written to %s successfully.\n", filename);
}

void add_last_column_f(gf16_t Gauss[80][81], gf16m_t Fvv_in_GF16Matrix[5]) {
    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < rank; ++j) {
            for (int k = 0; k < rank; ++k) {
                int index1 = i * lsq_SNOVA + j * rank + k;
                int index2 = m_SNOVA * lsq_SNOVA;
                Gauss[index1][index2] = gf16_get_add(Gauss[index1][index2], get_gf16m(Fvv_in_GF16Matrix[i], j, k));
            }
        }
    }
}

void clean_temp_f(gf16_t Temp[16][16]) {
    for (int ti = 0; ti < lsq_SNOVA; ++ti) {
        for (int tj = 0; tj < lsq_SNOVA; ++tj) {
            Temp[ti][tj] = 0;
        }
    }
}

void print_GF16Matrix_f(gf16m_t matrix) {
    for (int i = 0; i < sq_rank; i++) {
        printf("%d ", matrix[i]);  // Print in decimal format
        /*if ((i + 1) % rank == 0) {
            printf("\n");  // Newline after each row
        }*/
    }
    printf("\n");
}

void print_all_matrices(gf16m_t matrices[], int num_matrices) {
    for (int i = 0; i < num_matrices; i++) {
        print_GF16Matrix_f(matrices[i]);
    }
}

/**
 * Computes signature
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param Aalpha -
 * @param Balpha -
 * @param Qalpha1 -
 * @param Qalpha2 -
 * @param T12 -
 * @param F11 -
 * @param F12 -
 * @param F21 -
 * @param pt_public_key_seed - pointer to output public key seed.
 * @param pt_private_key_seed - pointer to output private key seed.
 */
void sign_digest_core_ref_f(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, uint8_t* array_salt, Aalpha_t Aalpha,
                          Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, T12_t T12, F11_t F11, F12_t F12, F21_t F21,
                          const uint8_t pt_public_key_seed[seed_length_public], const uint8_t pt_private_key_seed[seed_length_private]) {
    gf16_t Gauss[m_SNOVA * lsq_SNOVA][m_SNOVA * lsq_SNOVA + 1];
    gf16_t Temp[lsq_SNOVA][lsq_SNOVA];
    gf16_t t_GF16, solution[m_SNOVA * lsq_SNOVA];

    gf16m_t Left_X_tmp, Right_X_tmp;
    gf16_t *Left_X, *Right_X;

    gf16m_t Left[lsq_SNOVA][v_SNOVA], Right[lsq_SNOVA][v_SNOVA];
    gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
    //gf16m_t X_in_GF16Matrix_faulty[n_SNOVA] = {0};
    gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
    gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA] = {0};

    uint8_t signed_hash[bytes_hash];
    uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1];

    // temp
    gf16m_t gf16m_temp0;
    gf16m_t gf16m_temp1;
    gf16m_t gf16m_secret_temp0;

    Left_X = Left_X_tmp;
    Right_X = Right_X_tmp;
    int flag_redo = 1;
    uint8_t num_sign = 0;
    int counter;

    createSignedHash_f(digest, bytes_digest, pt_public_key_seed, array_salt, signed_hash);

    // put hash value in GF16 array
    convert_bytes_to_GF16s_f(signed_hash, hash_in_GF16, GF16s_hash);

    // printf("hash value in GF16: \n");
    // print_gf16(key->hash_in_GF16, m_SNOVA * lsq_SNOVA);
    // printf("===================\n");
    //--------------------------------------
    do {
        num_sign++;
        flag_redo = 0;
        // put hash value in the last column of Gauss matrix
        for (int index = 0; index < (m_SNOVA * lsq_SNOVA); index++) {
            Gauss[index][m_SNOVA * lsq_SNOVA] = hash_in_GF16[index];
        }
        // generate the vinegar value
        Keccak_HashInstance hashInstance;
        Keccak_HashInitialize_SHAKE256(&hashInstance);
        Keccak_HashUpdate(&hashInstance, pt_private_key_seed, 8 * seed_length_private);
        Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
        Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
        Keccak_HashUpdate(&hashInstance, &num_sign, 8);
        Keccak_HashFinal(&hashInstance, NULL);
        Keccak_HashSqueeze(&hashInstance, vinegar_in_byte, 8 * ((v_SNOVA * lsq_SNOVA + 1) >> 1));

        counter = 0;

        // Inject a fault by fixing a value in vinegar_in_byte
        // For example, fixing the first byte (or any other byte) to a constant value
        // So que have the vinegar values fixed (16 bits)
        vinegar_in_byte[0] = 0xFF;
        vinegar_in_byte[1] = 0xFF;
        for (int index = 0; index < v_SNOVA; index++) {
            for (int i = 0; i < rank; ++i) {
                for (int j = 0; j < rank; ++j) {
                    set_gf16m(X_in_GF16Matrix[index], i, j,
                              ((counter & 1) ? (vinegar_in_byte[counter >> 1] >> 4) : (vinegar_in_byte[counter >> 1] & 0xF)));
                    counter++;
                }
            }
        }

        // evaluate the vinegar part of central map
        evaluation_part_f(Aalpha, Balpha, Qalpha1, Qalpha2, F11, Left, Right, X_in_GF16Matrix, Fvv_in_GF16Matrix, gf16m_temp0,
                        gf16m_temp1);

        // add to the last column of Gauss matrix
        add_last_column_f(Gauss, Fvv_in_GF16Matrix);

        // compute the coefficients of Xo and put into Gauss matrix and compute
        // the coefficients of Xo^t and add into Gauss matrix
        for (int i = 0; i < m_SNOVA; ++i) {
            for (int index = 0; index < o_SNOVA; ++index) {
                clean_temp_f(Temp);

                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
                    for (int j = 0; j < v_SNOVA; ++j) {
                        gf16m_mul(Left[alpha][j], F12[i][j][index], gf16m_temp0);
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
                                temp3 = gf16_get_mul(get_gf16m(Left_X, ti / rank, tj / rank),
                                                     get_gf16m(Right_X, tj % rank, ti % rank));
                                Temp[ti][tj] = gf16_get_add(Temp[ti][tj], temp3);
                            }
                        }
                    }
                }

                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
                    for (int j = 0; j < v_SNOVA; ++j) {
                        Left_X = Aalpha[alpha];
                        gf16m_mul(Qalpha1[alpha], F21[i][index][j], gf16m_temp0);
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
                                temp2 = gf16_get_mul(get_gf16m(Left_X, ti / rank, tj % rank),
                                                     get_gf16m(Right_X, tj / rank, ti % rank));
                                Temp[ti][tj] = gf16_get_add(Temp[ti][tj], temp2);
                            }
                        }
                    }
                }
                for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                    for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                        Gauss[i * lsq_SNOVA + ti][index * lsq_SNOVA + tj] = Temp[ti][tj];
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
                        Gauss[j][k] = gf16_get_add(Gauss[j][k], gf16_get_mul(Gauss[i][k], t_GF16));
                    }
                }
            }
        }

        if (!flag_redo) {
            for (int i = m_SNOVA * lsq_SNOVA - 1; i >= 0; --i) {
                t_GF16 = 0;
                for (int k = i + 1; k < m_SNOVA * lsq_SNOVA; ++k) {
                    t_GF16 = gf16_get_add(t_GF16, gf16_get_mul(Gauss[i][k], solution[k]));
                }
                solution[i] = gf16_get_add(Gauss[i][m_SNOVA * lsq_SNOVA], t_GF16);
            }
        }

    } while (flag_redo);


    // printf("times of Gauss elimination : %d\n", num_sign);

    for (int index = 0; index < o_SNOVA; ++index) {
        for (int i = 0; i < rank; ++i) {
            for (int j = 0; j < rank; ++j) {
                set_gf16m(X_in_GF16Matrix[index + v_SNOVA], i, j, solution[index * lsq_SNOVA + i * rank + j]);
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
            gf16m_mul(T12[index][i], X_in_GF16Matrix[v_SNOVA + i], gf16m_secret_temp0);
            gf16m_add(signature_in_GF16Matrix[index], gf16m_secret_temp0, signature_in_GF16Matrix[index]);
            /*
            signature_in_GF16Matrix[index] = signature_in_GF16Matrix[index] +
            T12[index][i] * X_in_GF16Matrix[v_SNOVA + i];
            */
        }
    }
    for (int index = 0; index < o_SNOVA; ++index) {
        gf16m_clone(signature_in_GF16Matrix[v_SNOVA + index], X_in_GF16Matrix[v_SNOVA + index]);
    }
    printf("signature_in_GF16Matrix:\n");
    print_all_matrices(signature_in_GF16Matrix, n_SNOVA);


    // output signature
    for (int index = 0; index < n_SNOVA * lsq_SNOVA; ++index) {
        ((gf16_t*)signature_in_GF16Matrix)[index] =
            get_gf16m(signature_in_GF16Matrix[index / lsq_SNOVA], (index % lsq_SNOVA) / l_SNOVA, (index % lsq_SNOVA) % l_SNOVA);
    }
    convert_GF16s_to_bytes_f(pt_signature, (gf16_t*)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);
    for (int i = 0; i < bytes_salt; ++i) {
        pt_signature[bytes_signature + i] = array_salt[i];
    }

    SNOVA_CLEAR(gf16m_secret_temp0);
}


void sign_digest_core_ref_f_fault(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, uint8_t* array_salt, Aalpha_t Aalpha,
                          Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, T12_t T12, F11_t F11, F12_t F12, F21_t F21,
                          const uint8_t pt_public_key_seed[seed_length_public], const uint8_t pt_private_key_seed[seed_length_private]) {
    gf16_t Gauss[m_SNOVA * lsq_SNOVA][m_SNOVA * lsq_SNOVA + 1];
    gf16_t Temp[lsq_SNOVA][lsq_SNOVA];
    gf16_t t_GF16, solution[m_SNOVA * lsq_SNOVA];

    gf16m_t Left_X_tmp, Right_X_tmp;
    gf16_t *Left_X, *Right_X;

    gf16m_t Left[lsq_SNOVA][v_SNOVA], Right[lsq_SNOVA][v_SNOVA];
    gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
    //gf16m_t X_in_GF16Matrix_faulty[n_SNOVA] = {0};
    gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
    gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA] = {0};

    uint8_t signed_hash[bytes_hash];
    uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1];

    // temp
    gf16m_t gf16m_temp0;
    gf16m_t gf16m_temp1;
    gf16m_t gf16m_secret_temp0;

    Left_X = Left_X_tmp;
    Right_X = Right_X_tmp;
    int flag_redo = 1;
    uint8_t num_sign = 0;
    int counter;

    createSignedHash_f(digest, bytes_digest, pt_public_key_seed, array_salt, signed_hash);

    // put hash value in GF16 array
    convert_bytes_to_GF16s_f(signed_hash, hash_in_GF16, GF16s_hash);

    // printf("hash value in GF16: \n");
    // print_gf16(key->hash_in_GF16, m_SNOVA * lsq_SNOVA);
    // printf("===================\n");
    //--------------------------------------
    do {
        num_sign++;
        flag_redo = 0;
        // put hash value in the last column of Gauss matrix
        for (int index = 0; index < (m_SNOVA * lsq_SNOVA); index++) {
            Gauss[index][m_SNOVA * lsq_SNOVA] = hash_in_GF16[index];
        }
        // generate the vinegar value
        Keccak_HashInstance hashInstance;
           Keccak_HashInitialize_SHAKE256(&hashInstance);
           Keccak_HashUpdate(&hashInstance, pt_private_key_seed, 8 * seed_length_private);
           Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
           Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
           Keccak_HashUpdate(&hashInstance, &num_sign, 8);
           Keccak_HashFinal(&hashInstance, NULL);
           Keccak_HashSqueeze(&hashInstance, vinegar_in_byte, 8 * ((v_SNOVA * lsq_SNOVA + 1) >> 1));

           // Apply transient faults to specific parts of vinegar_in_byte based on fault injection model.
           // Here we fix the first byte of vinegar_in_byte as per pseudocode and instructions.
           vinegar_in_byte[0] = 0xFF; // Fixed byte to simulate a fault
           vinegar_in_byte[1] = 0xFF; // Another example of a fixed byte to simulate fault

           // Proceed with assigning vinegar values to the matrix elements
           counter = 0;
           for (int index = 0; index < v_SNOVA; index++) {
               for (int i = 0; i < rank; ++i) {
                   for (int j = 0; j < rank; ++j) {
                       set_gf16m(X_in_GF16Matrix[index], i, j,
                                 ((counter & 1) ? (vinegar_in_byte[counter >> 1] >> 4) : (vinegar_in_byte[counter >> 1] & 0xF)));
                       counter++;
                   }
               }
           }

        // evaluate the vinegar part of central map
        evaluation_part_f(Aalpha, Balpha, Qalpha1, Qalpha2, F11, Left, Right, X_in_GF16Matrix, Fvv_in_GF16Matrix, gf16m_temp0,
                        gf16m_temp1);

        // add to the last column of Gauss matrix
        add_last_column_f(Gauss, Fvv_in_GF16Matrix);

        // compute the coefficients of Xo and put into Gauss matrix and compute
        // the coefficients of Xo^t and add into Gauss matrix
        for (int i = 0; i < m_SNOVA; ++i) {
            for (int index = 0; index < o_SNOVA; ++index) {
                clean_temp_f(Temp);

                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
                    for (int j = 0; j < v_SNOVA; ++j) {
                        gf16m_mul(Left[alpha][j], F12[i][j][index], gf16m_temp0);
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
                                temp3 = gf16_get_mul(get_gf16m(Left_X, ti / rank, tj / rank),
                                                     get_gf16m(Right_X, tj % rank, ti % rank));
                                Temp[ti][tj] = gf16_get_add(Temp[ti][tj], temp3);
                            }
                        }
                    }
                }

                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
                    for (int j = 0; j < v_SNOVA; ++j) {
                        Left_X = Aalpha[alpha];
                        gf16m_mul(Qalpha1[alpha], F21[i][index][j], gf16m_temp0);
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
                                temp2 = gf16_get_mul(get_gf16m(Left_X, ti / rank, tj % rank),
                                                     get_gf16m(Right_X, tj / rank, ti % rank));
                                Temp[ti][tj] = gf16_get_add(Temp[ti][tj], temp2);
                            }
                        }
                    }
                }
                for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                    for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                        Gauss[i * lsq_SNOVA + ti][index * lsq_SNOVA + tj] = Temp[ti][tj];
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
                          Gauss[j][k] = gf16_get_add(Gauss[j][k], gf16_get_mul(Gauss[i][k], t_GF16));
                      }
                  }
              }
          }

        if (!flag_redo) {
            for (int i = m_SNOVA * lsq_SNOVA - 1; i >= 0; --i) {
                t_GF16 = 0;
                for (int k = i + 1; k < m_SNOVA * lsq_SNOVA; ++k) {
                    t_GF16 = gf16_get_add(t_GF16, gf16_get_mul(Gauss[i][k], solution[k]));
                }
                solution[i] = gf16_get_add(Gauss[i][m_SNOVA * lsq_SNOVA], t_GF16);
            }
        }

    } while (flag_redo);


    // printf("times of Gauss elimination : %d\n", num_sign);

    for (int index = 0; index < o_SNOVA; ++index) {
        for (int i = 0; i < rank; ++i) {
            for (int j = 0; j < rank; ++j) {
                set_gf16m(X_in_GF16Matrix[index + v_SNOVA], i, j, solution[index * lsq_SNOVA + i * rank + j]);
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
            gf16m_mul(T12[index][i], X_in_GF16Matrix[v_SNOVA + i], gf16m_secret_temp0);
            gf16m_add(signature_in_GF16Matrix[index], gf16m_secret_temp0, signature_in_GF16Matrix[index]);
            /*
            signature_in_GF16Matrix[index] = signature_in_GF16Matrix[index] +
            T12[index][i] * X_in_GF16Matrix[v_SNOVA + i];
            */
        }
    }
    for (int index = 0; index < o_SNOVA; ++index) {
        gf16m_clone(signature_in_GF16Matrix[v_SNOVA + index], X_in_GF16Matrix[v_SNOVA + index]);
    }
    printf("signature_in_GF16Matrix:\n");
    print_all_matrices(signature_in_GF16Matrix, n_SNOVA);


    // output signature
    for (int index = 0; index < n_SNOVA * lsq_SNOVA; ++index) {
        ((gf16_t*)signature_in_GF16Matrix)[index] =
            get_gf16m(signature_in_GF16Matrix[index / lsq_SNOVA], (index % lsq_SNOVA) / l_SNOVA, (index % lsq_SNOVA) % l_SNOVA);
    }
    convert_GF16s_to_bytes_f(pt_signature, (gf16_t*)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);
    for (int i = 0; i < bytes_salt; ++i) {
        pt_signature[bytes_signature + i] = array_salt[i];
    }

    SNOVA_CLEAR(gf16m_secret_temp0);
}

#define Q 16
#define NUM_RUNS 100
#define L 4 // Example for l; replace as appropriate for SNOVA
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
	} while (random_value <epsilon);
	return random_value;
}

// Fault injection with Bernoulli distribution
void inject_faults(uint8_t *V, uint8_t *vinegar_byte) {
	uint8_t fixed_value = random_F_q_excluding(0); // Example fixed fault value in F_q

	for (int i = 0; i < v_SNOVA; i++) {
		for (int j = 0; j < l_SNOVA; j++) {
			for (int k = 0; k < l_SNOVA; k++) {
				float res = (float) rand();
				if (res <= P[i][j][k]) {
					V[(i * l_SNOVA * l_SNOVA) + j * l_SNOVA + k] = fixed_value;
				} else {
					V[(i * l_SNOVA * l_SNOVA) + j * l_SNOVA + k] =vinegar_byte[i];
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


// Initialize P matrix with epsilon value
void initialize_P_matrix(float epsilon, int inf, int sup) {
	// By default, each entry can be chosen uniformly randomly.
    for (int i = 0; i < v_SNOVA; ++i) {
		for (int j = 0; j < l_SNOVA; ++j) {
			for (int k = 0; k < l_SNOVA; ++k) {
				P[i][j][k] = (float)1/Q;
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


bool generate_signature(uint8_t *pt_signature, const uint8_t *digest,
		uint64_t bytes_digest, uint8_t *array_salt, Aalpha_t Aalpha,
		Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, T12_t T12,
		F11_t F11, F12_t F12, F21_t F21,
		const uint8_t pt_public_key_seed[seed_length_public],
		const uint8_t pt_private_key_seed[seed_length_private], uint8_t *vinegar_in_byte) {

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

	int counter= 0;


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
		evaluation_part_f(Aalpha, Balpha, Qalpha1, Qalpha2, F11, Left, Right,
				X_in_GF16Matrix, Fvv_in_GF16Matrix, gf16m_temp0, gf16m_temp1);

		// add to the last column of Gauss matrix
		add_last_column_f(Gauss, Fvv_in_GF16Matrix);

		// compute the coefficients of Xo and put into Gauss matrix and compute
		// the coefficients of Xo^t and add into Gauss matrix
		for (int i = 0; i < m_SNOVA; ++i) {
			for (int index = 0; index < o_SNOVA; ++index) {
				clean_temp_f(Temp);

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
    //only one iteration
	} while (attemps <= 2);
	if (flag_redo) {
		// it means the Gauss elimnation has failed.
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
	convert_GF16s_to_bytes_f(pt_signature, (gf16_t*) signature_in_GF16Matrix,
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
		const uint8_t pt_private_key_seed[seed_length_private], uint8_t *V) {


	gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];

	uint8_t signed_hash[bytes_hash];


	bool signature_success = false;
	createSignedHash_f(digest, bytes_digest, pt_public_key_seed, array_salt,
			signed_hash);
	convert_bytes_to_GF16s_f(signed_hash, hash_in_GF16, GF16s_hash);



	do {

		signature_success = false;



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

#endif
#endif /* SNOVA_KERNEL_FAULT_H_ */
