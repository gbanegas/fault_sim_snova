/*
 * snova_utils.h
 *
 *  Created on: Oct 27, 2024
 *      Author: gustavo
 */

#ifndef SNOVA_UTILS_H_
#define SNOVA_UTILS_H_


#include "gf16_matrix_inline.h"
#include "snova.h"
#include "ct_functions.h"
#include "aes/snova_aes.h"

void convert_bytes_to_GF16s(const uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s);

#endif /* SNOVA_UTILS_H_ */
