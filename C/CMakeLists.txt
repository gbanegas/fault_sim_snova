cmake_minimum_required(VERSION 3.0)
project(snova C ASM)

set(CMAKE_C_STANDARD 99)

add_executable(snova_test test/test.c ct_functions.c gf16.h gf16_init.h gf16_matrix.h gf16_matrix_inline.h params.h
              shake/KeccakHash.c shake/KeccakSponge.c shake/SimpleFIPS202.c shake/KeccakP-1600-opt64.c shake/snova_shake_ref.c
            shake/snova_shake_ref.c
        aes/aes_c.c aes/aes_impl.c aes/aes_ossl.c
        aes/aes.c aes/snova_aes.c rng.c snova.c sign.c wasmapi.c)
#target_link_libraries(snova_test OpenSSL::Crypto)
#target_include_directories(snova_test PRIVATE .)