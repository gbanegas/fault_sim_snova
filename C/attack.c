#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "deriv_params.h"
#include "snova.h"
#include "util/util.h"


void save_to_file(const char *filename, uint8_t *pk, size_t nr_bytes_pk, uint8_t *sk, size_t nr_bytes_sk, uint8_t *seed_, size_t nr_seed_bytes) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }
    // Write pk
    fwrite(pk, sizeof(uint8_t), nr_bytes_pk, file);
    // Write sk
    fwrite(sk, sizeof(uint8_t), nr_bytes_sk, file);
    //seed
    fwrite(seed_, sizeof(uint8_t), nr_seed_bytes, file);

    fclose(file);
}


void save_signatures(const char *filename, uint8_t *signature, size_t nr_sig_plus_salt, uint8_t *message, size_t nr_message) {
    FILE *file = fopen(filename, "ab");  // Open file in append mode
    if (file == NULL) {
        perror("Error opening file for appending");
        exit(EXIT_FAILURE);
    }

    // Append the signature
    fwrite(signature, sizeof(uint8_t), nr_sig_plus_salt, file);
    // Append the corresponding message
    fwrite(message, sizeof(uint8_t), nr_message, file);

    fclose(file);

}

void read_from_file(const char *filename, uint8_t *pk, size_t nr_bytes_pk, uint8_t *sk, size_t nr_bytes_sk, uint8_t *seed_, size_t nr_seed_bytes) {
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

int main() {
    snova_init();
    uint8_t array_digest[64];
    uint8_t array_signature1[bytes_signature + bytes_salt];
    uint8_t array_signature2[bytes_signature + bytes_salt];

    uint8_t seed[seed_length];
    uint8_t* pt_private_key_seed;
    uint8_t* pt_public_key_seed;
    uint8_t pk[bytes_pk], sk[bytes_sk];
    uint8_t array_salt[bytes_salt];

    uint8_t entropy_input[48];
    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }
    randombytes_init(entropy_input, NULL, 256);
    randombytes(seed, seed_length);

    pt_public_key_seed = seed;
    pt_private_key_seed = seed + seed_length_public;

	randombytes(array_salt, bytes_salt);


    /**private key seed (32 bytes):
    767f2e24cc2bc479d09d86dc9abcfde7056a8c266f9ef97ed08541dbd2e1ffa1
    public key seed (16 bytes):
    061550234d158c5ec95595fe04ef7a25
    **/
    /*printf("generate_keys_pack\n");
    generate_keys_esk(pk, sk, pt_public_key_seed, pt_private_key_seed);*/

    printf("private key seed (%d bytes): \n", seed_length_private);
    print_byte(pt_private_key_seed, seed_length_private);
    printf("public key seed (%d bytes): \n", seed_length_public);
    print_byte(pt_public_key_seed, seed_length_public);

    read_from_file("keys.txt",pk, bytes_pk, sk, bytes_sk, seed, seed_length);

    /*printf("private key size: (%d bytes): \n", bytes_sk);
    print_byte(sk, bytes_sk);
    printf("public key size: (%d bytes): \n", bytes_pk);
    print_byte(pk, bytes_pk);*/



    for(int i = 0; i < 1;i++){
        printf("Message: \n");
        randombytes(array_digest, 64);
        print_byte(array_digest, 64);
        printf("=======================\n");
        sign_digest_esk(array_signature1, array_digest, 64, array_salt, sk);

        //save_signatures("signatures.bin", array_signature1, bytes_signature + bytes_salt, array_digest, 64);
        printf("signature (%d byte): \n", bytes_signature + bytes_salt);
        print_byte(array_signature1, bytes_signature + bytes_salt);
        printf("=======================\n");
    }


    int r = verify_signture(array_digest, 64, array_signature1, pk);

    if (r == 0) {
        printf("verification successful!\n");
    } else {
        printf("verification failed! err = %d\n", r);
    }

    /*printf("\nsign_digest_by_seed: \n");
    printf("=======================\n");
    sign_digest_ssk(array_signature2, array_digest, 64, array_salt, seed);
    printf("signature (%d byte): \n", bytes_signature + bytes_salt);
    print_byte(array_signature2, bytes_signature + bytes_salt);

    r = verify_signture(array_digest, 64, array_signature2, pk);
    if (r == 0) {
        printf("verification successful!\n");
    } else {
        printf("verification failed! err = %d\n", r);
    }*/

    return 0;
}
