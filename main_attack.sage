load("utils.sage")
load("snova_faulty.sage")
load("attack.sage")


q = 16
F = GF(q, 'x')
FE = []
dict_for_S = {}



# Store all the elements from GF(16) in the list FE.
for i in range(16):
    FE.append(F.from_integer(i))

# Construct the matrices S2,S3,S4 (l=2,3,4)
S2 = matrix(F, [[FE[8], FE[7]], [FE[7], FE[6]]])
S3 = matrix(F, [[FE[8], FE[7], FE[6]], [FE[7], FE[6], FE[5]], [FE[6], FE[5], FE[4]]])
S4 = matrix(F, [[FE[8], FE[7], FE[6], FE[5]], [FE[7], FE[6], FE[5], FE[4]], [FE[6], FE[5], FE[4], FE[3]],
                [FE[5], FE[4], FE[3], FE[2]]])
# Define a dictionary to easy access to S2,S3, or S4.
dict_for_S = {2: S2, 3: S3, 4: S4}
SL = [(28, 17, F, 2), (25, 8, F, 3), (24, 5, F, 4), (43, 25, F, 2), (49, 11, F, 3), (37, 8, F, 4), (61, 33, F, 2),
      (66, 15, F, 3), (60, 10, F, 4)]
snova_parameters = SL[0]

print(snova_parameters)
prob = 0.7
test_F16(snova_parameters, prob)

#Sprivate = generate_seed(32)
#Spublic = generate_seed(32)
#digest = generate_seed(32)
#salt = generate_seed(32)
#(Spublic, P22) = generate_public_key(snova_parameters, Spublic, Sprivate)
#(sig, salt) = sign_fault_F16(snova_parameters, Spublic, Sprivate, digest, salt)
#acc = signature_verification(snova_parameters, (Spublic, P22), digest, (sig, salt))
#print("Passed?", acc == True)