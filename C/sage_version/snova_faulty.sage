from hashlib import shake_256
import os
import random

load("utils.sage")


def generate_seed(number_bytes):
    return os.urandom(number_bytes)


def compute_the_vinegar_part_of_the_central_map(snova_parameters, F11, A, B, Q1, Q2, V):
    (v, o, F, l) = snova_parameters
    ring = MatrixSpace(F, l)
    # The following is the first step
    Left = []
    Right = []
    for a in range(l * l):
        L = []
        R = []
        for j in range(v):
            L.append(A[a] * V[j].transpose() * Q1[a])
            R.append(Q2[a] * V[j] * B[a])

        Left.append(L)
        Right.append(R)

    # The following is the second step
    m = o
    FVV = []
    for i in range(m):
        FVV.append(ring.zero_matrix())
        for a in range(l * l):
            for j in range(v):
                for k in range(v):
                    FVV[i] = FVV[i] + Left[a][j] * F11[i][j, k] * Right[a][k]

    # FVV is a list containing m values (from R), where each represents Fi,VV
    return FVV

def assign_values_to_vinegar_variables_fault_F16(snova_parameters, II):
    (v, o, field, l) = snova_parameters
    R = MatrixSpace(field, l)
    (x, I) = II
    V = []
    for i in range(v):
        Vi = []
        for r0 in range(l):
            for r1 in range(l):
                if x[i * l * l + r0 * l + r1] == 1:
                    Vi.append(I[i * l * l + r0 * l + r1])
                else:
                    Vi.append(field.random_element())

        V.append(R.matrix(Vi))

    return V

def compute_the_cofficient_matrix_of_the_oil_variable(snova_parameters, Fi_12, Fi_21, A, B, Q1, Q2, V, k):
    (v, o, F, l) = snova_parameters

    # first part
    Left = []
    Right = []
    for a in range(l * l):
        L = []
        R = []
        for j in range(v):
            L.append(A[a] * V[j].transpose() * Q1[a])
            R.append(Q2[a] * V[j] * B[a])

        Left.append(L)
        Right.append(R)

    R2 = MatrixSpace(F, l ** 2)
    Mik = R2.matrix()

    # second part
    for a in range(l ** 2):
        for j in range(v):
            LeftXk = Left[a][j] * Fi_12[j, k] * Q2[a]
            RightXk = B[a]
            for ti in range(l ** 2):
                for tj in range(l ** 2):
                    Mik[ti, tj] = Mik[ti, tj] + LeftXk[ti // l][tj // l] * RightXk[tj % l][ti % l]

    # At this point Mik should be M0ik
    # third part
    for a in range(l ** 2):
        for j in range(v):
            LeftXk = A[a]
            RightXk = Q1[a] * Fi_21[k, j] * Right[a][j]
            for ti in range(l ** 2):
                for tj in range(l ** 2):
                    Mik[ti, tj] = Mik[ti, tj] + LeftXk[ti // l][tj % l] * RightXk[tj // l][ti % l]

                    # At this point Mik should be M0ik+M1ik
    return Mik


def build_the_augmented_matrix_of_the_system(snova_parameters, FVV, M, digest, Spublic, salt):
    (v, o, F, l) = snova_parameters
    forder = F.order()
    m = o
    G = matrix(F, m * l * l, m * l * l + 1)

    # Generation of Yi and settting the vectorization of each Yi to the corresponding part of G's last column.
    byte_string = shake_256(Spublic + digest + salt).digest(m * l * l)
    for i in range(m * l * l):
        G[i, m * l * l] = F.from_integer(byte_string[i] % forder)

    # sum the vectorization of each FiVV to G's last column.
    for i in range(m):
        for j in range(l):
            for k in range(l):
                G[i * l * l + j * l + k, m * l * l] = G[i * l * l + j * l + k, m * l * l] + FVV[i][j, k]

    # Set each Mik at the submatrix G[i*l*l+ti,k*l*l+tj], for ti,tj in {0,1,...,l^2-1}
    for i in range(m):
        for k in range(m):
            for ti in range(l * l):
                for tj in range(l * l):
                    G[i * l * l + ti, k * l * l + tj] = M[i][k][ti, tj]

    return G


def generate_the_random_part_of_public_key(snova_parameters, Spublic):
    (v, o, field, l) = snova_parameters
    forder = field.order()
    m = o
    byte_string = shake_256(Spublic).digest(v * v * l * l * m + 2 * v * o * l * l * m + 2 * l ** 4 + 2 * (l ** 3))
    R = MatrixSpace(field, l)
    P11 = []  ## P11 is a list. It will contain m v*v matrices over R

    M_Over_R = MatrixSpace(R, v, v)

    offset = 0

    number_of_elements = m
    number_of_entries_R = l * l
    unit = v * v * number_of_entries_R
    total = unit * number_of_elements

    # obtain v*v*l*l*m bytes.
    byte_string_chunk = byte_string[offset:offset + total]

    for z in range(number_of_elements):
        byte_string_chunk_z = byte_string_chunk[
                              unit * z:unit * (z + 1)]  # a v*v matrix over R is derived from v*v*l*l bytes
        P = []
        for i in range(v):
            for j in range(v):
                C = []
                dis = (i * v + j) * number_of_entries_R

                for k in range(number_of_entries_R):
                    C.append(field.from_integer(byte_string_chunk_z[dis + k] % forder))
                # at this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix
                # Convert the list C to a matrix in R and add it to P
                P.append(R.matrix(C))
        # At this point, P should contain v*v R elements, which may be arranged to form a v*v matrix over R
        # Convert the list P to a matrix in M_Over_R and add it to P11
        P11.append(M_Over_R.matrix(P))

    # At this point, P11 should contain m v*v matrices over R

    P12 = []  ## P12 is a list. It will contain m v*o matrices over R

    offset = offset + total

    number_of_elements = m
    number_of_entries_R = l * l
    unit = v * o * number_of_entries_R
    total = unit * number_of_elements

    # obtain v*o*l*l*m bytes.
    byte_string_chunk = byte_string[offset:offset + total]
    M_Over_R = MatrixSpace(R, v, o)

    for z in range(number_of_elements):
        byte_string_chunk_z = byte_string_chunk[
                              unit * z:unit * (z + 1)]  # a v*o matrix over R is derived from v*o*l*l bytes
        P = []
        for i in range(v):
            for j in range(o):
                C = []
                dis = (i * o + j) * number_of_entries_R
                for k in range(number_of_entries_R):
                    C.append(field.from_integer(byte_string_chunk_z[dis + k] % forder))
                # At this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix
                # Convert the list C to a matrix in R and add it to P
                P.append(R.matrix(C))
        # At this point, P should contain v*v R elements, which may be arranged to form a v*v matrix over R
        # Convert the list P to a matrix in M_Over_R and add it to P12
        P12.append(M_Over_R.matrix(P))

    # At this point, P12 should contain m v*o matrices over R

    P21 = []  ## P21 is a list. It will contain m o*v matrices over R
    offset = offset + total

    number_of_elements = m
    number_of_entries_R = l * l
    unit = o * v * number_of_entries_R
    total = unit * number_of_elements

    # obtain o*v*l*l*m bytes.
    byte_string_chunk = byte_string[offset:offset + total]

    M_Over_R = MatrixSpace(R, o, v)

    for z in range(number_of_elements):
        byte_string_chunk_z = byte_string_chunk[
                              unit * z:unit * (z + 1)]  # a o*v matrix over R is derived from v*o*l*l bytes
        P = []
        for i in range(o):
            for j in range(v):
                C = []
                dis = (i * v + j) * number_of_entries_R
                for k in range(number_of_entries_R):
                    C.append(field.from_integer(byte_string_chunk_z[dis + k] % forder))
                # At this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix
                # Convert the list C to a matrix in R and add it to P
                P.append(R.matrix(C))
        # At this point, P should contain v*v R elements, which may be arranged to form a v*v matrix over R
        # Convert the list P to a matrix in M_Over_R and add it to P21
        P21.append(M_Over_R.matrix(P))

    # At this point, P21 should contain m o*v matrices over R

    A = []  # A is a list. It will contain l^2 invertible elements randomly chosen from R.

    offset = offset + total

    number_of_elements = l * l
    number_of_entries_R = l * l
    unit = number_of_entries_R
    total = unit * number_of_elements

    # obtain l^4 bytes.
    byte_string_chunk = byte_string[offset:offset + total]

    for z in range(number_of_elements):
        byte_string_chunk_z = byte_string_chunk[unit * z:unit * (z + 1)]  # a lxl matrix is derived from l*l bytes
        C = []
        for j in range(number_of_entries_R):
            C.append(field.from_integer(byte_string_chunk_z[j] % forder))
        # At this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix (R element)
        # Generate an invertible matrix in R from the derived from the list C, and add it to A.
        A.append(generate_invertible_matrix(snova_parameters, R.matrix(C)))

    # At this point, A should contain l^2 invertible elements randomly chosen from R

    B = []  # B is a list. It will contain l^2 invertible elements randomly chosen from R.

    offset = offset + total
    number_of_elements = l * l
    number_of_entries_R = l * l
    unit = number_of_entries_R
    total = unit * number_of_elements

    # obtain l^4 bytes.
    byte_string_chunk = byte_string[offset:offset + total]

    for z in range(number_of_elements):
        byte_string_chunk_z = byte_string_chunk[unit * z:unit * (z + 1)]  # a lxl matrix is derived from l*l bytes
        C = []
        for j in range(number_of_entries_R):
            C.append(F.from_integer(byte_string_chunk_z[j] % forder))

        # At this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix (R element)
        # Generate an invertible matrix in R from the derived from the list C, and add it to A.
        B.append(generate_invertible_matrix(snova_parameters, R.matrix(C)))

    # At this point, B should contain l^2 invertible elements randomly chosen from R

    Q1 = []  # Q1 is a list. It will contain l^2 invertible matrices randomly chosen from Fq[S].
    offset = offset + total
    number_of_elements = l * l
    number_of_entries_R = l * l
    unit = l
    total = unit * number_of_elements

    # obtain l^3 bytes.
    byte_string_chunk = byte_string[offset:offset + total]

    for z in range(number_of_elements):
        byte_string_chunk_z = byte_string_chunk[unit * z:unit * (z + 1)]  # a Fq[S] element is derived from l bytes
        C = []
        for i in range(l):
            C.append(F.from_integer(byte_string_chunk_z[i] % forder))
        # At this point, C should contain l GF(16) elements, which may be used to generate an element from Fq[S]
        # Generate an invertible matrix in Fq[S] from the derived from the list C, and add it to Q1.
        Q1.append(generate_invertible_matrix(snova_parameters, generate_element_in_F_S(snova_parameters, C)))

    # At this point, Q1 should contain l^2 invertible elements randomly chosen from Fq[S]

    Q2 = []  # Q2 is a list. It will contain l^2 invertible matrices randomly chosen from Fq[S].

    offset = offset + total
    number_of_elements = l * l
    number_of_entries_R = l * l
    unit = l
    total = unit * number_of_elements

    for z in range(number_of_elements):
        byte_string_chunk_z = byte_string_chunk[unit * z:unit * (z + 1)]
        C = []
        for i in range(l):
            C.append(F.from_integer(byte_string_chunk_z[i] % forder))
        # At this point, C should contain l GF(16) elements, which may be used to generate an element from Fq[S]
        # Generate an invertible matrix in Fq[S] from the derived from the list C, and add it to Q2.
        Q2.append(generate_invertible_matrix(snova_parameters, generate_element_in_F_S(snova_parameters, C)))

    # At this point, Q2 should contain l^2 invertible elements randomly chosen from Fq[S]

    return P11, P12, P21, A, B, Q1, Q2


def assign_values_to_vinegar_variables_orginal(snova_parameters, Sprivate, digest, salt, num_sig):
    (v, o, field, l) = snova_parameters
    forder = field.order()
    num = int(num_sig)
    # generate v*l*l bytes from Sprivate||digest||salt||num.to_bytes() to derive values for the vinegar variables
    byte_string = shake_256(Sprivate + digest + salt + num.to_bytes()).digest(v * l * l)

    R = MatrixSpace(field, l)
    # V is a list. It will contain v elements from R,
    # where each represents a value assigned to the corresponding vinegar variable
    V = []
    number_of_elements = v
    number_of_entries_R = l * l
    unit = number_of_entries_R

    for i in range(number_of_elements):
        # Get an string of l*l bytes to form an element in R and assign it to the vinegar variable Xi
        byte_string_chunk_i = byte_string[unit * i:unit * (i + 1)]
        C = []
        for j in range(number_of_entries_R):
            C.append(field.from_integer(byte_string_chunk_i[j] % forder))
        # At this point, C should contain l*l GF(16) elements, which may be used to generate an element in R
        # append it to V
        V.append(R.matrix(C))

    # At this point, V should contain v elements in R,
    # where each represents a value assigned to the corresponding variable.
    return V


def generate_public_key(snova_parameters, seed_public, Sprivate):
    (v, o, field, l) = snova_parameters
    m = o
    # Generate the v*o matrix T12 over R from Sprivate
    T12 = generate_linear_T12(snova_parameters, Sprivate)
    # generate the random part of public key from seed_public
    P11, P12, P21, A, B, Q1, Q2 = generate_the_random_part_of_public_key(snova_parameters, seed_public)
    P22 = []
    # Compute P22[i] = T12.t*(P11[i]*T12 + P12[i]) + P21[i]*T12 for all i in [m]
    for i in range(m):
        I1 = P11[i]._multiply_strassen(T12, 1) + P12[i]
        I2 = P21[i]._multiply_strassen(T12, 1)
        T12t = T12.transpose()
        P22.append(T12t._multiply_strassen(I1, 1) + I2)

    # At this point, the list P22 should contain m o*o matrices over R.
    return (seed_public, P22)


def generate_private_key(snova_parameters, Spublic, Sprivate):
    T12 = generate_linear_T12(snova_parameters, Sprivate)  # Generate the v*o matrix T12 over R from Sprivate
    (v, o, field, l) = snova_parameters
    m = o
    # generate the random part of public key from Spublic
    P11, P12, P21, A, B, Q1, Q2 = generate_the_random_part_of_public_key(snova_parameters, Spublic)
    F11 = []
    F12 = []
    F21 = []
    # The following for loop computes:
    # F11[i]=P11[i]
    # F12[i] =P11[i]*T12 +P12[i]
    # F21[i] =(T12).t*P11[i] +P21[i]
    for i in range(m):
        F11.append(P11[i])

        F12.append(P11[i]._multiply_strassen(T12, 1) + P12[i])
        F21.append(T12.transpose()._multiply_strassen(P11[i], 1) + P21[i])

    # At this point,
    # The list F11 should contain m v*v matrices over R.
    # The list F12 should contain m v*o matrices over R.
    # The list F21 should contain m o*v matrices over R.

    return T12, F11, F12, F21


def sign_fault_F16(snova_parameters, Spublic, Sprivate, digest, salt, I):
    (v, o, F, l) = snova_parameters
    m = o
    # Generate the random part of public key from Spublic
    _, _, _, A, B, Q1, Q2 = generate_the_random_part_of_public_key(snova_parameters, Spublic)
    # Generate the corresponding private_key from both Spublic and Sprivate
    T12, F11, F12, F21 = generate_private_key(snova_parameters, Spublic, Sprivate)
    # Build the augmented matrix [T]
    T = block_matrix_T(snova_parameters, T12)
    R = MatrixSpace(F, l)
    # Set numsig as 0
    numsig = 0
    # Set is_done as False
    is_done = False
    while (not is_done):
        ## Asume all vallues for i in I are fixed. I is a subset of {0,1,2,..., v-1}
        V = assign_values_to_vinegar_variables_fault_F16(snova_parameters, I)
        # V=assign_values_to_vinegar_variables(snova_parameters,Sprivate,digest,salt,numsig)

        # Evaluate the vinegar part for each Fi.
        FVV = compute_the_vinegar_part_of_the_central_map(snova_parameters, F11, A, B, Q1, Q2, V)
        # The following for loop computes Mik for all i,k in {0,1,...,m-1}.
        # It computes M, a m*m matrix where M[i][k] stores Mik
        M = []
        for i in range(m):
            Mi = []
            for k in range(o):
                Mik = compute_the_cofficient_matrix_of_the_oil_variable(snova_parameters, F12[i], F21[i], A, B, Q1, Q2,V, k)
                Mi.append(Mik)

            M.append(Mi)
        # build the augmented matrix of the system G, an (ml^2)*(ml^2+1) matrix over F16,from FVV,M,digest,Spublic,salt
        G = build_the_augmented_matrix_of_the_system(snova_parameters, FVV, M, digest, Spublic, salt)
        out, O = gauss(snova_parameters, G)
        # if out = true, the equation system represented by G is solved
        if out:
            is_done = True
            # update the list V to {V0,V1,...V{v-1},O0,O1,...,O{o-1}}
            V = V + O
            # obtain a column vector MV from V
            MR = MatrixSpace(R, 1, len(V))
            MV = MR.matrix(V)
            # compute signature = [T](MV).t
            signature = T._multiply_strassen(MV.transpose(), 1)
            # The following for loop converts signature to a list containing n elements in R.
            sig = []
            for i in range(v + o):
                sig.append(signature[i, 0])
        else:
            # if out = false, the equation system represented by G is not solved
            # numsig is increased by one to generate new random values V_0,V_1, ..., V_{v-1} in R in the next iteration
            is_done = False
            numsig = numsig + 1

    return (sig, salt)


def sign_fault_F2(snova_parameters, Spublic, Sprivate, digest, salt, I):
    (v, o, F, l) = snova_parameters
    m = o
    # Generate the random part of public key from Spublic
    _, _, _, A, B, Q1, Q2 = generate_the_random_part_of_public_key(snova_parameters, Spublic)
    # Generate the corresponding private_key from both Spublic and Sprivate
    T12, F11, F12, F21 = generate_private_key(snova_parameters, Spublic, Sprivate)
    # Build the augmented matrix [T]
    T = block_matrix_T(snova_parameters, T12)
    R = MatrixSpace(F, l)
    # Set numsig as 0
    numsig = 0
    # Set is_done as False
    is_done = False
    while (not is_done):
        ## Asume all vallues for i in I are fixed. I is a subset of {0,1,2,..., v-1}
        V = assign_values_to_vinegar_variables_fault_F2(snova_parameters, I)
        # V=assign_values_to_vinegar_variables(snova_parameters,Sprivate,digest,salt,numsig)

        # Evaluate the vinegar part for each Fi.
        FVV = compute_the_vinegar_part_of_the_central_map(snova_parameters, F11, A, B, Q1, Q2, V)
        # The following for loop computes Mik for all i,k in {0,1,...,m-1}.
        # It computes M, a m*m matrix where M[i][k] stores Mik
        M = []
        for i in range(m):
            Mi = []
            for k in range(o):
                Mik = compute_the_coeﬀicient_matrix_of_the_oil_variable(snova_parameters, F12[i], F21[i], A, B, Q1, Q2,
                                                                        V, k)
                Mi.append(Mik)

            M.append(Mi)
        # build the augmented matrix of the system G, an (ml^2)*(ml^2+1) matrix over F16,from FVV,M,digest,Spublic,salt
        G = build_the_augmented_matrix_of_the_system(snova_parameters, FVV, M, digest, Spublic, salt)
        out, O = gauss(snova_parameters, G)
        # if out = true, the equation system represented by G is solved
        if out:
            is_done = True
            # update the list V to {V0,V1,...V{v-1},O0,O1,...,O{o-1}}
            V = V + O
            # obtain a column vector MV from V
            MR = MatrixSpace(R, 1, len(V))
            MV = MR.matrix(V)
            # compute signature = [T](MV).t
            signature = T._multiply_strassen(MV.transpose(), 1)
            # The following for loop converts signature to a list containing n elements in R.
            sig = []
            for i in range(v + o):
                sig.append(signature[i, 0])
        else:
            # if out = false, the equation system represented by G is not solved
            # numsig is increased by one to generate new random values V_0,V_1, ..., V_{v-1} in R in the next iteration
            is_done = False
            numsig = numsig + 1

    return (sig, salt)


def signature_verification(snova_parameters, public_key, digest, signature):
    (Spublic, P22) = public_key
    (sig, salt) = signature
    (v, o, F, l) = snova_parameters
    forder = F.order()
    m = o
    ring = MatrixSpace(F, l)
    # Generate the random part of public key from the public seed
    P11, P12, P21, A, B, Q1, Q2 = generate_the_random_part_of_public_key(snova_parameters, Spublic)
    # Generation of a byte string from Spublic||digest|salt to derive the expected hash
    byte_string = shake_256(Spublic + digest + salt).digest(m * l * l)
    hash1 = []
    for i in range(m):
        hh = []
        for j in range(l * l):
            hh.append(F.from_integer(byte_string[i * l * l + j] % forder))
        hash1.append(ring.matrix(hh))
    # Evaluate the public map P on input sig, i.e. hash2 = P[P_0[sig],P_1[sig],...,P_{m−1}[sig]]
    hash2 = evaluate_the_public_map(snova_parameters, A, B, Q1, Q2, P11, P12, P21, P22, sig)

    ## Comparison between hash1 and hash2.
    acc = True
    for i in range(m):
        acc = acc and (hash1[i] == hash2[i])

    # acc will end up as True if and only if hash1[i]==hash2[i] for all i

    return acc
