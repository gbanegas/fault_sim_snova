from hashlib import shake_256
import os
import random


def sample_random_S(Fq, l, max_iter=10):
    iter_ = 1
    while iter_ < max_iter + 1:
        S = random_matrix(Fq, l, l)
        char_poly_S = S.characteristic_polynomial()
        if char_poly_S.is_irreducible():
            #print(“Matrix S found”)
            return S
        iter_ += 1
        #print(iter_)
    
    return None

def generate_element_in_F16_S(snova_parameters,scalars):
    (v,o,F,l)=snova_parameters
    
    if  len(scalars)!= l:
        return 
    
    zero=F.from_integer(0)
    if scalars[l-1]==zero:
        if scalars[0]!=zero:
            scalars[l-1]=F.from_integer(16-scalars[0].to_integer())
        else:
            scalars[l-1]=F.from_integer(15)
    
    R = MatrixSpace(F, l)    
    S=dict_for_S[l] 
    M = R.identity_matrix()
    O=R.matrix()
    
    for i in range(l):
        O=O+scalars[i]*M
        M=M*S
      
    
    return O


def generate_element_in_Fq_S(snova_parameters,scalars):
    (v,o,F,l)=snova_parameters
    forder=F.order()
    if  len(scalars)!= l:
        return 
    
    zero=F.from_integer(0)
    if scalars[l-1]==zero:
        scalars[l-1]=F.from_integer(forder-1)
    
    R = MatrixSpace(F, l)    
    S=dict_for_S[l] 
    M = R.identity_matrix()
    O=R.matrix()
    
    for i in range(l):
        O=O+scalars[i]*M
        M=M*S
      
    
    return O


def generate_element_in_F_S(snova_parameters,scalars):
    (v,o,F,l)=snova_parameters
    forder=F.order()
    if forder==16:
        return generate_element_in_F16_S(snova_parameters,scalars)
    else:
        return generate_element_in_Fq_S(snova_parameters,scalars)

# This algorithm generates an lxl invertible matrix
def generate_invertible_matriz (snova_parameters,M):
    (v,o,F,l)=snova_parameters
     
    S=dict_for_S[l]
    zero=F.from_integer(0)
    if M.determinant()==zero:
        for i in range(1,F.order()):
            if (M+F.from_integer(i)*S).determinant()!=zero:
                M=M+F.from_integer(i)*S
                break
    
    return M

def generate_seed(number_bytes):
 
    return os.urandom(number_bytes)

#Generates a v*o matrix T12 from Sprivate.
#T12 is a v × o matrix consisting of nonzero entries T12_{i,j} chosen randomly in Fq [S ]
def generate_linear_T12 (snova_parameters, Sprivate):
     
    (v,o,field,l)=snova_parameters
    forder=field.order()
    byte_string=shake_256(Sprivate).digest(v*o*l)#Each F_16[S] element may be generated from l scalars 
    T12=[]
    R = MatrixSpace(field, l)
    M_Over_R= MatrixSpace(R,v,o)
    for i in range(v):
        for j in range(o):
            scalars=[]
            for k in range(l):
                 scalars.append(field.from_integer(byte_string[(i*o +j)*l+k]%forder)) # Get a byte and convert it to a field element
            
            #At this point, scalars should contain l scalars in F_16
            entry=generate_element_in_F_S(snova_parameters,scalars) # generates the nonzero entry T12_{i,j} from scalars
            T12.append(entry)
        
    return M_Over_R.matrix(T12)

def generate_the_random_part_of_public_key(snova_parameters,Spublic):
        
        
        (v,o,field,l)=snova_parameters
        forder=field.order()
        m=o
        byte_string=shake_256(Spublic).digest(v*v*l*l*m + 2*v*o*l*l*m +2*l**4 + 2*(l**3)) 
        
       
        R = MatrixSpace(field, l)
        
        
        
               
        P11=[]## P11 is a list. It will contain m v*v matrices over R
        
        M_Over_R= MatrixSpace(R,v,v)
        
        offset=0
        
        number_of_elements=m
        number_of_entries_R=l*l
        unit=v*v*number_of_entries_R
        total=unit*number_of_elements
        
        # obtain v*v*l*l*m bytes. 
        byte_string_chunk=byte_string[offset:offset+total] 
        
        
        for z in range(number_of_elements):
            byte_string_chunk_z=byte_string_chunk[unit*z:unit*(z+1)] # a v*v matrix over R is derived from v*v*l*l bytes
            P=[]
            for i in range(v):
                for j in range(v):
                    C=[]
                    dis=(i*v+j)*number_of_entries_R
                    
                    for k in range(number_of_entries_R):
                        C.append(field.from_integer(byte_string_chunk_z[dis+k]%forder))
                    #at this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix
                    #Convert the list C to a matrix in R and add it to P
                    P.append(R.matrix(C))
            #At this point, P should contain v*v R elements, which may be arranged to form a v*v matrix over R
            #Convert the list P to a matrix in M_Over_R and add it to P11
            P11.append(M_Over_R.matrix(P))
        
        #At this point, P11 should contain m v*v matrices over R
    
    
        P12=[] ## P12 is a list. It will contain m v*o matrices over R
        
        offset=offset+total 
        
        number_of_elements=m
        number_of_entries_R=l*l
        unit=v*o*number_of_entries_R
        total=unit*number_of_elements
        
        # obtain v*o*l*l*m bytes. 
        byte_string_chunk=byte_string[offset:offset+total]
        M_Over_R= MatrixSpace(R,v,o)
            
        for z in range(number_of_elements):
            byte_string_chunk_z=byte_string_chunk[unit*z:unit*(z+1)]# a v*o matrix over R is derived from v*o*l*l bytes
            P=[]
            for i in range(v):
                for j in range(o):
                    C=[]
                    dis=(i*o+j)*number_of_entries_R
                    for k in range(number_of_entries_R):
                        C.append(field.from_integer(byte_string_chunk_z[dis+k]%forder))
                    #At this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix
                    #Convert the list C to a matrix in R and add it to P
                    P.append(R.matrix(C))
            #At this point, P should contain v*v R elements, which may be arranged to form a v*v matrix over R
            #Convert the list P to a matrix in M_Over_R and add it to P12
            P12.append(M_Over_R.matrix(P))
            
        #At this point, P12 should contain m v*o matrices over R
        
           
        P21=[] ## P21 is a list. It will contain m o*v matrices over R 
        offset=offset+total
        
        number_of_elements=m
        number_of_entries_R=l*l
        unit=o*v*number_of_entries_R
        total=unit*number_of_elements
        
        # obtain o*v*l*l*m bytes.
        byte_string_chunk=byte_string[offset:offset+ total]
        
        M_Over_R= MatrixSpace(R,o,v)
        
        for z in range(number_of_elements):
            byte_string_chunk_z=byte_string_chunk[unit*z:unit*(z+1)] # a o*v matrix over R is derived from v*o*l*l bytes
            P=[]
            for i in range(o):
                for j in range(v):
                    C=[]
                    dis=(i*v+j)*number_of_entries_R
                    for k in range(number_of_entries_R):
                        C.append(field.from_integer(byte_string_chunk_z[dis+k]%forder))
                    #At this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix
                    #Convert the list C to a matrix in R and add it to P
                    P.append(R.matrix(C))
            #At this point, P should contain v*v R elements, which may be arranged to form a v*v matrix over R
            #Convert the list P to a matrix in M_Over_R and add it to P21
            P21.append(M_Over_R.matrix(P))
            
        #At this point, P21 should contain m o*v matrices over R
        
        A=[]# A is a list. It will contain l^2 invertible elements randomly chosen from R.
        
        offset=offset+total
        
        number_of_elements=l*l
        number_of_entries_R=l*l
        unit=number_of_entries_R
        total=unit*number_of_elements
        
        # obtain l^4 bytes.
        byte_string_chunk=byte_string[offset:offset+total]
        
        for z in range(number_of_elements):
            byte_string_chunk_z=byte_string_chunk[unit*z:unit*(z+1)] # a lxl matrix is derived from l*l bytes
            C=[]
            for j in range(number_of_entries_R):
                 C.append(field.from_integer(byte_string_chunk_z[j]%forder))
            #At this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix (R element)
            #Generate an invertible matrix in R from the derived from the list C, and add it to A.   
            A.append(generate_invertible_matriz (snova_parameters,R.matrix(C)))
            
        #At this point, A should contain l^2 invertible elements randomly chosen from R
        
        B=[] # B is a list. It will contain l^2 invertible elements randomly chosen from R.
        
        offset=offset+ total
        number_of_elements=l*l
        number_of_entries_R=l*l
        unit=number_of_entries_R
        total=unit*number_of_elements
        
        # obtain l^4 bytes.
        byte_string_chunk=byte_string[offset:offset+total]
        
        for z in range(number_of_elements):
            byte_string_chunk_z=byte_string_chunk[unit*z:unit*(z+1)] # a lxl matrix is derived from l*l bytes
            C=[]
            for j in range(number_of_entries_R):
                 C.append(F.from_integer(byte_string_chunk_z[j]%forder))
            
            #At this point, C should contain l*l GF(16) elements, which may be arranged to form a l*l matrix (R element)
            #Generate an invertible matrix in R from the derived from the list C, and add it to A.  
            B.append(generate_invertible_matriz(snova_parameters, R.matrix(C)))
            
        #At this point, B should contain l^2 invertible elements randomly chosen from R
        
        Q1=[] # Q1 is a list. It will contain l^2 invertible matrices randomly chosen from Fq[S].
        offset=offset+ total
        number_of_elements=l*l
        number_of_entries_R=l*l
        unit=l
        total=unit*number_of_elements
        
        # obtain l^3 bytes.
        byte_string_chunk=byte_string[offset:offset+total]
        
        for z in range(number_of_elements):
            byte_string_chunk_z=byte_string_chunk[unit*z:unit*(z+1)]# a Fq[S] element is derived from l bytes
            C=[]
            for i in range(l):
                C.append(F.from_integer(byte_string_chunk_z[i]%forder))
            #At this point, C should contain l GF(16) elements, which may be used to generate an element from Fq[S]
            #Generate an invertible matrix in Fq[S] from the derived from the list C, and add it to Q1. 
            Q1.append(generate_invertible_matriz(snova_parameters,generate_element_in_F_S(snova_parameters,C)) )
        
        #At this point, Q1 should contain l^2 invertible elements randomly chosen from Fq[S]
        
        Q2=[]# Q2 is a list. It will contain l^2 invertible matrices randomly chosen from Fq[S].
        
        offset=offset+ total
        number_of_elements=l*l
        number_of_entries_R=l*l
        unit=l
        total=unit*number_of_elements
        
        for z in range(number_of_elements):
            byte_string_chunk_z=byte_string_chunk[unit*z:unit*(z+1)]
            C=[]
            for i in range(l):
                C.append(F.from_integer(byte_string_chunk_z[i]%forder))
            #At this point, C should contain l GF(16) elements, which may be used to generate an element from Fq[S]
            #Generate an invertible matrix in Fq[S] from the derived from the list C, and add it to Q2.     
            Q2.append(generate_invertible_matriz(snova_parameters,generate_element_in_F_S(snova_parameters,C)))
         
        #At this point, Q2 should contain l^2 invertible elements randomly chosen from Fq[S]
        
        
        return  P11,P12,P21,A,B,Q1,Q2

def generate_public_key (snova_parameters,Spublic, Sprivate):
        (v,o,field,l)=snova_parameters
        m=o
        #Generate the v*o matrix T12 over R from Sprivate
        T12=generate_linear_T12 (snova_parameters, Sprivate) 
        #generate the random part of public key from Spublic
        P11,P12,P21,A,B,Q1,Q2 =generate_the_random_part_of_public_key(snova_parameters,Spublic)
        P22=[]
        #Compute P22[i] = T12.t*(P11[i]*T12 + P12[i]) + P21[i]*T12 for all i in [m]
        for i in range(m):
            I1=P11[i]._multiply_strassen(T12,1)+P12[i]
            I2=P21[i]._multiply_strassen(T12,1)
            T12t=T12.transpose()
            P22.append(T12t._multiply_strassen(I1,1)+I2)
        
        #At this point, the list P22 should contain m o*o matrices over R.
        return (Spublic,P22)

def generate_private_key(snova_parameters,Spublic, Sprivate):
        
        T12=generate_linear_T12 (snova_parameters, Sprivate)#Generate the v*o matrix T12 over R from Sprivate
        (v,o,field,l)=snova_parameters
        m=o
        #generate the random part of public key from Spublic
        P11,P12,P21,A,B,Q1,Q2 =generate_the_random_part_of_public_key(snova_parameters,Spublic)
        F11=[]
        F12=[]
        F21=[]
        # The following for loop computes:
        # F11[i]=P11[i]
        # F12[i] =P11[i]*T12 +P12[i]
        # F21[i] =(T12).t*P11[i] +P21[i]
        for i in range(m):
            
            F11.append(P11[i])
            
            F12.append(P11[i]._multiply_strassen(T12,1)+P12[i])
            F21.append(T12.transpose()._multiply_strassen(P11[i],1)+P21[i])
        
        #At this point, 
        #The list F11 should contain m v*v matrices over R.
        #The list F12 should contain m v*o matrices over R.
        #The list F21 should contain m o*v matrices over R.

        return T12,F11,F12,F21

import random

def assign_values_to_vinegar_variables_fault_R(snova_parameters,I):
        (v,o,field,l)=snova_parameters
        R = MatrixSpace(field, l)
        V=[]
        
        for i in range(v):
            if i in I.keys():
                value=I[i]
            else:
                value=R.random_element()
            
            V.append(value)
            
        
        return V

def assign_values_to_vinegar_variables_orginal(snova_parameters,Sprivate,digest,salt,num_sig):
        (v,o,field,l)=snova_parameters
        forder=field.order()
        num=int(num_sig)
        #generate v*l*l bytes from Sprivate||digest||salt||num.to_bytes() to derive values for the vinegar variables
        byte_string=shake_256(Sprivate+digest+salt+num.to_bytes()).digest(v*l*l)
        
        R = MatrixSpace(field, l)
        # V is a list. It will contain v elements from R,
        #where each represents a value assigned to the corresponding vinegar variable
        V=[]
        number_of_elements=v
        number_of_entries_R=l*l
        unit=number_of_entries_R
        
        for i in range(number_of_elements):
            #Get an string of l*l bytes to form an element in R and assign it to the vinegar variable Xi
            byte_string_chunk_i=byte_string[unit*i:unit*(i+1)]
            C=[]
            for j in range(number_of_entries_R):
                 C.append(field.from_integer(byte_string_chunk_i[j]%forder))
            #At this point, C should contain l*l GF(16) elements, which may be used to generate an element in R
            #append it to V 
            V.append(R.matrix(C))
            
        # At this point, V should contain v elements in R, 
        #where each represents a value assigned to the corresponding variable.
        return V
 

  
def assign_values_to_vinegar_variables_fault_F16(snova_parameters,II):
        (v,o,field,l)=snova_parameters
        R = MatrixSpace(field, l)
        (x,I)=II
        V=[]
        for i in range(v):
            Vi=[]
            for r0 in range(l):
                for r1 in range(l):
                    if x[i*l*l+r0*l+r1]==1:
                        Vi.append(I[i*l*l+r0*l+r1])
                    else:
                        Vi.append(field.random_element())
            
            V.append(R.matrix(Vi))
        
        return V
    
    
def assign_values_to_vinegar_variables_fault_F2(snova_parameters,II):
        (v,o,field,l)=snova_parameters
        number_bits=int(log(field.order(), 2))
        R = MatrixSpace(field, l)
        (x,I)=II
        V=[]
        F2=field.base_ring()
        for i in range(v):
            Vi=[]
            for r0 in range(l):
                for r1 in range(l):
                    rand_element=[]
                    for r2 in range(number_bits):
                        if x[i*l*l+r0*l+r1*number_bits+r2]==1:
                            rand_element.append(I[i*l*l+r0*l+r1*number_bits+r2])
                        else:
                            rand_element.append(F2.random_element())
                    Vi.append(field(rand_element))
                    
            V.append(R.matrix(Vi))
        
        return V



def get_random_combination(u2,v):
    
    L=[i for i in range(v)]
    O=[]
    
    for _ in range(u2):
        random.shuffle(L)
        O.append(L[0])
        L=L[1:]
        
    return O
# get a dictionary such that i in I is a random value.   
def get_R(snova_parameters, u2):
    (v,o,field,l)=snova_parameters
    
    R = MatrixSpace(field, l)
    I={}
    
    keys=get_random_combination(u2,v)
    
    for i in keys:
        I[i]=R.random_element()
        
    return I

from numpy import random

def get_F16(snova_parameters,prob):
    (v,o,field,l)=snova_parameters
    
    
    I=[]
    for i in range(v*l*l):
         I.append(field.random_element())
            
    x = random.binomial(n=1, p=prob, size=v*l*l)
                
                
        
    return (x,I)
    


def get_F2(snova_parameters,prob):
    (v,o,field,l)=snova_parameters
    
    
    I=[]
    number_bits=int(log(field.order(), 2))
    F2=field.base_ring()
    for i in range(number_bits*v*l*l):
        I.append(F2.random_element())
            
    x = random.binomial(n=1, p=prob, size=number_bits*v*l*l)
                
                
        
    return (x,I)

def compute_the_vinegar_part_of_the_central_map (snova_parameters,F11,A,B,Q1,Q2,V):
    (v,o,F,l)=snova_parameters
    ring = MatrixSpace(F, l)
    # The following is the first step
    Left=[]
    Right=[]
    for a in range(l*l):
        L=[]
        R=[]
        for j in range(v):
            L.append(A[a]*V[j].transpose()*Q1[a])
            R.append(Q2[a]*V[j]*B[a])
            
        Left.append(L)
        Right.append(R)

    # The following is the second step   
    m=o
    FVV=[]
    for i in range(m):
        FVV.append(ring.zero_matrix())
        for a in range(l*l):
            for j in range(v):
                for k in range(v):
                    FVV[i]=FVV[i]+Left[a][j]*F11[i][j,k]*Right[a][k]
    
    #FVV is a list containing m values (from R), where each represents Fi,VV
    return FVV

def compute_the_coeﬀicient_matrix_of_the_oil_variable(snova_parameters,Fi_12,Fi_21,A,B,Q1,Q2,V,k):
    
    (v,o,F,l)=snova_parameters
   
    #first part
    Left=[]
    Right=[]
    for a in range(l*l):
        L=[]
        R=[]
        for j in range(v):
            L.append(A[a]*V[j].transpose()*Q1[a])
            R.append(Q2[a]*V[j]*B[a])
            
        Left.append(L)
        Right.append(R)
        
    R2=MatrixSpace(F, l**2)
    Mik=R2.matrix()
    
    #second part
    for a in range(l**2):
        for j in range(v):
            LeftXk=Left[a][j]*Fi_12[j,k]*Q2[a]
            RightXk=B[a]
            for ti in range(l**2):
                for tj in range(l**2): 
                    Mik[ti,tj]=Mik[ti,tj]+LeftXk[ti//l][tj//l]*RightXk[tj%l][ti%l]
     
    #At this point Mik should be M0ik  
    #third part               
    for a in range(l**2):
        for j in range(v):
            LeftXk=A[a]
            RightXk=Q1[a]*Fi_21[k,j]*Right[a][j]
            for ti in range(l**2):
                for tj in range(l**2): 
                    Mik[ti,tj]=Mik[ti,tj]+LeftXk[ti//l][tj%l]*RightXk[tj//l][ti%l]         
    
    #At this point Mik should be M0ik+M1ik 
    return Mik


def build_the_augmented_matrix_of_the_system(snova_parameters,FVV,M,digest,Spublic,salt):

    (v,o,F,l)=snova_parameters
    forder=F.order()
    m=o
    G=matrix(F,m*l*l,m*l*l+1)
    
    #Generation of Yi and settting the vectorization of each Yi to the corresponding part of G's last column.  
    byte_string=shake_256(Spublic+digest+salt).digest(m*l*l)
    for i in range(m*l*l):
        G[i,m*l*l]= F.from_integer(byte_string[i]%forder)
    
    #sum the vectorization of each FiVV to G's last column.
    for i in range(m):
        for j in range(l):
            for k in range(l):
                 G[i*l*l+j*l+k,m*l*l]=G[i*l*l+j*l+k,m*l*l]+FVV[i][j,k]
                    
                    
    #Set each Mik at the submatrix G[i*l*l+ti,k*l*l+tj], for ti,tj in {0,1,...,l^2-1}           
    for i in range(m):
        for k in range(m):
            for ti in range(l*l):
                for tj in range(l*l):
                    G[i*l*l+ti,k*l*l+tj]=M[i][k][ti,tj]
   

    return G

def block_matrix_T(snova_parameters, T12):
    
    (v,o,F,l)=snova_parameters
    n=v+o
    R=MatrixSpace(F,l)
    zero=R.zero_matrix()
    one=R.one()
    NN=MatrixSpace(R,n)
    
    
    T=[]
    # First step
    for i in range(v):
        Ti=[]
        for j in range(n):
            if(j<v):
                if i==j:
                    Ti.append(one)
                else:
                    Ti.append(zero)
                
            else:
                Ti.append(T12[i,j-v])
        
        T.append(Ti) 
    # Second step
    for i in range(o):
        Ti=[]
        for j in range(n):
            if(j<v):
                Ti.append(zero)
            else:
                if i==(j-v):
                    Ti.append(one)
                else:
                    Ti.append(zero)
        
        
        T.append(Ti) 

    
    
    return NN.matrix(T)

def gauss(snova_parameters,G):
   
    (v,o,F,l)=snova_parameters
    R=MatrixSpace(F,l)
    nrows=G.nrows()
    ncols=G.ncols()
    O=[]
    try:
        A=G[0:nrows,0:nrows]
        Y=G[0:nrows,nrows:ncols]
        #Attempt to solve the equation system AX=Y
        X= A\Y 
        num_entries=l*l
        #Since X is a (ml^2 x 1) matrix over F16, it is converted to a list contaning m elements in R.
        for k in range(o):
            MM=[]
            dis=k*num_entries
            for i in range(num_entries):
                MM.append(X[dis+i,0])
            
            O.append(R.matrix(MM))
        #If this line is reached, the equation system represented by G is solved.     
        out= True
    except Exception as ee:
        #If this line is reached, the equation system represented by G is not solved. 
        O=[]
        out= False
        
    return out,O


def sign_original(snova_parameters,Spublic, Sprivate,digest,salt):
    
    (v,o,F,l) =snova_parameters
    m=o
    #Generate the random part of public key from Spublic
    _,_,_,A,B,Q1,Q2 = generate_the_random_part_of_public_key(snova_parameters,Spublic)
    #Generate the corresponding private_key from both Spublic and Sprivate
    T12,F11,F12,F21 = generate_private_key(snova_parameters,Spublic, Sprivate)
    # Build the augmented matrix [T]
    T=block_matrix_T(snova_parameters, T12)
    R=MatrixSpace(F,l)
    #Set numsig as 0
    numsig=0
    #Set is_done as False
    is_done=False
    while (not is_done):
        ## Asume all vallues for i in I are fixed. I is a subset of {0,1,2,..., v-1}
        V=assign_values_to_vinegar_variables_orginal(snova_parameters,Sprivate,digest,salt,numsig)
        
                
            
        #V=assign_values_to_vinegar_variables(snova_parameters,Sprivate,digest,salt,numsig)
        
        #Evaluate the vinegar part for each Fi. 
        FVV= compute_the_vinegar_part_of_the_central_map (snova_parameters,F11,A,B,Q1,Q2,V)
        #The following for loop computes Mik for all i,k in {0,1,...,m-1}. 
        #It computes M, a m*m matrix where M[i][k] stores Mik
        M=[]
        for i in range(m):
            Mi=[]
            for k in range(o):
                Mik=compute_the_coeﬀicient_matrix_of_the_oil_variable(snova_parameters,F12[i],F21[i],A,B,Q1,Q2,V,k)
                Mi.append(Mik)
                
            M.append(Mi)
        #build the augmented matrix of the system G, an (ml^2)*(ml^2+1) matrix over F16,from FVV,M,digest,Spublic,salt
        G=build_the_augmented_matrix_of_the_system(snova_parameters,FVV,M,digest,Spublic,salt)    
        out,O=gauss(snova_parameters,G)
        #if out = true, the equation system represented by G is solved
        if out:
            is_done=True
            # update the list V to {V0,V1,...V{v-1},O0,O1,...,O{o-1}}
            V=V+O
            #obtain a column vector MV from V
            MR=MatrixSpace(R,1,len(V))
            MV=MR.matrix(V)
            #compute signature = [T](MV).t
            signature=T._multiply_strassen(MV.transpose(),1)
            #The following for loop converts signature to a list containing n elements in R.
            sig=[]
            for i in range(v+o):
                sig.append(signature[i,0])
        else:    
            #if out = false, the equation system represented by G is not solved 
            # numsig is increased by one to generate new random values V_0,V_1, ..., V_{v-1} in R in the next iteration
            is_done=False
            numsig=numsig+1 
            
         
    return (sig,salt)


def sign_leaked_vinegar_values(snova_parameters,Spublic, Sprivate,digest,salt, index):
    
    (v,o,F,l) =snova_parameters
    m=o
    #Generate the random part of public key from Spublic
    _,_,_,A,B,Q1,Q2 = generate_the_random_part_of_public_key(snova_parameters,Spublic)
    #Generate the corresponding private_key from both Spublic and Sprivate
    T12,F11,F12,F21 = generate_private_key(snova_parameters,Spublic, Sprivate)
    # Build the augmented matrix [T]
    T=block_matrix_T(snova_parameters, T12)
    R=MatrixSpace(F,l)
    #Set numsig as 0
    numsig=0
    #Set is_done as False
    is_done=False
    leaked_values=[]
    while (not is_done):
        ## Asume all vallues for i in I are fixed. I is a subset of {0,1,2,..., v-1}
        V=assign_values_to_vinegar_variables_orginal(snova_parameters,Sprivate,digest,salt,numsig)
        
                
            
        #V=assign_values_to_vinegar_variables(snova_parameters,Sprivate,digest,salt,numsig)
        
        #Evaluate the vinegar part for each Fi. 
        FVV= compute_the_vinegar_part_of_the_central_map (snova_parameters,F11,A,B,Q1,Q2,V)
        #The following for loop computes Mik for all i,k in {0,1,...,m-1}. 
        #It computes M, a m*m matrix where M[i][k] stores Mik
        M=[]
        for i in range(m):
            Mi=[]
            for k in range(o):
                Mik=compute_the_coeﬀicient_matrix_of_the_oil_variable(snova_parameters,F12[i],F21[i],A,B,Q1,Q2,V,k)
                Mi.append(Mik)
                
            M.append(Mi)
        #build the augmented matrix of the system G, an (ml^2)*(ml^2+1) matrix over F16,from FVV,M,digest,Spublic,salt
        G=build_the_augmented_matrix_of_the_system(snova_parameters,FVV,M,digest,Spublic,salt)    
        out,O=gauss(snova_parameters,G)
        #if out = true, the equation system represented by G is solved
        if out:
            for fi in range(v):
                for fj in range(l):
                    leaked_values.append(V[fi][fj, index])
                    
            is_done=True
            # update the list V to {V0,V1,...V{v-1},O0,O1,...,O{o-1}}
            V=V+O
            #obtain a column vector MV from V
            MR=MatrixSpace(R,1,len(V))
            MV=MR.matrix(V)
            #compute signature = [T](MV).t
            signature=T._multiply_strassen(MV.transpose(),1)
            #The following for loop converts signature to a list containing n elements in R.
            sig=[]
            for i in range(v+o):
                sig.append(signature[i,0])
        else:    
            #if out = false, the equation system represented by G is not solved 
            # numsig is increased by one to generate new random values V_0,V_1, ..., V_{v-1} in R in the next iteration
            is_done=False
            numsig=numsig+1 
            
         
    return (sig,salt, leaked_values)



def sign_fault_F16(snova_parameters,Spublic, Sprivate,digest,salt, I):
    
    (v,o,F,l) =snova_parameters
    m=o
    #Generate the random part of public key from Spublic
    _,_,_,A,B,Q1,Q2 = generate_the_random_part_of_public_key(snova_parameters,Spublic)
    #Generate the corresponding private_key from both Spublic and Sprivate
    T12,F11,F12,F21 = generate_private_key(snova_parameters,Spublic, Sprivate)
    # Build the augmented matrix [T]
    T=block_matrix_T(snova_parameters, T12)
    R=MatrixSpace(F,l)
    #Set numsig as 0
    numsig=0
    #Set is_done as False
    is_done=False
    while (not is_done):
        ## Asume all vallues for i in I are fixed. I is a subset of {0,1,2,..., v-1}
        V=assign_values_to_vinegar_variables_fault_F16(snova_parameters,I)
        #V=assign_values_to_vinegar_variables(snova_parameters,Sprivate,digest,salt,numsig)
        
        #Evaluate the vinegar part for each Fi. 
        FVV= compute_the_vinegar_part_of_the_central_map (snova_parameters,F11,A,B,Q1,Q2,V)
        #The following for loop computes Mik for all i,k in {0,1,...,m-1}. 
        #It computes M, a m*m matrix where M[i][k] stores Mik
        M=[]
        for i in range(m):
            Mi=[]
            for k in range(o):
                Mik=compute_the_coeﬀicient_matrix_of_the_oil_variable(snova_parameters,F12[i],F21[i],A,B,Q1,Q2,V,k)
                Mi.append(Mik)
                
            M.append(Mi)
        #build the augmented matrix of the system G, an (ml^2)*(ml^2+1) matrix over F16,from FVV,M,digest,Spublic,salt
        G=build_the_augmented_matrix_of_the_system(snova_parameters,FVV,M,digest,Spublic,salt)    
        out,O=gauss(snova_parameters,G)
        #if out = true, the equation system represented by G is solved
        if out:
            is_done=True
            # update the list V to {V0,V1,...V{v-1},O0,O1,...,O{o-1}}
            V=V+O
            #obtain a column vector MV from V
            MR=MatrixSpace(R,1,len(V))
            MV=MR.matrix(V)
            #compute signature = [T](MV).t
            signature=T._multiply_strassen(MV.transpose(),1)
            #The following for loop converts signature to a list containing n elements in R.
            sig=[]
            for i in range(v+o):
                sig.append(signature[i,0])
        else:    
            #if out = false, the equation system represented by G is not solved 
            # numsig is increased by one to generate new random values V_0,V_1, ..., V_{v-1} in R in the next iteration
            is_done=False
            numsig=numsig+1 
            
         
    return (sig,salt)


def sign_fault_F2(snova_parameters,Spublic, Sprivate,digest,salt, I):
    
    (v,o,F,l) =snova_parameters
    m=o
    #Generate the random part of public key from Spublic
    _,_,_,A,B,Q1,Q2 = generate_the_random_part_of_public_key(snova_parameters,Spublic)
    #Generate the corresponding private_key from both Spublic and Sprivate
    T12,F11,F12,F21 = generate_private_key(snova_parameters,Spublic, Sprivate)
    # Build the augmented matrix [T]
    T=block_matrix_T(snova_parameters, T12)
    R=MatrixSpace(F,l)
    #Set numsig as 0
    numsig=0
    #Set is_done as False
    is_done=False
    while (not is_done):
        ## Asume all vallues for i in I are fixed. I is a subset of {0,1,2,..., v-1}
        V=assign_values_to_vinegar_variables_fault_F2(snova_parameters,I)
        #V=assign_values_to_vinegar_variables(snova_parameters,Sprivate,digest,salt,numsig)
        
        #Evaluate the vinegar part for each Fi. 
        FVV= compute_the_vinegar_part_of_the_central_map (snova_parameters,F11,A,B,Q1,Q2,V)
        #The following for loop computes Mik for all i,k in {0,1,...,m-1}. 
        #It computes M, a m*m matrix where M[i][k] stores Mik
        M=[]
        for i in range(m):
            Mi=[]
            for k in range(o):
                Mik=compute_the_coeﬀicient_matrix_of_the_oil_variable(snova_parameters,F12[i],F21[i],A,B,Q1,Q2,V,k)
                Mi.append(Mik)
                
            M.append(Mi)
        #build the augmented matrix of the system G, an (ml^2)*(ml^2+1) matrix over F16,from FVV,M,digest,Spublic,salt
        G=build_the_augmented_matrix_of_the_system(snova_parameters,FVV,M,digest,Spublic,salt)    
        out,O=gauss(snova_parameters,G)
        #if out = true, the equation system represented by G is solved
        if out:
            is_done=True
            # update the list V to {V0,V1,...V{v-1},O0,O1,...,O{o-1}}
            V=V+O
            #obtain a column vector MV from V
            MR=MatrixSpace(R,1,len(V))
            MV=MR.matrix(V)
            #compute signature = [T](MV).t
            signature=T._multiply_strassen(MV.transpose(),1)
            #The following for loop converts signature to a list containing n elements in R.
            sig=[]
            for i in range(v+o):
                sig.append(signature[i,0])
        else:    
            #if out = false, the equation system represented by G is not solved 
            # numsig is increased by one to generate new random values V_0,V_1, ..., V_{v-1} in R in the next iteration
            is_done=False
            numsig=numsig+1 
            
         
    return (sig,salt)

def evaluate_the_public_map(snova_parameters,A,B,Q1,Q2,P11,P12,P21,P22,U):
    
    (v,o,F, l)=snova_parameters
    m=o
    n=v+o
    ring=MatrixSpace(F,l)
    left=[]
    right=[]
    # First step
    for a in range(l*l):
        L=[]
        R=[]
        for j in range(n):
            L.append(A[a]*U[j].transpose()*Q1[a])
            R.append(Q2[a]*U[j]*B[a])
        left.append(L)
        right.append(R)
    out=[]  
    # Second step
    for i in range(m):
        out.append(ring.zero_matrix())
        for a in range(l*l):
            
            # Compute the sum terms depending on P11
            for dj in range(v):
                for dk in range(v):
                     out[i]= out[i]+left[a][dj]*P11[i][dj,dk]*right[a][dk]
            
            # Compute the sum terms depending on P12      
            for dj in range(v):
                for dk in range(o):
                     out[i]= out[i]+left[a][dj]*P12[i][dj,dk]*right[a][v+dk]
                    
            # Compute the sum terms depending on P21        
            for dj in range(o):
                for dk in range(v):
                     out[i]= out[i]+left[a][v+dj]*P21[i][dj,dk]*right[a][dk]
                    
            # Compute the sum terms depending on P22
            for dj in range(o):
                for dk in range(o):
                     out[i]= out[i]+left[a][v+dj]*P22[i][dj,dk]*right[a][v+dk]
        

    return out 

def signature_verification(snova_parameters,public_key, digest, signature ):
            (Spublic,P22)=public_key
            (sig,salt)=signature
            (v,o,F, l)=snova_parameters
            forder=F.order()
            m=o
            ring=MatrixSpace(F,l)
            # Generate the random part of public key from the public seed
            P11,P12,P21,A,B,Q1,Q2=generate_the_random_part_of_public_key(snova_parameters,Spublic)
            #Generation of a byte string from Spublic||digest|salt to derive the expected hash 
            byte_string=shake_256(Spublic+digest+salt).digest(m*l*l)
            hash1=[]
            for i in range(m):
                hh=[]
                for j in range(l*l):
                    hh.append(F.from_integer(byte_string[i*l*l+j]%forder))
                hash1.append(ring.matrix(hh))
            # Evaluate the public map P on input sig, i.e. hash2 = P[P_0[sig],P_1[sig],...,P_{m−1}[sig]]
            hash2=evaluate_the_public_map(snova_parameters,A,B,Q1,Q2,P11,P12,P21,P22,sig)      
            
            ## Comparison between hash1 and hash2. 
            acc=True
            for i in range(m):
                acc= acc and (hash1[i]==hash2[i])
            
            #acc will end up as True if and only if hash1[i]==hash2[i] for all i   
            
            return acc


q = 16
F = GF(q,'x')
FE=[]
dict_for_S={}
#Store all the elements from GF(16) in the list FE.
for i in range(16):
    FE.append(F.from_integer(i)) 
    # Construct the matrices S2,S3,S4 (l=2,3,4)
S2 = matrix(F,[[FE[8],FE[7]],[FE[7],FE[6]]])
S3 = matrix(F,[[FE[8],FE[7], FE[6]],[FE[7],FE[6], FE[5]], [FE[6],FE[5], FE[4]]])
S4 = matrix(F,[[FE[8],FE[7], FE[6],FE[5]], [FE[7],FE[6], FE[5],FE[4]],[FE[6],FE[5], FE[4],FE[3]],[FE[5],FE[4], FE[3],FE[2]]] )
#Define a dictionary to easy access to S2,S3, or S4.
dict_for_S={2:S2,3:S3,4:S4}
SL=[(28,17,F,2),(25,8,F,3),(24,5,F,4), (43, 25, F, 2),(49, 11, F, 3), (37, 8, F, 4), (61, 33, F, 2), (66, 15, F, 3), (60, 10, F, 4)]
snova_parameters = SL[0]
print(snova_parameters)
Sprivate= generate_seed(32)
Spublic= generate_seed(32)
digest=generate_seed(32)
salt=generate_seed(32)
(Spublic,P22)= generate_public_key (snova_parameters,Spublic, Sprivate)            

print("Pub_key: {} -- {}".format(Spublic, P22))
(sig, salt)=sign_original(snova_parameters,Spublic, Sprivate,digest,salt)
acc=signature_verification(snova_parameters,(Spublic,P22), digest, (sig, salt))
print("Passed?",acc==True)
