load("utils.sage")
load("snova_faulty.sage")

from numpy import random


def get_F16(snova_parameters, prob):
    (v, o, field, l) = snova_parameters

    I = []
    for i in range(v * l * l):
        I.append(field.random_element())

    x = random.binomial(n=1, p=prob, size=v * l * l)

    return (x, I)


def get_F2(snova_parameters, prob):
    (v, o, field, l) = snova_parameters

    I = []
    number_bits = int(log(field.order(), 2))
    F2 = field.base_ring()
    for i in range(number_bits * v * l * l):
        I.append(F2.random_element())

    x = random.binomial(n=1, p=prob, size=number_bits * v * l * l)

    return (x, I)

def reconstruct_element_in_F16_S(snova_parameters,scalars):
    (v,o,F,l)=snova_parameters
    
    
    if  len(scalars)!= l:
        return 
    
    R = MatrixSpace(F, l)    
    S=dict_for_S[l] 
    M = R.identity_matrix()
    O=R.matrix()
    
    for i in range(l):
        O=O+scalars[i]*M
        M=M*S
    
    return O



def recover_partial_T_ring_element(snova_parameters, signatures):
    (v,o,F, l)=snova_parameters
    R=MatrixSpace(F,l)
    S=copy(dict_for_S[l])    
    PS=[R.identity_matrix(), S]
    for k3 in range(1,l-1,1):
        PS.append(PS[k3]*S)
    
    output={}
    
    Nmsg=len(signatures)
    
    for i in range(v):
        SS=matrix(F,(Nmsg-1)*l*l,l*o)
        YY=matrix(F,(Nmsg-1)*l*l,1)
        for k in range(1,Nmsg,1):
            A=signatures[0][i]-signatures[k][i]
            for r0 in range(l):
                for r1 in range(l):
                    for j in range(o):
                        for j1 in range(l):
                            B=PS[j1]*(signatures[k][v+j]-signatures[0][v+j])
                            SS[(k-1)*l*l+r0*l+r1,j*l+j1]=B[r0,r1]
                    
                    YY[(k-1)*l*l+r0*l+r1,0]=A[r0,r1]
        try:
            X= SS\YY
            Ti=[]
            for t in range(o):
                scalars=[X[i,0] for i in range(t*l, (t+1)*l,1)]
                Tij=reconstruct_element_in_F16_S(snova_parameters,scalars)
                Ti.append(Tij)
                
            output[i]=Ti
                    
            
        except Exception as ee:
            output[i]= None
            
        
    
    return output



def recover_partial_T_F16(snova_parameters, signatures):
    (v,o,F, l)=snova_parameters
    R=MatrixSpace(F,l)
    S=copy(dict_for_S[l])    
    PS=[R.identity_matrix(), S]
    for k3 in range(1,l-1,1):
        PS.append(PS[k3]*S)
    
    output={}
    
    Nmsg=len(signatures)
    for i in range(v):
        for r0 in range(l):
            for r1 in range(l):
                SS=matrix(F,(Nmsg-1),l*o)
                YY=matrix(F,(Nmsg-1),1)
                for k in range(1,Nmsg,1):
                    A=signatures[0][i]-signatures[k][i]
                    for j in range(o):
                        for j1 in range(l):
                            B=PS[j1]*(signatures[k][v+j]-signatures[0][v+j])
                            SS[(k-1),j*l+j1]=B[r0,r1]
                    
                    YY[(k-1),0]=A[r0,r1]
                
                try:
                    X= SS\YY
                    Ti=[]
                    for t in range(o):
                        scalars=[X[i,0] for i in range(t*l, (t+1)*l,1)]
                        Tij=reconstruct_element_in_F16_S(snova_parameters,scalars)
                        Ti.append(Tij)
                
                    output[(i,r0,r1)]=Ti
                    
            
                except Exception as ee:
                    output[(i,r0,r1)]= None
            
        
    
    return output


def recover_partial_T_F2(snova_parameters, signatures):
    (v,o,F, l)=snova_parameters
    R=MatrixSpace(F,l)
    S=copy(dict_for_S[l])    
    PS=[R.identity_matrix(), S]
    for k3 in range(1,l-1,1):
        PS.append(PS[k3]*S)
    
    output={}
    
    Nmsg=len(signatures)
    FF=F.base_ring()
    for i in range(v):
        for r0 in range(l):
            for r1 in range(l):
               
                for z in range(4):
                    SS=matrix(FF,(Nmsg-1),4*l*o)
                    YY=matrix( FF,(Nmsg-1),1)
                    for k in range(1,Nmsg,1):
                        C=signatures[0][i]-signatures[k][i]
                    
                        for j in range(o):
                             for j1 in range(l):
                                Ar=PS[j1]*(signatures[k][v+j]-signatures[0][v+j])
                                A=Ar[r0,r1]
                                
                                if z==0:
                                    #c0=a[0]*b[0]+a[1]*b[3]+a[2]*b[2]+a[3]*b[1]
                                    SS[(k-1),j*l+j1*4]=A[0]
                                    SS[(k-1),j*l+j1*4+1]=A[3]
                                    SS[(k-1),j*l+j1*4+2]=A[2]
                                    SS[(k-1),j*l+j1*4+3]=A[1]
                                    
                                elif z==1:
                                    #c1=(a[2]+a[1])*b[3]+(a[2]+a[3])*b[2]+(a[0]+a[3])*b[1]+a[1]*b[0]
                                    
                                    SS[(k-1),j*l+j1*4]=A[1]
                                    SS[(k-1),j*l+j1*4+1]=A[3]+A[0]
                                    SS[(k-1),j*l+j1*4+2]=A[3]+A[2]
                                    SS[(k-1),j*l+j1*4+3]=A[2]+A[1]
    
                                elif z==2:
                                    #c2=(a[3]+a[2])*b[3]+(a[0]+a[3])*b[2]+a[1]*b[1]+a[2]*b[0]
                                    SS[(k-1),j*l+j1*4]=A[2]
                                    SS[(k-1),j*l+j1*4+1]=A[1]
                                    SS[(k-1),j*l+j1*4+2]=A[0]+A[3]
                                    SS[(k-1),j*l+j1*4+3]=A[3]+A[2]
                                else:
                                    #c3=(a[3]+a[0])*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0]
                                    
                                    SS[(k-1),j*l+j1*4]=A[3]
                                    SS[(k-1),j*l+j1*4+1]=A[2]
                                    SS[(k-1),j*l+j1*4+2]=A[1]
                                    SS[(k-1),j*l+j1*4+3]=A[3]+A[0]
                    
                        YY[(k-1),0]=C[r0,r1][z]
                
                    try:
                        X= SS\YY
                        Ti=[]
                        for t1 in range(o):
                            scalars=[]
                            for t2 in range(l):
                                bits=[]
                                for t3 in range(4):
                                     bits.append(X[t1*l+t2*4+t3,0])
                                scalars.append(bits)
                            Tij=reconstruct_element_in_F16_S(snova_parameters,scalars)
                            Ti.append(Tij)
                
                        output[(i,r0,r1,z)]=Ti
                    
            
                    except Exception as ee:
                        output[(i,r0,r1,z)]= None
            
        
    
    return output

def test_F16(params, prob):

    (v,o,F, l)=snova_parameters
    print(snova_parameters)
    Sprivate= generate_seed(32)
    Spublic= generate_seed(32)
    (Spublic,P22)=generate_public_key (snova_parameters,Spublic, Sprivate)
    (x,I)=get_F16(snova_parameters,prob)
    print(x)
    #print("After the fault,",u2," vinegar values are fixed")
    #print("I:=", I.keys())
    signatures=[]
    Nmsg=o*l +2
    for pp in range(Nmsg):
        digest=generate_seed(32)
        salt=generate_seed(32)
        (sig, _) =sign_fault_F16(snova_parameters,Spublic, Sprivate,digest,salt, (x,I))
        signatures.append(sig)
        print("signature",pp, "Received")

    print("Num of signatures", len(signatures))
    output=recover_partial_T_F16(snova_parameters, signatures)
    T12,_,_,_ = generate_private_key(snova_parameters,Spublic, Sprivate)
    R=MatrixSpace(F,l)
    ST=MatrixSpace(R,1,o)
    LL=[]
    dd=0
    for key in output.keys():

        if output[key] != None:
            row=ST.matrix(output[key])
            if row==T12[key[0],:]:
                dd=dd+x[key[0]*l*l+key[1]*l+key[2]]

    print(dd, sum(x))
    for i in range(v):
        print(i, (sum(x[i*l*l:(i+1)*l*l])>=1))


def test_F2(params, prob):

    (v,o,F, l)=snova_parameters
    print(snova_parameters)
    Sprivate= generate_seed(32)
    Spublic= generate_seed(32)
    (Spublic,P22)=generate_public_key (snova_parameters,Spublic, Sprivate)
    (x,I)=get_F2(snova_parameters, prob)
    print(x)
    #print("After the fault,",u2," vinegar values are fixed")
    #print("I:=", I.keys())
    signatures=[]
    Nmsg=4*o*l +2
    for pp in range(Nmsg):
        digest=generate_seed(32)
        salt=generate_seed(32)
        (sig, _) =sign_fault_F2(snova_parameters,Spublic, Sprivate,digest,salt, (x,I))
        signatures.append(sig)
        print("signature",pp, "recieved")

    print("Num of signatures", len(signatures))
    output=recover_partial_T_F16(snova_parameters, signatures)
    T12,_,_,_ = generate_private_key(snova_parameters,Spublic, Sprivate)
    R=MatrixSpace(F,l)
    ST=MatrixSpace(R,1,o)
    LL=[]
    dd=0
    for key in output.keys():

        if output[key] != None:
            row=ST.matrix(output[key])
            if row==T12[key[0],:]:
                dd=dd+x[key[0]*l*l+key[1]*l+key[2]*4 +key[3]]

    print(dd, sum(x))
    for i in range(v):
        print(i, (sum(x[4*i*l*l:4*(i+1)*l*l])>=1))


def get_list_matrix_reconciliation(snova_parameters,public_key):
    (v,o,F, l)=snova_parameters
    (Spublic,P22)=public_key
    P11,P12,P21,A,B,Q1,Q2=generate_the_random_part_of_public_key(snova_parameters,Spublic)
    n=v+o
    R=MatrixSpace(F,l)
    S=copy(dict_for_S[l])
    PS=[R.identity_matrix(), S]
    for k3 in range(1,l-1,1):
        PS.append(PS[k3]*S)

    RR=MatrixSpace(R,n)
    MR=[]
    for i in range(o):
        Pi=block_matrix([[P11[i], P12[i]],[P21[i], P22[i]]], subdivide=False)
        for j in range(l):
            for k in range(l):
                rs=RR.matrix()
                for r0 in range(n):
                    for r1 in range(n):
                        rs[r0,r1]=PS[j]*Pi[r0,r1]*PS[k]

                MR.append(rs)


    Fln=MatrixSpace(F,l*n)
    output=[]
    for i in range(len(MR)):
        rs=MR[i]
        rsfln=Fln.matrix()
        for r0 in range(n):
            for j0 in range(l):
                for r1 in range(n):
                    for j1 in range(l):
                        rsfln[r0*l+j0, r1*l+j1]=rs[r0,r1][j0,j1]
        output.append(rsfln)

    return output

def check_vector(list_matrix_reconciliation, vector):

    output=True
    for mat in list_matrix_reconciliation:
        output= output and (vector.transpose()*mat*vector==0)

    return output


def get_secret_space(snova_parameters,list_matrix_reconciliation, vec):

    C=[]
    (v,o,F, l)=snova_parameters
    output=check_vector(list_matrix_reconciliation, vec)
    print("u_0",output)
    C.append(vector(vec[:, 0]))
    R=MatrixSpace(F,l)
    n=v+o
    Fln=MatrixSpace(F,l*n, 1)
    S=copy(dict_for_S[l])
    PS=[R.identity_matrix(), S]
    for k3 in range(1,l-1,1):
        PS.append(PS[k3]*S)

    print("constructing other elements")
    for k3 in range(1, l):
        rs=Fln.matrix()
        for i in range(n):
            rs[i*l:(i+1)*l, 0]=PS[k3]*vec[i*l:(i+1)*l, 0]

        output=check_vector(list_matrix_reconciliation, rs)
        print("u_"+str(k3), output)
        C.append(vector(rs[:, 0]))

    nrows=2*len(list_matrix_reconciliation)
    ncolumns=n*l
    MM=matrix(F,nrows,ncolumns)
    for i in range(len(list_matrix_reconciliation)):
        mat=list_matrix_reconciliation[i]
        col=mat*vec
        MM[2*i,:]=col.transpose()
        row=vec.transpose()*mat
        MM[2*i+1,:]=row

    kernel_space=MM.right_kernel()
    i=2
    while i<=l*o:

        C.append(kernel_space.random_element())
        W = span(C)
        if W.dimension()<i:
            C= C[:-1]
        else:
            output=check_vector(list_matrix_reconciliation, rs)
            print("u_"+str(i), output)
            i=i+1

    return span(C)

def random_vector_secret_space_snova(snova_parameters, O):
    (v,o,F, l)=snova_parameters
    Fql=VectorSpace(F,l)
    n=v+o
    vec_from_O=O.random_element()
    vec=Fql.random_element()
    R=MatrixSpace(F,l)
    out=[]
    for i in range(n):
        mm=R.matrix()
        for j in range(l):
            rowj=vec_from_O[i*l+j]*vec
            for k in range(l):
                mm[j,k]=rowj[k]
        out.append(mm)

    return out
