

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

def generate_invertible_matrix(snova_parameters,M):
    (v,o,F,l)=snova_parameters
     
    S=dict_for_S[l]
    zero=F.from_integer(0)
    if M.determinant()==zero:
        for i in range(1,F.order()):
            if (M+F.from_integer(i)*S).determinant()!=zero:
                M=M+F.from_integer(i)*S
                break
    
    return M

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
