import numpy as np

def unique(arr):
    un = []
    for elem in arr:
        if elem not in un:
            un.append(elem)
    return un

def transpose(matrix):
    # Calculate number of rows and columns
    rows = len(matrix)
    cols = len(matrix[0])  # Assuming all rows have the same number of columns
    
    # Initialize an empty matrix for the transpose
    transpose = [[0 for _ in range(rows)] for _ in range(cols)]
    
    # Compute the transpose
    for i in range(rows):
        for j in range(cols):
            transpose[j][i] = matrix[i][j]
    
    return transpose

# multiplies 2 matrices mod p and returns the resulting matrix
def multiplyModP(matrix1, matrix2, p):
    result = [[0 for _ in range(len(matrix2[0]))] for _ in range(len(matrix1))]
    for i in range(len(matrix1)):
        for j in range(len(matrix2[0])):
            for k in range(len(matrix2)):
                result[i][j] += (matrix1[i][k] * matrix2[k][j])
    for i in range(len(result)):
        for j in range(len(result[0])):
            result[i][j] = result[i][j] % p
    return result


def matrix_vector_multiply(matrix, vector):
    # Get dimensions of the matrix and vector
    m = len(matrix)
    n = len(matrix[0])  # Number of columns in the matrix
    vector_length = len(vector)
    
    # Check if dimensions are compatible for multiplication
    if n != vector_length:
        raise ValueError("Matrix columns must match vector length for multiplication")
    
    # Initialize result vector with zeros
    result = [0] * m
    
    # Perform matrix-vector multiplication
    for i in range(m):
        for j in range(n):
            result[i] += matrix[i][j] * vector[j]
    return result

def mult(c, v):
    result = []
    for i in range(len(v)):
        result.append(c*v[i])
    return result



def add(v1, v2):
    # Ensure both vectors have the same length
    if len(v1) != len(v2):
        raise ValueError("Vectors must have the same length for addition")
    
    # Perform vector addition
    result = [v1[i] + v2[i] for i in range(len(v1))]
    
    return result  

# takes in a matrix in SL2(Fp) and returns its inverse
def inverse(matrix, p):
    result = [[0,0],[0,0]]
    result[0][0] = matrix[1][1]
    result[0][1] = (-matrix[0][1] + p) % p
    result[1][0] = (-matrix[1][0] + p) % p
    result[1][1] = matrix[0][0]
    return result

def cosetsGmodB(p):
    cosetsGmodB = []
    for i in range(p):
        cosetsGmodB.append([[1, 0], [i, 1]])
    cosetsGmodB.append([[0, 1],[-1, 0]])    
    return cosetsGmodB

# given a p, returns an array (of size p(p-1)(p+1))containing all elements of SL2(Fp) 
def SL2Fp(p):
    result = []
    for a in range(p):
        for b in range(p):
            for c in range(p):
                for d in range(p):
                    # Create the 2x2 matrix
                    matrix = [[a % p, b % p],
                              [c % p, d % p]]
                    # Calculate determinant
                    determinant = (a*d - b*c) % p
                    # Check if determinant is 1 mod p
                    if determinant == 1:
                        result.append(matrix)
    return result


# given g in SL2(Fp), returns which coset G/B it lives in
def whichCoset(g, p):
    for coset in cosetsGmodB(p):
        temp = multiplyModP(inverse(coset, p), g, p)
        if temp[1][0] == 0:
            return coset

# inputs a vector v and writes it as a linCombo of the basis vectors for the steinberg rep
def linCombo(v, p):
    basis_vectors = []
    for i in range(p):
        basis_vectors.append([])
    for vector in basis_vectors:
        for i in range(p+1):
            vector.append(0)
    for i in range(p):
        basis_vectors[i][i] = 1
        basis_vectors[i][i+1] = -1
    
    
    #basis vectors holds the steinberg basis vectors
    
    # should return coeffs, an array of length p which will be the column of the pxp matrix rep 
    
    #TODO write v as lin combo of basis_vectors
    #v will always be a p+1 length array!
    # loop through all possible combinations of basis_vectors: 3p times because coeff either 0, 1, -1 and there are p vectors/coeffs
    # if sum == v then extract coeffs
    options = create_perms([0, 1, -1], p)
    for coeffs in options:
        temp = []
        for i in range(len(v)):
            temp.append(0)
        for c in range(len(coeffs)):
            temp = add(temp, mult(coeffs[c], basis_vectors[c]))
        if temp == v:
            return coeffs
        
     


def create_perms(options, p):
    def perm_recursive(perms, k):
        #print(perms)
        #print(k)
        if (k == 0):
            return perms
        new_perms = []
        #print(perms)
        for perm in perms:
            for c in options:
                #print(perm)
                new_perms.append(perm + [c])
        return perm_recursive(new_perms, k-1)

    return perm_recursive([[]], p)


#print(linCombo([0, 1, 0, -1], 3)) #should return 0, 1, 1

# returns the p x p matrix of the linear rep of g
def linearRep(g, p):
    result = []
    StBasisVectors = []
    for e in cosetsGmodB(p):
        transformedCoset = whichCoset(multiplyModP(g, e, p),p)
        StBasisVectors.append(transformedCoset[1][0])
    #now convert the weil -1 to a p
    for i in range(len(StBasisVectors)):
        if StBasisVectors[i] == -1:
            StBasisVectors[i] = p
    #StBasisVectors is the set of p vectors which form a basis for steinberg rep 
    transformedStBasisVectors = []
    for i in range(p): # there will be p steinberg basis vectors 
        transformedStBasisVectors.append([[StBasisVectors[i]], [StBasisVectors[i+1]]])
    #this array holds the transformed basis vectors as e_i - e_j


    linearRep = []
    
    for basisVector in transformedStBasisVectors:
        # now for each transformed basis vector we want to write it as a (p+1) column as coeffs of the og basis vectors
        # and then write that collumn as a lin combo of "transformedBasisVectors"
        tempColumn = [] # tempColumn will be the coeffs of the og, which we will then write as a linCombo
        #print(basisVector)
        for i in range(p+1):
            tempColumn.append(0)
            if basisVector[0][0] == i:
                tempColumn[i] = 1
            if basisVector[1][0] == i:
                tempColumn[i] = -1
        #print(tempColumn)
        coeffs = linCombo(tempColumn, p)
        linearRep.append(coeffs)
        #now temp column is the proper coeffs for the transformed basis vectors
        # need to write it as a lincombo of the og (StBasisVectors)

        #linearRep.append(coeffs)
                

    
    
    return(transpose(linearRep))

def print_matrix_from_columns(columns):
    # Transpose the columns to get rows
    rows = list(zip(*columns))
    
    # Print the matrix
    for row in rows:
        print(" ".join(map(str, row)))



#input a list of vectors in the orbit and will return a list of the unique lines spanned by these vectors
def unique_lines(orbit):
    # so need to check that each vector is not a scalar multiple of any vector in the orbit 
    return lines

    
p = 3
G = SL2Fp(p)
all_lin_reps = []

for g in G:
    for row in g:
        print(row)
    print("")
    linRep = linearRep(g, p)
    all_lin_reps.append(linRep)
    for row in linRep:
        print(row)
    print("")
print("There are " + str(len(all_lin_reps)) + " elements in G, and " + str(len(unique(all_lin_reps))) + " unique linear transformations.")

v = [1, 3, 2]
orbit = []

print("Orbit of v = "+ str(v))
for g in all_lin_reps:
    print(g)
    print(matrix_vector_multiply(g, v))
    orbit.append(matrix_vector_multiply(g, v))
    print("")
print("There are " + str(len(unique(orbit))) + " unique vectors in the orbit")

print(unique(orbit))
