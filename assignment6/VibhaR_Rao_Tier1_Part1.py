import sys
import numpy as np
from scipy.linalg import eig

def fibonacci_nth_term(a, b, n):
    # Define the Fibonacci matrix
    A = np.array([[1, 1], [1, 0]])
    
    # Calculate eigenvalues using scipy.linalg.eig and ignore eigenvectors
    eigenvalues, _ = eig(A)
    
    # Extract the real parts of the eigenvalues
    lambda1, lambda2 = eigenvalues.real

    # Use the known eigenvector structure for the Fibonacci matrix
    v1 = np.array([lambda1, 1])
    v2 = np.array([lambda2, 1])
    
    # Initial vector with custom starting values
    v0 = np.array([b, a])
    
    # Calculate coefficients for the closed-form solution using the initial vector

    # To find c1c_1c1​ and c2c_2c2​, we need to "project" v0 onto each eigenvector. This is done by taking the dot product of 
    #v0 with each eigenvector and normalizing by the dot product of the eigenvector with itself. 

    c1 = np.dot(v0, v1) / np.dot(v1, v1)
    c2 = np.dot(v0, v2) / np.dot(v2, v2)
    
    # Calculate the nth term using the derived formula
    nth_term = c1 * (lambda1 ** (n - 1)) + c2 * (lambda2 ** (n - 1))
    
    # Round to the nearest integer, as Fibonacci terms are integers
    return int(round(nth_term))

if __name__ == "__main__":
    # Take arguments from the command line
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    n = int(sys.argv[3])
    
    # Calculate and print the nth term
    print(fibonacci_nth_term(a, b, n))
