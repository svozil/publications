from fractions import Fraction
import numpy as np

def get_first_n_primes(n):
    """
    Generates the first n prime numbers.

    Args:
        n: The number of prime numbers to generate.

    Returns:
        A list of the first n prime numbers.
    """
    primes = []
    num = 2
    while len(primes) < n:
        is_prime = True
        # Check for factors from 2 up to the square root of num
        for i in range(2, int(num**0.5) + 1):
            if num % i == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(num)
        num += 1
    return primes

# (1) Input vectors
vectors = [
    (0, 0, 1, -1), (1, -1, 0, 0), (1, 1, -1, -1),
    (1, 1, 1, 1), (1, -1, 1, -1), (1, 0, -1, 0),
    (0, 1, 0, -1), (1, 0, 1, 0), (1, 1, -1, 1),
    (-1, 1, 1, 1), (1, 1, 1, -1), (1, 0, 0, 1),
    (0, 1, -1, 0), (0, 1, 1, 0), (0, 0, 0, 1),
    (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 1),
]

# (3) Get the first 18 prime numbers
primes = get_first_n_primes(18)

# Initialize a 4x4 numpy array to hold Fraction objects for the sum
sum_of_projectors = np.array([
    [Fraction(0) for _ in range(4)] for _ in range(4)
], dtype=object)

# (2) & (3) Calculate projectors, multiply by primes, and sum
for i, a_i in enumerate(vectors):
    # Get the corresponding prime number
    prime = primes[i]

    # Convert the vector components to Fraction objects
    v = [Fraction(x) for x in a_i]

    # Calculate the squared norm of the vector to keep arithmetic rational.
    squared_norm = sum(comp**2 for comp in v)

    if squared_norm == 0:
        continue # Skip the zero vector if it were present

    # Calculate the projector Pi = (v * v^T) / ||v||^2 and add its weighted value
    for row in range(4):
        for col in range(4):
            projector_element = (v[row] * v[col]) / squared_norm
            # Multiply by the prime number and add to the sum
            sum_of_projectors[row, col] += prime * projector_element

# Print the final resulting matrix with exact fractions
print("The sum of the projectors weighted by the first 18 prime numbers is:\n")
for row in sum_of_projectors:
    print([str(f) for f in row])

# (4) --- Eigensystem Calculation ---

# Convert the Fraction matrix to a standard float matrix for numerical computation
float_matrix = np.array(sum_of_projectors, dtype=float)

print("\n" + "="*50)
print("Calculating Eigensystem of the Final Matrix")
print("="*50)
print("\nFloating-point representation of the matrix for calculation:\n")
print(float_matrix)

# Calculate eigenvalues and eigenvectors
try:
    eigenvalues, eigenvectors = np.linalg.eig(float_matrix)

    print("\nCalculated Eigenvalues:\n")
    print(eigenvalues)

    print("\nCorresponding Eigenvectors (each column is an eigenvector):\n")
    print(eigenvectors)

    print("\n--- Eigenvalue and corresponding Eigenvector pairs ---\n")
    # The matrix is symmetric, so eigenvalues and eigenvectors will be real.
    for i in range(len(eigenvalues)):
        print(f"Eigenvalue {i+1}: {eigenvalues[i]:.4f}")
        print(f"Eigenvector {i+1}: {eigenvectors[:, i]}\n")

except np.linalg.LinAlgError as e:
    print(f"\nCould not compute the eigensystem: {e}")