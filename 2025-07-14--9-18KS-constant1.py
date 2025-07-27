from fractions import Fraction
import numpy as np

# (1) Input vectors
vectors = [
    (0, 0, 1, -1), (1, -1, 0, 0), (1, 1, -1, -1),
    (1, 1, 1, 1), (1, -1, 1, -1), (1, 0, -1, 0),
    (0, 1, 0, -1), (1, 0, 1, 0), (1, 1, -1, 1),
    (-1, 1, 1, 1), (1, 1, 1, -1), (1, 0, 0, 1),
    (0, 1, -1, 0), (0, 1, 1, 0), (0, 0, 0, 1),
    (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 1),
]

# Initialize a 4x4 numpy array to hold Fraction objects for the sum
sum_of_projectors = np.array([
    [Fraction(0) for _ in range(4)] for _ in range(4)
], dtype=object)

# (2) & (3) Calculate projectors and sum them (weighting by a constant 1)
for a_i in vectors:
    # Convert the vector components to Fraction objects
    v = [Fraction(x) for x in a_i]

    # Calculate the squared norm of the vector to keep arithmetic rational.
    squared_norm = sum(comp**2 for comp in v)

    if squared_norm == 0:
        continue # Skip the zero vector if it were present

    # Calculate the projector Pi = (v * v^T) / ||v||^2 and add its value to the sum
    for row in range(4):
        for col in range(4):
            projector_element = (v[row] * v[col]) / squared_norm
            # Add to the sum (implicitly multiplying by the constant factor 1)
            sum_of_projectors[row, col] += projector_element

# Print the final resulting matrix with exact fractions
print("The sum of the projectors with a constant factor of 1 is:\n")
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