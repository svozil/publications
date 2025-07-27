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

# Initialize the product matrix as a 4x4 identity matrix with Fraction objects
product_matrix = np.array([
    [Fraction(1), Fraction(0), Fraction(0), Fraction(0)],
    [Fraction(0), Fraction(1), Fraction(0), Fraction(0)],
    [Fraction(0), Fraction(0), Fraction(1), Fraction(0)],
    [Fraction(0), Fraction(0), Fraction(0), Fraction(1)],
], dtype=object)

# Loop through each vector to calculate its projector and multiply it into the product
for a_i in vectors:
    # --- (2) Calculate the orthogonal projector Pi ---
    v = [Fraction(x) for x in a_i]
    squared_norm = sum(comp**2 for comp in v)

    if squared_norm == 0:
        continue # Skip if vector is zero

    # Create the projector matrix for the current vector
    P_i = np.array([[Fraction(0)]*4 for _ in range(4)], dtype=object)
    for row in range(4):
        for col in range(4):
            P_i[row, col] = (v[row] * v[col]) / squared_norm

    # --- (3) Multiply the current projector into the overall product ---
    # The @ operator performs matrix multiplication
    product_matrix = product_matrix @ P_i

# Print the final product matrix with exact fractions
print("The final matrix from multiplying all 18 projectors is:\n")
for row in product_matrix:
    print([str(f) for f in row])

# --- (4) Compute the Eigensystem of the product matrix ---

# Convert the Fraction matrix to a standard float matrix for numerical computation
float_matrix = np.array(product_matrix, dtype=float)

print("\n" + "="*50)
print("Calculating Eigensystem of the Final Product Matrix")
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
    for i in range(len(eigenvalues)):
        # Use np.real_if_close to display complex numbers as real if they are close
        eigenvalue_display = np.real_if_close(eigenvalues[i])
        print(f"Eigenvalue {i+1}: {eigenvalue_display:.4f}")

        eigenvector_display = np.real_if_close(eigenvectors[:, i])
        print(f"Eigenvector {i+1}: {eigenvector_display}\n")

except np.linalg.LinAlgError as e:
    print(f"\nCould not compute the eigensystem: {e}")