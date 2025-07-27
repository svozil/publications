import numpy as np
from collections import defaultdict
import itertools

# (1) Define the "complex cross product" function
def complex_cross_product(u, v):
    """
    Calculates the complex cross product as defined:
    conj(standard_cross_product(u, v))
    """
    # The standard cross product in numpy
    standard_cross = np.cross(u, v)
    # Return the complex conjugate of the result
    return np.conj(standard_cross)

# Define a function to pretty-print vectors
def format_vector(v):
    """Formats a numpy vector for clean printing."""
    # Format each component to show real and imag parts cleanly
    # Handles cases like (1+0j) -> 1, (0+2j) -> 2j, (1+1j) -> 1+1j
    return f"({', '.join(f'{c.real if c.imag == 0 else c}' for c in v)})"

print("--- Python Program for Vector System Analysis ---")

# (2.1) Identify x = y = z = 1
# We use complex numbers to be rigorous, e.g., 1+0j
x = 1 + 0j
y = 1 + 0j
z = 1 + 0j

# Pre-calculate conjugates for readability
xc = np.conj(x)
yc = np.conj(y)
zc = np.conj(z)

# (2.2) Evaluate all vectors
# We store them in a dictionary for easy access
vectors = {}
dtype = np.complex128 # Use complex numbers

# --- Input the system of vectors ---
vectors['a1'] = np.array([1, 0, 0], dtype=dtype)
vectors['a2'] = np.array([0, 1, 0], dtype=dtype)
vectors['a3'] = np.array([0, 0, 1], dtype=dtype)
vectors['u'] = np.array([x, y, z], dtype=dtype)
vectors['b1'] = np.array([0, y, z], dtype=dtype)
vectors['b2'] = np.array([x, 0, z], dtype=dtype)
vectors['b3'] = np.array([x, y, 0], dtype=dtype)
vectors['c1'] = np.array([0, zc, -yc], dtype=dtype)
vectors['c2'] = np.array([zc, 0, -xc], dtype=dtype)
vectors['c3'] = np.array([yc, -xc, 0], dtype=dtype)
vectors['d1'] = np.array([-y*yc - z*zc, xc*y, xc*z], dtype=dtype)
vectors['d2'] = np.array([x*yc, -x*xc - z*zc, yc*z], dtype=dtype)
vectors['d3'] = np.array([x*zc, y*zc, -x*xc - y*yc], dtype=dtype)
vectors['b12'] = np.array([np.conj(y*z), np.conj(x*z), -np.conj(x*y)], dtype=dtype)
vectors['b112'] = np.array([x*y*yc + x*z*zc, -y*z*zc, y*yc*z], dtype=dtype)
vectors['b212'] = np.array([-x*z*zc, x*xc*y + y*z*zc, x*xc*z], dtype=dtype)
vectors['b13'] = np.array([np.conj(y*z), -np.conj(x*z), np.conj(x*y)], dtype=dtype)
vectors['b113'] = np.array([x*y*yc + x*z*zc, y*z*zc, -y*yc*z], dtype=dtype)
vectors['b313'] = np.array([-x*y*yc, x*xc*y, x*xc*z + y*yc*z], dtype=dtype)
vectors['b23'] = np.array([-np.conj(y*z), np.conj(x*z), np.conj(x*y)], dtype=dtype)
vectors['b223'] = np.array([x*z*zc, y*z*zc + x*xc*y, -x*xc*z], dtype=dtype)
vectors['b323'] = np.array([x*y*yc, -x*xc*y, x*xc*z + y*yc*z], dtype=dtype)

# Calculate the new vectors using the complex cross product
vectors['b12c3'] = complex_cross_product(vectors['b12'], vectors['c3'])
vectors['b13c2'] = complex_cross_product(vectors['b13'], vectors['c2'])
vectors['b23c1'] = complex_cross_product(vectors['b23'], vectors['c1'])

# --- Print the evaluated vectors ---
print("\n(2.2) Evaluated Vectors (for x=y=z=1):\n")
for name, vec in vectors.items():
    print(f"{name:<6} = {format_vector(vec)}")

# (2.4) Identify possible multiplicities
print("\n" + "-"*40)
print("\n(2.4) Vector Multiplicities (Identical Vectors):\n")
# Group vectors by their value (converting to tuple to be hashable)
value_to_names = defaultdict(list)
for name, vec in vectors.items():
    # Round to handle potential minor floating point differences
    vec_tuple = tuple(np.round(v, 8) for v in vec)
    value_to_names[vec_tuple].append(name)

# Create a clean list of unique vectors for the next step
unique_vectors = {
    value_to_names[val][0]: np.array(val, dtype=dtype)
    for val in value_to_names
}

# Print the groups of identical vectors
multiplicity_found = False
for val, names in value_to_names.items():
    if len(names) > 1:
        print(f"Value: {format_vector(val)}")
        print(f"  --> Names: {', '.join(sorted(names))}\n")
        multiplicity_found = True

if not multiplicity_found:
    print("No multiplicities found. All vectors are unique.")

# (2.3) List mutual orthogonalities
print("\n" + "-"*40)
print("\n(2.3) Mutual Orthogonality Cliques\n")
print("Orthogonality is defined by the Hermitian dot product u·v = Σ uᵢconj(vᵢ) = 0.")
print("We use np.vdot(v, u) for this calculation.\n")

names = list(unique_vectors.keys())
vecs = list(unique_vectors.values())
n = len(names)
adjacency = defaultdict(list)

# Build the orthogonality graph
for i in range(n):
    for j in range(i + 1, n):
        # The Hermitian dot product in numpy is vdot(u, v) = sum(conj(u)*v)
        # We check if u and v are orthogonal: vdot(v,u) == 0
        dot_product = np.vdot(vecs[j], vecs[i])
        if np.isclose(dot_product, 0):
            adjacency[names[i]].append(names[j])
            adjacency[names[j]].append(names[i])

# Bron-Kerbosch algorithm to find all maximal cliques
def find_cliques(potential_nodes, candidates, excluded, cliques):
    if not candidates and not excluded:
        cliques.append(potential_nodes)
        return

    # To avoid duplicate cliques, iterate over a copy of candidates
    for node in list(candidates):
        new_potential_nodes = potential_nodes + [node]
        new_candidates = [n for n in candidates if n in adjacency[node]]
        new_excluded = [n for n in excluded if n in adjacency[node]]
        find_cliques(new_potential_nodes, new_candidates, new_excluded, cliques)
        
        candidates.remove(node)
        excluded.append(node)

all_cliques = []
find_cliques([], names, [], all_cliques)

# Filter for cliques of size 2 or more and print
orthogonality_cliques = [c for c in all_cliques if len(c) >= 2]

if not orthogonality_cliques:
    print("No sets of mutually orthogonal vectors found.")
else:
    # Sort cliques by size (largest first) then alphabetically
    orthogonality_cliques.sort(key=lambda c: (-len(c), sorted(c)))
    print(f"Found {len(orthogonality_cliques)} mutual orthogonality cliques:\n")
    for i, clique in enumerate(orthogonality_cliques):
        print(f"Clique {i+1} (Size {len(clique)}):")
        print(f"  {{ {', '.join(sorted(clique))} }}\n")