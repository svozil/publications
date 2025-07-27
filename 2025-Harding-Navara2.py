# -*- coding: utf-8 -*-
import numpy as np
import itertools

# --- 0. Setup and Vector Definitions ---

np.set_printoptions(precision=3, suppress=True)
omega = np.exp(2 * np.pi * 1j / 3)
omega2 = omega**2

# Vector definitions remain the same as before
vector_definitions = {
    "a1": lambda x, y, z: np.array([1, 0, 0], dtype=complex), "a2": lambda x, y, z: np.array([0, 1, 0], dtype=complex), "a3": lambda x, y, z: np.array([0, 0, 1], dtype=complex),
    "u":  lambda x, y, z: np.array([x, y, z]), "b1": lambda x, y, z: np.array([0, y, z]), "b2": lambda x, y, z: np.array([x, 0, z]), "b3": lambda x, y, z: np.array([x, y, 0]),
    "c1": lambda x, y, z: np.array([0, z.conjugate(), -y.conjugate()]), "c2": lambda x, y, z: np.array([z.conjugate(), 0, -x.conjugate()]), "c3": lambda x, y, z: np.array([y.conjugate(), -x.conjugate(), 0]),
    "d1": lambda x, y, z: np.array([-y*y.conjugate() - z*z.conjugate(), x.conjugate()*y, x.conjugate()*z]), "d2": lambda x, y, z: np.array([x*y.conjugate(), -x*x.conjugate() - z*z.conjugate(), y.conjugate()*z]), "d3": lambda x, y, z: np.array([x*z.conjugate(), y*z.conjugate(), -x*x.conjugate() - y*y.conjugate()]),
    "b12": lambda x, y, z: np.array([(y*z).conjugate(), (x*z).conjugate(), -(x*y).conjugate()]), "b112": lambda x, y, z: np.array([x*y*y.conjugate() + x*z*z.conjugate(), -y*z*z.conjugate(), y*y.conjugate()*z]), "b212": lambda x, y, z: np.array([-x*z*z.conjugate(), x*x.conjugate()*y + y*z*z.conjugate(), x*x.conjugate()*z]),
    "b13": lambda x, y, z: np.array([(y*z).conjugate(), -(x*z).conjugate(), (x*y).conjugate()]), "b113": lambda x, y, z: np.array([x*y*y.conjugate() + x*z*z.conjugate(), y*z*z.conjugate(), -y*y.conjugate()*z]), "b313": lambda x, y, z: np.array([-x*y*y.conjugate(), x*x.conjugate()*y, x*x.conjugate()*z + y*y.conjugate()*z]),
    "b23": lambda x, y, z: np.array([-(y*z).conjugate(), (x*z).conjugate(), (x*y).conjugate()]), "b223": lambda x, y, z: np.array([x*z*z.conjugate(), y*z*z.conjugate() + x*x.conjugate()*y, -x*x.conjugate()*z]), "b323": lambda x, y, z: np.array([x*y*y.conjugate(), -x*x.conjugate()*y, x*x.conjugate()*z + y*y.conjugate()*z]),
}

def mesh_and_analyze():
    """
    Generates vectors for all cases and finds inter-case orthogonalities.
    """
    print(f"\n{'='*70}")
    print("ANALYSIS OF INTER-CASE ORTHOGONALITY")
    print("This program checks if a vector from one case is orthogonal to one from another case.")
    print(f"{'='*70}\n")
    
    cases = {
        "Case1_Real":    {'x': 1+0j, 'y': 1+0j, 'z': 1+0j, 'desc': "x=1, y=1, z=1"},
        "Case2_Omega":   {'x': 1+0j, 'y': omega, 'z': omega, 'desc': "x=1, y=w, z=w"},
        "Case3_OmegaSq": {'x': 1+0j, 'y': omega, 'z': omega2, 'desc': "x=1, y=w, z=w^2"}
    }
    
    # --- Step 1: Generate and store all vectors from all cases ---
    all_evaluated_vectors = {}
    for case_name, params in cases.items():
        x, y, z = params['x'], params['y'], params['z']
        for vec_name, func in vector_definitions.items():
            # The key is a tuple of (vector_name, case_name)
            key = (vec_name, case_name)
            all_evaluated_vectors[key] = func(x, y, z)
            
    # --- Step 2: Find vectors that are identical across different cases ---
    print("--- Section 1: Vectors That Are Identical Across Cases ---")
    identical_pairs = []
    vector_keys = list(all_evaluated_vectors.keys())
    
    for key1, key2 in itertools.combinations(vector_keys, 2):
        # We only care about the same vector name in different cases
        if key1[0] != key2[0] or key1[1] == key2[1]:
            continue
        
        v1 = all_evaluated_vectors[key1]
        v2 = all_evaluated_vectors[key2]
        
        # np.allclose checks if two arrays are element-wise equal within a tolerance
        if np.allclose(v1, v2):
            identical_pairs.append(f"  {key1[0]:<5} in {key1[1]} == {key2[0]:<5} in {key2[1]}")

    if not identical_pairs:
        print("  No vectors (besides a1,a2,a3) were found to be identical across different cases.")
    else:
        print("  The following vectors have identical values in different cases:")
        for pair in sorted(identical_pairs):
            print(pair)
            
    # --- Step 3: Find non-trivial inter-case orthogonalities ---
    print("\n--- Section 2: Non-Trivial Inter-Case Orthogonalities ---")
    
    inter_case_orthogonalities = []
    
    for key1, key2 in itertools.combinations(vector_keys, 2):
        # Condition 1: Must be from different cases
        if key1[1] == key2[1]:
            continue
            
        # Condition 2: Exclude the constant basis vectors a1, a2, a3
        if key1[0] in {'a1', 'a2', 'a3'} or key2[0] in {'a1', 'a2', 'a3'}:
            continue
            
        v1 = all_evaluated_vectors[key1]
        v2 = all_evaluated_vectors[key2]
        
        # Condition 3: Check for orthogonality
        if np.isclose(np.dot(v1, v2.conjugate()), 0):
            pair_string = f"  {key1[0]:<5} ({cases[key1[1]]['desc']})  _|_  {key2[0]:<5} ({cases[key2[1]]['desc']})"
            inter_case_orthogonalities.append(pair_string)
            
    if not inter_case_orthogonalities:
        print("  No non-trivial inter-case orthogonalities were found.")
    else:
        print("  Found the following orthogonal pairs between different cases:")
        for pair in sorted(inter_case_orthogonalities):
            print(pair)

# --- Main Execution ---
if __name__ == "__main__":
    mesh_and_analyze()