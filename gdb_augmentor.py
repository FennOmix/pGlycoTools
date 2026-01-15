import os
import numpy as np
from numba import njit
import argparse

ORD_A = 65
ORD_G = 71

@njit(cache=True)
def get_popcount(n):
    count = 0
    while n > 0:
        n &= (n - 1)
        count += 1
    return count

@njit(cache=True)
def generate_variants_core(base_arr, loc_A, loc_G, max_mutation_n_A, max_mutation_n_G):
    n_A = len(loc_A)
    n_G = len(loc_G)

    valid_masks_A = []
    limit_A = n_A if max_mutation_n_A < 0 else min(n_A, max_mutation_n_A)
    for mask in range(1 << n_A):
        if get_popcount(mask) <= limit_A:
            valid_masks_A.append(mask)

    valid_masks_G = []
    limit_G = n_G if max_mutation_n_G < 0 else min(n_G, max_mutation_n_G)
    for mask in range(1 << n_G):
        if get_popcount(mask) <= limit_G:
            valid_masks_G.append(mask)

    total_variants = len(valid_masks_A) * len(valid_masks_G)

    result_block = np.empty((total_variants, len(base_arr)), dtype=np.uint8)

    row_idx = 0
    for ma in valid_masks_A:
        for mg in valid_masks_G:
            result_block[row_idx] = base_arr
            for bit_idx in range(n_A):
                if (ma >> bit_idx) & 1:
                    result_block[row_idx][loc_A[bit_idx]] = ORD_G
            for bit_idx in range(n_G):
                if (mg >> bit_idx) & 1:
                    result_block[row_idx][loc_G[bit_idx]] = ORD_A
            row_idx += 1
            
    return result_block

def run_augmentation(input_file, output_file, max_mutation_n_A=None, max_mutation_n_G=None):
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    if lines:
        header = lines[0]
        data_lines = lines[1:]
    else:
        return

    unique_structs = set()
    limit_A = -1 if max_mutation_n_A is None else max_mutation_n_A
    limit_G = -1 if max_mutation_n_G is None else max_mutation_n_G

    for line in data_lines:
        base_arr = np.frombuffer(line.encode('ascii'), dtype=np.uint8)
        loc_A = np.where(base_arr == ORD_A)[0]
        loc_G = np.where(base_arr == ORD_G)[0]

        variants_block = generate_variants_core(base_arr, loc_A, loc_G, limit_A, limit_G)

        for row in variants_block:
            unique_structs.add(row.tobytes().decode('ascii'))
            
    sorted_result = sorted(list(unique_structs), key=lambda x: (len(x), x))

    with open(output_file, 'w') as f:
        f.write(header + "\n")
        f.write("\n".join(sorted_result) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="gdb_augmentor",
        description="Augment glycan structures in GDB file by swapping A/G residues.",
    )
    
    parser.add_argument("-i", "--input_gdb", type=str, required=True, help="Path to input .gdb file")
    parser.add_argument("-o", "--output_gdb", type=str, required=True, help="Path to output .gdb file")

    parser.add_argument("-ma", "--max_mutation_a", type=int, default=2, help="Max number of 'A' to swap (default: 2)")
    parser.add_argument("-mg", "--max_mutation_g", type=int, default=2, help="Max number of 'G' to swap (default: 2)")
    
    args = parser.parse_args()
    
    run_augmentation(
        input_file=args.input_gdb,
        output_file=args.output_gdb,
        max_mutation_n_A=args.max_mutation_a,
        max_mutation_n_G=args.max_mutation_g,
    )