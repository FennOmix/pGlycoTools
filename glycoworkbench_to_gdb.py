import copy
import os
import pandas as pd

import argparse

gwb_glyco_map = {
    "1D-GlcNAc": "N",
    "1D-GalNAc": "N",
    "1D-Man": "H",
    "1D-Glc": "H",
    "1D-Gal": "H",
    "1L-Fuc": "F",
    "1D-Xyl": "X",
    "2D-NeuAc": "A",
    "2D-NeuGc": "G",
}


def parse_gwp(
    input_gwp_file,
    output_gdb_file=None,
    max_glycan_len=100,
    use_old_gdb_format=True,
):
    print(f"Parsing {input_gwp_file} ...")
    gwb_struct_list = load_gwp(input_gwp_file)
    canons = gwb_structs_to_canons(gwb_struct_list)
    print(f"Loaded {len(gwb_struct_list)} gwb structures ...")

    canons = list(generate_subtrees_by_canons(canons))
    print(f"Generated {len(canons)} pGlyco sub-structures ...")
    canons.sort(key=lambda x: (len(x), x))

    print(f"Save as {os.path.abspath(output_gdb_file)}")
    df = pd.DataFrame(dict(structure_canon=canons))
    df["composition"] = df.structure_canon.apply(
        lambda x: composition_dict_to_str(get_glyco_compositions(x))
    )
    df["n_glyco"] = df.structure_canon.str.count("\\(")
    df = df.query(f"n_glyco <= {max_glycan_len} and n_glyco > 0")
    if use_old_gdb_format:
        _save_old_gdb(df, gdb_file=output_gdb_file)
    else:
        df.to_csv(output_gdb_file, index=False, sep="\t")
    return df


def _save_old_gdb(df, gdb_file):
    glyco_units = df.structure_canon.apply(get_glyco_units)
    glyco_units = set.union(*glyco_units.values)
    with open(gdb_file, "w") as f:
        f.write(",".join(list(glyco_units)) + "\n")
        for canon in df.structure_canon.values:
            f.write(canon + "\n")


def load_gwp(gwp):
    gwb_struct_list = []
    with open(gwp) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith("<Glycan structure="):
                glystruct = line[line.find('"') + 1 : line.rfind('"')]
                gwb_struct_list.append(glystruct)
    return gwb_struct_list


def gwb_structs_to_canons(gwb_struct_list):
    canons = set()
    for gwb_struct in gwb_struct_list:
        canons.update([gwb_struct_to_canon(gwb_struct)])
    return canons


def gwb_struct_to_canon(gwb_struct):
    find_dict = {}
    for key in gwb_glyco_map.keys():
        find_dict[key] = str_find_all_substr(gwb_struct, key + ",p")
    left_ = str_find_all_substr(gwb_struct, "(")
    right_ = str_find_all_substr(gwb_struct, ")")
    items = []
    for key, vals in find_dict.items():
        for val in vals:
            items.append((val, gwb_glyco_map[key]))
    for val in left_:
        items.append((val, "("))
    for val in right_:
        items.append((val, ")"))
    items.sort()
    items = [item[1] for item in items]
    canon = _gwb_items_to_canon(items)
    return canon


def _gwb_items_to_canon(items):
    root_items = []
    branches = []
    branch_count = 0
    for i in range(len(items)):
        if items[i] == "(":
            break
        root_items.append(items[i])

    for i in range(i, len(items)):
        if items[i] != "(":
            break
        branch_count += 1

    start = i
    left_parent_count = branch_count
    for i in range(start, len(items)):
        if items[i] == ")":
            if left_parent_count == branch_count:
                branches.append(items[start:i])
                start = i + 1
                branch_count -= 1
                if branch_count == 0:
                    branches.append(items[i + 1 :])
                    break
            left_parent_count -= 1
        elif items[i] == "(":
            left_parent_count += 1

    merged_codes = []
    full_root = "(" + "(".join(root_items)
    full_right_ = ")" * len(root_items)
    for branch in branches:
        merged_codes.append(_gwb_items_to_canon(branch))
    merged_codes.sort(key=lambda x: (len(x), x))
    canon = full_root + "".join(merged_codes) + full_right_
    return canon


def str_find_all_substr(s, sub):
    ret = []
    idx = s.find(sub)
    while idx != -1:
        ret.append(idx)
        idx = s.find(sub, idx + 1)
    return ret


# print(generate_subtree_by_canon("(N(H(A))(H(A)))"))


def get_glyco_compositions(canon):
    items = [item.strip(")") for item in canon[1:].split("(")]
    comp_dict = {}
    for item in items:
        if item in comp_dict:
            comp_dict[item] += 1
        else:
            comp_dict[item] = 1
    return comp_dict


def composition_dict_to_str(comp_dict):
    return "".join(f"{g}({n})" for g, n in comp_dict.items())


def get_glyco_units(canon):
    return set([item.strip(")") for item in canon[1:].split("(")])


def generate_subtrees_by_canons(canons):
    canon_set = set()
    for canon in canons:
        canon_set.update(generate_subtrees_by_canon(canon))
    return canon_set


def generate_subtrees_by_canon(canon):
    items = []
    start = 0
    for i in range(len(canon)):
        if canon[i] == "(" or canon[i] == ")":
            if start < i:
                items.append(canon[start:i])
            items.append(canon[i])
            start = i + 1
    return canon_items_to_subtree_canons(items)


def canon_items_to_subtree_canons(items):
    root = items[1]
    branches = []
    left_count = 0
    start = 2
    for i in range(2, len(items) - 1):
        if items[i] == "(":
            left_count += 1
        elif items[i] == ")":
            left_count -= 1
            if left_count == 0:
                branches.append(items[start : i + 1])
                start = i + 1

    branch_canons = set([""])
    for branch in branches:
        branch_codes = canon_items_to_subtree_canons(branch)
        tmp_set = copy.deepcopy(branch_canons)
        for subcode in branch_codes:
            for merge_canon in tmp_set:
                if (len(merge_canon), merge_canon) < (len(subcode), subcode):
                    branch_canons.add(merge_canon + subcode)
                else:
                    branch_canons.add(subcode + merge_canon)
    rooted_subcanons = set()
    for _code in branch_canons:
        rooted_subcanons.add("(" + root + _code + ")")
    rooted_subcanons.add("")
    return rooted_subcanons


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="glycoworkbench_to_gdb",
        description="Convert glycoworkbench gwp file into pGlyco's gdb file",
    )
    parser.add_argument("-i", "--gwp_file", type=str)
    parser.add_argument("-o", "--gdb_file", type=str)
    parser.add_argument("-m", "--max_glycan_len", type=int, default=100)
    parser.add_argument("--use_old_gdb_format", type=bool, default=True)
    args = parser.parse_args()
    parse_gwp(
        input_gwp_file=args.gwp_file,
        output_gdb_file=args.gdb_file,
        max_glycan_len=args.max_glycan_len,
        use_old_gdb_format=args.use_old_gdb_format,
    )
