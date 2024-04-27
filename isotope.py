from emass import EMassSimple
from pGlyco_config import pGlycoConfig
import numpy as np


class IsotopeSimple:
    def __init__(self, config, max_isotope_len=7):
        self.max_isotope_len = max_isotope_len
        self._emass = EMassSimple(max_isotope_len)
        self.config = config

        self.chem_list = ["C", "H", "N", "O", "S", "P"]
        self.chem_idx = dict(zip(self.chem_list, range(len(self.chem_list))))

        self.modlist = [
            self.config.mod_dict[mod] for mod in self.config.protein_var_mod if mod
        ]
        self.modlist += [
            self.config.mod_dict[mod] for mod in self.config.protein_fix_mod if mod
        ]

        self.replace_element = {}
        self.check_element()

    def check_element(self, replace_element=None):
        if replace_element is not None:
            self.replace_element = replace_element
        self._check_chem(self.modlist)
        self._check_chem(self.config.aa_dict.values())
        self._check_chem(self.config.gly_dict.values())

        self._chem_array_dict()

    def _chem_tuples(self, chem):
        items = chem.strip(")").split(")")
        items = [item.split("(") for item in items]
        return [
            (
                elem
                if elem not in self.replace_element
                else self.replace_element[elem],
                int(n),
            )
            for elem, n in items
        ]

    def _check_chem(self, obj_list):
        for obj in obj_list:
            vec = self._chem_tuples(obj.chemical)
            for elem, n in vec:
                if elem not in self.chem_list:
                    self.chem_idx[elem] = len(self.chem_list)
                    self.chem_list.append(elem)

    def _chem_array_dict(self):
        self.mod_dict = {}
        self.aa_dict = {}
        self.gly_dict = {}

        def _set_dict(xx_dict, obj_list):
            for obj in obj_list:
                arr = np.zeros(len(self.chem_list), dtype=int)
                vec = self._chem_tuples(obj.chemical)
                idx = [self.chem_idx[elem] for elem, n in vec]
                num = [n for elem, n in vec]
                arr[idx] = num
                xx_dict[obj.name] = arr

        _set_dict(self.mod_dict, self.modlist)
        _set_dict(self.aa_dict, self.config.aa_dict.values())
        _set_dict(self.gly_dict, self.config.gly_dict.values())

        self.chem_H2O = np.zeros(len(self.chem_list), dtype=int)
        self.chem_H2O[self.chem_idx["H"]] = 2
        self.chem_H2O[self.chem_idx["O"]] = 1

    def get_distribution(self, peptide, modifications, glycan_comp, glyco_names):
        chem = self.chem_H2O.copy()
        for aa in peptide:
            chem += self.aa_dict[aa]
        if modifications:
            modlist = modifications.strip(";").split(";")
            modlist = [mod.split(",")[1] for mod in modlist]
            for mod in modlist:
                chem += self.mod_dict[mod]

        for i in range(len(glycan_comp)):
            chem += self.gly_dict[glyco_names[i]] * glycan_comp[i]

        self.chem = chem

        return self._emass.get_distribution(list(zip(self.chem_list, chem)))


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    import sys

    argd = {}
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i + 1]

    mono = argd["-mono"]
    if mono == "13C":
        config = pGlycoConfig(r"d:\GlycoData\yeast_15N\pGlyco3aHx2\pGlycoXIC-13C.cfg")
    elif mono == "15N":
        config = pGlycoConfig(r"d:\GlycoData\yeast_15N\pGlyco3aHx2\pGlycoXIC-15N.cfg")
    else:
        config = pGlycoConfig(r"d:\GlycoData\yeast_15N\pGlyco3aHx2\pGlycoXIC-none.cfg")

    peptide = argd["-seq"]
    modification = argd["-mod"] if "-mod" in argd else ""
    glycan = argd["-glycan"]
    glycos = glycan.strip(")").split(")")
    glyco_names = []
    glycan_comp = []
    for glyco in glycos:
        name, comp = glyco.split("(")
        glyco_names.append(name)
        glycan_comp.append(int(comp))

    isocalc = IsotopeSimple(config, 10)
    yint, mono_idx = isocalc.get_distribution(
        peptide, modification, glycan_comp, glyco_names
    )

    plt.vlines(list(range(len(yint))), [0] * len(yint), yint)
    plt.vlines(mono_idx, [0], yint[mono_idx], color="red")
    plt.show()
