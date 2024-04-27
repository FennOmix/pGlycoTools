import numpy as np

ASCII = 128

class BaseMass:
    def __init__(self):
        self.mass_H = 1.0078250321
        self.mass_O = 15.9949146221
        self.mass_N = 14.0030740052
        self.mass_C = 12.00
        self.mass_isotope = 1.003
        self.mass_proton = 1.007276

        self.mass_H2O = self.mass_H * 2 + self.mass_O
        self.mass_CO = self.mass_C + self.mass_O
        self.mass_CO2 = self.mass_C + self.mass_O * 2
        self.mass_NH = self.mass_N + self.mass_H
        self.mass_NH3 = self.mass_N + self.mass_H * 3
        self.mass_HO = self.mass_H + self.mass_O


class GlyIonCalc:
    def __init__(self, glyco_names, gly_dict):
        self.glyco_mass = []
        self.glyco_name = []
        for glyco in glyco_names:
            self.glyco_name.append(glyco)
            self.glyco_mass.append(gly_dict[glyco].mass)
        self.rand_decoy = 1
        self.rand_lower = 1
        self.rand_upper = 30

    def calc_glycan_mass(self, glycan):
        mass = 0.0
        for i in range(len(glycan)):
            mass += glycan[i] * self.glyco_mass[i]
        return mass

    def calc_Y_ions_no_pep(self, glycan_fragment, glycan_mass, isdecoy):
        Yions = np.zeros(len(glycan_fragment), dtype=np.float64)
        for j in range(len(glycan_fragment)):
            for i in range(len(self.glyco_mass)):
                Yions[j] += self.glyco_mass[i] * glycan_fragment[j][i]
        if isdecoy:
            if self.rand_decoy:
                np.random.seed(int(glycan_mass * 1000))
                Yions += (
                    np.random.random(len(Yions)) * (self.rand_upper - self.rand_lower)
                    + self.rand_lower
                )
            else:
                Yions = Yions[::-1]
        return Yions


class PepIonCalc:
    def __init__(self, mod_dict, aa_mass_dict=None):
        self.base_mass = BaseMass()
        self.AAMass = np.zeros(ASCII)
        if aa_mass_dict:
            for i in range(ord("A"), ord("Z") + 1):
                if chr(i) in aa_mass_dict:
                    self.AAMass[i] = aa_mass_dict[chr(i)].mass
        else:
            self.AAMass[ord("A")] = 71.037114
            self.AAMass[ord("C")] = 103.009185
            self.AAMass[ord("D")] = 115.026943
            self.AAMass[ord("E")] = 129.042593
            self.AAMass[ord("F")] = 147.068414
            self.AAMass[ord("G")] = 57.021464
            self.AAMass[ord("H")] = 137.058912
            self.AAMass[ord("I")] = 113.084064
            self.AAMass[ord("J")] = 114.042927
            self.AAMass[ord("K")] = 128.094963
            self.AAMass[ord("L")] = 113.084064
            self.AAMass[ord("M")] = 131.040485
            self.AAMass[ord("N")] = 114.042927
            self.AAMass[ord("P")] = 97.052764
            self.AAMass[ord("Q")] = 128.058578
            self.AAMass[ord("R")] = 156.101111
            self.AAMass[ord("S")] = 87.032028
            self.AAMass[ord("T")] = 101.047679
            self.AAMass[ord("V")] = 99.068414
            self.AAMass[ord("W")] = 186.079313
            self.AAMass[ord("Y")] = 163.06332

        self.ModMass = {}
        for modname, modobj in mod_dict.items():
            self.ModMass[modname] = modobj.mass

    def aamass(self, aa):
        return self.AAMass[ord(aa)]

    def set_aamass(self, aa, mass):
        self.AAMass[ord(aa)] = mass

    def calc_mod_mass_list(self, peptide, modinfo):
        modmass = np.zeros(len(peptide) + 2)

        items = modinfo.strip(";").split(";")
        modlist = []

        if modinfo == "" or modinfo == "null":
            return modmass

        for mod in items:
            strSite, modname = mod.split(",")
            site = int(strSite)
            modlist.append((site, modname))
        modlist.sort()

        for site, modname in modlist:
            _mono = self.ModMass[modname]
            modmass[site] = _mono
        return modmass

    def calc_b_ions(self, peptide, modmass):
        b_ions = np.zeros(len(peptide) + 2, dtype=np.float64)
        b_ions[1:-1] = self.AAMass[np.array(peptide, "c").view(np.int8)]
        b_ions += modmass
        return b_ions[1:-2]

    def calc_pepmass_from_b(self, peptide, modmass, bions):
        return (
            bions[len(bions) - 1]
            + self.aamass(peptide[len(bions)])
            + modmass[len(peptide)]
            + modmass[len(peptide) + 1]
            + self.base_mass.mass_H2O
        )

    def calc_y_from_b(self, bions, pepmass):
        return pepmass - bions

    def calc_a_from_b(self, bions):
        return bions - self.base_mass.mass_CO

    def calc_c_from_b(self, bions):
        return bions + self.base_mass.mass_NH3

    def calc_z_from_b(self, bions, pepmass):
        return (pepmass - self.base_mass.mass_NH3 + self.base_mass.mass_H) - bions

    def calc_H2O_loss(self, ions):
        return ions - self.base_mass.mass_H2O

    def calc_NH3_loss(self, ions):
        return ions - self.base_mass.mass_NH3
