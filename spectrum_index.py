import numpy as np

BIN = 0.05


def PeakIndexing(peaklist, tol, isppm):
    MAX_BIN = int((peaklist[-1, 0] + 1) / BIN)
    if isppm:
        _tol = peaklist[:, 0] * tol * 1e-6
    else:
        _tol = np.full(len(peaklist), tol, dtype=np.float64)
    index_table = np.full(MAX_BIN, -1, dtype=np.int32)
    mass_lower_list = ((peaklist[:, 0] - _tol) / BIN).astype(np.int32)
    mass_upper_list = ((peaklist[:, 0] + _tol) / BIN + 1).astype(np.int32)
    for i in range(
        len(peaklist) - 1, -1, -1
    ):  # make sure index_table record the lower mass peak
        for j in range(mass_lower_list[i], mass_upper_list[i] + 1):
            index_table[j] = i
    return index_table


def MatchOneMass(peaklist, index_table, mass, tol, isppm):
    imass = int(mass / BIN)
    if imass >= len(index_table):
        return -1
    peakidx = index_table[imass]
    if peakidx == -1:
        return peakidx
    ret = -1
    min_delta = 1e100
    if isppm:
        tol = mass * tol * 1e-6
    while peakidx < len(peaklist):
        delta = peaklist[peakidx, 0] - mass
        if delta > tol:
            break
        elif delta < -tol:
            peakidx += 1
        else:
            if min_delta > abs(delta):
                ret = peakidx
                min_delta = abs(delta)
            peakidx += 1
    return ret


def Match(peaklist, index_table, masslist, tol, isppm):
    ret_idx = np.full(len(masslist), -1, dtype=int)
    ret_tol = np.full(len(masslist), tol * 2, dtype=np.float64)
    ret_inten = np.zeros(len(masslist), dtype=np.float64)

    for i in range(len(masslist)):
        ret_idx[i] = MatchOneMass(peaklist, index_table, masslist[i], tol, isppm)
        if ret_idx[i] != -1:
            if isppm:
                ret_tol[i] = (peaklist[ret_idx[i], 0] - masslist[i]) * 1e6 / masslist[i]
            else:
                ret_tol[i] = peaklist[ret_idx[i], 0] - masslist[i]
            ret_inten[i] = peaklist[ret_idx[i], 1]
    return ret_idx, ret_tol, ret_inten


def isMatch(peaklist, index_table, masslist, tol, isppm):
    ret = np.zeros(len(masslist), dtype=int)
    for i in range(len(masslist)):
        is_matched = MatchOneMass(peaklist, index_table, masslist[i], tol, isppm)
        if is_matched != -1:
            ret[i] = 1
    return ret
