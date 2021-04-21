""" Helper functions for GAP """

def relabels_dets_to_disjoint_cycles(relabels_dets, show_progress = 0):
    """ Converts the relabellings on deterministic points to disjoint cycles """
    all_cycles = []

    for i, relabel in enumerate(relabels_dets):
        if show_progress:
            print('Get cycle {} / {}'.format(i, len(relabels_dets)))
        # convert to dict
        d = {}
        for i in range(relabel.shape[0]):
            d[i] = relabel[i]
        # generate cycles
        cycles = []
        while len(d) > 0:
            cycle = []
            key = list(d)[0]
            val = d.pop(key)
            if val != key:
                cycle.append(key)
                cycle.append(val)
            while val in list(d):
                key = val
                val = d.pop(key)
                if val not in cycle:
                    cycle.append(val)
            if cycle:
                cycles.append(cycle)
        if cycles:
            all_cycles.append(cycles)

    # Store the cycles
    out = ""
    for cycles in all_cycles:
        for cycle in cycles:
            out = out + "("
            for val in cycle:
                out += str(val + 1) + ","
            out = out[:-1] + ")"
        out = out + ",\n"
    out = out[:-2]
    return out