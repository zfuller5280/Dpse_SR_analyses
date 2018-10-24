from optparse import OptionParser
import itertools
import operator as op
import random
import numpy as np

USAGE= """Usage: %prog [options] -i infile.txt"""
OPT_DEFAULTS={'infile':'-','bootstrap':1000}
DESCRIPTION="""Program description"""
EPILOG="""Requirements:"""

def get_options(defaults, usage, description='',epilog=''):
    """Get options, print usage text."""
    parser = OptionParser(usage=usage,description=description,epilog=epilog)
    parser.add_option("-i","--infile",action="store",dest="infile",type="string",
                      default=defaults.get('infile'),
                      help='Name of input file')
    parser.add_option("-N","--names",action="store",dest="names",type="string",
                      default=defaults.get('names'),
                      help='The name of each population, in order of the individuls in the population list')
    parser.add_option("-L","--locations",action="store",dest="locations",type="string",
                      default=defaults.get('locations'),
                      help='The start and end coordinates to estiamte the distance statistics over')
    parser.add_option("-w","--window_size",action="store",dest="window_size",type="int",
                      default=defaults.get('window'),
                      help='Size of window')
    parser.add_option("-b","--bootstrap",action="store",dest="bootstrap",type="int",
                      default=defaults.get('bootstrap'),
                      help='Number of bootstrap replicates')
    (options,args)=parser.parse_args()

    return (options, args)

def read_header(infile, names):
    with open(infile) as f:
        header = f.readline().strip("\r\n").split("\t")
    f.close()
    pos_col = header.index("Position")
    pop_idxs = [header.index(x) for x in names.split(",")]

    return pos_col, pop_idxs

def subset_data(infile, pos_col, locations):
    subset = []
    start, stop = locations
    with open(infile) as fp:
        fp.next()
        for line in fp:
            line = line.strip("\r\n").split("\t")
            position = int(float(line[pos_col]))
            if position >= start and position <= stop:
                subset.append(line)
    return subset

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom

def calc_fst(comb):
    pop1, pop2 = [int(x) for x in comb[0].split(",")], [int(x) for x in comb[1].split(",")]
    a1, a2, b1, b2 = pop1[1], pop2[1], pop1[0], pop2[0]
    n1, n2 = (a1 + b1), (a2 + b2)
    h1 = float(a1 * (n1 - a1))/(n1 * (n1 - 1))
    #
    h2_num, h2_denom = float(a2 * (n2 - a2)),float(n2 * (n2 - 1))
    #h2 = (float(a2 * (n2 - a2))) / float(n2 * (n2 - 1))
    #print h2_num, h2_denom, h2_num/h2_denom
    h2 = h2_num/h2_denom
    N = (((float(a1)/n1) - (float(a2)/n2))**2) - (h1/n1) - (h2/n2)
    if N == 0:
        return 0.0
    else:
        D = N + h1 + h2
        F = float(N)/D
        return F

def calc_dxy(comb):
    pop1, pop2 = [int(x) for x in comb[0].split(",")], [int(x) for x in comb[1].split(",")]
    a1, a2, b1, b2 = pop1[1], pop2[1], pop1[0], pop2[0]
    n1, n2 = (a1 + b1), (a2 + b2)
    x1 = float(b1)/n1
    x2 = float(a2)/n2
    dxy_i = x1 * x2
    dxy_j = (float(a1)/n1) * (float(b2)/n2)
    dxy = dxy_i + dxy_j
    return dxy

def get_freqs(window, pop_idxs, window_size):
    fst, dxy, diff = [], [], []

    for site in window:
        site_fsts, site_dxys = [], []
        alleles = [site[x] for x in pop_idxs]
        if "0,0" in alleles or "1,0" in alleles or "0,1" in alleles:
            continue
        #print alleles
        diffs = get_within_diffs(alleles)

        site_comb = list(itertools.combinations(alleles,2))
        for comb in site_comb:
            #print comb
            site_fsts.append(calc_fst(comb))
            site_dxys.append(calc_dxy(comb))
        fst.append(site_fsts)
        dxy.append(site_dxys)
        diff.append(diffs)
    window_diffs = [(float(x))/window_size for x in map(sum,zip(*diff))]
    window_fsts = [float(x)/len(window) for x in map(sum,zip(*fst))]
    window_dxy = [(float(x))/window_size for x in map(sum,zip(*dxy))]
    return window_fsts + window_dxy + window_diffs

def get_within_diffs(alleles):
    diffs = []
    for i in alleles:
        x, y = [int(j) for j in i.split(",")]
        p_diffs = x * y
        total = x + y
        #print
        diffs.append(float(p_diffs)/ncr(total,2))
    return diffs



def read_windows(infile, window_size, pos_col):
    start_pos = int(float(infile[0][pos_col]))
    window, start, stop, pos = [], start_pos , start_pos+window_size, None
    for line in infile:
        pos = int(float(line[pos_col]))
        if pos > stop:
            yield window
            start += window_size
            stop += window_size
            window = []
            window.append(line)
        else:
            window.append(line)

def bootstrap_window(window, pop_idxs, window_size, replicates):
    bootstrapped_stats = []
    for replicate in xrange(replicates):
        samp_window = [random.choice(window) for _ in window]
        window_fsts = get_freqs(samp_window, pop_idxs, window_size)
        bootstrapped_stats.append(window_fsts)
    cis = []
    #print bootstrapped_stats
    num_stats = len(bootstrapped_stats[0])
    for stat in xrange(num_stats):
        vals = [x[stat] for x in bootstrapped_stats]
        cis.append([str(ci) for ci in [np.percentile(vals,2.5),np.percentile(vals,97.5)]])
    return [item for sublist in cis for item in sublist]


def main():
    (options,args) = get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)
    infile = options.infile
    names = options.names
    locations = options.locations
    window_size = options.window_size
    bootstrap_replicates = options.bootstrap

    stats = ["fst","dxy"]
    comps = list(itertools.combinations(names.split(","),2))

    names_list = names.split(",")


    header = ["w_start","w_stop"]
    for stat in stats:
        for comp in comps:

            comp = "-".join([comp[0],comp[1]])
            lbl = "_".join([comp, stat])
            header.append(lbl)
    header += ["_".join([x,"diffs"]) for x in names_list]
    ci_header = []
    for i in header[2:]:
        ci_header.append("_".join([i,"ci_l"]))
        ci_header.append("_".join([i,"ci_u"]))
    header += ci_header
    locations = [int(x) for x in locations.split(",")]
    pos_col, pop_idxs = read_header(infile, names)

    subset = subset_data(infile, pos_col, locations)

    print "\t".join(header)
    for window in read_windows(subset, window_size, pos_col):
        if len(window) <= 2:
            continue
        start, stop = int(float(window[0][pos_col])), int(float(window[-1][pos_col]))
        window_fsts = get_freqs(window, pop_idxs, window_size)
        window_cis = bootstrap_window(window, pop_idxs, window_size,bootstrap_replicates)
        output = [str(x) for x in [start, stop]]
        #print output
        #print window_cis
        print "\t".join(output + [str(x) for x in window_fsts] + [str(x) for x in window_cis])
        #break


if __name__ == '__main__':
    main()
