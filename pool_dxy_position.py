from optparse import OptionParser
import pandas as pd
import random
import numpy as np

USAGE= """Usage: %prog [options]"""
OPT_DEFAULTS={'infile':'-'}
DESCRIPTION="""Program description: """
EPILOG="""Requirements:"""

random.seed(16)

def get_options(defaults, usage, description='',epilog=''):
    """Get options, print usage text."""
    parser=OptionParser(usage=usage,description=description,epilog=epilog)
    parser.add_option("-i","--infile",action="store",dest="infile",type="string",
                      default=defaults.get('infile'),
                      help='Name of input gene')
    parser.add_option("-m","--mpileup",action="store",dest="mpileup",type="string",
                      default=defaults.get('mpileup'),
                      help='Name of mpileup file containing per base coverage')
    parser.add_option("-w","--window_size",action="store",dest="window_size",type="int",
                      default=defaults.get('window'),
                      help='Size of window')
    parser.add_option("-p","--positions",action="store",dest="positions",type="string",
                      default=defaults.get('positions'),
                      help='Central position of the inversion bp to calculate dXY around')
    parser.add_option("-s","--samples",action="store",dest="samples",type="string",
                      default=defaults.get('samples'),
                      help='Samples, two pooled populations separated by comma')
    parser.add_option("-c","--coverages",action="store",dest="coverages",type="string",
                      default=defaults.get('coverages'),
                      help='Coverage cutoffs, the cutoffs of coverage for the two pooled populations separated by comma')
    (options,args)=parser.parse_args()

    return (options, args)

def read_vcf_header(infile, pop1, pop2):
    with open(infile) as fp:
        for line in fp:
            if line.startswith("#CHROM"):
                line = line.strip("\r\n").split("\t")
                for field in line:
                    if pop1 in field:
                        idx1 = line.index(field)
                    if pop2 in field:
                        idx2 = line.index(field)
                break

    return idx1, idx2

def read_vcf_windows(infile, window_size, position):
    if window_size % 2 == 0:
        window_size += 1
    window, pos_L = [], []
    position_low, position_high = position-(float(window_size)/2), position + (float(window_size)/2)
    start, stop = min(position_low, position_high), max(position_low, position_high)
    with open(infile) as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            line = line.strip("\r\n").split("\t")
            chrom, pos = line[0], int(line[1])

            if pos > stop:
                last_pos = int(window[-1][1])
                yield chrom, window, start, last_pos
                break
            elif pos > start:
                window.append(line)
                pos_L.append(pos)
            else:
                continue

def get_coverage(pop1, pop2, scaffold):
    infile1 = "%s.%s.mpileup"%(scaffold,pop1)
    infile2 = "%s.%s.mpileup"%(scaffold,pop2)
    coverage_df1 = pd.read_csv(infile1, sep="\t", header=None)
    coverage_df2 = pd.read_csv(infile2, sep="\t", header=None)
    coverage_df1.columns, coverage_df2.columns = ['chrom', 'pos', 'coverage_%s'%pop1], ['chrom', 'pos', 'coverage_%s'%pop2]
    coverage_df = coverage_df1.merge(coverage_df2, on=["chrom", "pos"], how='outer')
    coverage_df = coverage_df.sort_values('pos')
    coverage_df['coverage_%s'%pop1] = coverage_df['coverage_%s'%pop1].fillna(0)
    coverage_df['coverage_%s'%pop2] = coverage_df['coverage_%s'%pop2].fillna(0)

    return coverage_df

def get_callable_sites(coverage_df, start, stop, coverage1, coverage2, pop1, pop2):
    coverage_window = coverage_df[(coverage_df['pos']>=start) & (coverage_df['pos']<=stop)]
    coverage_window = coverage_window[(coverage_window['coverage_%s'%pop1]>=coverage1) & (coverage_window['coverage_%s'%pop2]>=coverage2)]
    callable_sites = len(coverage_window.index)

    return callable_sites

def make_consensus_call(geno, coverage, ref, alt):
    r = random.uniform(0.0,1.0)
    ad = geno.split(":")[1].split(",")
    ref_reads, alt_reads = [int(x) for x in ad]
    depth = ref_reads + alt_reads
    if depth < coverage or geno.split(":")[0]==".":
        call = None
    else:
        alt_freq = alt_reads/depth
        if r <= alt_freq:
            call = alt
        else:
            call = ref
    return call

def get_consensus_calls(window, idx1, idx2, coverage1, coverage2):
    matches, diffs, missing_sites = 0, 0, 0
    for site in window:
        ref, alt = site[3], site[4]
        geno1, geno2 = site[idx1], site[idx2]
        call1, call2 = make_consensus_call(geno1, coverage1, ref, alt), make_consensus_call(geno2, coverage2, ref, alt)
        if call1 != None and call2 != None:
            if call1 == call2:
                matches += 1
            elif call1 != call2:
                diffs += 1
            else:
                print("Something went wrong!")
        else:
            missing_sites += 1
    return matches, diffs, missing_sites

def read_positions(positions):
    position_list = []
    with open(positions) as fp:
        for line in fp:
            line = line.strip("\r\n").split("\t")
            breakpoint, infile, position = line
            position_list.append([breakpoint, infile, position])
    return position_list

def main():
    (options,args)=get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)
    infile = options.infile
    window_size = options.window_size
    mpileup = options.mpileup
    samples = options.samples
    coverages = options.coverages
    positions = options.positions


    positions = read_positions(positions)

    pop1, pop2 = samples.split(",")
    coverage1, coverage2 = [int(x) for x in coverages.split(",")]

    windows = []
    total_matches, total_diffs, total_missing_sites, total_raw_callable_sites, total_callable_sites = 0, 0, 0, 0, 0
    for position in positions:
        infile = "%s.merged.vcf.gz.indel_filtered.vcf"%position[1]
        scaffold = infile.split(".")[0]
        for chrom, window, start, stop in read_vcf_windows(infile, window_size, int(position[2])):

            coverage_df = get_coverage(pop1, pop2, scaffold)
            idx1, idx2 = read_vcf_header(infile, pop1, pop2)
            windows.append(window)
            callable_sites = get_callable_sites(coverage_df, start, stop, coverage1, coverage2, pop1, pop2)
            total_raw_callable_sites += callable_sites
            matches, diffs, missing_sites = get_consensus_calls(window, idx1, idx2, coverage1, coverage2)
            total_matches += matches
            total_diffs += diffs
            total_missing_sites += missing_sites
    breakpoint = positions[-1][0].split("_")[0]
    windows = [window for x in windows for window in x]
    #print(total_matches, total_diffs, total_missing_sites, total_raw_callable_sites)
    total_callable_sites = total_raw_callable_sites
    total_callable_sites -= total_missing_sites
    dxy = total_diffs/total_callable_sites

    bootstrapped_dxy = []
    for x in range(1000):
        window = [random.choice(windows) for _ in range(len(windows))]
        matches, diffs, missing_sites = get_consensus_calls(windows, idx1, idx2, coverage1, coverage2)
        rep_callable_sites = total_raw_callable_sites - missing_sites
        rep_dxy = diffs/rep_callable_sites
        bootstrapped_dxy.append(rep_dxy)
    low_ci = np.percentile(bootstrapped_dxy,2.5)
    hi_ci = np.percentile(bootstrapped_dxy,97.5)

    print(breakpoint, total_raw_callable_sites-total_missing_sites, diffs, dxy, low_ci, hi_ci)
        #break

if __name__ == '__main__':
    main()
