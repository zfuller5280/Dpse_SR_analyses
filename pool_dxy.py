from optparse import OptionParser
import pandas as pd
import random

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

def read_vcf_windows(infile, window_size):
    if window_size % 2 == 0:
        window_size += 1
    window, pos_L = [], []
    start, stop = 1, window_size
    with open(infile) as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            line = line.strip("\r\n").split("\t")
            chrom, pos = line[0], int(line[1])

            if pos > stop:
                last_pos = int(window[-1][1])
                yield chrom, window, start, last_pos
                window, pos_L = [], []
                window.append(line)
                pos_L.append(pos)
                start += window_size
                stop += window_size
            else:
                window.append(line)
                pos_L.append(pos)

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
    if depth < coverage:
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


def main():
    (options,args)=get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)
    infile = options.infile
    window_size = options.window_size
    mpileup = options.mpileup
    samples = options.samples
    coverages = options.coverages

    scaffold = infile.split(".")[0]

    pop1, pop2 = samples.split(",")
    coverage1, coverage2 = [int(x) for x in coverages.split(",")]

    idx1, idx2 = read_vcf_header(infile, pop1, pop2)

    coverage_df = get_coverage(pop1, pop2, scaffold)

    for chrom, window, start, stop in read_vcf_windows(infile, window_size):
        callable_sites = get_callable_sites(coverage_df, start, stop, coverage1, coverage2, pop1, pop2)
        matches, diffs, missing_sites = get_consensus_calls(window, idx1, idx2, coverage1, coverage2)
        callable_sites -= missing_sites
        dxy = diffs/callable_sites
        print(chrom, start, stop, callable_sites, diffs, dxy)
        #break

if __name__ == '__main__':
    main()
