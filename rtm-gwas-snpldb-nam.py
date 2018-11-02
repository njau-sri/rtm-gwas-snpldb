import sys

def read_blocks(infile):
    chr, start, stop = [], [], []
    with open(infile) as f:
        for line in f:
            v = line.split()
            chr.append(v[0])
            start.append(int(v[1]))
            stop.append(int(v[2]))
    return chr, start, stop

def read_snp_genotype(infile):
    loc, chr, pos = [], [], []
    geno = []
    with open(infile) as f:
        ind = f.readline().split()[4:]
        for line in f:
            v = line.split()
            loc.append(v[0])
            chr.append(v[1])
            pos.append(int(v[2]))
            geno.append(v[4:])
    return ind, loc, chr, pos, geno

def index_block_snps(blk_chr, blk_start, blk_stop, chr, pos):
    idx = []
    for i in range(len(chr)):
        if chr[i] == blk_chr and pos[i] >= blk_start and pos[i] <= blk_stop:
            idx.append(i)
    return idx

def index_family_rils(ns):
    idx = []
    idx.append(range(ns[0]))
    for i in range(1,len(ns)):
        a = sum(ns[:i])
        idx.append( [e + a for e in range(ns[i])] )
    return idx

# layout: cp p1 p2 p3 n11 n12 ... n21 n22 ... n31 n32 ...
def group_snps(geno, ns):
    n = len(geno[0])
    m = len(geno)
    np = n - sum(ns)
    nf = len(ns)

    phap = []
    for i in range(np):
        phap.append("".join([e[i] for e in geno]))

    p1hap = phap[0]
    ldb = [phap.index(e) for e in phap]

    rec = 0
    fams = index_family_rils(ns)
    for j in range(nf):
        p2hap = phap[j+1]
        for i in fams[j]:
            hap = "".join([e[i+np] for e in geno])
            code = -1
            if hap == p1hap:
                code = ldb[0]
            elif hap == p2hap:
                code = ldb[j+1]
            else:
                rec += 1
                s1 = [p1hap[k] == hap[k] for k in range(m)].count(True)
                s2 = [p2hap[k] == hap[k] for k in range(m)].count(True)
                code = ldb[0] if s1 >= s2 else ldb[j+1]
            ldb.append(code)

    print m, rec
    return ldb

def SNPLDB_NAM(blk_file, snp_file, ns):
    blk_chr, blk_start, blk_stop = read_blocks(blk_file)
    ind, loc, chr, pos, geno = read_snp_genotype(snp_file)
    for i in range(len(blk_chr)):
        snps = index_block_snps(blk_chr[i], blk_start[i], blk_stop[i], chr, pos)
        if not snps:
            continue
        g = group_snps([geno[j] for j in snps], ns)
        for j in snps:
            loc[j] = None
        loc[snps[0]] = "BLOCK_{0}_{1}_{2}".format(blk_chr[i], blk_start[i], blk_stop[i])
        chr[snps[0]] = blk_chr[i]
        pos[snps[0]] = blk_start[i]
        geno[snps[0]] = [str(e) for e in g]
    with open("NAM.snpldb.geno","w") as f:
        f.write("Locus\tChromosome\tPosition\tDistance\t")
        f.write("\t".join(ind))
        f.write("\n")
        for i in range(len(loc)):
            if not loc[i]:
                continue
            f.write("{0}\t{1}\t{2}\t0\t".format(loc[i], chr[i], pos[i]))
            f.write("\t".join(geno[i]))
            f.write("\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "SNPLDB_NAM.py blk_file snp_file n1,n2,n3,..."
        print "blk_file format: https://github.com/njau-sri/rtm-gwas/wiki/File-format"
        print "snp_file format: https://github.com/njau-sri/gconv"
        print "n1,n2,n3: the number of lines in each family"
        print "the first individual is the common parent, following with other parents"
    else:
        ns = [int(e) for e in sys.argv[3].split(",")]
        SNPLDB_NAM(sys.argv[1], sys.argv[2], ns)
