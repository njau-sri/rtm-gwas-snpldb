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

# layout: p1 p2 ind1 ind2 ...
def group_snps(geno):
    n = len(geno[0])
    m = len(geno)

    phap = []
    phap.append("".join([e[0] for e in geno]))
    phap.append("".join([e[1] for e in geno]))

    ldb = [0, 1]
    if phap[0] == phap[1]:
        ldb[1] = 0

    rec = 0
    for i in range(2,n):
        hap = "".join([e[i] for e in geno])
        code = -1
        if hap == phap[0]:
            code = ldb[0]
        elif hap == phap[1]:
            code = ldb[1]
        else:
            rec += 1
            s1 = [phap[0][j]==hap[j] for j in range(m)].count(True)
            s2 = [phap[1][j]==hap[j] for j in range(m)].count(True)
            code = ldb[0] if s1 >= s2 else ldb[1]
        ldb.append( code )
    print m, rec
    return ldb

def SNPLDB_RIL(blk_file, snp_file):
    blk_chr, blk_start, blk_stop = read_blocks(blk_file)
    ind, loc, chr, pos, geno = read_snp_genotype(snp_file)
    for i in range(len(blk_chr)):
        snps = index_block_snps(blk_chr[i], blk_start[i], blk_stop[i], chr, pos)
        if not snps:
            continue
        g = group_snps([geno[j] for j in snps])
        for j in snps:
            loc[j] = None
        loc[snps[0]] = "BLOCK_{0}_{1}_{2}".format(blk_chr[i], blk_start[i], blk_stop[i])
        chr[snps[0]] = blk_chr[i]
        pos[snps[0]] = blk_start[i]
        geno[snps[0]] = [str(e) for e in g]
    with open("RIL.snpldb.geno","w") as f:
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
    if len(sys.argv) != 3:
        print "SNPLDB_RIL.py blk_file snp_file"
        print "blk_file format: https://github.com/njau-sri/rtm-gwas/wiki/File-format"
        print "snp_file format: https://github.com/njau-sri/gconv"
        print "the first two individuals must be the parents"
    else:
        SNPLDB_RIL(sys.argv[1], sys.argv[2])
