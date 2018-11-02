#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "cmdline.h"
#include "vcf.h"
#include "split.h"
#include "util.h"


using std::size_t;


namespace {


struct Parameter
{
    std::string vcf;
    std::string blk;
    std::string out;
    std::string nam;
    double maf = 0.01;
    double inform = 0.95;
    int maxlen = 200000;
    int llim = 70;
    int ulim = 98;
    int recomb = 90;
    int batch = 20000;
    bool ril = false;
} par ;


struct Block
{
    int first;
    int last;
    int length;
    float inform;

    Block(int fi, int la, int le, float in)
        : first(fi), last(la), length(le),inform(in)
    {}
};


// Gabriel, S.B. et al. The structure of haplotype blocks in the human genome. Science, 2002, 296(5576): 2225-9.
// Barrett, J.C. et al. Haploview: analysis and visualization of LD and haplotype maps. Bioinformatics, 2005, 21: 263-5.

// Estimates haplotype frequencies via the EM algorithm
// AB, Ab, aB, ab, AaBb
void calc_hap_prob_EM(int n11, int n12, int n21, int n22, int ndh, double &p11, double &p12, double &p21, double &p22)
{
    static const int maxit = 1000;
    static const double tol = 1e-10;

    double n = n11 + n12 + n21 + n22 + ndh * 2;
    p11 = n11 / n;
    p12 = n12 / n;
    p21 = n21 / n;
    p22 = n22 / n;

    if (ndh == 0)
        return;

    auto cp11 = p11;
    auto cp12 = p12;
    auto cp21 = p21;
    auto cp22 = p22;

    auto h = ndh / n;
    auto x = h / 2;
    auto y = h - x;

    for (int i = 0; i < maxit; ++i) {
        p11 = cp11 + x;
        p12 = cp12 + y;
        p21 = cp21 + y;
        p22 = cp22 + x;
        auto z = h * p11 * p22 / (p11 * p22 + p12 * p21);
        if (std::fabs(x - z) < tol)
            break;
        x = z;
        y = h - x;
    }
}

void count_freq_3x3(const std::vector<char> &x, const std::vector<char> &y, int freq[3][3])
{
    auto n = x.size();

    freq[0][0] = freq[0][1] = freq[0][2] = 0;
    freq[1][0] = freq[1][1] = freq[1][2] = 0;
    freq[2][0] = freq[2][1] = freq[2][2] = 0;

    for (size_t i = 0; i < n; ++i) {
        if (x[i] >= 0 && y[i] >= 0) {
            auto xi = static_cast<size_t>(x[i]);
            auto yi = static_cast<size_t>(y[i]);
            ++freq[xi][yi];
        }
    }
}

// D' 95% confidence interval estimate
int calc_dprime_CI(int freq[3][3], int &lower, int &upper)
{
    int n11 = freq[0][0] * 2 + freq[0][1] + freq[1][0];
    int n12 = freq[0][2] * 2 + freq[0][1] + freq[1][2];
    int n21 = freq[2][0] * 2 + freq[1][0] + freq[2][1];
    int n22 = freq[2][2] * 2 + freq[2][1] + freq[1][2];
    int ndh = freq[1][1];

    double nn = n11 + n12 + n21 + n22 + ndh * 2;
    if (nn < 4)
        return 1;

    double p11 = 0.0;
    double p12 = 0.0;
    double p21 = 0.0;
    double p22 = 0.0;

    calc_hap_prob_EM(n11, n12, n21, n22, ndh, p11, p12, p21, p22);

    if (n11 > n22) {
        std::swap(n11, n22);
        std::swap(n12, n21);
        std::swap(p11, p22);
        std::swap(p12, p21);
    }

    if (n11 > n12 || n11 > n21) {
        if (n12 < n21) {
            std::swap(n11, n12);
            std::swap(n21, n22);
            std::swap(p11, p12);
            std::swap(p21, p22);
        }
        else {
            std::swap(n11, n21);
            std::swap(n12, n22);
            std::swap(p11, p21);
            std::swap(p12, p22);
        }
    }

    auto p1x = (n11 + n12 + ndh) / nn;
    auto p2x = 1 - p1x;
    auto px1 = (n11 + n21 + ndh) / nn;
    auto px2 = 1 - px1;

    auto D = p11 - p1x*px1;
    auto Dmax = D < 0.0 ? std::min(p1x*px1, p2x*px2) : std::min(p1x*px2, p2x*px1);

    if (p11 < 1e-10) p11 = 1e-10;
    if (p12 < 1e-10) p12 = 1e-10;
    if (p21 < 1e-10) p21 = 1e-10;
    if (p22 < 1e-10) p22 = 1e-10;
    auto LL1 = n11*log(p11) + n12*log(p12) + n21*log(p21) + n22*log(p22) + ndh*log(p11*p22 + p12*p21);

    if (D < 0.0) {
        std::swap(p1x, p2x);
        std::swap(n11, n21);
        std::swap(n12, n22);
    }

    double tp = 0.0;
    double ls[101];
    auto Dstep = Dmax / 100;

    for (int i = 0; i <= 100; ++i) {
        auto q11 = i*Dstep + p1x*px1;
        auto q12 = p1x - q11;
        auto q21 = px1 - q11;
        auto q22 = p2x - q21;
        if (i == 100) {
            if (q11 < 1e-10) q11 = 1e-10;
            if (q12 < 1e-10) q12 = 1e-10;
            if (q21 < 1e-10) q21 = 1e-10;
            if (q22 < 1e-10) q22 = 1e-10;
        }
        auto LL2 = n11*log(q11) + n12*log(q12) + n21*log(q21) + n22*log(q22) + ndh*log(q11*q22 + q12*q21);
        auto prob = std::exp(LL2 - LL1);
        ls[i] = prob;
        tp += prob;
    }

    double sp = 0.0;
    auto tp5 = tp * 0.05;
    for (int i = 0; i <= 100; ++i) {
        sp += ls[i];
        if (sp > tp5 && sp - ls[i] < tp5) {
            lower = i - 1;
            break;
        }
    }

    sp = 0.0;
    for (int i = 100; i >= 0; --i) {
        sp += ls[i];
        if (sp > tp5 && sp - ls[i] < tp5) {
            upper = i + 1;
            break;
        }
    }

    return 0;
}

void classify_dprime_CI(int llim, int ulim, bool &strong, bool &recomb, bool &strong2, bool &strong34)
{
    if (ulim >= par.ulim) {
        if (llim > par.llim)
            strong = true;
        if (llim > 80)
            strong2 = strong34 = true;
        else if (llim > 50)
            strong34 = true;
    }

    if (ulim < par.recomb)
        recomb = true;
}

double test_block_Gabriel(size_t x, size_t y, size_t n,
                          const std::vector< std::vector<char> > &sr,
                          const std::vector< std::vector<char> > &ss)
{
    int s = 0, r = 0;
    for (size_t i = x; i <= y; ++i) {
        for (auto j = i + 1; j <= y; ++j) {
            r += sr[j][i];
            if (n > 4)
                s += sr[i][j];
            else if (n == 2)
                s += ss[i][j];
            else if (n == 3 || n == 4)
                s += ss[j][i];
        }
    }

    int t = s + r;
    if (n == 2) {
        if (t < 1)
            return -1;
    }
    else if (n == 3) {
        if (t < 3)
            return -2;
    }
    else {
        if (t < 6)
            return -3;
    }

    static const double eps = std::numeric_limits<double>::epsilon();
    double inform =  static_cast<double>(s) / t;
    if (inform - par.inform > eps)
        return inform;

    return -4;
}

int find_block_Gabriel(const Genotype &gt, const std::vector<size_t> &snps, std::vector< std::pair<int,int> > &ppos)
{
    auto n = gt.ind.size();

    std::vector<int> pos;
    std::vector< std::vector<char> > geno;

    for (auto j : snps) {
        pos.push_back(gt.pos[j]);

        std::vector<char> jgeno(n, -1);  // AA,Aa,aa -> 0,1,2

        if (gt.ploidy == 1) {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i];
                if (a == 1 || a == 2)
                    jgeno[i] = (a - 1) * 2;
            }
        }
        else {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][2*i];
                auto b = gt.dat[j][2*i+1];
                if ((a == 1 || a == 2) && (b == 1 || b == 2))
                    jgeno[i] = a + b - 2;
            }
        }

        geno.push_back(jgeno);
    }

    auto w = static_cast<size_t>(par.batch);
    auto m = pos.size();

    for (size_t i = 0; i < m; ++i) {
        auto pos1 = pos[i];
        for (size_t j = i + 1; j < m; ++j) {
            auto pos2 = pos[j];
            if (pos2 <= pos1) {
                std::cerr << "ERROR: chromosome positions must be in ascending order: " << pos1 << " " << pos2 << "\n";
                return 1;
            }
            auto dist = pos2 - pos1;
            if (dist <= par.maxlen && w < (j - i + 1))
                w = j - i + 1;
        }
    }

    if (w > m)
        w = m;

    int freq[3][3];

    std::vector< std::vector<char> > ci(w, std::vector<char>(w, -1));
    std::vector< std::vector<char> > sr(w, std::vector<char>(w, 0));
    std::vector< std::vector<char> > ss(w, std::vector<char>(w, 0));

    int llim = -1, ulim = -1;
    bool strong = false, recomb = false, strong2 = false, strong34 = false;

    std::vector<Block> blks;

    size_t start = 0;
    std::cerr << "INFO: " << start + 1 << " - " << start + w << "\n";

    for (;;) {
        for (size_t i = 0; i < w; ++i) {
            auto x = start + i;
            auto pos1 = pos[x];
            for (auto j = i + 1; j < w; ++j) {
                auto y = start + j;
                auto pos2 = pos[y];
                if (pos2 - pos1 > par.maxlen)
                    break;

                count_freq_3x3(geno[x], geno[y], freq);

                llim = ulim = -1;
                calc_dprime_CI(freq, llim, ulim);

                strong = recomb = strong2 = strong34 = false;
                classify_dprime_CI(llim, ulim, strong, recomb, strong2, strong34);

                ci[j][i] = llim;
                ci[i][j] = ulim;
                sr[i][j] = strong;
                sr[j][i] = recomb;
                ss[i][j] = strong2;
                ss[j][i] = strong34;
            }
        }

        for (size_t i = 0; i < w; ++i) {
            auto pos1 = pos[start + i];
            for (auto j = i + 1; j < w; ++j) {
                auto pos2 = pos[start + j];
                auto dist = pos2 - pos1;
                if (dist > par.maxlen)
                    break;
                if (ci[j][i] < par.llim || ci[i][j] < par.ulim)
                    continue;
                auto q = j - i + 1;
                if ((q == 2 && dist > 20000) || (q == 3 && dist > 30000))
                    continue;
                auto inform = test_block_Gabriel(i, j, q, sr, ss);
                if (inform > 0)
                    blks.emplace_back(start + i, start + j, dist, inform);
            }
        }

        if (start + w >= m)
            break;

        auto pos2 = pos[start + w];
        for (size_t i = 1; i < w; ++i) {
            auto pos1 = pos[++start];
            auto dist = pos2 - pos1;
            if (dist <= par.maxlen)
                break;
        }

        if (start + w >= m)
            w = m - start;

        std::cerr << "INFO: " << start + 1 << " - " << start + w << "\n";
    }

    if (blks.empty())
        return 0;

    auto cmp = [&](const Block &a, const Block &b) {
        if (a.length > b.length) return true;
        if (a.length < b.length) return false;
        if ((b.first > a.first && b.first < a.last) || (b.last > a.first && b.last < a.last)) {
            if (a.inform > b.inform) return true;
            if (a.inform < b.inform) return false;
            int a1 = -1, a2 = -1;
            count_freq_3x3(geno[a.first], geno[a.last], freq);
            calc_dprime_CI(freq, a1, a2);
            int b1 = -1, b2 = -1;
            count_freq_3x3(geno[b.first], geno[b.last], freq);
            calc_dprime_CI(freq, b1, b2);
            if (a1 > b1) return true;
            if (a1 < b1) return false;
        }
        return a.first < b.first;
    };

    std::sort(blks.begin(), blks.end(), cmp);

    std::vector<char> inblock(m+1,0);

    for (auto &e : blks) {
        auto first = e.first, last = e.last;
        if (inblock[first] || inblock[last])
            continue;
        ppos.emplace_back(pos[first], pos[last]);
        for (auto i = first; i <= last; ++i)
            inblock[i] = 1;
    }

    std::sort(ppos.begin(), ppos.end());

    return 0;
}

int read_block(const std::string &filename, std::vector<std::string> &chr,
               std::vector<int> &pos1, std::vector<int> &pos2)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file: " << filename << "\n";
        return 1;
    }

    size_t ln = 0;
    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if ( vs.empty() )
            continue;

        if (vs.size() < 3) {
            std::cerr << "ERROR: expected at least 3 columns for each line of block file\n";
            return 1;
        }

        auto start = std::stoi(vs[1]);
        auto stop = std::stoi(vs[2]);

        if (stop <= start) {
            std::cerr << "ERROR: invalid block: " << vs[0] << " " << vs[1] << " " << vs[2] << "\n";
            return 1;
        }

        chr.push_back(vs[0]);
        pos1.push_back(start);
        pos2.push_back(stop);
    }

    return 0;
}

std::vector< std::vector<size_t> > index_snps(const Genotype &gt, const std::vector<std::string> &chr)
{
    auto nchr = chr.size();

    std::map<std::string,size_t> mapchr;
    for (size_t i = 0; i < nchr; ++i)
        mapchr[chr[i]] = i;

    auto m = gt.loc.size();

    std::vector< std::vector<size_t> > idx(nchr);
    for (size_t i = 0; i < m; ++i) {
        auto j = mapchr[gt.chr[i]];
        idx[j].push_back(i);
    }

    return idx;
}

size_t count_match(const std::vector<allele_t> &x, const std::vector<allele_t> &y)
{
    size_t c = 0;

    auto n = x.size();
    for (size_t i = 0; i < n; ++i) {
        if (x[i] == y[i])
            ++c;
    }

    return c;
}

void group_snps(const Genotype &gt, std::vector<size_t> &snps, Genotype &ogt)
{
    auto n = gt.ind.size();
    std::vector< std::vector<allele_t> > dat;

    for (size_t i = 0; i < n; ++i) {
        std::vector<allele_t> v1, v2;
        for (auto j : snps) {
            if (gt.ploidy == 1) {
                v1.push_back(gt.dat[j][i]);
            }
            else {
                v1.push_back(gt.dat[j][i*2]);
                v2.push_back(gt.dat[j][i*2+1]);
            }
        }
        dat.push_back(v1);
        if ( ! v2.empty() )
            dat.push_back(v2);
    }

    auto haps = unique(dat);
    auto nhap = haps.size();

    std::vector<size_t> freq;
    for (auto &e : haps)
        freq.push_back( count(dat,e) );

    auto ord = order(freq);
    std::reverse(ord.begin(), ord.end());
    subset(haps,ord).swap(haps);
    subset(freq,ord).swap(freq);

    auto maf = static_cast<size_t>(std::ceil(par.maf * n * gt.ploidy));

    auto na = count_if(freq, [maf](size_t a) { return a >= maf; });

    // TODO: use cluster analysis if na = 1
    // TODO: maf_new may less than maf_threshold

    std::vector<size_t> hc(nhap);
    std::iota(hc.begin(), hc.end(), 0);
    for (auto i = na; i < nhap; ++i) {
        std::vector<size_t> v;
        for (size_t j = 0; j < na; ++j) {
            auto s = count_match(haps[i], haps[j]);
            v.push_back(s);
        }
        hc[i] = index(v, * std::max_element(v.begin(), v.end()));
    }

    auto start = gt.pos[snps[0]];
    auto stop = start;
    for (auto j : snps) {
        auto pos = gt.pos[j];
        if (pos < start)
            start = pos;
        else if (pos > stop)
            stop = pos;
    }

    auto loc = "LDB_" + gt.chr[snps[0]] + "_" + std::to_string(start) + "_" + std::to_string(stop);
    ogt.loc.push_back(loc);
    ogt.chr.push_back(gt.chr[snps[0]]);
    ogt.pos.push_back(start);

    std::vector<allele_t> v;

    if (gt.ploidy == 1) {
        for (size_t i = 0; i < n; ++i) {
            auto wh = index(haps, dat[i]);
            auto code = static_cast<allele_t>(hc[wh]+1);
            v.push_back(code);
        }
    }
    else {
        for (size_t i = 0; i < n; ++i) {
            auto wh = index(haps, dat[i*2]);
            auto code = static_cast<allele_t>(hc[wh]+1);
            v.push_back(code);
            wh = index(haps, dat[i*2+1]);
            code = static_cast<allele_t>(hc[wh]+1);
            v.push_back(code);
        }
    }

    if (na == 0)
        std::fill(v.begin(), v.end(), 0);

    ogt.dat.push_back(v);

    std::vector<std::string> allele;
    for (size_t i = 0; i < na; ++i) {
        std::string si;
        auto p = snps.size();
        for (size_t j = 0; j < p; ++j) {
            auto jj = snps[j];
            auto a = haps[i][j];
            if (a)
                si.append(gt.allele[jj][a-1]);
            else
                si.push_back('N');
        }
        allele.push_back(si);
    }

    ogt.allele.push_back(allele);
}

void recode_save_allele(Genotype &gt)
{
    std::ofstream ofs(par.out + ".allele");

    if ( ofs )
        ofs << "Locus\tCode\tAllele\n";
    else
        std::cerr << "ERROR: can't open file for writing: " << par.out << ".allele\n";

    auto m = gt.loc.size();
    for (size_t i = 0; i < m; ++i) {
        auto n = gt.allele[i].size();
        for (size_t j = 0; j < n; ++j) {
            if ( ofs )
                ofs << gt.loc[i] << "\t" << j << "\t" << gt.allele[i][j] << "\n";
            gt.allele[i][j] = std::to_string(j);
        }
    }
}

// RIL population
// layout: p1 p2 ind1 ind2 ...
int rtm_gwas_snpldb_ril()
{
    std::cerr << "ERROR: RIL function is not implemented\n";

    return 1;
}

// NAM population
// layout: cp p1 p2 p3 f1-ind1 f1-ind2 ... f2-ind1 f2-ind2 ... f3-ind1 f3-ind2 ...
int rtm_gwas_snpldb_nam()
{
    std::vector<std::string> vs;
    split(par.nam, ",", vs);

    std::cerr << "INFO: " << vs.size() << " families in NAM population\n";

    std::vector<size_t> vn;
    for (auto &e : vs) {
        auto a = std::stoi(e);
        if (a < 1) {
            std::cerr << "ERROR: invalid family size of NAM population: " << e << "\n";
            return 1;
        }
        vn.push_back( static_cast<size_t>(a) );
    }

    std::cerr << "ERROR: NAM function is not implemented\n";

    return 1;
}


} // namespace


int rtm_gwas_snpldb(int argc, char *argv[])
{
    std::cerr << "RTM-GWAS v1.5 SNPLDB (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd;

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--blk", "predefined block file", "");
    cmd.add("--out", "output file", "rtm-gwas-snpldb.out");
    cmd.add("--maf", "minimum minor haplotype frequency", "0.01");
    cmd.add("--maxlen", "maximum length of blocks", "200000");

    cmd.add("--llim", "lower limit CI for strong LD", "70");
    cmd.add("--ulim", "upper limit CI for string LD", "98");
    cmd.add("--recomb", "upper limit CI for strong recombination", "90");
    cmd.add("--inform", "minimum fraction of informative strong LD", "0.95");

    cmd.add("--batch", "number of SNPs in a batch", "10000");

    cmd.add("--nam", "nested association mapping population, n1 n2 n3 ...", "");
    cmd.add("--ril", "recombinant inbred line population");

    cmd.parse(argc, argv);

    if (argc < 2) {
        cmd.show();
        return 1;
    }

    par.vcf = cmd.get("--vcf");
    par.blk = cmd.get("--blk");
    par.out = cmd.get("--out");

    par.maf = std::stod(cmd.get("--maf"));
    par.maxlen = std::stoi(cmd.get("--maxlen"));

    par.llim = std::stoi(cmd.get("--llim"));
    par.ulim = std::stoi(cmd.get("--ulim"));
    par.recomb = std::stoi(cmd.get("--recomb"));
    par.inform = std::stod(cmd.get("--inform"));

    par.batch = std::stoi(cmd.get("--batch"));
    if (par.batch < 0)
        par.batch = 10000;

    par.nam = cmd.get("--nam");
    par.ril = cmd.has("--ril");

    if ( par.ril )
        return rtm_gwas_snpldb_ril();

    if ( ! par.nam.empty() )
        return  rtm_gwas_snpldb_nam();

    Genotype gt;

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::vector<std::string> blk_chr;
    std::vector<int> blk_pos1, blk_pos2;

    if ( ! par.blk.empty() )
        read_block(par.blk, blk_chr, blk_pos1, blk_pos2);

    auto chrid = stable_unique(gt.chr);
    auto nchr = chrid.size();
    auto snps = index_snps(gt, chrid);

    if ( blk_chr.empty() ) {
        for (size_t i = 0; i < nchr; ++i) {
            std::cerr << "INFO: finding blocks on chromosome " << chrid[i] << "\n";
            std::vector< std::pair<int,int> > ppos;
            if (find_block_Gabriel(gt, snps[i], ppos) != 0)
                return 1;
            blk_chr.insert(blk_chr.end(), ppos.size(), chrid[i]);
            for (auto &e : ppos) {
                blk_pos1.push_back(e.first);
                blk_pos2.push_back(e.second);
            }
        }
    }

    auto m = gt.loc.size();
    auto nb = blk_chr.size();

    Genotype ldb;
    std::vector<char> inblock(m,0);

    std::vector<int> blk_length(nb,0);
    std::vector<size_t> blk_size(nb,0);

    for (size_t i = 0; i < nchr; ++i) {
        Genotype gtchr;
        for (size_t k = 0; k < nb; ++k) {
            if (blk_chr[k] != chrid[i])
                continue;
            std::vector<size_t> idx;
            for (auto j : snps[i]) {
                if (gt.pos[j] < blk_pos1[k] || gt.pos[j] > blk_pos2[k])
                    continue;
                inblock[j] = true;
                idx.push_back(j);
            }
            blk_length[k] = blk_pos2[k] - blk_pos1[k];
            blk_size[k] = idx.size();
            if ( idx.empty() ) {
                std::cerr << "WARNING: no SNPs were found in block: " << chrid[i] << " "
                          << blk_pos1[k] << " " << blk_pos2[k] << "\n";
                continue;
            }
            group_snps(gt, idx, gtchr);
        }

        for (auto j : snps[i]) {
            if ( inblock[j] ) {
                auto k = index(gtchr.pos, gt.pos[j]);
                if (k != gtchr.pos.size()) {
                    ldb.loc.push_back(gtchr.loc[k]);
                    ldb.chr.push_back(gtchr.chr[k]);
                    ldb.pos.push_back(gtchr.pos[k]);
                    ldb.dat.push_back(gtchr.dat[k]);
                    ldb.allele.push_back(gtchr.allele[k]);
                }
            }
            else {
                ldb.loc.push_back("LDB_" + gt.loc[j]);
                ldb.chr.push_back(gt.chr[j]);
                ldb.pos.push_back(gt.pos[j]);
                ldb.dat.push_back(gt.dat[j]);
                ldb.allele.push_back(gt.allele[j]);
            }
        }
    }

    std::ofstream ofs(par.out + ".block");

    if ( ! ofs )
        std::cerr << "ERROR: can't open file for writing: " << par.out << ".block\n";
    else {
        ofs << "Chromosome\tStart\tStop\tLength\tSNPs\n";
        for (size_t i = 0; i < nb; ++i)
            ofs << blk_chr[i] << "\t" << blk_pos1[i] << "\t" << blk_pos2[i] << "\t"
                << blk_length[i] << "\t" << blk_size[i] << "\n";
    }

    ldb.ind = gt.ind;
    ldb.ploidy = gt.ploidy;

    recode_save_allele(ldb);

    if (write_vcf(ldb, par.out + ".vcf") != 0)
        return 3;

    return 0;
}
