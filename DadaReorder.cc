#include "DadaReorder.h"
#include <sstream>
#include <stdexcept>
#include <cstdio>
#include <complex>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cstring>

namespace dada {

DadaReorder::DadaReorder(int nAnt, int nFreq, int nPol, int nCorr) :
    mApplyCal(false), mAutosOnly(false), mIndexIsValid(false),
    mNAnt(nAnt), mNFreq(nFreq), mNPol(nPol), mNCorr(nCorr),
    mNBaseline((nAnt + 1) * nAnt/2), mGpuBaselines(nAnt * (nAnt / 2 + 1)),
    mGpuHalfBlock(mGpuBaselines * nFreq * nCorr),
    mLineMap(new int[nAnt*nPol]),
    mBaselineIndex(nAnt * nPol, std::vector<int>(nAnt * nPol)),
    mConjBaseline(nAnt * nPol, std::vector<int>(nAnt * nPol, -1))
{
    // Initialize to nominal mapping
    for (int i=0; i<mNAnt*mNPol; ++i)
        mLineMap[i] = i;
}

DadaReorder::~DadaReorder()
{
    delete[] mLineMap;
}

void DadaReorder::setLineMapping(int corrInput, int cable)
{
    // Lines as integers indexed from 0.
    int nLines(mNAnt*mNPol);
    if (corrInput >= nLines || cable>= nLines) {
        std::stringstream message;
        message << "Line mapping " << corrInput << " -> " << cable << " is not valid (> Num Lines)";
        throw std::runtime_error(message.str());
    }
    mLineMap[corrInput] = cable;
    mIndexIsValid = false;
}

void DadaReorder::setLineMapping(const char *corrInput, const char *cable)
{
    // Lines in the format "1A", "256B" (indexed from 1).
    setLineMapping(simpleLineNum(corrInput), simpleLineNum(cable));
}

int DadaReorder::setLineMappingFromFile(const char *filename)
{
    int corrInput, cable, read(0);
    std::string ant1, ant2;
    std::ifstream inf(filename);
    while (inf.good()) {
        char c = inf.peek();
        while ((c == '#') || (c == '\n')) {
            // Ignore lines starting with # and empty lines
            std::string line;
            std::getline(inf, line);
            c = inf.peek();
        }
        inf >> ant1 >> ant2;
        corrInput = simpleLineNum(ant1.c_str());
        cable = simpleLineNum(ant2.c_str());
        setLineMapping(corrInput, cable);
        read++;
    }
    if (!inf.eof())
        throw std::runtime_error("Error reading linemap file");
    return read;
}

void DadaReorder::buildIndex()
{
    for (int i=0; i<mNAnt/2; i++) {
        for (int rx=0; rx<2; rx++) {
            for (int j=0; j<=i; j++) {
                for (int ry=0; ry<2; ry++) {
                    int ant1 = 2*i+rx;
                    int ant2 = 2*j+ry;
                    int l = (2*ry+rx) * mGpuBaselines / 4 + i*(i+1)/2 + j;
                    for (int pol1=0; pol1<mNPol; pol1++) {
                        for (int pol2=0; pol2<mNPol; pol2++) {
                            int reg_index = (l*mNPol+pol1)*mNPol+pol2;
                            int line1, line2;
                            line1 = mLineMap[2 * ant1 + pol1];
                            line2 = mLineMap[2 * ant2 + pol2];
                            if (line1 < line2) {
                                int tmp = line1;
                                line1 = line2;
                                line2 = tmp;
                                mConjBaseline[line2][line1] = 1;
                            }
                            if (line1 / 2 == line2 / 2)
                                mBaselineIndex[line1][line2] = reg_index;
                            mBaselineIndex[line2][line1] = reg_index;
                        }
                    }
                }
            }
        }
    }
}

void DadaReorder::applyGains(std::complex<float> *gains, char *gainFlags, char *outVisFlags)
{
    mApplyCal = true;
    mGains = gains;
    mGainFlags = gainFlags;
    mOutVisFlags = outVisFlags;
}

void DadaReorder::applyJones(std::complex<float> *jones, char *jonesFlags, char *outVisFlags)
{
    if (mNPol != 2)
        throw std::logic_error("must have fully polarized visibilities to apply a Jones matrix calibration");
    mApplyJones = true;
    mJones = jones;
    mJonesFlags = jonesFlags;
    mOutVisFlags = outVisFlags;
}

void DadaReorder::sortData(float *dadaArr, float *outArr)
{
	if (!mIndexIsValid)
        buildIndex();
    int freqs_x_pols = mNFreq * mNPol;
    int out_offset=0, baseline=-1;
    for (int ant1=0; ant1<mNAnt; ant1++) {
        for (int ant2=ant1; ant2<mNAnt; ant2++) {
            if (mAutosOnly && ant1 != ant2)
                continue;
            baseline++;
            for (int f=0; f<mNFreq; f++) {
                int dada_freq_offset = f * mGpuBaselines * mNCorr;
                for (int pol1=0; pol1<mNPol; pol1++) {
                    int line1 = 2*ant1 + pol1;
                    for (int pol2=0; pol2<mNPol; pol2++) {
                        int line2 = 2*ant2 + pol2;
                        int dada_index = dada_freq_offset + mBaselineIndex[line1][line2];
                        if (mApplyCal) {
                            float vr = dadaArr[dada_index];
                            float vi = mConjBaseline[line1][line2] * dadaArr[dada_index+mGpuHalfBlock];
                            int l0_index = ant1 * freqs_x_pols + f * mNPol + pol1;
                            int l1_index = ant2 * freqs_x_pols + f * mNPol + pol2;
                            std::complex<float> g0 = mGains[l0_index];
                            float g0r = g0.real();
                            float g0i = g0.imag();
                            std::complex<float> g1 = mGains[l1_index];
                            float g1r = g1.real();
                            float g1i = g1.imag();
                            // (G0)(V)(G1)* factored out:
                            outArr[out_offset] = g0r*vr*g1r - g0i*vi*g1r + g0i*vr*g1i + g0r*vi*g1i;
                            outArr[out_offset+1] = g0i*vr*g1r + g0r*vi*g1r - g0r*vr*g1i + g0i*vi*g1i;
                            int corr = pol1 * mNPol + pol2;
                            bool gainFlag0 = static_cast<bool>(mGainFlags[l0_index]);
                            bool gainFlag1 = static_cast<bool>(mGainFlags[l1_index]);
                            mOutVisFlags[baseline * mNFreq * mNCorr + f * mNCorr + corr] = static_cast<char>(gainFlag0 || gainFlag1);
                        } else {
                            outArr[out_offset] = dadaArr[dada_index];
                            outArr[out_offset+1] = mConjBaseline[line1][line2] * dadaArr[dada_index+mGpuHalfBlock];
                        }
                        out_offset += 2;
                    }
                }
                if (mApplyJones) {
                    // mApplyJones is only true if mNPol is 2.
                    // Therefore we can specialize on the case of 2x2 correlation products.

                    // Start be rewinding the offset counter so that we can take another pass
                    // over the same set of visibilities.
                    int offset = out_offset - 2*mNPol*mNPol;
                    int j0_offset = ant1*mNFreq*mNPol*mNPol + f*mNPol*mNPol;
                    int j1_offset = ant2*mNFreq*mNPol*mNPol + f*mNPol*mNPol;

                    // Now compute the matrix product (J0)(V)(J1)*
                    std::complex<float> vxx = std::complex<float>(outArr[offset+0],outArr[offset+1]);
                    std::complex<float> vxy = std::complex<float>(outArr[offset+2],outArr[offset+3]);
                    std::complex<float> vyx = std::complex<float>(outArr[offset+4],outArr[offset+5]);
                    std::complex<float> vyy = std::complex<float>(outArr[offset+6],outArr[offset+7]);

                    std::complex<float> a0 = mJones[j0_offset+0];
                    std::complex<float> b0 = mJones[j0_offset+1];
                    std::complex<float> c0 = mJones[j0_offset+2];
                    std::complex<float> d0 = mJones[j0_offset+3];

                    std::complex<float> a1 = std::conj(mJones[j1_offset+0]);
                    std::complex<float> b1 = std::conj(mJones[j1_offset+1]);
                    std::complex<float> c1 = std::conj(mJones[j1_offset+2]);
                    std::complex<float> d1 = std::conj(mJones[j1_offset+3]);

                    std::complex<float> vxx_ = a0*vxx*a1 + b0*vyx*a1 + a0*vxy*b1 + b0*vyy*b1;
                    std::complex<float> vxy_ = a0*vxx*c1 + b0*vyx*c1 + a0*vxy*d1 + b0*vyy*d1;
                    std::complex<float> vyx_ = c0*vxx*a1 + d0*vyx*a1 + c0*vxy*b1 + d0*vyy*b1;
                    std::complex<float> vyy_ = c0*vxx*c1 + d0*vyx*c1 + c0*vxy*d1 + d0*vyy*d1;

                    outArr[offset+0] = std::real(vxx_);
                    outArr[offset+1] = std::imag(vxx_);
                    outArr[offset+2] = std::real(vxy_);
                    outArr[offset+3] = std::imag(vxy_);
                    outArr[offset+4] = std::real(vyx_);
                    outArr[offset+5] = std::imag(vyx_);
                    outArr[offset+6] = std::real(vyy_);
                    outArr[offset+7] = std::imag(vyy_);

                    // Apply flags
                    if (static_cast<bool>(mJonesFlags[ant1*mNFreq + f]) || static_cast<bool>(mJonesFlags[ant2*mNFreq + f])) {
                        for (int i = 0; i < mNPol*mNPol; ++i) {
                            mOutVisFlags[baseline*mNFreq*mNCorr + f*mNCorr + i] = static_cast<char>(true);
                        }
                    }
                }
            }
        }
    }
}

int DadaReorder::simpleLineNum(const char *antName)
{
    // Take a 1-indexed string, eg "256B" and return the integer line, eg 511.

    int ant_num, line, read;
    char pol;

    // TODO: More robust string parsing
    read = sscanf(antName, "%d%c", &ant_num, &pol);
    if ((read != 2) || (ant_num > 256 || (pol != 'A' && pol != 'B'))) {
        std::stringstream message;
        message << "Invalid Ant string: " << antName << " encountered in DadaReorder::simpleLineNum";
        throw std::runtime_error(message.str());
    }

    line = (ant_num-1) * 2;
    if (pol == 'B')
        line++;

    return line;
}

} // namespace dada
