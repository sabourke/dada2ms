#ifndef DADAREORDER_H_
#define DADAREORDER_H_

#include <vector>
#include <complex>

namespace dada {

class DadaReorder
{
public:
    DadaReorder(int nAnt, int nFreq, int nPol, int nCorr);
    ~DadaReorder();
    int nAnt() const {return mNAnt;};
    int nFreq() const {return mNFreq;};
    int nBaseline() const {return mNBaseline;};
    int nPol() const {return mNPol;};
    int nCorr() const {return mNCorr;};
    int inputSize() const {return mGpuBaselines * mNFreq * mNCorr * 2;};
    int outputSize() const {return mNBaseline * mNFreq * mNCorr;};
    void setLineMapping(int corrInput, int cable);
    void setLineMapping(const char *corrInput, const char *cable);
    int setLineMappingFromFile(const char *filename);
    void applyGains(std::complex<float> *gains, char *gainFlags, char *outVisFlags);
    void resetGains() {mApplyCal = false;};
    void sortData(float *inArr, float *outArr);
    static int simpleLineNum(const char *antName);

private:
    bool mApplyCal, mAutosOnly, mIndexIsValid;
    const int mNAnt, mNFreq, mNPol, mNCorr, mNBaseline;
    const int mGpuBaselines, mGpuHalfBlock; // Larger than nBaselines due to alignment
    int * const mLineMap; // Defines physically remapped lines. Eg, line X could be connected to correlator input Y
    std::vector<std::vector<int> > mBaselineIndex;
    std::vector<std::vector<int> > mConjBaseline;
    std::complex<float> *mGains;      // size MUST be nAnt * nFreq * nPol
    // These are char instead of bool to be compatible with std::vector
    char *mGainFlags;                 // size MUST be nAnt * nFreq * nPol
    char *mOutVisFlags;               // size MUST be nBaseline * nFreq * nCorr
    void buildIndex();
};

} // namespace dada

#endif // DADAREORDER_H_
