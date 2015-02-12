#ifndef SORTEDDADA_H_
#define SORTEDDADA_H_

#include "DadaHeader.h"
#include "DadaReorder.h"
#include <complex>
#include <string>
#include <fstream>
#include <vector>

namespace dada {

class SortedDada
{
public:
    SortedDada(const char *dadaFilename);
    const DadaHeader header;
    int inputSize() const {return mOrder.inputSize();};
    int outputSize() const {return mOrder.outputSize();};
    void setLineMapping(int corrInput, int cable) {return mOrder.setLineMapping(corrInput, cable);};
    void setLineMapping(const char *corrInput, const char *cable) {return mOrder.setLineMapping(corrInput, cable);};
    int setLineMappingFromFile(const char *filename) {return mOrder.setLineMappingFromFile(filename);};
    void applyGains(const std::vector<std::complex<float> > &gains, const std::vector<char> &gainFlags);
    void resetGains() {mOrder.resetGains();};
    std::vector<float> &rRawChunk(int index);
    std::vector<std::complex<float> > &rGetChunk(int index);
    std::vector<std::complex<float> > &rNextChunk();
    std::vector<char> &rCurrentVisFlags() {return mOutVisFlags;};
    int prevChunkIndex() const {return mPrevChunk;};
    const char *filename() const {return mFileName.c_str();};
    void rewind() {mPrevChunk=-1;};
private:
    std::string mFileName;
    std::ifstream mDadaFile;
    int mPrevChunk;
    DadaReorder mOrder;
    int mInChunkBytes;
    std::vector<float> mRawData;
    std::vector<std::complex<float> > mSortedData;
    std::vector<std::complex<float> > mGains;
    // std::vector<bool> is bit packed so we'll use chars
    std::vector<char> mGainFlags;
    std::vector<char> mOutVisFlags;
};

} // namespace dada

#endif // SORTEDDADA_H_
