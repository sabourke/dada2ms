#include "SortedDada.h"
#include <stdexcept>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace dada {

SortedDada::SortedDada(const char *dadaFilename) :
    header(dadaFilename),
    mFileName(dadaFilename),
    mDadaFile(dadaFilename, std::ifstream::in | std::ifstream::binary),
    mPrevChunk(-1),
    mOrder(header.nAnt(), header.nFreq(), header.nPol(), header.nCorr()),
    mInChunkBytes(mOrder.inputSize() * sizeof(float)),
    mRawData(mOrder.inputSize()),
    mSortedData(mOrder.outputSize()),
    mOutVisFlags(outputSize(), static_cast<char>(false))
{
}

std::vector<float> &SortedDada::rRawChunk(int index)
{
    if (index < 0 || index >= header.nTime())
        throw std::out_of_range("SortedDada::rRawChunk() Invalid index");
    mDadaFile.seekg(header.headerSize() + index * mInChunkBytes, std::ios_base::beg);
    if (!mDadaFile.good())
        throw std::runtime_error("Seek Error in SortedDada::rRawChunk()");
    mDadaFile.read(reinterpret_cast<char*>(mRawData.data()), mInChunkBytes);
    if (!mDadaFile.good())
        throw std::runtime_error("Read Error in SortedDada::rRawChunk()");
    mPrevChunk = index;
    return mRawData;
}

std::vector<std::complex<float> > &SortedDada::rGetChunk(int index)
{
    rRawChunk(index);
    mOrder.sortData(mRawData.data(), reinterpret_cast<float*>(mSortedData.data()));
    return mSortedData;
}

std::vector<std::complex<float> > &SortedDada::rNextChunk()
{
    return rGetChunk(mPrevChunk+1);
}

void SortedDada::applyGains(const std::vector<std::complex<float> > &gains, const std::vector<char> &gainFlags)
{
    int num_gains = header.nAnt() * header.nFreq() * header.nPol();
    if (gains.size() != num_gains)
        throw std::length_error("gains vector size error in SortedDada::applyGains");
    if (gainFlags.size() != num_gains)
        throw std::length_error("gainFlags vector size error in SortedDada::applyGains");
    mGains = gains;
    // Inputs are assumed to be CASA style gains, so invert them
    for (int i=0; i<mGains.size(); ++i)
        mGains[i] = std::complex<float>(1) / mGains[i];
    mGainFlags = gainFlags;
    mOrder.applyGains(mGains.data(), mGainFlags.data(), mOutVisFlags.data());
}

} // namespace dada
