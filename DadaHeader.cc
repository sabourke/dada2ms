#include "DadaHeader.h"
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <ctime>
#include <iostream>

namespace dada {

DadaHeader::DadaHeader(const char *dadaFilename)
{
    loadDadaHeader(dadaFilename);
    mNAvg = atoi(mHeaderMap.at("NAVG").c_str());
    mTSamp = atof(mHeaderMap.at("TSAMP").c_str());
    mIntTime = mNAvg * mTSamp / 1e6;
    mNAnt = atoi(mHeaderMap["NSTATION"].c_str());
    mNBaseline = (mNAnt + 1) * mNAnt / 2;
    mNFreq = atoi(mHeaderMap["NCHAN"].c_str());
    mNPol = atoi(mHeaderMap["NPOL"].c_str());
    mNCorr = mNPol * mNPol;
    mCFreq = atof(mHeaderMap["CFREQ"].c_str()) * 1e6;
    mBandwidth = atof(mHeaderMap["BW"].c_str()) * 1e6;
    mChannelBandwidth = mBandwidth / mNFreq;
    mHeaderSize = atoi(mHeaderMap["HDR_SIZE"].c_str());
    mFileSize = atol(mHeaderMap["FILE_SIZE"].c_str());
    mBytesPerAvg = atol(mHeaderMap["BYTES_PER_AVG"].c_str());
    mNTime = mFileSize / mBytesPerAvg;
    mObsOffset = atol(mHeaderMap["OBS_OFFSET"].c_str()); // OBS_OFFSET is in bytes
    mObsOffsetSec = static_cast<double>(mObsOffset) / mBytesPerAvg * mIntTime;
    mStartTime = str2epoch(mHeaderMap["UTC_START"].c_str(), mObsOffsetSec);
    mFinishTime = mStartTime + mNTime * mIntTime;
}

void DadaHeader::loadDadaHeader(const char *dadaFilename)
{
    // Read DADA header into map
    std::ifstream inf(dadaFilename);
    if (inf.fail())
	    throw std::invalid_argument("Error opening dada file in DadaHeader::loadDadaHeader()");
    std::string key, value;
    while (inf >> key >> value) {
	if (mHeaderMap.count("HDR_SIZE") > 0 && inf.tellg() > atol(mHeaderMap["HDR_SIZE"].c_str()))
		break;
        mHeaderMap[key] = value;
    }
}

double DadaHeader::channelFreq(int channelNum)
{
	return mCFreq;
	return 0.0;
}

double DadaHeader::str2epoch(const char *time, double offsetSec)
{
    struct tm t;
    if (sscanf(time, "%d-%d-%d-%d:%d:%d", &t.tm_year, &t.tm_mon, &t.tm_mday, &t.tm_hour, &t.tm_min, &t.tm_sec) != 6)
        throw std::runtime_error("Invalid String enountered in str2MEpoch()");
    --t.tm_mon;
    t.tm_year -= 1900;
    t.tm_isdst = 0;

    time_t t0 = mktime(&t) - timezone;
    return static_cast<double>(t0) + offsetSec;
}

double DadaHeader::unix2mjd(double time)
{
	// Difference between 17.11.1858 and 1.1.1970 in seconds
	return time + 3506716800.0;
}


} // namespace dada

