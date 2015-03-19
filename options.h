/*
 * options.h
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>
#include <vector>

namespace dada2ms {

struct options {
	double latitude;
	double longitude;
	double altitude;

	bool append;       // Append to an existing MS
	bool firstOnly;    // Take only first integration
	bool autosOnly;    // Take only auto-correlations
	bool azel;         // MS should AZ-EL for coordinates (default is J2000)
	bool addWtSpec;    // Write a WEIGHT_SPECTRUM column
	bool addSPW;       // Add a new SPW/DATA_DESC
	bool applyCal;     // Apply existing calibrations during conversion
    bool applyTTCalBandpass; // Apply existing TTCal bandpass calibration
    bool applyTTCalPolcal;   // Apply existing TTCal polcal calibration
	bool antsAreITRF;  // Antenna positions are ITRF (default is relative to array position)

	int dataDescID;
	int startScan;

	std::string configFile;
	std::string remapFile;
	std::string calTable;
	std::string bcalTable; // TTCal bandpass calibration
	std::string jcalTable; // TTCal polcal calibration
	std::string antFile;
	std::string msName;

	std::vector<int> integrations;
	std::vector<std::string> dadaFile;

	options(int argc, char *argv[]);
};

template <typename T>
T StringToNumber (const std::string &Text);

template <typename T>
std::vector<T> split(const std::string &s, char delim);

} // namespace dada2ms

#endif /* OPTIONS_H_ */
