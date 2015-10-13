/*
 * ms_funcs.h
 */

#ifndef MS_FUNCS_H_
#define MS_FUNCS_H_

#include <string>

// casacore headers
#include <casa/Arrays.h>
#include <tables/Tables.h>
#include <ms/MeasurementSets.h>

#define FIELD_NAME "Zenith"
// Maybe these should be configurable:
#define ANT_NAME "ANT"
#define STATION_NAME "LWA-OVRO"
#define CORRELATOR_NAME "LEDA512"
#define OBSERVER "Default"
#define PROJECT "Default"
#define TELESCOPE_NAME "OVRO_MMA" // Used by CASA as reference position
#define ANTENNA_TYPE "GROUND-BASED"
#define ANTENNA_MOUNT "X-Y"
#define ANTENNA_DISH_DIAMETER 2.0
#define DEFAULT_INT_TIME 8.33

casa::MEpoch str2MEpoch(const char *time, double offset);
casa::MDirection getZenith(const casa::MPosition &pos, const casa::MEpoch &epoch);
inline double radians(double degrees);
double seaLevel(double latitude);
casa::Matrix<double> readAnts(const char *filename, int nAnt);
casa::Matrix<double> zenithUVWs(casa::Matrix<double> antPos);
casa::Matrix<double> itrfAnts(casa::Matrix<double> antPos, int utmzone);
void utm2latlong(int utmzone, double northing, double easting, double* latitude, double* longitude);
int addSourceTab(casa::MeasurementSet &ms);
int fillAntTab(casa::MSAntenna &ant, int nAnt, casa::Matrix<double> itrfPos);
int fillFeedTab(casa::MSFeed &feed, int nAnt);
int fillFieldTab(casa::MSField &field, const casa::MDirection *dir);
int addField(casa::MSField &field, const casa::String &name, casa::MDirection *dir);
int fillObservationTab(casa::MSObservation &observation, double startTime, double finishTime);
int fillPointingTab(casa::MSPointing &pointing, int nAnt, double time, const casa::MDirection *dir);
int fillPolarizationTab(casa::MSPolarization &polarization);
int fillProcessorTab(casa::MSProcessor &processor);
int fillSpWindowTab(casa::MSSpectralWindow &spw, int nFreq, double cFreq, double bw);
int fillSourceTab(casa::MSSource &source, double startTime, double finishTime, const casa::MDirection *dir);
int updateSourceTab(casa::MSSource &source, double startTime, double finishTime);
int updateObservationTab(casa::MSObservation &observation, double startTime, double finishTime);
void boolArray2charVector(casa::Array<casa::Bool> &boolArr, std::vector<char> &charVec);
void charVector2boolArray(std::vector<char> &charVec, casa::Array<casa::Bool> &boolArr);
void readCalTable(const char *calName, std::vector<std::complex<float> > &gain, std::vector<char> &flag);

#endif /* MS_FUNCS_H_ */
