#include "BCalTable.h"

#include <fstream>
#include <stdexcept>

using namespace std;

BCalTable::BCalTable(char const* filename) {
    ifstream file(filename,ios::in|ios::binary);
    if (file.is_open()) {
        // Read in the calibration type
        char T;
        file.read(&T,1);

        // Throw if this is not the correct type of calibration table
        if (T != 'B') throw invalid_argument("Supplied calibration table is not a bandpass calibration.");

        // Read in the number of antennas and channels
        file.read(reinterpret_cast<char*>(&_Nant), sizeof(int));
        file.read(reinterpret_cast<char*>(&_Nchan),sizeof(int));

        // Read in the flags
        bool bool_flags[_Nant*_Nchan*2];
        _flags = vector<char>(_Nant*_Nchan*2);
        file.read(reinterpret_cast<char*>(bool_flags), _Nant*_Nchan*2*sizeof(bool));
        for (int i = 0; i < _Nant*_Nchan*2; ++i) {
            _flags[i] = static_cast<char>(bool_flags[i]);
        }

        // Read in the gains
        float float_gains[_Nant*_Nchan*2*2];
        _gains = vector<complex<float> >(_Nant*_Nchan*2);
        file.read(reinterpret_cast<char*>(float_gains), _Nant*_Nchan*2*2*sizeof(float));
        for (int i = 0; i < _Nant*_Nchan*2; ++i) {
            _gains[i] = complex<float>(float_gains[2*i],float_gains[2*i+1]);
        }
    }
    else {
        throw invalid_argument("Could not open calibration table.");
    }
}

