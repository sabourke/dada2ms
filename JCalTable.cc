#include "JCalTable.h"

#include <fstream>
#include <stdexcept>

using namespace std;

JCalTable::JCalTable(char const* filename) {
    ifstream file(filename,ios::in|ios::binary);
    if (file.is_open()) {
        // Read in the calibration type
        char T;
        file.read(&T,1);

        // Throw if this is not the correct type of calibration table
        if (T != 'J') throw invalid_argument("Supplied calibration table is not a polarized calibration.");

        // Read in the number of antennas and channels
        file.read(reinterpret_cast<char*>(&_Nant), sizeof(int));
        file.read(reinterpret_cast<char*>(&_Nchan),sizeof(int));

        // Read in the flags
        bool bool_flags[_Nant*_Nchan];
        _flags = vector<char>(_Nant*_Nchan);
        file.read(reinterpret_cast<char*>(bool_flags), _Nant*_Nchan*sizeof(bool));
        for (int i = 0; i < _Nant*_Nchan; ++i) {
            _flags[i] = static_cast<char>(bool_flags[i]);
        }

        // Read in the gains
        float float_gains[_Nant*_Nchan*2*2*2];
        _gains = vector<complex<float> >(_Nant*_Nchan*2*2);
        file.read(reinterpret_cast<char*>(float_gains), _Nant*_Nchan*2*2*2*sizeof(float));
        for (int i = 0; i < _Nant*_Nchan*2*2; ++i) {
            _gains[i] = complex<float>(float_gains[2*i],float_gains[2*i+1]);
        }
    }
    else {
        throw invalid_argument("Could not open calibration table.");
    }
}

int JCalTable::index(int ant, int chan, int pol1, int pol2) {
    return ant*_Nchan*2*2 + chan*2*2 + pol1*2 + pol2;
}

