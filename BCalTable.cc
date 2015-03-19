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
        _flags = new bool[_Nant*_Nchan*2];
        file.read(reinterpret_cast<char*>(_flags), _Nant*_Nchan*2*sizeof(bool));

        // Read in the gains
        float float_gains[_Nant*_Nchan*2*2];
        _gains = new complex<float>[_Nant*_Nchan*2];
        file.read(reinterpret_cast<char*>(float_gains), _Nant*_Nchan*2*2*sizeof(float));
        for (int i = 0; i < _Nant; ++i) {
            for (int j = 0; j < _Nchan; ++j) {
                for (int k = 0; k < 2; ++k) {
                    int idx = index(i,j,k);
                    _gains[idx] = complex<float>(float_gains[2*idx],float_gains[2*idx+1]);
                }
            }
        }
    }
    else {
        throw invalid_argument("Could not open calibration table.");
    }
}

BCalTable::~BCalTable() {
    delete[] _flags;
    delete[] _gains;
}

int BCalTable::index(int ant, int chan, int pol) {
    return ant*_Nchan*2 + chan*2 + pol;
}

