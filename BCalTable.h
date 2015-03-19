#ifndef BCALTABLE_H_
#define BCALTABLE_H_

#include <complex>

class BCalTable {
public:
    BCalTable(char const* filename);
    ~BCalTable();

    int Nant() {return _Nant;}
    int Nchan() {return _Nchan;}
    bool flag(int ant, int chan, int pol) {return _flags[index(ant,chan,pol)];}
    std::complex<float> gain(int ant, int chan, int pol) {return _gains[index(ant,chan,pol)];}

private:
    int _Nant, _Nchan;

    bool* _flags; // Nant x Nchan x 2 array
    std::complex<float>* _gains; // Nant x Nchan x 2 array

    int index(int ant, int chan, int pol);
};

#endif // BCALTABLE_H_

