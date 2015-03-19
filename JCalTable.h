#ifndef JCALTABLE_H_
#define JCALTABLE_H_

#include <complex>

class JCalTable {
public:
    JCalTable(char const* filename);
    ~JCalTable();

    int Nant() {return _Nant;}
    int Nchan() {return _Nchan;}
    bool flag(int ant, int chan) {return _flags[ant*_Nchan + chan];}
    std::complex<float> gain(int ant, int chan, int pol1, int pol2) {return _gains[index(ant,chan,pol1,pol2)];}

private:
    int _Nant, _Nchan;

    bool* _flags; // Nant x Nchan array
    std::complex<float>* _gains; // Nant x Nchan x 2 x 2 array

    int index(int ant, int chan, int pol1, int pol2);
};

#endif // JCALTABLE_H_

