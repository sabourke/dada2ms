#ifndef JCALTABLE_H_
#define JCALTABLE_H_

#include <complex>
#include <vector>

class JCalTable {
public:
    JCalTable(char const* filename);
    ~JCalTable() {}

    int Nant() const {return _Nant;}
    int Nchan() const {return _Nchan;}
    std::vector<char> const& flags() const {return _flags;}
    std::vector<std::complex<float> > const& gains() const {return _gains;}

private:
    int _Nant, _Nchan;
    std::vector<char> _flags; // Nant x Nchan
    std::vector<std::complex<float> > _gains; // Nant x Nchan x 2 x 2

    int index(int ant, int chan, int pol1, int pol2);
};

#endif // JCALTABLE_H_

