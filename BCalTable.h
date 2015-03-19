#ifndef BCALTABLE_H_
#define BCALTABLE_H_

#include <complex>
#include <vector>

class BCalTable {
public:
    BCalTable(char const* filename);
    ~BCalTable() {}

    int Nant() const {return _Nant;}
    int Nchan() const {return _Nchan;}
    std::vector<char> const& flags() const {return _flags;}
    std::vector<std::complex<float> > const& gains() const {return _gains;}

private:
    int _Nant, _Nchan;
    std::vector<char> _flags; // Nant x Nchan x 2
    std::vector<std::complex<float> > _gains; // Nant x Nchan x 2
};

#endif // BCALTABLE_H_

