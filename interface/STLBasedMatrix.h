#ifndef __STLBasedMatrix_h__
#define __STLBasedMatrix_h__

//#include "TString.h"
#include <string>
#include <cstring>
#include <vector>
#include <fstream>

//#define N_STATES 13

using namespace std;

// Creates a matrix with nRows x nCols initialized to 0

template <class T>
class STLBasedMatrix {
  vector<vector<T> > m;
public:
  STLBasedMatrix(unsigned int x, unsigned int y) {
    m.resize(x, vector<T>(y,0));
  }
  STLBasedMatrix(STLBasedMatrix &other) {
    for(T::iterator i=other.begin(); i < other.end(); ++i) i->
  }
  class matrix_row {
    vector<T>& row;
  public:
    matrix_row(vector<T>& r) : row(r) {
    }
    T& operator[](unsigned int y) {
      return row.at(y);
    }
  };
  matrix_row& operator[](unsigned int x) {
    return matrix_row(m.at(x));
  }
  STLBasedMatrix<T>& operator =(const STLBasedMatrix<T>& rhs);
  friend STLBasedMatrix<T>& operator +(const STLBasedMatrix<T>& lhs, const STLBasedMatrix<T>& rhs);
  friend STLBasedMatrix<T>& operator *(const STLBasedMatrix<T>& lhs, const STLBasedMatrix<T>& rhs);
};

#endif
