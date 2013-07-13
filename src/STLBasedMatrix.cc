#include "../interface/STLBasedMatrix.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <fstream>

#include <TMath.h>

using namespace std;

template <class T>
STLBasedMatrix<T>::STLBasedMatrix(STLBasedMatrix<T> &other) {
    for(vector<vector<T> >::iterator i=other.begin(); i < other.end(); ++i) i->assign(i->begin(), i->end());
}

template <class T>
STLBasedMatrix<T>& STLBasedMatrix<T>::operator =(const STLBasedMatrix<T>& rhs){
	if( m.size() != rhs.m.size() || m[0].size() != rhs.m[0].size() ){
		cerr << "Error: trying to assign matrices of different size!" << endl;
		exit(1);
	}
	if(this != rhs){
		int row=0, col;
		for(vector<vector<T> >::iterator i=rhs.begin; i < rhs.end(); ++i){
			col=0;
			for(vector<T>::iterator j=i.begin; j < i.end(); ++j){
				m[row][col]=rhs.m[row][col];
				++col;
			}
			++row;
		}
	}
	return *this;
}


template <class T>
STLBasedMatrix<T>& operator +(const STLBasedMatrix<T>& lhs, const STLBasedMatrix<T>& rhs){
	STLBasedMatrix<T> *result = new STLBasedMatrix<T>(lhs.m.size(), lhs.m[0].size());
	int row=0, col=0;
	if(m.size() != rhs.m.size() || m[0].size() != rhs.m[0].size() ){
		cerr << "Error: trying to sum matrices of incompatible sizes!" << endl;
		exit(1);
	}
	col=0;
	int row=0, col;
	for(vector<vector<T> >::iterator i=rhs.begin; i < rhs.end(); ++i){
		col=0;
		for(vector<T>::iterator j=i.begin; j < i.end(); ++j){
			result->m[row][col]=lhs.m[row][col]+rhs.m[row][col];
			++col;
		}
		++row;
	}
	return *result;
}


template <class T>
STLBasedMatrix<T>& operator *(const STLBasedMatrix<T>& lhs, const STLBasedMatrix<T>& rhs){
	STLBasedMatrix<T> *result = new STLBasedMatrix<T>(lhs.m.size(), rhs.m[0].size());
	int row=0, col=0;
	if( lhs.m[0].size() != rhs.m.size() ){
		cerr << "Error: trying to multiply matrices of incompatible sizes!" << endl;
		exit(1);
	}
	col=0;
	int row=0, col;
	for(vector<vector<T> >::iterator i=lhs.begin; i < lhs.end(); ++i){
		col=0;
		for(vector<T>::iterator j=i.begin; j < i.end(); ++j){
			int k=0;
			for(vector<T>::iterator pivot=i.begin; pivot < i.end(); ++pivot){
				result->m[row][col]+=lhs.m[row][k]*rhs.m[k][col];
				++k;
			}
			++col;
		}
		++row;
	}
	return *result;
}
