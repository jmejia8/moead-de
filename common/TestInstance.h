// TestInstance.h: interface for the CTestInstance class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include <math.h>

#define PI  3.1415926535897932384626433832795

using namespace std;

class CTestInstance  
{
public:
	CTestInstance();
	virtual ~CTestInstance();
	

	void DTLZ1(vector<double> &x, vector<double> &f, const unsigned int nx);
	void DTLZ2(vector<double> &x, vector<double> &f, const unsigned int nx);


	void ZDT3(vector<double> &x, vector<double> &f, const unsigned int nx);
	void ZDT2(vector<double> &x, vector<double> &f, const unsigned int nx);
	void ZDT1(vector<double> &x, vector<double> &f, const unsigned int nx);


    void TEC09_LZ1(vector<double> &x, vector<double> &f, const unsigned int nx);
    void TEC09_LZ2(vector<double> &x, vector<double> &f, const unsigned int nx);
    void TEC09_LZ3(vector<double> &x, vector<double> &f, const unsigned int nx);
	void TEC09_LZ4(vector<double> &x, vector<double> &f, const unsigned int nx);
    void TEC09_LZ5(vector<double> &x, vector<double> &f, const unsigned int nx);
    void TEC09_LZ6(vector<double> &x, vector<double> &f, const unsigned int nx);
    void TEC09_LZ7(vector<double> &x, vector<double> &f, const unsigned int nx);
    void TEC09_LZ8(vector<double> &x, vector<double> &f, const unsigned int nx);
    void TEC09_LZ9(vector<double> &x, vector<double> &f, const unsigned int nx);

};


