#pragma once

#ifndef FFTPROPAGATOR_H
#define FFTPROPAGATOR_H

#include <string>
#include <complex>
#include <vector>
#include "../Util/FieldBase.h"
#include "../Util/GraphTrans.h"

using namespace std;
class _declspec(dllexport) FFTPropagator
{
public:
	FFTPropagator(int i = 1);
	~FFTPropagator();
	void setInput(FieldBase* input);
	void setInputFile(const std::string & file);
	bool calculate(double fre, GraphTrans GTt, int mode = 0);//0 先转再传，1 先传再转
	void SetReturnFloat(void(*returnFloat)(int, void*), void*_user);// 注册回调函数
	void getCalculatedFieldBase(FieldBase &_out);

private:
	void Allocate(int Nu, int Nv);
	void FreeCal();

private:
	complex<double> ** Eu1;
	complex<double> ** Ev1;
	complex<double> ** En1;
	complex<double> ** Hu1;
	complex<double> ** Hv1;
	complex<double> ** Hn1;
	FieldBase fieldout;
	void* inputptr;
	bool inputfile;
	string inputFieldFile;
	void(*returnFloat)(int, void*);
	void *user; // 回调函数的类指针

};

#endif
