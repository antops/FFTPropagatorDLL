#include "FFTPropagator.h"
#include "../FFTDIDLL/FFTDI.h"
#include "../FFTASRotationDLL/FFTASRotation.h"
#include "../Util/Matrix4D.h"
#include "../Util/Vector3D.h"
#include "../Util/Position3D.h"
#include "../Util/FieldFile.h"



FFTPropagator::FFTPropagator(int i)
{
	Eu1 = NULL;	Ev1 = NULL;	En1 = NULL;
	Hu1 = NULL;	Hv1 = NULL;	Hn1 = NULL;
	this->returnFloat = NULL;
	this->user = NULL;
	inputfile = true;
}

FFTPropagator::~FFTPropagator()
{
	FreeCal();
}

void FFTPropagator::setInput(FieldBase * input)
{
	inputptr = input;
	inputfile = false;
}

void FFTPropagator::Allocate(int Nu, int Nv) {
	Eu1 = Allocate_2D(Eu1, Nu, Nv);
	Ev1 = Allocate_2D(Ev1, Nu, Nv);
	En1 = Allocate_2D(En1, Nu, Nv);
	Hu1 = Allocate_2D(Hu1, Nu, Nv);
	Hv1 = Allocate_2D(Hv1, Nu, Nv);
	Hn1 = Allocate_2D(Hn1, Nu, Nv);
}

void FFTPropagator::FreeCal() {
	if(Eu1 != NULL)	Free_2D(Eu1);
	if(Ev1 != NULL)	Free_2D(Ev1);
	if(En1 != NULL)	Free_2D(En1);
	if(Hu1 != NULL) Free_2D(Hu1);
	if(Hv1 != NULL)	Free_2D(Hv1);
	if(Hn1 != NULL)	Free_2D(Hn1);
	Eu1 = NULL;		Ev1 = NULL;		En1 = NULL;
	Hu1 = NULL;		Hv1 = NULL;		Hn1 = NULL;
}

void FFTPropagator::setInputFile(const std::string & file)
{
	inputFieldFile = file;
	inputfile = true;
}

bool FFTPropagator::calculate(double fre, GraphTrans GTt, int mode)
{
	//mode 是控制先传再转 还是先转再传
	//mode == 0 先转再传 //适合镜前场传播
	//mode == 1 先传再转 //适合镜后场传播
	

	double theta;
	double phi;
 	GraphTrans GT0;
	Vector3 U0, V0, N0;
	Vector3 Ut, Vt, Nt;
	double X0, Y0, Z0;
	double Xt, Yt, Zt;
	double Up, Vp, Np;
	int Nu, Nv;
	double du,dv;

	FieldBase inputField;
	//读入入射场
	if (inputfile) {
		fromFileBinary(inputFieldFile, inputField);
		GT0 = inputField.graphTransField;
		Nu = inputField.Ex.size();
		Nv = inputField.Ex[0].size();
		du = inputField.ds_x;
		dv = inputField.ds_y;
		Allocate(Nu, Nv);
		for (int i = 0; i < Nu; i++) {
			for (int j = 0; j < Nv; j++) {
				Eu1[i][j] = inputField.Ex[i][j];
				Ev1[i][j] = inputField.Ey[i][j];
			}
		}
	}
	else {
		FieldBase* inFieldPtr = (FieldBase *)inputptr;
		GT0 = inFieldPtr->graphTransField;
		cout << GT0.getTrans_x() << " " << GT0.getTrans_y() << " " << GT0.getTrans_z() << " "
			<< GT0.getRotate_x() << " " << GT0.getRotate_y() << " " << GT0.getRotate_z() << " "
			<< GT0.getRotate_theta() << endl;

		Nu = inFieldPtr->N_width;
		Nv = inFieldPtr->M_depth;
		du = inFieldPtr->ds_x;
		dv = inFieldPtr->ds_y;

		Allocate(Nu, Nv);
		for (int i = 0; i < Nu; i++) {
			for (int j = 0; j < Nv; j++) {
				Eu1[i][j] = inFieldPtr->Ex[i][j];
				Ev1[i][j] = inFieldPtr->Ey[i][j];
			}
		}

	}

	//为输出做准备
	fieldout.ds_x = du;
	fieldout.ds_y = dv;
	fieldout.N_width = Nu;
	fieldout.M_depth = Nv;
	fieldout.graphTransField = GTt;
	

	cout << fieldout.graphTransField.getTrans_x() << " "
		<< fieldout.graphTransField.getTrans_y() << " "
		<< fieldout.graphTransField.getTrans_z() << " "
		<< fieldout.graphTransField.getRotate_x() << " "
		<< fieldout.graphTransField.getRotate_y() << " "
		<< fieldout.graphTransField.getRotate_z() << " "
		<< fieldout.graphTransField.getRotate_theta() << endl;

	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(10, user);
	}

	X0 = GT0.getTrans_x();	Y0 = GT0.getTrans_y();	Z0 = GT0.getTrans_z();
	U0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z()) * Vector3(1.0, 0.0, 0.0);
	V0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z()) * Vector3(0.0, 1.0, 0.0);
	N0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z()) * Vector3(0.0, 0.0, 1.0);

	Xt = GTt.getTrans_x();	Yt = GTt.getTrans_y();	Zt = GTt.getTrans_z();
	Ut = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z()) * Vector3(1.0, 0.0, 0.0);
	Vt = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z()) * Vector3(0.0, 1.0, 0.0);
	Nt = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z()) * Vector3(0.0, 0.0, 1.0);

	if (N0.Dot(Nt) < 0) return false;	//转大发了，不算了
	//计算传播时相对的位移偏移
	theta = acos(N0.Dot(Nt));
	phi = atan2((Nt - N0*(Nt.Dot(N0))).Dot(V0), (Nt - N0*(Nt.Dot(N0))).Dot(U0));
	Up = Xt - X0;	Vp = Yt - Y0;	Np = Zt - Z0;
	Vector3 UVNt = Matrix4D::getRotateMatrix(-GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z())*Vector3(Up, Vp, Np);
	Up = UVNt.x;	Vp = UVNt.y;	Np = UVNt.z;
	
	Vector3 tU, tV, tN;
	tU = Matrix4D::getInvRotateMatrixRad(theta,phi)*U0;
	tV = Matrix4D::getInvRotateMatrixRad(theta, phi)*V0; 
	tN = Matrix4D::getInvRotateMatrixRad(theta, phi)*N0;

	//旋转和传播所需的变量都写好了
	if (mode == 0) {//先转再传
		//旋转部分
		if (Nt.Dot(N0) > 0.9999) {}//不旋转了
		else {//旋转-要插值
			FFTASRotation* ASR = new FFTASRotation;
			ASR->SetParas(fre, Nu, Nv, du,dv);
			ASR->SetRotationParas(GT0, GTt);
			ASR->SetInput(Eu1, Ev1);

			if (returnFloat) // 如果没有注册则不回调
			{
				returnFloat(20, user);
			}

			ASR->PerformRotate();

			if (returnFloat) // 如果没有注册则不回调
			{
				returnFloat(40, user);
			}

			ASR->output(Eu1, Ev1, En1);
			//调试检查用
			GT0.setGraphTransPar(X0, Y0, Z0, GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z(), GTt.getRotate_theta());
			//inputField.setPlane(GT0, ds);
			//inputField.setField(Eu1, Ev1, Eu1, Ev1, Eu1, Ev1);
			//inputField.saveall(string("./PVVD/out.txt"));
			delete ASR;
		}

		complex<double> test1 = Eu1[(Nu - 1) / 2][(Nv - 1) / 2];


		if (returnFloat) // 如果没有注册则不回调
		{
			returnFloat(45, user);
		}

		//传播部分
		FFTDI * FDI = new FFTDI(fre, Up, Vp, Np, Nu, Nv);
		cout << fre << " "<< du<<" "<<dv<< endl;
		FDI->Setds(du,dv);
		FDI->SetInput(Eu1, Ev1);
		FDI->StartCal();

		FDI->output(Eu1, Ev1, En1, Hu1, Hv1, Hn1);
		delete FDI;

		if (returnFloat) // 如果没有注册则不回调
		{
			returnFloat(80, user);
		}

	}
	else if (mode == 1) {//先传再转

		//传播部分 //这里的偏移不对！
		FFTDI * FDI = new FFTDI(fre, Up, Vp, Np, Nu, Nv);
		FDI->Setds(du,dv);
		FDI->SetInput(Eu1, Ev1);
		FDI->StartCal();

		FDI->output(Eu1, Ev1, En1, Hu1, Hv1, Hn1);
		delete FDI;

		if (returnFloat) // 如果没有注册则不回调
		{
			returnFloat(45, user);
		}
		//旋转部分
		if (Nt.Dot(N0) > 0.9999) {}//不旋转了
		else {//旋转-要插值
			FFTASRotation* ASR = new FFTASRotation;
			ASR->SetParas(fre, Nu, Nv, du,dv);
			ASR->SetRotationParas(GT0, GTt);
			ASR->SetInput(Eu1, Ev1);

			if (returnFloat) // 如果没有注册则不回调
			{
				returnFloat(20, user);
			}

			ASR->PerformRotate();

			if (returnFloat) // 如果没有注册则不回调
			{
				returnFloat(40, user);
			}

			ASR->output(Eu1, Ev1, En1);
			//调试检查用
			delete ASR;
		}


		if (returnFloat) // 如果没有注册则不回调
		{
			returnFloat(80, user);
		}


	}

	fieldout.Ex.resize(Nu);	fieldout.Ey.resize(Nu);	fieldout.Ez.resize(Nu);
	fieldout.Hx.resize(Nu);	fieldout.Hy.resize(Nu);	fieldout.Hz.resize(Nu);
	for (int i = 0; i < Nu; i++) {
		fieldout.Ex[i].resize(Nv);	fieldout.Ey[i].resize(Nv);	fieldout.Ez[i].resize(Nv);
		fieldout.Hx[i].resize(Nv);	fieldout.Hy[i].resize(Nv);	fieldout.Hz[i].resize(Nv);
	}
	for (int i = 0; i < Nu; i++) {
		for (int j = 0; j < Nv; j++) {
			fieldout.Ex[i][j] = Eu1[i][j];
			fieldout.Ey[i][j] = Ev1[i][j];
			fieldout.Ez[i][j] = En1[i][j];
			fieldout.Hx[i][j] = Hu1[i][j];
			fieldout.Hy[i][j] = Hv1[i][j];
			fieldout.Hz[i][j] = Hn1[i][j];
		}
	}

	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(100, user);
	}
	return true;
	//pvva.getField(&inputField);
}

void FFTPropagator::SetReturnFloat(void(*returnFloat)(int, void *), void * _user)
{
	this->returnFloat = returnFloat;
	this->user = _user;
}

void FFTPropagator::getCalculatedFieldBase(FieldBase &_out) {
	_out = fieldout;
}



