#pragma once
#include "math.h"
#include "Vector.h"
#include <cstdlib>
#define PI 3.14159265358979324

class CPolyhedron//��������
{
public:
	CPolyhedron(void);
	CPolyhedron(double ppara[3], double rrr[2]);
	~CPolyhedron(void);

public:
	double pA[3];//x,y,z�᳤�ȣ�a=b;  p:ָ�� |x/a|^2p+|y/b|^2p+|z/c|^2p=1;
	double area;
	double rout, rin;//�����������뾶
	CVector center;	//���ģ�ŷ������ϵ��
	CVector center_L;//����	������������ϵ��
	CVector e[2];//���᷽��e[2]Ϊ����Z������������Ϊ��ֱ�������ƽ���ڵĻ��ഹֱ����������������

	CVector center_L_old;
	CVector e_old[2];
	int ix[2], icell;
	int kind;

public:
	double sgnp(double x);
	void Getroutrin();
	void ChangeSize(double ratio);
	void ChangePosations(CVector center0);
	void ChangeOritations(CVector e0[2]);
	void ChangePosAndOris(CVector center0, CVector e0[2]);

	void Update_EtoL(double L_inv[3][3]);
	void Update_LtoE(double lambda[3][3]);
	void Copyfrom(CPolyhedron * pp);
	void Jump(CVector TV);//process periodic boundary condition	


	void TranslateRandom(double a, double lambda[3][3]);
	void RotateRandom(double a);
	void RetainTrans(double lambda[3][3]);
	void RetainRotate();
};

