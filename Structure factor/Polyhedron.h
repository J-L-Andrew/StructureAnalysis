#pragma once
#include "math.h"
#include "Vector.h"
#include <cstdlib>
#define PI 3.14159265358979324

class CPolyhedron//��������
{
public:
	CPolyhedron(void);
	CPolyhedron(double ppara[5], double rrr[2]);
	~CPolyhedron(void);

public:
	double pA[3], pP[2];//x,y,z�᳤�ȣ�a=b;  p:ָ�� |x/a|^2p+|y/b|^2p+|z/c|^2p=1;
	double volume;
	double rout, rin;//�����������뾶
	CVector center;	//���ģ�ŷ������ϵ��
	CVector center_L;//����	������������ϵ��
	CVector e[3];//���᷽��e[2]Ϊ����Z������������Ϊ��ֱ�������ƽ���ڵĻ��ഹֱ����������������

	CVector center_L_old;
	CVector e_old[3];
	int ix[3], icell;
	int kind;

public:
	double sgnp(double x);
	void Getroutrin();
	void ChangeSize(double ratio);
	void ChangePosations(CVector center0);
	void ChangeOritations(CVector e0[3]);
	void ChangePosAndOris(CVector center0, CVector e0[3]);

	void Update_EtoL(double L_inv[3][3]);
	void Update_LtoE(double lambda[3][3]);
	void Copyfrom(CPolyhedron * pp);
	void Jump(CVector TV);//process periodic boundary condition	

	void Rotate(CVector V, double a);//��V������ת�Ƕ�a
	void Randomize(CVector a[3]);
	void RandomizeCenter(CVector a[3]);
	void RandomizeOrient();

	void TranslateRandom(double a, double lambda[3][3]);
	void RotateRandom(double a);
	void RetainTrans(double lambda[3][3]);
	void RetainRotate();
};

