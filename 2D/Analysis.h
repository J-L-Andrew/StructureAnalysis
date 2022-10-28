#include <iostream>
#include "fstream"
#include <sstream>
#include "iomanip"
#include "string" 
#include <algorithm>

#include "Vector.h"
#include "Matrix.h"
#include "Superball.h"

void ReadAddress();
void Initialize2();//��ȡ�ۺϴ��ļ�

void InitialParaAndSpace();//��ʼ��һЩ��������+ΪPPD����ռ�
void askforPPDspace();//ΪPPD����ռ�
void releaseallspace();//�ͷ����пռ�
void askforpolyspace();//������Ԫ�ռ�����
void releasepolyspace();//������Ԫ�ռ��ͷ�
void askfororderspace();
void releaseorderspace();
void Generatepolyhedras();//���ɿ���
void Getboundary();//�߽��ʼ��
double GetDensity();//�������
void Randomall();//���п�����λ�á�������Ϊ���
void Randomallposition();//���п�����λ����Ϊ���
void Releasetononoverlap(double release_ratio);//���ṹ�ɳڵ����ص�״̬
void OutPackinginfor();//���������Ϣ����Ļ
void Lambda_inv();
void Rescale(double bilv);
void PeriodCheckP(CVector &pp);
void PeriodicalCheck(CPolyhedron * pp);

//�������״��Ϣ
void OutPolyInfor();
void InputPolyInfor();//ͨ���ֶ����룬�õ���Ԫ�����������
void GetPolyInfor();//�õ���Ԫ�����������
void Getroutrin(double ppara[3], double rrout[2]);

double GetArea(double ppara[3]);

double Ff(double ppE[5], double cita, double fai);
double calpow(double x, double y);
double GetERY(double pP[2]);
double GetERZ(double pP[2]);
void GetEllipsoidity(double ppara[5]);

//�󽻣�����λ��
void contact(double bilv);//�������߳�������ԭ����bilv��������λ��
void ErrorAnalysis1();//�����ݲ����log
void ErrorAnalysis2();//�����ݲ����
int GlobalOverlapCheck();//�������Ƿ����ص����з���1���޷���0��
int SingleOverlapCheck(int Mi);//��ⵥ�������Ƿ��������������ص����з����ص���ţ��޷���0��
int doubleOverlapCheck(int Mi, int Ni);

//����ṹ�ļ�
void OutputPAC(int Counter);
void OutputPAC();
void OutputPACsample(int Counter);

//���SCR�ļ�//���pov-ray�ļ�
void SCRmn(int mp, int np);
void SCRmn(int mp, int np, int aabb, CPolyhedron * psp0, CPolyhedron * psp1);
void SCR();
void SCRcenter();//����
void SCRmodel(double ppara[5]);
void povray();
void povrayeach();

//λ������ԣ������Լ��,RDF
void IsotropyTest();//�����Լ��
void RadialDistributionFunction();//RDF������ֲ�����
void OutputPairDis();//���ľ������

					 //��Χ������ϵ
void getppd();//particle particle deitence
void updateppd(int m);//��m�������ƶ�������ppd������
void getpairdismin(int m, int i);//�ҵ�m��������ĵ�i������
void outinfor(int Numm);//��Ϣ������

//�������
void CalculateS();
double CalculateP4global();
double CalculateP6global();



//��״����������
double sgn(double x);
double makeMnonsing(double eps, double scale, double M[3][3]);//����Ԫʹ������棬�����ظı�������ʽ
double SuperballFunction(double ppara[5], CVector X);//x^2p+y^2p+z^2p;
void SuperballFunctionDall(double ppara[5], CVector X, CVector &xn, double ddx[3][3]);//һ�׵������������Ͷ��׵���������
int Inspect(CPolyhedron * psp0, CPolyhedron * psp1);//�����������Ƿ��ཻ���ཻ����1�����ཻ����0.ǰ�᣺֮ǰ�Ѿ����й����Ե���������жϡ�
int Inspectnum(CPolyhedron * psp0, CPolyhedron * psp1);//������ж�
double Inspectvalue(CPolyhedron * pa, CPolyhedron * pb);//�������������ص�����

														//voronoi�ʷ�
void getsurpoint(double ppara[5]);//��������㣬����ռ䣬��ֲ�������������
void OutputVoroPoint(int num);//�������voronoi�ʷֵı�����ɢ��
int VoroVolume(int num);//ͳ��vorovolume
void OutputVoroPointC(); //�������voronoi�ʷֵ����ĵ�
int VoroVolumeC();//ͳ��vorovolume
