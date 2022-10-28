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
void Initialize2();//读取综合大文件

void InitialParaAndSpace();//初始化一些参数常量+为PPD申请空间
void askforPPDspace();//为PPD申请空间
void releaseallspace();//释放所有空间
void askforpolyspace();//颗粒多元空间申请
void releasepolyspace();//颗粒多元空间释放
void askfororderspace();
void releaseorderspace();
void Generatepolyhedras();//生成颗粒
void Getboundary();//边界初始化
double GetDensity();//求填充率
void Randomall();//所有颗粒的位置、方向设为随机
void Randomallposition();//所有颗粒的位置设为随机
void Releasetononoverlap(double release_ratio);//将结构松弛到无重叠状态
void OutPackinginfor();//输出基本信息到屏幕
void Lambda_inv();
void Rescale(double bilv);
void PeriodCheckP(CVector &pp);
void PeriodicalCheck(CPolyhedron * pp);

//求颗粒形状信息
void OutPolyInfor();
void InputPolyInfor();//通过手动输入，得到多元颗粒具体情况
void GetPolyInfor();//得到多元颗粒具体情况
void Getroutrin(double ppara[3], double rrout[2]);

double GetArea(double ppara[3]);

double Ff(double ppE[5], double cita, double fai);
double calpow(double x, double y);
double GetERY(double pP[2]);
double GetERZ(double pP[2]);
void GetEllipsoidity(double ppara[5]);

//求交，求配位数
void contact(double bilv);//将颗粒边长扩大至原来的bilv倍，求配位数
void ErrorAnalysis1();//单个容差分析log
void ErrorAnalysis2();//单个容差分析
int GlobalOverlapCheck();//检测填充是否有重叠，有返回1，无返回0；
int SingleOverlapCheck(int Mi);//检测单个颗粒是否与其它颗粒有重叠，有返回重叠编号，无返回0；
int doubleOverlapCheck(int Mi, int Ni);

//输出结构文件
void OutputPAC(int Counter);
void OutputPAC();
void OutputPACsample(int Counter);

//输出SCR文件//输出pov-ray文件
void SCRmn(int mp, int np);
void SCRmn(int mp, int np, int aabb, CPolyhedron * psp0, CPolyhedron * psp1);
void SCR();
void SCRcenter();//质心
void SCRmodel(double ppara[5]);
void povray();
void povrayeach();

//位置随机性，均匀性监测,RDF
void IsotropyTest();//均匀性监测
void RadialDistributionFunction();//RDF，径向分布函数
void OutputPairDis();//质心距离输出

					 //周围颗粒建系
void getppd();//particle particle deitence
void updateppd(int m);//第m个颗粒移动，更新ppd和所有
void getpairdismin(int m, int i);//找到m颗粒最近的第i个距离
void outinfor(int Numm);//信息输出检测

//有序参数
void CalculateS();
double CalculateP4global();
double CalculateP6global();



//形状函数求交所用
double sgn(double x);
double makeMnonsing(double eps, double scale, double M[3][3]);//加主元使矩阵可逆，并返回改变后的行列式
double SuperballFunction(double ppara[5], CVector X);//x^2p+y^2p+z^2p;
void SuperballFunctionDall(double ppara[5], CVector X, CVector &xn, double ddx[3][3]);//一阶导数（向量）和二阶导数（矩阵）
int Inspect(CPolyhedron * psp0, CPolyhedron * psp1);//返回两颗粒是否相交，相交返回1，不相交返回0.前提：之前已经进行过粗略的内切外接判断。
int Inspectnum(CPolyhedron * psp0, CPolyhedron * psp1);//多次求交判断
double Inspectvalue(CPolyhedron * pa, CPolyhedron * pb);//返回两颗粒的重叠势能

														//voronoi剖分
void getsurpoint(double ppara[5]);//颗粒表面点，申请空间，求局部拉格朗日坐标
void OutputVoroPoint(int num);//输出用于voronoi剖分的表面离散点
int VoroVolume(int num);//统计vorovolume
void OutputVoroPointC(); //输出用于voronoi剖分的质心点
int VoroVolumeC();//统计vorovolume
