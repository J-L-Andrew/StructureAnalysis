#include "StructureAnalysis.h"
using namespace std;
#define ERROR1		1.000E-6
#define ERROR2		1.000E-12
#define PI			3.14159265358979324
#define Maxlocal	26

CPolyhedron **polyhedra;
CPolyhedron *pTmp[2];

int	   Npolyhedron;//颗粒总数
int    NKind;//组分种数
int	   *NpEach;//各组分颗粒数
double *VolfEach;//各组分体积分数
double *SizeEach;//各组分相对大小
double *VEach;//各组分体积
double **PPara;//形状参数
double	Pcdmax, Pidmin; //最大外接球直径,最小内切球半径

//int	   Npolyhedron;//颗粒总数
//int    NKind;//组分种数
//int	   *NpEach;//各组分颗粒数
//double *VolfEach;//各组分体积分数
//double **PPara;//形状参数
//double	Pcdmax, Pidmin; //最大外接球直径,最小内切球半径

double SumVolume;//颗粒总体积
double V_box;//盒子总体积
double PackingDensity;
double Lambda[3][3], L_inv[3][3];
CVector att[3];//att[0],att[1],att[2]表示三个边界向量

double sblengthmin;//边界的最短长度的平方，最近颗粒距离的最大取值
double DISMAX;//=sblengthmin*10000；
double S4random, S4length;
double Q6localrandom[Maxlocal], Q6locallength[Maxlocal];
double S4localrandom[Maxlocal], S4locallength[Maxlocal];
double **pairdis;//颗粒间距离
double **pairdismin;//颗粒间距离
int **pairnumber;//颗粒间编号
int **pairnumbermin;//颗粒间编号
CVector **pairdcenter;//颗粒间质心向量

CVector *SurfPot;//颗粒表面的点，求交用
int NSuP; //颗粒表面点数

int Counterfile = 0;
ifstream input;
ofstream result("结构因子.txt");

double* S_k;
int* Num;
double delta_k = 0.04;
int N_sam = int(6.0 / delta_k);


/*读取综合大文件分析*/
void main()
{
	srand((unsigned)time(NULL));
	int STRUCTURENUM, iteam, Ponum;
	string pline;
	result << fixed << setprecision(8);

	Ponum = 20; //cout << "please input surface point number:\t"; cin >> Ponum;
	for (iteam = 0; iteam < 1; iteam++)
	{
		cout << "iteam:\t" << iteam << endl;
		input.open("H:\\totalstructure.txt");	//input.open("totalstructure.txt");////input.open("packingsample.txt");
		if (!input.good()) { cout << "totalstructure.txt不存在" << endl;	break; }	
		result << "Counter\tNsp\tPD\tNKind\t";
		result << "Sratio\tVolfL\tNumfL\tSize3\tVolf3\tS1\ts2\ts3\tv1\tv2\tv3\tn1\tn2\tn3\n";



		//result << "Counter\tNsp\tPD\tNKind\tSupsp\tSupvolf\tSupsize\tNmix\tmixratio\tmixsp\tmixVolF\tmixsize\t";
		//result << "S4\tS4super\tS4superN\tQ6local\tS4local" << endl;
		//result << "num\tNSuP\tavVlocal\tdetaVlocal\tPDlocal\tdetaPDlocal\tNsuper\tNmix\tVsuper\tVSuperd\tVmix\tVmixd\tPDsuper\tPDsuperd\tPDmix\tPDmixd\n";
		Counterfile = -1;
		do
		{
			getline(input, pline);
			if (pline == "NEW")
			{
				input >> STRUCTURENUM;
				Counterfile++;
				cout << Counterfile << endl;
				Initialize2(); //cout << "yes!\n";
				
				{					
					
					result << Counterfile << "\t" << Npolyhedron << "\t" << PackingDensity << "\t" << NKind << "\t";					
					OutPolyInfor();				
					Structure_factor(S_k);
					
					for (int j = 0; j < N_sam; j++)
					{
						result << (j+0.5) * delta_k << "\t" << S_k[j] << endl;
					}	
				}
				releaseallspace();
				releasepolyspace();				
			}
		} while (!input.eof());
		input.close();
	}
	
	result.close();
	cout << "press any key to continue!\n";
	getchar(); //getchar(); //getchar();
}

//初始化
void Initialize2()
{
	int i, j;
	string line;
	double nonesize;

	do
	{
		getline(input, line);
		if (line == "PackingSpaceInformation")
		{
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					input >> Lambda[j][i];
		}
		if (line == "Particle&Structure")
		{
			input >> nonesize >> nonesize >> nonesize >> nonesize;
			do
			{
				getline(input, line);
				if (line == "Superellipsoidevery")
				{					
					input >> Npolyhedron >> NKind;
					askforpolyspace();
					for (int i = 0; i < NKind; i++)
					{
						input >> NpEach[i] >> VolfEach[i] >> SizeEach[i];
						//input >> NpEach[i] >> VolfEach[i];
						for (int j = 0; j < 5; j++)
							input >> PPara[i][j];
					}
					InitialParaAndSpace();
					for (int i = 0; i < Npolyhedron; i++)
					{
						for (int j = 0; j < 3; j++)
							input >> polyhedra[i]->center[j];
						for (int j = 0; j < 3; j++)
							for (int k = 0; k < 3; k++)
								input >> polyhedra[i]->e[j][k];
					}
				}
			} while (line != "END");
		}
	} while (line != "ENDOFFILE");

	Getboundary();
	GetPolyInfor();
	//Releasetononoverlap(1.000000001);
	PackingDensity = GetDensity();
}

//初始化辅助函数
void InitialParaAndSpace()
{
	//填充基本信息
	Generatepolyhedras();//申请颗粒空间	
	askforPPDspace();//申请其他空间	

	//有序度参数阈值
	int i;
	double lkav, lksigma;
	S4random = 0.97980 / sqrt(Npolyhedron);
	S4length = 0.27048 / sqrt(Npolyhedron);
	lkav = 0.30972;
	lksigma = 0.38720 / sqrt(Npolyhedron);
	for (i = 0; i < Maxlocal; i++)
	{
		S4localrandom[i] = lkav / sqrt(i + 1.0);
		S4locallength[i] = lksigma / sqrt(i + 1.0);
	}
	lkav = 0.98123 / sqrt(Npolyhedron);
	lksigma = 0.19379 / sqrt(Npolyhedron);
	for (i = 0; i < Maxlocal; i++)
	{
		Q6localrandom[i] = lkav / sqrt(i + 1.0);
		Q6locallength[i] = lksigma / sqrt(i + 1.0);
	}	
}
void askforPPDspace()
{	
	pairdis = new double*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		pairdis[i] = new double[Npolyhedron];
	pairdismin = new double*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		pairdismin[i] = new double[Maxlocal];
	pairnumber = new int*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		pairnumber[i] = new int[Npolyhedron];
	pairnumbermin = new int*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		pairnumbermin[i] = new int[Maxlocal];
	pairdcenter = new CVector*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		pairdcenter[i] = new CVector[Npolyhedron];
}
void releaseallspace()//释放所有空间
{
	for (int i = 0; i < 2; i++) delete pTmp[i];
	for (int i = 0; i < Npolyhedron; i++) delete polyhedra[i];  delete[] polyhedra;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairdis[i]; delete[] pairdis;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairnumber[i]; delete[] pairnumber;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairdismin[i]; delete[] pairdismin;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairnumbermin[i]; delete[] pairnumbermin;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairdcenter[i]; delete[] pairdcenter;
}
void askforpolyspace()//颗粒多元空间申请
{
	NpEach = new int[NKind];
	VolfEach = new double[NKind];
	SizeEach = new double[NKind];
	VEach = new double[NKind];
	PPara = new double*[NKind];

	S_k = new double[N_sam];
	Num = new int[N_sam];

	for (int i = 0; i < NKind; i++)
		PPara[i] = new double[5];
}
void releasepolyspace()//颗粒多元空间释放
{
	delete[] NpEach; delete[] VolfEach;
	delete[] SizeEach; delete[] VEach;
	delete[] S_k; delete[] Num;
	for (int i = 0; i < NKind; i++)delete[] PPara[i]; delete[] PPara;

}
void Generatepolyhedras()//生成颗粒
{
	int i, m, mmi;
	double rrou[2];
	//CPolyhedron * pPolyhedron;

	for (int i = 0; i < 2; i++)
		pTmp[i] = new CPolyhedron();

	polyhedra = new CPolyhedron*[Npolyhedron];
	SumVolume = 0.0; mmi = -1;
	for (m = 0; m < NKind; m++)
	{
		Getroutrin(PPara[m], rrou);
		VolfEach[m] = 0.0;
		for (i = 0; i < NpEach[m]; i++)
		{
			mmi++;
			polyhedra[mmi] = new CPolyhedron(PPara[m], rrou);
			polyhedra[mmi]->kind = m;
			VolfEach[m] += polyhedra[mmi]->volume;
		}
		SumVolume += VolfEach[m];
	}
	for (m = 0; m < NKind; m++)
		VolfEach[m] = VolfEach[m] / SumVolume;
	//颗粒最大外接球直径
	Pcdmax = polyhedra[0]->rout; Pidmin = Pcdmax;
	for (i = 1; i < Npolyhedron; i++)
	{
		if ((polyhedra[i]->rout) > Pcdmax)
			Pcdmax = polyhedra[i]->rout;
		if ((polyhedra[i]->rin) < Pidmin)
			Pidmin = polyhedra[i]->rin;
	}
	Pcdmax = Pcdmax*2.0;
	Pidmin = Pidmin*2.0;
}
void Getboundary()//边界初始化
{
	for (int k = 0; k < 3; k++)
		for (int m = 0; m < 3; m++)
			att[k][m] = Lambda[m][k];
	Lambda_inv();
	for (int i = 0; i < Npolyhedron; i++)
	{
		polyhedra[i]->Update_EtoL(L_inv);
		PeriodicalCheck(polyhedra[i]);
	}
}
double GetDensity()
{
	CVector a_temp = att[0].Cross(att[1]);
	V_box = fabs(a_temp*att[2]);
	return SumVolume / V_box;
}
void Randomall()//所有颗粒的位置、方向设为随机
{
	for (int i = 0; i < Npolyhedron; i++)
	{
		polyhedra[i]->Randomize(att);
		polyhedra[i]->Update_EtoL(L_inv);
	}
}
void Randomallposition()//所有颗粒的位置设为随机
{
	for (int i = 0; i < Npolyhedron; i++)
	{
		polyhedra[i]->RandomizeCenter(att);
		polyhedra[i]->Update_EtoL(L_inv);
	}
}
void Releasetononoverlap(double release__ratio)//将结构松弛到无重叠状态
{
	int xxi = 0;
	while (GlobalOverlapCheck() == 1)
	{
		xxi++;
		printf("%d The Initial Configuration Contains Overlapping Pairs!\n", xxi);
		Rescale(release__ratio);
	}
}
void OutPackinginfor()//输出基本信息到屏幕
{
	cout << "Npolyhedron and Nkind:\t" << Npolyhedron << "\t" << NKind << endl;
	cout << "Number\tVolf\tSize\tVeach\tratio\t" << endl;
	for (int i = 0; i < NKind; i++)
	{
		cout << NpEach[i] << "\t" << VolfEach[i] << "\t" << SizeEach[i] << "\t" << VEach[i] << "\t" << PPara[i][2] / PPara[i][0] << "\t\t\t";
		for (int j = 0; j < 5; j++)cout << PPara[i][j] << "\t";
		cout << endl;
	}

	cout << "Lambda:\n";
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			cout << Lambda[i][j] << " ";
		}
		cout << endl;
	}
	cout << "PackingDensity:\t" << PackingDensity << endl;
}
void Lambda_inv()
{
	int i, j;
	double Linv_m;
	Linv_m = Lambda[0][0] * Lambda[1][1] * Lambda[2][2] + Lambda[0][1] * Lambda[1][2] * Lambda[2][0] + Lambda[0][2] * Lambda[1][0] * Lambda[2][1] - Lambda[0][0] * Lambda[1][2] * Lambda[2][1] - Lambda[0][1] * Lambda[1][0] * Lambda[2][2] - Lambda[0][2] * Lambda[1][1] * Lambda[2][0];
	L_inv[0][0] = Lambda[1][1] * Lambda[2][2] - Lambda[1][2] * Lambda[2][1];
	L_inv[1][0] = -Lambda[1][0] * Lambda[2][2] + Lambda[1][2] * Lambda[2][0];
	L_inv[2][0] = Lambda[1][0] * Lambda[2][1] - Lambda[1][1] * Lambda[2][0];
	L_inv[0][1] = -Lambda[0][1] * Lambda[2][2] + Lambda[0][2] * Lambda[2][1];
	L_inv[1][1] = Lambda[0][0] * Lambda[2][2] - Lambda[2][0] * Lambda[0][2];
	L_inv[2][1] = -Lambda[0][0] * Lambda[2][1] + Lambda[2][0] * Lambda[0][1];
	L_inv[0][2] = Lambda[0][1] * Lambda[1][2] - Lambda[1][1] * Lambda[0][2];
	L_inv[1][2] = -Lambda[0][0] * Lambda[1][2] + Lambda[1][0] * Lambda[0][2];
	L_inv[2][2] = Lambda[0][0] * Lambda[1][1] - Lambda[0][1] * Lambda[1][0];
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			L_inv[i][j] = L_inv[i][j] / Linv_m;
	//	printf("L_inv%lf %lf %lf\n",L_inv[0][0],L_inv[1][1],L_inv[1][1]);
}
void Rescale(double bilv) //if packing containing overlapping pairs, rescale the packing
{
	int i, j;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			Lambda[i][j] = bilv*Lambda[i][j];
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			att[i][j] = Lambda[j][i];
	Lambda_inv();
	for (i = 0; i < Npolyhedron; i++)
		polyhedra[i]->Update_LtoE(Lambda);
}
void PeriodCheckP(CVector &pp)
{
	for (int i = 0; i < 3; i++)
	{
		while (pp[i] >= att[i][i])
		{
			pp[i] = pp[i] - att[i][i];
		}
		while (pp[i] < 0)
		{
			pp[i] = pp[i] + att[i][i];
		}
	}
}
void PeriodicalCheck(CPolyhedron * pp)
{
	for (int i = 0; i < 3; i++)
	{
		while (pp->center_L[i] >= 1.0)
		{
			pp->center_L[i] -= 1.0;
			pp->center -= att[i];
		}
		while (pp->center_L[i] < 0)
		{
			pp->center_L[i] += 1.0;
			pp->center += att[i];
		}
	}
}

//求颗粒形状信息
void OutPolyInfor()
{

	result << SizeEach[1]/ SizeEach[0] << "\t";
	result << VolfEach[1] / (VolfEach[0]+ VolfEach[1]) << "\t";
	result << NpEach[1]*1.0 / (NpEach[0] + NpEach[1]) << "\t";
	result << SizeEach[2] << "\t" << VolfEach[2] << "\t";
	for (int i = 0; i < NKind; i++) result << SizeEach[i] << "\t";
	for (int i = 0; i < NKind; i++) result << VolfEach[i] << "\t";
	for (int i = 0; i < NKind; i++) result << NpEach[i] << "\t";
}
void InputPolyInfor()
{
	int i, j;
	int AllN;//初始颗粒总数，并不是最终个数
	int NKEll;//椭球种类的个数
	double *Ellsratio, *Ellsize, *Ellvolf;//椭球的长细比，大小和体积分数；
	double Supsp, Supsize, Supvolf;//超球的形状，大小和体积分数
	double v0, v1, Aratio, parasph0[5];//大小变化所需

	//input poly infor
	AllN = 420; //cout << "Please input AllN of particles：\t"; cin >> AllN;	 
	Supsp = 4.0; //cout << "Please input Supsp of superball：\t"; cin >> Supsp;
	Supsize = 1.0;// cout << "Please input size of superball：\t"; cin >> Supsize;
	Supvolf = 0.5; //cout << "Please input vol frac of superball：\t"; cin >> Supvolf;
	NKEll = 0;// cout << "Please input kinds of ellipsoids：\t"; cin >> NKEll;
	NKind = NKEll + 1;
	askforpolyspace();
	Ellsratio = new double[NKind];
	Ellsize = new double[NKind];
	Ellvolf = new double[NKind];
	for (i = 0; i < NKEll; i++)
	{
		Ellsratio[i] = 1.0; //cout << "Please input Ratio2 of ellipsoid：\t"; cin >> Ellsratio[i];
		Ellsize[i] = 1.0; //cout << "Please input size of ellipsoid：\t"; cin >> Ellsize[i];
		if (NKEll > 0) Ellvolf[i] = (1.0 - Supvolf) / double(NKEll); //cout << "Please input vol frac of ellipsoid：\t"; cin >> Ellvolf[i];			
	}

	//get shapes, sizes and volumes	
	PPara[0][0] = 1.0; PPara[0][1] = 1.0; PPara[0][2] = 1.0; PPara[0][3] = Supsp; PPara[0][4] = Supsp;
	SizeEach[0] = Supsize;
	for (i = 0; i < NKEll; i++)
	{
		PPara[i + 1][0] = 1.0; PPara[i + 1][1] = 1.0; PPara[i + 1][2] = Ellsratio[i]; PPara[i + 1][3] = 1.0; PPara[i + 1][4] = 1.0;
		SizeEach[i + 1] = Ellsize[i];
	}
	//大小变化，与单位球等体积的等效直径相比	
	parasph0[0] = 1.0; parasph0[1] = 1.0; parasph0[2] = 1.0; parasph0[3] = 1.0; parasph0[4] = 1.0;
	v0 = GetVolume(parasph0); cout << "000\t" << v0 << endl;
	for (i = 0; i < NKind; i++)
	{
		v1 = GetVolume(PPara[i]); cout << i << "\t" << v1 << endl;
		Aratio = pow(v0 / v1, 1.0 / 3.0);
		Aratio = Aratio*SizeEach[i];
		for (j = 0; j < 3; j++)
			PPara[i][j] = PPara[i][j] * Aratio;
		VEach[i] = GetVolume(PPara[i]); cout << i << "\t" << VEach[i] << endl;
	}

	////get  amounts 方案2：总体积一定
	//SumVolume = AllN*v0;
	//NpEach[0] = int(round(SumVolume*Supvolf / VEach[0])); Npolyhedron = NpEach[0];
	//for (i = 0; i < NKEll; i++)
	//{
	//	NpEach[i + 1] =int(round(SumVolume*Ellvolf[i] / VEach[i + 1])); Npolyhedron += NpEach[i + 1];
	//}

	//get  amounts 方案1：总颗粒数一定	
	VolfEach[0] = Supvolf / VEach[0]; SumVolume = VolfEach[0];
	for (i = 0; i < NKEll; i++)
	{
		VolfEach[i + 1] = Ellvolf[i] / VEach[i + 1]; SumVolume += VolfEach[i + 1];
	}
	Npolyhedron = 0;
	for (i = 0; i < NKind; i++)
	{
		NpEach[i] = int(round(AllN*VolfEach[i] / SumVolume)); Npolyhedron += NpEach[i];
	}

	delete[] Ellsratio; delete[] Ellsize; delete[] Ellvolf;
}
void GetPolyInfor()//得到多元颗粒具体情况
{
	int i;
	double v0, v1, Aratio1;
	double parasph0[5], parasph1[5];

	//get sizes	and volumes
	parasph0[0] = 1.0; parasph0[1] = 1.0; parasph0[2] = 1.0; parasph0[3] = 1.0; parasph0[4] = 1.0;
	v0 = GetVolume(parasph0); //cout << "000\t" << v0 << endl;
	for (i = 0; i < NKind; i++)
	{
		parasph1[0] = 1.0; parasph1[1] = PPara[i][1] / PPara[i][0]; parasph1[2] = PPara[i][2] / PPara[i][0]; parasph1[3] = PPara[i][3]; parasph1[4] = PPara[i][4];
		v1 = GetVolume(parasph1); //cout << i << "\t" << v1 << endl;
		Aratio1 = pow(v0 / v1, 1.0 / 3.0);
		SizeEach[i] = PPara[i][0] / Aratio1;
		VEach[i] = GetVolume(PPara[i]);
	}
}
void Getroutrin(double ppara[5], double rrout[2])
{
	int num, numi;//剖分精度
	int i, j;
	double dv, v, u, cos1, sin1, cos2, sin2;
	double *pcos1, *psin1, *pcos2, *psin2;
	double *sgncos1, *sgnsin1, *sgnsin2;
	double rrmax, rrmin, rrr;

	num = 20;
	dv = PI / (double)num;
	NSuP = 2 * (num*num - num + 1);
	SurfPot = new CVector[NSuP];
	pcos1 = new double[2 * num];
	psin1 = new double[2 * num];
	pcos2 = new double[num - 1];
	psin2 = new double[num - 1];
	sgncos1 = new double[2 * num];
	sgnsin1 = new double[2 * num];
	sgnsin2 = new double[num - 1];

	//v,cita (-0.5*PI,0.5*PI)
	for (i = 1; i < num; i++)
	{
		v = -0.5*PI + dv*i;
		cos2 = cos(v);
		sin2 = sin(v);
		pcos2[i - 1] = pow(fabs(cos2), 1.0 / ppara[4]);
		psin2[i - 1] = pow(fabs(sin2), 1.0 / ppara[4]);
		sgnsin2[i - 1] = sgn(sin2);
	}
	//u,fai [-PI,PI)
	for (j = 0; j < 2 * num; j++)
	{
		u = -PI + dv*j;
		cos1 = cos(u);
		sin1 = sin(u);
		pcos1[j] = pow(fabs(cos1), 1.0 / ppara[3]);
		psin1[j] = pow(fabs(sin1), 1.0 / ppara[3]);
		sgncos1[j] = sgn(cos1);
		sgnsin1[j] = sgn(sin1);
	}
	SurfPot[0][0] = 0.0;
	SurfPot[0][1] = 0.0;
	SurfPot[0][2] = -ppara[2];
	numi = 1;
	for (i = 1; i < num; i++)
		for (j = 0; j < 2 * num; j++)
		{
			SurfPot[numi][0] = ppara[0] * sgncos1[j] * pcos1[j] * pcos2[i - 1];
			SurfPot[numi][1] = ppara[1] * sgnsin1[j] * psin1[j] * pcos2[i - 1];
			SurfPot[numi][2] = ppara[2] * sgnsin2[i - 1] * psin2[i - 1];
			numi++;
		}
	SurfPot[numi][0] = 0.0;
	SurfPot[numi][1] = 0.0;
	SurfPot[numi][2] = ppara[2];

	rrmax = SurfPot[0].Length(); rrmin = rrmax;
	for (i = 1; i < NSuP; i++)
	{
		rrr = SurfPot[i].Length();
		if (rrr > rrmax) rrmax = rrr;
		if (rrr < rrmin) rrmin = rrr;
	}
	rrout[0] = rrmax;
	rrout[1] = rrmin;

	delete[] pcos1; delete[] psin1; delete[] pcos2; delete[] psin2;
	delete[] sgncos1; delete[] sgnsin1; delete[] sgnsin2;
	delete[] SurfPot;
}
void GetPInfo(double ppara[5])
{
	double a1, a2, b1, b2, b3;
	double vp, sp, rs, spericity, csxy, csxz, csyz;
	double rry, rrz;

	a1 = tgamma(0.5 / ppara[3]);
	a2 = tgamma(1.0 / ppara[3]);
	b1 = tgamma(0.5 / ppara[4]);
	b2 = tgamma(1.0 / ppara[4]);
	b3 = tgamma(1.5 / ppara[4]);
	vp = 2.0*ppara[0] * ppara[1] * ppara[2] / 3.0 / ppara[3] / ppara[4] * a1*a1*b1*b2 / a2 / b3; //result << vp << "\t";
	sp = GetSurface(ppara); //result << sp << "\t";
	rs = pow(3.0*vp / 4.0 / PI, 1.0 / 3.0);
	spericity = 4.0*PI*rs*rs / sp; //result << spericity << "\t";
	csxy = 1.0*ppara[0] * ppara[1] / ppara[3] * a1*a1 / a2;// result << csxy << "\t";
	csxz = 1.0*ppara[0] * ppara[2] / ppara[4] * b1*b1 / b2;//result << csxz << "\t";
	csyz = 1.0*ppara[1] * ppara[2] / ppara[4] * b1*b1 / b2;//result << csyz << "\t";
	rry = sqrt(csxy*csyz) / csxz; //result << rry << "\t";
	rrz = sqrt(csxz*csyz) / csxy; //result << rrz << "\t";
}
double GetVolume(double ppara[5])
{
	double vp, a1, a2, b1, b2, b3;
	a1 = tgamma(0.5 / ppara[3]);
	a2 = tgamma(1.0 / ppara[3]);
	b1 = tgamma(0.5 / ppara[4]);
	b2 = tgamma(1.0 / ppara[4]);
	b3 = tgamma(1.5 / ppara[4]);
	vp = 2.0*ppara[0] * ppara[1] * ppara[2] / 3.0 / ppara[3] / ppara[4] * a1*a1*b1*b2 / a2 / b3;
	return vp;
}
double GetSurface(double ppara[5])
{
	int i, maxi, M, N;
	double surf0, surf1;
	double eps;

	maxi = 20;
	M = 1000; N = M;
	surf0 = GetSurfaceini(ppara, M, N);
	for (i = 0; i < maxi; i++)
	{
		if (M > 64000) break;
		//cout << "i nd M" << i << "\t" << M << endl;
		M = M * 2; N = M;
		surf1 = GetSurfaceini(ppara, M, N);
		eps = fabs((surf1 - surf0) / surf1);
		cout << "i and M:\t" << i << "\t" << M << "\t" << eps << "\t" << surf0 << "\t" << surf1 << endl;
		if (eps < 5E-3)break;
		surf0 = surf1;
	}
	if (i == maxi)
	{
		cout << "ERROR! surface area not converge!" << endl;
		getchar(); getchar();
	}
	//result << M << "\t" << surf1 << "\t";
	return surf1;
}
double GetSurfaceini(double ppara[5], int M, int N)
{
	int i, j;
	double h, k;
	double pE[5];
	double *x;
	double *y;
	double surf, temp;

	for (i = 0; i < 3; i++) pE[i] = ppara[i];
	for (i = 3; i < 5; i++) pE[i] = 1.0 / ppara[i];
	x = new double[2 * M + 1]; y = new double[2 * N + 1];
	h = PI / 4.0 / M; k = PI / 4.0 / N;
	for (i = 0; i < 2 * M + 1; i++) x[i] = i*h;
	for (j = 0; j < 2 * N + 1; j++) y[j] = j*k;
	x[2 * M] = 10.0; y[2 * N] = 10.0;

	surf = 0.0;
	surf += Ff(pE, x[0], y[0]) + Ff(pE, x[0], y[2 * N]) + Ff(pE, x[2 * M], y[0]) + Ff(pE, x[2 * M], y[2 * N]);

	temp = 0.0;
	for (i = 1; i < M; i++)
		temp += Ff(pE, x[2 * i], y[0]) + Ff(pE, x[2 * i], y[2 * N]);
	for (j = 1; j < N; j++)
		temp += Ff(pE, x[0], y[2 * j]) + Ff(pE, x[2 * M], y[2 * j]);
	surf += 2.0*temp;

	temp = 0.0;
	for (i = 1; i < M + 1; i++)
		temp += Ff(pE, x[2 * i - 1], y[0]) + Ff(pE, x[2 * i - 1], y[2 * N]);
	for (j = 1; j < N + 1; j++)
		temp += Ff(pE, x[0], y[2 * j - 1]) + Ff(pE, x[2 * M], y[2 * j - 1]);
	surf += 4.0*temp;

	temp = 0.0;
	for (i = 1; i < M + 1; i++)
		for (j = 1; j < N; j++)
			temp += Ff(pE, x[2 * i - 1], y[2 * j]);
	surf += 8.0*temp;

	temp = 0.0;
	for (i = 1; i < M; i++)
		for (j = 1; j < N + 1; j++)
			temp += Ff(pE, x[2 * i], y[2 * j - 1]);
	surf += 8.0*temp;

	temp = 0.0;
	for (i = 1; i < M; i++)
		for (j = 1; j < N; j++)
			temp += Ff(pE, x[2 * i], y[2 * j]);
	surf += 4.0*temp;

	temp = 0.0;
	for (i = 1; i < M + 1; i++)
		for (j = 1; j < N + 1; j++)
			temp += Ff(pE, x[2 * i - 1], y[2 * j - 1]);
	surf += 16.0*temp;

	surf = h*k *surf / 9.0;
	surf = surf*8.0;
	delete[]x; delete[]y;
	if (fabs(surf) > 10000000) cout << surf << "\tfW Wrong!\n", getchar(), getchar();
	return surf;
}
double Ff(double ppE[5], double cita, double fai)
{
	CVector Rcita, Rfai;
	double fE, fF, fG, fW;
	double sincita, coscita, sinfai, cosfai;

	if (cita == 10.0) sincita = 1.0, coscita = 0.0;
	else sincita = sin(cita), coscita = cos(cita);
	if (fai == 10.0) sinfai = 1.0, cosfai = 0.0;
	else sinfai = sin(fai), cosfai = cos(fai);

	Rcita[0] = -ppE[0] * ppE[3] * sincita*calpow(coscita, ppE[3] - 1.0)*calpow(cosfai, ppE[4]);
	Rcita[1] = ppE[1] * ppE[3] * coscita*calpow(sincita, ppE[3] - 1.0)*calpow(cosfai, ppE[4]);
	Rcita[2] = 0.0;

	Rfai[0] = -ppE[0] * ppE[4] * calpow(coscita, ppE[3])*sinfai*calpow(cosfai, ppE[4] - 1.0);
	Rfai[1] = -ppE[1] * ppE[4] * calpow(sincita, ppE[3])*sinfai*calpow(cosfai, ppE[4] - 1.0);
	Rfai[2] = ppE[2] * ppE[4] * cosfai*calpow(sinfai, ppE[4] - 1.0);

	fE = Rcita[0] * Rcita[0] + Rcita[1] * Rcita[1] + Rcita[2] * Rcita[2];
	fF = Rcita[0] * Rfai[0] + Rcita[1] * Rfai[1] + Rcita[2] * Rfai[2];
	fG = Rfai[0] * Rfai[0] + Rfai[1] * Rfai[1] + Rfai[2] * Rfai[2];
	fW = fE*fG - fF*fF;	//cout << cita << "\t" << fai << "\t" << fW << endl;
	fW = sqrt(fW);
	return (fW);
}
double calpow(double x, double y)
{
	if (x == 0.0) return 0.0;
	else return pow(x, y);
}
double GetERY(double pP[2])
{
	double a1, a2, b1, b2;
	double csxy, csxz, csyz;
	a1 = tgamma(0.5 / pP[0]);
	a2 = tgamma(1.0 / pP[0]);
	b1 = tgamma(0.5 / pP[1]);
	b2 = tgamma(1.0 / pP[1]);
	csxy = a1*a1 / a2 / pP[0];
	csxz = b1*b1 / b2 / pP[1];
	csyz = b1*b1 / b2 / pP[1];
	return sqrt(csxy*csyz) / csxz;
}
double GetERZ(double pP[2])
{
	double a1, a2, b1, b2;
	double csxy, csxz, csyz;
	a1 = tgamma(0.5 / pP[0]);
	a2 = tgamma(1.0 / pP[0]);
	b1 = tgamma(0.5 / pP[1]);
	b2 = tgamma(1.0 / pP[1]);
	csxy = a1*a1 / a2 / pP[0];
	csxz = b1*b1 / b2 / pP[1];
	csyz = b1*b1 / b2 / pP[1];
	return sqrt(csxz*csyz) / csxy;
}
void GetEllipsoidity(double ppara[5])
{
	int num, numi;//剖分精度
	int i, j;
	double dv, v, u, cos1, sin1, cos2, sin2;
	double *pcos1, *psin1, *pcos2, *psin2;
	double *sgncos1, *sgnsin1, *sgnsin2;
	double rrmax, rrmin, rrr;	
	double vp, vcd, vid;
	double chong1, chong2, chong3;

	num = 20;
	dv = PI / (double)num;
	NSuP = 2 * (num*num - num + 1);
	SurfPot = new CVector[NSuP];
	pcos1 = new double[2 * num];
	psin1 = new double[2 * num];
	pcos2 = new double[num - 1];
	psin2 = new double[num - 1];
	sgncos1 = new double[2 * num];
	sgnsin1 = new double[2 * num];
	sgnsin2 = new double[num - 1];

	//v,cita (-0.5*PI,0.5*PI)
	for (i = 1; i < num; i++)
	{
		v = -0.5*PI + dv*i;
		cos2 = cos(v);
		sin2 = sin(v);
		pcos2[i - 1] = pow(fabs(cos2), 1.0 / ppara[4]);
		psin2[i - 1] = pow(fabs(sin2), 1.0 / ppara[4]);
		sgnsin2[i - 1] = sgn(sin2);
	}
	//u,fai [-PI,PI)
	for (j = 0; j < 2 * num; j++)
	{
		u = -PI + dv*j;
		cos1 = cos(u);
		sin1 = sin(u);
		pcos1[j] = pow(fabs(cos1), 1.0 / ppara[3]);
		psin1[j] = pow(fabs(sin1), 1.0 / ppara[3]);
		sgncos1[j] = sgn(cos1);
		sgnsin1[j] = sgn(sin1);
	}
	SurfPot[0][0] = 0.0;
	SurfPot[0][1] = 0.0;
	SurfPot[0][2] = -ppara[2];
	numi = 1;
	for (i = 1; i < num; i++)
		for (j = 0; j < 2 * num; j++)
		{
			SurfPot[numi][0] = ppara[0] * sgncos1[j] * pcos1[j] * pcos2[i - 1];
			SurfPot[numi][1] = ppara[1] * sgnsin1[j] * psin1[j] * pcos2[i - 1];
			SurfPot[numi][2] = ppara[2] * sgnsin2[i - 1] * psin2[i - 1];
			numi++;
		}
	SurfPot[numi][0] = 0.0;
	SurfPot[numi][1] = 0.0;
	SurfPot[numi][2] = ppara[2];

	rrmax = SurfPot[0].Length(); rrmin = rrmax;
	for (i = 1; i < NSuP; i++)
	{
		rrr = SurfPot[i].Length();
		if (rrr > rrmax) rrmax = rrr;
		if (rrr < rrmin) rrmin = rrr;
	}	

	vcd = 4.0 / 3.0*PI*rrmax*rrmax*rrmax;
	vid = 4.0 / 3.0*PI*rrmin*rrmin*rrmin;
	vp = GetVolume(ppara);

	chong1 = vp / vcd;
	chong2 = vid / vp;
	chong3 = vid / vcd;
	result << chong1 << "\t" << chong2 << "\t" << chong3 << "\t";

	delete[] pcos1; delete[] psin1; delete[] pcos2; delete[] psin2;
	delete[] sgncos1; delete[] sgnsin1; delete[] sgnsin2;
	delete[] SurfPot;
}

//求交，求配位数
void contact(double bilv)//将颗粒边长扩大至原来的bilv倍，求配位数
{
	int m, n;
	int i, j, k;
	CVector V_PBC, centerm, dcenter, V;
	double AB_dis;
	double totalnumber;
	double CCinerror, CCouterror;//变化之后的内切球外接球直径平方。	

	totalnumber = 0;
	for (m = 0; m < Npolyhedron; m++)
	{
		centerm = polyhedra[m]->center;
		pTmp[0]->Copyfrom(polyhedra[m]);
		pTmp[0]->ChangeSize(bilv);
		for (n = m; n < Npolyhedron; n++)
		{
			dcenter = polyhedra[n]->center.operator- (centerm);
			CCinerror = (polyhedra[m]->rin + polyhedra[n]->rin)*bilv - ERROR1*2.0; CCinerror = CCinerror*CCinerror;
			CCouterror = (polyhedra[m]->rout + polyhedra[n]->rout) *bilv + ERROR1*2.0; CCouterror = CCouterror*CCouterror;
			for (i = -1; i <= 1; i++)
				for (j = -1; j <= 1; j++)
					for (k = -1; k <= 1; k++)
					{
						if (m == n && i == 0 && j == 0 && k == 0)
							continue;
						V_PBC = att[0] * i + att[1] * j + att[2] * k;
						V = dcenter + V_PBC;
						AB_dis = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
						if (AB_dis < CCinerror)//两颗粒质心距离小于内接球直径，必定相交
						{
							totalnumber = totalnumber + 2;
							continue;
						}
						if (AB_dis > CCouterror)//两颗粒质心距离大于外接球直径，必定相离
							continue;

						//否则需仔细判断
						pTmp[1]->Copyfrom(polyhedra[n]);
						pTmp[1]->ChangeSize(bilv);
						pTmp[1]->Jump(V_PBC);
						if (Inspect(pTmp[0], pTmp[1]) == 1)
							totalnumber = totalnumber + 2;
					}
		}
	}
	totalnumber = totalnumber / double(Npolyhedron);
	result << totalnumber << "\t";
	//midinfor << totalnumber << "\t";
}
void ErrorAnalysis1()//单个容差分析log
{
	ofstream midinfor("midinfor.txt");
	int i;
	double er1, loger1;
	midinfor << "i\tloger1\ter1\tCN\tDOF" << endl;
	for (i = 0; i < 40; i++)
	{
		loger1 = -8.0 + 0.2 * i;
		er1 = pow(10.0, loger1);
		midinfor << i << "\t" << loger1 << "\t" << er1 << "\t";
		contact(1 + er1);
		midinfor << endl;
	}
	midinfor.close();
}
void ErrorAnalysis2()//单个容差分析
{
	ofstream midinfor("midinfor.txt");
	int i;
	double er1;
	midinfor << "i\ter1\tCN\tDOF" << endl;
	for (i = 0; i < 50; i++)
	{
		//er1 = i*0.01 + 0.01;
		er1 = i*0.001 + 0.001;
		//er1 = i*5.0 + 5.0;

		midinfor << i << "\t" << er1 << "\t";
		contact(1 + er1);
		midinfor << endl;
	}
	midinfor.close();
}
int GlobalOverlapCheck()
{
	int m, n;
	int i, j, k;
	CVector V_PBC, centerm, dcenter, V;
	double AB_dis;
	double CCinerror, CCouterror;//变化之后的内切球外接球直径平方。	

	for (m = 0; m < Npolyhedron; m++)
	{
		centerm = polyhedra[m]->center;
		pTmp[0]->Copyfrom(polyhedra[m]);
		for (n = 0; n < Npolyhedron; n++)
		{
			dcenter = polyhedra[n]->center.operator- (centerm);
			CCinerror = polyhedra[m]->rin + polyhedra[n]->rin; CCinerror = CCinerror*CCinerror*0.9;
			CCouterror = polyhedra[m]->rout + polyhedra[n]->rout; CCouterror = CCouterror*CCouterror*1.1;
			for (i = -1; i <= 1; i++)
				for (j = -1; j <= 1; j++)
					for (k = -1; k <= 1; k++)
					{
						if (m == n && i == 0 && j == 0 && k == 0)
							continue;
						V_PBC = att[0] * i + att[1] * j + att[2] * k;
						V = dcenter + V_PBC;
						AB_dis = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
						if (AB_dis < CCinerror)//两颗粒质心距离小于内接球直径，必定相交
							return 1;
						/*{
							cout << polyhedra[m]->rin << "\t" << polyhedra[n]->rin << "\t" << CCinerror << "\t" << CCouterror << "\t" << AB_dis << endl;
							SCRmn(m, n, 0, pTmp[0], pTmp[1]); getchar(); getchar();
							return 1;
						}*/
						if (AB_dis > CCouterror)//两颗粒质心距离大于外接球直径，必定相离
							continue;

						//否则需仔细判断
						pTmp[1]->Copyfrom(polyhedra[n]);
						pTmp[1]->Jump(V_PBC);
						if (Inspect(pTmp[0], pTmp[1]) == 1)
						{
							//SCRmn(m, n, 0, pTmp[0], pTmp[1]);
							return 1;
						}
						if (Inspect(pTmp[1], pTmp[0]) == 1)
						{
							//SCRmn(m, n, 1, pTmp[0], pTmp[1]);
							return 1;
						}
					}
		}
	}
	return 0;
}
int SingleOverlapCheck(int Mi)
{
	int m, n;
	int i, j, k;
	CVector V_PBC, centerm, dcenter, V;
	double AB_dis;
	double CCinerror, CCouterror;//变化之后的内切球外接球直径平方。	

	//for (m = 0; m < Npolyhedron; m++)
	m = Mi;
	{
		centerm = polyhedra[m]->center;
		pTmp[0]->Copyfrom(polyhedra[m]);
		for (n = 0; n < Npolyhedron; n++)
		{
			dcenter = polyhedra[n]->center.operator- (centerm);
			CCinerror = polyhedra[m]->rin + polyhedra[n]->rin - ERROR1*2.0; CCinerror = CCinerror*CCinerror;
			CCouterror = polyhedra[m]->rout + polyhedra[n]->rout + ERROR1*2.0; CCouterror = CCouterror*CCouterror;
			for (i = -1; i <= 1; i++)
				for (j = -1; j <= 1; j++)
					for (k = -1; k <= 1; k++)
					{
						if (m == n && i == 0 && j == 0 && k == 0)
							continue;
						V_PBC = att[0] * i + att[1] * j + att[2] * k;
						V = dcenter + V_PBC;
						AB_dis = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
						if (AB_dis < CCinerror)//两颗粒质心距离小于内接球直径，必定相交
							return n;
						if (AB_dis > CCouterror)//两颗粒质心距离大于外接球直径，必定相离
							continue;
						//否则需仔细判断
						pTmp[1]->Copyfrom(polyhedra[n]);
						pTmp[1]->Jump(V_PBC);
						if (Inspect(pTmp[0], pTmp[1]) == 1)
						{
							SCRmn(m, n, 0, pTmp[0], pTmp[1]);
							return n;
						}
						if (Inspect(pTmp[1], pTmp[0]) == 1)
						{
							SCRmn(m, n, 1, pTmp[0], pTmp[1]);
							return n;
						}
					}
		}
	}
	return -1;
}
int doubleOverlapCheck(int Mi, int Ni)
{
	int m, n;
	int ab, ba, lapmark;
	int i, j, k;
	CVector V_PBC, centerm, dcenter, V;
	double AB_dis, avalue, bvalue;
	double CCinerror, CCouterror;//变化之后的内切球外接球直径平方。	

	lapmark = 0;
	//for (m = 0; m < Npolyhedron; m++)
	m = Mi;
	{
		centerm = polyhedra[m]->center;
		pTmp[0]->Copyfrom(polyhedra[m]);
		//for (n = 0; n < Npolyhedron; n++)
		n = Ni;
		{
			dcenter = polyhedra[n]->center.operator- (centerm);
			CCinerror = polyhedra[m]->rin + polyhedra[n]->rin - ERROR1*2.0; CCinerror = CCinerror*CCinerror;
			CCouterror = polyhedra[m]->rout + polyhedra[n]->rout + ERROR1*2.0; CCouterror = CCouterror*CCouterror;
			for (i = -1; i <= 1; i++)
				for (j = -1; j <= 1; j++)
					for (k = -1; k <= 1; k++)
					{
						if (m == n && i == 0 && j == 0 && k == 0)
							continue;
						V_PBC = att[0] * i + att[1] * j + att[2] * k;

						V = dcenter + V_PBC;
						AB_dis = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
						if (AB_dis < CCinerror)//两颗粒质心距离小于内接球直径，必定相交
						{
							lapmark = 1;
							continue;
						}
						if (AB_dis > CCouterror)//两颗粒质心距离大于外接球直径，必定相离
							continue;
						//否则需仔细判断
						pTmp[1]->Copyfrom(polyhedra[n]);
						pTmp[1]->Jump(V_PBC);
						ab = 0, ba = 0; avalue = 0; bvalue = 0;
						ab = Inspect(pTmp[0], pTmp[1]);
						ba = Inspect(pTmp[1], pTmp[0]);
						//midinfor << m << "\t" << n << "\t";
						avalue = Inspectvalue(pTmp[0], pTmp[1]);
						//midinfor << m << "\t" << n << "\t";
						bvalue = Inspectvalue(pTmp[1], pTmp[0]);
						SCRmn(m, n, ab + ba, pTmp[0], pTmp[1]);
						//midinfor << m << "\t" << n << "\t" << ab << "\t" << ba << "\t" << avalue << "\t" << avalue << endl;
					}
		}
	}
	return lapmark;
}

//输出结构文件
void OutputPAC(int Counter)
{
	char cha[100];
	sprintf_s(cha, "%dpackingstructure.txt", Counter);
	ofstream pac(cha);
	pac << fixed << setprecision(16);
	
	int i, j;
	CVector V;
	
	pac << "NEW" << endl;
	pac << Counter << endl;
	pac << "PackingSpaceInformation" << endl;
	for (i = 0; i < 3; i++)
		pac << att[i][0] << "\t" << att[i][1] << "\t" << att[i][2] << endl;
	pac << "Particle&Structure" << endl;
	pac << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << PackingDensity << endl;
	pac << "Superellipsoidevery" << endl;
	pac << Npolyhedron << "\t" << NKind << endl;
	for (i = 0; i < NKind; i++)
	{
		pac << NpEach[i] << "\t" << VolfEach[i] << "\t" << SizeEach[i];
		for (j = 0; j < 5; j++)
			pac << "\t" << PPara[i][j];
		pac << endl;
	}
	for (i = 0; i < Npolyhedron; i++)
	{
		V = polyhedra[i]->center;
		pac << V[0] << "\t" << V[1] << "\t" << V[2];
		for (j = 0; j < 3; j++)
		{
			V = polyhedra[i]->e[j];
			pac << "\t" << V[0] << "\t" << V[1] << "\t" << V[2];
		}
		pac << endl;
	}
	pac << "END" << endl;
	pac << "ENDOFFILE" << endl;
	pac.close();	
}
void OutputPAC()
{
	ofstream pac;
	pac.open("totalstructurenew.txt", ios::app);

	int i, j;
	CVector V;	

	pac << fixed << setprecision(16);
	pac << "NEW" << endl;
	pac << Counterfile << endl;
	pac << "PackingSpaceInformation" << endl;
	for (i = 0; i < 3; i++)
		pac << att[i][0] << "\t" << att[i][1] << "\t" << att[i][2] << endl;
	pac << "Particle&Structure" << endl;
	pac << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << PackingDensity << endl;
	pac << "Superellipsoidevery" << endl;
	pac << Npolyhedron << "\t" << NKind << endl;
	for (i = 0; i < NKind; i++)
	{
		pac << NpEach[i] << "\t" << VolfEach[i] << "\t" << SizeEach[i];
		//pac << NpEach[i] << "\t" << VolfEach[i];
		for (j = 0; j < 5; j++)
			pac << "\t" << PPara[i][j];
		pac << endl;
	}
	for (i = 0; i < Npolyhedron; i++)
	{
		V = polyhedra[i]->center;
		pac << V[0] << "\t" << V[1] << "\t" << V[2];
		for (j = 0; j < 3; j++)
		{
			V = polyhedra[i]->e[j];
			pac << "\t" << V[0] << "\t" << V[1] << "\t" << V[2];
		}
		pac << endl;
	}
	pac << "END" << endl;
	pac << "ENDOFFILE" << endl;
	pac.close();
}
void OutputPACsample(int Counter)
{
	char cha[100];
	sprintf_s(cha, "packingsample%d.txt", Counter);
	ofstream pac(cha);
	pac << fixed << setprecision(16);

	int i, j;
	CVector V;

	for (i = 0; i < 3; i++)
		pac << att[i][0] << "\t" << att[i][1] << "\t" << att[i][2] << endl;
	pac << Npolyhedron << "\t" << NKind << endl;
	for (i = 0; i < NKind; i++)
	{
		pac << NpEach[i] << "\t" << VolfEach[i] << "\t" << SizeEach[i];
		for (j = 0; j < 5; j++)
			pac << "\t" << PPara[i][j];
		pac << endl;
	}
	for (i = 0; i < Npolyhedron; i++)
	{
		V = polyhedra[i]->center;
		pac << V[0] << "\t" << V[1] << "\t" << V[2];
		for (j = 0; j < 3; j++)
		{
			V = polyhedra[i]->e[j];
			pac << "\t" << V[0] << "\t" << V[1] << "\t" << V[2];
		}
		pac << endl;
	}
	pac.close();
}

//输出SCR文件
void SCRmn(int mp, int np)
{
	char cha[100];
	sprintf_s(cha, "%d %d %d.scr", Counterfile, mp, np);
	ofstream scr(cha);
	scr << fixed << setprecision(15);

	int mnp[2];
	int num = 20;	
	CVector ce, xt, xtt;
	CVector *points;
	double pP[2], pA[3];
	double A[3][3];
	double v, u, cos1, sin1, cos2, sin2;
	double x, y, z;
	double dv = PI / (double)num;	
	
	points = new CVector[num + 1];
	mnp[0] = mp; mnp[1] = np;
	for (int k = 0; k < 2; k++)
	{
		pP[0] = polyhedra[mnp[k]]->pP[0];
		pP[1] = polyhedra[mnp[k]]->pP[1];
		pA[0] = polyhedra[mnp[k]]->pA[0];
		pA[1] = polyhedra[mnp[k]]->pA[1];
		pA[2] = polyhedra[mnp[k]]->pA[2];
		ce = polyhedra[mnp[k]]->center;
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++)
				A[jj][ii] = polyhedra[mnp[k]]->e[ii][jj];
		//3dmesh in AutoCAD
		scr << "3dmesh" << endl;
		int m = 2 * num + 1;
		int n = num + 1;
		scr << m << endl;
		scr << n << endl;
		for (int j = 0; j < 2 * num; j++)
		{
			u = -PI + dv*j;
			cos1 = cos(u);
			sin1 = sin(u);

			for (int i = 0; i <= num; i++)
			{
				v = -0.5*PI + dv*i;
				cos2 = cos(v);
				sin2 = sin(v);

				x = pA[0] * sgn(cos1)*pow(fabs(cos1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
				y = pA[1] * sgn(sin1)*pow(fabs(sin1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
				z = pA[2] * sgn(sin2)*pow(fabs(sin2), 1.0 / pP[1]);

				xt.Set(x, y, z);
				MatrixMult(A, xt, xtt);
				xtt = xtt + ce;
				scr << xtt[0] << "," << xtt[1] << "," << xtt[2] << endl;
				if (j == 0) points[i].Set(xtt);
			}
		}
		for (int i = 0; i <= num; i++)
			scr << points[i][0] << "," << points[i][1] << "," << points[i][2] << endl;
	}

	/*CVector V, vv[3];
	double rrr = Pidmin*0.05;
	for (int j = 0; j < 2; j++)
		for (int k = 0; k < 2; k++)
		{
			vv[0] = att[1] * j + att[2] * k;
			vv[1] = att[0] * j + att[2] * k;
			vv[2] = att[0] * j + att[1] * k;

			for (int kk = 0; kk < 3; kk++)
			{
				V = att[kk] + vv[kk];
				scr << "cylinder " << vv[kk][0] << "," << vv[kk][1] << "," << vv[kk][2] << " " << rrr << " a " << V[0] << "," << V[1] << "," << V[2] << endl;
			}
		}*/
	scr << "grid off\nVSCURRENT c\nzoom e\n";
	scr.close();
	delete[] points;	
}
void SCRmn(int mp, int np, int aabb, CPolyhedron * psp0, CPolyhedron * psp1)
{
	char cha[100];
	sprintf_s(cha, "%d,%d,%d,%d.scr", Counterfile, mp, np, aabb);
	ofstream scr(cha);
	scr << fixed << setprecision(15);
	
	CPolyhedron* mnp[2];
	int num = 20;	
	CVector ce, xt, xtt;
	CVector *points;
	double pP[2], pA[3];
	double A[3][3];
	double v, u, cos1, sin1, cos2, sin2;
	double x, y, z;
	double dv = PI / (double)num;

	points = new CVector[num + 1];
	mnp[0] = psp0; mnp[1] = psp1;
	for (int k = 0; k < 2; k++)
	{
		pP[0] = mnp[k]->pP[0];
		pP[1] = mnp[k]->pP[1];
		pA[0] = mnp[k]->pA[0];
		pA[1] = mnp[k]->pA[1];
		pA[2] = mnp[k]->pA[2];
		ce = mnp[k]->center;
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++)
				A[jj][ii] = mnp[k]->e[ii][jj];
		//3dmesh in AutoCAD
		scr << "3dmesh" << endl;
		int m = 2 * num + 1;
		int n = num + 1;
		scr << m << endl;
		scr << n << endl;
		for (int j = 0; j < 2 * num; j++)
		{
			u = -PI + dv*j;
			cos1 = cos(u);
			sin1 = sin(u);

			for (int i = 0; i <= num; i++)
			{
				v = -0.5*PI + dv*i;
				cos2 = cos(v);
				sin2 = sin(v);

				x = pA[0] * sgn(cos1)*pow(fabs(cos1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
				y = pA[1] * sgn(sin1)*pow(fabs(sin1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
				z = pA[2] * sgn(sin2)*pow(fabs(sin2), 1.0 / pP[1]);

				xt.Set(x, y, z);
				MatrixMult(A, xt, xtt);
				xtt = xtt + ce;
				scr << xtt[0] << "," << xtt[1] << "," << xtt[2] << endl;
				if (j == 0) points[i].Set(xtt);
			}
		}
		for (int i = 0; i <= num; i++)
			scr << points[i][0] << "," << points[i][1] << "," << points[i][2] << endl;
	}

	/*CVector V, vv[3];
	double rrr = Pidmin*0.05;
	for (int j = 0; j < 2; j++)
		for (int k = 0; k < 2; k++)
		{
			vv[0] = att[1] * j + att[2] * k;
			vv[1] = att[0] * j + att[2] * k;
			vv[2] = att[0] * j + att[1] * k;

			for (int kk = 0; kk < 3; kk++)
			{
				V = att[kk] + vv[kk];
				scr << "cylinder " << vv[kk][0] << "," << vv[kk][1] << "," << vv[kk][2] << " " << rrr << " a " << V[0] << "," << V[1] << "," << V[2] << endl;
			}
		}*/
	scr << "grid off\nVSCURRENT c\nzoom e\n";
	scr.close();
	delete[] points;	
}
void SCR()
{
	char cha[100];
	sprintf_s(cha, "%d.scr", Counterfile);
	ofstream scr(cha);
	scr << fixed << setprecision(15);

	int num = 20;
	CVector ce, xt, xtt;
	CVector *points;
	double pP[2], pA[3];
	double A[3][3];
	double v, u, cos1, sin1, cos2, sin2;
	double x, y, z;
	double dv = PI / (double)num;

	points = new CVector[num + 1];
	for (int k = 0; k < Npolyhedron; k++)
	{
		pP[0] = polyhedra[k]->pP[0];
		pP[1] = polyhedra[k]->pP[1];
		pA[0] = polyhedra[k]->pA[0];
		pA[1] = polyhedra[k]->pA[1];
		pA[2] = polyhedra[k]->pA[2];
		ce = polyhedra[k]->center;
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++)
				A[jj][ii] = polyhedra[k]->e[ii][jj];
		//3dmesh in AutoCAD
		scr << "3dmesh" << endl;
		int m = 2 * num + 1;
		int n = num + 1;
		scr << m << endl;
		scr << n << endl;
		for (int j = 0; j < 2 * num; j++)
		{
			u = -PI + dv*j;
			cos1 = cos(u);
			sin1 = sin(u);

			for (int i = 0; i <= num; i++)
			{
				v = -0.5*PI + dv*i;
				cos2 = cos(v);
				sin2 = sin(v);

				x = pA[0] * sgn(cos1)*pow(fabs(cos1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
				y = pA[1] * sgn(sin1)*pow(fabs(sin1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
				z = pA[2] * sgn(sin2)*pow(fabs(sin2), 1.0 / pP[1]);

				xt.Set(x, y, z);
				MatrixMult(A, xt, xtt);
				xtt = xtt + ce;
				scr << xtt[0] << "," << xtt[1] << "," << xtt[2] << endl;
				if (j == 0) points[i].Set(xtt);
			}
		}
		for (int i = 0; i <= num; i++)
			scr << points[i][0] << "," << points[i][1] << "," << points[i][2] << endl;
	}

	CVector V, vv[3];
	double rrr = Pidmin*0.05;
	for (int j = 0; j < 2; j++)
		for (int k = 0; k < 2; k++)
		{
			vv[0] = att[1] * j + att[2] * k;
			vv[1] = att[0] * j + att[2] * k;
			vv[2] = att[0] * j + att[1] * k;

			for (int kk = 0; kk < 3; kk++)
			{
				V = att[kk] + vv[kk];
				scr << "cylinder " << vv[kk][0] << "," << vv[kk][1] << "," << vv[kk][2] << " " << rrr << " a " << V[0] << "," << V[1] << "," << V[2] << endl;
			}
		}
	scr << "grid off\nVSCURRENT c\nzoom e\n";
	scr.close();
	delete[] points;
}
void SCRcenter()
{
	char cha[100];
	sprintf_s(cha, "%d SPNum_%d center.scr", Counterfile, Npolyhedron);
	ofstream scr(cha);
	scr << fixed << setprecision(15);
	
	int m;
	CVector  v[2];
	double radius = Pidmin*1.0;	

	for (m = 0; m < Npolyhedron; m++)
	{
		if (radius > 0)
		{
			v[0] = polyhedra[m]->center;
			//if (v[0][0]<0.15*a[0][0])
			//scr << "sphere " << v[0][0] << "," << v[0][1] << "," << v[0][2] << " " << radius << endl;
			scr << "sphere " << v[0][0] << "," << v[0][1] << "," << v[0][2] << " " << polyhedra[m]->pA[0] << endl;
		}
	}
	scr << "grid off\nVSCURRENT c\nzoom e\n";
	scr.close();	
}
void SCRmodel(double ppara[5])
{
	ofstream scr("model.scr");
	scr << fixed << setprecision(15);

	int num = 30;
	CVector ce, xt, xtt;
	CVector *points;
	double pP[2], pA[3];
	double A[3][3];
	double v, u, cos1, sin1, cos2, sin2;
	double x, y, z;
	double dv = PI / (double)num;

	points = new CVector[num + 1];
	{
		pA[0] = ppara[0];
		pA[1] = ppara[1];
		pA[2] = ppara[2];
		pP[0] = ppara[3];
		pP[1] = ppara[4];
		ce.SetZero();
		for (int ii = 0; ii < 3; ii++)
		{
			for (int jj = 0; jj < 3; jj++)
				A[ii][jj] = 0.0;
			A[ii][ii] = 1.0;
		}
		//3dmesh in AutoCAD
		scr << "3dmesh" << endl;
		int m = 2 * num + 1;
		int n = num + 1;
		scr << m << endl;
		scr << n << endl;
		for (int j = 0; j < 2 * num; j++)
		{
			u = -PI + dv*j;
			cos1 = cos(u);
			sin1 = sin(u);

			for (int i = 0; i <= num; i++)
			{
				v = -0.5*PI + dv*i;
				cos2 = cos(v);
				sin2 = sin(v);

				x = pA[0] * sgn(cos1)*pow(fabs(cos1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
				y = pA[1] * sgn(sin1)*pow(fabs(sin1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
				z = pA[2] * sgn(sin2)*pow(fabs(sin2), 1.0 / pP[1]);

				xt.Set(x, y, z);
				MatrixMult(A, xt, xtt);
				xtt = xtt + ce;
				scr << xtt[0] << "," << xtt[1] << "," << xtt[2] << endl;
				if (j == 0) points[i].Set(xtt);
			}
		}
		for (int i = 0; i <= num; i++)
			scr << points[i][0] << "," << points[i][1] << "," << points[i][2] << endl;
	}

	scr << "grid off\nVSCURRENT c\nzoom e\n";
	scr.close();
	delete[] points;
}
void povrayrandomcolor()//颗粒
{
	char cha[100];
	sprintf_s(cha, "%d povfile_superball.pov", Counterfile);
	ofstream scr(cha);
	scr << fixed << setprecision(15);

	int i, j, k, kk, m;
	CVector  V, v[4];
	int mark;
	double SBmax;
	double rrr = 0.005;
	double rgbb[3], rgb[8][3];

	SBmax = 0.0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			if (att[i][j]>SBmax)
				SBmax = att[i][j];
		}
	rgb[0][0] = 0.000, rgb[0][1] = 0.353, rgb[0][2] = 0.671;
	rgb[1][0] = 0.251, rgb[1][1] = 0.455, rgb[1][2] = 0.204;
	rgb[2][0] = 0.263, rgb[2][1] = 0.231, rgb[2][2] = 0.271;
	rgb[3][0] = 0.282, rgb[3][1] = 0.608, rgb[3][2] = 0.263;
	rgb[4][0] = 0.286, rgb[4][1] = 0.294, rgb[4][2] = 0.588;
	rgb[5][0] = 0.620, rgb[5][1] = 0.294, rgb[5][2] = 0.278;
	rgb[6][0] = 0.773, rgb[6][1] = 0.122, rgb[6][2] = 0.122;
	rgb[7][0] = 0.863, rgb[7][1] = 0.341, rgb[7][2] = 0.071;
	
	scr << "#include \"colors.inc\"" << endl;
	scr << "#include \"textures.inc\"" << endl;
	scr << "camera{ location <8, 6, 8>\nright 0.15*x*image_width / image_height\nup 0.15*x\nlook_at <0, 0, 0>}" << endl;
	scr << "background {White}" << endl;
	scr << "light_source{ <8,20,30> color rgb <0.77,0.75,0.75> }" << endl;
	scr << "light_source{ <25,12,12> color rgb <0.63,0.65,0.65> }" << endl;
	scr << "light_source{ <12,25,12> color rgb <0.63,0.65,0.65> }" << endl;

	scr << "union {" << endl;
	for (m = 0; m < Npolyhedron; m++)
		//for (m = 0; m < 1; m++)
	{
		mark = rand() % 8;
		for (kk = 0; kk < 3; kk++)
			rgbb[kk] = rgb[mark][kk];
		//mark = 0;
		scr << "object{superellipsoid{<" << 1.0 / polyhedra[m]->pP[0] << "," << 1.0 / polyhedra[m]->pP[1] << ">";
		scr << "scale<" << polyhedra[m]->pA[0] / SBmax << "," << polyhedra[m]->pA[1] / SBmax << "," << polyhedra[m]->pA[2] / SBmax << ">";
		scr << "texture{pigment{color rgb<" << rgbb[0] << "," << rgbb[1] << "," << rgbb[2] << ">}finish{reflection 0.06 specular 0.13}}}";
		scr << "matrix<";
		for (int j = 0; j < 3; j++)
			scr << polyhedra[m]->e[j][0] << "," << polyhedra[m]->e[j][1] << "," << polyhedra[m]->e[j][2] << ",";
		scr << polyhedra[m]->center[0] / SBmax << "," << polyhedra[m]->center[1] / SBmax << "," << polyhedra[m]->center[2] / SBmax << ">}" << endl;
	}
	rgbb[0] = 0.0, rgbb[1] = 0.0, rgbb[2] = 1.0;
	for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++)
		{
			v[0] = att[1] * j + att[2] * k;
			v[1] = att[0] * j + att[2] * k;
			v[2] = att[0] * j + att[1] * k;

			for (kk = 0; kk < 3; kk++)
			{
				V = att[kk] + v[kk];
				scr << "object{cylinder{"
					<< "<" << v[kk][0] / SBmax << "," << v[kk][1] / SBmax << "," << v[kk][2] / SBmax << ">,"
					<< "<" << V[0] / SBmax << "," << V[1] / SBmax << "," << V[2] / SBmax << ">,"
					<< rrr
					<< " pigment{rgb<" << rgbb[0] << "," << rgbb[1] << "," << rgbb[2] << ">" << "}"
					<< "finish{reflection 0.06 specular 0.13}}}"
					<< endl;
			}
		}
	scr << "translate <0.2,-0.1,-0.2>\nrotate -22.5*y\n}" << endl;
	scr.close();	
}
void povray()//颗粒
{
	char cha[100];
	sprintf_s(cha, "%d povfile_superball.pov", Counterfile);
	ofstream scr(cha);
	scr << fixed << setprecision(15);

	int i, j, k, kk, m;
	CVector  V, v[4];
	int mark;
	double SBmax;
	double rrr = 0.005;
	double rgbb[3], rgb[9][3];

	SBmax = 0.0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			if (att[i][j]>SBmax)
				SBmax = att[i][j];
		}

	rgb[0][0] = 0.00000, rgb[0][1] = 0.50196, rgb[0][2] = 0.00000;
	rgb[1][0] = 1.00000, rgb[1][1] = 1.00000, rgb[1][2] = 0.00000;
	rgb[2][0] = 1.00000, rgb[2][1] = 0.38824, rgb[2][2] = 0.27843;
	rgb[3][0] = 0.00000, rgb[3][1] = 0.00000, rgb[3][2] = 1.00000;
	rgb[4][0] = 1.00000, rgb[4][1] = 0.00000, rgb[4][2] = 1.00000;
	rgb[5][0] = 0.60392, rgb[5][1] = 0.80392, rgb[5][2] = 0.19608;
	rgb[6][0] = 0.50196, rgb[6][1] = 0.00000, rgb[6][2] = 0.50196;
	rgb[7][0] = 0.50196, rgb[7][1] = 0.50196, rgb[7][2] = 0.00000;
	rgb[8][0] = 0.60000, rgb[8][1] = 0.19608, rgb[8][2] = 0.80000;

	rgb[0][0] = 0.00000, rgb[0][1] = 0.00000, rgb[0][2] = 1.00000;
	rgb[1][0] = 1.00000, rgb[1][1] = 1.00000, rgb[1][2] = 0.00000;

	scr << "#include \"colors.inc\"" << endl;
	scr << "#include \"textures.inc\"" << endl;
	scr << "camera{ location <8, 6, 8>\nright 0.15*x*image_width / image_height\nup 0.15*x\nlook_at <0, 0, 0>}" << endl;
	scr << "background {White}" << endl;
	scr << "light_source{ <8,20,30> color rgb <0.77,0.75,0.75> }" << endl;
	scr << "light_source{ <25,12,12> color rgb <0.63,0.65,0.65> }" << endl;
	scr << "light_source{ <12,25,12> color rgb <0.63,0.65,0.65> }" << endl;

	scr << "union {" << endl;
	for (m = 0; m < Npolyhedron; m++)
		//for (m = 0; m < 1; m++)
	{
		mark = (polyhedra[m]->kind) % 9;
		for (kk = 0; kk < 3; kk++)
			rgbb[kk] = rgb[mark][kk];
		//mark = 0;
		scr << "object{superellipsoid{<" << 1.0 / polyhedra[m]->pP[0] << "," << 1.0 / polyhedra[m]->pP[1] << ">";
		scr << "scale<" << polyhedra[m]->pA[0] / SBmax << "," << polyhedra[m]->pA[1] / SBmax << "," << polyhedra[m]->pA[2] / SBmax << ">";
		scr << "texture{pigment{color rgb<" << rgbb[0] << "," << rgbb[1] << "," << rgbb[2] << ">}finish{reflection 0.1 specular 0.13}}}";
		scr << "matrix<";
		for (int j = 0; j < 3; j++)
			scr << polyhedra[m]->e[j][0] << "," << polyhedra[m]->e[j][1] << "," << polyhedra[m]->e[j][2] << ",";
		scr << polyhedra[m]->center[0] / SBmax << "," << polyhedra[m]->center[1] / SBmax << "," << polyhedra[m]->center[2] / SBmax << ">}" << endl;
	}
	rgbb[0] = 0.0, rgbb[1] = 0.0, rgbb[2] = 0.0;
	for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++)
		{
			v[0] = att[1] * j + att[2] * k;
			v[1] = att[0] * j + att[2] * k;
			v[2] = att[0] * j + att[1] * k;

			for (kk = 0; kk < 3; kk++)
			{
				V = att[kk] + v[kk];
				scr << "object{cylinder{"
					<< "<" << v[kk][0] / SBmax << "," << v[kk][1] / SBmax << "," << v[kk][2] / SBmax << ">,"
					<< "<" << V[0] / SBmax << "," << V[1] / SBmax << "," << V[2] / SBmax << ">,"
					<< rrr
					<< " pigment{rgb<" << rgbb[0] << "," << rgbb[1] << "," << rgbb[2] << ">" << "}"
					<< "finish{reflection 0.06 specular 0.13}}}"
					<< endl;
			}
		}
	scr << "translate <0.2,-0.1,-0.2>\nrotate -22.5*y\n}" << endl;
	scr.close();
}

//位置随机性:均匀性监测,RDF
void IsotropyTest()//均匀性监测
{
	//cout<<"颗粒分布的各向同性分析...";
	if (Npolyhedron < 5)
	{
		cout << "颗粒过少，本程序拒绝进行各向同性分析." << endl;
		return;
	}
	double * Dis_x, *Dis_y, *Dis_z;
	int N = 10;
	int i;
	Dis_x = new double[N];
	Dis_y = new double[N];
	Dis_z = new double[N];
	double d = 1.0 / double(N);
	int ix, iy, iz;
	for (int i = 0; i < N; i++) { Dis_x[i] = 0.0, Dis_y[i] = 0.0, Dis_z[i] = 0.0; }

	for (i = 0; i < Npolyhedron; i++)
	{
		ix = int(polyhedra[i]->center_L[0] / d);
		iy = int(polyhedra[i]->center_L[1] / d);
		iz = int(polyhedra[i]->center_L[2] / d);
		if (ix < 0) ix = 0;	else if (ix >= N) ix = N - 1;
		if (iy < 0) iy = 0;	else if (iy >= N) iy = N - 1;
		if (iz < 0) iz = 0;	else if (iz >= N) iz = N - 1;
		Dis_x[ix] += 1.0;
		Dis_y[iy] += 1.0;
		Dis_z[iz] += 1.0;
	}
	d = double(Npolyhedron) / double(N);
	for (i = 0; i < N; i++) { Dis_x[i] /= d, Dis_y[i] /= d, Dis_z[i] /= d; }

	char cha[100];
	sprintf_s(cha, "%d ParticleDistributionOn_x_y_z.txt", Counterfile);
	ofstream it(cha);
	it << fixed << setprecision(8);	
	it << "i\tx\ty\tz\n";
	for (i = 0; i < N; i++)
		it << i + 1 << "\t" << Dis_x[i] << "\t" << Dis_y[i] << "\t" << Dis_z[i] << endl;
	it.close();
	delete[] Dis_x;
	delete[] Dis_y;
	delete[] Dis_z;	
	//cout<<"完成."<<endl;
}
void RadialDistributionFunction()//RDF，径向分布函数
{
	int i, j, k, imax, imin, p1, p2;
	double Dmin, Dmax, dR, V, PRInradius, t;
	int N, n, nmax;
	double * RDF;
	double * x;
	CVector P2;
	double a0, b0, dx, dy, dz;

	Dmin = Pidmin;
	N = 10;
	n = 20;
	dR = Dmin / double(n);
	nmax = n*N + 1;
	RDF = new double[nmax];
	x = new double[nmax];
	Dmax = Dmin*double(N);

	a0 = 0;
	for (i = 0; i < 3; i++)
	{
		b0 = att[i % 3].Cross(att[(i + 1) % 3]).Length();
		if (b0 > a0) a0 = b0;
	}
	V = fabs(att[0].Cross(att[1])*att[2]);
	PRInradius = V / a0;//填充区域内接球直径	
	imax = int(Dmax / PRInradius) + 2;
	imin = -imax + 1;
	Dmax *= Dmax;
	for (i = 0; i < nmax; i++) RDF[i] = 0.0;
	for (p1 = 0; p1 < Npolyhedron; p1++)
	{
		for (p2 = p1 + 1; p2 < Npolyhedron; p2++)
		{
			for (i = imin; i < imax; i++)
			{
				for (j = imin; j < imax; j++)
				{
					for (k = imin; k < imax; k++)
					{
						P2 = polyhedra[p2]->center + att[0] * double(i) + att[1] * double(j) + att[2] * double(k);
						dx = P2[0] - polyhedra[p1]->center[0];
						dy = P2[1] - polyhedra[p1]->center[1];
						dz = P2[2] - polyhedra[p1]->center[2];
						a0 = dx*dx + dy*dy + dz*dz;
						if (a0 < Dmax)
						{
							a0 = sqrt(a0);
							n = int(a0 / dR);
							if (n<nmax&&n>-1) RDF[n] += 1.0;
						}
					}
				}
			}
		}
	}
	a0 = 4.0*PI*dR*double(Npolyhedron);
	a0 *= double(Npolyhedron) / fabs(att[0].Cross(att[1])*att[2]);
	RDF[0] = 0.0;
	x[0] = 0.0;
	for (i = 1; i < nmax; i++)
	{
		x[i] = double(i)*dR;
		RDF[i] /= x[i] * x[i] * a0;
		RDF[i] *= 2.0;
	}
	t = 0.0;//评价RDF的波动程度
	a0 = pow(V / double(Npolyhedron), 1.0 / 3.0);
	a0 = 3.5*a0;
	j = int(a0 / dR) + 2;
	if (j < nmax)
	{
		for (k = 0; k<nmax; k++) if (RDF[k]>2.0E-6) break;
		b0 = double(k)*dR;
		for (i = k; i < j; i++) t += (fabs(RDF[i] - 1.0) *dR);
		if (a0 > b0) t /= (a0 - b0);
		else t = -1.0;
	}
	else t = -1.0;

	//文件输出
	char cha[100];
	sprintf_s(cha, "%d RDFfile.txt", Counterfile);
	ofstream rdf(cha);
	rdf << fixed << setprecision(8);	
	rdf << "x\tx/Dmin\tRDF\t";
	//rdf << Dmin;
	rdf << endl;
	for (i = 0; i < nmax; i++)
		rdf << x[i] << "\t" << x[i] / Dmin << "\t" << RDF[i] << endl;
	//rdf<<"\nT : \t"<<t<<endl;
	delete[]x;
	delete[]RDF;
	rdf.close();	
}
void OutputPairDis()//质心距离输出
{
	char cha[100];
	sprintf_s(cha, "%d PairDisfile.txt", Counterfile);
	ofstream rdf(cha);
	rdf << fixed << setprecision(8);

	int m, n;
	int i, j, k;
	CVector V_PBC, centerm, dcenter, V;
	double AB_dis;
	
	for (m = 0; m < Npolyhedron; m++)
	{
		centerm = polyhedra[m]->center;
		pTmp[0]->Copyfrom(polyhedra[m]);
		for (n = 0; n < Npolyhedron; n++)
		{
			dcenter = polyhedra[n]->center.operator- (centerm);
			for (i = -1; i <= 1; i++)
				for (j = -1; j <= 1; j++)
					for (k = -1; k <= 1; k++)
					{
						if (m == n)
							//if (m == n && i == 0 && j == 0 && k == 0)
							continue;
						V_PBC = att[0] * i + att[1] * j + att[2] * k;

						V = dcenter + V_PBC;
						AB_dis = sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
						rdf << AB_dis << endl;
					}
		}
	}
	rdf.close();	
}

//周围颗粒建系
void getppd()//particle particle deitence
{
	int i, j, k, m, n;
	double dmin, dd;
	CVector PBC, center0, dcenter, V, Vdmin;

	sblengthmin = att[0][0] * att[0][0] + att[0][1] * att[0][1] + att[0][2] * att[0][2];
	for (i = 1; i < 3; i++)
	{
		dmin = att[i][0] * att[i][0] + att[i][1] * att[i][1] + att[i][2] * att[i][2];
		if (dmin < sblengthmin)
			sblengthmin = dmin;
	}
	DISMAX = sblengthmin*10000.0;
	for (m = 0; m < Npolyhedron; m++)
	{
		pairdis[m][m] = sblengthmin / 4.0;
		pairnumber[m][m] = -1;
		pairdcenter[m][m].Set(0.0, 0.0, 0.0);
		center0 = polyhedra[m]->center;
		for (n = m + 1; n < Npolyhedron; n++)
		{
			dmin = DISMAX;
			dcenter = polyhedra[n]->center.operator- (center0);
			for (i = -1; i < 2; i++)
				for (j = -1; j < 2; j++)
					for (k = -1; k < 2; k++)
					{
						PBC = att[0] * i + att[1] * j + att[2] * k;
						V = dcenter + PBC;
						dd = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
						if (dd < dmin)
						{
							dmin = dd;
							Vdmin = V;
						}

					}
			pairdis[m][n] = dmin;
			pairdis[n][m] = dmin;
			pairnumber[m][n] = -1;
			pairnumber[n][m] = -1;
			pairdcenter[m][n] = Vdmin;
			pairdcenter[n][m] = Vdmin*(-1.0);
		}
	}

	for (m = 0; m < Npolyhedron; m++)
	{
		for (i = 0; i < Maxlocal; i++)
			getpairdismin(m, i);
	}
}
void updateppd(int m)
{
	int i, j, k, n;
	double dmin, dd;
	CVector PBC, center0, dcenter, V, Vdmin;

	int t, midnumber, tempnumber;
	double middis, tempdis;

	pairdis[m][m] = sblengthmin / 2.0;
	pairdcenter[m][m].Set(0.0, 0.0, 0.0);
	//if (pairnumber[m][m] >= 0){ cout << "error!! 颗粒本身统计" << endl; getchar(); }
	pairnumber[m][m] = -1;
	center0 = polyhedra[m]->center;
	for (n = 0; n < Npolyhedron; n++)
	{
		if (n == m)
		{
			pairdis[m][m] = sblengthmin / 2.0;
			pairnumber[m][m] = -1;
			continue;
		}
		dmin = DISMAX;
		dcenter = polyhedra[n]->center.operator- (center0);
		for (i = -1; i < 2; i++)
			for (j = -1; j < 2; j++)
				for (k = -1; k < 2; k++)
				{
					PBC = att[0] * i + att[1] * j + att[2] * k;
					V = dcenter + PBC;
					dd = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
					if (dd < dmin)
					{
						dmin = dd;
						Vdmin = V;
					}
				}
		pairdcenter[m][n] = Vdmin;
		pairdcenter[n][m] = Vdmin*(-1.0);
		pairdis[m][n] = dmin;
		pairnumber[m][n] = -1;

		pairdis[n][m] = dmin;

		//方案一：全都更新
		//if ((pairnumber[n][m] >= 0) || (dmin <= pairdismin[n][Maxlocal - 1]))//说明：原来的m已经在n的min序列当中或者原来未在序列中，现在要在序列中，此时需要更新min
		//{
		//	for (i = 0; i < Npolyhedron; i++) pairnumber[n][i] = -1;
		//	for (i = 0; i < Maxlocal; i++)
		//		getpairdismin(n, i);
		//}

		//方案二：局部更新
		if (pairnumber[n][m] >= 0)//说明：原来的m已经在n的min序列当中，需要更新min序列
		{
			t = pairnumber[n][m];
			pairdismin[n][t] = dmin;//更新此时t位置的min值

			if ((t != 0) && (t != (Maxlocal - 1)) && (dmin <= pairdismin[n][t + 1]) && (dmin >= pairdismin[n][t - 1]))
				continue;
			if ((t == 0) && (dmin <= pairdismin[n][t + 1]))
				continue;
			if ((t == (Maxlocal - 1)) && (dmin >= pairdismin[n][t - 1]))
			{
				pairnumber[n][m] = -1;
				getpairdismin(n, Maxlocal - 1);
				continue;
			}

			for (i = 0; i < t; i++)
			{
				if (pairdismin[n][i]>dmin)//距离：t之前的i比t还大，从i开始往后更新,直到t
				{
					middis = dmin;
					midnumber = m;
					for (j = i; j < t + 1; j++)
					{
						pairnumber[n][midnumber] = j;

						tempdis = pairdismin[n][j];
						pairdismin[n][j] = middis;
						middis = tempdis;

						tempnumber = pairnumbermin[n][j];
						pairnumbermin[n][j] = midnumber;
						midnumber = tempnumber;
					}
					break;
				}
			}

			if (i == t)//此时说明：t变化之后依然比t之前小，t之前无需变化。或者t=0；
			{
				for (i = Maxlocal - 1; i > t; i--)
				{
					if (pairdismin[n][i] < dmin)//距离：t之后的i比t还小，从i开始往前更新, 直到t
					{
						middis = dmin;
						midnumber = m;
						for (j = i; j > t - 1; j--)
						{
							pairnumber[n][midnumber] = j;

							tempdis = pairdismin[n][j];
							pairdismin[n][j] = middis;
							middis = tempdis;

							tempnumber = pairnumbermin[n][j];
							pairnumbermin[n][j] = midnumber;
							midnumber = tempnumber;
						}
						break;
					}
				}
			}

			if (pairnumbermin[n][Maxlocal - 1] == m)//t变化之后在最后，全局更新最后一个
			{
				pairnumber[n][m] = -1;
				getpairdismin(n, Maxlocal - 1);
			}
			continue;
		}

		else if (dmin <= pairdismin[n][Maxlocal - 1])//原来未在min序列当中，现在要在min序列当中；此时需要更新min
		{
			for (i = 0; i < Maxlocal; i++)
			{
				if (pairdismin[n][i]>dmin)//距离：i比dmin还大，从i开始往后更新,直到Maxlocal-1
				{
					middis = dmin;
					midnumber = m;
					for (j = i; j < Maxlocal; j++)
					{
						pairnumber[n][midnumber] = j;

						tempdis = pairdismin[n][j];
						pairdismin[n][j] = middis;
						middis = tempdis;

						tempnumber = pairnumbermin[n][j];
						pairnumbermin[n][j] = midnumber;
						midnumber = tempnumber;
					}
					pairnumber[n][midnumber] = -1;//此时old版本的最后一个的颗粒编号已经记为midnumber已经排除在最小序列中，类型置为-1
					break;
				}
			}
		}
	}//end of looping all particle;			
	for (i = 0; i < Maxlocal; i++)
		getpairdismin(m, i);
}
void getpairdismin(int M, int I)
{
	int j, markn;
	double minvalue;
	minvalue = DISMAX; markn = -1;
	for (j = 0; j < Npolyhedron; j++)
	{
		if (pairdis[M][j] < minvalue)
		{
			if (pairnumber[M][j] < 0)
			{
				minvalue = pairdis[M][j];
				markn = j;
			}
		}
	}
	pairnumber[M][markn] = I;
	pairdismin[M][I] = minvalue;
	pairnumbermin[M][I] = markn;
}
void outinfor(int Numm)//信息输出检测
{
	char cha[100];
	sprintf_s(cha, "%d ppdinformation.txt", Numm);
	ofstream pac(cha);
	pac << fixed << setprecision(16);

	int m, n;
	pac << "\t";
	for (m = 0; m < Npolyhedron; m++)
		pac << m << "\t";
	pac << endl;
	for (m = 0; m < Npolyhedron; m++)
	{
		pac << m << "\t";
		for (n = 0; n < Npolyhedron; n++)
			pac << pairdis[m][n] << "\t";
		pac << endl;
	}
	pac << "*************************************************************\n";
	pac << "\t";
	for (m = 0; m < Npolyhedron; m++)
		pac << m << "\t";
	pac << endl;
	for (m = 0; m < Npolyhedron; m++)
	{
		pac << m << "\t";
		for (n = 0; n < Npolyhedron; n++)
			pac << pairnumber[m][n] << "\t";
		pac << endl;
	}

	pac << "*************************************************************\n";
	pac << "\t";
	for (m = 0; m < Npolyhedron; m++)
		pac << m << "\t";
	pac << endl;
	for (m = 0; m < Maxlocal; m++)
	{
		pac << m << "\t";
		for (n = 0; n < Npolyhedron; n++)
			pac << pairdismin[n][m] << "\t";
		pac << endl;
	}
	pac << "*************************************************************\n";
	pac << "\t";
	for (m = 0; m < Npolyhedron; m++)
		pac << m << "\t";
	pac << endl;
	for (m = 0; m < Maxlocal; m++)
	{
		pac << m << "\t";
		for (n = 0; n < Npolyhedron; n++)
			pac << pairnumbermin[n][m] << "\t";
		pac << endl;
	}	
	pac.close();
}

//有序参数
double CalculateS4()
{
	int i, j, ii, jj;
	double p, p4, coscita, c;
	CVector trialvector;

	p4 = 0.0;
	for (i = 0; i < Npolyhedron; i++)
		for (j = 0; j < 3; j++)
		{
			trialvector = polyhedra[i]->e[j];
			p = 0.0;
			for (ii = 0; ii < Npolyhedron; ii++)
				for (jj = 0; jj < 3; jj++)
				{
					coscita = trialvector* polyhedra[ii]->e[jj];
					c = coscita*coscita;
					p = p + (35.0*c*c - 30.0*c + 3.0) / 14.0;
				}
			p = p / double(Npolyhedron);
			if (p > p4) p4 = p;
		}
	//p4 = fabs((p4 - S4random) / S4length);
	result << p4 << "\t";
	return p4;
}
double CalculateQ6local()
{
	int i, j, ii, jj, m;
	double t, q6[Maxlocal], Q6ev[Maxlocal], Q6av;
	double x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6, z1, z2, z4, r2, r4, r6ni;
	double AB[Maxlocal][13];//开头为0的实部，0的虚部为0；之后从1到4的实部和虚部；
	CVector dcenter;

	for (i = 0; i < Maxlocal; i++)
		for (j = 0; j < 13; j++)
			AB[i][j] = 0.0;
	for (ii = 0; ii < Npolyhedron; ii++)
	{
		for (jj = 0; jj < Maxlocal; jj++)
		{
			m = pairnumbermin[ii][jj];
			//if (m == ii) { cout << "error!! 颗粒本身统计" << endl; getchar(); }
			dcenter = pairdcenter[ii][m];
			x1 = dcenter[0]; y1 = dcenter[1]; z1 = dcenter[2];
			x2 = x1*x1; y2 = y1*y1; z2 = z1*z1;
			x3 = x2*x1; y3 = y2*y1;
			x4 = x3*x1; y4 = y3*y1; z4 = z2*z2;
			x5 = x4*x1; y5 = y4*y1;
			x6 = x5*x1; y6 = y5*y1;
			r2 = x2 + y2 + z2; r4 = r2*r2; r6ni = 1.0 / (r2*r4);

			//m=0;
			t = (231.0*z2*z4 - 315.0*z4*r2 + 105.0*z2*r4)*r6ni - 5.0;
			AB[jj][0] = AB[jj][0] + t;
			//m=1;
			t = z1*(33.0*z4 - 30.0*z2*r2 + 5.0*r4)*r6ni;
			AB[jj][1] = AB[jj][1] + t*x1;
			AB[jj][2] = AB[jj][2] + t*y1;
			//m=2;
			t = (33.0*z4 - 18.0*z2*r2 + r4)*r6ni;
			AB[jj][3] = AB[jj][3] + t*(x2 - y2);
			AB[jj][4] = AB[jj][4] + t*2.0*x1*y1;
			//m=3;
			t = z1*(11.0*z2 - 3.0*r2)*r6ni;
			AB[jj][5] = AB[jj][5] + t*(x3 - 3.0*x1*y2);
			AB[jj][6] = AB[jj][6] + t*(3.0*x2*y1 - y3);
			//m=4;
			t = (11.0*z2 - r2) * r6ni;
			AB[jj][7] = AB[jj][7] + t*(x4 - 6.0*x2*y2 + y4);
			AB[jj][8] = AB[jj][8] + t*4.0*(x3*y1 - x1*y3);
			//m=5;
			t = z1 * r6ni;
			AB[jj][9] = AB[jj][9] + t*(x5 - 10.0*x3*y2 + 5.0*x1*y4);
			AB[jj][10] = AB[jj][10] + t*(5.0*x4*y1 - 10.0*x2*y3 + y5);
			//m=6;
			t = r6ni;
			AB[jj][11] = AB[jj][11] + t*(x6 - 15.0*x4*y2 + 15.0*x2*y4 - y6);
			AB[jj][12] = AB[jj][12] + t*(6.0*x5*y1 - 20.0*x3*y3 + 6.0*x1*y5);
		}
	}
	for (i = 1; i < Maxlocal; i++)
		for (j = 0; j < 13; j++)
			AB[i][j] = AB[i][j] + AB[i - 1][j];

	////2范数
	//Q6av = 0.0;
	//for (i = 0; i < Maxlocal; i++)
	//{
	//	q6[i] = 4.0*AB[i][0] * AB[i][0];
	//	q6[i] = q6[i] + 336.0*(AB[i][1] * AB[i][1] + AB[i][2] * AB[i][2]);
	//	q6[i] = q6[i] + 210.0*(AB[i][3] * AB[i][3] + AB[i][4] * AB[i][4]);
	//	q6[i] = q6[i] + 840.0*(AB[i][5] * AB[i][5] + AB[i][6] * AB[i][6]);
	//	q6[i] = q6[i] + 252.0*(AB[i][7] * AB[i][7] + AB[i][8] * AB[i][8]);
	//	q6[i] = q6[i] + 5544.0*(AB[i][9] * AB[i][9] + AB[i][10] * AB[i][10]);
	//	q6[i] = q6[i] + 462.0*(AB[i][11] * AB[i][11] + AB[i][12] * AB[i][12]);
	//	q6[i] = sqrt(q6[i]);
	//	Q6ev[i] = q6[i] / (32.0*double(Npolyhedron*(i + 1)));
	//	Q6ev[i] = fabs((Q6ev[i] - Q6random[i]) / Q6length[i]);
	//	Q6av = Q6av + Q6ev[i] * Q6ev[i];
	//}
	//Q6av = sqrt(Q6av / double(Maxlocal));

	//无穷范数
	Q6av = 0.0;
	for (i = 0; i < Maxlocal; i++)
	{
		q6[i] = 4.0*AB[i][0] * AB[i][0];
		q6[i] = q6[i] + 336.0*(AB[i][1] * AB[i][1] + AB[i][2] * AB[i][2]);
		q6[i] = q6[i] + 210.0*(AB[i][3] * AB[i][3] + AB[i][4] * AB[i][4]);
		q6[i] = q6[i] + 840.0*(AB[i][5] * AB[i][5] + AB[i][6] * AB[i][6]);
		q6[i] = q6[i] + 252.0*(AB[i][7] * AB[i][7] + AB[i][8] * AB[i][8]);
		q6[i] = q6[i] + 5544.0*(AB[i][9] * AB[i][9] + AB[i][10] * AB[i][10]);
		q6[i] = q6[i] + 462.0*(AB[i][11] * AB[i][11] + AB[i][12] * AB[i][12]);
		q6[i] = sqrt(q6[i]);
		Q6ev[i] = q6[i] / (32.0*double(Npolyhedron*(i + 1)));
		Q6ev[i] = fabs((Q6ev[i] - Q6localrandom[i]) / Q6locallength[i]);
		if (Q6av < Q6ev[i]) Q6av = Q6ev[i];
	}
	result << Q6av << "\t";
	return Q6av;
}
double CalculateS4local()
{
	int i, ii, jj, kk, ll, m;
	double coscita, c;
	double pmax, p[3];
	CVector pe[3], pf[3];
	double S4localmax;
	double S4localev[Maxlocal];

	for (i = 0; i < Maxlocal; i++) S4localev[i] = 0.0;
	for (ii = 0; ii < Npolyhedron; ii++)
	{
		for (i = 0; i < 3; i++)
		{
			pe[i] = polyhedra[ii]->e[i];
			p[i] = 0.0;
		}
		for (jj = 0; jj < Maxlocal; jj++)
		{
			m = pairnumbermin[ii][jj];
			//if (m == ii){ cout << "error!! 颗粒本身统计" << endl; getchar(); }
			for (kk = 0; kk < 3; kk++)
			{
				pf[kk] = polyhedra[m]->e[kk];
				for (ll = 0; ll < 3; ll++)
				{
					coscita = pe[ll] * pf[kk];
					c = coscita*coscita;
					p[ll] = p[ll] + (35.0*c*c - 30.0*c + 3.0) / 14.0;
				}
			}
			pmax = p[0];
			if (p[1] > pmax) pmax = p[1];
			if (p[2] > pmax) pmax = p[2];
			pmax = pmax / double(jj + 1);
			S4localev[jj] = S4localev[jj] + pmax;
		}
	}

	////2范数
	//S4localmax = 0.0;
	//for (i = 0; i < Maxlocal; i++)
	//{
	//	S4localev[i] = S4localev[i] / double(Npolyhedron);//result << i << "\t" << S4local[i] << endl;
	//	S4localev[i] = fabs((S4localev[i] - S4localrandom[i]) / S4locallength[i]); //result << i << "\t" << S4local[i] << endl;
	//	S4localmax = S4localmax + S4localev[i] * S4localev[i];		
	//}
	//S4localmax = sqrt(S4localmax / double(Maxlocal));

	//无穷范数
	S4localmax = 0.0;
	for (i = 0; i < Maxlocal; i++)
	{
		S4localev[i] = S4localev[i] / double(Npolyhedron);//result << i << "\t" << S4local[i] << endl;
		S4localev[i] = fabs((S4localev[i] - S4localrandom[i]) / S4locallength[i]); //result << i << "\t" << S4local[i] << endl;
		if (S4localmax < S4localev[i]) S4localmax = S4localev[i];
	}
	result << S4localmax << "\t";
	return S4localmax;
}
double CalculateS4Super()
{
	//if (NKind == 2 && PPara[0][4]>1.0)
	{
		int i, j, ii, jj;
		double p, p4, coscita, c;
		CVector trialvector;

		p4 = 0.0;
		for (i = 0; i < NpEach[0]; i++)
			for (j = 0; j < 3; j++)
			{
				trialvector = polyhedra[i]->e[j];
				p = 0.0;
				for (ii = 0; ii < NpEach[0]; ii++)
					for (jj = 0; jj < 3; jj++)
					{
						coscita = trialvector* polyhedra[ii]->e[jj];
						c = coscita*coscita;
						p = p + (35.0*c*c - 30.0*c + 3.0) / 14.0;
					}
				p = p / double(NpEach[0]);
				if (p > p4) p4 = p;
			}
		result << p4 << "\t";
		S4random = 0.97980 / sqrt(NpEach[0]);
		S4length = 0.27048 / sqrt(NpEach[0]);
		p4 = fabs((p4 - S4random) / S4length);
		result << p4 << "\t";
		return p4;
	}
	//else return -1;
}

//求交函数,Newton 迭代
double sgn(double x)
{
	double t;
	if (x == 0) t = 0.0;
	else if (x > 0) t = 1.0;
	else t = -1.0;
	return t;
}
double makeMnonsing(double eps, double scale, double M[3][3])//加主元使矩阵可逆，并返回改变后的行列式
{
	int i, j;
	double Mtemp1[3][3];
	double maxv, DetM;
	maxv = MatrixMax(M);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			Mtemp1[i][j] = 0.0;
		Mtemp1[i][i] = scale*maxv;
	}
	MatrixAdd(M, Mtemp1);
	for (i = 0; i < 10000; i++)
	{
		//cout << "DetM\t" << i << endl;
		DetM = MatrixDet(M);
		if (fabs(DetM) > eps)//M可逆
			break;
		else
			MatrixAdd(M, Mtemp1);
	}
	if (i == 10000)//主对角元处理不给力
	{
		cout << "too much trials for DetM!\t" << i << endl; getchar(); getchar(); getchar();
		DetM = 0.0;
	}
	return DetM;
}
double SuperballFunction(double ppara[5], CVector X)//x^2p+y^2p+z^2p
{
	int i;
	double logxx[3];
	for (i = 0; i < 3; i++)
		logxx[i] = log(fabs(X[i] / ppara[i]));
	return exp(ppara[4] / ppara[3] * log(exp(2.0*ppara[3] * logxx[0]) + exp(2.0*ppara[3] * logxx[1]))) + exp(2.0*ppara[4] * logxx[2]);
}
void SuperballFunctionDall(double ppara[5], CVector X, CVector &xn, double ddx[3][3])//一阶导数（向量）和二阶导数（矩阵）
{
	int i, j;
	double logPP, logxx[3];
	for (i = 0; i < 3; i++)
		logxx[i] = log(fabs(X[i] / ppara[i]));
	logPP = log(exp(2.0*ppara[3] * logxx[0]) + exp(2.0*ppara[3] * logxx[1]));
	for (i = 0; i < 2; i++)
		xn[i] = 2.0*ppara[4] * sgn(X[i]) / ppara[i] * exp((2.0*ppara[3] - 1.0)*logxx[i] + (ppara[4] / ppara[3] - 1.0)*logPP);
	xn[i] = 2.0*ppara[4] * sgn(X[i]) / ppara[i] * exp((2.0*ppara[4] - 1.0)*logxx[i]);

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			ddx[i][j] = 0;
	for (i = 0; i < 2; i++)
	{
		ddx[i][i] = (2.0*ppara[3] - 1.0)*exp((2.0*ppara[3] - 2.0)*logxx[i] + (ppara[4] / ppara[3] - 1.0)*logPP);
		ddx[i][i] += 2.0*(ppara[4] - ppara[3])*exp((4.0*ppara[3] - 2.0)*logxx[i] + (ppara[4] / ppara[3] - 2.0)*logPP);
		ddx[i][i] *= 2.0*ppara[4] / ppara[i] / ppara[i];
	}
	ddx[i][i] = 2.0*ppara[4] * (2.0*ppara[4] - 1.0) / ppara[i] / ppara[i] * exp((2.0*ppara[4] - 2.0)*logxx[i]);
	ddx[0][1] = 4.0*ppara[4] * (ppara[4] - ppara[3]) * sgn(X[0]) * sgn(X[1]) / ppara[0] / ppara[1] * exp((2.0*ppara[3] - 1.0)*(logxx[0] + logxx[1]) + (ppara[4] / ppara[3] - 2.0)*logPP);
	ddx[1][0] = ddx[0][1];
}
int Inspect(CPolyhedron * pa, CPolyhedron * pb)

{
	int i, j;
	double parp[2], ppara[2][5];//particle_p[1];		
	CVector ra, rb, r_ab, r_mid, ee0, ee1;//center of a,b;deta of ra, rb;mid of ra, rb; verticle to r_ab unit vector
	double Aa[3][3], AaT[3][3], Ab[3][3], AbT[3][3];//transformation matrix

													/**///目标函数值和自变量相关
	double XiAB;
	double lamda0, lamda, dlamda;
	CVector r_c0, dr_c0, r_c, drc, ra_l, rb_l;//local coordinates of c in a,b

											  /**///特定指代变量
	double xa, xb;//function value
	CVector dxa, dxb;//normal vector
	double ddxa[3][3], ddxb[3][3];//second order derivative of f
	double Xia, Xib;//Xi A,B的函数值 
	double ta, tb, ta1, tb1;//Xi  A,B 一阶导数局部和二阶导数局部
	CVector dxa2, dxb2;//Xi  A,B 一阶导数全局
	double ddxa2[3][3], ddxb2[3][3];//Xi  A,B 一阶导数全局
	double DetM, XiLamda;
	double M[3][3], Minv[3][3];
	CVector dg, dn;
	double dXi, dlength;

	/**///中间代换量
	int ab0, ab1, ab2;
	CVector rt, dxt, Vtemp;
	double ddxt[3][3], Mtemp1[3][3], Mtemp2[3][3];

	/**///迭代控制量
	double eps, eps1;//误差
	int step; //iteration steps in single repli	
	int repli; //replication time in this Newton method
	int replimax, numcita, numradius;//初值重置次数
	double cita, dcita, dradius;//初值移动相关
	double scale, coeff, distmax; //零处理，牛顿法控制步长	

								  /**///颗粒基本信息copy
	parp[0] = pa->pP[1];
	parp[1] = pb->pP[1];
	for (i = 0; i < 3; i++)
	{
		ppara[0][i] = pa->pA[i];
		ppara[1][i] = pb->pA[i];
	}
	for (i = 0; i < 2; i++)
	{
		ppara[0][i + 3] = pa->pP[i];
		ppara[1][i + 3] = pb->pP[i];
	}
	ra = pa->center;
	rb = pb->center;
	r_ab = rb.operator-(ra);
	r_mid = ra + r_ab*0.5;//middle point	
	for (i = 0; i < 3; i++)//get ee0,ee1;
	{
		if (fabs(r_ab[i])>ERROR1)
		{
			ab0 = i; ab1 = (i + 1) % 3; ab2 = (i + 2) % 3;
			ee0[ab1] = 1.0;	ee0[ab2] = 1.0; ee0[ab0] = -(r_ab[ab1] + r_ab[ab2]) / r_ab[ab0];
			ee1 = r_ab.Cross(ee0);
			ee0.Normalize(); ee1.Normalize();
			break;
		}
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			Aa[j][i] = pa->e[i][j];//a局部坐标的旋转矩阵
			AaT[i][j] = pa->e[i][j];//a旋转矩阵的逆矩阵
			Ab[j][i] = pb->e[i][j];//b局部坐标的旋转矩阵
			AbT[i][j] = pb->e[i][j];//b旋转矩阵的逆矩阵
		}

	/**///迭代控制
	eps = 1.0E-15;//收敛误差控制
	eps1 = 1.0E-26;//中间为零控制
	numcita = 36;//角度变化次数
	numradius = 20;//半径变化次数
	replimax = numcita *numradius + 1;
	dcita = 2.0*PI / double(numcita);
	dradius = 0.1;//半径加倍尺度
	scale = 0.01;//零处理

				 /**///迭代开始！！！！
	repli = -1;
	lamda0 = 0.5;
	r_c0 = r_mid;
reinitialize: //as New method fails
	repli++;
	if (repli > replimax)
	{
		cout << "too much replication! failed!!\t" << repli << endl; //outErr << "too much replication! failed!!\t" << repli << endl;
																	 //OutputPACerr(-1000);
																	 //getchar(); getchar();
		return 1;
	}
	if (repli > 0)
	{
		cita = dcita*double(repli % numcita);
		dr_c0 = ee0*cos(cita) + ee1*sin(cita);
		r_c0 = r_mid + dr_c0* (r_mid.Length()*dradius*(int(repli / numcita) + 1.0));
		lamda0 = (r_c0.operator-(ra)).Length() / ((r_c0.operator-(ra)).Length() + (r_c0.operator-(rb)).Length());
	}

	lamda = lamda0;
	r_c = r_c0;
	//local coordinates of r_c and f
	rt = r_c.operator-(ra);
	MatrixMult(AaT, rt, ra_l);
	rt = r_c.operator-(rb);
	MatrixMult(AbT, rt, rb_l);
	xa = SuperballFunction(ppara[0], ra_l);//中点在a局部坐标中的函数值
	xb = SuperballFunction(ppara[1], rb_l);//中点在b局部坐标中的函数值
	if ((xa < 1.0) && (xb < 1.0))	return 1;
	XiAB = lamda*(exp(1.0 / parp[0] * log(xa)) - 1.0) + (1.0 - lamda)*(exp(1.0 / parp[1] * log(xb)) - 1.0);
	for (step = 0; step < 50000; step++)
	{
		//get df, ddf 
		SuperballFunctionDall(ppara[0], ra_l, dxt, ddxt);//中点在a局部坐标中的一阶导数值和二阶导数值
		MatrixMult(Aa, dxt, dxa);//中点在颗粒a上的全局坐标中的一阶导数值
		MatrixMult(Aa, ddxt, ddxa);
		MatrixMult(ddxa, AaT, ddxt);
		MatrixCopy(ddxt, ddxa);//中点在颗粒a上的全局坐标中的二阶导数值

		SuperballFunctionDall(ppara[1], rb_l, dxt, ddxt);//中点在b局部坐标中的一阶导数值和二阶导数值
		MatrixMult(Ab, dxt, dxb);//中点在颗粒b上的全局坐标中的一阶导数值
		MatrixMult(Ab, ddxt, ddxb);
		MatrixMult(ddxb, AbT, ddxt);
		MatrixCopy(ddxt, ddxb);//中点在颗粒b上的全局坐标中的二阶导数值

							   /**///get Xi, dXi, ddXi  
		Xia = exp(1.0 / parp[0] * log(xa)) - 1.0;//Xi  A //Xia = pow(xa, 1.0 / p1) - 1.0; 
		ta = 1.0 / parp[0] * exp((double)(1.0 / parp[0] - 1.0)*log(xa));//Xi  A 一阶导数局部 //double ta = (1.0 / p1)*pow(xa, 1.0 / p1 - 1.0);			
		ta1 = (1.0 / parp[0] - 1.0)*ta / xa;//Xi A 二阶导数局部		
		dxa2 = dxa*ta; //Xi  A 一阶导数全局
		MatrixMult(ddxa, ta, ddxa2);
		Dyadic(dxa, dxa, Mtemp1);
		MatrixMult(Mtemp1, ta1, Mtemp2);
		MatrixAdd(ddxa2, Mtemp2);//Xi A 二阶导数全局

		Xib = exp(1.0 / parp[1] * log(xb)) - 1.0;//Xi  B //Xib = pow(xb, 1.0 / p2) - 1.0;
		tb = 1.0 / parp[1] * exp((double)(1.0 / parp[1] - 1.0)*log(xb));//Xi  B 一阶导数局部   //double tb = (1.0 / p2)*pow(xb, 1.0 / p2 - 1.0);
		tb1 = (1.0 / parp[1] - 1.0)*tb / xb;//caita B 二阶导数局部
		dxb2 = dxb*tb; //Xi  B 一阶导数全局
		MatrixMult(ddxb, tb, ddxb2);
		Dyadic(dxb, dxb, Mtemp1);
		MatrixMult(Mtemp1, tb1, Mtemp2);
		MatrixAdd(ddxb2, Mtemp2);//Xi  B 二阶导数全局		

								 /**/// get dlamda, drc;
		dXi = Xia - Xib;//XiAB 关于lamda的一阶导数
		dn = dxa2*lamda + dxb2*(1.0 - lamda); //XiAB 关于r_c的一阶导数 
		dlength = dXi*dXi + dn.LengthSquare();
		if (dlength < eps) break;
		MatrixMult(ddxa2, lamda, M);
		MatrixMult(ddxb2, 1.0 - lamda, Mtemp1);
		MatrixAdd(M, Mtemp1);//M = lamda*ddxa2 + (1-lamda)*ddxb2;	
		DetM = MatrixDet(M);
		if (fabs(DetM) < eps1)//M不可逆，主对角元处理
		{
			//cout << "DetM equal to zero!\t" << step << endl; outErr << "DetM equal to zero!\t" << step << endl;//cout <<"before\t"<< DetM <<"\t"<<M[0][0]<< endl;
			DetM = makeMnonsing(eps1, scale, M); //cout << "after\t" << DetM << "\t" << M[0][0] << endl;
		}
		MatrixInv(M, Minv, DetM); //if (step == 1) midinfor << DetM << "\t" << i << "\t";		
		dg = dxb2.operator-(dxa2);//dg = dxb2.operator-(dxa2);
		MatrixMult(Minv, dg, Vtemp);
		XiLamda = dg*Vtemp; //Xilamda = (dg')*M1*dg;
		if (fabs(XiLamda) < eps1)
		{
			//cout << "XiLamda eaqual to zero!!!" << endl; outErr << "XiLamda eaqual to zero!!!" << endl;
			XiLamda = -scale;
		}
		MatrixMult(Minv, dn, Vtemp);
		dlamda = 1.0 / XiLamda * (dXi + dg*Vtemp);
		Vtemp = (dg*dlamda).operator-(dn);
		MatrixMult(Minv, Vtemp, drc);//drc = M1 * (dlamda*dg-dn);	

									 /**///get coeff  for NR method		
		coeff = 0.5;
		distmax = drc.Length();
		if (distmax < fabs(dlamda)) distmax = fabs(dlamda);
		if (distmax > 1.0) coeff /= distmax;
		lamda = lamda + dlamda*coeff;
		r_c = r_c + drc*coeff;
		rt = r_c.operator-(ra);
		MatrixMult(AaT, rt, ra_l);
		rt = r_c.operator-(rb);
		MatrixMult(AbT, rt, rb_l);
		xa = SuperballFunction(ppara[0], ra_l);//中点在a局部坐标中的函数值
		xb = SuperballFunction(ppara[1], rb_l);//中点在b局部坐标中的函数值
		if ((xa < 1.0) && (xb < 1.0))	return 1;
		XiAB = lamda*(exp(1.0 / parp[0] * log(xa)) - 1.0) + (1.0 - lamda)*(exp(1.0 / parp[1] * log(xb)) - 1.0);
	}
	if ((lamda < 0) || (lamda>1))
	{
		//cout << "Lamda is not in [0,1], repli!!!\t" << lamda << endl; outErr << "Lamda is not in [0,1], repli!!!\t" << lamda << endl;
		goto reinitialize;
	}
	if (step == 50000)//cannot converge
	{
		cout << "too much step,repli!!\t" << step << "\t" << lamda << endl; //outErr << "too much step,repli!!\t" << step << "\t" << lamda << endl;
		goto reinitialize;
	}
	//midinfor << repli << "\t" << step << "\t" << repli * 2000 + step << "\t" << XiAB << "\t" << lamda << "\t" << r_c[0] << "\t" << r_c[1] << "\t" << r_c[2] << endl;
	//return XiAB;
	if (XiAB < 0.0) return 1;
	else return 0;
}
int Inspectnum(CPolyhedron * pa, CPolyhedron * pb)
{
	int ii;
	double XiAB0, XiAB1;
	XiAB0 = Inspectvalue(pa, pb); cout << XiAB0 << "\n";
	for (ii = 0; ii < 5; ii++)
	{
		XiAB1 = Inspectvalue(pa, pb); cout << ii << "\t" << XiAB1 << "\n";
		if (fabs((XiAB1 - XiAB0) / XiAB1) > ERROR1)
		{
			cout << "inspect not converge!\n";
			break;
		}
	}
	if (XiAB0 < 0.0) return 1;
	else return 0;
}
double Inspectvalue(CPolyhedron * pa, CPolyhedron * pb)
{
	int i, j;
	double parp[2], ppara[2][5];//particle_p[1];		
	CVector ra, rb, r_ab, r_mid, ee0, ee1;//center of a,b;deta of ra, rb;mid of ra, rb; verticle to r_ab unit vector
	double Aa[3][3], AaT[3][3], Ab[3][3], AbT[3][3];//transformation matrix

													/**///目标函数值和自变量相关
	double XiAB;
	double lamda0, lamda, dlamda;
	CVector r_c0, dr_c0, r_c, drc, ra_l, rb_l;//local coordinates of c in a,b

											  /**///特定指代变量
	double xa, xb;//function value
	CVector dxa, dxb;//normal vector
	double ddxa[3][3], ddxb[3][3];//second order derivative of f
	double Xia, Xib;//Xi A,B的函数值 
	double ta, tb, ta1, tb1;//Xi  A,B 一阶导数局部和二阶导数局部
	CVector dxa2, dxb2;//Xi  A,B 一阶导数全局
	double ddxa2[3][3], ddxb2[3][3];//Xi  A,B 一阶导数全局
	double DetM, XiLamda;
	double M[3][3], Minv[3][3];
	CVector dg, dn;
	double dXi, dlength;

	/**///中间代换量
	int ab0, ab1, ab2;
	CVector rt, dxt, Vtemp;
	double ddxt[3][3], Mtemp1[3][3], Mtemp2[3][3];

	/**///迭代控制量
	double eps, eps1;//误差
	int step; //iteration steps in single repli	
	int repli; //replication time in this Newton method
	int replimax, numcita, numradius;//初值重置次数
	double cita, dcita, dradius;//初值移动相关
	double scale, coeff, distmax; //零处理，牛顿法控制步长	

								  /**///颗粒基本信息copy
	parp[0] = pa->pP[1];
	parp[1] = pb->pP[1];
	for (i = 0; i < 3; i++)
	{
		ppara[0][i] = pa->pA[i];
		ppara[1][i] = pb->pA[i];
	}
	for (i = 0; i < 2; i++)
	{
		ppara[0][i + 3] = pa->pP[i];
		ppara[1][i + 3] = pb->pP[i];
	}
	ra = pa->center;
	rb = pb->center;
	r_ab = rb.operator-(ra);
	r_mid = ra + r_ab*0.5;//middle point	
	for (i = 0; i < 3; i++)//get ee0,ee1;
	{
		if (fabs(r_ab[i])>ERROR1)
		{
			ab0 = i; ab1 = (i + 1) % 3; ab2 = (i + 2) % 3;
			ee0[ab1] = 1.0;	ee0[ab2] = 1.0; ee0[ab0] = -(r_ab[ab1] + r_ab[ab2]) / r_ab[ab0];
			ee1 = r_ab.Cross(ee0);
			ee0.Normalize(); ee1.Normalize();
			break;
		}
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			Aa[j][i] = pa->e[i][j];//a局部坐标的旋转矩阵
			AaT[i][j] = pa->e[i][j];//a旋转矩阵的逆矩阵
			Ab[j][i] = pb->e[i][j];//b局部坐标的旋转矩阵
			AbT[i][j] = pb->e[i][j];//b旋转矩阵的逆矩阵
		}

	/**///迭代控制
	eps = 1.0E-15;//收敛误差控制
	eps1 = 1.0E-26;//中间为零控制
	numcita = 36;//角度变化次数
	numradius = 20;//半径变化次数
	replimax = numcita *numradius + 1;
	dcita = 2.0*PI / double(numcita);
	dradius = 0.1;//半径加倍尺度
	scale = 0.01;//零处理

				 /**///迭代开始！！！！
	repli = -1;
	lamda0 = 0.5;
	r_c0 = r_mid;
reinitialize: //as New method fails
	repli++;
	if (repli > replimax)
	{
		cout << "too much replication! failed!!\t" << repli << endl; //outErr << "too much replication! failed!!\t" << repli << endl;
																	 //OutputPACerr(-1000);
																	 //getchar(); getchar();
		return 1;
	}
	if (repli > 0)
	{
		cita = dcita*double(repli % numcita);
		dr_c0 = ee0*cos(cita) + ee1*sin(cita);
		r_c0 = r_mid + dr_c0* (r_mid.Length()*dradius*(int(repli / numcita) + 1.0));
		lamda0 = (r_c0.operator-(ra)).Length() / ((r_c0.operator-(ra)).Length() + (r_c0.operator-(rb)).Length());
	}

	lamda = lamda0;
	r_c = r_c0;
	//local coordinates of r_c and f
	rt = r_c.operator-(ra);
	MatrixMult(AaT, rt, ra_l);
	rt = r_c.operator-(rb);
	MatrixMult(AbT, rt, rb_l);
	xa = SuperballFunction(ppara[0], ra_l);//中点在a局部坐标中的函数值
	xb = SuperballFunction(ppara[1], rb_l);//中点在b局部坐标中的函数值
	if ((xa < 1.0) && (xb < 1.0))	return 1;
	XiAB = lamda*(exp(1.0 / parp[0] * log(xa)) - 1.0) + (1.0 - lamda)*(exp(1.0 / parp[1] * log(xb)) - 1.0);
	for (step = 0; step < 50000; step++)
	{
		//get df, ddf 
		SuperballFunctionDall(ppara[0], ra_l, dxt, ddxt);//中点在a局部坐标中的一阶导数值和二阶导数值
		MatrixMult(Aa, dxt, dxa);//中点在颗粒a上的全局坐标中的一阶导数值
		MatrixMult(Aa, ddxt, ddxa);
		MatrixMult(ddxa, AaT, ddxt);
		MatrixCopy(ddxt, ddxa);//中点在颗粒a上的全局坐标中的二阶导数值

		SuperballFunctionDall(ppara[1], rb_l, dxt, ddxt);//中点在b局部坐标中的一阶导数值和二阶导数值
		MatrixMult(Ab, dxt, dxb);//中点在颗粒b上的全局坐标中的一阶导数值
		MatrixMult(Ab, ddxt, ddxb);
		MatrixMult(ddxb, AbT, ddxt);
		MatrixCopy(ddxt, ddxb);//中点在颗粒b上的全局坐标中的二阶导数值

							   /**///get Xi, dXi, ddXi  
		Xia = exp(1.0 / parp[0] * log(xa)) - 1.0;//Xi  A //Xia = pow(xa, 1.0 / p1) - 1.0; 
		ta = 1.0 / parp[0] * exp((double)(1.0 / parp[0] - 1.0)*log(xa));//Xi  A 一阶导数局部 //double ta = (1.0 / p1)*pow(xa, 1.0 / p1 - 1.0);			
		ta1 = (1.0 / parp[0] - 1.0)*ta / xa;//Xi A 二阶导数局部		
		dxa2 = dxa*ta; //Xi  A 一阶导数全局
		MatrixMult(ddxa, ta, ddxa2);
		Dyadic(dxa, dxa, Mtemp1);
		MatrixMult(Mtemp1, ta1, Mtemp2);
		MatrixAdd(ddxa2, Mtemp2);//Xi A 二阶导数全局

		Xib = exp(1.0 / parp[1] * log(xb)) - 1.0;//Xi  B //Xib = pow(xb, 1.0 / p2) - 1.0;
		tb = 1.0 / parp[1] * exp((double)(1.0 / parp[1] - 1.0)*log(xb));//Xi  B 一阶导数局部   //double tb = (1.0 / p2)*pow(xb, 1.0 / p2 - 1.0);
		tb1 = (1.0 / parp[1] - 1.0)*tb / xb;//caita B 二阶导数局部
		dxb2 = dxb*tb; //Xi  B 一阶导数全局
		MatrixMult(ddxb, tb, ddxb2);
		Dyadic(dxb, dxb, Mtemp1);
		MatrixMult(Mtemp1, tb1, Mtemp2);
		MatrixAdd(ddxb2, Mtemp2);//Xi  B 二阶导数全局		

								 /**/// get dlamda, drc;
		dXi = Xia - Xib;//XiAB 关于lamda的一阶导数
		dn = dxa2*lamda + dxb2*(1.0 - lamda); //XiAB 关于r_c的一阶导数 
		dlength = dXi*dXi + dn.LengthSquare();
		if (dlength < eps) break;
		MatrixMult(ddxa2, lamda, M);
		MatrixMult(ddxb2, 1.0 - lamda, Mtemp1);
		MatrixAdd(M, Mtemp1);//M = lamda*ddxa2 + (1-lamda)*ddxb2;	
		DetM = MatrixDet(M);
		if (fabs(DetM) < eps1)//M不可逆，主对角元处理
		{
			//cout << "DetM equal to zero!\t" << step << endl; outErr << "DetM equal to zero!\t" << step << endl;//cout <<"before\t"<< DetM <<"\t"<<M[0][0]<< endl;
			DetM = makeMnonsing(eps1, scale, M); //cout << "after\t" << DetM << "\t" << M[0][0] << endl;
		}
		MatrixInv(M, Minv, DetM); //if (step == 1) midinfor << DetM << "\t" << i << "\t";		
		dg = dxb2.operator-(dxa2);//dg = dxb2.operator-(dxa2);
		MatrixMult(Minv, dg, Vtemp);
		XiLamda = dg*Vtemp; //Xilamda = (dg')*M1*dg;
		if (fabs(XiLamda) < eps1)
		{
			//cout << "XiLamda eaqual to zero!!!" << endl; outErr << "XiLamda eaqual to zero!!!" << endl;
			XiLamda = -scale;
		}
		MatrixMult(Minv, dn, Vtemp);
		dlamda = 1.0 / XiLamda * (dXi + dg*Vtemp);
		Vtemp = (dg*dlamda).operator-(dn);
		MatrixMult(Minv, Vtemp, drc);//drc = M1 * (dlamda*dg-dn);	

									 /**///get coeff  for NR method		
		coeff = 0.5;
		distmax = drc.Length();
		if (distmax < fabs(dlamda)) distmax = fabs(dlamda);
		if (distmax > 1.0) coeff /= distmax;
		lamda = lamda + dlamda*coeff;
		r_c = r_c + drc*coeff;
		rt = r_c.operator-(ra);
		MatrixMult(AaT, rt, ra_l);
		rt = r_c.operator-(rb);
		MatrixMult(AbT, rt, rb_l);
		xa = SuperballFunction(ppara[0], ra_l);//中点在a局部坐标中的函数值
		xb = SuperballFunction(ppara[1], rb_l);//中点在b局部坐标中的函数值
		if ((xa < 1.0) && (xb < 1.0))	return 1;
		XiAB = lamda*(exp(1.0 / parp[0] * log(xa)) - 1.0) + (1.0 - lamda)*(exp(1.0 / parp[1] * log(xb)) - 1.0);
	}
	if ((lamda < 0) || (lamda>1))
	{
		//cout << "Lamda is not in [0,1], repli!!!\t" << lamda << endl; outErr << "Lamda is not in [0,1], repli!!!\t" << lamda << endl;
		goto reinitialize;
	}
	if (step == 50000)//cannot converge
	{
		cout << "too much step,repli!!\t" << step << "\t" << lamda << endl; //outErr << "too much step,repli!!\t" << step << "\t" << lamda << endl;
		goto reinitialize;
	}
	//midinfor << repli << "\t" << step << "\t" << repli * 2000 + step << "\t" << XiAB << "\t" << lamda << "\t" << r_c[0] << "\t" << r_c[1] << "\t" << r_c[2] << endl;
	return XiAB;
	//if (XiAB < 0.0) return 1;
	//else return 0;
}

//voronoi剖分
void getsurpoint(double ppara[5])//颗粒表面点，申请空间，求局部拉格朗日坐标
{
	int num, numi;//剖分精度
	int i, j;
	double dv, v, u, cos1, sin1, cos2, sin2;
	double *pcos1, *psin1, *pcos2, *psin2;
	double *sgncos1, *sgnsin1, *sgnsin2;

	num = 20;
	dv = PI / (double)num;
	NSuP = 2 * (num*num - num + 1);
	SurfPot = new CVector[NSuP];
	pcos1 = new double[2 * num];
	psin1 = new double[2 * num];
	pcos2 = new double[num - 1];
	psin2 = new double[num - 1];
	sgncos1 = new double[2 * num];
	sgnsin1 = new double[2 * num];
	sgnsin2 = new double[num - 1];

	//v,cita (-0.5*PI,0.5*PI)
	for (i = 1; i < num; i++)
	{
		v = -0.5*PI + dv*i;
		cos2 = cos(v);
		sin2 = sin(v);
		pcos2[i - 1] = pow(fabs(cos2), 1.0 / ppara[4]);
		psin2[i - 1] = pow(fabs(sin2), 1.0 / ppara[4]);
		sgnsin2[i - 1] = sgn(sin2);
	}
	//u,fai [-PI,PI)
	for (j = 0; j < 2 * num; j++)
	{
		u = -PI + dv*j;
		cos1 = cos(u);
		sin1 = sin(u);
		pcos1[j] = pow(fabs(cos1), 1.0 / ppara[3]);
		psin1[j] = pow(fabs(sin1), 1.0 / ppara[3]);
		sgncos1[j] = sgn(cos1);
		sgnsin1[j] = sgn(sin1);
	}
	SurfPot[0][0] = 0.0;
	SurfPot[0][1] = 0.0;
	SurfPot[0][2] = -ppara[2];
	numi = 1;
	for (i = 1; i < num; i++)
		for (j = 0; j < 2 * num; j++)
		{
			SurfPot[numi][0] = ppara[0] * sgncos1[j] * pcos1[j] * pcos2[i - 1];
			SurfPot[numi][1] = ppara[1] * sgnsin1[j] * psin1[j] * pcos2[i - 1];
			SurfPot[numi][2] = ppara[2] * sgnsin2[i - 1] * psin2[i - 1];
			numi++;
		}
	SurfPot[numi][0] = 0.0;
	SurfPot[numi][1] = 0.0;
	SurfPot[numi][2] = ppara[2];

	delete[] pcos1; delete[] psin1; delete[] pcos2; delete[] psin2;
	delete[] sgncos1; delete[] sgnsin1; delete[] sgnsin2;

	/*cout << numi << "\t" << NSuP << endl;
	double radius = particle_id*0.01;
	ofstream scr("surfacepoint.scr");
	scr << fixed << setprecision(15);
	for (i = 0; i <NSuP; i++)
	{
	if (radius > 0)
	{
	scr << "sphere " << SurfPot[i][0] << "," << SurfPot[i][1] << "," << SurfPot[i][2] << " " << radius << endl;
	}
	}
	scr << "grid off\nVSCURRENT c\nzoom e\n";
	scr.close();*/
}
void OutputVoroPoint(int num)//输出用于voronoi剖分的表面离散点
{
	char cha[100];
	sprintf_s(cha, "surfpoint %d.txt", Counterfile);
	ofstream scr(cha);
	scr << fixed << setprecision(10);

	int i, j, m, n, ii, jj;
	CVector surfpl0, surfpe;
	CVector center;
	double A[3][3];

	//剖分相关
	int numi;//剖分精度	
	double dv, v, u, cos1, sin1, cos2, sin2;
	double *pcos1, *psin1, *pcos2, *psin2;
	double *sgncos1, *sgnsin1, *sgnsin2;
	
	dv = PI / (double)num;
	NSuP = 2 * (num*num - num + 1); cout << "num and NSuP:\t" << num << "\t" << NSuP << endl;
	SurfPot = new CVector[NSuP];
	pcos1 = new double[2 * num];
	psin1 = new double[2 * num];
	pcos2 = new double[num - 1];
	psin2 = new double[num - 1];
	sgncos1 = new double[2 * num];
	sgnsin1 = new double[2 * num];
	sgnsin2 = new double[num - 1];

	m = -1;
	for (ii = 0; ii < NKind; ii++)
	{
		//表面点离散//v,cita (-0.5*PI,0.5*PI)
		for (i = 1; i < num; i++)
		{
			v = -0.5*PI + dv*i;
			cos2 = cos(v);
			sin2 = sin(v);
			pcos2[i - 1] = pow(fabs(cos2), 1.0 / PPara[ii][4]);
			psin2[i - 1] = pow(fabs(sin2), 1.0 / PPara[ii][4]);
			sgnsin2[i - 1] = sgn(sin2);
		}
		//u,fai [-PI,PI)
		for (j = 0; j < 2 * num; j++)
		{
			u = -PI + dv*j;
			cos1 = cos(u);
			sin1 = sin(u);
			pcos1[j] = pow(fabs(cos1), 1.0 / PPara[ii][3]);
			psin1[j] = pow(fabs(sin1), 1.0 / PPara[ii][3]);
			sgncos1[j] = sgn(cos1);
			sgnsin1[j] = sgn(sin1);
		}
		SurfPot[0][0] = 0.0;
		SurfPot[0][1] = 0.0;
		SurfPot[0][2] = -PPara[ii][2];
		numi = 1;
		for (i = 1; i < num; i++)
			for (j = 0; j < 2 * num; j++)
			{
				SurfPot[numi][0] = PPara[ii][0] * sgncos1[j] * pcos1[j] * pcos2[i - 1];
				SurfPot[numi][1] = PPara[ii][1] * sgnsin1[j] * psin1[j] * pcos2[i - 1];
				SurfPot[numi][2] = PPara[ii][2] * sgnsin2[i - 1] * psin2[i - 1];
				numi++;
			}
		SurfPot[numi][0] = 0.0;
		SurfPot[numi][1] = 0.0;
		SurfPot[numi][2] = PPara[ii][2];
		//表面点输出
		for (jj = 0; jj < NpEach[ii]; jj++)
		{
			m++;
			center = polyhedra[m]->center;
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
				{
					A[j][i] = polyhedra[m]->e[i][j];//a局部坐标的旋转矩阵				
				}
			for (i = 0; i < NSuP; i++)
			{
				surfpl0 = SurfPot[i];
				MatrixMult(A, surfpl0, surfpe);
				surfpe = surfpe + center;
				PeriodCheckP(surfpe);
				n = m*NSuP + i;
				scr << n << "\t" << surfpe[0] << "\t" << surfpe[1] << "\t" << surfpe[2] << endl;
			}
		}
	}
	scr.close();	
	delete[] pcos1; delete[] psin1; delete[] pcos2; delete[] psin2;
	delete[] sgncos1; delete[] sgnsin1; delete[] sgnsin2;
	delete[] SurfPot;
}
int VoroVolume000(int num)//统计vorovolume
{
	char cha[100];
	sprintf_s(cha, "vorovol %d.txt", Counterfile);
	ifstream voro;
	voro.open(cha);
	if (!voro.good()) { cout << cha << "不存在" << endl; return 1; }

	sprintf_s(cha, "vororesults %d.txt", Counterfile);
	ofstream scr(cha);
	scr << fixed << setprecision(15);

	int i,n;
	double x, y, z, v;
	double *Vlocal, *PDlocal;
	double avVlocal, avPDlocal, detaVlocal, detaPDlocal;

	NSuP = 2 * (num*num - num + 1);
	PDlocal = new double[Npolyhedron];
	Vlocal = new double[Npolyhedron];
	for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
	voro >> n >> x >> y >> z >> v;
	while (1)
	{
		if (voro.eof()) break;
		//scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
		i = int(n / NSuP);
		Vlocal[i] += v;
		voro >> n >> x >> y >> z >> v;
	}

	avVlocal = 0.0; avPDlocal = 0.0;
	scr << "i\tVorov\tPDlocal" << endl;
	for (i = 0; i < Npolyhedron; i++)
	{
		Vlocal[i] /= polyhedra[i]->volume;
		avVlocal += Vlocal[i];
		PDlocal[i] = 1.0 / Vlocal[i];
		avPDlocal += PDlocal[i];
		scr << i << "\t" << Vlocal[i] << "\t" << PDlocal[i] << endl;
	}

	avVlocal /= double(Npolyhedron);
	avPDlocal /= double(Npolyhedron);

	detaVlocal = 0.0; detaPDlocal = 0.0;
	for (i = 0; i < Npolyhedron; i++)
	{
		v = Vlocal[i] - avVlocal;
		detaVlocal += v*v;
		v = PDlocal[i] - avPDlocal;
		detaPDlocal += v*v;
	}
	detaVlocal = sqrt(detaVlocal / double(Npolyhedron));
	detaPDlocal = sqrt(detaPDlocal / double(Npolyhedron));

	result << avVlocal << "\t" << avPDlocal << "\t" << detaVlocal << "\t" << detaPDlocal << "\t";
	voro.close(); scr.close();
	delete[] Vlocal; delete[] PDlocal;
	return 0;
}
int VoroVolume(int num)//统计vorovolume
{
	char cha[100];
	sprintf_s(cha, "vorovol %d.txt", Counterfile);
	ifstream voro;
	voro.open(cha);
	if (!voro.good()) { cout << cha << "不存在" << endl; return 1; }

	sprintf_s(cha, "vororesults %d.txt", Counterfile);
	ofstream scr(cha);
	scr << fixed << setprecision(15);

	int i, n, ii, jj;
	double x, y, z, v;
	double *Vlocal, *PDlocal;
	double *avVlocal, *avPDlocal, *detaVlocal, *detaPDlocal;
	avVlocal = new double[NKind + 1];
	avPDlocal = new double[NKind + 1];
	detaVlocal = new double[NKind + 1];
	detaPDlocal = new double[NKind + 1];


	NSuP = 2 * (num*num - num + 1); 
	cout << "num and NSuP:\t" << num << "\t" << NSuP << endl;
	result << num << "\t" << NSuP<< "\t";
	PDlocal = new double[Npolyhedron];
	Vlocal = new double[Npolyhedron];
	for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
	voro >> n >> x >> y >> z >> v;
	while (1)
	{
		if (voro.eof()) break;
		//scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
		i = int(n / NSuP);
		Vlocal[i] += v;
		voro >> n >> x >> y >> z >> v;
	}
	scr << "i\tkind\tVorov\tPDlocal" << endl;
	for (i = 0; i < Npolyhedron; i++)
	{
		Vlocal[i] /= polyhedra[i]->volume;
		PDlocal[i] = 1.0 / Vlocal[i];
		scr << i << "\t" << polyhedra[i]->kind << "\t" << Vlocal[i] << "\t" << PDlocal[i] << endl;
	}

	i = -1; 
	avVlocal[NKind] = 0.0;
	avPDlocal[NKind] = 0.0;
	for (ii = 0; ii < NKind; ii++)//种类数
	{
		avVlocal[ii] = 0.0;
		avPDlocal[ii] = 0.0;
		for (jj = 0; jj < NpEach[ii]; jj++)//每种的个数
		{
			i++;			
			avVlocal[ii] += Vlocal[i];
			avPDlocal[ii] += PDlocal[i];
		}
		avVlocal[NKind] += avVlocal[ii];
		avPDlocal[NKind] += avPDlocal[ii];
		if (NpEach[ii] > 0)
		{
			avVlocal[ii] /= double(NpEach[ii]);
			avPDlocal[ii] /= double(NpEach[ii]);
		}
	}
	avVlocal[NKind] /= double(Npolyhedron);
	avPDlocal[NKind] /= double(Npolyhedron);

	i = -1;
	detaVlocal[NKind] = 0.0;
	detaPDlocal[NKind] = 0.0;
	for (ii = 0; ii < NKind; ii++)//种类数
	{
		detaVlocal[ii] = 0.0;
		detaPDlocal[ii] = 0.0;
		for (jj = 0; jj < NpEach[ii]; jj++)//每种的个数
		{
			i++;
			v = Vlocal[i] - avVlocal[NKind];
			detaVlocal[NKind] += v*v;
			v = PDlocal[i] - avPDlocal[NKind];
			detaPDlocal[NKind] += v*v;
			v = Vlocal[i] - avVlocal[ii];
			detaVlocal[ii] += v*v;
			v = PDlocal[i] - avPDlocal[ii];
			detaPDlocal[ii] += v*v;
		}
		if (NpEach[ii] > 0)
		{
			detaVlocal[ii] = sqrt(detaVlocal[ii] / double(NpEach[ii]));
			detaPDlocal[ii] = sqrt(detaPDlocal[ii] / double(NpEach[ii]));
		}
	}
	detaVlocal[NKind] = sqrt(detaVlocal[NKind] / double(Npolyhedron));
	detaPDlocal[NKind] = sqrt(detaPDlocal[NKind] / double(Npolyhedron));

	result << avVlocal[NKind] << "\t" << detaVlocal[NKind] << "\t" << avPDlocal[NKind] << "\t" << detaPDlocal[NKind] << "\t";
	for (ii = 0; ii < NKind; ii++)	result << NpEach[ii] << "\t";
	for (ii = 0; ii < NKind; ii++)	result << avVlocal[ii] << "\t" << detaVlocal[ii] << "\t";
	for (ii = 0; ii < NKind; ii++)	result << avPDlocal[ii] << "\t" << detaPDlocal[ii] << "\t";
	
	voro.close(); scr.close();
	delete[] Vlocal; delete[] PDlocal;
	delete[]avVlocal; delete[]avPDlocal; delete[]detaVlocal; delete[]detaPDlocal;
	return 0;
}
void OutputVoroPointC()//输出用于voronoi剖分的表面离散点
{
	char cha[100];
	sprintf_s(cha, "surfpoint %d.txt", Counterfile);	
	ofstream scr(cha);
	scr << fixed << setprecision(15);

	int m;
	CVector center;
	for (m = 0; m < Npolyhedron; m++)
	{
		center = polyhedra[m]->center;
		scr << m << "\t" << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
	}
	scr.close();
}
int VoroVolumeC()//统计vorovolume
{
	char cha[100];
	sprintf_s(cha, "vorovol %d.txt", Counterfile);
	ifstream voro;
	voro.open(cha);
	if (!voro.good()) { cout << cha << "不存在" << endl; return 1; }

	sprintf_s(cha, "vororesults %d.txt", Counterfile);
	ofstream scr(cha);
	scr << fixed << setprecision(15);

	
	int i,n;
	double x, y, z, v;
	double *Vlocal, *PDlocal;
	double avVlocal, avPDlocal, detaVlocal, detaPDlocal;

	PDlocal = new double[Npolyhedron];
	Vlocal = new double[Npolyhedron];
	for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
	voro >> n >> x >> y >> z >> v;
	while (1)
	{
		if (voro.eof()) break;
		//scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
		//i = int(n / NSuP);
		i = n;
		Vlocal[i] += v;
		voro >> n >> x >> y >> z >> v;
	}

	avVlocal = 0.0; avPDlocal = 0.0;
	scr << "i\tVorov\tPDlocal" << endl;
	for (i = 0; i < Npolyhedron; i++)
	{
		Vlocal[i] /= polyhedra[i]->volume;
		avVlocal += Vlocal[i];
		PDlocal[i] = 1.0 / Vlocal[i];
		avPDlocal += PDlocal[i];
		scr << i << "\t" << Vlocal[i] << "\t" << PDlocal[i] << endl;
	}

	avVlocal /= double(Npolyhedron);
	avPDlocal /= double(Npolyhedron);

	detaVlocal = 0.0; detaPDlocal = 0.0;
	for (i = 0; i < Npolyhedron; i++)
	{
		v = Vlocal[i] - avVlocal;
		detaVlocal += v*v;
		v = PDlocal[i] - avPDlocal;
		detaPDlocal += v*v;
	}
	detaVlocal = sqrt(detaVlocal / double(Npolyhedron));
	detaPDlocal = sqrt(detaPDlocal / double(Npolyhedron));

	result << avVlocal << "\t" << avPDlocal << "\t" << detaVlocal << "\t" << detaPDlocal << "\t";
	voro.close(); scr.close();
	delete[] Vlocal; delete[] PDlocal;
	return 0;
}


/*************************************************************************************************************/

//Hyperuniform检测
void Structure_factor(double *Sk)
{
	int j;
	int nx, ny, nz;
	double kx, ky, kz;
	double K;
	double k;
	double re, im;
	double temp = 2 * PI / att[0][0];
	double mol;
	double sk;
	int x;
	
	
		
	for (j = 0; j < N_sam; j++)
	{
		Sk[j] = 0.0;
		Num[j] = 0;
	}
	for (nx = -42; nx < 43; nx++)
	{
		kx = temp * nx;
		for (ny = -42; ny < 43; ny++)
		{
			ky = temp * ny;
			for (nz = -42; nz < 43; nz++)
			{
				kz = temp * nz;
				if (nx == ny == nz == 0){continue;}

				k = sqrt(kx * kx + ky * ky + kz * kz);
				//cout << kx << "\t" << ky << "\t" << kz << "\t" << k << endl;
				K = k * PPara[1][0] / (2 * PI); //cout << K << endl;
				x = int(K / delta_k);
				//cout << x << endl;
				re = im = 0.0;
				for (j = 0; j < Npolyhedron; j++)
				{
					if (polyhedra[j]->kind == 0)
					{
						mol = kx * polyhedra[j]->center[0] + ky * polyhedra[j]->center[1] + kz * polyhedra[j]->center[2];
						re += cos(mol);
						im += sin(mol);
					}
					
				}
				sk = (re * re + im * im) / NpEach[0];//Npolyhedron;
				Sk[x] += sk;
				//cout << Sk[x] << endl;
				Num[x] += 1;
			}
		}
	}

	for (j = 0; j < N_sam; j++)
	{
		if (Num[j] == 0)
		{
			continue;
		}
		Sk[j] /= Num[j];
	}
}