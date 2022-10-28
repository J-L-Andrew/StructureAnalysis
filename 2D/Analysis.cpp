#include "Analysis.h"
using namespace std;

#define ERROR1		1.000E-6
#define ERROR2		1.000E-12
#define PI			3.14159265358979324
#define Maxlocal	8

CPolyhedron **polyhedra;
CPolyhedron *pTmp[2];

int	   Npolyhedron;//��������
int    NKind;//�������
int	   *NpEach;//����ֿ�����
double *AerfEach;//������������
double *SEach;//����ֵ������������
double **PPara;//��״����
double	Pcdmax, Pidmin; //��������ֱ��,��С������뾶


double SumArea;//���������
double V_box;//���������
double PackingDensity;
double Lambda[2][2], L_inv[2][2];
CVector att[2];//att[0],att[1]��ʾ�����߽�����

double sblengthmin;//�߽����̳��ȵ�ƽ�������������������ȡֵ
double DISMAX;//=sblengthmin*10000��

double S2;//the nematic order
double S4;//tetratic order
double P4g;
double P6g;

//Pglobal
double **Gre4max;
double Gre4[Maxlocal];
double **Gre6max;
double Gre6[Maxlocal];
double **Gim4max;
double Gim4[Maxlocal];
double **Gim6max;
double Gim6[Maxlocal];

double GRe2, GRe4;
double GIm2, GIm4;

double **cita2;
double **cita4;

double meancn_1, meancn_2;


double **pairdis;//���������
double **pairdismin;//���������
int **pairnumber;//��������
int **pairnumbermin;//��������
CVector **pairdcenter;//��������������

CVector *SurfPot;//��������ĵ㣬����
int NSuP; //�����������

int Counterfile = 0;
ifstream myfile;
ofstream result("�������.txt");
//ofstream result("thresholdget.txt");
//ofstream midinfor("midinfor.txt");
ofstream midinfor("boundary0.txt");


/*��ȡ�ۺϴ��ļ�����*/
void main()
{
	srand((unsigned)time(NULL));
	int STRUCTURENUM, index, iteam, Ponum;
	string a2, a3;
	string pline;
	result << fixed << setprecision(8);
	midinfor << fixed << setprecision(15);

	Ponum = 30; //cout << "please input surface point number:\t"; cin >> Ponum;
	//for (index = 0; index < 11; index++)
	//{

	int num = 9;
	a2 = "packingstructure.txt";
	double s_2[30], s_4[30], g_4[30], g_6[30], mean1[30], mean2[30];
	double sm_2=0, sm_4=0, gm_4=0, gm_6=0, m1=0, m2=0;
	double sv_2=0, sv_4=0, gv_4=0, gv_6=0;
	
	for (iteam = 0; iteam < num; iteam++)
	{
		//cout << "iteam:\t" << iteam << endl;
		string a1;
		stringstream transfer;
		int a = iteam;
		transfer << a;
		transfer >> a1;
		a3 = a1 + a2;
		a3 = "E:\\Ell 3.0\\2.0\\" + a3;
		//if (index < 10)
		//{
			myfile.open(a3);
			
		//}
		//else
			//myfile.open("E:\\disk\\2.0\\%dpackingstructure.txt", iteam);
		
		if (!myfile.good()) { cout << "�ļ�������" << endl;	break; }
		
		result << "S2\tS4\tP4\tP6\t" << endl;


		//result << "Counter\tNsp\tPD\tNKind\tSupsp\tSupvolf\tSupsize\tNmix\tmixratio\tmixsp\tmixVolF\tmixsize\t";
		//result << "S4\tS4super\tS4superN\tQ6local\tS4local" << endl;
		//result << "num\tNSuP\tavVlocal\tdetaVlocal\tPDlocal\tdetaPDlocal\tNsuper\tNmix\tVsuper\tVSuperd\tVmix\tVmixd\tPDsuper\tPDsuperd\tPDmix\tPDmixd\n";
		Counterfile = -1;
		do
		{
			getline(myfile, pline);
			if (pline == "NEW")
			{
				myfile >> STRUCTURENUM;
				Counterfile++;
				//cout << Counterfile << endl;
				Initialize2(); //cout << "yes!\n";
				askfororderspace();
				{
					midinfor << Counterfile << "\t" << att[0][0] << endl;
					//midinfor << Counterfile << "\t" << att[0][0] << "\t" << polyhedra[0]->volume << endl;
					//result << Counterfile << "\t" << Npolyhedron << "\t" << PackingDensity << "\t" << PPara[0][2] / PPara[0][0] << "\t";
										
					contact(1.005);
					mean1[iteam] = meancn_1; m1 += meancn_1 / num;
					mean2[iteam] = meancn_2; m2 += meancn_2 / num;
					CalculateS(); 
					s_2[iteam] = S2; sm_2 += S2 / num;
					s_4[iteam] = S4; sm_4 += S4 / num;
					getppd();
					g_4[iteam] = CalculateP4global(); gm_4 += g_4[iteam] / num;
					g_6[iteam] = CalculateP6global(); gm_6 += g_6[iteam] / num;
					//cout << S2 << "\t" << S4 << "\t" << P4g << "\t" << P6g << endl;
					//if (NKind==1) CalculateS4local();
					
					result << Counterfile << "\t" << S2 << "\t" << S4 << "\t" << P4g << "\t" << P6g << endl;


					//�ļ����	
					//povray(); SCRcenter();
					//povrayeach();
					//SCR(); //SCRcenter(); povray();
					//OutputPAC();
					//OutputPAC(Counterfile);
					//OutputPACsample(Counterfile);

					//΢�۷���							
					//contact(1.001);														
					//ErrorAnalysis1();
					//ErrorAnalysis2();
					//IsotropyTest();
					//RadialDistributionFunction();	
					//OutputVoroPoint(Ponum);							
					//VoroVolume(Ponum);
					//OutputVoroPointC();
					//VoroVolumeC();
					
				}
				releaseallspace();
				releasepolyspace();
				releaseorderspace();
			}
		} while (!myfile.eof());
		myfile.close();
	//}
	result << "\n" << endl;
	result.close();
	midinfor.close();
	}
	for (int i = 0; i < num; i++)
	{
		sv_2 += (s_2[i] - sm_2)*(s_2[i] - sm_2) / num;
		sv_4 += (s_4[i] - sm_4)*(s_4[i] - sm_4) / num;
		gv_4 += (g_4[i] - gm_4)*(g_4[i] - gm_4) / num;
		gv_6 += (g_6[i] - gm_6)*(g_6[i] - gm_6) / num;
	}
	sv_2 = sqrt(sv_2);
	sv_4 = sqrt(sv_4);
	gv_4 = sqrt(gv_4);
	gv_6 = sqrt(gv_6);


	cout << sm_2 << "\t" << sm_4 << "\t" << gm_4 << "\t" << gm_6 << "\t" << meancn_1 << "\t" << meancn_2 << endl;


	cout << "press any key to continue!\n";
	getchar(); //getchar(); //getchar();
}



//��ʼ��
void Initialize2()
{
	int i, j;
	string line;
	double nonesize;

	do
	{
		getline(myfile, line);
		if (line == "PackingSpaceInformation")
		{
			for (i = 0; i < 2; i++)
				for (j = 0; j < 2; j++)
					myfile >> Lambda[j][i];
		}
		if (line == "Particle&Structure")
		{
			myfile >> nonesize >> nonesize >> nonesize >> nonesize;
			do
			{
				getline(myfile, line);
				if (line == "Superellipsoidevery")
				{
					myfile >> Npolyhedron >> NKind;
					askforpolyspace();
					for (int i = 0; i < NKind; i++)
					{
						myfile >> NpEach[i] >> AerfEach[i];
						//input >> NpEach[i] >> VolfEach[i];
						for (int j = 0; j < 3; j++)
							myfile >> PPara[i][j];
					}
					InitialParaAndSpace();
					for (int i = 0; i < Npolyhedron; i++)
					{
						for (int j = 0; j < 2; j++)
							myfile >> polyhedra[i]->center[j];
						for (int j = 0; j < 2; j++)
							for (int k = 0; k < 2; k++)
								myfile >> polyhedra[i]->e[j][k];
					}
				}
			} while (line != "END");
		}
	} while (line != "ENDOFFILE");


	Getboundary();
	PackingDensity = GetDensity();
}

//��ʼ����������
void InitialParaAndSpace()
{
	//��������Ϣ
	Generatepolyhedras();//��������ռ�	
	askforPPDspace();//���������ռ�	
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
void releaseallspace()//�ͷ����пռ�
{
	for (int i = 0; i < 2; i++) delete pTmp[i];
	for (int i = 0; i < Npolyhedron; i++) delete polyhedra[i];  delete[] polyhedra;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairdis[i]; delete[] pairdis;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairnumber[i]; delete[] pairnumber;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairdismin[i]; delete[] pairdismin;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairnumbermin[i]; delete[] pairnumbermin;
	for (int i = 0; i < Npolyhedron; i++) delete[] pairdcenter[i]; delete[] pairdcenter;
}
void askforpolyspace()//������Ԫ�ռ�����
{
	NpEach = new int[NKind];
	AerfEach = new double[NKind];
	SEach = new double[NKind];
	PPara = new double*[NKind];
	for (int i = 0; i < NKind; i++)
		PPara[i] = new double[3];
}
void releasepolyspace()//������Ԫ�ռ��ͷ�
{
	delete[] NpEach; delete[] AerfEach;
	delete[] SEach;
	for (int i = 0; i < NKind; i++)delete[] PPara[i]; delete[] PPara;

}
void releaseorderspace()//�ͷ�order�ռ�
{
	for (int i = 0; i < Npolyhedron; i++) delete cita2[i];  delete[] cita2;
	for (int i = 0; i < Npolyhedron; i++) delete cita4[i];  delete[] cita4;

	for (int i = 0; i < Npolyhedron; i++)delete[] Gre4max[i];
	delete[] Gre4max;
	for (int i = 0; i < Npolyhedron; i++)delete[] Gre6max[i];
	delete[] Gre6max;

	for (int i = 0; i < Npolyhedron; i++)delete[] Gim4max[i];
	delete[] Gim4max;
	for (int i = 0; i < Npolyhedron; i++)delete[] Gim6max[i];
	delete[] Gim6max;
}
void askfororderspace()//Ϊorder����ռ�
{
	cita2 = new double*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		cita2[i] = new double[2];
	cita4 = new double*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		cita4[i] = new double[2];

	Gre4max = new double*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		Gre4max[i] = new double[Maxlocal];
	Gre6max = new double*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		Gre6max[i] = new double[Maxlocal];

	Gim4max = new double*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		Gim4max[i] = new double[Maxlocal];
	Gim6max = new double*[Npolyhedron];
	for (int i = 0; i < Npolyhedron; i++)
		Gim6max[i] = new double[Maxlocal];
}
void Getboundary()//�߽��ʼ��
{
	for (int k = 0; k < 2; k++)
		for (int m = 0; m < 2; m++)
			att[k][m] = Lambda[m][k];
}
void Generatepolyhedras()//���ɿ���
{
	int i, m, mmi;
	double rrou[2];
	//CPolyhedron * pPolyhedron;

	for (int i = 0; i < 2; i++)
		pTmp[i] = new CPolyhedron();

	polyhedra = new CPolyhedron*[Npolyhedron];
	SumArea = 0.0; mmi = -1;
	for (m = 0; m < NKind; m++)
	{
		Getroutrin(PPara[m], rrou);
		AerfEach[m] = 0.0;
		for (i = 0; i < NpEach[m]; i++)
		{
			mmi++;
			polyhedra[mmi] = new CPolyhedron(PPara[m], rrou);
			polyhedra[mmi]->kind = m;
			AerfEach[m] += polyhedra[mmi]->area;
		}
		SumArea += AerfEach[m];
	}
	for (m = 0; m < NKind; m++)
		AerfEach[m] = AerfEach[m] / SumArea;
	//������������ֱ��
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
double GetDensity()
{
	double V_box = double(fabs(att[0][0] * att[1][1]));//�߽��������
	return SumArea / V_box;
}

void PeriodicalCheck(CPolyhedron * pp)
{
	for (int i = 0; i < 2; i++)
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

//�������״��Ϣ
void GetPolyInfor()//�õ���Ԫ�����������
{
	int i;

	for (i = 0; i < NKind; i++)
	{
		SEach[i] = GetArea(PPara[i]);
	}
}
void Getroutrin(double ppara[3], double rrout[2])
{
	

	if (ppara[2] > 0.0&&ppara[2] <= 1.0)
	{
		rrout[0] = ppara[0];
		rrout[1] = ppara[1];
	}
	else if (ppara[2]>1.0)
	{
        double k = ppara[2] / (ppara[2] - 1.0);
	    double m = pow(pow(ppara[0], 2.0*k) + pow(ppara[1], 2.0*k), 0.5 / k);
		rrout[0] = m;
		rrout[1] = ppara[1];
	}
	else if (ppara[2] <= 0.0)
	{
		cout << "p <= 0\tWrong!!!" << endl;
		exit(1);
	}
}
double GetArea(double ppara[3])
{
	double s = 0;
	double a1 = tgamma(0.5 / ppara[2]);
	double a2 = tgamma(1.0 / ppara[2]);
	s = ppara[0] * ppara[1] * a1 * a1 / (a2 * ppara[2]);
	return s;
}

//��λ��
void contact(double bilv)//�������߳�������ԭ����bilv��������λ��
{
	int m, n;
	int i, j;
	CVector V_PBC, centerm, dcenter, V;
	double AB_dis;
	double CCinerror, CCouterror;//�仯֮��������������ֱ��ƽ����	

	meancn_1 = 0;
	meancn_2 = 0;
	int component_1 = NpEach[0] + NpEach[1];
	for (m = 0; m < component_1; m++)
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
				{
					if (m == n && i == 0 && j == 0)
						continue;
					V_PBC = att[0] * i + att[1] * j;
					V = dcenter + V_PBC;
					AB_dis = V[0] * V[0] + V[1] * V[1];
					if (AB_dis < CCinerror)//���������ľ���С���ڽ���ֱ�����ض��ཻ
					{
						if (n < component_1)
							meancn_1 += 2;
						else if (n >= component_1)
						{
							meancn_1 += 1; meancn_2 += 1;
						}
						continue;
					}
					if (AB_dis > CCouterror)//���������ľ�����������ֱ�����ض�����
						continue;

					//��������ϸ�ж�
					pTmp[1]->Copyfrom(polyhedra[n]);
					pTmp[1]->ChangeSize(bilv);
					pTmp[1]->Jump(V_PBC);
					if (Inspect(pTmp[0], pTmp[1]) == 1)
					{
						if (n < component_1)
							meancn_1 += 2;
						else if (n >= component_1)
						{
							meancn_1 += 1; meancn_2 += 1;
						}
					}
				}
		}
	}

	for (m = component_1; m < Npolyhedron; m++)
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
				{
					if (m == n && i == 0 && j == 0)
						continue;
					V_PBC = att[0] * i + att[1] * j;
					V = dcenter + V_PBC;
					AB_dis = V[0] * V[0] + V[1] * V[1];
					if (AB_dis < CCinerror)//���������ľ���С���ڽ���ֱ�����ض��ཻ
					{
						meancn_2 += 2;
						continue;
					}
					if (AB_dis > CCouterror)//���������ľ�����������ֱ�����ض�����
						continue;

					//��������ϸ�ж�
					pTmp[1]->Copyfrom(polyhedra[n]);
					pTmp[1]->ChangeSize(bilv);
					pTmp[1]->Jump(V_PBC);
					if (Inspect(pTmp[0], pTmp[1]) == 1)
					{
						meancn_2 += 2;
					}
				}
		}
	}
	meancn_1 = meancn_1 / double(component_1);
	meancn_2 = meancn_2 / double(Npolyhedron - component_1);
}


//��Χ������ϵ
void getppd()//particle particle deitence
{
	int i, j, k, m, n;
	double dmin, dd;
	CVector PBC, center0, dcenter, V, Vdmin;

	sblengthmin = att[0][0] * att[0][0] + att[0][1] * att[0][1];
	for (i = 1; i < 2; i++)
	{
		dmin = att[i][0] * att[i][0] + att[i][1] * att[i][1];
		if (dmin < sblengthmin)
			sblengthmin = dmin;
	}
	DISMAX = sblengthmin*10000.0;
	for (m = 0; m < Npolyhedron; m++)
	{
		pairdis[m][m] = sblengthmin / 4.0;
		pairnumber[m][m] = -1;
		pairdcenter[m][m].Set(0.0, 0.0);
		center0 = polyhedra[m]->center;
		for (n = m + 1; n < Npolyhedron; n++)
		{
			dmin = DISMAX;
			dcenter = polyhedra[n]->center.operator- (center0);
			for (i = -1; i < 2; i++)
				for (j = -1; j < 2; j++)
				{
					PBC = att[0] * i + att[1] * j;
					V = dcenter + PBC;
					dd = V[0] * V[0] + V[1] * V[1];
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
	int i, j, n;
	double dmin, dd;
	CVector PBC, center0, dcenter, V, Vdmin;

	int t, midnumber, tempnumber;
	double middis, tempdis;

	pairdis[m][m] = sblengthmin / 2.0;
	pairdcenter[m][m].Set(0.0, 0.0);
	//if (pairnumber[m][m] >= 0){ cout << "error!! ��������ͳ��" << endl; getchar(); }
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
			{
				PBC = att[0] * i + att[1] * j;
				V = dcenter + PBC;
				dd = V[0] * V[0] + V[1] * V[1];
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

		//����һ��ȫ������
		//if ((pairnumber[n][m] >= 0) || (dmin <= pairdismin[n][Maxlocal - 1]))//˵����ԭ����m�Ѿ���n��min���е��л���ԭ��δ�������У�����Ҫ�������У���ʱ��Ҫ����min
		//{
		//	for (i = 0; i < Npolyhedron; i++) pairnumber[n][i] = -1;
		//	for (i = 0; i < Maxlocal; i++)
		//		getpairdismin(n, i);
		//}

		//���������ֲ�����
		if (pairnumber[n][m] >= 0)//˵����ԭ����m�Ѿ���n��min���е��У���Ҫ����min����
		{
			t = pairnumber[n][m];
			pairdismin[n][t] = dmin;//���´�ʱtλ�õ�minֵ

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
				if (pairdismin[n][i]>dmin)//���룺t֮ǰ��i��t���󣬴�i��ʼ�������,ֱ��t
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

			if (i == t)//��ʱ˵����t�仯֮����Ȼ��t֮ǰС��t֮ǰ����仯������t=0��
			{
				for (i = Maxlocal - 1; i > t; i--)
				{
					if (pairdismin[n][i] < dmin)//���룺t֮���i��t��С����i��ʼ��ǰ����, ֱ��t
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

			if (pairnumbermin[n][Maxlocal - 1] == m)//t�仯֮�������ȫ�ָ������һ��
			{
				pairnumber[n][m] = -1;
				getpairdismin(n, Maxlocal - 1);
			}
			continue;
		}

		else if (dmin <= pairdismin[n][Maxlocal - 1])//ԭ��δ��min���е��У�����Ҫ��min���е��У���ʱ��Ҫ����min
		{
			for (i = 0; i < Maxlocal; i++)
			{
				if (pairdismin[n][i]>dmin)//���룺i��dmin���󣬴�i��ʼ�������,ֱ��Maxlocal-1
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
					pairnumber[n][midnumber] = -1;//��ʱold�汾�����һ���Ŀ�������Ѿ���Ϊmidnumber�Ѿ��ų�����С�����У�������Ϊ-1
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

//�������
void CalculateS()
{
	double s = 0;
	double cita;
	double cos1 = 0, sin1 = 0;
	GRe2 = GIm2 = GRe4 = GIm4 = 0;
	int component_1 = NpEach[0] + NpEach[1];
	for (int i = 0; i<component_1; i++)
	{
		cos1 = polyhedra[i]->e[0][0];
		sin1 = polyhedra[i]->e[0][1];

		//Ell
		if (sin1 > 0)
		{
			cita = acos(cos1);
		}
		else
			cita = 2 * PI - acos(cos1);
		cos1 = cos(cita);
		sin1 = sin(cita);
		
		//Hyperciircle
		/*if (sin1 > 0)
		{
			cita = acos(cos1);
		}
		else
			cita = 2 * PI - acos(cos1);
		do
		{
			cita -= 0.5*PI;
		} while (cita > 0 && cita < 0.5*PI);
		cos1 = cos(cita);
		sin1 = sin(cita);*/

		cita2[i][1] = 2.0*cos1*sin1;
		GIm2 += cita2[i][1];
		cita2[i][0] = 2.0*cos1*cos1 - 1.0;
		GRe2 += cita2[i][0];
		cita4[i][1] = 2.0*cita2[i][1] * cita2[i][0];
		GIm4 += cita4[i][1];
		cita4[i][0] = 2.0*cita2[i][0] * cita2[i][0] - 1.0;
		GRe4 += cita4[i][0];
	}
	s = GRe2*GRe2 + GIm2*GIm2;
	S2 = sqrt(s) / component_1;
	s = GRe4*GRe4 + GIm4*GIm4;
	S4 = sqrt(s) / component_1;
}

double CalculateP4global()
{
	int i, ii, jj, m;
	double cos1, sin1;
	double cos2, sin2;
	double Re, Im;
	double L;
	CVector dcenter, center0;
	double P4globalmax;
	double P4reev[Maxlocal];
	double P4imev[Maxlocal];
	double P4globalev[Maxlocal];

	for (i = 0; i < Maxlocal; i++)
	{
		P4reev[i] = 0.0;
		P4imev[i] = 0.0;
	}
	int component_1 = NpEach[0] + NpEach[1];
	for (ii = 0; ii < component_1; ii++)
	{
		center0 = polyhedra[ii]->center;
		Re = 0.0; Im = 0.0;
		for (jj = 0; jj < Maxlocal; jj++)
		{
			m = pairnumbermin[ii][jj];
			if (m == ii) { cout << "error!! ��������ͳ��" << endl; getchar(); }
			dcenter = pairdcenter[ii][m];
			L = sqrt(dcenter[0] * dcenter[0] + dcenter[1] * dcenter[1]);
			//cout << L << endl;
			cos1 = dcenter[0] / L;



			sin1 = dcenter[1] / L;
			sin2 = 2.0*sin1*cos1;
			cos2 = 2.0*cos1*cos1 - 1.0;
			sin2 = 2.0*sin2*cos2;
			cos2 = 2.0*cos2*cos2 - 1.0;
			Re += cos2;
			Im += sin2;
			P4reev[jj] += Re / double(jj + 1);
			P4imev[jj] += Im / double(jj + 1);
			//cout << s << endl;
			//cout << pmax << endl;
		}
	}

	//�����
	P4globalmax = 0.0;
	for (i = 0; i < Maxlocal; i++)
	{
		P4globalev[i] = sqrt(P4reev[i] * P4reev[i] + P4imev[i] * P4imev[i]);
		P4globalev[i] /= double(component_1);
		//P4globalev[i] = fabs(P4globalev[i] - P4globalrandom[i]) / P4globallength[i];
		if (P4globalmax < P4globalev[i]) P4globalmax = P4globalev[i];
		P4globalmax = P4globalev[3];
	}
	return P4globalmax;
}

double CalculateP6global()
{
	int i, ii, jj, m;
	double cos1, sin1;
	double cos3, sin3;//������
	double Re, Im;
	double L;
	CVector dcenter, center0;
	double P6globalmax;
	double P6reev[Maxlocal];
	double P6imev[Maxlocal];
	double P6globalev[Maxlocal];

	for (i = 0; i < Maxlocal; i++)
	{
		P6reev[i] = 0.0;
		P6imev[i] = 0.0;
	}
	int component_1 = NpEach[0] + NpEach[1];
	for (ii = 0; ii < component_1; ii++)
	{
		center0 = polyhedra[ii]->center;
		Re = 0.0, Im = 0.0;
		for (jj = 0; jj < Maxlocal; jj++)
		{
			m = pairnumbermin[ii][jj];
			if (m == ii) { cout << "error!! ��������ͳ��" << endl; getchar(); }
			dcenter = pairdcenter[ii][m];
			L = sqrt(dcenter[0] * dcenter[0] + dcenter[1] * dcenter[1]);
			//cout << L << endl;
			cos1 = dcenter[0] / L;
			sin1 = dcenter[1] / L;
			sin3 = 3.0*sin1 - 4.0*sin1*sin1*sin1;
			cos3 = -3.0*cos1 + 4.0*cos1*cos1*cos1;
			sin3 = 2.0*sin3*cos3;
			cos3 = 2.0*cos3*cos3 - 1.0;
			Re += cos3;
			Im += sin3;
			P6reev[jj] += Re / double(jj + 1);
			P6imev[jj] += Im / double(jj + 1);
		}
	}

	//�����
	P6globalmax = 0.0;
	for (i = 0; i < Maxlocal; i++)
	{
		P6globalev[i] = sqrt(P6reev[i] * P6reev[i] + P6imev[i] * P6imev[i]);
		P6globalev[i] = P6globalev[i] / double(component_1);
		//P6globalev[i] = fabs(P6globalev[i] - P6globalrandom[i]) / P6globallength[i];
		if (P6globalmax < P6globalev[i]) P6globalmax = P6globalev[i];
		P6globalmax = P6globalev[5];
	}
	return P6globalmax;
}

double sgn(double x)
{
	double t;
	if (x == 0) t = 0.0;
	else if (x > 0) t = 1.0;
	else t = -1.0;
	return t;
}
double makeMnonsing(double eps, double scale, double M[2][2])//����Ԫʹ������棬�����ظı�������ʽ
{
	int i, j;
	double Mtemp1[2][2];
	double maxv, DetM;
	maxv = MatrixMax(M);
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
			Mtemp1[i][j] = 0.0;
		Mtemp1[i][i] = scale*maxv;
	}
	MatrixAdd(M, Mtemp1);
	for (i = 0; i < 10000; i++)
	{
		//cout << "DetM\t" << i << endl;
		DetM = MatrixDet(M);
		if (fabs(DetM) > eps)//M����
			break;
		else
			MatrixAdd(M, Mtemp1);
	}
	if (i == 10000)//���Խ�Ԫ��������
	{
		cout << "too much trials for DetM!\t" << i << endl; getchar(); getchar(); getchar();
		DetM = 0.0;
	}
	return DetM;
}
double EllipseFunction(double ra, double rb, double pa, CVector X)//
{
	double logxx[2];
	logxx[0] = log(fabs(X[0] / ra));
	logxx[1] = log(fabs(X[1] / rb));
	return exp(2 * pa * logxx[0]) + exp(2 * pa * logxx[1]);
}
void EllipseFunctionDall(double ra, double rb, double pa, CVector X, CVector &xn, double ddx[2][2])//һ�׵������������Ͷ��׵���������
{
	int i, j;
	double logxx[2];
	logxx[0] = log(fabs(X[0] / ra));
	logxx[1] = log(fabs(X[1] / rb));
	xn[0] = 2.0 * pa * sgn(X[0]) / ra * exp((2.0*pa - 1.0)*logxx[0]);
	xn[1] = 2.0 * pa * sgn(X[1]) / rb * exp((2.0*pa - 1.0)*logxx[1]);

	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			ddx[i][j] = 0;
	ddx[0][0] = 2.0*pa*(2.0*pa - 1.0)* exp((2.0*pa - 2.0)*logxx[0]) / (ra*ra);
	ddx[1][1] = 2.0*pa*(2.0*pa - 1.0)* exp((2.0*pa - 2.0)*logxx[1]) / (rb*rb);
}
int Inspect(CPolyhedron * pa, CPolyhedron * pb)

{
	int i, j;
	double parp[2], radiusa[2], radiusb[2];
	CVector ra, rb, r_ab, r_mid, ee0;//center of a,b;deta of ra, rb;mid of ra, rb; verticle to r_ab unit vector
	double Aa[2][2], AaT[2][2], Ab[2][2], AbT[2][2];//�ֲ�����;ֲ������ת�ã��������󣬵����棩

													/**///Ŀ�꺯��ֵ���Ա������
	double CtAB;//����AB
	double lambda0, lambda, dlambda;//��ʼ����ǰ�������Lambda
	CVector r_c0, dr_c0, r_c, drc, ra_l, rb_l;//ע����L��local�ֲ�����

											  /**///�ض�ָ������
	double xa, xb;//function�ĺ���ֵ
	CVector dxa, dxb;//һ�׵���
	double ddxa[2][2], ddxb[2][2];//���׵���
	double Cta, Ctb;//A,B�ĺ���ֵ 
	double ta, tb, ta1, tb1;//A,B һ�׵����ֲ��Ͷ��׵����ֲ�
	CVector dxa2, dxb2;//A,B ���׵���ȫ��
	double ddxa2[2][2], ddxb2[2][2];//A,B ���׵���ȫ��
	double DetM, CtLambda;
	double M[2][2], Minv[2][2];
	CVector dg, dn;
	double dXi, dlength;

	/**///�м������
	int ab0, ab1;
	CVector rt, dxt, Vtemp;
	double ddxt[2][2], Mtemp1[2][2], Mtemp2[2][2];

	/**///����������
	double eps, eps1;//���
	int step; //iteration steps in single repli	
	int repli; //replication time in this Newton method
	int replimax, numcita, numradius;//��ֵ���ô���
	double cita, dcita, dradius;//��ֵ�ƶ����
	double scale, coeff, distmax; //�㴦��ţ�ٷ����Ʋ���	

								  /**///����������Ϣcopy
	parp[0] = pa->pA[2];
	parp[1] = pb->pA[2];
	radiusa[0] = pa->pA[0];
	radiusb[0] = pa->pA[1];
	radiusa[1] = pb->pA[0];
	radiusb[1] = pb->pA[1];
	ra = pa->center;
	rb = pb->center;
	r_ab = rb.operator-(ra);
	r_mid = ra + r_ab*0.5;//middle point	
	for (i = 0; i < 2; i++)//�õ��������ߴ�ֱ����;
	{
		if (fabs(r_ab[i])>ERROR1)
		{
			ab0 = i; ab1 = (i + 1) % 2;
			ee0[ab1] = 1.0; ee0[ab0] = -r_ab[ab1] / r_ab[ab0];
			ee0.Normalize();
			break;
		}
	}
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
			Aa[j][i] = pa->e[i][j];//a�ֲ��������ת����e�Ǻ���洢��
			AaT[i][j] = pa->e[i][j];//a��ת����������
			Ab[j][i] = pb->e[i][j];//b�ֲ��������ת����
			AbT[i][j] = pb->e[i][j];//b��ת����������
		}

	/**///��������
	eps = 1.0E-9;//����������
	eps1 = 1.0E-23;//�м�Ϊ�����
	numcita = 2;//�Ƕȱ仯����
	numradius = 20;//�뾶�仯����
	replimax = 2 * numradius + 1;
	dcita = PI;
	dradius = 0.1;//�뾶�ӱ��߶�
	scale = 0.01;//�㴦��

				 /**///������ʼ��������
	repli = -1;
	lambda0 = 0.5;
	r_c0 = r_mid;
reinitialize: //as New method fails
	repli++;
	if (repli > replimax)
	{
		//cout << "too much replication! failed!!\t" << repli << endl; //outErr << "too much replication! failed!!\t" << repli << endl;
		//OutputPACerr(-1000);
		//getchar(); getchar();
		return 1;
	}
	if (repli > 0)
	{
		cita = dcita*double(repli % numcita);
		dr_c0 = ee0*cos(cita);
		r_c0 = r_mid + dr_c0* (r_mid.Length()*dradius*(int(repli / numcita) + 1.0));
		lambda0 = (r_c0.operator-(ra)).Length() / ((r_c0.operator-(ra)).Length() + (r_c0.operator-(rb)).Length());
	}

	lambda = lambda0;
	r_c = r_c0;
	//local coordinates of r_c and f
	rt = r_c.operator-(ra);
	MatrixMult(AaT, rt, ra_l);//ra_l���е���ڿ���a���ĵľֲ�����
	rt = r_c.operator-(rb);
	MatrixMult(AbT, rt, rb_l);//rb_l���е���ڿ���b���ĵľֲ�����
	xa = EllipseFunction(radiusa[0], radiusb[0], parp[0], ra_l);//�е���a�ֲ������еĺ���ֵ
	xb = EllipseFunction(radiusa[1], radiusb[1], parp[1], rb_l);//�е���b�ֲ������еĺ���ֵ
	if ((xa < 1.0) && (xb < 1.0))
		return 1;
	CtAB = lambda*(xa - 1.0) + (1.0 - lambda)*(xb - 1.0);
	for (step = 0; step < 50000; step++)
	{

		//�õ�df, ddf 
		EllipseFunctionDall(radiusa[0], radiusb[0], parp[0], ra_l, dxt, ddxt);//�е���a�ֲ������е�һ�׵���ֵ�Ͷ��׵���ֵ
		MatrixMult(Aa, dxt, dxa);//�е��ڿ���a�ϵ�ȫ�������е�һ�׵���ֵ
		MatrixMult(Aa, ddxt, ddxa);//
		MatrixMult(ddxa, AaT, ddxt);//ddxt���м��������������Ҫ������ҳ�
		MatrixCopy(ddxt, ddxa);//�е��ڿ���a�ϵ�ȫ�������еĶ��׵���ֵ

		EllipseFunctionDall(radiusa[1], radiusb[1], parp[1], rb_l, dxt, ddxt);//�е���a�ֲ������е�һ�׵���ֵ�Ͷ��׵���ֵ
		MatrixMult(Ab, dxt, dxb);//�е��ڿ���b�ϵ�ȫ�������е�һ�׵���ֵ
		MatrixMult(Ab, ddxt, ddxb);
		MatrixMult(ddxb, AbT, ddxt);
		MatrixCopy(ddxt, ddxb);//�е��ڿ���b�ϵ�ȫ�������еĶ��׵���ֵ


		/**///�õ�Xi, dXi, ddXi  
		Cta = exp(1.0 / parp[0] * log(xa)) - 1.0;//Xi  A //Xia = pow(xa, 1.0 / p1) - 1.0; 
		ta = 1.0 / parp[0] * exp((double)(1.0 / parp[0] - 1.0)*log(xa));//Xi  A һ�׵����ֲ� //double ta = (1.0 / p1)*pow(xa, 1.0 / p1 - 1.0);			
		ta1 = (1.0 / parp[0] - 1.0)*ta / xa;//Xi A ���׵����ֲ�		
		dxa2 = dxa*ta; //Xi  A һ�׵���ȫ��
		MatrixMult(ddxa, ta, ddxa2);
		Dyadic(dxa, dxa, Mtemp1);
		MatrixMult(Mtemp1, ta1, Mtemp2);
		MatrixAdd(ddxa2, Mtemp2);//Xi A ���׵���ȫ��

		Ctb = exp(1.0 / parp[1] * log(xb)) - 1.0;//Xi  B //Xib = pow(xb, 1.0 / p2) - 1.0;
		tb = 1.0 / parp[1] * exp((double)(1.0 / parp[1] - 1.0)*log(xb));//Xi  B һ�׵����ֲ�   //double tb = (1.0 / p2)*pow(xb, 1.0 / p2 - 1.0);
		tb1 = (1.0 / parp[1] - 1.0)*tb / xb;//caita B ���׵����ֲ�
		dxb2 = dxb*tb; //Xi  B һ�׵���ȫ��
		MatrixMult(ddxb, tb, ddxb2);
		Dyadic(dxb, dxb, Mtemp1);
		MatrixMult(Mtemp1, tb1, Mtemp2);
		MatrixAdd(ddxb2, Mtemp2);//Xi  B ���׵���ȫ��			

								 /**/// get dlamda, drc;
		dXi = Cta - Ctb;//XiAB ����lamda��һ�׵���
		dn = dxa2*lambda + dxb2*(1.0 - lambda); //XiAB ����r_c��һ�׵��� 
		dlength = dXi*dXi + dn.LengthSquare();
		if (dlength < eps) break;
		MatrixMult(ddxa2, lambda, M);
		MatrixMult(ddxb2, 1.0 - lambda, Mtemp1);
		MatrixAdd(M, Mtemp1);//M = lamda*ddxa2 + (1-lamda)*ddxb2;	
		DetM = MatrixDet(M);
		if (fabs(DetM) < eps1)//M�����棬���Խ�Ԫ����
		{
			//cout << "DetM equal to zero!\t" << step << endl; outErr << "DetM equal to zero!\t" << step << endl;//cout <<"before\t"<< DetM <<"\t"<<M[0][0]<< endl;
			DetM = makeMnonsing(eps1, scale, M); //cout << "after\t" << DetM << "\t" << M[0][0] << endl;
		}
		MatrixInv(M, Minv, DetM); //if (step == 1) midinfor << DetM << "\t" << i << "\t";		
		dg = dxb2.operator-(dxa2);//dg = dxb2.operator-(dxa2);
		MatrixMult(Minv, dg, Vtemp);
		CtLambda = dg*Vtemp; //Xilamda = (dg')*M1*dg;
		if (fabs(CtLambda) < eps1)
		{
			//cout << "CtLambda eaqual to zero!!!" << endl; outErr << "CtLamda eaqual to zero!!!" << endl;
			CtLambda = -scale;
		}
		MatrixMult(Minv, dn, Vtemp);
		dlambda = 1.0 / CtLambda * (dXi + dg*Vtemp);
		Vtemp = (dg*dlambda).operator-(dn);
		MatrixMult(Minv, Vtemp, drc);//drc = M1 * (dlamda*dg-dn);	

									 /**///get coeff  for NR method		
		coeff = 0.5;
		distmax = drc.Length();
		if (distmax < fabs(dlambda)) distmax = fabs(dlambda);
		if (distmax > 1.0) coeff /= distmax;
		lambda = lambda + dlambda*coeff;
		r_c = r_c + drc*coeff;
		rt = r_c.operator-(ra);
		MatrixMult(AaT, rt, ra_l);
		rt = r_c.operator-(rb);
		MatrixMult(AbT, rt, rb_l);
		xa = EllipseFunction(radiusa[0], radiusb[0], parp[0], ra_l);//�е���a�ֲ������еĺ���ֵ
		xb = EllipseFunction(radiusa[1], radiusb[1], parp[1], rb_l);//�е���b�ֲ������еĺ���ֵ
		if ((xa < 1.0) && (xb < 1.0))	return 1;
		CtAB = lambda*(xa - 1.0) + (1.0 - lambda)*(xb - 1.0);
	}
	if ((lambda < 0) || (lambda>1))
	{
		//cout << "Lambda is not in [0,1], repli!!!\t" << lambda << endl; outErr << "Lambda is not in [0,1], repli!!!\t" << lambda << endl;
		goto reinitialize;
	}
	if (step == 50000)//cannot converge
	{
		//cout << "too much step,repli!!\t" << step << "\t" << lamda << endl; //outErr << "too much step,repli!!\t" << step << "\t" << lamda << endl;
		goto reinitialize;
	}
	//midinfor << repli << "\t" << step << "\t" << repli * 2000 + step << "\t" << XiAB << "\t" << lamda << "\t" << r_c[0] << "\t" << r_c[1] << "\t" << r_c[2] << endl;
	//return XiAB;
	if (CtAB < 0.0) return 1;
	else return 0;
}

////voronoi�ʷ�
//void getsurpoint(double ppara[3])//��������㣬����ռ䣬��ֲ�������������
//{
//	int num, numi;//�ʷ־���
//	int i, j;
//	double dv, v, cos1, sin1;
//	double *pcos1, *psin1;
//	double *sgncos1, *sgnsin1, *sgnsin2;
//
//	num = 20;
//	dv = PI / (double)num;
//	SurfPot = new CVector[num];
//	pcos1 = new double[num - 1];
//	psin1 = new double[num - 1];
//
//	SurfPot[0][0] = -ppara[0];
//	SurfPot[0][1] = 0.0;
//
//	//v,cita (-0.5*PI,0.5*PI)
//	for (i = 1; i < num; i++)
//	{
//		v = -0.5*PI + dv*i;
//		cos1 = cos(v);
//		sin1 = sin(v);
//		pcos1[i - 1] = pow(fabs(cos1), 1.0 / ppara[2]);
//		psin1[i - 1] = pow(fabs(sin1), 1.0 / ppara[2]);
//
//		SurfPot[i][0] = ppara[0] * pcos1[i - 1];
//		SurfPot[i][1] = ppara[1] * psin1[i - 1];
//	}
//
//	delete[] pcos1; delete[] psin1;
//
//	/*cout << numi << "\t" << NSuP << endl;
//	double radius = particle_id*0.01;
//	ofstream scr("surfacepoint.scr");
//	scr << fixed << setprecision(15);
//	for (i = 0; i <NSuP; i++)
//	{
//	if (radius > 0)
//	{
//	scr << "sphere " << SurfPot[i][0] << "," << SurfPot[i][1] << "," << SurfPot[i][2] << " " << radius << endl;
//	}
//	}
//	scr << "grid off\nVSCURRENT c\nzoom e\n";
//	scr.close();*/
//}
//void OutputVoroPoint(int num)//�������voronoi�ʷֵı�����ɢ��
//{
//	char cha[100];
//	sprintf_s(cha, "surfpoint %d.txt", Counterfile);
//	ofstream scr(cha);
//	scr << fixed << setprecision(10);
//
//	int i, j, m, n, ii, jj;
//	CVector surfpl0, surfpe;
//	CVector center;
//	double A[3][3];
//
//	//�ʷ����
//	int numi;//�ʷ־���	
//	double dv, v, u, cos1, sin1, cos2, sin2;
//	double *pcos1, *psin1, *pcos2, *psin2;
//	double *sgncos1, *sgnsin1, *sgnsin2;
//
//	dv = PI / (double)num;
//	NSuP = 2 * (num*num - num + 1); cout << "num and NSuP:\t" << num << "\t" << NSuP << endl;
//	SurfPot = new CVector[NSuP];
//	pcos1 = new double[2 * num];
//	psin1 = new double[2 * num];
//	pcos2 = new double[num - 1];
//	psin2 = new double[num - 1];
//	sgncos1 = new double[2 * num];
//	sgnsin1 = new double[2 * num];
//	sgnsin2 = new double[num - 1];
//
//	m = -1;
//	for (ii = 0; ii < NKind; ii++)
//	{
//		//�������ɢ//v,cita (-0.5*PI,0.5*PI)
//		for (i = 1; i < num; i++)
//		{
//			v = -0.5*PI + dv*i;
//			cos2 = cos(v);
//			sin2 = sin(v);
//			pcos2[i - 1] = pow(fabs(cos2), 1.0 / PPara[ii][4]);
//			psin2[i - 1] = pow(fabs(sin2), 1.0 / PPara[ii][4]);
//			sgnsin2[i - 1] = sgn(sin2);
//		}
//		//u,fai [-PI,PI)
//		for (j = 0; j < 2 * num; j++)
//		{
//			u = -PI + dv*j;
//			cos1 = cos(u);
//			sin1 = sin(u);
//			pcos1[j] = pow(fabs(cos1), 1.0 / PPara[ii][3]);
//			psin1[j] = pow(fabs(sin1), 1.0 / PPara[ii][3]);
//			sgncos1[j] = sgn(cos1);
//			sgnsin1[j] = sgn(sin1);
//		}
//		SurfPot[0][0] = 0.0;
//		SurfPot[0][1] = 0.0;
//		SurfPot[0][2] = -PPara[ii][2];
//		numi = 1;
//		for (i = 1; i < num; i++)
//			for (j = 0; j < 2 * num; j++)
//			{
//				SurfPot[numi][0] = PPara[ii][0] * sgncos1[j] * pcos1[j] * pcos2[i - 1];
//				SurfPot[numi][1] = PPara[ii][1] * sgnsin1[j] * psin1[j] * pcos2[i - 1];
//				SurfPot[numi][2] = PPara[ii][2] * sgnsin2[i - 1] * psin2[i - 1];
//				numi++;
//			}
//		SurfPot[numi][0] = 0.0;
//		SurfPot[numi][1] = 0.0;
//		SurfPot[numi][2] = PPara[ii][2];
//		//��������
//		for (jj = 0; jj < NpEach[ii]; jj++)
//		{
//			m++;
//			center = polyhedra[m]->center;
//			for (i = 0; i < 3; i++)
//				for (j = 0; j < 3; j++)
//				{
//					A[j][i] = polyhedra[m]->e[i][j];//a�ֲ��������ת����				
//				}
//			for (i = 0; i < NSuP; i++)
//			{
//				surfpl0 = SurfPot[i];
//				MatrixMult(A, surfpl0, surfpe);
//				surfpe = surfpe + center;
//				PeriodCheckP(surfpe);
//				n = m*NSuP + i;
//				scr << n << "\t" << surfpe[0] << "\t" << surfpe[1] << "\t" << surfpe[2] << endl;
//			}
//		}
//	}
//	scr.close();
//	delete[] pcos1; delete[] psin1; delete[] pcos2; delete[] psin2;
//	delete[] sgncos1; delete[] sgnsin1; delete[] sgnsin2;
//	delete[] SurfPot;
//}
//int VoroVolume000(int num)//ͳ��vorovolume
//{
//	char cha[100];
//	sprintf_s(cha, "vorovol %d.txt", Counterfile);
//	ifstream voro;
//	voro.open(cha);
//	if (!voro.good()) { cout << cha << "������" << endl; return 1; }
//
//	sprintf_s(cha, "vororesults %d.txt", Counterfile);
//	ofstream scr(cha);
//	scr << fixed << setprecision(15);
//
//	int i, n;
//	double x, y, z, v;
//	double *Vlocal, *PDlocal;
//	double avVlocal, avPDlocal, detaVlocal, detaPDlocal;
//
//	NSuP = 2 * (num*num - num + 1);
//	PDlocal = new double[Npolyhedron];
//	Vlocal = new double[Npolyhedron];
//	for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
//	voro >> n >> x >> y >> z >> v;
//	while (1)
//	{
//		if (voro.eof()) break;
//		//scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
//		i = int(n / NSuP);
//		Vlocal[i] += v;
//		voro >> n >> x >> y >> z >> v;
//	}
//
//	avVlocal = 0.0; avPDlocal = 0.0;
//	scr << "i\tVorov\tPDlocal" << endl;
//	for (i = 0; i < Npolyhedron; i++)
//	{
//		Vlocal[i] /= polyhedra[i]->volume;
//		avVlocal += Vlocal[i];
//		PDlocal[i] = 1.0 / Vlocal[i];
//		avPDlocal += PDlocal[i];
//		scr << i << "\t" << Vlocal[i] << "\t" << PDlocal[i] << endl;
//	}
//
//	avVlocal /= double(Npolyhedron);
//	avPDlocal /= double(Npolyhedron);
//
//	detaVlocal = 0.0; detaPDlocal = 0.0;
//	for (i = 0; i < Npolyhedron; i++)
//	{
//		v = Vlocal[i] - avVlocal;
//		detaVlocal += v*v;
//		v = PDlocal[i] - avPDlocal;
//		detaPDlocal += v*v;
//	}
//	detaVlocal = sqrt(detaVlocal / double(Npolyhedron));
//	detaPDlocal = sqrt(detaPDlocal / double(Npolyhedron));
//
//	result << avVlocal << "\t" << avPDlocal << "\t" << detaVlocal << "\t" << detaPDlocal << "\t";
//	voro.close(); scr.close();
//	delete[] Vlocal; delete[] PDlocal;
//	return 0;
//}
//int VoroVolume(int num)//ͳ��vorovolume
//{
//	char cha[100];
//	sprintf_s(cha, "vorovol %d.txt", Counterfile);
//	ifstream voro;
//	voro.open(cha);
//	if (!voro.good()) { cout << cha << "������" << endl; return 1; }
//
//	sprintf_s(cha, "vororesults %d.txt", Counterfile);
//	ofstream scr(cha);
//	scr << fixed << setprecision(15);
//
//	int i, n, ii, jj;
//	double x, y, z, v;
//	double *Vlocal, *PDlocal;
//	double *avVlocal, *avPDlocal, *detaVlocal, *detaPDlocal;
//	avVlocal = new double[NKind + 1];
//	avPDlocal = new double[NKind + 1];
//	detaVlocal = new double[NKind + 1];
//	detaPDlocal = new double[NKind + 1];
//
//
//	NSuP = 2 * (num*num - num + 1);
//	cout << "num and NSuP:\t" << num << "\t" << NSuP << endl;
//	result << num << "\t" << NSuP << "\t";
//	PDlocal = new double[Npolyhedron];
//	Vlocal = new double[Npolyhedron];
//	for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
//	voro >> n >> x >> y >> z >> v;
//	while (1)
//	{
//		if (voro.eof()) break;
//		//scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
//		i = int(n / NSuP);
//		Vlocal[i] += v;
//		voro >> n >> x >> y >> z >> v;
//	}
//	scr << "i\tkind\tVorov\tPDlocal" << endl;
//	for (i = 0; i < Npolyhedron; i++)
//	{
//		Vlocal[i] /= polyhedra[i]->volume;
//		PDlocal[i] = 1.0 / Vlocal[i];
//		scr << i << "\t" << polyhedra[i]->kind << "\t" << Vlocal[i] << "\t" << PDlocal[i] << endl;
//	}
//
//	i = -1;
//	avVlocal[NKind] = 0.0;
//	avPDlocal[NKind] = 0.0;
//	for (ii = 0; ii < NKind; ii++)//������
//	{
//		avVlocal[ii] = 0.0;
//		avPDlocal[ii] = 0.0;
//		for (jj = 0; jj < NpEach[ii]; jj++)//ÿ�ֵĸ���
//		{
//			i++;
//			avVlocal[ii] += Vlocal[i];
//			avPDlocal[ii] += PDlocal[i];
//		}
//		avVlocal[NKind] += avVlocal[ii];
//		avPDlocal[NKind] += avPDlocal[ii];
//		if (NpEach[ii] > 0)
//		{
//			avVlocal[ii] /= double(NpEach[ii]);
//			avPDlocal[ii] /= double(NpEach[ii]);
//		}
//	}
//	avVlocal[NKind] /= double(Npolyhedron);
//	avPDlocal[NKind] /= double(Npolyhedron);
//
//	i = -1;
//	detaVlocal[NKind] = 0.0;
//	detaPDlocal[NKind] = 0.0;
//	for (ii = 0; ii < NKind; ii++)//������
//	{
//		detaVlocal[ii] = 0.0;
//		detaPDlocal[ii] = 0.0;
//		for (jj = 0; jj < NpEach[ii]; jj++)//ÿ�ֵĸ���
//		{
//			i++;
//			v = Vlocal[i] - avVlocal[NKind];
//			detaVlocal[NKind] += v*v;
//			v = PDlocal[i] - avPDlocal[NKind];
//			detaPDlocal[NKind] += v*v;
//			v = Vlocal[i] - avVlocal[ii];
//			detaVlocal[ii] += v*v;
//			v = PDlocal[i] - avPDlocal[ii];
//			detaPDlocal[ii] += v*v;
//		}
//		if (NpEach[ii] > 0)
//		{
//			detaVlocal[ii] = sqrt(detaVlocal[ii] / double(NpEach[ii]));
//			detaPDlocal[ii] = sqrt(detaPDlocal[ii] / double(NpEach[ii]));
//		}
//	}
//	detaVlocal[NKind] = sqrt(detaVlocal[NKind] / double(Npolyhedron));
//	detaPDlocal[NKind] = sqrt(detaPDlocal[NKind] / double(Npolyhedron));
//
//	result << avVlocal[NKind] << "\t" << detaVlocal[NKind] << "\t" << avPDlocal[NKind] << "\t" << detaPDlocal[NKind] << "\t";
//	for (ii = 0; ii < NKind; ii++)	result << NpEach[ii] << "\t";
//	for (ii = 0; ii < NKind; ii++)	result << avVlocal[ii] << "\t" << detaVlocal[ii] << "\t";
//	for (ii = 0; ii < NKind; ii++)	result << avPDlocal[ii] << "\t" << detaPDlocal[ii] << "\t";
//
//	voro.close(); scr.close();
//	delete[] Vlocal; delete[] PDlocal;
//	delete[]avVlocal; delete[]avPDlocal; delete[]detaVlocal; delete[]detaPDlocal;
//	return 0;
//}
//void OutputVoroPointC()//�������voronoi�ʷֵı�����ɢ��
//{
//	char cha[100];
//	sprintf_s(cha, "surfpoint %d.txt", Counterfile);
//	ofstream scr(cha);
//	scr << fixed << setprecision(15);
//
//	int m;
//	CVector center;
//	for (m = 0; m < Npolyhedron; m++)
//	{
//		center = polyhedra[m]->center;
//		scr << m << "\t" << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
//	}
//	scr.close();
//}
//int VoroVolumeC()//ͳ��vorovolume
//{
//	char cha[100];
//	sprintf_s(cha, "vorovol %d.txt", Counterfile);
//	ifstream voro;
//	voro.open(cha);
//	if (!voro.good()) { cout << cha << "������" << endl; return 1; }
//
//	sprintf_s(cha, "vororesults %d.txt", Counterfile);
//	ofstream scr(cha);
//	scr << fixed << setprecision(15);
//
//
//	int i, n;
//	double x, y, z, v;
//	double *Vlocal, *PDlocal;
//	double avVlocal, avPDlocal, detaVlocal, detaPDlocal;
//
//	PDlocal = new double[Npolyhedron];
//	Vlocal = new double[Npolyhedron];
//	for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
//	voro >> n >> x >> y >> z >> v;
//	while (1)
//	{
//		if (voro.eof()) break;
//		//scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
//		//i = int(n / NSuP);
//		i = n;
//		Vlocal[i] += v;
//		voro >> n >> x >> y >> z >> v;
//	}
//
//	avVlocal = 0.0; avPDlocal = 0.0;
//	scr << "i\tVorov\tPDlocal" << endl;
//	for (i = 0; i < Npolyhedron; i++)
//	{
//		Vlocal[i] /= polyhedra[i]->volume;
//		avVlocal += Vlocal[i];
//		PDlocal[i] = 1.0 / Vlocal[i];
//		avPDlocal += PDlocal[i];
//		scr << i << "\t" << Vlocal[i] << "\t" << PDlocal[i] << endl;
//	}
//
//	avVlocal /= double(Npolyhedron);
//	avPDlocal /= double(Npolyhedron);
//
//	detaVlocal = 0.0; detaPDlocal = 0.0;
//	for (i = 0; i < Npolyhedron; i++)
//	{
//		v = Vlocal[i] - avVlocal;
//		detaVlocal += v*v;
//		v = PDlocal[i] - avPDlocal;
//		detaPDlocal += v*v;
//	}
//	detaVlocal = sqrt(detaVlocal / double(Npolyhedron));
//	detaPDlocal = sqrt(detaPDlocal / double(Npolyhedron));
//
//	result << avVlocal << "\t" << avPDlocal << "\t" << detaVlocal << "\t" << detaPDlocal << "\t";
//	voro.close(); scr.close();
//	delete[] Vlocal; delete[] PDlocal;
//	return 0;
//}
