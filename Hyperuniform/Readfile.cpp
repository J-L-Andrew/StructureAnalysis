#include "Analysis.h"
using namespace std;



void Load_file()
{
	srand((unsigned)time(NULL));
	int STRUCTURENUM, iteam, Ponum;
	string pline;

	for (iteam = 0; iteam < 1; iteam++)
	{
		cout << "iteam:\t" << iteam << endl;
		input.open("totalstructure.txt");
		if (!input.good()) { cout << "totalstructure.txt不存在" << endl;	break; }



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

							   //if (Counterfile==526)
							   //if (Counterfile == 18 || Counterfile == 13)
							   //if (Counterfile == 1 || Counterfile == 150 || Counterfile == 162 || Counterfile == 218)// || Counterfile == 170 || Counterfile == 191 || Counterfile == 351 || Counterfile == 496 || Counterfile == 520 || Counterfile == 529 || Counterfile == 534 || Counterfile == 541)
							   //if ((Counterfile >29 && Counterfile <35) || (Counterfile >54 && Counterfile <60) || (Counterfile >209 && Counterfile <215) || (Counterfile >234 && Counterfile <240))
							   //if (Counterfile >184 && Counterfile <190)
				{
					midinfor << Counterfile << "\t" << att[0][0] << endl;
					//midinfor << Counterfile << "\t" << att[0][0] << "\t" << polyhedra[0]->volume << endl;
					//result << Counterfile << "\t" << Npolyhedron << "\t" << PackingDensity << "\t" << PPara[0][2] / PPara[0][0] << "\t";
					result << Counterfile << "\t" << Npolyhedron << "\t" << PackingDensity << "\t" << NKind << "\t";
					OutPolyInfor();
					//GetPInfo(PPara[0]);							
					//GetEllipsoidity(PPara[0]);

					/*CalculateS4();
					CalculateS4Super();
					getppd();
					CalculateQ6local();
					if (NKind==1) CalculateS4local();*/


					//文件输出	
					//povray(); SCRcenter();
					//povrayeach();
					//SCR(); //SCRcenter(); povray();
					OutputPAC();
					//OutputPAC(Counterfile);
					//OutputPACsample(Counterfile);

					//微观分析							
					//contact(1.001);														
					//ErrorAnalysis1();
					//ErrorAnalysis2();
					//IsotropyTest();
					//RadialDistributionFunction();	
					//OutputVoroPoint(Ponum);							
					//VoroVolume(Ponum);
					//OutputVoroPointC();
					//VoroVolumeC();
					result << endl;
				}
				releaseallspace();
				releasepolyspace();
			}
		} while (!input.eof());
		input.close();
	}

	result.close();
	midinfor.close();
	cout << "press any key to continue!\n";
	getchar(); //getchar(); //getchar();
}

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
	Releasetononoverlap(1.000000001);
	PackingDensity = GetDensity();
}