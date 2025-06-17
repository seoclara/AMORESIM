#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "TVector3.h"
#include "TTree.h"
#include <iostream>
#include <getopt.h>

using namespace std;

void PrintHelp();

int main(int argc, char ** argv) {

    if (argc < 2) {
        PrintHelp();
        return -1;
    }

	// int n1;
    int opt;
	int nevt = 0;
    double e_start = 0;
    double e_end = 0;
	bool kSave = false;

	// if (argc!=2) {
		// std::cout <<"./NeutGen [entry]\n";
		// exit (-1);
	// }
	// n1 = atoi(argv[1]);

	// if(n1>0) nevt = n1;
	// else {
		// std::cout <<"Invalid entry\n";
		// exit (-1);
	// }

    while ((opt = getopt(argc, argv, "n:s:e:t")) != -1){
        switch (opt) {
            case 'n': nevt = atoi(optarg); break;
            case 's': e_start = atof(optarg); break;
            case 'e': e_end = atof(optarg); break;
			case 't': kSave = true; break;
            default: PrintHelp();
        }
    }

	// particle, energy, momentum
	double NeutMass = 0.9395654133; // GeV
	double TEnergy = 0; // GeV
	double KEnergy = 0; // GeV 
	double momentum_GeVx = 0;
	double momentum_GeVy = 0;
	double momentum_GeVz = 0;
	int PDGcode = 2112; // neutron

	// position selection
	double FlRadius = 6000;
	//double pospick = 0;	
	double xpos = 0;
	double ypos = 0;
	double zpos = 0;

	TRandom enernd = TRandom(time(0));
	TRandom3 *rnd3 = new TRandom3(0);

	TFile *savefile = nullptr;
	TTree *respara = nullptr;
	if (kSave){
		savefile = new TFile("neutgen.root", "RECREATE");
		respara = new TTree("gen_tree", "generated neutrons");
		respara->Branch("xpos", &xpos, "xpos/D");
		respara->Branch("ypos", &ypos, "ypos/D");
		respara->Branch("zpos", &zpos, "zpos/D");
		respara->Branch("px", &momentum_GeVx, "px/D");
		respara->Branch("py", &momentum_GeVy, "py/D");
		respara->Branch("pz", &momentum_GeVz, "pz/D");
		respara->Branch("ke", &KEnergy, "ke/D");
		respara->Branch("totene", &TEnergy, "totene/D");
	}

	for (int  i = 0; i < nevt; i++) {

		rnd3->Sphere(xpos, ypos, zpos, FlRadius);
		TVector3 spos = TVector3(xpos, ypos, zpos);

		TVector3 dir;
		rnd3->Sphere(dir[0], dir[1], dir[2], 1.0);

		if(spos.Dot(dir) > 0) {
			dir = -dir;
		}

		//KEnergy = h1->GetRandom(); // MeV
		KEnergy = enernd.Uniform(e_start, e_end);
		//TEnergy = KEnergy + NeutMass;
		TEnergy = sqrt(pow((KEnergy/1000. + NeutMass), 2) - pow(NeutMass, 2));

		momentum_GeVx = dir.X() * TEnergy;
		momentum_GeVy = dir.Y() * TEnergy;
		momentum_GeVz = dir.Z() * TEnergy;

		//momrnd.Sphere(momentum_GeVx, momentum_GeVy, momentum_GeVz, sqrt(TEnergy*TEnergy - NeutMass*NeutMass));
		//posrnd3.Sphere(xpos, ypos, zpos, FlRadius);
//*
		printf("1\n1 %d 0 0 %.10le %.10le %.10le %.10lf 0 %.10le %.10le %.10le\n", 
				PDGcode, momentum_GeVx, momentum_GeVy, momentum_GeVz, NeutMass, xpos, ypos, zpos);
//*/
		if (kSave) {
			respara->Fill();
		}
	}
	if (kSave) {
		savefile->cd();
		respara->Write();
		savefile->Close();
	}

	return 0;
}

void PrintHelp(){
    cout << endl;
    cout << "Usage: neutgen [-n # of event] [-s start energy] [-e end energy] [-t (optional for making output rootfile)]" << endl;
}
