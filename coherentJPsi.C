#include <iostream>
#include <string>

#include "Lcore.h"
#include "TFile.h"
#include "TTree.h"

struct event{
  double x;
  double y;
  double Q2;
  double W;
  double W2;
  double t;
  double x_true;
  double y_true;
  double Q2_true;
  double W_true;
  double W2_true;
  double t_true;
  double tmin;
  double tmax;
};

int main(int argc, char *argv[])
/*

  ./coherentJPsi <photo/electro> <events> <electronBeamEnergy> <time (days)> <kmin> <kmax> <filePrefix> <batchID>

 */
  
{

  //-----------------------------------------------------------------------------------------------------------
  // Parse through input
  //-----------------------------------------------------------------------------------------------------------

  const int MAX_ARG=9;
  int processID=0; // 0 = photo, 1 = electro
  string processStr="photo";
  int events=1000;  // Number of events
  double beamE=11;  // Incoming electron beam energy
  double days =50;  // Number of days of beamtime
  string fprefix = "test"; // output file name prefix
  int batchID =0;   // ID for batch
  double kmin = 0;  // minimum photon energy
  double kmax =10;  // maximum photon energy
  
  if(argc>MAX_ARG)
    {
      std::cout << "Too many arguments, " << argc << " > " << MAX_ARG << " , ending." << std::endl;
      return -1;
    }
  if(argc==MAX_ARG)
    {
      if(string(argv[1]).compare("photo")==0)
	{
	  processID=0;
	  processStr="photo";
	}
      else if(string(argv[1]).compare("electro")==0)
	{
	  processID=1;
	  processStr="electro";
	}
      else
	std::cout << "Unrecognized process, defaulting to photoproduction" << std::endl;

      events=std::stoi(argv[2]);
      beamE=std::stod(argv[3]);
      days=std::stod(argv[4]);
      kmin=std::stod(argv[5]);
      kmax=std::stod(argv[6]);
      fprefix=string(argv[7]);
      batchID=std::stoi(argv[8]);
    }
  else
    {
      std::cout << "Missing arguments, assuming photoproduction at beamE="<<beamE<<",events="<<events << std::endl;
    }

  //-----------------------------------------------------------------------------------------------------------
  // File location
  //-----------------------------------------------------------------------------------------------------------

  string fname = Form("./data/%s_beamE_%.2f_evts_%d_process_%s_batch_%d.root",fprefix.c_str(),beamE,events,processStr.c_str(),batchID);
  TFile *f = new TFile(TString(fname).Data(),"RECREATE");
  TTree *t = new TTree("tree","");
  t->SetDirectory(f);
  
  //-----------------------------------------------------------------------------------------------------------
  // Calculate luminosity and print event information
  //-----------------------------------------------------------------------------------------------------------

  const double lumi = 1.2e37 * 0.5 * 1.0e-24 * 1.0e-9 * 3600.0 * 24 * days;
  std::cout << "Coherent J/Psi Production Off Deuteron\n (gamma + D --> J/Psi + D)" << std::endl;
  std::cout << "Simulated Events = " << events << std::endl;
  std::cout << "Incoming Electron Beam Energy = " << beamE << std::endl;
  if(processID==0)
    {
      std::cout << "processType = Photoproduction" << std::endl;
      std::cout << "Real Photon Energy (GeV) " << kmin << " < Egamma < " << kmax << std::endl;
    }
  else
    {
      std::cout << "processType = Electroproduction" << std::endl;
      std::cout << "Virtual Photon Energy (GeV) " << kmin << " < Egamma < " << kmax << std::endl;
    }
  std::cout << "Luminosity = " << lumi/3600.0/24/days << " (events/nb/s)" << std::endl;
  std::cout << "Integrated Luminosity = " << lumi << " (events/nb) " << std::endl;
  std::cout << "File name = " << fname << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;

  //-----------------------------------------------------------------------------------------------------------
  // Initial parameters
  //-----------------------------------------------------------------------------------------------------------
  
  double mJPsi = 3.0969; // GeV/c^2
  double mD    = 2.014101778 * 0.9315; // GeV/c^2
  event myEvent;
  TLorentzVector ki[2],kf[4],q;
  // ki[0] --> Incoming electron
  // ki[1] --> Target Deuteron
  // kf[0] --> Scattered Electron
  // kf[1] --> Scattered Deuteron
  // kf[2] --> e-
  // kf[3] --> e+
  // q --> Real/Virtual Photon
  ki[0].SetXYZM(0,0,beamE,PARTICLE::e.M());
  ki[1].SetXYZM(0,0,0,mD);

  double weight = 0.0;
  
  JPSIMODEL::SetModel();  
  //-----------------------------------------------------------------------------------------------------------
  // Setting constraints
  //-----------------------------------------------------------------------------------------------------------
  
  if(processID==1)
    // minimum and maximum scattered electron energy for electroproduction
    {
      if(beamE-kmax<0)
	GENERATE::perange[0] = 0;
      else
	GENERATE::perange[0] = beamE-kmax;
      GENERATE::perange[1] = beamE-kmin;
    }
  //-----------------------------------------------------------------------------------------------------------
  // Tree Branching
  //-----------------------------------------------------------------------------------------------------------
  // t->Branch("event", &myEvent.x,
  //	       "x/D:y/D:Q2/D:W/D:W2/D:t/D:x_true/D:y_true/D:Q2_true/D:W_true/D:W2_true/D:t_true/D:tmin/D:tmax/D");
  t->Branch("eOut","TLorentzVector",&kf[0]);
  t->Branch("dOut","TLorentzVector",&kf[1]);
  t->Branch("ePlusOut","TLorentzVector",&kf[3]);
  t->Branch("eMinusOut","TLorentzVector",&kf[2]);
  t->Branch("gamma","TLorentzVector",&q);
  t->Branch("weight",&weight,"Double/D");
  //-----------------------------------------------------------------------------------------------------------
  // Event Generation
  //-----------------------------------------------------------------------------------------------------------
  for(int i = 0 ; i < events ; i++)
    {
      if(i%(events/100)==0)
	std::cout << "Finished " << i << " of " << events << " : (" << i/(events/100) << "% remaining)" << std::endl;
      if(processID==0)
	{
	  TLorentzVector eIn = ki[0];
	  GENERATE::SetBremsstrahlung();
	  weight = GENERATE::BremsstrahlungPhoton(&ki[0], kmin, kmax, beamE) * 1.95 / 2;//15cm LD2 target
	  weight *= GENERATE::Event_gD2Dee_Jpsi(ki, kf);
	  q=ki[0];
	  kf[0]=eIn-q;
	}      
      else
	{
	  weight = GENERATE::Event_eD2eDee_Jpsi(ki,kf);
	  q=ki[0]-kf[0];
	}

      t->Fill();
   
    }
  
  f->Write();
  f->Close();
  return 0;
}
