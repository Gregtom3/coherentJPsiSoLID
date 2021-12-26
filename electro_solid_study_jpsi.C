#include "Lcore.h"
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

double tmin(double W, double Q2)
{
  const double m1_2 = -Q2;
  const double m2_2 = pow(0.938272,2);
  const double m3_2 = pow(3.0969,2);
  const double m4_2 = pow(0.938272,2);
  
  double p1 = sqrt(pow(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
  double p2 = sqrt(pow(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);

  return  m2_2 + m4_2 - 2 * sqrt(m2_2 + p1*p1) * sqrt(m4_2 + p2*p2) + 2 * p1 * p2;
}

double tmax(double W, double Q2)
{
  const double m1_2 = -Q2;
  const double m2_2 = pow(0.938272,2);
  const double m3_2 = pow(3.0969,2);
  const double m4_2 = pow(0.938272,2);
  
  double p1 = sqrt(pow(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
  double p2 = sqrt(pow(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);

  return  m2_2 + m4_2 - 2 * sqrt(m2_2 + p1*p1) * sqrt(m4_2 + p2*p2) - 2 * p1 * p2;
}

event myEvent;

int CheckAcceptance(const TLorentzVector P, const double thmin, const double thmax){
  if (P.Theta() < thmin * M_PI / 180.0) return 0;
  if (P.Theta() > thmax * M_PI / 180.0) return 0;
  if (P.P() < 0.3) return 0;
  return 1;
}
double CalcEg(const TLorentzVector P){
  double Md = 1.8756;
  double Mn = 0.93957;
  double Eg = (pow(P.E() - Md, 2) - pow(P.P(), 2) - Mn * Mn) / (2.0 * (P.E() - P.Pz() - Md));
  return Eg;
}
int electro_solid_study_jpsi(double Ebeam = 17.2, bool do_deuterium = false, Long64_t Nsim = 100000000, TString output_dir = ""){

  
  TString stringA;
  if(do_deuterium)
    stringA = "D";
  else
    stringA = "p";
   
  // Set simulation
  gRandom->SetSeed(0);
  



  // Electron beam energy and luminosity
  //  double Ebeam = 11;//GeV
  double lumi = 1.2e37; // events / cm^2 / s
  if(do_deuterium)
    lumi = lumi * 0.5;
  lumi = lumi * 1.0e-24; // 1e24 barn = 1 cm^2
  lumi = lumi * 1.0e-9; // 1 barn = 10^9 nano barn
  // events / nb / s
  double time = 3600.0 * 24 * 50;//s

  // Set nuclear
  if(do_deuterium)
    {
      NUCLEAR::SetNuclear("D");
      GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
      GENERATE::TF_fMomentum->SetNpx(1000);
    }
  else
    {
      NUCLEAR::SetNuclear("p");
    }

  // Set Psi2S production model
  JPSIMODEL::SetModel("23g");

  // Set scattered electron range
  double degtorad = M_PI / 180.0;
  //  GENERATE::cthrange[0] = ;
  //  GENERATE::cthrange[1] = 1;
  GENERATE::perange[0] = 0.0;//GeV
  GENERATE::perange[1] = Ebeam-7;//GeV

  // Set detector
  DETECTOR::SetDetector("SoLID");
  

  // TTree Info
  TString output_dir_save;
  if(output_dir.IsWhitespace())
    output_dir_save = Form("result-electro-jpsi/%s_solid_electro_%.1fGeV.root",stringA.Data(),Ebeam);
  else
    output_dir_save = Form("%s/%s_solid_electro_%.1fGeV.root",output_dir.Data(),stringA.Data(),Ebeam);
  
  TFile * fall2 = new TFile(output_dir_save, "RECREATE");
  TTree * tree = new TTree("tree","");
  tree->SetDirectory(fall2);

  TLorentzVector ki[2], kf[4], q;
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());

  double weight_jacobian = 0.0;
  double weight_nucleon_flux = 0.0;
  double weight_photon_flux = 0.0;
  double weight_dodt = 0.0;
  double weight_decay = 0.0;
  double weight_total = 0.0;
  double weight_smear = 0.0;
  
  double Mpsi = 3.0969;
  int is_sub;

  double Eg=0.0;
  double Eg_smear=0.0;
  double Eg_true=0.0;
  
  TLorentzVector *eIn = new TLorentzVector();
  TLorentzVector *pIn = new TLorentzVector();
  TLorentzVector *eOut = new TLorentzVector();
  TLorentzVector *pOut = new TLorentzVector();
  TLorentzVector *ePlusOut = new TLorentzVector();
  TLorentzVector *eMinusOut = new TLorentzVector();
  TLorentzVector _eOutSmear;
  TLorentzVector _pOutSmear;
  TLorentzVector _ePlusOutSmear;
  TLorentzVector _eMinusOutSmear;
  TLorentzVector *eOutSmear = new TLorentzVector();
  TLorentzVector *pOutSmear = new TLorentzVector();
  TLorentzVector *ePlusOutSmear = new TLorentzVector();
  TLorentzVector *eMinusOutSmear = new TLorentzVector();
  TLorentzVector vm;
  TLorentzVector vmSmear;
  TLorentzVector gamma;
  double vm_p = 0;
  double vm_theta = 0;
  double vm_mass_true = 0;
  double vm_mass_reco = 0;
  double is_accept_eOut = 0;
  double is_accept_pOut = 0;
  double is_accept_ePlusOut = 0;
  double is_accept_eMinusOut = 0;
  // Add TTree Branches
  tree->Branch("event", &myEvent.x,
	       "x/D:y/D:Q2/D:W/D:W2/D:t/D:x_true/D:y_true/D:Q2_true/D:W_true/D:W2_true/D:t_true/D:tmin/D:tmax/D");
  /* tree->Branch("eIn","TLorentzVector",&eIn);
  tree->Branch("pIn","TLorentzVector",&pIn);
  tree->Branch("eOut","TLorentzVector",&eOut);
  tree->Branch("pOut","TLorentzVector",&pOut);
  tree->Branch("ePlusOut","TLorentzVector",&ePlusOut);
  tree->Branch("eMinusOut","TLorentzVector",&eMinusOut);
  tree->Branch("eOutSmear","TLorentzVector",&eOutSmear);
  tree->Branch("pOutSmear","TLorentzVector",&pOutSmear);
  tree->Branch("ePlusOutSmear","TLorentzVector",&ePlusOutSmear);
  tree->Branch("eMinusOutSmear","TLorentzVector",&eMinusOutSmear);
  tree->Branch("vm","TLorentzVector",&vm);
  tree->Branch("gamma","TLorentzVector",&gamma);
  tree->Branch("vmSmear","TLorentzVector",&vmSmear);*/
   tree->Branch("weight_jacobian",&weight_jacobian,"Double/D");
  tree->Branch("weight_nucleon_flux",&weight_nucleon_flux,"Double/D");
  tree->Branch("weight_photon_flux",&weight_photon_flux,"Double/D");
  tree->Branch("weight_dodt",&weight_dodt,"Double/D");
  tree->Branch("weight_decay",&weight_decay,"Double/D");
  tree->Branch("weight_total",&weight_total,"Double/D");
  tree->Branch("weight_smear",&weight_smear,"Double/D");
  tree->Branch("is_accept_eOut",&is_accept_eOut,"Double/D");
  tree->Branch("is_accept_pOut",&is_accept_pOut,"Double/D");
  tree->Branch("is_accept_ePlusOut",&is_accept_ePlusOut,"Double/D");
  tree->Branch("is_accept_eMinusOut",&is_accept_eMinusOut,"Double/D");
  tree->Branch("Eg",&Eg,"Double/D");
  tree->Branch("Eg_smear",&Eg_smear,"Double/D");
  tree->Branch("Eg_true",&Eg_true,"Double/D");
  tree->Branch("is_sub",&is_sub,"Sub/I");
  tree->Branch("Nsim",&Nsim,"Sims/I");
  tree->Branch("lumi",&lumi,"Luminosity/D");
  tree->Branch("time",&time,"Time/D");
  tree->Branch("vm_p",&vm_p,"Double/D");
  tree->Branch("vm_theta",&vm_theta,"Double/D");
  tree->Branch("vm_mass_true",&vm_mass_true,"Double/D");
  tree->Branch("vm_mass_reco",&vm_mass_reco,"Double/D");


  TLorentzVector target_proton;
  target_proton.SetXYZM(0,0,0,0.9382721);
  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight_total = GENERATE::GetNucleon(&ki[1]);

    weight_total *= GENERATE::Event_eN2eNee_Jpsi(ki, kf); 

    if (weight_total > 0.0){
    
      q = ki[0] - kf[0];
      gamma = q;

      myEvent.x = - (q*q) / (2.0 * target_proton * q);
      myEvent.y = (target_proton * q ) / (target_proton * ki[0]);
      myEvent.Q2 = - (q * q);
      myEvent.W = sqrt( (target_proton + q) * (target_proton + q) );
      myEvent.W2 = (target_proton + q) * (target_proton + q);
      myEvent.t = (target_proton - kf[1]) * (target_proton - kf[1]);

      myEvent.x_true = - (q*q) / (2.0 * ki[1] * q);
      myEvent.y_true = (ki[1] * q ) / (ki[1] * ki[0]);
      myEvent.Q2_true = - (q * q);
      myEvent.W_true = sqrt( (ki[1] + q) * (ki[1] + q) );
      myEvent.W2_true = (ki[1] + q) * (ki[1] + q);
      myEvent.t_true = (ki[1] - kf[1]) * (ki[1] - kf[1]);


      _eOutSmear=kf[0];
      _pOutSmear=kf[1];
      _ePlusOutSmear=kf[2];
      _eMinusOutSmear=kf[3];
     
      // Grab all event weights
      weight_jacobian = GENERATE::WEIGHT_JACOBIAN;
      weight_nucleon_flux = GENERATE::WEIGHT_NUCLEON_FLUX;
      weight_photon_flux = GENERATE::WEIGHT_PHOTON_FLUX;
      weight_dodt = GENERATE::WEIGHT_DODT;
      weight_decay = GENERATE::WEIGHT_DECAY;

      // Calculate if the process was "subthreshold"
      is_sub = 0;
      if (q.E() < Mpsi + (Mpsi * Mpsi - q * q) / (2.0 * Mp)){
	is_sub = 1;
      }

      // Calculate unsmeared quantities
      Eg = CalcEg(kf[2]+kf[3]+kf[1]);
      Eg_true = q.E();
      vm=kf[2]+kf[3];

      // Calculate smeared quantities
      DETECTOR::SmearSoLID(_eOutSmear, "e-"); // smear final state electron, but its detection is unimportant for event reco
      is_accept_eOut = DETECTOR::AcceptanceSoLID(kf[0],"e-");
      is_accept_pOut = DETECTOR::AcceptanceSoLID(kf[1],"p");
      is_accept_ePlusOut = DETECTOR::AcceptanceSoLID(kf[2],"e+");
      is_accept_eMinusOut = DETECTOR::AcceptanceSoLID(kf[3],"e-");
      
      weight_smear = weight_total*DETECTOR::SmearSoLID(_pOutSmear, "p")*DETECTOR::SmearSoLID(_ePlusOutSmear, "e+")*DETECTOR::SmearSoLID(_eMinusOutSmear, "e-");
   
      Eg_smear = CalcEg(_ePlusOutSmear+_eMinusOutSmear+_pOutSmear);
      vmSmear=_ePlusOutSmear+_eMinusOutSmear;

      myEvent.tmin = tmin(myEvent.W_true, myEvent.Q2_true);
      myEvent.tmax = tmax(myEvent.W_true, myEvent.Q2_true);

      vm_p=vm.P();
      vm_theta=vm.Theta();
      vm_mass_true=vm.M();
      vm_mass_reco=(_ePlusOutSmear+_eMinusOutSmear).M();

      tree->Fill();

    }
  }

  //2
  fall2->Write();
  fall2->Close();


  return 0;
}
