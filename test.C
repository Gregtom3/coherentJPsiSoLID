#include "TFile.h"
double tmin(double W, double Q2)
{
  const double m1_2 = -Q2;
  const double m2_2 = pow(0.938272*2,2);
  const double m3_2 = pow(3.0969,2);
  const double m4_2 = pow(0.938272*2,2);
  
  double p1 = sqrt(pow(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
  double p2 = sqrt(pow(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);

  return  m2_2 + m4_2 - 2 * sqrt(m2_2 + p1*p1) * sqrt(m4_2 + p2*p2) + 2 * p1 * p2;
}

double tmax(double W, double Q2)
{
  const double m1_2 = -Q2;
  const double m2_2 = pow(0.938272*2,2);
  const double m3_2 = pow(3.0969,2);
  const double m4_2 = pow(0.938272*2,2);
  
  double p1 = sqrt(pow(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
  double p2 = sqrt(pow(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);

  return  m2_2 + m4_2 - 2 * sqrt(m2_2 + p1*p1) * sqrt(m4_2 + p2*p2) - 2 * p1 * p2;
}
TGraph2D *tg = new TGraph2D(16*9);
double dsdt(Int_t nDim, Double_t *Xarg)
{
  return tg->Interpolate(Xarg[0],Xarg[1]);
}
int test()
{
  int E[] = {6,7,8,9,10,12,14,16,18};
  for(int i = 0 ; i < 9 ; i++)
    {
      ifstream infile(Form("tables/relwf/dsdt-v18-fold-rel-%dgev.out",E[i])); 
      char tmp[100];
      infile.getline(tmp,100);
      double beamE, t, tmtmin, dsdt;
      for(int j = 0*i ; j < 16*i ; j++)
	{
	  infile >> beamE >> t >> tmtmin >> dsdt;
	  //	  cout << beamE << " " << t << " " << dsdt << endl;
	  tg->SetPoint(j,beamE,t,dsdt);
	  if(infile.eof())
	    break;
	}
    }
  //cout << tg->Interpolate(5.6,0.15) << endl;
  double W = 5.94;
  double Q2=0;
  //  cout << tmin(W,Q2) << " " << tmax(W,Q2) << endl;
  cout<<"--- kanwa started ---"<<endl;
  TH2D  *hst_xy = new TH2D("hst_xy" ,  "x-y plot", 50,0,1.0, 50,0,1.0);
  Double_t MCvect[2]; // 2-dim vector generated in the MC run
  //
  TRandom     *PseRan   = new TRandom3();  // Create random number generator
  PseRan->SetSeed(4357);
  TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
  FoamX->SetkDim(2);         // No. of dimensions, obligatory!
  FoamX->SetnCells(500);     // Optionally No. of cells, default=2000
  FoamX->SetRhoInt(dsdt);  // Set 2-dim distribution, included below
  FoamX->SetPseRan(PseRan);  // Set random number generator
  FoamX->Initialize();       // Initialize simulator, may take time...
  //
  // visualising generated distribution
  TCanvas *cKanwa = new TCanvas("cKanwa","Canvas for plotting",600,600);
  cout << "H" << endl;
  cKanwa->cd();
  // From now on FoamX is ready to generate events
  int nshow=5000;
  for(long loop=0; loop<100000; loop++){
    FoamX->MakeEvent();            // generate MC event
    FoamX->GetMCvect( MCvect);     // get generated vector (x,y)
    Double_t x=MCvect[0];
    Double_t y=MCvect[1];
    if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
    hst_xy->Fill(x,y);
    // live plot
    if(loop == nshow){
      nshow += 5000;
      hst_xy->Draw("lego2");
      cKanwa->Update();
    }
  }// loop
   //
  hst_xy->Draw("lego2");  // final plot
  cKanwa->Update();
  //
  Double_t MCresult, MCerror;
  FoamX->GetIntegMC( MCresult, MCerror);  // get MC integral, should be one
  cout << " MCresult= " << MCresult << " +- " << MCerror <<endl;
  cout<<"--- kanwa ended ---"<<endl;
 
  return 0;
}
