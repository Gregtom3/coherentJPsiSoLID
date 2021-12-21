/* A class defining particles
   Declare some particles
*/

#ifndef _LPARTICLE_H_
#define _LPARTICLE_H_

#include <cmath>

#include "TRandom.h"
#include "TRandom3.h"

/* physics constant */
namespace Phys{
  const double hbar = 0.1973269718;//GeV fm
  const double c = 299792458.0;// m/s
}

/* class */
class Lparticle{
 private:
  double kmass;//Breit-Wigner mass
  double kwidth;//Breit-Wigner width
  double klife;//mean life
  int kpid;//particle id
  double kj;//spin
  TRandom3 * kran;//a pointer to random number
 public:
  Lparticle();
  Lparticle(const double mass, const double width);
  Lparticle(const double mass, const double life, const int ctr);
  double Gamma();
  double GetLife();
  double GetMass();
  int GetPID();
  double GetSpin();
  double GetWidth();
  double J();
  double M();
  double RandomM();
  double RandomM(const double Mmin, const double Mmax);
  int SetMass(const double mass);
  int SetPID(const int pid);
  int SetSpin(const double j);
  int SetWidth(const double width);
  double Tau();
  Lparticle& operator=(const Lparticle& part); 
};

/* member functions */
Lparticle::Lparticle(){
  kmass = 0;
  kwidth = 0;
  klife = 1.0;
  kran = 0;
  kj = 0;
}

Lparticle::Lparticle(const double mass, const double width){
  kmass = mass;
  kwidth = width;
  if (width > 0)
    klife = Phys::hbar * 1.0e-15 / width / Phys::c;
  else
    klife = 1.0;
  kran = 0;
  kj = 0;
}

Lparticle::Lparticle(const double mass, const double life, const int ctr){
  kmass = mass;
  klife = life;
  kwidth = Phys::hbar * 1.0e-15 / life / Phys::c;
  kran = 0;
  kj = 0;
}

double Lparticle::Gamma(){
  return GetWidth();
}

double Lparticle::GetLife(){
  return klife;
}

double Lparticle::GetMass(){
  return kmass;
}

int Lparticle::GetPID(){
  return kpid;
}

double Lparticle::GetSpin(){
  return kj;
}

double Lparticle::GetWidth(){
  return kwidth;
}

double Lparticle::J(){
  return GetSpin();
}

double Lparticle::M(){
  return GetMass();
}

double Lparticle::RandomM(){
  double Mmin = kmass - 5.0 * kwidth;
  double Mmax = kmass + 5.0 * kwidth;
  return RandomM(Mmin, Mmax);
}

double Lparticle::RandomM(const double Mmin, const double Mmax){
  if (kwidth == 0)
    return kmass;
  if (kran == 0)
    kran = new TRandom3(0);
  double mass = kmass;
  do {
    mass = kran->Uniform(Mmin, Mmax);
  } while (kran->Uniform(0.0, 1.0) > pow(kmass * kwidth, 2) / (pow(mass * mass - kmass * kmass, 2) + pow(kmass * kwidth, 2)));
  return mass;
}

int Lparticle::SetMass(const double mass){
  kmass = mass;
  return 1;
}

int Lparticle::SetPID(const int pid){
  kpid = pid;
  return 1;
}

int Lparticle::SetSpin(const double j){
  kj = j;
  return 1;
}

int Lparticle::SetWidth(const double width){
  kwidth = width;
  return 1;
}

double Lparticle::Tau(){
  return GetLife();
}

Lparticle& Lparticle::operator=(const Lparticle& part){
  kmass = part.kmass;
  kwidth = part.kwidth;
  klife = part.klife;
  kran = part.kran;
  kj = part.kj;
  return *this;
}



/* particles */
namespace PARTICLE{
  /* Gauge bosons */
  Lparticle photon(0.0, 0.0);
  Lparticle gluon(0.0, 0.0);
  /* Leptons */
  // charged leptons
  Lparticle e(0.5109989461e-3, 0.0);
  Lparticle e$m(0.5109989461e-3, 0.0);
  Lparticle e$p(0.5109989461e-3, 0.0);
  Lparticle mu(0.1056583745, 2.1969811e-6, 1);
  Lparticle mu$m(0.1056583745, 2.1969811e-6, 1);
  Lparticle mu$p(0.1056583745, 2.1969811e-6, 1);
  Lparticle tau(1.77686, 2.903e-13, 1);
  Lparticle tau$m(1.77686, 2.903e-13, 1);
  Lparticle tau$p(1.77686, 2.903e-13, 1);
  /* Mesons */
  // light mesons
  Lparticle pi(0.13957018, 2.6033e-8, 1);//charged
  Lparticle pi$p(0.13957018, 2.6033e-8, 1);
  Lparticle pi$m(0.13957018, 2.6033e-8, 1);
  Lparticle pi$0(0.1349766, 8.52e-17, 1);
  Lparticle eta(0.547862, 1.31e-6);
  Lparticle eta958(0.95778, 0.197e-3);
  Lparticle rho(0.77511, 0.1491);//charged
  Lparticle rho$p(0.77511, 0.1491);
  Lparticle rho$m(0.77511, 0.1491);
  Lparticle rho$0(0.77526, 0.1478);
  Lparticle phi(1.019461, 4.266e-3);
  Lparticle a1320(1.3183, 0.110);
  Lparticle K(0.493677, 1.2380e-8, 1);//charged
  Lparticle K$p(0.493677, 1.2380e-8, 1);
  Lparticle K$m(0.493677, 1.2380e-8, 1);
  Lparticle K$S(0.497614, 0.8954e-10, 1);
  Lparticle K$L(0.497614, 5.116e-8, 1);
  Lparticle K892(0.89166, 0.0508);//charged
  Lparticle K892$p(0.89166, 0.0508);
  Lparticle K892$m(0.89166, 0.0508);
  Lparticle K892$0(0.89581, 0.0474);
  // charmed mesons
  Lparticle D(1.86958, 1.040e-12, 1);//charged
  Lparticle D$p(1.86958, 1.040e-12, 1);
  Lparticle D$m(1.86958, 1.040e-12, 1);
  Lparticle D$0(1.86483, 0.4101e-12, 1);
  Lparticle D2010(2.01026, 83.4e-6);//charged
  Lparticle D2010$p(2.01026, 83.4e-6);
  Lparticle D2010$m(2.01026, 83.4e-6);
  Lparticle D2007(2.00685, 2.1e-3);
  // charmonium
  Lparticle etac(2.9834, 31.8e-3);
  Lparticle Jpsi(3.096900, 92.9e-6);
  Lparticle psi2S(3.68609, 286.0e-6);
  /* Baryons */
  // light baryons
  Lparticle proton(0.938272081, 0);
  Lparticle neutron(0.939565413, 880.2, 1);
  Lparticle Delta(1.232, 0.117);
  Lparticle Lambda(1.115683, 2.632e-10, 1);
  Lparticle Lambda1405(1.4051, 0.0505);
  Lparticle Lambda1520(1.5195, 0.0156);
  Lparticle Sigma$p(1.18937, 0.8018e-10, 1);
  Lparticle Sigma$0(1.192642, 7.4e-20, 1);
  Lparticle Sigma$m(1.197449, 1.479e-10, 1);
  Lparticle Sigma1385$p(1.38280, 0.0360);
  Lparticle Sigma1385$0(1.3837, 0.036);
  Lparticle Sigma1385$m(1.3872, 0.0394);
  // charmed baryons
  Lparticle Lambdac(2.28646, 0.2e-12, 1);
  Lparticle Sigmac$pp(2.45397, 1.89e-3);
  Lparticle Sigmac$p(2.4529, 4.6e-3);
  Lparticle Sigmac$0(2.45375, 1.83e-3);
  Lparticle Sigmac2520$pp(2.51841, 14.78e-3);
  Lparticle Sigmac2520$p(2.5175, 17.0e-3);
  Lparticle Sigmac2520$0(2.51848, 15.3e-3);
  //bottomonium
  Lparticle upsilon1S(9.46030, 0.00005402);
}


#endif
