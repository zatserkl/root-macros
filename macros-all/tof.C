#include <TROOT.h>
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>

using std::cout;     using std::endl;

Double_t pdelay(Double_t T, Double_t L)    // kinetic energy T in MeV, base L in m
{
   // delay of a proton with kinetic energy T wrt positron on base L

   const Double_t c = 0.299792458;  // m/ns
   const Double_t M = 938.27;       // MeV

   Double_t beta_p = TMath::Sqrt(1 + 2*M/T) / (1 + M/T);
   Double_t beta_e = 1;

   Double_t t_p = L / (beta_p*c);   // ns
   Double_t t_e = L / (beta_e*c);   // ns

   Double_t pdelay = t_p - t_e;
   return pdelay;
}

void pdelay_plot(Double_t L, Double_t T1=300, Double_t T2=900)    // base L in m, kinetic energy T in MeV
{
   Double_t dT = 1;     // MeV
   Int_t np = 0;

   TGraph* gr_delay_e = new TGraph(10000);
   gr_delay_e->SetNameTitle("gr_pdelay", Form("Delay of protons wrt electrons for base L = %0.2f m;proton's kinetic energy T, MeV;delay, ns", L));
   gr_delay_e->SetMarkerStyle(7);
   gr_delay_e->SetMarkerColor(2);

   TGraph* gr_delay_pi = new TGraph(10000);
   gr_delay_pi->SetNameTitle("gr_pdelay", Form("Delay of protons wrt pions for base L = %0.2f m;proton's kinetic energy T, MeV;delay, ns", L));
   gr_delay_pi->SetMarkerStyle(7);
   gr_delay_pi->SetMarkerColor(4);

   Double_t T = T1;
   while (T < T2) {
      gr_delay_e->SetPoint(np, T, pdelay(T, L));
      gr_delay_e->Set(np+1);
      gr_delay_pi->SetPoint(np, T, pdelay(T, L));
      gr_delay_pi->Set(np+1);
      np++;
      T += dT;
   }

   new TCanvas;
   gr_delay_e->Draw("ap");

   new TCanvas;
   gr_delay_pi->Draw("ap");
}

////////////////////////////////////////////////
//
// new approach:
//
// 1) find beta of the proton with kinetic energy T
// 2) find momentum of the proton with kinetic energy T
// 3) find beta of the beam particle with momentum p and mass m (has the same momentum as the proton)
//
////////////////////////////////////////////////

//-- kinematics

Double_t beta_alpha(Double_t T)   // beta of a proton with kinetic energy T, MeV
{
   const Double_t M = 4*938.27;             // MeV
   Double_t beta = TMath::Sqrt(2.*T/M) * TMath::Sqrt(1.+0.5*T/M) / (1.+T/M);
   return beta;
}

Double_t beta_proton(Double_t T)   // beta of a proton with kinetic energy T, MeV
{
   const Double_t M = 938.27;             // MeV
   Double_t beta = TMath::Sqrt(2.*T/M) * TMath::Sqrt(1.+0.5*T/M) / (1.+T/M);
   return beta;
}

Double_t Pproton(Double_t T)       // momentum of a proton with kinetic energy T
{
   const Double_t M = 938.27;             // MeV
   Double_t p = TMath::Sqrt(2.*M*T) * TMath::Sqrt(1.+0.5*T/M);
   return p;
}

void plot_Proton(Double_t T1=300, Double_t T2=900)    // kinetic energy and momentum in MeV
{
   // TGraph* gr_p_T = new TGraph(10000);
   TGraph* gr_p_T = new TGraph(10);
   gr_p_T->SetNameTitle("gr_Tp", "Momentum vs kinetic energy for proton;proton's kinetic energy T, MeV;proton's momentum p, MeV/c");
   gr_p_T->SetMarkerStyle(6);
   gr_p_T->SetMarkerColor(2);
   gr_p_T->SetLineColor(2);

   Double_t T = T1;
   Double_t dT = 1.;    // MeV
   Int_t np = 0;

   while (T < T2) {
      Double_t p = Pproton(T);
      gr_p_T->SetPoint(np, T, p);
      np++;
      T += dT;
   }

   new TCanvas;
   gr_p_T->Draw("ap");
}

Double_t Tproton(Double_t p)       // proton kinetic energy from the proton momentum
{
   const Double_t M = 938.27;             // MeV
   Double_t E = TMath::Sqrt(p*p + M*M);
   Double_t T = E - M;
   return T;
}

Double_t beta_particle(Double_t p, Double_t m)  // beta of a particle with momentum p and mass m
{
   Double_t beta = p / TMath::Sqrt(p*p + m*m);
   return beta;
}

void delay_plot(Double_t L, Double_t T1=300, Double_t T2=900)    // base L in m, kinetic energy T in MeV
{
   // delay of a proton wrt electron and pion

   const Double_t c = 0.299792458;  // m/ns
   const Double_t m_e = 0.511;      // MeV
   const Double_t m_pi = 206.0;     // MeV

   Double_t dT = 1;     // MeV
   Int_t np = 0;

   TGraph* gr_delay_e = new TGraph(10000);
   gr_delay_e->SetNameTitle("gr_pdelay", Form("Delay of protons wrt electrons for base L = %0.2f m;proton's kinetic energy T, MeV;delay, ns", L));
   gr_delay_e->SetMarkerStyle(7);
   gr_delay_e->SetMarkerColor(2);

   TGraph* gr_delay_pi = new TGraph(10000);
   gr_delay_pi->SetNameTitle("gr_pdelay", Form("Delay of protons wrt pions for base L = %0.2f m;proton's kinetic energy T, MeV;delay, ns", L));
   gr_delay_pi->SetMarkerStyle(7);
   gr_delay_pi->SetMarkerColor(4);

   Double_t T = T1;
   while (T < T2) {
      Double_t p = Pproton(T);                     // proton momentum
      Double_t beta_p = beta_proton(T);            // beta of proton

      Double_t beta_e = beta_particle(p, m_e);     // beta of electron with momentum equal to the proton momentum
      Double_t delay_e = (L/c) * (1./beta_p - 1./beta_e);

      Double_t beta_pi = beta_particle(p, m_pi);   // beta of pion with momentum equal to the proton momentum
      Double_t delay_pi = (L/c) * (1./beta_p - 1./beta_pi);

      gr_delay_e->SetPoint(np, T, delay_e);
      gr_delay_e->Set(np+1);
      gr_delay_pi->SetPoint(np, T, delay_pi);
      gr_delay_pi->Set(np+1);
      np++;
      T += dT;
   }

   new TCanvas;
   gr_delay_e->Draw("ap");

   new TCanvas;
   gr_delay_pi->Draw("ap");
}

void delay_plot_p(Double_t L=8.70, Double_t p1=4000, Double_t p2=8000)    // base L in m, proton momentum p in MeV/c
{
   // delay of a proton wrt electron and pion vs proton momentum in MeV/c. Fermilab beam test: L = 8.7 m

   const Double_t c = 0.299792458;  // m/ns
   const Double_t m_e = 0.511;      // MeV
   const Double_t m_pi = 206.0;     // MeV
   const Double_t M = 938.27;       // MeV

   Double_t dp = 1;     // MeV/c
   Int_t np = 0;

   TGraph* gr_delay_e = new TGraph(10000);
   gr_delay_e->SetNameTitle("gr_pdelay", Form("Delay of protons wrt electrons for base L = %0.2f m;proton's momentum, MeV/c;delay, ns", L));
   gr_delay_e->SetMarkerStyle(7);
   gr_delay_e->SetMarkerColor(2);

   TGraph* gr_delay_pi = new TGraph(10000);
   gr_delay_pi->SetNameTitle("gr_pdelay", Form("Delay of protons wrt pions for base L = %0.2f m;proton's momentum, MeV/c;delay, ns", L));
   gr_delay_pi->SetMarkerStyle(7);
   gr_delay_pi->SetMarkerColor(4);

   Double_t p = p1;
   while (p < p2) {
      Double_t beta_p = beta_particle(p, M);             // beta of proton

      Double_t beta_e = beta_particle(p, m_e);           // beta of an electron with momentum equal to the proton momentum
      Double_t delay_e = (L/c) * (1./beta_p - 1./beta_e);

      Double_t beta_pi = beta_particle(p, m_pi);         // beta of a pion with momentum equal to the proton momentum
      Double_t delay_pi = (L/c) * (1./beta_p - 1./beta_pi);

      gr_delay_e->SetPoint(np, p, delay_e);
      gr_delay_e->Set(np+1);
      gr_delay_pi->SetPoint(np, p, delay_pi);
      gr_delay_pi->Set(np+1);
      np++;
      p += dp;
   }

   new TCanvas;
   gr_delay_e->Draw("ap");

   new TCanvas;
   gr_delay_pi->Draw("ap");
}

void delay_plot_p_check(Double_t L=8.70, Double_t p1=4000, Double_t p2=8000)    // base L in m, proton momentum p in MeV/c
{
   // delay of a proton wrt electron and pion vs proton momentum in MeV/c. Fermilab beam test: L = 8.7 m

   const Double_t c = 0.299792458;  // m/ns
   const Double_t m_e = 0.511;      // MeV
   const Double_t m_pi = 206.0;     // MeV
   const Double_t M = 938.27;       // MeV

   Double_t dp = 1;     // MeV/c
   Int_t np = 0;

   TGraph* gr_delay_e = new TGraph(10000);
   gr_delay_e->SetNameTitle("gr_pdelay", Form("Delay of protons wrt electrons for base L = %0.2f m;proton's momentum, MeV/c;delay, ns", L));
   gr_delay_e->SetMarkerStyle(7);
   gr_delay_e->SetMarkerColor(2);

   TGraph* gr_delay_pi = new TGraph(10000);
   gr_delay_pi->SetNameTitle("gr_pdelay", Form("Delay of protons wrt pions for base L = %0.2f m;proton's momentum, MeV/c;delay, ns", L));
   gr_delay_pi->SetMarkerStyle(7);
   gr_delay_pi->SetMarkerColor(4);

   Double_t p = p1;
   while (p < p2) {
      // calculate proton total energy

      Double_t E_total = TMath::Sqrt(p*p + M*M);
      Double_t T = E_total - M;

      // proton momentum from the kinetic energy (to check the formulae)

      Double_t proton_momentum = Pproton(T);

      Double_t beta_p = beta_particle(proton_momentum, M);             // beta of proton

      Double_t beta_e = beta_particle(proton_momentum, m_e);           // beta of an electron with momentum equal to the proton momentum
      Double_t delay_e = (L/c) * (1./beta_p - 1./beta_e);

      Double_t beta_pi = beta_particle(proton_momentum, m_pi);         // beta of a pion with momentum equal to the proton momentum
      Double_t delay_pi = (L/c) * (1./beta_p - 1./beta_pi);

      gr_delay_e->SetPoint(np, proton_momentum, delay_e);
      gr_delay_e->Set(np+1);
      gr_delay_pi->SetPoint(np, proton_momentum, delay_pi);
      gr_delay_pi->Set(np+1);
      np++;
      p += dp;
   }

   new TCanvas;
   gr_delay_e->Draw("ap");

   new TCanvas;
   gr_delay_pi->Draw("ap");
}

//-- quartz (fused silica) and Photek240

class Quartz {
public:
   TGraph* gr_trans;
   Int_t n_trans;
   Quartz() {
      gr_trans = new TGraph(1000);
      n_trans = 0;
      // data: transparency of 1 cm of quartz
      gr_trans->SetPoint(n_trans++, 150, 0);
      gr_trans->SetPoint(n_trans++, 160, 0.30);
      gr_trans->SetPoint(n_trans++, 170, 0.57);
      gr_trans->SetPoint(n_trans++, 180, 0.74);
      gr_trans->SetPoint(n_trans++, 190, 0.83);
      gr_trans->SetPoint(n_trans++, 200, 0.87);
      gr_trans->SetPoint(n_trans++, 210, 0.89);
      gr_trans->SetPoint(n_trans++, 220, 0.90);
      gr_trans->SetPoint(n_trans++, 230, 0.90);
      gr_trans->SetPoint(n_trans++, 240, 0.90);
      gr_trans->SetPoint(n_trans++, 250, 0.900);
      gr_trans->SetPoint(n_trans++, 260, 0.900);
      gr_trans->SetPoint(n_trans++, 270, 0.901);
      gr_trans->SetPoint(n_trans++, 280, 0.901);
      gr_trans->SetPoint(n_trans++, 290, 0.902);
      gr_trans->SetPoint(n_trans++, 300, 0.902);
      gr_trans->SetPoint(n_trans++, 310, 0.902);
      gr_trans->SetPoint(n_trans++, 320, 0.903);
      gr_trans->SetPoint(n_trans++, 330, 0.903);
      gr_trans->SetPoint(n_trans++, 340, 0.904);
      gr_trans->SetPoint(n_trans++, 350, 0.904);
      gr_trans->SetPoint(n_trans++, 360, 0.904);
      gr_trans->SetPoint(n_trans++, 370, 0.905);
      gr_trans->SetPoint(n_trans++, 380, 0.905);
      gr_trans->SetPoint(n_trans++, 390, 0.906);
      gr_trans->SetPoint(n_trans++, 400, 0.906);
      gr_trans->SetPoint(n_trans++, 410, 0.906);
      gr_trans->SetPoint(n_trans++, 420, 0.907);
      gr_trans->SetPoint(n_trans++, 430, 0.907);
      gr_trans->SetPoint(n_trans++, 440, 0.908);
      gr_trans->SetPoint(n_trans++, 450, 0.908);
      gr_trans->SetPoint(n_trans++, 460, 0.908);
      gr_trans->SetPoint(n_trans++, 470, 0.909);
      gr_trans->SetPoint(n_trans++, 480, 0.909);
      gr_trans->SetPoint(n_trans++, 490, 0.910);
      gr_trans->SetPoint(n_trans++, 500, 0.910);
      gr_trans->SetPoint(n_trans++, 510, 0.910);
      gr_trans->SetPoint(n_trans++, 520, 0.911);
      gr_trans->SetPoint(n_trans++, 530, 0.911);
      gr_trans->SetPoint(n_trans++, 540, 0.912);
      gr_trans->SetPoint(n_trans++, 550, 0.912);
      gr_trans->SetPoint(n_trans++, 560, 0.912);
      gr_trans->SetPoint(n_trans++, 570, 0.913);
      gr_trans->SetPoint(n_trans++, 580, 0.913);
      gr_trans->SetPoint(n_trans++, 590, 0.914);
      gr_trans->SetPoint(n_trans++, 600, 0.914);
      gr_trans->SetPoint(n_trans++, 610, 0.914);
      gr_trans->SetPoint(n_trans++, 620, 0.915);
      gr_trans->SetPoint(n_trans++, 630, 0.915);
      gr_trans->SetPoint(n_trans++, 640, 0.916);
      gr_trans->SetPoint(n_trans++, 650, 0.916);
      gr_trans->SetPoint(n_trans++, 660, 0.916);
      gr_trans->SetPoint(n_trans++, 670, 0.917);
      gr_trans->SetPoint(n_trans++, 680, 0.917);
      gr_trans->SetPoint(n_trans++, 690, 0.918);
      gr_trans->SetPoint(n_trans++, 700, 0.918);
      gr_trans->SetPoint(n_trans++, 710, 0.918);
      gr_trans->SetPoint(n_trans++, 720, 0.919);
      gr_trans->SetPoint(n_trans++, 730, 0.919);
      gr_trans->SetPoint(n_trans++, 740, 0.920);
      gr_trans->SetPoint(n_trans++, 750, 0.920);
      gr_trans->SetPoint(n_trans++, 760, 0.920);
      gr_trans->SetPoint(n_trans++, 770, 0.921);
      gr_trans->SetPoint(n_trans++, 780, 0.921);
      gr_trans->SetPoint(n_trans++, 790, 0.922);
      gr_trans->SetPoint(n_trans++, 800, 0.922);
      gr_trans->SetPoint(n_trans++, 810, 0.922);
      gr_trans->SetPoint(n_trans++, 820, 0.923);
      gr_trans->SetPoint(n_trans++, 830, 0.923);
      gr_trans->SetPoint(n_trans++, 840, 0.924);
      gr_trans->SetPoint(n_trans++, 850, 0.924);
      gr_trans->SetPoint(n_trans++, 860, 0.924);
      gr_trans->SetPoint(n_trans++, 870, 0.925);
      gr_trans->SetPoint(n_trans++, 880, 0.925);
      gr_trans->SetPoint(n_trans++, 890, 0.926);
      gr_trans->SetPoint(n_trans++, 900, 0.926);
      gr_trans->SetPoint(n_trans++, 910, 0.926);
      gr_trans->SetPoint(n_trans++, 920, 0.927);
      gr_trans->SetPoint(n_trans++, 930, 0.927);
      gr_trans->SetPoint(n_trans++, 940, 0.928);
      gr_trans->SetPoint(n_trans++, 950, 0.928);
      gr_trans->SetPoint(n_trans++, 960, 0.928);
      gr_trans->SetPoint(n_trans++, 970, 0.929);
      gr_trans->SetPoint(n_trans++, 980, 0.929);
      gr_trans->SetPoint(n_trans++, 990, 0.930);
      gr_trans->SetPoint(n_trans++, 1000,0.930);
      gr_trans->Set(n_trans);
   }

   static Double_t nref(Double_t lambda)
   {
      lambda *= 1e-3;   // convert to um

      if (lambda < 0.15 or lambda > 1) return 0;

      const Double_t lambda1 = 0.0684043;
      const Double_t lambda2 = 0.1162414;
      const Double_t lambda3 = 9.896161;

      // NB: the n for the wavelength close to the Sellmeier coefficients should be calculated using other model, like Helmholtz

      Double_t denom1 = lambda*lambda - lambda1*lambda1;
      Double_t denom2 = lambda*lambda - lambda2*lambda2;
      Double_t denom3 = lambda*lambda - lambda3*lambda3;

      // const Double_t eps = 1e-7;
      // if (denom1 < eps or denom2 < eps or denom3 < eps) return 0;

      Double_t n2 = 1
         + 0.6961663*lambda*lambda / denom1
         + 0.4079426*lambda*lambda / denom2
         + 0.8974794*lambda*lambda / denom3
         ;

      Double_t n = n2 > 0? TMath::Sqrt(n2): 0;
      return n;
   }
   static Double_t nref(Double_t* xx, Double_t*)
   {
      Double_t lambda = *xx;
      return nref(lambda);
   }
   ~Quartz() {
      delete gr_trans;     gr_trans = 0;
   }
};

class Quartz_trans: public Quartz {
public:
   Quartz_trans(): Quartz() {}
   double operator ()(double* xx, double*) {
      double& x = *xx;
      if (x < gr_trans->GetX()[0] or x > gr_trans->GetX()[gr_trans->GetN()-1]) return 0;
      return gr_trans->Eval(x);
   }
};

void Quartz_nref()
{
   TF1* fnref = new TF1("fnref", Quartz::nref, 200, 1000, 0);
   fnref->SetTitle("Quartz index of refraction;#lambda, nm");

   new TCanvas;
   fnref->Draw();
}

void test_quartz()
{
   Quartz_trans* quartz_trans = new Quartz_trans;
   TF1* ftrans = new TF1("ftrans", quartz_trans, 100, 1000, 0, "Quartz_trans");
   ftrans->SetTitle("Quartz transparency;#lambda, nm;transparency, %");

   new TCanvas;
   ftrans->Draw();
}

class Photek240 {
public:
   Double_t integral_sensitivity;   // uA/lm
   Int_t np;
   TGraph* gr_sens;
   TGraph* gr_QE;
   Photek240() {
      integral_sensitivity = 145.;
      // NB: two starting zeros and two ending zeros provide correct interpolation value, 0 out of range
      Double_t lambda[] =        {100.,   200.,   214.,    254.,    270.,    280.,    290.,    300.,    350.,    400.,    450.,    500.,    532.,    550.,    650.,    750.,    800.,    850.,    900.,    925.,  950.};
      Double_t sensitivity[] =   {0.,     0.,     30.39,   42.29,   47.07,   45.61,   46.81,   46.32,   52.54,   44.57,   51.725,  37.821,  38.532,  36.482,  38.265,  4.828,   0.78,    0.10,    0.03,    0.,    0.};
      Double_t QE[] =            {0.,     0.,     17.61,   20.64,   21.62,   20.20,   20.02,   19.14,   18.61,   13.82,   14.253,  9.380,   8.981,   8.225,   7.300,   0.80,    0.12,    0.01,    0.00,    0.,    0.};
      np = sizeof(lambda) / sizeof(Double_t);
      gr_sens = new TGraph(np, lambda, sensitivity);
      gr_sens->SetNameTitle("gr_sens","Photek240 Spectral Response;Wavelength, nm;Sensitivity, mA/W");
      gr_sens->SetMarkerStyle(7);
      gr_sens->SetMarkerColor(2);
      gr_sens->SetLineColor(2);

      gr_QE = new TGraph(np, lambda, QE);
      gr_QE->SetNameTitle("gr_QE", "Photek240 Quantum Efficiency;Wavelength, nm;QE, %");
      gr_QE->SetMarkerStyle(7);
      gr_QE->SetMarkerColor(4);
      gr_QE->SetLineColor(4);
   }
   ~Photek240() {
      delete gr_sens;   gr_sens = 0;
      delete gr_QE;     gr_QE = 0;
   }
   void plot_Photek240_gr_QE(const char* options="ap") const {
      gr_QE->Draw(options);
   }
};

class Photek240_sensitivity: public Photek240 {
public:
   Photek240_sensitivity(): Photek240() {}
   double operator ()(double* xx, double*) {
      double& x = *xx;
      if (x < gr_sens->GetX()[0] or x > gr_sens->GetX()[gr_sens->GetN()-1]) return 0;
      return gr_sens->Eval(x);
   }
};

class Photek240_QE: public Photek240 {
public:
   Photek240_QE(): Photek240() {}
   double operator ()(double* xx, double*) {
      double& x = *xx;
      if (x < gr_QE->GetX()[0] or x > gr_QE->GetX()[gr_QE->GetN()-1]) return 0;
      return gr_QE->Eval(x);
   }
};

void test_Photek240()
{
   Photek240_sensitivity* photek240_sensitivity = new Photek240_sensitivity;
   TF1* fsens = new TF1("fsens", photek240_sensitivity, 100, 1000, 0, "Photek240_sensitivity");
   fsens->SetTitle("Photek 240 Sensitivity;#lambda, nm;sensitivity, mA/W");

   new TCanvas;
   fsens->Draw();

   Photek240_QE* photek240_QE = new Photek240_QE;
   TF1* fQE = new TF1("fQE", photek240_QE, 100, 1000, 0, "Photek240_QE");
   fQE->SetTitle("Photek 240 Quantum Efficiency;#lambda, nm;QE, %");
   fQE->SetLineColor(4);

   new TCanvas;
   fQE->Draw();
}

class Nphotons: public Photek240_QE, public Quartz {
public:
   Double_t beta;
   Double_t length;
   Nphotons(Double_t beta_, Double_t length_): Photek240_QE(), Quartz()
      , beta(beta_)
      , length(length_)
   {}
   double operator ()(double* xx, double*) {
      double& lambda = *xx;

      const double alpha = 1./137.;
      const double nm2cm = 1e9 * 1e-2;    // to use the lambda in nm and the length in cm
      const double eps = 1e-7;

      double beta_n = beta * nref(*xx);
      //cout<< beta_n <<endl;
      double Nphot = 0;
      if (beta_n > eps) Nphot = length * (TMath::TwoPi() * alpha) * (1. - 1./(beta_n*beta_n)) * 1./(lambda*lambda) * nm2cm;
      return Nphot;
   }
};

void test_Nphotons_heap(Double_t beta=1., Double_t length=0.9, Double_t lambda1=220, Double_t lambda2=650)
{
   Nphotons* nphot = new Nphotons(beta, length);
   TF1* fnphot = new TF1("fnphot", nphot, lambda1, lambda2, 0, "Nphotons");
   fnphot->SetTitle("Integrand for Photek 240;#lambda, nm;#frac{d^{2}N}{dxd#lambda}");

   new TCanvas;
   fnphot->DrawCopy();

   Double_t Nphot = fnphot->Integral(lambda1, lambda2);
   cout<< "Nphot = " << Nphot <<endl;
}

void test_Nphotons(Double_t beta=1., Double_t length=0.9, Double_t lambda1=220, Double_t lambda2=650)
{
   Nphotons nphot(beta, length);
   TF1 fnphot("fnphot", &nphot, lambda1, lambda2, 0, "Nphotons");
   fnphot.SetTitle("Integrand for Photek 240;#lambda, nm;d^{2}N/dxd#lambda");

   new TCanvas;
   fnphot.DrawCopy();

   Double_t Nphot = fnphot.Integral(lambda1, lambda2);
   cout<< "Nphot = " << Nphot <<endl;
}

class Nphotoelectrons: public Photek240, public Quartz {
public:
   Double_t beta;
   Double_t length;
   Nphotoelectrons(Double_t beta_, Double_t length_): Photek240(), Quartz()
      , beta(beta_)
      , length(length_)
   {}
   double operator ()(double* xx, double*) {
      double& lambda = *xx;

      const double alpha = 1./137.;
      const double nm2cm = 1e9 * 1e-2;    // to use the lambda in nm and the length in cm
      const double eps = 1e-7;

      double beta_n = beta * nref(*xx);
      //cout<< beta_n <<endl;

      double Npe = 0;
      if (beta_n > eps) {
         double Nphot = 0;
         Nphot = length * (TMath::TwoPi() * alpha) * (1. - 1./(beta_n*beta_n)) * 1./(lambda*lambda) * nm2cm;
         double qe = 0;
         qe = 0.01 * gr_QE->Eval(lambda);    // 0.01: convert percents to fraction
         Npe = Nphot * qe;
      }

      return Npe;
   }
};

void test_Nphotoelectrons(Double_t beta=1., Double_t length=0.9, Double_t lambda1=220, Double_t lambda2=650)
{
   Nphotoelectrons* npe = new Nphotoelectrons(beta, length);
   TF1 fnpe("fnpe", npe, lambda1, lambda2, 0, "Nphotoelectrons");
   fnpe.SetTitle("Integrand for Photek 240;#lambda, nm;d^{2}N/dxd#lambda");

   new TCanvas;
   fnpe.DrawCopy();

   Double_t Npe = fnpe.Integral(lambda1, lambda2);
   cout<< "Npe = " << Npe <<endl;
}

Double_t Npe_quartz(Double_t beta=1., Double_t length=0.9, Double_t lambda1=220, Double_t lambda2=650)
{
   Nphotoelectrons* npe = new Nphotoelectrons(beta, length);
   TF1 fnpe("fnpe", npe, lambda1, lambda2, 0, "Nphotoelectrons");
   fnpe.SetTitle("Integrand for Photek 240;#lambda, nm;d^{2}N/dxd#lambda");

   Double_t Npe = fnpe.Integral(lambda1, lambda2);
   return Npe;
}

Double_t Npe_quartz_momentum(Double_t p, Double_t mass, Double_t length=0.9, Double_t lambda1=220, Double_t lambda2=650)
{
   // length in cm

   Double_t beta = beta_particle(p, mass);

   Nphotoelectrons* npe = new Nphotoelectrons(beta, length);
   TF1 fnpe("fnpe", npe, lambda1, lambda2, 0, "Nphotoelectrons");
   fnpe.SetTitle("Integrand for Photek 240;#lambda, nm;d^{2}N/dxd#lambda");

   Double_t Npe = fnpe.Integral(lambda1, lambda2);
   return Npe;
}

Double_t Npe_Photek240(Double_t p, Double_t mass, Double_t lambda1=220, Double_t lambda2=650)
{
   Double_t length = 0.9;     // cm

   Double_t beta = beta_particle(p, mass);

   Nphotoelectrons* npe = new Nphotoelectrons(beta, length);
   TF1 fnpe("fnpe", npe, lambda1, lambda2, 0, "Nphotoelectrons");
   fnpe.SetTitle("Integrand for Photek 240;#lambda, nm;d^{2}N/dxd#lambda");

   Double_t Npe = fnpe.Integral(lambda1, lambda2);
   return Npe;
}

/////////////////////////////////////////////////////////////////
//
// Cherenkov light
//
/////////////////////////////////////////////////////////////////

Double_t Nphotons_nref(Double_t beta, Double_t n, Double_t lambda1, Double_t lambda2, Double_t length)
{
   // length in cm

   const Double_t alpha = 1./137.;
   const Double_t nm2cm = 1e7;

   Double_t beta_n = beta*n;
   Double_t term = 1. - 1./(beta_n*beta_n);
   if (term < 0) term = 0;
   Double_t nphotons = TMath::TwoPi()*alpha*length*term*(1./lambda1 - 1./lambda2)*nm2cm;
   return nphotons;
}

Double_t Nphotons_nref_momentum(Double_t p, Double_t mass, Double_t n, Double_t lambda1, Double_t lambda2, Double_t length)
{
   // length in cm

   Double_t beta = beta_particle(p, mass);
   return Nphotons_nref(beta, n, lambda1, lambda2, length);
}

void Nphotons_momentum_plot(Double_t nref=1.5, Double_t lambda1=200, Double_t lambda2=650, Double_t length=1.)  // lambda in nm, length in cm
{
   // NB: Hanz Wenzel plot corresponds to approx. lambda range 300-450 (n = 1.5)

   const Double_t m_p = 938.27;        // MeV
   const Double_t m_e = 0.511;         // MeV
   const Double_t m_pi = 139.6;        // MeV
   const Double_t m_mu = 105.7;        // MeV
   const Double_t m_K = 493.7;         // MeV

   Int_t np_p = 0;
   TGraph* gr_p = new TGraph(0);
   gr_p->SetNameTitle("gr_p", Form("Cherenkov radiation of proton. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_p->SetFillColor(0);
   gr_p->SetMarkerStyle(1);
   gr_p->SetMarkerColor(2);
   gr_p->SetLineColor(2);

   Int_t np_e = 0;
   TGraph* gr_e = new TGraph(0);
   gr_e->SetNameTitle("gr_e", Form("Cherenkov radiation of electron. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_e->SetFillColor(0);
   gr_e->SetMarkerStyle(1);
   gr_e->SetMarkerColor(4);
   gr_e->SetLineColor(4);

   Int_t np_mu = 0;
   TGraph* gr_mu = new TGraph(0);
   gr_mu->SetNameTitle("gr_mu", Form("Cherenkov radiation of #mu. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_mu->SetFillColor(0);
   gr_mu->SetMarkerStyle(1);
   gr_mu->SetMarkerColor(6);
   gr_mu->SetLineColor(6);

   Int_t np_pi = 0;
   TGraph* gr_pi = new TGraph(0);
   gr_pi->SetNameTitle("gr_pi", Form("Cherenkov radiation of #pi. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_pi->SetFillColor(0);
   gr_pi->SetMarkerStyle(1);
   gr_pi->SetMarkerColor(1);
   gr_pi->SetLineColor(1);

   Int_t np_K = 0;
   TGraph* gr_K = new TGraph(0);
   gr_K->SetNameTitle("gr_K", Form("Cherenkov radiation of #font[12]{K}. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_K->SetFillColor(0);
   gr_K->SetMarkerStyle(1);
   gr_K->SetMarkerColor(8);
   gr_K->SetLineColor(8);

   Double_t p1 = 10.;   // MeV/c
   Double_t p2 = 5000.; // MeV/c
   Double_t p = p1;
   Double_t dp = 10.;

   while (p < p2)
   {
      Double_t nphotons_p = Nphotons_nref_momentum(p, m_p, nref, lambda1, lambda2, length);
      gr_p->SetPoint(np_p, p, nphotons_p);
      gr_p->Set(++np_p);

      Double_t nphotons_mu = Nphotons_nref_momentum(p, m_mu, nref, lambda1, lambda2, length);
      gr_mu->SetPoint(np_mu, p, nphotons_mu);
      gr_mu->Set(++np_mu);

      Double_t nphotons_pi = Nphotons_nref_momentum(p, m_pi, nref, lambda1, lambda2, length);
      gr_pi->SetPoint(np_pi, p, nphotons_pi);
      gr_pi->Set(++np_pi);

      Double_t nphotons_e = Nphotons_nref_momentum(p, m_e, nref, lambda1, lambda2, length);
      gr_e->SetPoint(np_e, p, nphotons_e);
      gr_e->Set(++np_e);

      Double_t nphotons_K = Nphotons_nref_momentum(p, m_K, nref, lambda1, lambda2, length);
      gr_K->SetPoint(np_K, p, nphotons_K);
      gr_K->Set(++np_K);

      p += dp;
   }

   TCanvas* can = new TCanvas;
   can->DrawFrame(0,1, 1.10*p2, 1.10*gr_e->GetY()[np_e-1], Form("Cherenkov radiation n=%0.2f, #lambda=%0.0f-%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", nref,lambda1,lambda2,length));
   gr_p->Draw("pl");
   gr_pi->Draw("pl");
   gr_mu->Draw("pl");
   gr_e->Draw("pl");
   gr_K->Draw("pl");

   TLegend* legend = new TLegend(0.70,0.20, 0.78,0.70);
   legend->SetFillColor(0);
   legend->SetFillStyle(1001);
   legend->AddEntry(gr_e, " #font[12]{e}");
   legend->AddEntry(gr_mu, " #font[12]{#mu}");
   legend->AddEntry(gr_pi, " #font[12]{#pi}");
   legend->AddEntry(gr_K, " #font[12]{K}");
   legend->AddEntry(gr_p, " #font[12]{p}");

   legend->Draw();
}

// Double_t Nphotons_quartz(Double_t beta, Double_t lambda1, Double_t lambda2, Double_t length)
// {
//    const Double_t alpha = 1./137.;
//    const Double_t nm2cm = 1e7;
// 
//    //---------------------------------------Double_t beta_n = beta * Quartz::nref();
//    Double_t term = 1. - 1./(beta_n*beta_n);
//    if (term < 0) term = 0;
//    Double_t nphotons = TMath::TwoPi()*alpha*length*term*(1./lambda1 - 1./lambda2)*nm2cm;
//    return nphotons;
// }

void Nphotons_quartz_momentum_plot(Double_t nref=1.5, Double_t lambda1=200, Double_t lambda2=650, Double_t length=1.)  // lambda in nm, length in cm
{
   // NB: Hanz Wenzel plot corresponds to approx. lambda range 300-450 (n = 1.5)

   const Double_t m_p = 938.27;        // MeV
   const Double_t m_e = 0.511;         // MeV
   const Double_t m_pi = 139.6;        // MeV
   const Double_t m_mu = 105.7;        // MeV
   const Double_t m_K = 493.7;         // MeV

   Int_t np_p = 0;
   TGraph* gr_p = new TGraph(0);
   gr_p->SetNameTitle("gr_p", Form("Cherenkov radiation of proton. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_p->SetFillColor(0);
   gr_p->SetMarkerStyle(1);
   gr_p->SetMarkerColor(2);
   gr_p->SetLineColor(2);

   Int_t np_e = 0;
   TGraph* gr_e = new TGraph(0);
   gr_e->SetNameTitle("gr_e", Form("Cherenkov radiation of electron. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_e->SetFillColor(0);
   gr_e->SetMarkerStyle(1);
   gr_e->SetMarkerColor(4);
   gr_e->SetLineColor(4);

   Int_t np_mu = 0;
   TGraph* gr_mu = new TGraph(0);
   gr_mu->SetNameTitle("gr_mu", Form("Cherenkov radiation of #mu. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_mu->SetFillColor(0);
   gr_mu->SetMarkerStyle(1);
   gr_mu->SetMarkerColor(6);
   gr_mu->SetLineColor(6);

   Int_t np_pi = 0;
   TGraph* gr_pi = new TGraph(0);
   gr_pi->SetNameTitle("gr_pi", Form("Cherenkov radiation of #pi. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_pi->SetFillColor(0);
   gr_pi->SetMarkerStyle(1);
   gr_pi->SetMarkerColor(1);
   gr_pi->SetLineColor(1);

   Int_t np_K = 0;
   TGraph* gr_K = new TGraph(0);
   gr_K->SetNameTitle("gr_K", Form("Cherenkov radiation of #font[12]{K}. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", lambda1,lambda2,length));
   gr_K->SetFillColor(0);
   gr_K->SetMarkerStyle(1);
   gr_K->SetMarkerColor(8);
   gr_K->SetLineColor(8);

   Double_t p1 = 10.;   // MeV/c
   Double_t p2 = 5000.; // MeV/c
   Double_t p = p1;
   Double_t dp = 10.;

   while (p < p2)
   {
      Double_t nphotons_p = Nphotons_nref_momentum(p, m_p, nref, lambda1, lambda2, length);
      gr_p->SetPoint(np_p, p, nphotons_p);
      gr_p->Set(++np_p);

      Double_t nphotons_mu = Nphotons_nref_momentum(p, m_mu, nref, lambda1, lambda2, length);
      gr_mu->SetPoint(np_mu, p, nphotons_mu);
      gr_mu->Set(++np_mu);

      Double_t nphotons_pi = Nphotons_nref_momentum(p, m_pi, nref, lambda1, lambda2, length);
      gr_pi->SetPoint(np_pi, p, nphotons_pi);
      gr_pi->Set(++np_pi);

      Double_t nphotons_e = Nphotons_nref_momentum(p, m_e, nref, lambda1, lambda2, length);
      gr_e->SetPoint(np_e, p, nphotons_e);
      gr_e->Set(++np_e);

      Double_t nphotons_K = Nphotons_nref_momentum(p, m_K, nref, lambda1, lambda2, length);
      gr_K->SetPoint(np_K, p, nphotons_K);
      gr_K->Set(++np_K);

      p += dp;
   }

   TCanvas* can = new TCanvas;
   can->DrawFrame(0,1, 1.10*p2, 1.10*gr_e->GetY()[np_e-1], Form("Cherenkov radiation n=%0.2f, #lambda=%0.0f-%0.0f nm, #font[12]{l}=%0.0f cm;p, MeV/c;N_{#gamma}", nref,lambda1,lambda2,length));
   gr_p->Draw("pl");
   gr_pi->Draw("pl");
   gr_mu->Draw("pl");
   gr_e->Draw("pl");
   gr_K->Draw("pl");

   TLegend* legend = new TLegend(0.70,0.20, 0.78,0.70);
   legend->SetFillColor(0);
   legend->SetFillStyle(1001);
   legend->AddEntry(gr_e, " #font[12]{e}");
   legend->AddEntry(gr_mu, " #font[12]{#mu}");
   legend->AddEntry(gr_pi, " #font[12]{#pi}");
   legend->AddEntry(gr_K, " #font[12]{K}");
   legend->AddEntry(gr_p, " #font[12]{p}");

   legend->Draw();
}

void Npe_Photek240_momentum_plot(Double_t lambda1=250, Double_t lambda2=650)  // lambda in nm
{
   const Double_t m_p = 938.27;        // MeV
   const Double_t m_e = 0.511;         // MeV
   const Double_t m_pi = 139.6;        // MeV
   const Double_t m_mu = 105.7;        // MeV
   const Double_t m_K = 493.7;         // MeV

   Int_t np_p = 0;
   TGraph* gr_p = new TGraph(0);
   gr_p->SetNameTitle("gr_p", Form("Cherenkov radiation of proton in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's momentum p, MeV/c;Npe", lambda1,lambda2));
   gr_p->SetFillColor(0);
   gr_p->SetMarkerStyle(1);
   gr_p->SetMarkerColor(2);
   gr_p->SetLineColor(2);

   Int_t np_e = 0;
   TGraph* gr_e = new TGraph(0);
   gr_e->SetNameTitle("gr_e", Form("Cherenkov radiation of electron in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's momentum p, MeV/c;Npe", lambda1,lambda2));
   gr_e->SetFillColor(0);
   gr_e->SetMarkerStyle(1);
   gr_e->SetMarkerColor(4);
   gr_e->SetLineColor(4);

   Int_t np_mu = 0;
   TGraph* gr_mu = new TGraph(0);
   gr_mu->SetNameTitle("gr_mu", Form("Cherenkov radiation of #mu in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's momentum p, MeV/c;Npe", lambda1,lambda2));
   gr_mu->SetFillColor(0);
   gr_mu->SetMarkerStyle(1);
   gr_mu->SetMarkerColor(6);
   gr_mu->SetLineColor(6);

   Int_t np_pi = 0;
   TGraph* gr_pi = new TGraph(0);
   gr_pi->SetNameTitle("gr_pi", Form("Cherenkov radiation of #pi in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's momentum p, MeV/c;Npe", lambda1,lambda2));
   gr_pi->SetFillColor(0);
   gr_pi->SetMarkerStyle(1);
   gr_pi->SetMarkerColor(1);
   gr_pi->SetLineColor(1);

   Int_t np_K = 0;
   TGraph* gr_K = new TGraph(0);
   gr_K->SetNameTitle("gr_K", Form("Cherenkov radiation of #font[12]{K} in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's momentum p, MeV/c;Npe", lambda1,lambda2));
   gr_K->SetFillColor(0);
   gr_K->SetMarkerStyle(1);
   gr_K->SetMarkerColor(8);
   gr_K->SetLineColor(8);

   Double_t p1 = 10.;   // MeV/c
   Double_t p2 = 5000.; // MeV/c
   Double_t p = p1;
   Double_t dp = 10.;

   while (p < p2)
   {
      Double_t npe_p = Npe_Photek240(p, m_p, lambda1, lambda2);
      gr_p->SetPoint(np_p, p, npe_p);
      gr_p->Set(++np_p);

      Double_t npe_mu = Npe_Photek240(p, m_mu, lambda1, lambda2);
      gr_mu->SetPoint(np_mu, p, npe_mu);
      gr_mu->Set(++np_mu);

      Double_t npe_pi = Npe_Photek240(p, m_pi, lambda1, lambda2);
      gr_pi->SetPoint(np_pi, p, npe_pi);
      gr_pi->Set(++np_pi);

      Double_t npe_e = Npe_Photek240(p, m_e, lambda1, lambda2);
      gr_e->SetPoint(np_e, p, npe_e);
      gr_e->Set(++np_e);

      Double_t npe_K = Npe_Photek240(p, m_K, lambda1, lambda2);
      gr_K->SetPoint(np_K, p, npe_K);
      gr_K->Set(++np_K);

      p += dp;
   }

   TCanvas* can = new TCanvas;
   can->DrawFrame(0,0.001, 1.10*p2, 1.10*gr_e->GetY()[np_e-1], Form("Cherenkov radiation in Photek240, #lambda=%0.0f-%0.0f nm;proton's momentum p, MeV/c;Npe",lambda1,lambda2));
   gr_p->Draw("pl");
   gr_pi->Draw("pl");
   gr_mu->Draw("pl");
   gr_e->Draw("pl");
   gr_K->Draw("pl");

   TLegend* legend = new TLegend(0.70,0.20, 0.78,0.70);
   legend->SetFillColor(0);
   legend->SetFillStyle(1001);
   legend->AddEntry(gr_e, " #font[12]{e}");
   legend->AddEntry(gr_mu, " #font[12]{#mu}");
   legend->AddEntry(gr_pi, " #font[12]{#pi}");
   legend->AddEntry(gr_K, " #font[12]{K}");
   legend->AddEntry(gr_p, " #font[12]{p}");

   legend->Draw();
}

void Npe_Photek240_T_plot(Double_t lambda1=250, Double_t lambda2=650)  // lambda in nm
{
   // Npe vs kinetic energy T
   //
   // Find proton momentum from the kinetic energy T
   // calculate Npe from other particles with this momentum

   const Double_t m_p = 938.27;        // MeV
   const Double_t m_e = 0.511;         // MeV
   const Double_t m_pi = 139.6;        // MeV
   const Double_t m_mu = 105.7;        // MeV
   const Double_t m_K = 493.7;         // MeV

   Int_t np_p = 0;
   TGraph* gr_p = new TGraph(0);
   gr_p->SetNameTitle("gr_p", Form("Cherenkov radiation of proton in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's kinetic energy T, MeV;Npe", lambda1,lambda2));
   gr_p->SetFillColor(0);
   gr_p->SetMarkerStyle(1);
   gr_p->SetMarkerColor(2);
   gr_p->SetLineColor(2);

   Int_t np_e = 0;
   TGraph* gr_e = new TGraph(0);
   gr_e->SetNameTitle("gr_e", Form("Cherenkov radiation of electron in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's kinetic energy T, MeV;Npe", lambda1,lambda2));
   gr_e->SetFillColor(0);
   gr_e->SetMarkerStyle(1);
   gr_e->SetMarkerColor(4);
   gr_e->SetLineColor(4);

   Int_t np_mu = 0;
   TGraph* gr_mu = new TGraph(0);
   gr_mu->SetNameTitle("gr_mu", Form("Cherenkov radiation of #mu in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's kinetic energy T, MeV;Npe", lambda1,lambda2));
   gr_mu->SetFillColor(0);
   gr_mu->SetMarkerStyle(1);
   gr_mu->SetMarkerColor(6);
   gr_mu->SetLineColor(6);

   Int_t np_pi = 0;
   TGraph* gr_pi = new TGraph(0);
   gr_pi->SetNameTitle("gr_pi", Form("Cherenkov radiation of #pi in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's kinetic energy T, MeV;Npe", lambda1,lambda2));
   gr_pi->SetFillColor(0);
   gr_pi->SetMarkerStyle(1);
   gr_pi->SetMarkerColor(1);
   gr_pi->SetLineColor(1);

   Int_t np_K = 0;
   TGraph* gr_K = new TGraph(0);
   gr_K->SetNameTitle("gr_K", Form("Cherenkov radiation of #font[12]{K} in Photek240. #lambda_{1}=%0.0f nm, #lambda_{2}=%0.0f nm;proton's kinetic energy T, MeV;Npe", lambda1,lambda2));
   gr_K->SetFillColor(0);
   gr_K->SetMarkerStyle(1);
   gr_K->SetMarkerColor(8);
   gr_K->SetLineColor(8);

   Double_t T1 = 10.;   // MeV/c
   Double_t T2 = 1200.; // MeV/c
   Double_t T = T1;
   Double_t dT = 1.;

   while (T < T2)
   {
      Double_t p = Pproton(T);   // momentum of a proton with the kinetic energy T

      Double_t npe_p = Npe_Photek240(p, m_p, lambda1, lambda2);
      gr_p->SetPoint(np_p, T, npe_p);
      gr_p->Set(++np_p);

      Double_t npe_mu = Npe_Photek240(p, m_mu, lambda1, lambda2);
      gr_mu->SetPoint(np_mu, T, npe_mu);
      gr_mu->Set(++np_mu);

      Double_t npe_pi = Npe_Photek240(p, m_pi, lambda1, lambda2);
      gr_pi->SetPoint(np_pi, T, npe_pi);
      gr_pi->Set(++np_pi);

      Double_t npe_e = Npe_Photek240(p, m_e, lambda1, lambda2);
      gr_e->SetPoint(np_e, T, npe_e);
      gr_e->Set(++np_e);

      Double_t npe_K = Npe_Photek240(p, m_K, lambda1, lambda2);
      gr_K->SetPoint(np_K, T, npe_K);
      gr_K->Set(++np_K);

      T += dT;
   }

   TCanvas* can = new TCanvas;
   can->DrawFrame(0,0.001, 1.10*T2, 1.10*gr_e->GetY()[np_e-1], Form("Cherenkov radiation in Photek240, #lambda=%0.0f-%0.0f nm;proton's kinetic energy T, MeV;Npe",lambda1,lambda2));
   gr_p->Draw("pl");
   gr_pi->Draw("pl");
   gr_mu->Draw("pl");
   gr_e->Draw("pl");
   gr_K->Draw("pl");

   TLegend* legend = new TLegend(0.80,0.14, 0.88,0.64);
   legend->SetFillColor(0);
   legend->SetFillStyle(1001);
   legend->AddEntry(gr_e, " #font[12]{e}");
   legend->AddEntry(gr_mu, " #font[12]{#mu}");
   legend->AddEntry(gr_pi, " #font[12]{#pi}");
   legend->AddEntry(gr_K, " #font[12]{K}");
   legend->AddEntry(gr_p, " #font[12]{p}");

   legend->Draw();
}
