#include <TROOT.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <iostream>
#include <cassert>
#include <vector>

using std::cout;     using std::endl;

class GaussNewtonMinimizer {
private:
   Int_t npar_;
   TVectorD par_;       // vector of parameters, dimension npar_
   TVectorD Jrsum_;     // vector with dimension npar_ to keep product of Jacobian matrix with vector of residuals
   TMatrixD JJsum_;     // matrix (npar_,npar) to keep square of the Jacobian matrix: J(transposed)*J
   Double_t r2sum_;
public:
   GaussNewtonMinimizer(Int_t npar, Double_t* vpar0=0): npar_(npar), r2sum_(0)
   {
      par_.ResizeTo(npar_);
      Jrsum_.ResizeTo(npar_);
      JJsum_.ResizeTo(npar_,npar_);
      if (vpar0) for (int ipar=0; ipar<npar_; ++ipar) par_(ipar) = vpar0[ipar];  // initial values
      InitIteration();
   }
   void InitIteration() {
      for (int irow=0; irow<npar_; ++irow) {
         Jrsum_(irow) = 0;
         for (int icol=0; icol<npar_; ++icol) {
            JJsum_(irow,icol) = 0;
         }
      }
      r2sum_ = 0;
   }
   virtual ~GaussNewtonMinimizer() {}
   Double_t OldSumResiduals() const {return r2sum_;}
   Double_t GetParameter(Int_t ipar) const {return par_(ipar);}
   const TVectorD& GetParameters() const {return par_;}
   void SetParameter(Int_t ipar, Double_t val) {
      assert(ipar < npar_);
      par_(ipar) = val;
   }
   void SetParameters(const Double_t* vpar0) {
      for (int ipar=0; ipar<npar_; ++ipar) par_(ipar) = vpar0[ipar];
   }
   void Feed(Double_t x, Double_t y, Double_t* factors=0) {
      Double_t jacobian[100];
      for (int ipar=0; ipar<npar_; ++ipar) jacobian[ipar] = Jacobian(x,ipar,factors);
      Double_t r = Residual(x,y,factors);
      r2sum_ += r*r;
      for (int ipar=0; ipar<npar_; ++ipar) {
         // Jrsum_(ipar) += jacobian[ipar]*Residual(x,y,factors);
         Jrsum_(ipar) += jacobian[ipar]*r;
      }
      for (int ipar=0; ipar<npar_; ++ipar) {
         for (int jpar=0; jpar<npar_; ++jpar) {
            JJsum_(ipar,jpar) += jacobian[ipar]*jacobian[jpar];
         }
      }
   }
   virtual Double_t Residual(Double_t x, Double_t y, Double_t* factors=0) const = 0;
   // virtual Double_t Residual(Double_t x, Double_t y, Double_t* factors=0) {
   //    if (!factors) return 0;
   //    Double_t r = y - (par_(0) + par_(1)*x + par_(2)*pow(x,3));
   //    return r;
   // }
   virtual Double_t Jacobian(Double_t x, Int_t ipar, Double_t* factors=0) const = 0;
   // virtual Double_t Jacobian(Double_t x, Int_t ipar, Double_t* factors=0) {
   //    if (!factors) return 0;
   //    Double_t drdp = 0;
   //    switch (ipar) {
   //       case 0:  drdp = -1;        break;
   //       case 1:  drdp = -x;        break;
   //       case 2:  drdp = -pow(x,3); break;
   //       default: assert(ipar >= 0 && ipar < npar_);
   //    }
   //    return drdp;
   // }
   void Correct()
   {
      // cout<< "before inversion: JJsum_.Print()";
      // JJsum_.Print();
      Double_t det;
      JJsum_.Invert(&det);
      // cout<< "after inversion: JJsum_.Print(), det = " << det;
      // JJsum_.Print();

      // cout<< "mult JJ by Jr" <<endl;
      TVectorD corr = JJsum_*Jrsum_;

      cout<< "JJsum_(0,0) = " << JJsum_(0,0) << " Jrsum_(0) = " << Jrsum_(0) << "   Correction JJsum_*Jrsum_ to be subtracted:";    corr.Print();

      par_ -= corr;

      //cout<< "par_(0) = " << par_(0) << " par_(1) = " << par_(1) << " par_(2) = " << par_(2) <<endl;
   }
};

class GaussNewtonMinimizerGaussian: public GaussNewtonMinimizer {
public:
   GaussNewtonMinimizerGaussian(Int_t npar, Double_t* vpar0=0): GaussNewtonMinimizer(npar,vpar0) {}
   virtual Double_t Residual(Double_t x, Double_t y, Double_t* factors=0) const {
      if (!factors) return 0;
      Double_t arg = -0.5*(x - par_(1))*(x - par_(1)) / par_(2);
      Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;
      Double_t r = par_(0)*exp - y;
      return r;
   }
   virtual Double_t Jacobian(Double_t x, Int_t ipar, Double_t* factors=0) const {
      if (!factors) return 0;
      Double_t drdp = 0;
      Double_t arg = -0.5*(x - par_(1))*(x - par_(1)) / par_(2);
      Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;
      switch (ipar) {
         case 0:  drdp = exp;
                  break;
         case 1:  drdp = (x - par_(1))/(par_(2)*par_(2)) * par_(0)*exp;
                  break;
         case 2:  drdp = (x-par_(1))*(x-par_(1)) / pow(par_(2),3) * par_(0)*exp;
                  break;
         default: assert(ipar >= 0 && ipar < npar_);
      }
      return drdp;
   }
};

void minimizeGaussian()
{
   GaussNewtonMinimizerGaussian gnmin(3);

   const Double_t xmin = -5;
   const Double_t xmax = 5;

   TH1F* hgaus = new TH1F("hgaus", "hgaus", 100, xmin, xmax);
   hgaus->FillRandom("gaus", 1000);

   new TCanvas;
   hgaus->Draw();

   //
   // fit the parameters
   //

   Double_t par0[3];
   par0[0] = hgaus->GetMaximum();
   par0[1] = hgaus->GetMean();
   par0[2] = hgaus->GetRMS();

   gnmin.SetParameters(par0);
   const TVectorD& par = gnmin.GetParameters();

   Double_t res2sum = 0;
   for (int i=1; i<=hgaus->GetNbinsX(); ++i) {
      Double_t r = gnmin.Residual(hgaus->GetBinCenter(i), hgaus->GetBinContent(i));
      res2sum += r*r;
   }
   cout<< "init: par(0) = " << par(0) << " par(1) = " << par(1) << " par(2) = " << par(2) << " res2sum = " << res2sum <<endl;

   for (int iter=0; iter<3; ++iter)
   {
      cout<< "\niteration #" << iter <<endl;
      gnmin.InitIteration();
      for (int i=1; i<=hgaus->GetNbinsX(); ++i) {
         gnmin.Feed(hgaus->GetBinCenter(i), hgaus->GetBinContent(i));
      }
      gnmin.Correct();

      res2sum = 0;
      for (int i=1; i<=hgaus->GetNbinsX(); ++i) {
         Double_t r = gnmin.Residual(hgaus->GetBinCenter(i), hgaus->GetBinContent(i));
         res2sum += r*r;
      }
      cout<< "par(0) = " << par(0) << " par(1) = " << par(1) << " par(2) = " << par(2) << " res2sum = " << res2sum << " gnmin.OldSumResiduals() = " << gnmin.OldSumResiduals() <<endl;
   }

   TF1* fgaus = new TF1("fgaus","gaus",xmin,xmax);
   fgaus->SetParameters(par(0),par(1),par(2));
   //fgaus->SetLineColor(4);
   fgaus->Draw("same");
}

class Hit: public TObject {
   friend std::ostream& operator<<(std::ostream&, const Hit&);
public:
   Float_t x;
   Float_t y;
   Float_t z;
public:
   void clear() {
      x = y = z = 0;
   }
   Hit(): TObject() {
      clear();
   }
   Hit(Double_t x_, Double_t y_, Double_t z_): TObject(), x(x_), y(y_), z(z_) {}
   Float_t dHit(const Hit& hit) const {
      Float_t dx = x - hit.x;
      Float_t dy = y - hit.y;
      Float_t dz = z - hit.z;
      return TMath::Sqrt(dx*dx + dy*dy + dz*dz);
   }
   ClassDef(Hit, 1)
};

std::ostream& operator << (std::ostream& os, const Hit& hit) {
   return os << "hit.x " << hit.x << " hit.y " << hit.y << " hit.z " << hit.z;
}

class Ray: public TObject {
   friend std::ostream& operator<<(std::ostream&, const Ray&);
public:
   Float_t x, y, z;	        // radiant
   Float_t cx, cy, cz;	    // direction cosines
public:
   void clear() {
      x = y = z = 0;
      cx = cy = cz = 0;
   }
   Ray(): TObject() {
      clear();
   }
   Ray(const Ray& ray): TObject() {
      x = ray.x;
      y = ray.y;
      z = ray.z;
      cx = ray.cx;
      cy = ray.cy;
      cz = ray.cz;
   }
   Ray(Double_t x1, Double_t y1, Double_t z1, Double_t x2, Double_t y2, Double_t z2): TObject() {
      x = x1;
      y = y1;
      z = z1;
      cx = x2 - x1;
      cy = y2 - y1;
      cz = z2 - z1;
      Double_t length = cx*cx + cy*cy + cz*cz;
      length = TMath::Sqrt(length);
      if (length > 0) {
      	 cx /= length;
      	 cy /= length;
      	 cz /= length;
      }
      else cx = cy = cz = 0;
   }
   Ray(Double_t x1, Double_t y1, Double_t z1, Double_t costheta, Double_t phi): TObject() {
      x = x1;
      y = y1;
      z = z1;
      cz = costheta;    // NB: the parameter is costheta
      Double_t sintheta = TMath::Sqrt(1.-costheta*costheta);
      cx = sintheta*TMath::Cos(phi);
      cy = sintheta*TMath::Sin(phi);
      Double_t length = cx*cx + cy*cy + cz*cz;
      length = TMath::Sqrt(length);
      if (length > 0) {
      	 cx /= length;
      	 cy /= length;
      	 cz /= length;
      }
      else cx = cy = cz = 0;
   }
   Ray(const Hit& hit1, const Hit& hit2) {
      x = hit1.x;
      y = hit1.y;
      z = hit1.z;
      cx = hit2.x - hit1.x;
      cy = hit2.y - hit1.y;
      cz = hit2.z - hit1.z;
      Double_t length = cx*cx + cy*cy + cz*cz;
      length = TMath::Sqrt(length);
      if (length > 0) {
      	 cx /= length;
      	 cy /= length;
      	 cz /= length;
      }
      else cx = cy = cz = 0;
   }
   void Normalize() {
      Double_t length = Length();
      cx /= length;
      cy /= length;
      cz /= length;
   }
   // void NormalizeXY() {
   //    Double_t length = Length();
   //    cx /= length;
   //    cy /= length;
   //    cz /= length;
   // }
   Float_t static Angle(const Ray& ray1, const Ray& ray2) {
      return TMath::ACos(ray1.cx*ray2.cx + ray1.cy*ray2.cy + ray1.cz*ray2.cz);
   }
   Float_t Angle(const Ray& ray) const {
      return TMath::ACos(cx*ray.cx + cy*ray.cy + cz*ray.cz);
   }
   Float_t Length() const {return TMath::Sqrt(cx*cx + cy*cy + cz*cz);}    // should be 1. Good to check.
   // Float_t Length() const {
   //    Double_t sum = 0;
   //    Double_t min = cx;
   //    if (TMath::Abs(cy) < TMath::Abs(cx)) {
   //       min = cy;
   //       sum += cx*cx;
   //    }
   //      TODO
   //    return TMath::Sqrt(sum);
   // }
   Float_t Theta() const {return TMath::ACos(cz);}
   Float_t xat(Float_t zplane) const {
      const Float_t eps = 1e-7;
      Float_t xproj = 1./eps;
      if (TMath::Abs(cz) > eps) {
      	 Double_t t = (zplane - z)/cz;
      	 xproj = x + cx*t;
      }
      return xproj;
   }
   Float_t yat(Float_t zplane) const {
      const Float_t eps = 1e-7;
      Float_t yproj = 1./eps;
      if (TMath::Abs(cz) > eps) {
      	 Double_t t = (zplane - z)/cz;
      	 yproj = y + cy*t;
      }
      return yproj;
   }
   Float_t dHit(const Hit& hit) const {
      // distance to the hit
      Ray ray(x,y,z, hit.x,hit.y,hit.z);
      Float_t dx = hit.x - x;
      Float_t dy = hit.y - y;
      Float_t dz = hit.z - z;
      Float_t cos = cx*ray.cx + cy*ray.cy + cz*ray.cz;
      Float_t sin = 1. - cos*cos;
      sin = sin >= 0? TMath::Sqrt(sin): 0;
      return sin * TMath::Sqrt(dx*dx + dy*dy + dz*dz);
   }
   //void Print() {
   //}

   ClassDef(Ray, 3)
};

std::ostream& operator << (std::ostream& os, const Ray& ray) {
   os << "x = " << ray.x << " y = " << ray.y << " z = " << ray.z << " cx = " << ray.cx << " cy = " << ray.cy << " cz = " << ray.cz;
   return os;
}

class Sensor {
public:
   // q = R*(r - r0) - dq, where r is global
   TMatrixD R_;       // rotation matrix (unit matrix initially)
   TVectorD r0_;      // system origin in global system
   TVectorD dq_;      // local coordinates correction
public:
   const TVectorD& GetPosition() const {return r0_;}
   Sensor() {
      R_.ResizeTo(3,3);
      r0_.ResizeTo(3);
      dq_.ResizeTo(3);

      for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) R_(i,j) = 0;
      R_(0,0) = R_(1,1) = R_(2,2) = 1;

      dq_(0) = dq_(1) = dq_(2) = 0;

      r0_(0) = 0;
      r0_(1) = 0;
      r0_(2) = 0;
   }
   Sensor(Double_t z) {
      R_.ResizeTo(3,3);
      r0_.ResizeTo(3);
      dq_.ResizeTo(3);

      for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) R_(i,j) = 0;
      R_(0,0) = R_(1,1) = R_(2,2) = 1;

      dq_(0) = dq_(1) = dq_(2) = 0;

      r0_(0) = 0;
      r0_(1) = 0;
      r0_(2) = z;
   }
   void SetPosition(Double_t x, Double_t y, Double_t z) {
      // position of the sensor origin in global system
      r0_(0) = x;
      r0_(1) = y;
      r0_(2) = z;
   }
   void Print() const {
      cout<< "r0_: " << r0_(0) << " " << r0_(1) << " " << r0_(2) << " dq_: " << dq_(0) <<" "<< dq_(1) <<" "<< dq_(2) <<endl;
      cout<< "R_"; R_.Print();
   }
   void SetRotation(Double_t alpha, Double_t beta, Double_t gamma) {
      // rotations around these axis in this order:
      // around x -- alpha, around new y -- beta, around new z -- gamma
      R_(0,0) = cos(beta)*cos(gamma);
      R_(0,1) = -sin(alpha)*sin(beta)*cos(gamma) + cos(alpha)*sin(gamma);
      R_(0,2) = cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma);
      R_(1,0) = -cos(beta)*sin(gamma);
      R_(1,1) = sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma);
      R_(1,2) = -cos(alpha)*sin(beta)*sin(gamma) + sin(alpha)*cos(gamma);
      R_(2,0) = -sin(beta);
      R_(2,1) = -sin(alpha)*cos(beta);
      R_(2,2) = cos(alpha)*cos(beta);
      cout<< "--- Sensor::SetRotation: R_.Print()"; R_.Print();
   }
   void AddAlpha(Double_t alpha) {
      // rotation around the x-axis
      TMatrixD dR(3,3);
      for (int i=0; i<3; ++i) for (int j=0; j<2; ++j) dR(i,j) = 0;
      dR(1,1) = cos(alpha);
      dR(1,2) = sin(alpha);
      dR(2,1) = -sin(alpha);
      dR(2,2) = cos(alpha);
      dR(0,0) = 1.;
      R_ *= dR;
   }
   void AddBeta(Double_t beta) {
      // rotation around the y-axis
      TMatrixD dR(3,3);
      for (int i=0; i<3; ++i) for (int j=0; j<2; ++j) dR(i,j) = 0;
      dR(0,0) = cos(beta);
      dR(0,2) = sin(beta);
      dR(2,0) = -sin(beta);
      dR(2,2) = cos(beta);
      dR(1,1) = 1.;
      R_ *= dR;
   }
   const TMatrixD& GetRotation() const {return R_;}
   TVectorD ToLocal(const TVectorD& r) const {
      TVectorD q = R_*(r - r0_) - dq_;
      return q;
   }
   Ray ToLocal(const Ray& ray) const {
      // cout<< "ToLocal: r0_: " << r0_(0) << " " << r0_(1) << " " << r0_(2) << " dq_: " << dq_(0) <<" "<< dq_(1) <<" "<< dq_(2) <<endl;
      Ray oray;
      TVectorD r(3);
      r(0) = ray.x;
      r(1) = ray.y;
      r(2) = ray.z;
      // cout<< "ToLocal: R_.Print():"; R_.Print();
      TVectorD q = R_*(r - r0_) - dq_;
      oray.x = q(0);
      oray.y = q(1);
      oray.z = q(2);
      // direction cosines
      TVectorD r1(3);            // point along the r on distance 1 from the radiant
      r1(0) = ray.x + ray.cx;
      r1(1) = ray.y + ray.cy;
      r1(2) = ray.z + ray.cz;
      q = R_*(r1 - r0_) - dq_;   // end point in the local system
      // cout<< "ToLocal: oray.x = " << oray.x << " oray.y = " << oray.y << " oray.z = " << oray.z << " q(0) = " << q(0) << " q(1) = " << q(1) << " q(2) = " << q(2) <<endl;
      oray.cx = q(0) - oray.x;
      oray.cy = q(1) - oray.y;
      oray.cz = q(2) - oray.z;
      // cout<< "ToLocal: oray.cx = " << oray.cx << " oray.cy = " << oray.cy << " oray.cz = " << oray.cz <<endl;
      oray.Normalize();
      // cout<< "ray.Length() = " << ray.Length() << " oray.Length() = " << oray.Length() <<endl;
      return oray;
   }
   TVectorD ToGlobal(const TVectorD& q) const {
      TMatrixD Rinv(TMatrixD::kInverted, R_);
      TVectorD r = r0_ + Rinv*(q + dq_);
      return r;
   }
   Ray ToGlobal(const Ray& ray) const {
      TMatrixD Rinv(TMatrixD::kInverted, R_);
      Ray oray;
      TVectorD q(3);
      q(0) = ray.x;
      q(1) = ray.y;
      q(2) = ray.z;
      TVectorD r = r0_ + Rinv*(q + dq_);
      oray.x = r(0);
      oray.y = r(1);
      oray.z = r(2);
      // direction cosines
      TVectorD dq(3);
      dq(0) = ray.x + ray.cx;
      dq(1) = ray.y + ray.cy;
      dq(2) = ray.z + ray.cz;
      r = r0_ + Rinv*(dq + dq_);
      oray.cx = r(0) - oray.x;
      oray.cy = r(1) - oray.y;
      oray.cz = r(2) - oray.z;
      oray.Normalize();
      // cout<< "ray.Length() = " << ray.Length() << " oray.Length() = " << oray.Length() <<endl;
      return oray;
   }
   TVectorD CrossLoc(const Ray& rayloc) const {
      // intersection with XY plane
      // input: ray in the local system
      // output: cross point in the local system
      const Double_t eps = 1e-15;
      TVectorD cross(3);
      cross(0) = cross(1) = cross(2) = 0;
      if (TMath::Abs(rayloc.cz) > eps) {
         Double_t t = -rayloc.z/rayloc.cz;
         //cout<< "CrossLoc: t = " << t << " rayloc.cz*t = " << rayloc.cz*t <<endl;
         //printf("CrossLoc: t = %0.16lg rayloc.cz*t = %0.16lg\n", t, rayloc.cz*t);
         if (t > 0) {                           // correct direction
            cross(0) = rayloc.x + rayloc.cx*t;
            cross(1) = rayloc.y + rayloc.cy*t;
            cross(2) = rayloc.z + rayloc.cz*t;
         }
      }
      // cout<< "CrossLoc: in local: cross(0) = " << cross(0) << " cross(1) = " << cross(1) << " cross(2) = " << cross(2) <<endl;
      return cross;
   }
   TVectorD CrossGlo(const Ray& rayglo) const {
      // intersection with XY plane
      // input: ray in the global system
      // output: cross point in the global system
      const Double_t eps = 1e-15;
      Ray rayloc = ToLocal(rayglo);
      //cout<< "CrossGlo: rayloc.x = " << rayloc.x << " rayloc.y = " << rayloc.y << " rayloc.z = " << rayloc.z << "   rayloc.cx = " << rayloc.cx << " rayloc.cy = " << rayloc.cy << " rayloc.cz = " << rayloc.cz <<endl;
      TVectorD cross(3);
      cross(0) = cross(1) = cross(2) = 0;
      if (TMath::Abs(rayloc.cz) > eps) {
         Double_t t = -rayloc.z/rayloc.cz;
         if (t > 0) {                           // correct direction
            cross(0) = rayloc.x + rayloc.cx*t;
            cross(1) = rayloc.y + rayloc.cy*t;
            cross(2) = rayloc.z + rayloc.cz*t;
         }
      }
      //cout<< "CrossGlo: in local: cross(0) = " << cross(0) << " cross(1) = " << cross(1) << " cross(2) = " << cross(2) <<endl;
      return ToGlobal(cross);
   }
   TVectorD LocalCross(const Ray& rayglo) const {
      // intersection with XY plane
      // input: ray in the global system
      // output: cross point in the local system
      const Double_t eps = 1e-15;
      Ray rayloc = ToLocal(rayglo);
      //cout<< "LocalCross: rayloc.x = " << rayloc.x << " rayloc.y = " << rayloc.y << " rayloc.z = " << rayloc.z << "   rayloc.cx = " << rayloc.cx << " rayloc.cy = " << rayloc.cy << " rayloc.cz = " << rayloc.cz <<endl;
      TVectorD cross(3);
      cross(0) = cross(1) = cross(2) = 0;
      if (TMath::Abs(rayloc.cz) > eps) {
         Double_t t = -rayloc.z/rayloc.cz;
         if (t > 0) {                           // correct direction
            cross(0) = rayloc.x + rayloc.cx*t;
            cross(1) = rayloc.y + rayloc.cy*t;
            cross(2) = rayloc.z + rayloc.cz*t;
         }
      }
      return cross;
   }
   // TVectorD Cross(const Ray& rayglo) const {
   //    // intersection with XY plane
   //    Ray rayloc = ToLocal(rayglo);
   //    // cout<< "Cross: rayloc.x = " << rayloc.x << " rayloc.y = " << rayloc.y << " rayloc.z = " << rayloc.z << "   rayloc.cx = " << rayloc.cx << " rayloc.cy = " << rayloc.cy << " rayloc.cz = " << rayloc.cz <<endl;
   //    TVectorD crossloc(3);
   //    crossloc = CrossLoc(rayloc);
   //    return ToGlobal(crossloc);
   // }
};

void align_test()
{
   const Double_t F = 4000;                  // distance to the foil, mm
   const Double_t pitch = 0.228;             // mm
   const Double_t ssd_side_length = 89.5;    // mm
   //-- const Double_t plane_separation = 98.0;   // mm

   const Double_t costheta_span = 5e-5;

   Int_t xstrip[4][384];

   Sensor sensor[4];
   sensor[0].SetPosition(0,0,-200);
   sensor[1].SetPosition(0,0,-100);
   sensor[2].SetPosition(0,0,100);
   sensor[3].SetPosition(0,0,200);
   // sensor[0].SetRotation(0,0,0);
   // sensor[1].SetRotation(0,0,0);
   // sensor[2].SetRotation(0,0,0);
   // sensor[3].SetRotation(0,0,0);

   Double_t xpos[4][1000];
   Double_t ypos[4][1000];

   TRandom3 rand;

   // Int_t nevt = 1000000;
   Int_t nevt = 1000;

   TH2F* h2hits = new TH2F("h2hits","hits in the first plane", 100,-ssd_side_length/2,ssd_side_length/2, 100,-ssd_side_length/2,ssd_side_length/2);
   h2hits->SetMarkerColor(2);

   TH1F* h0 = new TH1F("h0","fired strips in the first plane", 384,0,384);

   Double_t beamspot_x = 0;
   Double_t beamspot_y = 0;
   Double_t beamspot_z = -F;

   for (int ievt=0; ievt<nevt; ++ievt)
   {
      // generate points uniformly distributed on the sphere

      // generate cos(theta) in range (1-5e-5, 1)
      Double_t delta = costheta_span*rand.Rndm(ievt);   // 1 - costheta
      Double_t costheta = 1 - delta;

      // generate phi
      Double_t phi = 2.*TMath::Pi()*rand.Rndm(nevt+ievt);

      Ray ray(beamspot_x,beamspot_y,beamspot_z, costheta,phi);

      //-- Double_t sintheta = TMath::Sqrt(2.*delta - delta*delta);  // for small delta
      Double_t sintheta = TMath::Sqrt(1.-costheta*costheta);
      Double_t cx = sintheta*TMath::Cos(phi);
      Double_t cy = sintheta*TMath::Sin(phi);
      // Double_t cz = costheta;    // NB cx*cx + cy*cy + cz*cz = 1

      // ray length
      Double_t t = (F - sensor[0].GetPosition()(2)) / costheta;
      Double_t x = cx*t;
      Double_t y = cy*t;

      xpos[0][ievt] = x;
      ypos[0][ievt] = y;

      // cout<< "cx " << cx << " cy " << cy << " t " << t << " x " << x << " y " << y <<endl;

      h2hits->Fill(x,y);

      Int_t istrip = Int_t((x + 0.5*ssd_side_length) / pitch);
      if (istrip < 384) ++xstrip[0][istrip];
      h0->Fill(istrip);
   }

   new TCanvas;
   h2hits->Draw();

   new TCanvas;
   h0->Draw();

   // play with last event

   cout<< "test ToLocal and ToGlobal back" <<endl;

   TVectorD vglobal0(3);
   TVectorD vglobal(3);
   TVectorD vlocal(3);
   vglobal0(0) = xpos[0][nevt-1];
   vglobal0(1) = ypos[0][nevt-1];
   vglobal0(2) = sensor[0].GetPosition()(2);
   vlocal = sensor[0].ToLocal(vglobal0);
   vglobal = sensor[0].ToGlobal(vlocal);

   cout<< "vglobal0";   vglobal0.Print();
   cout<< "vglobal ";   vglobal.Print();

   //-- rotate sensor[0]

   const Double_t pi = TMath::Pi();
   // sensor[0].SetRotation(pi/12, pi/9, pi/6);
   // // sensor[0].SetRotation(1, 1, 0.01);
   // cout<< "sensor[0].GetRotation().Print()"; sensor[0].GetRotation().Print();
   //-- sensor[0].SetRotation(pi/6, 0, 0);
   sensor[0].SetRotation(0, 0, pi/6);        // ok with theta=0, phi=0

   TMatrixD RT(TMatrixD::kTransposed, sensor[0].GetRotation());
   TMatrixD unity = RT*sensor[0].GetRotation();
   cout<< "rotation matrix is unitary"; unity.Print();


   cout<< "**********************************\nTest Ray" <<endl;
   //sensor[0].SetRotation(0, 0, 0);        // ok with theta=0, phi=0
   //sensor[0].SetRotation(0, 0, pi/6);
   sensor[0].SetRotation(pi/6, 0, 0);
   cout<< "sensor[0] is: " <<endl;
   sensor[0].Print();

   // generate cos(theta) in range (1-5e-5, 1)
   Double_t delta = costheta_span*rand.Rndm(0);   // 1 - costheta
   Double_t costheta = 1 - delta;
   costheta = 1;                             //-- set explicit zero
   costheta = TMath::Cos(TMath::ATan(19./3800.));

   // generate phi
   Double_t phi = 2.*TMath::Pi()*rand.Rndm(1);
   phi = 0;                                  //-- set explicit zero

   Ray ray(beamspot_x,beamspot_y,beamspot_z, costheta,phi);

   // ray
   cout<< "-- call ToLocal" <<endl;
   Ray rloc = sensor[0].ToLocal(ray);
   cout<< "-- call ToGlobal" <<endl;
   Ray rglo = sensor[0].ToGlobal(rloc);
   cout<< "ray:  " << ray <<endl;
   cout<< "rloc: " << rloc <<endl;
   cout<< "rglo: " << rglo <<endl;

   TVectorD crossLoc(3);
   crossLoc = sensor[0].CrossLoc(rloc);
   cout<< "Local intersection point: crossLoc(0) = " << crossLoc(0) << " crossLoc(1) = " << crossLoc(1) << " crossLoc(2) = " << crossLoc(2) <<endl;

   TVectorD cross(3);
   cross = sensor[0].CrossGlo(ray);
   cout<< "Global intersection point: cross(0) = " << cross(0) << " cross(1) = " << cross(1) << " cross(2) = " << cross(2) <<endl;
}

class GaussNewtonMinimizerAlign: public GaussNewtonMinimizer {
public:
   GaussNewtonMinimizerAlign(Int_t npar, Double_t* vpar0=0): GaussNewtonMinimizer(npar,vpar0) {}
   virtual Double_t Residual(Double_t x, Double_t y, Double_t* factors=0) const {
      //-- Double_t r = x - par_(0) - y;
      //Double_t r = x - par_(0)*factors[0] + x*par_[1]*factors[1]  - y;
      Double_t r = x - par_(0)*factors[0] - y;
      if (npar_ > 1) r += x*par_(1)*factors[1];
      return r;
   }
   virtual Double_t Jacobian(Double_t x, Int_t ipar, Double_t* factors=0) const {
      //x += 0;                 // to avoid a message about the unused parameter
      Double_t drdp = 0;
      switch (ipar) {
         case 0:  drdp = -1*factors[0];
                  break;
         case 1:  if (npar_ > 1) drdp = x*factors[1];
                  break;
         default: assert(ipar >= 0 && ipar < npar_);
      }
      return drdp;
   }
};

void align()
{
   const Double_t F = 4000;                  // distance to the foil, mm
   const Double_t pitch = 0.228;             // mm
   const Double_t ssd_side_length = 89.5;    // mm
   //-- const Double_t plane_separation = 98.0;   // mm

   //-- const Double_t costheta_span = 5e-5;
   const Double_t costheta_span = 5e-5;

   Int_t xstrip[4][384];

   Sensor sensorReal[4];
   //-- sensorReal[0].SetPosition(0,0,-200); // perfect alignment
   //sensorReal[0].SetPosition(0.1,   0.2,  -200);
   //sensorReal[1].SetPosition(0.1,   0.1,  -100);
   //sensorReal[2].SetPosition(-0.1,  -0.1, 100);
   //sensorReal[3].SetPosition(0.1,   -0.1, 200);

   // sensorReal[0].SetPosition(0.1,   0.2,  -200);
   sensorReal[0].SetPosition(0.0,   0.0,  -200);
   sensorReal[1].SetPosition(0.0,   0.0,  -100);
   sensorReal[2].SetPosition(-0.0,  -0.0,  100);
   sensorReal[3].SetPosition(0.0,   -0.0,  200);
   // rotation
   // sensorReal[0].SetRotation(0.01, 0.02, 0);
   sensorReal[0].SetRotation(0.00, 0.00, 0);
   // sensorReal[1].SetRotation(0,0,0);
   // sensorReal[2].SetRotation(0,0,0);
   // sensorReal[3].SetRotation(0,0,0);

   Sensor sensor[4];
   sensor[0].SetPosition(0,0,-200);
   sensor[1].SetPosition(0,0,-100);
   sensor[2].SetPosition(0,0, 100);
   sensor[3].SetPosition(0,0, 200);

   Double_t xpos[4][100000];
   Double_t ypos[4][100000];
   Double_t zpos[4];
   for (int isensor=0; isensor<4; ++isensor) zpos[isensor] = sensorReal[isensor].r0_(2);

   TRandom3 rand;

   // Int_t nevt = 100000;
   Int_t nevt = 1000;

   TH2F* h2hits = new TH2F("h2hits","hits in the first plane", 100,-ssd_side_length/2,ssd_side_length/2, 100,-ssd_side_length/2,ssd_side_length/2);
   h2hits->SetMarkerColor(2);

   TH2F* h2foil = new TH2F("h2foil","foil hits", 1000,-ssd_side_length/2,ssd_side_length/2, 1000,-ssd_side_length/2,ssd_side_length/2);
   h2foil->SetMarkerColor(4);

   TH2F* h2sensor[4];
   for (int isensor=0; isensor<4; ++isensor) {
      h2sensor[isensor] = new TH2F(Form("h2sensor%d",isensor),Form("sensor%d hits (in the first iteration)",isensor), 1000,-ssd_side_length/2,ssd_side_length/2, 1000,-ssd_side_length/2,ssd_side_length/2);
      h2sensor[isensor]->SetMarkerColor(8);
   }

   TH1F* h0 = new TH1F("h0","fired strips in the first plane", 384,0,384);

   TH1F* hr[4];
   for (int isensor=0; isensor<4; ++isensor) {
      hr[isensor] = new TH1F(Form("hr%d",isensor), Form("distance from fit in the sensor%d (in the first iteration)",isensor), 1000,-50,50);
   }
   TH1F* hxpos[4];
   for (int isensor=0; isensor<4; ++isensor) {
      hxpos[isensor] = new TH1F(Form("hxpos%d",isensor), Form("position of the real hit in the sensor%d",isensor), 1000,-50,50);
   }
   TH1F* hLocalCross[4];
   for (int isensor=0; isensor<4; ++isensor) {
      hLocalCross[isensor] = new TH1F(Form("hLocalCross%d",isensor), Form("LocalCross in the sensor%d (in the first iteration)",isensor), 1000,-50,50);
   }

   Double_t beamspot_x = 0;
   Double_t beamspot_y = 0;
   Double_t beamspot_z = -F;

   // just to make sure ...
   Ray ray0(beamspot_x, beamspot_y, beamspot_z, 0, 0, 1);
   TVectorD cross_point(3);
   cross_point = sensorReal[0].LocalCross(ray0);
   cout<< "sensorReal: cross_point(0) = " << cross_point(0) <<endl;
   cross_point = sensor[0].LocalCross(ray0);
   cout<< "sensor:     cross_point(0) = " << cross_point(0) <<endl;

   for (int ievt=0; ievt<nevt; ++ievt)
   {
      // generate points uniformly distributed on the sphere

      // generate cos(theta) in range (1-5e-5, 1)
      Double_t delta = costheta_span*rand.Rndm(ievt);   // 1 - costheta
      Double_t costheta = 1 - delta;

      // generate phi
      Double_t phi = 2.*TMath::Pi()*rand.Rndm(nevt+ievt);

      Ray ray(beamspot_x,beamspot_y,beamspot_z, costheta,phi);
      //-- Ray ray(beamspot_x,beamspot_y,beamspot_z, 1,0);

      for (int isensor=0; isensor<4; ++isensor) {
         TVectorD xpoint(3);
         xpoint = sensorReal[isensor].LocalCross(ray);
         xpos[isensor][ievt] = xpoint(0);
         ypos[isensor][ievt] = xpoint(1);
         if (isensor == 0) {
            h2hits->Fill(xpoint(0),xpoint(1));

            Int_t istrip = Int_t((xpoint(0) + 0.5*ssd_side_length) / pitch);
            if (istrip < 384) ++xstrip[0][istrip];
            h0->Fill(istrip);
         }
      }
   }

   //
   // fit the parameters
   //

   GaussNewtonMinimizerAlign gnmin_x(2);
   GaussNewtonMinimizerAlign gnmin_y(2);

   Double_t par0_x[10];
   par0_x[0] = sensor[0].r0_(0);
   par0_x[1] = 0;

   Double_t par0_y[10];
   par0_y[0] = sensor[0].r0_(1);
   par0_y[1] = 0;

   Double_t factors[10];
   factors[0] = 1;
   factors[1] = 0;

   gnmin_x.SetParameters(par0_x);
   const TVectorD& par_x = gnmin_x.GetParameters();    // pointer to parameters
   gnmin_y.SetParameters(par0_y);
   const TVectorD& par_y = gnmin_y.GetParameters();    // pointer to parameters

   Double_t res2sum_x = 0;
   Double_t res2sum_y = 0;
   // for (int ievt=0; ievt<nevt; ++ievt) {
   //    Double_t r = gnmin.Residual(hgaus->GetBinCenter(i), hgaus->GetBinContent(i));
   //    res2sum += r*r;
   // }

   TGraph* gx = new TGraph(10);
   TGraph* gy = new TGraph(10);
   Double_t xline[2];
   Double_t yline[2];

   for (int ialign=0; ialign<4; ++ialign)
   {
      cout<< "\nAlign sensor #" << ialign <<endl;

      par0_x[0] = sensor[ialign].r0_(0);
      par0_x[1] = sensor[ialign].R_(0,2);
      gnmin_x.SetParameters(par0_x);

      par0_y[0] = sensor[ialign].r0_(1);
      par0_y[1] = sensor[ialign].R_(1,2);
      gnmin_y.SetParameters(par0_x);

      cout<< "init: par_x(0) = " << par_x(0) << " par_x(1) = " << par_x(1) << " res2sum_x = " << res2sum_x <<endl;
      cout<< "init: par_y(0) = " << par_y(0) << " par_y(1) = " << par_y(1) << " res2sum_y = " << res2sum_y <<endl;

      gnmin_x.InitIteration();
      gnmin_y.InitIteration();

      for (int ievt=0; ievt<nevt; ++ievt)
      {
         Int_t np = 0;
         gx->SetPoint(np, beamspot_z, beamspot_x);
         gy->SetPoint(np, beamspot_z, beamspot_y);
         np++;

         for (int isensor=0; isensor<4; ++isensor) {
            if (isensor == ialign) continue;
            gx->SetPoint(np, zpos[isensor], xpos[isensor][ievt]);
            gy->SetPoint(np, zpos[isensor], ypos[isensor][ievt]);
            np++;
         }
         gx->Set(np);
         gy->Set(np);

         // fit the straight line
         gx->Fit("pol1", "Q", "goff");
         gx->GetFunction("pol1")->GetParameters(xline);
         gy->Fit("pol1", "Q", "goff");
         gy->GetFunction("pol1")->GetParameters(yline);
         // get the fitted ray
         Ray rayfit(xline[0]+xline[1]*beamspot_z, yline[0]+yline[1]*beamspot_z, beamspot_z, xline[0]+xline[1]*zpos[0], yline[0]+yline[1]*zpos[0], zpos[0]);
         if (ialign == 0) {
            h2foil->Fill(rayfit.xat(beamspot_z), rayfit.yat(beamspot_z));
         }
         h2sensor[ialign]->Fill(rayfit.xat(sensor[ialign].r0_(2)), rayfit.yat(sensor[ialign].r0_(2)));

         //correct only sensor[0] for now

         // apply the fitted ray to the sensor with current iteration to get the ux
         TVectorD xpoint(3);
         xpoint = sensor[ialign].LocalCross(rayfit);
         Double_t ucross = xpoint(0);
         Double_t vcross = xpoint(1);

         factors[0] = 1.;
         factors[1] = rayfit.cx/rayfit.cz;
         gnmin_x.Feed(ucross, xpos[ialign][ievt], factors);
         factors[1] = rayfit.cy/rayfit.cz;
         gnmin_y.Feed(vcross, ypos[ialign][ievt], factors);

         hr[ialign]->Fill(ucross - xpos[ialign][ievt]);
         hLocalCross[ialign]->Fill(ucross);
         hxpos[ialign]->Fill(xpos[ialign][ievt]);
      }

      gnmin_x.Correct();
      gnmin_y.Correct();

      res2sum_x = 0;
      res2sum_y = 0;
      // for (int ievt=0; ievt<nevt; ++ievt) {
      //    Double_t r = gnmin.Residual(hgaus->GetBinCenter(i), hgaus->GetBinContent(i));
      //    res2sum += r*r;
      // }
      cout<< "par_x(0) = " << par_x(0) << " par_x(1) = " << par_x(1) << " res2sum_x = " << res2sum_x << " gnmin_x.OldSumResiduals() = " << gnmin_x.OldSumResiduals() <<endl;
      cout<< "par_y(0) = " << par_y(0) << " par_y(1) = " << par_y(1) << " res2sum_y = " << res2sum_y << " gnmin.OldSumResiduals() = " << gnmin_y.OldSumResiduals() <<endl;

      // correct current sensor position
      sensor[ialign].r0_(0) += par_x(0);     //-- NB sign: +
      sensor[ialign].AddBeta(par_x(1));      //-- NB: Beta for x
      sensor[ialign].r0_(1) += par_y(0);     //-- NB sign: +
      sensor[ialign].AddAlpha(par_y(1));     //-- NB: Alpha for y

      cout<< "sensor[" << ialign << "].r0_(0) = " << sensor[ialign].r0_(0) << " sensor[" << ialign << "].R_(0,2) = " << sensor[ialign].R_(0,2) <<endl;
      cout<< "sensor[" << ialign << "].r0_(1) = " << sensor[ialign].r0_(1) << " sensor[" << ialign << "].R_(1,2) = " << sensor[ialign].R_(1,2) <<endl;

      // correct the data values -------- not sure now :-)
      cout<< "correct data for position of the sensor #" << ialign << " sensor[" << ialign << "].r0_(0) = " << sensor[ialign].r0_(0) << " sensor[" << ialign << "].r0_(1) = " << sensor[ialign].r0_(1) <<endl;
      for (int ievt=0; ievt<nevt; ++ievt) {
         xpos[ialign][ievt] -= sensor[ialign].r0_(0);
         ypos[ialign][ievt] -= sensor[ialign].r0_(1);
      }
   }

   new TCanvas;
   h2hits->Draw();

   new TCanvas;
   h2foil->Draw();

   new TCanvas;
   h2sensor[0]->Draw();

   new TCanvas;
   h0->Draw();

   new TCanvas;
   hr[0]->Draw();

   new TCanvas;
   hLocalCross[0]->Draw();

   new TCanvas;
   hxpos[0]->Draw();

   // // play with last event

   // cout<< "test ToLocal and ToGlobal back" <<endl;

   // TVectorD vglobal0(3);
   // TVectorD vglobal(3);
   // TVectorD vlocal(3);
   // vglobal0(0) = xpos[0][nevt-1];
   // vglobal0(1) = ypos[0][nevt-1];
   // vglobal0(2) = sensorReal[0].GetPosition()(2);
   // vlocal = sensorReal[0].ToLocal(vglobal0);
   // vglobal = sensorReal[0].ToGlobal(vlocal);

   // cout<< "vglobal0";   vglobal0.Print();
   // cout<< "vglobal ";   vglobal.Print();

   // //-- rotate sensorReal[0]

   // const Double_t pi = TMath::Pi();
   // // sensorReal[0].SetRotation(pi/12, pi/9, pi/6);
   // // // sensorReal[0].SetRotation(1, 1, 0.01);
   // // cout<< "sensorReal[0].GetRotation().Print()"; sensorReal[0].GetRotation().Print();
   // //-- sensorReal[0].SetRotation(pi/6, 0, 0);
   // sensorReal[0].SetRotation(0, 0, pi/6);        // ok with theta=0, phi=0

   // TMatrixD RT(TMatrixD::kTransposed, sensorReal[0].GetRotation());
   // TMatrixD unity = RT*sensorReal[0].GetRotation();
   // cout<< "rotation matrix is unitary"; unity.Print();


   // cout<< "**********************************\nTest Ray" <<endl;
   // sensorReal[0].SetRotation(0, 0, 0);        // ok with theta=0, phi=0
   // //sensorReal[0].SetRotation(0, 0, pi/6);
   // // sensorReal[0].SetRotation(pi/6, 0, 0);
   // // sensorReal[0].SetRotation(0, pi/6, 0);
   // cout<< "sensorReal[0] is: " <<endl;
   // sensorReal[0].Print();

   // // generate cos(theta) in range (1-5e-5, 1)
   // Double_t delta = costheta_span*rand.Rndm(0);   // 1 - costheta
   // Double_t costheta = 1 - delta;
   // costheta = 1;                             //-- set explicit zero
   // costheta = TMath::Cos(TMath::ATan(19./3800.));

   // // generate phi
   // Double_t phi = 2.*TMath::Pi()*rand.Rndm(1);
   // phi = 0;                                  //-- set explicit zero

   // Ray ray(beamspot_x,beamspot_y,beamspot_z, costheta,phi);

   // // ray
   // cout<< "-- call ToLocal" <<endl;
   // Ray rloc = sensorReal[0].ToLocal(ray);
   // cout<< "-- call ToGlobal" <<endl;
   // Ray rglo = sensorReal[0].ToGlobal(rloc);
   // cout<< "ray:  " << ray <<endl;
   // cout<< "rloc: " << rloc <<endl;
   // cout<< "rglo: " << rglo <<endl;

   // TVectorD crossLoc(3);
   // crossLoc = sensorReal[0].CrossLoc(rloc);
   // cout<< "Local intersection point: crossLoc(0) = " << crossLoc(0) << " crossLoc(1) = " << crossLoc(1) << " crossLoc(2) = " << crossLoc(2) <<endl;

   // TVectorD cross(3);
   // cross = sensorReal[0].CrossGlo(ray);
   // cout<< "Global intersection point: cross(0) = " << cross(0) << " cross(1) = " << cross(1) << " cross(2) = " << cross(2) <<endl;

   // TVectorD localCross(3);
   // localCross = sensorReal[0].LocalCross(ray);
   // cout<< "Cross in the Local system: localCross(0) = " << localCross(0) << " localCross(1) = " << localCross(1) << " localCross(2) = " << localCross(2) <<endl;
}
