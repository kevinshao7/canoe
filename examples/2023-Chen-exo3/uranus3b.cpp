//========================================================================================                                        
// Athena++ astrophysical MHD code                                                                                                
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code                                                        
// contributors Licensed under the 3-clause BSD License, see LICENSE file for                                                     
// details                                                                                                                        
//========================================================================================                                        
//uranus.2a added realistic solar flux, start with high flux and decrease towards equilibrium

// C++ headers                                                                                                                    
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

// athena                                                                                                                         
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application                                                                                                                    
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe                                                                                                                          
#include <air_parcel.hpp>
#include <configure.hpp>
#include <impl.hpp>

// exo3                                                                                                                           
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// snap                                                                                                                           
#include <snap/thermodynamics/thermodynamics.hpp>

#define _sqr(x) ((x) * (x))
#define _qur(x) ((x) * (x) * (x) * (x))

using namespace std;

static Real p0, Rd, cp, Ts, Rpol,Req, grav, z_stra, eq_heat_flux,sday,syear, pi, dday, Omega, sigma,emissivity;
Real current_time;
std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);                                                                       

//REMEMBER VARIABLE GRAV IS STORED AS NEGATIVE!!!
// Calculate daily and annual cycle
Real Insolation(MeshBlock *pmb,int k,int j,Real time){
  auto pexo3 = pmb->pimpl->pexo3;
  int is = pmb->is;
  Real lat, lon;
  pexo3->GetLatLon(&lat, &lon, k, j, is);
  Real latsun = (pi*82.2/180)* cos(2*pi*time /syear); //in radians
  Real dayn = std::floor(time/dday);
  Real longsun = -2*pi*(time-dayn*dday)/dday; //in radians
  Real cosomega = sin(lat)*sin(latsun)+cos(lat)*cos(latsun)*cos(longsun-lon);
  if (cosomega>0){
    Real fluxmax = eq_heat_flux + 3000;
    return fluxmax*cosomega;
  }
  else return 0;
}

Real IsobarLevel(Real lat){//input lat in radian
  return Rpol*Req/sqrt(_sqr(Rpol)*_sqr(cos(lat))+_sqr(Req)*_sqr(sin(lat))); //chord of ellipse
}
Real BottomTemp(Real lat){
  return Ts + (-grav/cp)*(IsobarLevel(lat)-Rpol);
}

Real BottomPres(Real lat){//input lat in radians
  Real z1 = IsobarLevel(lat);
  Real lnbotpres = log(p0) - (cp/Rd)*log(1-(z1-Rpol)/(-cp*Ts/grav+z1-Rpol));
  return exp(lnbotpres);
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  Real om_earth = Omega;
  current_time=time;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        // coriolis force                                                                                                         
        Real f = 2. * om_earth * sin(lat);
        Real f2 = 2. * om_earth * cos(lat);
        Real U, V;

        pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);

        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * U;
        Real m3 = w(IDN, k, j, i) * V;

        // du(IM1, k, j, i) += dt * f * m2;                                                                                       
        Real ll_acc_U = f * m3;  //- f2 * m1;                                                                                     
        Real ll_acc_V = -f * m2;
        Real acc1, acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i);
        pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
        du(IM2, k, j, i) += dt * acc2;
        du(IM3, k, j, i) += dt * acc3;
      }

    // Heat flux at the boundaries                                                                                                
  for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) {
          // Energy input from bottom at 1 bar level
          for (int i = pmb->is; i<= pmb->ie; ++i){
            if (pmb->phydro->w(IPR, k, j, i) <= 1e5){
              Real dz = pmb->pcoord->dx1f(pmb->is);
              du(IEN, k, j, pmb->is) += Insolation(pmb,k,j,time) * dt / dz;
              break;
            }
          }
          int coolinglayers = 0;
          // Energy taken away from the top at 0.001 bar level                                                                      
          for (int i = pmb->ie; i>=0; --i){
            if (pmb->phydro->w(IPR, k, j, i) >= 100){
              Real dz = pmb->pcoord->dx1f(i);
              Real temp = pmb->phydro->w(IPR, k, j, pmb->ie) / pmb->phydro->w(IDN, k, j, pmb->ie) / Rd;
              if (temp >= 40){
                du(IEN, k, j, i) -= emissivity*sigma*_qur(temp) * dt / dz / 10; //emissivity of 0.25
                coolinglayers++;
              }
            }
            if (coolinglayers>=10){
              break;
            }
          }
          //Dirchlet bottom pressure condition
          Real lat, lon;
          pexo3->GetLatLon(&lat, &lon, k, j, pmb->is);
          pmb->phydro->w(IPR, k, j, pmb->is) = BottomPres(lat);
  }
}

Real AngularMomentum(MeshBlock *pmb, int iout) {
  auto pexo3 = pmb->pimpl->pexo3;
  Real AMz = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks,
      ke = pmb->ke;

  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        Real x1l = pmb->pcoord->x1f(i);
        Real x1u = pmb->pcoord->x1f(i + 1);
        Real U, V;
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        pexo3->GetUV(&U, &V, pmb->phydro->w(IVY, k, j, i),
                     pmb->phydro->w(IVZ, k, j, i), k, j, i);

        Real xt = tan(pmb->pcoord->x2v(j));
        Real yt = tan(pmb->pcoord->x3v(k));
        Real sin_theta =
            sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));

        Real x1 = tan(pmb->pcoord->x2f(j));
        Real x2 = tan(pmb->pcoord->x2f(j + 1));
        Real y = tan(pmb->pcoord->x3v(k));
        Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
        Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
        Real dx2_ang = acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));

        Real x = tan(pmb->pcoord->x2v(j));
        Real y1 = tan(pmb->pcoord->x3f(k));
        Real y2 = tan(pmb->pcoord->x3f(k + 1));
        delta1 = sqrt(1.0 + x * x + y1 * y1);
        delta2 = sqrt(1.0 + x * x + y2 * y2);
        Real dx3_ang = acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));

        Real vol = pmb->pcoord->dx1f(i) * dx2_ang * dx3_ang * sin_theta;

        // Originally here used cos(lat), which is x2v-pi, strange                                                                
        AMz += pmb->phydro->w(IDN, k, j, i) * vol *
               sqrt((_sqr(x1l) + _sqr(x1u)) / 2.) * cos(lat) *
               (Omega * sqrt(0.5 * (_sqr(x1l) + _sqr(x1u))) * cos(lat) + U);
      }
    }
  }

  return AMz;
}

//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)                                                                          
void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real day_to_s = 8.64E4;
  // forcing parameters                                                                                                           
  // thermodynamic parameters                                                                                                     
  Real gamma = pin->GetReal("hydro", "gamma");
  grav = pin->GetReal("hydro", "grav_acc1");
  Ts = pin->GetReal("problem", "Ts");
  p0 = pin->GetReal("problem", "p0");
  z_stra = pin->GetReal("thermodynamics", "z_stra");
  Rd = pin->GetReal("thermodynamics", "Rd");
  eq_heat_flux = pin->GetReal("problem", "eq_heat_flux");
  sday = pin->GetReal("problem", "sday");
  syear = pin->GetReal("problem", "syear");
  emissivity = pin->GetReal("problem", "emissivity");
  Rpol = pin->GetReal("problem", "Rpol");
  Req = pin->GetReal("problem","Req");
  cp = gamma / (gamma - 1.) * Rd;
  pi = 3.14159265;
  dday = syear/(1+syear/sday); //diurnal day inseconds
  Omega = -2*pi/sday; //retrograde
  sigma = 5.670374419E-8;
  current_time=0.;
  // forcing function                                                                                                             
  EnrollUserExplicitSourceFunction(Forcing);
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, AngularMomentum, "z-angular-mom");
}

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)                                                                     
//  \brief Held-Suarez problem generator                               
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pexo3 = pimpl->pexo3;
  Real grav = -pin->GetReal("hydro", "grav_acc1");
  Real gamma = pin->GetReal("hydro", "gamma");
  p0 = pin->GetReal("problem", "p0");
  Ts = pin->GetReal("problem", "Ts");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma / (gamma - 1.) * Rd;

  // construct an isothermal atmosphere                                                                                           
  auto pthermo = Thermodynamics::GetInstance();
  AirParcel air(AirParcel::Type::MoleFrac);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real lat, lon;
      pexo3->GetLatLon(&lat, &lon, k, j, is);
      air.w[IDN] = BottomTemp(lat);
      air.w[IPR] = BottomPres(lat);
      int i = is;
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "dry", grav, -grav/cp);
        // add noise                                                                                                              
        air.w[IVY] = 100* distribution(generator);
        air.w[IVZ] = 100* distribution(generator);
      }
    }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(7);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "lat");
  SetUserOutputVariableName(3, "lon");
  SetUserOutputVariableName(4, "vlat");
  SetUserOutputVariableName(5, "vlon");
  SetUserOutputVariableName(6, "solflux");
}


  // \brif Output distributions of temperature and potential temperature.                                                           
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pexo3 = pimpl->pexo3;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real prim[NHYDRO];
        for (int n = 0; n < NHYDRO; ++n) prim[n] = phydro->w(n, j, i);
        Real temp = phydro->w(IPR, k, j, i) / phydro->w(IDN, k, j, i) / Rd;
        user_out_var(0, k, j, i) = temp;
        user_out_var(1, k, j, i) =
            temp * pow(p0 / phydro->w(IPR, k, j, i), Rd / cp);
        Real lat, lon;
        Real U, V;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        pexo3->GetUV(&U, &V, phydro->w(IVY, k, j, i), phydro->w(IVZ, k, j, i),
                     k, j, i);
        user_out_var(2, k, j, i) = lat;
        user_out_var(3, k, j, i) = lon;
        user_out_var(4, k, j, i) = U;
        user_out_var(5, k, j, i) = V;
        Real latsun = (pi*82.2/180)* cos(2*pi* current_time /syear); //in radians
        Real dayn = std::floor(current_time/dday);
        Real longsun = -2*pi*(current_time-dayn*dday)/dday; //in radians
        Real cosomega = sin(lat)*sin(latsun)+cos(lat)*cos(latsun)*cos(longsun-lon);
        if (cosomega>0){
          Real fluxmax = eq_heat_flux + 3000;
          user_out_var(6,k,j,i) = fluxmax*cosomega;
        }
        else user_out_var(6,k,j,i) = 0;
      }
}