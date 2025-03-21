//========================================================================================                                        
// Athena++ astrophysical MHD code                                                                                                
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code                                                        
// contributors Licensed under the 3-clause BSD License, see LICENSE file for                                                     
// details                                                                                                                        
//========================================================================================                                        
//uranus2p Newtonian cooling to specified temperature profile, including latitudinal variation

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
#define _cube(x) ((x) * (x) * (x))

using namespace std;

static Real p0, Rd, cp, Ts, Rp, grav, eq_heat_flux,sday,syear, pi, dday, Omega, sigma, emissivity,adtdz,cv,spinupflux,flux_ratio,Tint,initheatdecay,Kt,sponge_tau,spongeheight;
std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);                                                                       

// Real Insolation(MeshBlock *pmb,int k,int j,Real time){
//   auto pexo3 = pmb->pimpl->pexo3;
//   int is = pmb->is;
//   Real lat, lon;
//   Real shifttime = time + (syear/4);
//   pexo3->GetLatLon(&lat, &lon, k, j, is);
//   Real latsun = (pi*82.2/180)* cos(2*pi*shifttime /syear); //in radians
//   Real dayn = std::floor(shifttime/dday);
//   Real longsun = 2*pi*(shifttime-dayn*dday)/dday; //in radians
//   Real cosomega = sin(lat)*sin(latsun)+cos(lat)*cos(latsun)*cos(longsun-lon);
//   if (cosomega>0){
//     if (time > 5.E9){
//       return eq_heat_flux*cosomega;
//     }
//     else{
//       Real fluxmax = eq_heat_flux + spinupflux/exp( time /initheatdecay);
//       return fluxmax*cosomega;
//     }
    
//   }
//   else return 0;
// }
Real power(Real base, int exp){
  if (exp == 0) return 1;
  Real result = base;
  for (int k = 1; k < exp; k++){
    result *= base;
  }
  return result;
}

Real Tempprof(Real logp,Real lat,Real coeffs[]){
  Real logbar = logp-5;
  //log p in bar, lat in rad
  Real bestfit = coeffs[0];
  int index = 1;
  for (int i= 1; i<=6; i++){
    for (int j = 0; j<= i; j++){
      bestfit += coeffs[index]*power(lat,j)*power(logbar,(i-j));
      index++;
    }
  }
  return bestfit;
}


void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  Real om_earth = Omega;
  Real coeffs[28] = {77.85794651631285,63.05354059921049,-0.6579587841707587,48.82531413927198,-0.2797132026545594,-4.841788384525834,5.17455075960518,-0.24855224651931832,0.758499268400661,1.0035081120559388,-11.442993252368204,0.2562222422845137,0.5665583712287425,0.5280296497013807,5.591761024224716,-6.1371175259807575,0.2651917494349661,0.29129966012989605,-0.26069951559528914,-0.583875113022543,-0.19537113253113358,-0.9062494030895678,0.051886033275524344,0.06898914727838346,-0.0873511826991891,-0.22860085803293506,-0.1581405250484245,-1.5798139118321535};


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

    // Heating   
  for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i){
          Real logpress = log10(pmb->phydro->w(IPR, k, j,i));
          Real lat, lon;
          pexo3->GetLatLon(&lat, &lon, k, j, i);
          Real Teq = Tempprof(logpress,lat,coeffs);
          Real Temp =  pmb->phydro->w(IPR, k, j,i) / pmb->phydro->w(IDN, k, j,i) / Rd;
          du(IEN, k, j, i) += dt*(cp - Rd) * w(IDN, k, j, i) * Kt * (Temp-Teq);
          Real z = pmb->pcoord->x1v(i) - Rp;
          if (z > spongeheight) {  // sponge layer at top and bottom
            Real tau = sponge_tau;
            du(IVX, k, j, i) -= w(IVX, k, j, i) * (dt / tau) * w(IDN, k, j, i);
            du(IVY, k, j, i) -= w(IVY, k, j, i) * (dt / tau) * w(IDN, k, j, i);
            du(IVZ, k, j, i) -= w(IVZ, k, j, i) * (dt / tau) * w(IDN, k, j, i);
        }
      }
        

    // if(u(IEN, k, j, pmb->ie)<0) u(IEN, k, j, pmb->ie) = 0;                                                                      
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
  Rd = pin->GetReal("thermodynamics", "Rd");
  eq_heat_flux = pin->GetReal("problem", "eq_heat_flux");
  sday = pin->GetReal("problem", "sday");
  syear = pin->GetReal("problem", "syear");
  Tint = pin->GetReal("problem", "Tint");
  cp = gamma / (gamma - 1.) * Rd;
  cv = cp/gamma;
  pi = 3.14159265;
  dday = syear/(1+syear/sday); //diurnal day inseconds
  Omega = -2*pi/sday; //retrograde
  sigma = 5.670374419E-8;
  // flux_ratio = 0.25*1.14*eq_heat_flux/(sigma*_qur(Ts));
  // int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  // Real currentforcing[ke-ks][je-js][ie-is] = {};
  Kt = pin->GetReal("problem", "Kt"); //specific cooling/heating rate
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
  Rp = pin->GetReal("problem", "Rp");
  emissivity = pin->GetReal("problem", "emissivity");
  eq_heat_flux = pin->GetReal("problem", "eq_heat_flux");
  adtdz = pin->GetReal("problem", "adtdz");
  spinupflux = pin->GetReal("problem", "spinupflux");
  initheatdecay = pin->GetReal("problem", "initheatdecay");
  sponge_tau = pin->GetReal("problem", "sponge_tau");
  spongeheight = pin->GetReal("problem", "spongeheight");
  Kt = pin->GetReal("problem", "Kt"); //specific cooling/heating rate
  cp = gamma / (gamma - 1.) * Rd;
  pi = 3.14159265;
  dday = syear/(1+syear/sday); //diurnal day inseconds
  Omega = -2*pi/sday; //retrograde
  sigma = 5.670374419E-8;

  // construct an isothermal atmosphere                                                                                           
  auto pthermo = Thermodynamics::GetInstance();
  AirParcel air(AirParcel::Type::MoleFrac);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.w[IPR] = p0;
      air.w[IDN] = Tint;

      int i = is;
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "dry", grav, adtdz);
        // add noise                                                                                                              
        air.w[IVY] = 0.1* distribution(generator);
        air.w[IVZ] = 0.1* distribution(generator);
      }
    }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(8);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "lat");
  SetUserOutputVariableName(3, "lon");
  SetUserOutputVariableName(4, "vlat");
  SetUserOutputVariableName(5, "vlon");
  SetUserOutputVariableName(6, "forcing");
  SetUserOutputVariableName(7, "goalprof");
}


  // \brif Output distributions of temperature and potential temperature.                                                           
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  Real coeffs[28] = {77.85794651631285,63.05354059921049,-0.6579587841707587,48.82531413927198,-0.2797132026545594,-4.841788384525834,5.17455075960518,-0.24855224651931832,0.758499268400661,1.0035081120559388,-11.442993252368204,0.2562222422845137,0.5665583712287425,0.5280296497013807,5.591761024224716,-6.1371175259807575,0.2651917494349661,0.29129966012989605,-0.26069951559528914,-0.583875113022543,-0.19537113253113358,-0.9062494030895678,0.051886033275524344,0.06898914727838346,-0.0873511826991891,-0.22860085803293506,-0.1581405250484245,-1.5798139118321535};
  auto pexo3 = pimpl->pexo3;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real prim[NHYDRO];
        for (int n = 0; n < NHYDRO; ++n) prim[n] = phydro->w(n, j, i);
        Real temp = phydro->w(IPR, k, j, i) / phydro->w(IDN, k, j, i) / Rd; //assume hydrostatic??
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
        Real logpress = log10(phydro->w(IPR, k, j,i));
        Real Teq = Tempprof(logpress,lat,coeffs);
        Real Temp =  phydro->w(IPR, k, j,i) / phydro->w(IDN, k, j,i) / Rd;
        Real parcelmass = pcoord->GetCellVolume(k,j,i)*phydro->w(IDN, k, j,i);
        user_out_var(6,k,j,i) = (cp - Rd) * phydro->w(IDN, k, j, i) * Kt * (Temp-Teq);
        user_out_var(7,k,j,i) = Teq;
        //user_out_var(6,k,j,i) = Tempprof(logpress,lat);
    }
}

