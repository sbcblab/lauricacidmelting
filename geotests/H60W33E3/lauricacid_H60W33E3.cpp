/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Alexandre Tacques, Tim N. Bingert
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* The solution for the melting problem (solid-liquid phase change)
coupled with natural convection is found using the lattice Boltzmann
method after Rongzong Huang and Huiying Wu (2015)[1]. The equilibrium
distribution function for the temperature is modified in order to deal
with the latent-heat source term. That way, iteration steps or solving
a group of linear equations is avoided, which results in enhanced efficiency.
The phase interface is located by the current total enthalpy, and
its movement is considered by the immersed moving boundary scheme after
Noble and Torczynski (1998)[2]. This method was validated by comparison
with experimental values (e.g. Gau and Viskanta (1986) [3]).

[1] Maximilian Gaedtke, S. Abishek, Ryan Mead-Hunter, Andrew J.  C. King, Benjamin J. Mullins, Hermann Nirschl, Mathias J. Krause,
Total enthalpy-based lattice Boltzmann simulations of melting in paraffin/metal foam composite phase change materials,
International Journal of Heat and Mass Transfer, Volume 155, 2020, 119870, ISSN 0017-9310.

[2] Rongzong Huang, Huiying Wu, Phase interface effects in the total enthalpy-based lattice
Boltzmann model for solid–liquid phase change, Journal of Computational Physics 294 (2015) 346–362.

[3] D. Noble, J. Torczynski, A lattice-Boltzmann method for partially saturated
computational cells, Int. J. Modern Phys. C 9 (8) (1998) 1189–1202.

[4] C. Gau, R. Viskanta, Melting and Solidification of a Pure Metal on a
Vertikal Wall, Journal of Heat Transfer  108(1) (1986): 174–181.
 */

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using T = FLOATING_POINT_TYPE;
using NSDESCRIPTOR = D2Q9<POROSITY,VELOCITY_SOLID,FORCE>;
using TDESCRIPTOR  = D2Q5<VELOCITY,TEMPERATURE>;
using TotalEnthalpyAdvectionDiffusionDynamics = TotalEnthalpyAdvectionDiffusionBGKdynamics<T,TDESCRIPTOR>;

// Parameters for the simulation setup
const T lx  = 50e-3;          // length of the cavity
const T ly  = 120e-3;         // height of the cavity
const T finlength = 33.33e-3;    // length of the fin
const T finheight = 3e-3;    // length of the fin 
int     N = 100;              // resolution of the model
T       tau = 0.51;           // relaxation time
const T maxPhysT = 13000.;    // simulation time
const T density = 885.;       // density [kg / m³]
const T dyn_viscosity = 7.6e-3;    // dynamic viscosity [kg / m s] 
const T thermalexpcoeff = 6.15e-4; // thermal expansion coefficient [1 / K]
const T gravity = 9.81;            // gravity [m / s²]
const T Tcold = 25.;               // Temperature of heated wall [Celsius]
const T Tmelt = 47.5;              // Temperature of heated wall [Celsius]
const T Thot = 70.0;               // Temperature of heated wall [Celsius]
const T lambda_s = 0.158;          // Thermal conductivity [W / m K]
const T lambda_l = 0.145;          // Thermal conductivity [W / m K]
const T lambda_ref = (2 * lambda_s * lambda_l / (lambda_s + lambda_l)); // Harmonic mean
const T R_lambda =  lambda_s/lambda_l;
// Specific heat
const T cp_s = 2180.;                               // [J / kg K] OR [J / kg C] OR [m² / s² K]
const T cp_l = 2390.;                               // [J / kg K] OR [J / kg C] OR [m² / s² K]
const T cp_ref = 2.0 * cp_s * cp_l / (cp_s + cp_l); // Harmonic mean
const T L_phys = (187210.);                         // Latent Heat of Lauric Acid 187210 [J / kg] = [m² / s²]
// Hydraulic diameter for rectangular cavity 
// https://www.sciencedirect.com/topics/engineering/hydraulic-diameter
const T charac_length = (4*((lx*ly)-(finlength*finheight)))/(lx+lx+ly+ly+finlength+finlength); // charac_length = 4 x Area of PCM / Perimeter with fin included
const T thermaldiff = ( lambda_ref / (density * cp_ref) );
const T kin_viscosity = dyn_viscosity / density;
// Ra - "Strength of thermal buoyancy against the viscous and thermal diffusion"  
const T Ra = (gravity  * thermalexpcoeff * (Thot - Tmelt) * util::pow(charac_length, 3)) / (thermaldiff * kin_viscosity);
// Pr - Ratio of momentum diffusivity to thermal diffusivity
const T Pr = kin_viscosity / thermaldiff;
T char_lattice_u = 0.03422;   // characteristical lattice velocity (obtained empirically)
const T char_u = util::sqrt( gravity  * thermalexpcoeff * (Thot - Tmelt) * density ); // square root of F buoyancy force?
const T conversion_u = char_u / char_lattice_u;
T Hcold, Hhot, Hcold_lattice, Hhot_lattice;
T physDeltaX, physDeltaT;

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,2>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,3);
  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2,T());
  origin[0] = converter.getPhysLength(1);
  origin[1] = 0.5*converter.getPhysLength(1);
  IndicatorCuboid2D<T> cuboid2(extend, origin);
  superGeometry.rename(3,1,cuboid2);
  std::vector<T> extendwallright(2,T(0));
  extendwallright[0] = converter.getPhysLength(1);
  extendwallright[1] = ly+2*converter.getPhysLength(1);
  std::vector<T> originwallright(2,T(0));
  originwallright[0] = lx+2*converter.getPhysLength(1);
  originwallright[1] = 0.0;
  IndicatorCuboid2D<T> wallright(extendwallright, originwallright);
  std::vector<T> extendfin(2,T(0));
  extendfin[0] = 34.33e-3; 
  extendfin[1] = 3e-3;
  std::vector<T> originfin(2,T(0));
  originfin[0] = 16.67e-3;
  originfin[1] = 58.5e-3;
  IndicatorCuboid2D<T> fin(extendfin, originfin);
  superGeometry.rename(3,2,wallright);
  superGeometry.rename(1,2,fin);
  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;

}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice<T, TDESCRIPTOR>& ADlattice,
                     SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  
  //////////////////////////////////////////////// Convert from physical to lattice units ///////////////////////////////////////////////////////////////
  // Adimensional parameters that regulate the physics behaviour
  clout << "charac_length=" << charac_length << std::endl;
  clout << "Ra=" << Ra << " Pr=" << Pr << std::endl;
  // Temperatures
  clout << "Tcold=" << Tcold << " Thot=" << Thot << " Tmelt=" << Tmelt << std::endl;
  T Tcold_lattice = Tcold / converter.getConversionFactorTemperature();  
  T Thot_lattice = Thot / converter.getConversionFactorTemperature();   
  T Tmelt_lattice = Tmelt / converter.getConversionFactorTemperature();
  clout << "Tcold_lattice=" << Tcold_lattice << " Thot_lattice=" << Thot_lattice << " Tmelt_lattice=" << Tmelt_lattice << std::endl;
  // Density
  T density_lattice = density / converter.getConversionFactorDensity();
  clout << "density=" << density << " density_lattice=" << density_lattice << std::endl;
  // Specific heat      
  clout << "cp_s=" << cp_s << " cp_l=" << cp_l << " cp_ref=" << cp_ref << std::endl;
  T cp_s_lattice = cp_s / (conversion_u*conversion_u/converter.getConversionFactorTemperature());
  T cp_l_lattice = cp_l / (conversion_u*conversion_u/converter.getConversionFactorTemperature());
  T cp_ref_lattice = cp_ref / (conversion_u*conversion_u/converter.getConversionFactorTemperature());
  clout << "char_lattice_u=" << char_lattice_u << " char_u=" << char_u << " conversion_u=" << conversion_u << std::endl;  
  clout << "cp_s_lattice=" << cp_s_lattice << " cp_l_lattice=" << cp_l_lattice << " cp_ref_lattice=" << cp_ref_lattice << std::endl; 
  // Enthalpy 
  Hcold = cp_s * Tcold;
  Hhot = cp_l * Thot;
  clout << "H_cold=" << Hcold << " H_hot=" << Hhot << std::endl;
  Hcold_lattice = cp_s_lattice * Tcold_lattice;
  Hhot_lattice = cp_l_lattice * Thot_lattice;
  clout << "Hcold_lattice=" << Hcold_lattice << " Hhot_lattice=" << Hhot_lattice << std::endl;
  // Latent heat 
  T L_lattice = L_phys / conversion_u / conversion_u;  
  clout << "latentheat= " << L_phys << " latentheatinlatticeunits= " << L_lattice << std::endl;
  //////////////////////////////////////////////// Convert from physical to lattice units ///////////////////////////////////////////////////////////////

  T omega  = converter.getLatticeRelaxationFrequency();
  T Tomega = converter.getLatticeThermalRelaxationFrequency();
  /// sets dynamics in each lattice and bounceback dynamic in the adiabatic wall  
  NSlattice.defineDynamics<ForcedPSMBGKdynamics>(superGeometry.getMaterialIndicator({1, 2, 3}));
  ADlattice.defineDynamics<TotalEnthalpyAdvectionDiffusionDynamics>(superGeometry.getMaterialIndicator({1, 2, 3}));
  setBounceBackBoundary(ADlattice, superGeometry, 3);
  /// sets boundary
  setRegularizedTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry.getMaterialIndicator({2}));
  setInterpolatedVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry.getMaterialIndicator({2, 3}));
  NSlattice.setParameter<descriptors::OMEGA>(omega);
  ADlattice.setParameter<descriptors::OMEGA>(Tomega);
  ADlattice.setParameter<TotalEnthalpy::T_S>(Tmelt_lattice);
  ADlattice.setParameter<TotalEnthalpy::T_L>(Tmelt_lattice);
  ADlattice.setParameter<TotalEnthalpy::CP_S>(cp_s_lattice);
  ADlattice.setParameter<TotalEnthalpy::CP_L>(cp_l_lattice);   
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_S>(cp_ref_lattice / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5) * R_lambda);
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_L>(cp_ref_lattice / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5));
  ADlattice.setParameter<TotalEnthalpy::L>(L_lattice);
  /// define initial conditions  
  AnalyticalConst2D<T,T> rho(density_lattice);
  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> H_cold(Hcold_lattice);
  AnalyticalConst2D<T,T> H_hot(Hhot_lattice);
  AnalyticalConst2D<T,T> tauAD(1./converter.getLatticeThermalRelaxationFrequency());
  ADlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3}), tauAD );
  AnalyticalConst2D<T,T> tauNS(1./converter.getLatticeRelaxationFrequency());
  NSlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3}), tauNS );
  
  AnalyticalConst2D<T,T> T_hot_(Thot_lattice);  
  ADlattice.defineField<descriptors::TEMPERATURE>( superGeometry.getMaterialIndicator({2}), T_hot_ );
    
  /// for each material set Rho, U and the Equilibrium
  NSlattice.defineRhoU(superGeometry.getMaterialIndicator({1, 2, 3}), rho, u0);
  NSlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 3}), rho, u0);

  ADlattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({1, 2}), u0);
  ADlattice.defineRho(superGeometry.getMaterialIndicator({1}), H_cold);
  ADlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1}), H_cold, u0);
  ADlattice.defineRho(superGeometry, 2, H_hot);
  ADlattice.iniEquilibrium(superGeometry, 2, H_hot, u0);

  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice, int iT,
                 SuperGeometry<T,2>& superGeometry,
                 util::Timer<T>& timer,
                 bool converged)
{

  OstreamManager clout(std::cout,"getResults");

    SuperVTMwriter2D<T> vtkWriter("lauricacid");
    SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);

    SuperLatticeField2D<T, TDESCRIPTOR, VELOCITY> velocity(ADlattice);
    velocity.getName() = "velocity";
    SuperLatticePhysPressure2D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
    pressure.getName() = "pressure";
    SuperLatticeDensity2D<T, TDESCRIPTOR> enthalpy(ADlattice);
    enthalpy.getName() = "enthalpy";
    SuperLatticeField2D<T, NSDESCRIPTOR, POROSITY> liquid_frac(NSlattice);
    liquid_frac.getName() = "liquid fraction";
    SuperLatticeField2D<T, TDESCRIPTOR, TEMPERATURE> temperature(ADlattice);
    temperature.getName() = "temperature";
    SuperLatticeField2D<T, NSDESCRIPTOR, FORCE> force(NSlattice);
    force.getName() = "force";    
    SuperLatticePhysTemperature2D<T, NSDESCRIPTOR, TDESCRIPTOR> PhysTemperature_(ADlattice, converter);
    AnalyticalConst2D<T,T> ConversionTemperature_(converter.getConversionFactorTemperature());
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> ConversionTemperature(ConversionTemperature_, ADlattice);
    SuperIdentity2D<T,T> PhysTemperature(temperature * ConversionTemperature);
    PhysTemperature.getName()="PhysTemperature";    
    SuperLatticePhysPressure2D<T, NSDESCRIPTOR> PhysPressure(NSlattice, converter);
    PhysPressure.getName() = "PhysPressure";
    SuperLatticeField2D<T, TDESCRIPTOR, TAU_EFF> tauAD(ADlattice);
    tauAD.getName() = "AD tau";
    SuperLatticeField2D<T, NSDESCRIPTOR, TAU_EFF> tauNS(NSlattice);
    tauNS.getName() = "NS tau";
    vtkWriter.addFunctor( geometry );
    vtkWriter.addFunctor( pressure );
    vtkWriter.addFunctor( enthalpy );
    vtkWriter.addFunctor( liquid_frac );
    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( force );
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( tauAD );
    vtkWriter.addFunctor( tauNS );
    vtkWriter.addFunctor( PhysTemperature );
    vtkWriter.addFunctor( PhysPressure );    
    
    const int vtkIter = 1000000.;
    //const int vtkIter = 1000.;

    if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
        SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
        SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
        SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
        vtkWriter.write(geometry);
        vtkWriter.write(cuboid);
        vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT % vtkIter == 0 || converged) {

    timer.update(iT);
    timer.printStep();

    /// NSlattice statistics console output
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    /// ADlattice statistics console output
    ADlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    /*if ( NSlattice.getStatistics().getAverageRho() != NSlattice.getStatistics().getAverageRho() or ADlattice.getStatistics().getAverageRho() != ADlattice.getStatistics().getAverageRho() ) {
      clout << "simulation diverged! stopping now." << std::endl;
      exit(1);
    }*/

    vtkWriter.write(converter.getPhysTime(iT)*100);
  }
}

int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
    
  if (argc >= 2) {
    N = atoi(argv[1]);
  }
  if (argc >= 3) {
    tau = atof(argv[2]);
  }
  if (argc >= 4) {
    char_lattice_u = atof(argv[3]);
  }

  physDeltaX = lx / N;
  physDeltaT = density / dyn_viscosity / descriptors::invCs2<T,NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX;

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const converter(
    (T) physDeltaX, // physDeltaX
    (T) physDeltaT, // physDeltaT
    (T) charac_length, // charPhysLength
    (T) char_u, // charPhysVelocity
    (T) dyn_viscosity / density, // physViscosity
    (T) density, // physDensity
    (T) lambda_l, // physThermalConductivity
    (T) cp_l, // physSpecificHeatCapacity
    (T) thermalexpcoeff, // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );
  converter.print();
  
  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx + 2*converter.getPhysLength(1);
  extend[1] = ly + converter.getPhysLength(1);
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidGeometry
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry(cuboidGeometry, loadBalancer);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);
  SuperLattice<T, NSDESCRIPTOR> NSlattice(superGeometry);

  T cp_ref_lattice = cp_ref / (conversion_u*conversion_u/converter.getConversionFactorTemperature());
    
  std::vector<T> dir{0.0, 1.0};
  T boussinesqForcePrefactor = Ra / util::pow(T(N),3) * Pr * util::pow(cp_ref_lattice / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5), 2);
  //clout << "boussinesq " << Ra / util::pow(T(N), 3) * Pr * lambda_l * lambda_l << std::endl;
  clout << "boussinesqForcePrefactor " << boussinesqForcePrefactor << std::endl;
  
  T Tcold_lattice = Tcold / converter.getConversionFactorTemperature();    
    
  TotalEnthalpyPhaseChangeCouplingGenerator2D<T,NSDESCRIPTOR,TotalEnthalpyAdvectionDiffusionDynamics>
  coupling(0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(ly),
           boussinesqForcePrefactor, Tcold_lattice, 1., dir);

  NSlattice.addLatticeCoupling(superGeometry, 1, coupling, ADlattice);

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  /// === 5th Step: Definition of Initial and Boundary Conditions ===
  // prepareLattice and setBoundaryConditions
  prepareLattice(converter, NSlattice, ADlattice, superGeometry);

  SuperLatticeCoupling couplingSmago(
    EnthalpySmagorinskyCoupling{},
    names::NavierStokes{}, NSlattice,
    names::Temperature{}, ADlattice);
  couplingSmago.setParameter<EnthalpySmagorinskyCoupling::SMAGORINSKY_PREFACTOR>(0.2);
  couplingSmago.setParameter<EnthalpySmagorinskyCoupling::PR_TURB>(Pr);
  couplingSmago.setParameter<EnthalpySmagorinskyCoupling::OMEGA_NSE>(1./converter.getLatticeRelaxationTime());
  couplingSmago.setParameter<EnthalpySmagorinskyCoupling::OMEGA_ADE>(1./converter.getLatticeThermalRelaxationTime());

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT)+1; ++iT) {
    couplingSmago.execute();

    /// === 6th Step: Collide and Stream Execution ===
    NSlattice.executeCoupling();
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, false);
   // if (iT % converter.getLatticeTime(0.5) == 0) {
   //   SuperLatticeExternalScalarField2D<T, NSDESCRIPTOR, POROSITY> liquid_frac(NSlattice);
   //   SuperAverage2D<T> liquidFrac_total_(liquid_frac, superGeometry, 1);
   //   int in[2];
   //   T out[1];
   //   liquidFrac_total_(out, in);
   //   std::cout << "liquid fraction: " << out[0] << std::endl;
   // }
  }

  timer.stop();
  timer.printSummary();

}
