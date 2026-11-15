// ----------------------------------------------------------------------------
// nexus | ItacaChamber.cc
//
// Gas cylinder for ITACA track studies: enriched xenon with a helium
// admixture, no internal structures. Axis along z, centred at the origin.
//
// NAUSICAA addition to the soleti/nexus "precdr" branch.
// ----------------------------------------------------------------------------

#include "ItacaChamber.h"

#include "CylinderPointSampler.h"
#include "IonizationSD.h"
#include "XenonProperties.h"
#include "FactoryBase.h"

#include <G4GenericMessenger.hh>
#include <G4UserLimits.hh>
#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4Element.hh>
#include <G4Isotope.hh>
#include <G4VisAttributes.hh>
#include <G4SDManager.hh>
#include <G4Exception.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(ItacaChamber, GeometryBase)

namespace nexus {

  ItacaChamber::ItacaChamber():
    GeometryBase(),
    radius_(180.*cm), length_(130.*cm),
    pressure_(16.*bar), temperature_(290.*kelvin),
    helium_mole_fraction_(0.15), step_limit_(1.*mm),
    gas_vertex_gen_(nullptr)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/ItacaChamber/",
      "Control commands of geometry ItacaChamber.");

    G4GenericMessenger::Command& radius_cmd =
      msg_->DeclareProperty("radius", radius_, "Radius of the gas cylinder.");
    radius_cmd.SetUnitCategory("Length");
    radius_cmd.SetParameterName("radius", false);
    radius_cmd.SetRange("radius>0.");

    G4GenericMessenger::Command& length_cmd =
      msg_->DeclareProperty("length", length_, "Length of the gas cylinder.");
    length_cmd.SetUnitCategory("Length");
    length_cmd.SetParameterName("length", false);
    length_cmd.SetRange("length>0.");

    G4GenericMessenger::Command& pressure_cmd =
      msg_->DeclareProperty("pressure", pressure_, "Total gas pressure.");
    pressure_cmd.SetUnitCategory("Pressure");
    pressure_cmd.SetParameterName("pressure", false);
    pressure_cmd.SetRange("pressure>0.");

    G4GenericMessenger::Command& temperature_cmd =
      msg_->DeclareProperty("temperature", temperature_, "Gas temperature.");
    temperature_cmd.SetUnitCategory("Temperature");
    temperature_cmd.SetParameterName("temperature", false);
    temperature_cmd.SetRange("temperature>0.");

    G4GenericMessenger::Command& helium_cmd =
      msg_->DeclareProperty("helium_mole_fraction", helium_mole_fraction_,
                            "Helium mole fraction, 0 to 1 (0 gives pure xenon).");
    helium_cmd.SetParameterName("helium_mole_fraction", false);
    helium_cmd.SetRange("helium_mole_fraction>=0. && helium_mole_fraction<1.");

    G4GenericMessenger::Command& step_cmd =
      msg_->DeclareProperty("step_limit", step_limit_, "Maximum step in the gas.");
    step_cmd.SetUnitCategory("Length");
    step_cmd.SetParameterName("step_limit", false);
    step_cmd.SetRange("step_limit>0.");
  }


  ItacaChamber::~ItacaChamber()
  {
    delete gas_vertex_gen_;
    delete msg_;
  }


  G4Material* ItacaChamber::XenonHelium(G4double pressure, G4double temperature,
                                        G4double helium_mole_fraction)
  {
    G4String name = "GXeHe_ITACA";
    G4Material* mat = G4Material::GetMaterial(name, false);
    if (mat) return mat;

    // Enriched xenon: isotopic abundances of materials::GXeEnriched
    const std::vector<std::pair<G4int, G4double>> isotopes = {
      {129, 0.0656392*perCent}, {130, 0.0656392*perCent}, {131, 0.234361*perCent},
      {132, 0.708251*perCent},  {134, 8.6645*perCent},    {136, 90.2616*perCent}
    };

    G4double xenon_molar_mass = 0.;
    G4Element* xenon = new G4Element("EnrichedXe_ITACA", "Xe", isotopes.size());
    for (const auto& iso : isotopes) {
      G4double molar_mass = XenonMassPerMole(iso.first);
      xenon_molar_mass += iso.second * molar_mass;
      xenon->AddIsotope(new G4Isotope("Xe" + std::to_string(iso.first) + "_ITACA",
                                      54, iso.first, molar_mass), iso.second);
    }

    const G4double helium_molar_mass = 4.0026*g/mole;
    G4Element* helium = new G4Element("Helium_ITACA", "He", 1);
    helium->AddIsotope(new G4Isotope("He4_ITACA", 2, 4, helium_molar_mass), 100.*perCent);

    // Xenon compressibility: second virial coefficient, as in
    // materials::CompressibilityFactor (NEXT DocDB 1183). Helium: ideal gas.
    const G4double gas_constant = k_Boltzmann * Avogadro;
    const G4double T_r = temperature / (289.733*kelvin);
    const G4double P_r = pressure / (58.420*bar);
    const G4double z_xenon = 1. + (0.083 - 0.422/std::pow(T_r, 1.6)) * (P_r/T_r);

    const G4double xenon_density  = pressure * xenon_molar_mass  / (z_xenon * gas_constant * temperature);
    const G4double helium_density = pressure * helium_molar_mass / (gas_constant * temperature);

    // Additive partial volumes (Amagat): one mole of mixture occupies
    // (x_xe z_xe + x_he) R T / P, so the density of a mixture with mole
    // fractions x_xe and x_he is
    //   rho = P (x_xe M_xe + x_he M_he) / ((x_xe z_xe + x_he) R T).
    const G4double x_he = helium_mole_fraction;
    const G4double x_xe = 1. - x_he;
    const G4double density =
      pressure * (x_xe * xenon_molar_mass + x_he * helium_molar_mass) /
      ((x_xe * z_xenon + x_he) * gas_constant * temperature);

    // Geant4 takes MASS fractions when the fraction is a G4double.
    const G4double xenon_mass  = x_xe * xenon_molar_mass;
    const G4double helium_mass = x_he * helium_molar_mass;
    const G4double w_xe = xenon_mass / (xenon_mass + helium_mass);
    const G4double w_he = 1. - w_xe;

    const G4int components = (x_he > 0.) ? 2 : 1;
    mat = new G4Material(name, density, components, kStateGas, temperature, pressure);
    mat->AddElement(xenon, w_xe);
    if (x_he > 0.) mat->AddElement(helium, w_he);

    G4cout << "[ItacaChamber] gas: He mole fraction " << x_he
           << ", He mass fraction " << w_he
           << ", density " << density/(kg/m3) << " kg/m3"
           << " (Xe " << xenon_density/(kg/m3) << ", He " << helium_density/(kg/m3) << " kg/m3"
           << " at " << pressure/bar << " bar, " << temperature/kelvin << " K)" << G4endl;

    return mat;
  }


  void ItacaChamber::Construct()
  {
    const G4String name = "ITACA";
    const G4double half_length = length_ / 2.;

    // World: vacuum box. Particles leaving the gas are not returned.
    const G4double lab_half = 1.2 * std::max(radius_, half_length);
    G4Box* lab_solid = new G4Box(name + "_LAB", lab_half, lab_half, lab_half);
    G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4LogicalVolume* lab_logic = new G4LogicalVolume(lab_solid, vacuum, name + "_LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    GeometryBase::SetLogicalVolume(lab_logic);

    G4Tubs* gas_solid = new G4Tubs(name + "_GAS", 0., radius_, half_length, 0., twopi);
    G4Material* gas = XenonHelium(pressure_, temperature_, helium_mole_fraction_);
    G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, gas, name + "_GAS");
    new G4PVPlacement(nullptr, G4ThreeVector(), gas_logic, name + "_GAS", lab_logic, false, 0, true);

    IonizationSD* ionizsd = new IonizationSD("/ITACA_GAS");
    G4SDManager::GetSDMpointer()->AddNewDetector(ionizsd);
    gas_logic->SetSensitiveDetector(ionizsd);

    gas_logic->SetUserLimits(new G4UserLimits(step_limit_));

    // Built here, after the macro has set the dimensions.
    delete gas_vertex_gen_;
    gas_vertex_gen_ = new CylinderPointSampler(0., radius_, half_length, 0., twopi);

    G4cout << "[ItacaChamber] radius " << radius_/cm << " cm, length " << length_/cm
           << " cm, step limit " << step_limit_/mm << " mm" << G4endl;
  }


  G4ThreeVector ItacaChamber::GenerateVertex(const G4String& region) const
  {
    if (region == "GAS" || region == "VOLUME" || region == "XENON")
      return gas_vertex_gen_->GenerateVertex(VOLUME);
    if (region == "CENTER")
      return G4ThreeVector(0., 0., 0.);

    G4Exception("[ItacaChamber]", "GenerateVertex()", FatalException,
                ("Unknown vertex region: " + region).c_str());
    return G4ThreeVector();
  }

} // end namespace nexus
