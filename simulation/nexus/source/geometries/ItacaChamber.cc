// ----------------------------------------------------------------------------
// nexus | ItacaChamber.cc
//
// Gas cylinder for ITACA track studies: enriched xenon with a helium
// admixture, no internal structures. Axis along z; the gas spans z = 0
// (TopMetal pad plane) to z = length (EL entry), centred on x = y = 0.
//
// With shield_thickness > 0 the gas is surrounded by the inner copper shield
// (ICS) of the radiogenic background runs (see ItacaChamber.h). Only the
// fiducial gas is sensitive: the energy saved is the energy deposited there.
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
#include <Randomize.hh>

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
    shield_thickness_(0.), shield_sampling_length_(0.), shield_material_("G4_Cu"),
    barrel_gap_(0.), barrel_layer_thickness_(0.),
    anode_gap_(0.), anode_layer_thickness_(0.),
    cathode_gap_(0.), cathode_layer_thickness_(0.),
    barrel_layer_material_("G4_Galactic"), anode_layer_material_("G4_Galactic"),
    cathode_layer_material_("G4_Galactic"),
    cavity_radius_(0.), cavity_zmin_(0.), cavity_zmax_(0.),
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

    DeclareShieldingCommands();
  }


  ItacaChamber::~ItacaChamber()
  {
    delete gas_vertex_gen_;
    delete msg_;
  }


  void ItacaChamber::DeclareShieldingCommands()
  {
    auto length = [this](const char* name, G4double& value, const char* doc) {
      G4GenericMessenger::Command& cmd = msg_->DeclareProperty(name, value, doc);
      cmd.SetUnitCategory("Length");
      cmd.SetParameterName(name, false);
      cmd.SetRange((std::string(name) + ">=0.").c_str());
    };
    length("shield_thickness", shield_thickness_,
           "Copper thickness on every face of the cavity; 0 (default) builds no shield.");
    length("shield_sampling_length", shield_sampling_length_,
           "Region ICS: depth from the inner face sampled as exp(-depth/length); 0 gives a uniform depth.");
    msg_->DeclareProperty("shield_material", shield_material_, "NIST material of the shield (G4_Cu).");
    length("barrel_gap", barrel_gap_, "Xenon between the gas and the barrel layer.");
    length("barrel_layer_thickness", barrel_layer_thickness_, "Passive layer between the barrel gap and the shield; 0 for none.");
    msg_->DeclareProperty("barrel_layer_material", barrel_layer_material_, "NIST material of the barrel layer.");
    length("anode_gap", anode_gap_, "Xenon between the gas (z = length) and the anode layer.");
    length("anode_layer_thickness", anode_layer_thickness_, "Passive layer between the anode gap and the shield; 0 for none.");
    msg_->DeclareProperty("anode_layer_material", anode_layer_material_, "NIST material of the anode layer.");
    length("cathode_gap", cathode_gap_, "Xenon between the gas (z = 0) and the cathode layer.");
    length("cathode_layer_thickness", cathode_layer_thickness_, "Passive layer between the cathode gap and the shield; 0 for none.");
    msg_->DeclareProperty("cathode_layer_material", cathode_layer_material_, "NIST material of the cathode layer.");
  }


  G4Material* ItacaChamber::FindMaterial(const G4String& name)
  {
    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(name);
    if (!mat)
      G4Exception("[ItacaChamber]", "FindMaterial()", FatalException,
                  ("Unknown NIST material: " + name).c_str());
    return mat;
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


  G4LogicalVolume* ItacaChamber::ConstructShielding(G4LogicalVolume* lab_logic, G4Material* gas)
  {
    const G4String name = "ITACA";
    const G4double t = shield_thickness_;

    // Cavity: the same gas as the fiducial volume, from the cathode layer to
    // the anode layer. Its radius and z range are the inner faces of the copper.
    cavity_radius_ = radius_ + barrel_gap_ + barrel_layer_thickness_;
    cavity_zmin_ = -(cathode_gap_ + cathode_layer_thickness_);
    cavity_zmax_ = length_ + anode_gap_ + anode_layer_thickness_;
    const G4double cavity_length = cavity_zmax_ - cavity_zmin_;
    const G4ThreeVector cavity_centre(0., 0., (cavity_zmin_ + cavity_zmax_) / 2.);

    G4Tubs* cavity_solid = new G4Tubs(name + "_CAVITY", 0., cavity_radius_, cavity_length / 2., 0., twopi);
    G4LogicalVolume* cavity_logic = new G4LogicalVolume(cavity_solid, gas, name + "_CAVITY");
    new G4PVPlacement(nullptr, cavity_centre, cavity_logic, name + "_CAVITY", lab_logic, false, 0, true);

    // Passive layers inside the cavity, positions relative to the cavity centre.
    if (barrel_layer_thickness_ > 0.) {
      G4Tubs* solid = new G4Tubs(name + "_BARREL_LAYER", radius_ + barrel_gap_, cavity_radius_, length_ / 2., 0., twopi);
      G4LogicalVolume* logic = new G4LogicalVolume(solid, FindMaterial(barrel_layer_material_), name + "_BARREL_LAYER");
      new G4PVPlacement(nullptr, G4ThreeVector(0., 0., length_ / 2.) - cavity_centre, logic, name + "_BARREL_LAYER", cavity_logic, false, 0, true);
    }
    if (anode_layer_thickness_ > 0.) {
      G4Tubs* solid = new G4Tubs(name + "_ANODE_LAYER", 0., cavity_radius_, anode_layer_thickness_ / 2., 0., twopi);
      G4LogicalVolume* logic = new G4LogicalVolume(solid, FindMaterial(anode_layer_material_), name + "_ANODE_LAYER");
      new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cavity_zmax_ - anode_layer_thickness_ / 2.) - cavity_centre, logic, name + "_ANODE_LAYER", cavity_logic, false, 0, true);
    }
    if (cathode_layer_thickness_ > 0.) {
      G4Tubs* solid = new G4Tubs(name + "_CATHODE_LAYER", 0., cavity_radius_, cathode_layer_thickness_ / 2., 0., twopi);
      G4LogicalVolume* logic = new G4LogicalVolume(solid, FindMaterial(cathode_layer_material_), name + "_CATHODE_LAYER");
      new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cavity_zmin_ + cathode_layer_thickness_ / 2.) - cavity_centre, logic, name + "_CATHODE_LAYER", cavity_logic, false, 0, true);
    }

    // Copper: barrel shell over the cavity length, and two disks closing it.
    G4Material* copper = FindMaterial(shield_material_);
    const G4double outer = cavity_radius_ + t;
    G4Tubs* barrel_solid = new G4Tubs(name + "_ICS_BARREL", cavity_radius_, outer, cavity_length / 2., 0., twopi);
    G4LogicalVolume* barrel_logic = new G4LogicalVolume(barrel_solid, copper, name + "_ICS_BARREL");
    new G4PVPlacement(nullptr, cavity_centre, barrel_logic, name + "_ICS_BARREL", lab_logic, false, 0, true);
    G4Tubs* disk_solid = new G4Tubs(name + "_ICS_DISK", 0., outer, t / 2., 0., twopi);
    G4LogicalVolume* cathode_disk = new G4LogicalVolume(disk_solid, copper, name + "_ICS_CATHODE");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cavity_zmin_ - t / 2.), cathode_disk, name + "_ICS_CATHODE", lab_logic, false, 0, true);
    G4LogicalVolume* anode_disk = new G4LogicalVolume(disk_solid, copper, name + "_ICS_ANODE");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cavity_zmax_ + t / 2.), anode_disk, name + "_ICS_ANODE", lab_logic, false, 0, true);

    const G4double volume = pi * (outer*outer - cavity_radius_*cavity_radius_) * cavity_length + 2. * pi * outer*outer * t;
    G4cout << "[ItacaChamber] shield " << shield_material_ << ": thickness " << t/cm << " cm, inner radius "
           << cavity_radius_/cm << " cm, cavity z " << cavity_zmin_/cm << " to " << cavity_zmax_/cm
           << " cm, mass " << volume * copper->GetDensity() / kg << " kg, region ICS sampling length "
           << shield_sampling_length_/cm << " cm" << G4endl;
    return cavity_logic;
  }


  void ItacaChamber::Construct()
  {
    const G4String name = "ITACA";
    const G4double half_length = length_ / 2.;
    const G4bool shielded = shield_thickness_ > 0.;

    // World: vacuum box. Particles leaving the gas (or the shield) are not returned.
    G4double lab_half = 1.2 * std::max(radius_, length_);
    if (shielded) {
      const G4double r = radius_ + barrel_gap_ + barrel_layer_thickness_ + shield_thickness_;
      const G4double z = std::max(cathode_gap_ + cathode_layer_thickness_,
                                  length_ + anode_gap_ + anode_layer_thickness_) + shield_thickness_;
      lab_half = 1.2 * std::max(r, z);
    }
    G4Box* lab_solid = new G4Box(name + "_LAB", lab_half, lab_half, lab_half);
    G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4LogicalVolume* lab_logic = new G4LogicalVolume(lab_solid, vacuum, name + "_LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    GeometryBase::SetLogicalVolume(lab_logic);

    G4Material* gas = XenonHelium(pressure_, temperature_, helium_mole_fraction_);

    // Detector coordinates: z = 0 at the TopMetal pad plane.
    const G4ThreeVector gas_centre(0., 0., half_length);
    G4LogicalVolume* mother = lab_logic;
    G4ThreeVector gas_position = gas_centre;
    if (shielded) {
      mother = ConstructShielding(lab_logic, gas);
      gas_position = gas_centre - G4ThreeVector(0., 0., (cavity_zmin_ + cavity_zmax_) / 2.);
    }

    G4Tubs* gas_solid = new G4Tubs(name + "_GAS", 0., radius_, half_length, 0., twopi);
    G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, gas, name + "_GAS");
    new G4PVPlacement(nullptr, gas_position, gas_logic, name + "_GAS", mother, false, 0, true);

    IonizationSD* ionizsd = new IonizationSD("/ITACA_GAS");
    G4SDManager::GetSDMpointer()->AddNewDetector(ionizsd);
    gas_logic->SetSensitiveDetector(ionizsd);

    gas_logic->SetUserLimits(new G4UserLimits(step_limit_));

    // Built here, after the macro has set the dimensions. World coordinates.
    delete gas_vertex_gen_;
    gas_vertex_gen_ = new CylinderPointSampler(0., radius_, half_length, 0., twopi,
                                               nullptr, gas_centre);

    G4cout << "[ItacaChamber] radius " << radius_/cm << " cm, length " << length_/cm
           << " cm (z = 0 to " << length_/cm << " cm), step limit " << step_limit_/mm << " mm" << G4endl;
  }


  G4ThreeVector ItacaChamber::GenerateVertex(const G4String& region) const
  {
    if (region == "GAS" || region == "VOLUME" || region == "XENON")
      return gas_vertex_gen_->GenerateVertex(VOLUME);
    if (region == "CENTER")
      return G4ThreeVector(0., 0., 0.);

    // Uniform on the inner surface of the gas cylinder (barrel and both end
    // faces, weighted by area), one micrometre inside the gas: the starting
    // point of gammas emitted by the materials around the gas.
    if (region == "INNER_SURFACE") {
      const G4double inset = 1.*micrometer;
      const G4double phi = twopi * G4UniformRand();
      G4double r, z;
      if (G4UniformRand() < length_ / (length_ + radius_)) {   // barrel: 2 pi R L of 2 pi R (L + R)
        r = radius_ - inset;
        z = inset + (length_ - 2.*inset) * G4UniformRand();
      } else {                                                 // one of the two end faces
        r = (radius_ - inset) * std::sqrt(G4UniformRand());
        z = (G4UniformRand() < 0.5) ? inset : length_ - inset;
      }
      return G4ThreeVector(r * std::cos(phi), r * std::sin(phi), z);
    }

    // In the copper shield. The piece (barrel, cathode disk, anode disk) is
    // chosen by volume; the depth from its inner face is sampled with density
    // exp(-depth/lambda) truncated at the thickness (uniform when lambda = 0);
    // the position on the surface at that depth is uniform. The weight that
    // restores a uniform distribution in the copper is computed from the saved
    // vertex (src/background/Background.jl, weight).
    if (region == "ICS") {
      if (shield_thickness_ <= 0.)
        G4Exception("[ItacaChamber]", "GenerateVertex()", FatalException,
                    "Region ICS needs shield_thickness > 0.");
      const G4double t = shield_thickness_, inner = cavity_radius_, outer = inner + t;
      const G4double cavity_length = cavity_zmax_ - cavity_zmin_;
      const G4double v_barrel = pi * (outer*outer - inner*inner) * cavity_length;
      const G4double v_disk = pi * outer*outer * t;
      const G4double lambda = shield_sampling_length_;
      const G4double depth = (lambda > 0.)
        ? -lambda * std::log(1. - G4UniformRand() * (1. - std::exp(-t / lambda)))
        : t * G4UniformRand();
      const G4double u = G4UniformRand() * (v_barrel + 2. * v_disk);
      const G4double phi = twopi * G4UniformRand();
      G4double r, z;
      if (u < v_barrel) {
        r = inner + depth;
        z = cavity_zmin_ + cavity_length * G4UniformRand();
      } else {
        r = outer * std::sqrt(G4UniformRand());
        z = (u < v_barrel + v_disk) ? cavity_zmin_ - depth : cavity_zmax_ + depth;
      }
      return G4ThreeVector(r * std::cos(phi), r * std::sin(phi), z);
    }

    G4Exception("[ItacaChamber]", "GenerateVertex()", FatalException,
                ("Unknown vertex region: " + region).c_str());
    return G4ThreeVector();
  }

} // end namespace nexus
