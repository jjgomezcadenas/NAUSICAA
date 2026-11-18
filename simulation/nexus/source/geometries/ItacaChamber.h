// ----------------------------------------------------------------------------
// nexus | ItacaChamber.h
//
// Gas cylinder for ITACA track studies: enriched xenon with a helium
// admixture, no internal structures. Axis along z; the gas spans z = 0
// (TopMetal pad plane) to z = length (EL entry), centred on x = y = 0.
//
// With shield_thickness > 0 the gas is surrounded by the inner copper shield
// (ICS) of the radiogenic background runs: a xenon cavity with one passive
// layer on each face (barrel, anode end at z = length, cathode end at z = 0),
// closed by a copper barrel and two copper disks of that thickness. Region ICS
// generates vertices in the copper with a depth sampled as exp(-depth/lambda).
//
// NAUSICAA addition to the soleti/nexus "precdr" branch.
// ----------------------------------------------------------------------------

#ifndef ITACA_CHAMBER_H
#define ITACA_CHAMBER_H

#include "GeometryBase.h"

class G4Material;
class G4GenericMessenger;
class G4LogicalVolume;
namespace nexus { class CylinderPointSampler; }

namespace nexus {

  class ItacaChamber: public GeometryBase
  {
  public:
    ItacaChamber();
    ~ItacaChamber();

    /// Vertex uniform in the gas for region GAS (also VOLUME or XENON),
    /// uniform on the inner surface of the gas cylinder for INNER_SURFACE,
    /// or in the copper shield for ICS (exponential depth, see the .cc)
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

    /// Xe-He mixture: mole fractions converted to the mass fractions
    /// Geant4 requires; density from additive partial volumes (Amagat).
    static G4Material* XenonHelium(G4double pressure, G4double temperature,
                                   G4double helium_mole_fraction);

  private:
    /// Cavity, passive layers and copper around the gas; returns the cavity
    /// logical volume, the mother of the gas
    G4LogicalVolume* ConstructShielding(G4LogicalVolume* lab_logic, G4Material* gas);
    void DeclareShieldingCommands();
    static G4Material* FindMaterial(const G4String& name);

    G4double radius_;               ///< Gas radius (ITACA: 180 cm)
    G4double length_;               ///< Gas length along z (ITACA: 130 cm)
    G4double pressure_;             ///< Total pressure (16 bar)
    G4double temperature_;          ///< Gas temperature (290 K)
    G4double helium_mole_fraction_; ///< Helium mole fraction (0.15)
    G4double step_limit_;           ///< Maximum step in the gas (1 mm)

    G4double shield_thickness_;        ///< Copper thickness on every face; 0 = no shield (default)
    G4double shield_sampling_length_;  ///< Depth sampling length of region ICS; 0 = uniform depth
    G4String shield_material_;         ///< Shield material (G4_Cu)
    G4double barrel_gap_, barrel_layer_thickness_;    ///< Xenon gap from the gas, then the layer
    G4double anode_gap_, anode_layer_thickness_;
    G4double cathode_gap_, cathode_layer_thickness_;
    G4String barrel_layer_material_, anode_layer_material_, cathode_layer_material_;

    G4double cavity_radius_;           ///< Inner radius of the copper, set in Construct()
    G4double cavity_zmin_, cavity_zmax_; ///< Inner faces of the copper disks, world z

    CylinderPointSampler* gas_vertex_gen_; ///< Built in Construct()

    G4GenericMessenger* msg_;
  };

} // end namespace nexus

#endif
