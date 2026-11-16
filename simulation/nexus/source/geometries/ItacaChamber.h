// ----------------------------------------------------------------------------
// nexus | ItacaChamber.h
//
// Gas cylinder for ITACA track studies: enriched xenon with a helium
// admixture, no internal structures. Axis along z; the gas spans z = 0
// (TopMetal pad plane) to z = length (EL entry), centred on x = y = 0.
//
// NAUSICAA addition to the soleti/nexus "precdr" branch.
// ----------------------------------------------------------------------------

#ifndef ITACA_CHAMBER_H
#define ITACA_CHAMBER_H

#include "GeometryBase.h"

class G4Material;
class G4GenericMessenger;
namespace nexus { class CylinderPointSampler; }

namespace nexus {

  class ItacaChamber: public GeometryBase
  {
  public:
    ItacaChamber();
    ~ItacaChamber();

    /// Vertex uniform in the gas for region GAS (also VOLUME or XENON), or
    /// uniform on the inner surface of the gas cylinder for INNER_SURFACE
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

    /// Xe-He mixture: mole fractions converted to the mass fractions
    /// Geant4 requires; density from additive partial volumes (Amagat).
    static G4Material* XenonHelium(G4double pressure, G4double temperature,
                                   G4double helium_mole_fraction);

  private:
    G4double radius_;               ///< Gas radius (ITACA: 180 cm)
    G4double length_;               ///< Gas length along z (ITACA: 130 cm)
    G4double pressure_;             ///< Total pressure (16 bar)
    G4double temperature_;          ///< Gas temperature (290 K)
    G4double helium_mole_fraction_; ///< Helium mole fraction (0.15)
    G4double step_limit_;           ///< Maximum step in the gas (1 mm)

    CylinderPointSampler* gas_vertex_gen_; ///< Built in Construct()

    G4GenericMessenger* msg_;
  };

} // end namespace nexus

#endif
