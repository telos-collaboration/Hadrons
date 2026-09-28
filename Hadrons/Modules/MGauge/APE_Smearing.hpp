#ifndef Hadrons_MGauge_APE_Smearing_hpp_
#define Hadrons_MGauge_APE_Smearing_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         APE_Smearing                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class APE_SmearingPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(APE_SmearingPar, std::string, gauge, double,
                                  alpha, unsigned int, steps);
};

template <typename GImpl> class TAPE_Smearing : public Module<APE_SmearingPar> {
public:
  GAUGE_TYPE_ALIASES(GImpl, );

public:
  // constructor
  TAPE_Smearing(const std::string name);
  // destructor
  virtual ~TAPE_Smearing(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(APE_Smearing, TAPE_Smearing<FIMPL>, MGauge);

/******************************************************************************
 *                 TAPE_Smearing implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TAPE_Smearing<GImpl>::TAPE_Smearing(const std::string name)
    : Module<APE_SmearingPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TAPE_Smearing<GImpl>::getInput(void) {
  std::vector<std::string> in = {par().gauge};

  return in;
}

template <typename GImpl>
std::vector<std::string> TAPE_Smearing<GImpl>::getOutput(void) {
  std::vector<std::string> out = {getName()};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl> void TAPE_Smearing<GImpl>::setup(void) {
  envCreateLat(GaugeField, getName());
  envTmpLat(GaugeField, "buf");
}

// execution ///////////////////////////////////////////////////////////////////
/*! @brief APE Smearing.
 *
 * The implementation is based on the definition given in
 * E. Bennett et al., "Meson spectroscopy from spectral densities in lattice
 * gauge theories", Phys. Rev. D 110, 074509. See Eq. (22) in section IIIA:
 * APE and Wuppertal smearing algorithms.
 *
 * The staple operator is defined as
 * \f[S_\mu(x) \equiv \sum_{\nu \neq \mu} U_\nu(x), U_\mu(x+\hat\nu),
 * U_\nu^\dagger(x+\hat\mu).\f]
 *
 * The APE smearing iteration is
 * \f[ U_\mu^{(m)}(x) = P\left\{(1-\alpha_{\mathrm{APE}})U_\mu^{(m-1)}(x) +
 * \frac{\alpha_{\mathrm{APE}}}{2(N_{\mathrm{smear}}-1)} S_\mu^{(m-1)}(x)
 * \right\},\f]
 *
 * where \f$P\f$ denotes the projection back onto the gauge group,
 * \f$\alpha_{\mathrm{APE}}\f$ is the APE-smearing step size, and
 * \f$N_{\mathrm{smear}}\f$ is the number of spacetime directions included
 * in the smearing procedure.
 */
template <typename GImpl> void TAPE_Smearing<GImpl>::execute(void) {
  auto &U = envGet(GaugeField, par().gauge);
  auto &Usmr = envGet(GaugeField, getName());

  double a =
      par().alpha /
      4.; // alpha/(2*(Nd_smear - 1)) where Nd_smear is the number of dimensions
          // that the field is being smeared in. For eg., for spatial
          // smearing, Nd_smear = 3 after excluding temporal dimension
  std::vector<double> rho = {0, a, a, 0, a, 0, a, 0, a, a, 0, 0, 0, 0, 0, 0};

  Smear_APE<GImpl> smearer(rho);

  envGetTmp(GaugeField, buf);
  Usmr = U;

  LOG(Message) << "APE Smearing '" << par().gauge << "' for " << par().steps
               << " steps and alpha = " << par().alpha
               << " with initial plaquette = "
               << WilsonLoops<GImpl>::avgPlaquette(U) << std::endl;

  // Repeat the smearing nsteps times
  for (int i = 1; i <= par().steps; i++) {
    smearer.smear(buf, Usmr); // buf = rho * staples

    for (int mu = 0; mu < Nd - 1;
         mu++) // updates only the spatial components, leaves temporal
               // components untouched change 'a' and 'rho' matrix for general
               // smearing
    {
      auto U_mu = PeekIndex<LorentzIndex>(Usmr, mu);
      auto B_mu = PeekIndex<LorentzIndex>(buf, mu);

      auto new_mu = (1.0 - par().alpha) * U_mu + B_mu;

      PokeIndex<LorentzIndex>(Usmr, new_mu, mu);
    }

    auto U_t = PeekIndex<LorentzIndex>(Usmr, Nd - 1);

    startTimer("Projection timer");
    GImpl::GaugeGroup::ProjectOnSpecialGroup(Usmr);
    stopTimer("Projection timer");

    PokeIndex<LorentzIndex>(Usmr, U_t, Nd - 1);
  }

  LOG(Message) << "After " << par().steps
               << " smearing steps, total plaquette = "
               << WilsonLoops<GImpl>::avgPlaquette(Usmr) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_APE_Smearing_hpp_
