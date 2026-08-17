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

class APE_SmearingPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(APE_SmearingPar,
                                    std::string, gauge,
                                    double, alpha,
                                    unsigned int, steps);
};

template <typename GImpl>
class TAPE_Smearing: public Module<APE_SmearingPar>
{
public: 
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TAPE_Smearing(const std::string name);
    // destructor
    virtual ~TAPE_Smearing(void) {};
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
: Module<APE_SmearingPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TAPE_Smearing<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TAPE_Smearing<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TAPE_Smearing<GImpl>::setup(void)
{
   envCreateLat(GaugeField, getName());
   envTmpLat(GaugeField, "buf"); 
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TAPE_Smearing<GImpl>::execute(void)
{
    auto &U = envGet(GaugeField, par().gauge);
    auto &Usmr = envGet(GaugeField, getName());

    double a = par().alpha / 4.;
    std::vector<double> rho = {0, a, a, 0,
                               a, 0, a, 0,
                               a, a, 0, 0,
                               0, 0, 0, 0};

    Smear_APE<GImpl> smearer(rho);

    envGetTmp(GaugeField, buf);
    Usmr = U;

    std::vector<GaugeLinkField> Uinit(Nd, U.Grid());

    for(int mu = 0; mu < Nd; mu++) {
                Uinit[mu] = PeekIndex<LorentzIndex>(U, mu);
    }

    for(int mu = 0; mu < Nd; mu++){
          for (int nu = 0; nu < mu; nu++){

                        typename GImpl::ComplexField Plaq_init(U.Grid());

                        WilsonLoops<GImpl>::traceDirPlaquette(Plaq_init, Uinit, mu, nu);

                        auto Tp = sum(Plaq_init);
                        auto p = TensorRemove(Tp);

                        RealD avg_mu_nu = p.real() / U.Grid()->gSites() / Nc;

                        LOG(Debug) << "Initial directional plaquettes: " << " Plaquette[" << mu << "," << nu << "] = " << avg_mu_nu << std::endl;
          }
    }
    
    startTimer("Plaquette timer");
    LOG(Message) << "APE Smearing '" << par().gauge << "' for " << par().steps << " steps and alpha = " << par().alpha << " with initial plaquette = " << WilsonLoops<GImpl>::avgPlaquette(U) << std::endl;
    stopTimer("Plaquette timer");    

    // Repeat the smearing nsteps times
    for (int i=0; i<par().steps; i++) {
        smearer.smear(buf, Usmr); // buf = rho * staples

        for(int mu = 0; mu < Nd-1; mu++) //updates only the spatial components, leaves temporal components untouched
        {
                auto U_mu = PeekIndex<LorentzIndex>(Usmr, mu);
                auto B_mu = PeekIndex<LorentzIndex>(buf, mu);

                auto new_mu = (1.0 - par().alpha) * U_mu + B_mu;

                PokeIndex<LorentzIndex>(Usmr, new_mu, mu);
        }

        //Usmr = (1-par().alpha) * Usmr + buf;
	
	auto U_t = PeekIndex<LorentzIndex>(Usmr, Nd - 1);
        
        startTimer("Projection timer");
        GImpl::GaugeGroup::ProjectOnSpecialGroup(Usmr);
        stopTimer("Projection timer");

	PokeIndex<LorentzIndex>(Usmr, U_t, Nd - 1);

        LOG(Debug) << "Smearing step: " << i + 1 << " total plaquette = " << WilsonLoops<GImpl>::avgPlaquette(Usmr) << std::endl;
        LOG(Debug) << "Smearing step: " << i + 1 << " spatial plaquette = " << WilsonLoops<GImpl>::timesliceAvgSpatialPlaquette(Usmr)[0] << std::endl;

        std::vector<GaugeLinkField> Udir(Nd, Usmr.Grid());

        for(int mu = 0; mu < Nd; mu++) {
                Udir[mu] = PeekIndex<LorentzIndex>(Usmr, mu);
        }

        for(int mu = 0; mu < Nd; mu++){
                for (int nu = 0; nu < mu; nu++){

                        typename GImpl::ComplexField Plaq(Usmr.Grid());

                        WilsonLoops<GImpl>::traceDirPlaquette(Plaq, Udir, mu, nu);

                        auto Tp = sum(Plaq);
                        auto p = TensorRemove(Tp);

                        RealD avg_mu_nu = p.real() / Usmr.Grid()->gSites() / Nc;

                        LOG(Debug) << "Smearing step: " << i + 1 << " Plaquette[" << mu << "," << nu << "] = " << avg_mu_nu << std::endl;
                }
        }

    }

    LOG(Message) << "After " << par().steps << " smearing steps, total plaquette = " << WilsonLoops<GImpl>::avgPlaquette(Usmr) << std::endl;
    //LOG(Message) << "After " << par().steps << " smearing steps, spatial plaquette = " << spatial[0] << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_APE_Smearing_hpp_
