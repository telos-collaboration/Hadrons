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

    double a = par().alpha / 6.;
    std::vector<double> rho = {0, 0, 0, 0,
                               0, 0, a, a,
                               0, a, 0, a,
                               0, a, a, 0};

    Smear_APE<GImpl> smearer(rho);

    envGetTmp(GaugeField, buf);
    Usmr = U;

    LOG(Message) << "APE Smearing '" << par().gauge << "' for " << par().steps << " steps and alpha = " << par().alpha << " with initial plaquette = " << WilsonLoops<GImpl>::avgPlaquette(U) << std::endl;
    
    // Repeat the smearing nsteps times
    for (int i=0; i<par().steps; i++) {
        smearer.smear(buf, Usmr); // buf = rho * staples
        Usmr = (1-par().alpha) * Usmr + buf;
        
        GImpl::GaugeGroup::ProjectOnGeneralGroup(Usmr);

        LOG(Debug) << "Smearing step: " << i << " plaquette = " << WilsonLoops<GImpl>::avgPlaquette(Usmr) << std::endl;
    }

    LOG(Message) << "After " << par().steps << " smearing steps plaquette = " << WilsonLoops<GImpl>::avgPlaquette(Usmr) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_APE_Smearing_hpp_
