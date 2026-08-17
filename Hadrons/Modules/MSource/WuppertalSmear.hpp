#ifndef Hadrons_MSource_WuppertalSmear_hpp_
#define Hadrons_MSource_WuppertalSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         WuppertalSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class WuppertalSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WuppertalSmearPar,
                                    std::string, gauge,
                                    double, step,
                                    int, iterations,
                                    int, orthog,
                                    std::string, source);
};

template <typename FImpl>
class TWuppertalSmear: public Module<WuppertalSmearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TWuppertalSmear(const std::string name);
    // destructor
    virtual ~TWuppertalSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WuppertalSmear, TWuppertalSmear<FIMPL>, MSource);

/******************************************************************************
 *                 TWuppertalSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWuppertalSmear<FImpl>::TWuppertalSmear(const std::string name)
: Module<WuppertalSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWuppertalSmear<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWuppertalSmear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWuppertalSmear<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envTmp(std::vector<GaugeMat>, "Umu", 1, 4, envGetGrid(LatticeColourMatrix));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWuppertalSmear<FImpl>::execute(void)
{
    LOG(Message) << "Wuppertal Smearing starting..." << std::endl;
    auto &out = envGet(PropagatorField, getName());
    auto &src = envGet(PropagatorField, par().source);
    auto &U = envGet(GaugeField, par().gauge);
    envGetTmp(std::vector<GaugeMat>, Umu);
    for(int mu=0; mu<Nd; mu++)
    {
       Umu.at(mu)=peekLorentz(U,mu);
    }
    WuppertalSmearing<FImpl> wupsmear;
    out=src;
    startTimer("Wuppertal iteration");
    wupsmear.WuppertalSmear(Umu, out, par().step, par().iterations, par().orthog);
    stopTimer("Wuppertal iteration");
    LOG(Message) << "Wuppertal Smearing ending..." << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_WuppertalSmear_hpp_
