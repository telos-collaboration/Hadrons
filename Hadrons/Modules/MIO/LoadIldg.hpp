#ifndef Hadrons_MIO_LoadIldg_hpp_
#define Hadrons_MIO_LoadIldg_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadIldg                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadIldgPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadIldgPar,
                                    std::string, fileStem,
                                    bool,        waitForSave);
};

template <typename FImpl>
class TLoadIldg: public Module<LoadIldgPar>
{
public:
    // constructor
    TLoadIldg(const std::string name);
    // destructor
    virtual ~TLoadIldg(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadIldg, TLoadIldg<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadIldg implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadIldg<FImpl>::TLoadIldg(const std::string name)
: Module<LoadIldgPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadIldg<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    // if cfg not already on disk a MIO::SaveIldg module must first write it
    if( par().waitForSave )
    {
        in = {par().fileStem};
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadIldg<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadIldg<FImpl>::setup(void)
{
    envCreateLat(LatticeGaugeField, getName()); 
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadIldg<FImpl>::execute(void)
{
    FieldMetaData header;
    std::string   fileName = par().fileStem + "."
                             + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading ILDG gauge field from file '" << fileName
                 << "'" << std::endl;

    auto &U = envGet(LatticeGaugeField, getName());
    IldgReader _IldgReader;
   _IldgReader.open(fileName);
   _IldgReader.readConfiguration(U,header);
   _IldgReader.close();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadIldg_hpp_
