/*
 * SaveIldg.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Gaurav Ray <gsr95@pm.me>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MIO_SaveIldg_hpp_
#define Hadrons_MIO_SaveIldg_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Save an ILDG configuration

 gauge         Name of the gauge field object to write
 fileStem      Namestem of the file to write the gauge field to
 ensembleLabel Label of the ensemble. Recommended this is gauge info.
 ensembleId    Collaboration Name
 gaugeGroup    type of field - either su or sp
 reducedFormat save gauge field in reduced or full format
 precision     save gauge field in single or double precision
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MIO)

class SaveIldgPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveIldgPar,
                                    std::string, gauge,
                                    std::string, fileStem,
                                    std::string, ensembleLabel,
                                    std::string, ensembleId,
                                    std::string, gaugeGroup,
                                    bool,        reducedFormat,
                                    std::string, precision);
};

template <typename GImpl>
class TSaveIldg: public Module<SaveIldgPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TSaveIldg(const std::string name);
    // destructor
    virtual ~TSaveIldg(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveIldg,  TSaveIldg<GIMPL>,  MIO);

/******************************************************************************
*                       TSaveIldg implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TSaveIldg<GImpl>::TSaveIldg(const std::string name)
: Module<SaveIldgPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TSaveIldg<GImpl>::getInput(void)
{
    return { par().gauge };
}

template <typename GImpl>
std::vector<std::string> TSaveIldg<GImpl>::getOutput(void)
{
    return { par().fileStem };
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TSaveIldg<GImpl>::setup(void)
{
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TSaveIldg<GImpl>::execute(void)
{
    std::string fileName = par().fileStem + "." + std::to_string(vm().getTrajectory());
    LOG(Message) << "DEBUG:: Saving ILDG configuration to file '" << fileName
                 << "'" << std::endl;

    auto &U = envGet(GaugeField, par().gauge);
    using stats = PeriodicGaugeStatistics;
    makeFileDir(fileName, U.Grid());
    IldgWriter _IldgWriter(U.Grid()->IsBoss());
    _IldgWriter.open(fileName);

    if( par().gaugeGroup == "su" ) {
        if( par().precision == "single" ) {
            if( par().reducedFormat ) {
                _IldgWriter.writeConfiguration<stats,GroupName::SU,MatrixFormat::REDUCED,FloatingPointFormat::IEEE32BIG>(U, vm().getTrajectory(), par().ensembleId, par().ensembleLabel);
            } else if ( !par().reducedFormat) {
                _IldgWriter.writeConfiguration<stats,GroupName::SU,MatrixFormat::FULL,FloatingPointFormat::IEEE32BIG>(U, vm().getTrajectory(), par().ensembleId, par().ensembleLabel);
            }
        } else if ( par().precision == "double" ) {
             if( par().reducedFormat ) {
                _IldgWriter.writeConfiguration<stats,GroupName::SU,MatrixFormat::REDUCED,FloatingPointFormat::IEEE64BIG>(U, vm().getTrajectory(), par().ensembleId, par().ensembleLabel);
            } else if ( !par().reducedFormat) {
                _IldgWriter.writeConfiguration<stats,GroupName::SU,MatrixFormat::FULL,FloatingPointFormat::IEEE64BIG>(U, vm().getTrajectory(), par().ensembleId, par().ensembleLabel);
            }
        }
    } 

#if Sp2n_config == 1
    if( par().gaugeGroup == "sp" ) {
        if( par().precision == "single" ) {
            if( par().reducedFormat ) {
                _IldgWriter.writeConfiguration<stats,GroupName::Sp,MatrixFormat::REDUCED,FloatingPointFormat::IEEE32BIG>(U, vm().getTrajectory(), par().ensembleId, par().ensembleLabel);
            } else if ( !par().reducedFormat) {
                _IldgWriter.writeConfiguration<stats,GroupName::Sp,MatrixFormat::FULL,FloatingPointFormat::IEEE32BIG>(U, vm().getTrajectory(), par().ensembleId, par().ensembleLabel);
            }
        } else if ( par().precision == "double" ) {
             if( par().reducedFormat ) {
                _IldgWriter.writeConfiguration<stats,GroupName::Sp,MatrixFormat::REDUCED,FloatingPointFormat::IEEE64BIG>(U, vm().getTrajectory(), par().ensembleId, par().ensembleLabel);
            } else if ( !par().reducedFormat) {
                _IldgWriter.writeConfiguration<stats,GroupName::Sp,MatrixFormat::FULL,FloatingPointFormat::IEEE64BIG>(U, vm().getTrajectory(), par().ensembleId, par().ensembleLabel);
            }
        }
    }
#endif

    _IldgWriter.close();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_SaveIldg_hpp_
