/*
 * Test_field_io.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Gaurav Ray       <gsr95@pm.me>
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
#include <Hadrons/Modules.hpp>
#include <Hadrons/Application.hpp>

using namespace Grid;
using namespace Hadrons;

/*
    This test demonstrates the functionality of
    the ILDG IO modules in Hadrons.

    It does 4 things:
      1) Generates a cfg using MGauge::Unit/Random
      2) Saves to disk that cfg using MIO::SaveIldg
      3) Loads the saved cfg back using MIO::LoadIldg
      4) Performs a computation on the loaded cfg
         using MGradientFlow::WilsonFlow

    Except for the inclusion of steps 2 and 3 in the
    same programme, this layout is intended to parallel the
    layout of a typical Hadrons application.

    MIO::SaveIldgPar::gauge is a string that should
    match the output string of the Module producing 
    the cfg that is to be saved to disk.

    MIO::SaveIldgPar::gaugeGroup is a string that can
    only be either su or sp. Take care when setting this
    because if set incorrectly saved fields can be corrupted
    and data irretrievably lost.

    MIO::LoadIldgPar::waitForSave is a string that lets
    MIO::LoadIldg know whether it needs to wait for
    the cfg to be saved to disk before attempting to load
    it into memory. If not set Hadrons will
    expect the cfgs to already be on disk.
*/

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    Application application;

    Application::GlobalPar          globalPar; 
    MIO::SaveIldgPar                saveIldgPar;
    MIO::LoadIldgPar                loadIldgPar;
    MGradientFlow::WilsonFlow::Par  gfPar;

    globalPar.runId             = "TEST_ILDG_IO";
    globalPar.trajCounter.start = 1;
    globalPar.trajCounter.end   = 2;
    globalPar.trajCounter.step  = 1;

    application.setPar(globalPar);

    std::vector<std::string> groups = {"su"};  // only su and sp valid

    // generate cfgs
    application.createModule<MGauge::Random>("su-test-cfg"); // SU(Nc)

#if Sp2n_config == 1
    application.createModule<MGauge::Unit>("sp-test-cfg");
    groups.push_back("sp");
#endif 

    // demonstrate a save-load-measure workflow
    for(auto &g: groups) {
      saveIldgPar.gauge         = g + "-test-cfg";    // name of gauge field
      saveIldgPar.ensembleId    = "telos";
      saveIldgPar.ensembleLabel = g + std::to_string(Nc) + "hadrons_test";
      saveIldgPar.gaugeGroup    = g;

      saveIldgPar.precision     = "double"; // default value
      saveIldgPar.reducedFormat = false;    // default value
      saveIldgPar.fileStem      = g+"_full_double";
      // save
      application.createModule<MIO::SaveIldg>("save-"+g+"-full-double", saveIldgPar);

      loadIldgPar.file        = saveIldgPar.fileStem;
      loadIldgPar.waitForSave = "Y";        // don't set if cfg already on disk
      // load
      application.createModule<MIO::LoadIldg>("load-"+g+"-full-double", loadIldgPar);

      // gfPar.gauge must match the name of the 
      // corresponding LoadIldg module
      gfPar.gauge         = "load-"+g+"-full-double";
      gfPar.steps         = 6;
      gfPar.step_size     = 0.01;
      gfPar.meas_interval = 2;
      // measure
      application.createModule<MGradientFlow::WilsonFlow>(g+"-wilson-flow",gfPar);

      // demonstrate 3 other write/read combinations
      // that Grid's IldgWriter/Reader can handle
      saveIldgPar.precision     = "single";
      saveIldgPar.reducedFormat = false;
      saveIldgPar.fileStem      = g+"_full_single";
      application.createModule<MIO::SaveIldg>("save-"+g+"-full-single", saveIldgPar);

      loadIldgPar.file      = saveIldgPar.fileStem;
      application.createModule<MIO::LoadIldg>("load-"+g+"-full-single", loadIldgPar);

      saveIldgPar.precision     = "double";
      saveIldgPar.reducedFormat = true;
      saveIldgPar.fileStem      = g+"_red_double";
      application.createModule<MIO::SaveIldg>("save-"+g+"-reduced-double", saveIldgPar);

      loadIldgPar.file      = saveIldgPar.fileStem;
      application.createModule<MIO::LoadIldg>("load-"+g+"-reduced-double", loadIldgPar);

      saveIldgPar.precision     = "single";
      saveIldgPar.fileStem      = g+"_red_single";
      application.createModule<MIO::SaveIldg>("save-"+g+"-reduced-single", saveIldgPar);

      loadIldgPar.file      = saveIldgPar.fileStem;
      application.createModule<MIO::LoadIldg>("load-"+g+"-reduced-single", loadIldgPar);

    }

    application.run();

    Grid_finalize();
    
    return EXIT_SUCCESS;
}
