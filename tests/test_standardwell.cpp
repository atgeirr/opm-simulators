/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil Petroleum AS.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE StandardWellTest

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>





#include <ewoms/common/propertysystem.hh>


namespace Ewoms {
namespace Properties {

    // NEW_TYPE_TAG(TestProblem);
    // NEW_TYPE_TAG(DerivedTestProblem, INHERITS_FROM(TestProblem));

    NEW_PROP_TAG(Grid);
    NEW_PROP_TAG(Simulator);
    NEW_PROP_TAG(FluidSystem);
    NEW_PROP_TAG(Indices);
    NEW_PROP_TAG(IntensiveQuantities);
    NEW_PROP_TAG(MaterialLaw);

}
}


#include <opm/autodiff/WellInterface.hpp>

/*
#include <vector>
#include <unordered_set>
#include <memory>
#include <array>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>

#include <opm/core/grid.h>
#include <opm/core/props/satfunc/SaturationPropsFromDeck.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/wells/DynamicListEconLimited.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/createGlobalCellArray.hpp>
#include <opm/autodiff/GridInit.hpp>

#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/MultisegmentWells.hpp>
#include <opm/autodiff/WellInterface.hpp>
//#include <opm/autodiff/StandardWell.hpp>


// NEW_TYPE_TAG(TestProblem);


struct SetupStandardWell {

    using Grid = UnstructuredGrid;
    using GridInit = Opm::GridInit<Grid>;
    using FluidProps = Opm::BlackoilPropsAdFromDeck;
    using MaterialLawManager = FluidProps::MaterialLawManager;

    SetupStandardWell()
    {
        Opm::ParseContext parse_context;
        Opm::Parser parser;
        auto deck = parser.parseFile("msw.data", parse_context);
        Opm::EclipseState ecl_state(deck , parse_context);

        // Create grid.
        const std::vector<double>& porv =
                            ecl_state.get3DProperties().getDoubleGridProperty("PORV").getData();

        std::unique_ptr<GridInit> grid_init(new GridInit(ecl_state, porv));
        const Grid& grid = grid_init->grid();

        const size_t current_timestep = 0;

        // dummy_dynamic_list_econ_lmited
        const Opm::DynamicListEconLimited dummy_dynamic_list;

        // Create wells.
        Opm::WellsManager wells_manager(ecl_state,
                                        current_timestep,
                                        Opm::UgGridHelpers::numCells(grid),
                                        Opm::UgGridHelpers::globalCell(grid),
                                        Opm::UgGridHelpers::cartDims(grid),
                                        Opm::UgGridHelpers::dimensions(grid),
                                        Opm::UgGridHelpers::cell2Faces(grid),
                                        Opm::UgGridHelpers::beginFaceCentroids(grid),
                                        dummy_dynamic_list,
                                        false,
                                        // We need to pass the optionaly arguments
                                        // as we get the following error otherwise
                                        // with c++ (Debian 4.9.2-10) 4.9.2 and -std=c++11
                                        // converting to ‘const std::unordered_set<std::basic_string<char> >’ from initializer list would use explicit constructor
                                        std::unordered_set<std::string>());

        const Wells* wells = wells_manager.c_wells();
        const auto& wells_ecl = ecl_state.getSchedule().getWells(current_timestep);

        // ms_wells.reset(new Opm::MultisegmentWells(wells, &(wells_manager.wellCollection()), wells_ecl, current_timestep));
    };

    // StandardWell<TestProblem> standard_well;
};

// number of wells for this case
const int nw = 2;
// number of segments for this case
const int nseg = 7;
// number of perforations for this case
const int nperf = 8;

BOOST_AUTO_TEST_CASE(testConstruction)
{
    SetupStandardWell setup;
}
*/

BOOST_AUTO_TEST_CASE(testConstruction)
{
    int x;
}
