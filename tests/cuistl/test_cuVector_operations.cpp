/*
  Copyright 2022-2023 SINTEF AS

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

#define BOOST_TEST_MODULE TestCuVectorOperations
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <dune/istl/bcrsmatrix.hh>
#include <opm/simulators/linalg/cuistl/CuJac.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/cuistl/detail/vector_operations.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>

using NumericTypes = boost::mpl::list<double, float>;

BOOST_AUTO_TEST_CASE_TEMPLATE(ElementWiseMultiplicationOf3By3BlockVectorAndVectorVector, T, NumericTypes)
{
    /*
        Example in the test for multiplying by element a blockvector with a vector of vectors
        | |1 2 3| |   | |3| |   | |10| |
        | |5 2 3| | X | |2| | = | |22| |
        | |2 1 1| |   | |1| |   | |10| |
    */

    const size_t blocksize = 3;
    const size_t N = 1;

    std::vector<T> h_blockVector({1.0,2.0,3.0,5.0,2.0,3.0,2.0,1.0,2.0});
    std::vector<T> h_vecVector({3.0,2.0,1.0});
    Opm::cuistl::CuVector<T> d_blockVector(h_blockVector);
    Opm::cuistl::CuVector<T> d_vecVector(h_vecVector);

    Opm::cuistl::detail::blockVectorMultiplicationAtAllIndices(d_blockVector.data(), N, blocksize, d_vecVector.data());

    std::vector<T> expected_vec{10.0,22.0,10.0};
    std::vector<T> computed_vec = d_vecVector.asStdVector();

    BOOST_REQUIRE_EQUAL(expected_vec.size(), computed_vec.size());
    for (size_t i = 0; i < expected_vec.size(); i++){
        BOOST_CHECK_CLOSE(expected_vec[i], computed_vec[i], 1e-7);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ElementWiseMultiplicationOf2By2BlockVectorAndVectorVector, T, NumericTypes)
{
    /*
        Example in the test for multiplying by element a blockvector with a vector of vectors
        | |1 2| |   | |1| |   | | 7| |
        | |3 4| | X | |3| | = | |15| |
        |       |   |     |   |      |
        | |4 3| |   | |2| |   | |20| |
        | |2 1| |   | |4| |   | | 8| |
    */

    const size_t blocksize = 2;
    const size_t N = 2;

    std::vector<T> h_blockVector({1.0,2.0,3.0,4.0,4.0,3.0,2.0,1.0});
    std::vector<T> h_vecVector({1.0,3.0,2.0,4.0});
    Opm::cuistl::CuVector<T> d_blockVector(h_blockVector);
    Opm::cuistl::CuVector<T> d_vecVector(h_vecVector);

    Opm::cuistl::detail::blockVectorMultiplicationAtAllIndices(d_blockVector.data(), N, blocksize, d_vecVector.data());

    std::vector<T> expected_vec{7.0,15.0,20.0,8.0};
    std::vector<T> computed_vec = d_vecVector.asStdVector();

    BOOST_REQUIRE_EQUAL(expected_vec.size(), computed_vec.size());
    for (size_t i = 0; i < expected_vec.size(); i++){
        BOOST_CHECK_CLOSE(expected_vec[i], computed_vec[i], 1e-7);
    }
}