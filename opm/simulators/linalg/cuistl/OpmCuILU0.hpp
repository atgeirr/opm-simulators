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
#ifndef OPM_CUILU0_OPM_Impl_HPP
#define OPM_CUILU0_OPM_Impl_HPP

#include <memory>
#include <opm/grid/utility/SparseTable.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <optional>
#include <type_traits>
#include <vector>


namespace Opm::cuistl
{
//! \brief ILU0 preconditioner on the GPU.
//!
//! \tparam M The matrix type to operate on
//! \tparam X Type of the update
//! \tparam Y Type of the defect
//! \tparam l Ignored. Just there to have the same number of template arguments
//!    as other preconditioners.
//!
//! \note We assume X and Y are both CuVector<real_type>, but we leave them as template
//! arguments in case of future additions.
template <class M, class X, class Y, int l = 1>
class OpmCuILU0 : public Dune::PreconditionerWithUpdate<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::remove_const<M>::type;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief The GPU matrix type
    using CuMat = CuSparseMatrix<field_type>;

    //! \brief Constructor.
    //!
    //!  Constructor gets all parameters to operate the prec.
    //! \param A The matrix to operate on.
    //! \param w The relaxation factor.
    //!
    explicit OpmCuILU0(const M& A, bool splitMatrix, bool tuneKernels);

    //! \brief Prepare the preconditioner.
    //! \note Does nothing at the time being.
    void pre(X& x, Y& b) override;

    //! \brief Apply the preconditoner.
    void apply(X& v, const Y& d) override;

    //! \brief Post processing
    //! \note Does nothing at the moment
    void post(X& x) override;

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override;

    //! \brief Updates the matrix data.
    void update() final;

    //! \brief Compute LU factorization, and update the data of the reordered matrix
    void LUFactorizeAndMoveData();

    //! \brief function that will experimentally tune the thread block sizes of the important cuda kernels
    void tuneThreadBlockSizes();

    //! \returns false
    static constexpr bool shouldCallPre()
    {
        return false;
    }

    //! \returns false
    static constexpr bool shouldCallPost()
    {
        return false;
    }

    virtual bool hasPerfectUpdate() const override {
        return true;
    }


private:
    //! \brief Reference to the underlying matrix
    const M& m_cpuMatrix;
    //! \brief size_t describing the dimensions of the square block elements
    static constexpr const size_t blocksize_ = matrix_type::block_type::cols;
    //! \brief SparseTable storing each row by level
    Opm::SparseTable<size_t> m_levelSets;
    //! \brief converts from index in reordered structure to index natural ordered structure
    std::vector<int> m_reorderedToNatural;
    //! \brief converts from index in natural ordered structure to index reordered strucutre
    std::vector<int> m_naturalToReordered;
    //! \brief The A matrix stored on the gpu, and its reordred version
    CuMat m_gpuMatrix;
    std::unique_ptr<CuMat> m_gpuReorderedLU;
    //! \brief If matrix splitting is enabled, then we store the lower and upper part separately
    std::unique_ptr<CuMat> m_gpuMatrixReorderedLower;
    std::unique_ptr<CuMat> m_gpuMatrixReorderedUpper;
    //! \brief If matrix splitting is enabled, we also store the diagonal separately
    std::optional<CuVector<field_type>> m_gpuMatrixReorderedDiag;
    //! row conversion from natural to reordered matrix indices stored on the GPU
    CuVector<int> m_gpuNaturalToReorder;
    //! row conversion from reordered to natural matrix indices stored on the GPU
    CuVector<int> m_gpuReorderToNatural;
    //! \brief Stores the inverted diagonal that we use in ILU0
    CuVector<field_type> m_gpuDInv;
    //! \brief Bool storing whether or not we should store matrices in a split format
    bool m_splitMatrix;
    //! \brief Bool storing whether or not we will tune the threadblock sizes. Only used for AMD cards
    bool m_tuneThreadBlockSizes;
    //! \brief variables storing the threadblocksizes to use if using the tuned sizes and AMD cards
    //! The default value of -1 indicates that we have not calibrated and selected a value yet
    int m_applyThreadBlockSize = -1;
    int m_updateThreadBlockSize = -1;
};
} // end namespace Opm::cuistl

#endif