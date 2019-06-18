/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED
#define OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED

#include <ewoms/linear/matrixblock.hh>
#include <opm/simulators/linalg/findOverlapRowsAndColumns.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <memory>
#include <utility>

BEGIN_PROPERTIES

NEW_TYPE_TAG(FlowIstlSolverFlexible, INHERITS_FROM(FlowIstlSolverParams));

NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(Simulator);

END_PROPERTIES


namespace Opm
{

//=====================================================================
// Implementation for ISTL-matrix based operator
//=====================================================================
/// This class solves the fully implicit black-oil system by
/// solving the reduced system (after eliminating well variables)
/// as a block-structured matrix (one block for all cell variables) for a fixed
/// number of cell variables.
///
/// The solvers and preconditioners used are run-time configurable.
template <class TypeTag>
class ISTLSolverEbosFlexible
{
    using SparseMatrixAdapter = typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter);
    using VectorType = typename GET_PROP_TYPE(TypeTag, GlobalEqVector);
    using Simulator = typename GET_PROP_TYPE(TypeTag, Simulator);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MatrixType = typename SparseMatrixAdapter::IstlMatrix;
#if HAVE_MPI
    using Communication = Dune::OwnerOverlapCopyCommunication<int, int>;
#endif
    using SolverType = Dune::FlexibleSolver<MatrixType, VectorType>;
    using FloatMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<float, 3, 3>>;
    using FloatVectorType = Dune::BlockVector<Dune::FieldVector<float, 3>>;
    using FloatSolverType = Dune::FlexibleSolver<FloatMatrixType, FloatVectorType>;


public:
    static void registerParameters()
    {
        FlowLinearSolverParameters::registerParameters<TypeTag>();
    }

    explicit ISTLSolverEbosFlexible(const Simulator& simulator)
        : simulator_(simulator)
    {
        parameters_.template init<TypeTag>();
        prm_ = setupPropertyTree(parameters_);
        extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);
        detail::findOverlapRowsAndColumns(simulator_.vanguard().grid(), overlapRowAndColumns_);
#if HAVE_MPI
        if (parallelInformation_.type() == typeid(ParallelISTLInformation)) {
            // Parallel case.
            const ParallelISTLInformation* parinfo = boost::any_cast<ParallelISTLInformation>(&parallelInformation_);
            assert(parinfo);
            comm_.reset(new Communication(parinfo->communicator()));
        }
#endif
    }

    void eraseMatrix()
    {
    }

    void prepare(SparseMatrixAdapter& mat, VectorType& b)
    {
#if HAVE_MPI
        static bool firstcall = true;
        if (firstcall && parallelInformation_.type() == typeid(ParallelISTLInformation)) {
            // Parallel case.
            const ParallelISTLInformation* parinfo = boost::any_cast<ParallelISTLInformation>(&parallelInformation_);
            assert(parinfo);
            const size_t size = mat.istlMatrix().N();
            parinfo->copyValuesTo(comm_->indexSet(), comm_->remoteIndices(), size, 1);
            firstcall = false;
        }
        makeOverlapRowsInvalid(mat.istlMatrix());
#endif
        if (floatA_.N() == 0) {
            buildFloatMatrixStructure(mat.istlMatrix());
        }

        // Decide if we should recreate the solver or just do
        // a minimal preconditioner update.
        const int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
        bool recreate_solver = false;
        if (this->parameters_.cpr_reuse_setup_ == 0) {
            // Always recreate solver.
            recreate_solver = true;
        } else if (this->parameters_.cpr_reuse_setup_ == 1) {
            // Recreate solver on the first iteration of every timestep.
            if (newton_iteration == 0) {
                recreate_solver = true;
            }
        } else if (this->parameters_.cpr_reuse_setup_ == 2) {
            // Recreate solver if the last solve used more than 10 iterations.
            if (this->iterations() > 10) {
                recreate_solver = true;
            }
        } else {
            assert(this->parameters_.cpr_reuse_setup_ == 3);
            assert(recreate_solver == false);
            // Never recreate solver.
        }

        const bool useFloat = true;
        if (useFloat) {
            copyValuesToFloatMatrix(mat.istlMatrix());
            if (recreate_solver || !floatSolver_) {
                if (isParallel()) {
#if HAVE_MPI
                    floatSolver_.reset(new FloatSolverType(prm_, floatA_, *comm_));
#endif
                } else {
                    floatSolver_.reset(new FloatSolverType(prm_, floatA_));
                }
            } else {
                floatSolver_->preconditioner().update();
            }
            const size_t n = b.size();
            floatRhs_.resize(n);
            for (size_t i = 0; i < n; ++i) {
                floatRhs_[i] = b[i];
            }
            floatX_.resize(n);
        } else {
            if (recreate_solver || !solver_) {
                if (isParallel()) {
#if HAVE_MPI
                    solver_.reset(new SolverType(prm_, mat.istlMatrix(), *comm_));
#endif
                } else {
                    solver_.reset(new SolverType(prm_, mat.istlMatrix()));
                }
                rhs_ = b;
            } else {
                solver_->preconditioner().update();
                rhs_ = b;
            }
        }
    }

    bool solve(VectorType& x)
    {
        const bool useFloat = true;
        if (useFloat) {
            const size_t n = floatRhs_.size();
            floatX_.resize(n);
            for (size_t i = 0; i < n; ++i) {
                floatX_[i] = x[i];
            }
            floatSolver_->apply(floatX_, floatRhs_, res_);
            for (size_t i = 0; i < n; ++i) {
                x[i] = floatX_[i];
            }
        } else {
            solver_->apply(x, rhs_, res_);
        }
        return res_.converged;
    }

    bool isParallel() const
    {
#if HAVE_MPI
        return parallelInformation_.type() == typeid(ParallelISTLInformation);
#else
        return false;
#endif
    }

    int iterations() const
    {
        return res_.iterations;
    }

    void setResidual(VectorType& /* b */)
    {
        // rhs_ = &b; // Must be handled in prepare() instead.
    }

    void setMatrix(const SparseMatrixAdapter& /* M */)
    {
        // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
    }

protected:
    /// Zero out off-diagonal blocks on rows corresponding to overlap cells
    /// Diagonal blocks on ovelap rows are set to diag(1e100).
    void makeOverlapRowsInvalid(MatrixType& matrix) const
    {
        // Value to set on diagonal
        const int numEq = MatrixType::block_type::rows;
        typename MatrixType::block_type diag_block(0.0);
        for (int eq = 0; eq < numEq; ++eq)
            diag_block[eq][eq] = 1.0e100;

        // loop over precalculated overlap rows and columns
        for (auto row = overlapRowAndColumns_.begin(); row != overlapRowAndColumns_.end(); row++) {
            int lcell = row->first;
            // diagonal block set to large value diagonal
            matrix[lcell][lcell] = diag_block;

            // loop over off diagonal blocks in overlap row
            for (auto col = row->second.begin(); col != row->second.end(); ++col) {
                int ncell = *col;
                // zero out block
                matrix[lcell][ncell] = 0.0;
            }
        }
    }

    template <class M>
    void buildFloatMatrixStructure(const M& m)
    {
        // Create ISTL matrix with interleaved rows and columns (block structured).
        floatA_.setSize(m.N(), m.M(), m.nonzeroes());
        floatA_.setBuildMode(FloatMatrixType::row_wise);
        auto frow = floatA_.createbegin();
        auto row = m.begin();
        auto endrow = m.end();
        for (; row != endrow; ++row, ++frow) {
            const int ri = row.index();
            auto col = row->begin();
            auto endcol = row->end();
            for (; col != endcol; ++col) {
                frow.insert(col.index());
            }
        }
    }

    template <class M>
    void copyValuesToFloatMatrix(const M& m)
    {
        auto frow = floatA_.begin();
        auto row = m.begin();
        auto endrow = m.end();
        for (; row != endrow; ++row, ++frow) {
            auto fcol = frow->begin();
            auto col = row->begin();
            auto endcol = row->end();
            for (; col != endcol; ++col, ++fcol) {
                *fcol = *col;
            }
        }
    }

    const Simulator& simulator_;

    std::unique_ptr<SolverType> solver_;
    std::unique_ptr<FloatSolverType> floatSolver_;
    FlowLinearSolverParameters parameters_;
    boost::property_tree::ptree prm_;
    VectorType rhs_;
    FloatMatrixType floatA_;
    FloatVectorType floatRhs_;
    FloatVectorType floatX_;
    Dune::InverseOperatorResult res_;
    boost::any parallelInformation_;
#if HAVE_MPI
    std::unique_ptr<Communication> comm_;
#endif
    std::vector<std::pair<int, std::vector<int>>> overlapRowAndColumns_;
}; // end ISTLSolverEbosFlexible

} // namespace Opm

#endif // OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED
