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


#ifndef OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED
#define OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED

#include <opm/simulators/linalg/makePreconditioner.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>

#include <boost/property_tree/ptree.hpp>

namespace Dune
{



template <class MatrixTypeT, class VectorTypeT, class Communication>
class FlexibleSolver : Dune::InverseOperator<VectorTypeT, VectorTypeT>
{
public:
    using MatrixType = MatrixTypeT;
    using VectorType = VectorTypeT;

    FlexibleSolver(const boost::property_tree::ptree& prm, const MatrixType& matrix)
    {
        init(prm, matrix, Communication());
    }

    FlexibleSolver(const boost::property_tree::ptree& prm, const MatrixType& matrix, const Communication& comm)
    {
        init(prm, matrix, comm);
    }

    virtual void apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res) override
    {
        linsolver_->apply(x, rhs, res);
    }

    virtual void apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res) override
    {
        linsolver_->apply(x, rhs, reduction, res);
    }

    void updatePreconditioner()
    {
        preconditioner_->update();
    }

    virtual Dune::SolverCategory::Category category() const override
    {
        return IsSequential ? Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
    }

private:

    using SeqOperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
    using ParOperatorType = Dune::OverlappingSchwarzOperator<MatrixType, VectorType, VectorType, Communication>;
    static constexpr bool IsSequential = std::is_same<Communication, Dune::Amg::SequentialInformation>::value;
    using OperatorType = std::conditional_t<IsSequential, SeqOperatorType, ParOperatorType>;

    // Machinery for making sequential or parallel operator.
    template <class Comm>
    std::shared_ptr<OperatorType> makeOperator(const MatrixType& matrix, const Comm& comm)
    {
        return std::make_shared<OperatorType>(matrix, comm); // Parallel case.
    }
    template <>
    std::shared_ptr<OperatorType> makeOperator<Dune::Amg::SequentialInformation>(const MatrixType& matrix, const Dune::Amg::SequentialInformation&)
    {
        return std::make_shared<OperatorType>(matrix); // Sequential case.
    }

    // Machinery for making sequential or parallel preconditioner.
    template <class Comm>
    std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
    makePrecond(OperatorType& op,
                const boost::property_tree::ptree& prm,
                const Comm& comm)
    {
        return Dune::makePreconditioner<MatrixType, VectorType>(op, prm, comm); // Parallel case.
    }
    template <>
    std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
    makePrecond<Dune::Amg::SequentialInformation>(OperatorType& op,
                                                  const boost::property_tree::ptree& prm,
                                                  const Dune::Amg::SequentialInformation&)
    {
        return Dune::makePreconditioner<MatrixType, VectorType>(op, prm); // Sequential case.
    }

    void init(const boost::property_tree::ptree& prm, const MatrixType& matrix, const Communication& comm)
    {
        linearoperator_ = makeOperator(matrix, comm);
        preconditioner_ = makePrecond(*linearoperator_, prm, comm);
        scalarproduct_ = createScalarProduct<VectorType, Communication>(comm, category());
        const double tol = prm.get<double>("tol");
        const int maxiter = prm.get<int>("maxiter");
        const int verbosity = prm.get<int>("verbosity");
        const std::string solver_type = prm.get<std::string>("solver");
        if (solver_type == "bicgstab") {
            linsolver_.reset(new Dune::BiCGSTABSolver<VectorType>(*linearoperator_,
                                                                  *scalarproduct_,
                                                                  *preconditioner_,
                                                                  tol, // desired residual reduction factor
                                                                  maxiter, // maximum number of iterations
                                                                  verbosity));
        } else if (solver_type == "loopsolver") {
            linsolver_.reset(new Dune::LoopSolver<VectorType>(*linearoperator_,
                                                              *scalarproduct_,
                                                              *preconditioner_,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity));
        } else if (solver_type == "gmres") {
            int restart = prm.get<int>("restart");
            linsolver_.reset(new Dune::RestartedGMResSolver<VectorType>(*linearoperator_,
                                                                        *scalarproduct_,
                                                                        *preconditioner_,
                                                                        tol,
                                                                        restart, // desired residual reduction factor
                                                                        maxiter, // maximum number of iterations
                                                                        verbosity));
#if HAVE_SUITESPARSE_UMFPACK
        } else if (solver_type == "umfpack") {
            bool dummy = false;
            linsolver_.reset(new Dune::UMFPack<MatrixType>(linearoperator_->getmat(), verbosity, dummy));
#endif
        } else {
            std::string msg("Solver not known ");
            msg += solver_type;
            throw std::runtime_error(msg);
        }
    }

    std::shared_ptr<OperatorType> linearoperator_;
    std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>> preconditioner_;
    std::shared_ptr<Dune::ScalarProduct<VectorType>> scalarproduct_;
    std::shared_ptr<Dune::InverseOperator<VectorType, VectorType>> linsolver_;
};

} // namespace Dune



#endif // OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED
