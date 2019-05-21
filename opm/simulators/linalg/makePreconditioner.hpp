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

#ifndef OPM_MAKEPRECONDITIONER_HEADER_INCLUDED
#define OPM_MAKEPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/amgcpr.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/preconditioners.hh>

#include <boost/property_tree/ptree.hpp>

namespace Dune
{

template <class MatrixType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeSeqPreconditioner(const Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                      const boost::property_tree::ptree& prm)
{
    auto& matrix = linearoperator.getmat();
    double w = prm.get<double>("w");
    int n = prm.get<int>("n");
    std::string precond(prm.get<std::string>("preconditioner"));
    if (precond == "ILU0") {
        return wrapPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>>(matrix, w);
    } else if (precond == "ParOverILU0") {
        return wrapPreconditioner<Opm::ParallelOverlappingILU0<MatrixType, VectorType, VectorType>>(matrix, n, w, Opm::MILU_VARIANT::ILU);
    } else if (precond == "Jac") {
        return wrapPreconditioner<Dune::SeqJac<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else if (precond == "GS") {
        return wrapPreconditioner<Dune::SeqGS<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else if (precond == "SOR") {
        return wrapPreconditioner<Dune::SeqSOR<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else if (precond == "SSOR") {
        return wrapPreconditioner<Dune::SeqSSOR<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else if (precond == "ILUn") {
        return wrapPreconditioner<Dune::SeqILUn<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else {
        std::string msg("No such preconditioner : ");
        msg += precond;
        throw std::runtime_error(msg);
    }
}

template <class Smoother, class MatrixType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeAmgPreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                      const boost::property_tree::ptree& global_prm)
{
    boost::property_tree::ptree prm = global_prm.get_child("amg");
    using OperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
    using CriterionBase = Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<MatrixType, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;
    int coarsenTarget = prm.get<int>("coarsenTarget");
    int ml = prm.get<int>("maxlevel");
    Criterion criterion(15, coarsenTarget);
    criterion.setDefaultValuesIsotropic(2);
    criterion.setAlpha(prm.get<double>("alpha"));
    criterion.setBeta(prm.get<double>("beta"));
    criterion.setMaxLevel(ml);
    criterion.setSkipIsolated(false);
    criterion.setDebugLevel(prm.get<int>("verbosity"));
    if (global_prm.get<std::string>("preconditioner") == "famg") {
        Dune::Amg::Parameters parms;
        parms.setNoPreSmoothSteps(1);
        parms.setNoPostSmoothSteps(1);
        return wrapPreconditioner<Dune::Amg::FastAMG<OperatorType, VectorType>>(linearoperator, criterion, parms);
    } else {
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = prm.get<int>("n");
        // smootherArgs.overlap=SmootherArgs::vertex;
        // smootherArgs.overlap=SmootherArgs::none;
        // smootherArgs.overlap=SmootherArgs::aggregate;
        smootherArgs.relaxationFactor = prm.get<double>("w");
        return wrapPreconditioner<Dune::Amg::AMGCPR<OperatorType, VectorType, Smoother>>(linearoperator, criterion, smootherArgs);
    }
}




template <class MatrixType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeAmgPreconditioners(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                       const boost::property_tree::ptree& prm)
{
    if (prm.get<std::string>("preconditioner") == "famg") {
        // smoother type should not be used
        return makeAmgPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    }

    std::string precond = prm.get<std::string>("amg.smoother");
    if (precond == "ILU0") {
        return makeAmgPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    } else if (precond == "Jac") {
        return makeAmgPreconditioner<Dune::SeqJac<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
        // } else if (precond == "GS") {
        //   return makeAmgPreconditioner<
        // 	Dune::SeqGS<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(linearoperator, prm);
    } else if (precond == "SOR") {
        return makeAmgPreconditioner<Dune::SeqSOR<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    } else if (precond == "SSOR") {
        return makeAmgPreconditioner<Dune::SeqSSOR<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    } else if (precond == "ILUn") {
        return makeAmgPreconditioner<Dune::SeqILUn<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    } else {
        std::string msg("No such seq preconditioner : ");
        msg += precond;
        throw std::runtime_error(msg);
    }
}

template <class MatrixType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeTwoLevelPreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                           const boost::property_tree::ptree& global_prm)
{
    boost::property_tree::ptree prm = global_prm.get_child("cpr");
    if (global_prm.get<std::string>("preconditioner") == "cpr") {
        return std::make_shared<OwningTwoLevelPreconditioner<MatrixType, VectorType, false>>(linearoperator, prm);
    } else if (global_prm.get<std::string>("preconditioner") == "cprt") {
        return std::make_shared<OwningTwoLevelPreconditioner<MatrixType, VectorType, true>>(linearoperator, prm);
    } else {
        std::string msg("Wrong cpr Should not happen");
        throw std::runtime_error(msg);
    }
}

template <class MatrixType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makePreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                   const boost::property_tree::ptree& prm)
{
    if ((prm.get<std::string>("preconditioner") == "famg") or (prm.get<std::string>("preconditioner") == "amg")) {
        return makeAmgPreconditioners<MatrixType, VectorType>(linearoperator, prm);
    } else if ((prm.get<std::string>("preconditioner") == "cpr")
               or (prm.get<std::string>("preconditioner") == "cprt")) {
        return makeTwoLevelPreconditioner<MatrixType, VectorType>(linearoperator, prm);
    } else {
        return makeSeqPreconditioner<MatrixType, VectorType>(linearoperator, prm);
    }
}


} // namespace Dune


#endif // OPM_MAKEPRECONDITIONER_HEADER_INCLUDED
