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

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/matrixredistribute.hh>
#include <dune/istl/paamg/graph.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solvers.hh>

int
main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    auto world_comm = Dune::MPIHelper::getCollectiveCommunication();
    const int BS = 3; // block size, sparse scalar matrix
    typedef Dune::FieldMatrix<double, BS, BS> MatrixBlock; // matrix block type
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat; // sparse matrix type
    typedef Dune::FieldVector<double, BS> VectorBlock; // vector block type
    typedef Dune::BlockVector<VectorBlock> BVector; // vector type
    BCRSMat A;

    // read matrix on rank 0
    if (world_comm.rank() == 0) {
        loadMatrixMarket(A, std::string("sherman5.mtx"));
    }

    typedef std::size_t GlobalId; // The type for the global index
    typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;

    Communication comm(world_comm);
    // No need to add any indices to comm.indexSet()

    // Now we can create a parallel representation of the matrix and use PT-Scotch or ParMETIS (ParMETIS is non-free for
    // non-academic use!) to distribute the sequential matrix to all processes. (See this post on how to install and use
    // PT-Scotch.)

    Communication* comm_redist;
    BCRSMat parallel_A;
    typedef Dune::Amg::MatrixGraph<BCRSMat> MatrixGraph;
    Dune::RedistributeInformation<Communication> rinfo;
    bool hasDofs = Dune::graphRepartition(
        MatrixGraph(A), comm, static_cast<int>(world_comm.size()), comm_redist, rinfo.getInterface(), true);
    rinfo.setSetup();
    redistributeMatrix(A, parallel_A, comm, *comm_redist, rinfo);
    // Note that you have to make sure that DUNE detected either ParMETIS or PT-Scotch on your system. Otherwise the
    // above code will not work. As the first step we created the empty parallel matrix parallel_A and the data
    // structure rinfo. The latter encapsulates the information how to redistribute the matrix. Then we use the graph of
    // the matrix to load balance it and setup the communication interface. This happens in the call to function
    // graphRepartition. In the last step we redistribute the matrix according to the information of the load balancer.

    // Of course we might also want to distribute a vector repesenting the right hand side of the linear system. We can
    // do that using the data structure rinfo:

    BVector b(A.N());
    BVector parallel_b(parallel_A.N());
    BVector parallel_x(parallel_A.M());
    b = 100.0;
    rinfo.redistribute(b, parallel_b);
    // Now we have the complete linear system and can solve it in parallel. To do this we need to set up the parallel
    // solver components and call the apply method of the solver. Below we use the conjugate gradient method
    // preconditioned with hybrid SSOR:

    if (hasDofs) // if hasDofs is false we do not compute.
    {
        // the index set has changed. Rebuild the remote information
        comm_redist->remoteIndices().rebuild<false>();
        typedef Dune::SeqSSOR<BCRSMat, BVector, BVector> Prec;
        typedef Dune::BlockPreconditioner<BVector, BVector, Communication, Prec>
            ParPrec; // type of parallel preconditioner
        typedef Dune::OverlappingSchwarzScalarProduct<BVector, Communication>
            ScalarProduct; // type of parallel scalar product
        typedef Dune::OverlappingSchwarzOperator<BCRSMat, BVector, BVector, Communication>
            Operator; // type of parallel linear operator

        ScalarProduct sp(*comm_redist);
        Operator op(parallel_A, *comm_redist);
        Prec prec(parallel_A, 1, 1.0);
        ParPrec pprec(prec, *comm_redist);
        Dune::InverseOperatorResult r;
        Dune::BiCGSTABSolver<BVector> cg(op, sp, pprec, 10e-8, 80, world_comm.rank() == 0 ? 2 : 0);
        cg.apply(parallel_x, parallel_b, r);
    }
    // If you really need to then you can also gather all the information on the master process afterwards:

    BVector x(A.M());
    rinfo.redistributeBackward(x, parallel_x);
}
