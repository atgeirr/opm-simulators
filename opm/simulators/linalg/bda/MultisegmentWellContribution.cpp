/*
  Copyright 2020 Equinor ASA

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


#include <cstdlib>
#include <cstring>
#include <config.h> // CMake
#include <fstream>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#include <opm/simulators/linalg/bda/cuda_header.hpp>
#endif // HAVE_UMFPACK

#include <opm/simulators/linalg/bda/MultisegmentWellContribution.hpp>

namespace Opm
{

        MultisegmentWellContribution::MultisegmentWellContribution(unsigned int dim_, unsigned int dim_wells_,
            unsigned int Nb_, unsigned int Mb_,
            unsigned int BnumBlocks_, double *Bvalues, unsigned int *BcolIndices, unsigned int *BrowPointers,
            unsigned int DnumBlocks_, double *Dvalues, int *DcolPointers, int *DrowIndices,
            double *Cvalues)
        :
            dim(dim_),
            dim_wells(dim_wells_),
            N(Nb_*dim),
            Nb(Nb_),
            M(Mb_*dim_wells),
            Mb(Mb_),
            DnumBlocks(DnumBlocks_),
            BnumBlocks(BnumBlocks_)
        {
            Cvals.insert(Cvals.end(), Cvalues, Cvalues + BnumBlocks*dim*dim_wells);
            Dvals.insert(Dvals.end(), Dvalues, Dvalues + DnumBlocks*dim_wells*dim_wells);
            Bvals.insert(Bvals.end(), Bvalues, Bvalues + BnumBlocks*dim*dim_wells);
            Bcols.insert(Bcols.end(), BcolIndices, BcolIndices + BnumBlocks);
            Brows.insert(Brows.end(), BrowPointers, BrowPointers + M + 1);

            Dcols.insert(Dcols.end(), DcolPointers, DcolPointers + M + 1);
            Drows.insert(Drows.end(), DrowIndices, DrowIndices + DnumBlocks*dim_wells*dim_wells);

            z1.resize(Mb * dim_wells);
            z2.resize(Mb * dim_wells);

            // allocate pinned memory on host
            cudaMallocHost(&h_x, sizeof(double) * N);
            cudaMallocHost(&h_y, sizeof(double) * N);

            umfpack_di_symbolic(M, M, Dcols.data(), Drows.data(), Dvals.data(), &UMFPACK_Symbolic, nullptr, nullptr);
            umfpack_di_numeric(Dcols.data(), Drows.data(), Dvals.data(), UMFPACK_Symbolic, &UMFPACK_Numeric, nullptr, nullptr);
        }

        MultisegmentWellContribution::~MultisegmentWellContribution()
        {
            cudaFreeHost(h_x);
            cudaFreeHost(h_y);

            umfpack_di_free_symbolic(&UMFPACK_Symbolic);
            umfpack_di_free_numeric(&UMFPACK_Numeric);
        }


        // Apply the MultisegmentWellContribution, similar to MultisegmentWell::apply()
        // y -= (C^T * (D^-1 * (B * x)))
        void MultisegmentWellContribution::apply(double *d_x, double *d_y)
        {
            // copy vectors x and y from GPU to CPU
            cudaMemcpyAsync(h_x, d_x, sizeof(double) * N, cudaMemcpyDeviceToHost, stream);
            cudaMemcpyAsync(h_y, d_y, sizeof(double) * N, cudaMemcpyDeviceToHost, stream);
            cudaStreamSynchronize(stream);

            // reset z1 and z2
            std::fill(z1.begin(), z1.end(), 0.0);
            std::fill(z2.begin(), z2.end(), 0.0);

            // z1 = B * x
            for (unsigned int row = 0; row < Mb; ++row) {
                // for every block in the row
                for (unsigned int blockID = Brows[row]; blockID < Brows[row+1]; ++blockID) {
                    unsigned int colIdx = Bcols[blockID];
                    for (unsigned int j = 0; j < dim_wells; ++j) {
                        double temp = 0.0;
                        for (unsigned int k = 0; k < dim; ++k) {
                            temp += Bvals[blockID * dim * dim_wells + j * dim + k] * h_x[colIdx * dim + k];
                        }
                        z1[row * dim_wells + j] += temp;
                    }
                }
            }

            // z2 = D^-1 * (B * x)
            // umfpack
            umfpack_di_solve(UMFPACK_A, Dcols.data(), Drows.data(), Dvals.data(), z2.data(), z1.data(), UMFPACK_Numeric, nullptr, nullptr);

            // y -= (C^T * z2)
            // y -= (C^T * (D^-1 * (B * x)))
            for (unsigned int row = 0; row < Mb; ++row) {
                // for every block in the row
                for (unsigned int blockID = Brows[row]; blockID < Brows[row+1]; ++blockID) {
                    unsigned int colIdx = Bcols[blockID];
                    for (unsigned int j = 0; j < dim; ++j) {
                        double temp = 0.0;
                        for (unsigned int k = 0; k < dim_wells; ++k) {
                            temp += Cvals[blockID * dim * dim_wells + j + k * dim] * z2[row * dim_wells + k];
                        }
                        h_y[colIdx * dim + j] -= temp;
                    }
                }
            }

            // copy vector y from CPU to GPU
            cudaMemcpyAsync(d_y, h_y, sizeof(double) * N, cudaMemcpyHostToDevice, stream);
            cudaStreamSynchronize(stream);
        }

        void MultisegmentWellContribution::setCudaStream(cudaStream_t stream_)
        {
            stream = stream_;
        }


} //namespace Opm

