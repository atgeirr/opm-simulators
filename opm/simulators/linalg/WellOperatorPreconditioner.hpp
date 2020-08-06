/*
  Copyright 2020 SINTEF Digital, Mathematics and Cybernetics.

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


#ifndef OPM_WELL_OPERATOR_PRECONDITIONER_HEADER_INCLUDED
#define OPM_WELL_OPERATOR_PRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

namespace Opm
{

template <class ReservoirPreconditioner>
class WellOperatorPreconditioner : public PreconditionerWithUpdate<typename ReservoirPreconditioner::domain_type,
                                                                   typename ReservoirPreconditioner::range_type>
{
public:
    template <class... Args>
    WellOperatorPreconditioner(const Dune::LinearOperator<X, Y>& well_op,
                               Args&&... args)
        : well_op_(well_op)
        , res_precond_(std::forward<Args>(args)...)
    {
    }

    using X = typename ReservoirPreconditioner::domain_type;
    using Y = typename ReservoirPreconditioner::range_type;

    virtual void pre(X& x, Y& b) override
    {
        // TODO effect of wells here?
        res_precond_.pre(x, b);
    }

    virtual void apply(X& v, const Y& d) override
    {
        X temp = v;
        well_op_.apply(d, temp);     // temp = (I - CD^{-1}B)d
        res_precond_.apply(v, temp); // ~= A^{-1}(I -CD^{-1}B)d
        // (because applying this preconditioner is approximating A^{-1})
    }

    virtual void post(X& x) override
    {
        // TODO effect of wells here?
        res_precond_.post(x);
    }

    virtual SolverCategory::Category category() const override
    {
        return res_precond_.category();
    }

    virtual void update() override
    {
        res_precond_.update();
    }

private:
    const Dune::LinearOperator<X, Y>& well_op_;
    ReservoirPreconditioner res_precond_;
};

} // namespace Opm

#endif // OPM_WELL_OPERATOR_PRECONDITIONER_HEADER_INCLUDED
