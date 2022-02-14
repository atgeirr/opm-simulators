/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Equinor.

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

#ifndef OPM_CONVERGENCEREPORT_HEADER_INCLUDED
#define OPM_CONVERGENCEREPORT_HEADER_INCLUDED

#include <cassert>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace Opm
{

    /// Represents the convergence status of the whole simulator, to
    /// make it possible to query and store the reasons for
    /// convergence failures.
    class ConvergenceReport
    {
    public:

        // ----------- Types -----------

        enum Status { AllGood            = 0,
                      ReservoirFailed    = 1 << 0,
                      WellFailed         = 1 << 1 };
        enum struct Severity { None       = 0,
                               Normal     = 1,
                               TooLarge   = 2,
                               NotANumber = 3 };
        static std::string stringify(const Severity s)
        {
            switch (s) {
            case Severity::None: return "None";
            case Severity::Normal: return "Normal";
            case Severity::TooLarge: return "TooLarge";
            case Severity::NotANumber: return "NotANumber";
            }
        }

        class ReservoirFailure
        {
        public:
            enum struct Type { Invalid, MassBalance, Cnv };
            static std::string stringify(const Type t)
            {
                switch (t) {
                case Type::Invalid: return "Invalid";
                case Type::MassBalance: return "MassBalance";
                case Type::Cnv: return "Cnv";
                }
            }

            ReservoirFailure(Type t, Severity s, int phase, int cell_number, double magnitude)
                : type_(t), severity_(s), phase_(phase), cell_number_(cell_number), magnitude_(magnitude)
            {
            }
            Type type() const { return type_; }
            Severity severity() const { return severity_; }
            int phase() const { return phase_; }
            int cellNumber() const { return cell_number_; }
            double magnitude() const { return magnitude_; }
        private:
            Type type_;
            Severity severity_;
            int phase_;
            int cell_number_;
            double magnitude_;
        };

        class WellFailure
        {
        public:
            enum struct Type { Invalid, MassBalance, Pressure, ControlBHP, ControlTHP, ControlRate };
            static std::string stringify(const Type t)
            {
                switch (t) {
                case Type::Invalid: return "Invalid";
                case Type::MassBalance: return "MassBalance";
                case Type::Pressure: return "Pressure";
                case Type::ControlBHP: return "ControlBHP";
                case Type::ControlTHP: return "ControlTHP";
                case Type::ControlRate: return "ControlRate";
                }
            }
            WellFailure(Type t, Severity s, int phase, const std::string& well_name, double magnitude)
                : type_(t), severity_(s), phase_(phase), well_name_(well_name), magnitude_(magnitude)
            {
            }
            Type type() const { return type_; }
            Severity severity() const { return severity_; }
            int phase() const { return phase_; }
            const std::string& wellName() const { return well_name_; }
            double magnitude() const { return magnitude_; }
        private:
            Type type_;
            Severity severity_;
            int phase_;
            std::string well_name_;
            double magnitude_;
        };

        // ----------- Mutating member functions -----------

        ConvergenceReport()
            : status_{AllGood}
            , res_failures_{}
            , well_failures_{}
            , groupConverged_(true)
        {
        }

        void clear()
        {
            status_ = AllGood;
            res_failures_.clear();
            well_failures_.clear();
            groupConverged_ = true;
        }

        void setReservoirFailed(const ReservoirFailure& rf)
        {
            status_ = static_cast<Status>(status_ | ReservoirFailed);
            res_failures_.push_back(rf);
        }

        void setWellFailed(const WellFailure& wf)
        {
            status_ = static_cast<Status>(status_ | WellFailed);
            well_failures_.push_back(wf);
        }

        void setGroupConverged(const bool groupConverged)
        {
            groupConverged_ = groupConverged;
        }

        ConvergenceReport& operator+=(const ConvergenceReport& other)
        {
            status_ = static_cast<Status>(status_ | other.status_);
            res_failures_.insert(res_failures_.end(), other.res_failures_.begin(), other.res_failures_.end());
            well_failures_.insert(well_failures_.end(), other.well_failures_.begin(), other.well_failures_.end());
            assert(reservoirFailed() != res_failures_.empty());
            assert(wellFailed() != well_failures_.empty());
            return *this;
        }

        // ----------- Const member functions (queries) -----------

        bool converged() const
        {
            return status_ == AllGood && groupConverged_;
        }

        bool reservoirFailed() const
        {
            return status_ & ReservoirFailed;
        }

        bool wellFailed() const
        {
            return status_ & WellFailed;
        }

        bool groupFailed() const
        {
            return !groupConverged_;
        }

        const std::vector<ReservoirFailure>& reservoirFailures() const
        {
            return res_failures_;
        }

        const std::vector<WellFailure>& wellFailures() const
        {
            return well_failures_;
        }

        Severity severityOfWorstFailure() const
        {
            // A function to get the worst of two severities.
            auto smax = [](Severity s1, Severity s2) {
                return s1 < s2 ? s2 : s1;
            };
            auto s = Severity::None;
            for (const auto& f : res_failures_) {
                s = smax(s, f.severity());
            }
            for (const auto& f : well_failures_) {
                s = smax(s, f.severity());
            }
            return s;
        }

    private:

        // ----------- Member variables -----------
        Status status_;
        std::vector<ReservoirFailure> res_failures_;
        std::vector<WellFailure> well_failures_;
        bool groupConverged_;
    };


    inline std::ostream& operator<<(std::ostream& os, const ConvergenceReport::WellFailure& wf)
    {
        os << "||| Failure in well " << wf.wellName()
           << " of type " << ConvergenceReport::WellFailure::stringify(wf.type())
           << " and severity " << ConvergenceReport::stringify(wf.severity())
           << " for phase " << wf.phase()
           << " with magnitude " << wf.magnitude() << "\n";
        return os;
    }


    inline std::ostream& operator<<(std::ostream& os, const ConvergenceReport::ReservoirFailure& rf)
    {
        os << "||| Failure of type " << ConvergenceReport::ReservoirFailure::stringify(rf.type())
           << " and severity " << ConvergenceReport::stringify(rf.severity())
           << " for phase " << rf.phase()
           << " in cell " << rf.cellNumber()
           << " with magnitude " << rf.magnitude() << "\n";
        return os;
    }




    inline std::ostream& operator<<(std::ostream& os, const ConvergenceReport& cr)
    {
        os << "||| Convergence report: ";
        if (cr.converged()) {
            os << "All good, converged\n";
        } else {
            os << "Not converged\n";
        }
        if (cr.reservoirFailed()) {
            for (const auto& rf : cr.reservoirFailures()) {
                os << rf;
            }
        }
        if (cr.wellFailed()) {
            for (const auto& wf : cr.wellFailures()) {
                os << wf;
            }
        }
        if (cr.groupFailed()) {
            os << "||| Well groups not converged.";
        }
        return os;
    }


} // namespace Opm

#endif // OPM_CONVERGENCEREPORT_HEADER_INCLUDED
