/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_TIMEDISCRETIZATION_HEADER_INCLUDED
#define OPM_TIMEDISCRETIZATION_HEADER_INCLUDED

#include <opm/simulators/timestepping/SimulatorTimer.hpp>

namespace Opm
{

    class ForwardingSimulatorTimer : public SimulatorTimerInterface
    {
    public:
        explicit ForwardingSimulatorTimer(SimulatorTimerInterface& outer_timer)
            : outer_timer_(outer_timer)
        {
        }
        /// Current step number. This is the number of timesteps that
        /// has been completed from the start of the run. The time
        /// after initialization but before the simulation has started
        /// is timestep number zero.
        virtual int currentStepNum() const
        {
            return outer_timer_.currentStepNum();
        }

        /// Current report step number. This might differ from currentStepNum in case of sub stepping
        virtual int reportStepNum() const
        {
            return outer_timer_.reportStepNum();
        }

        /// Current step length. This is the length of the step
        /// the simulator will take in the next iteration.
        ///
        /// @note if done(), it is an error to call currentStepLength().
        virtual double currentStepLength() const
        {
            return outer_timer_.currentStepLength();
        }

        /// Previous step length. This is the length of the step that
        /// was taken to arrive at this time.
        ///
        /// @note if no increments have been done (i.e. the timer is
        /// still in its constructed state and currentStepNum() == 0),
        /// it is an error to call stepLengthTaken().
        virtual double stepLengthTaken() const
        {
            return outer_timer_.stepLengthTaken();
        }

        /// Previous report step length. This is the length of the step that
        /// was taken to arrive at this report step time.
        ///
        /// @note if no increments have been done (i.e. the timer is
        /// still in its constructed state and reportStepNum() == 0),
        /// it is an error to call stepLengthTaken().
        virtual double reportStepLengthTaken() const
        {
            return outer_timer_.reportStepLengthTaken();
        }

        /// Time elapsed since the start of the simulation until the
        /// beginning of the current time step [s].
        virtual double simulationTimeElapsed() const
        {
            return outer_timer_.simulationTimeElapsed();
        }

        /// advance time by currentStepLength
        virtual void advance()
        {
            outer_timer_.advance();
        }

        /// Return true if timer indicates that simulation of timer interval is finished
        virtual bool done() const
        {
            return outer_timer_.done();
        }

        /// Return start date of simulation
        virtual boost::posix_time::ptime startDateTime() const
        {
            return outer_timer_.startDateTime();
        }

        /// Return the current time as a posix time object.
        virtual boost::posix_time::ptime currentDateTime() const
        {
            return outer_timer_.currentDateTime();
        }

        /// Time elapsed since the start of the POSIX epoch (Jan 1st,
        /// 1970) until the current time step begins [s].
        virtual time_t currentPosixTime() const
        {
            return outer_timer_.currentPosixTime();
        }

        /// Return true if last time step failed
        virtual bool lastStepFailed() const
        {
            return outer_timer_.lastStepFailed();
        }

        /// return copy of current timer instance
        virtual std::unique_ptr< SimulatorTimerInterface > clone() const
        {
            return std::make_unique<ForwardingSimulatorTimer>(*this);
        }

    private:
        SimulatorTimerInterface& outer_timer_;
    };





    class SubStepSimulatorTimer : public ForwardingSimulatorTimer
    {
    public:
        explicit SubStepSimulatorTimer(SimulatorTimerInterface& outer_timer)
            : ForwardingSimulatorTimer(outer_timer)
        {
        }

        void setSteps(std::vector<double> substep_lengths)
        {
            substep_lengths_ = std::move(substep_lengths);
        }

        virtual double currentStepLength() const
        {
            return substep_lengths_[substep_];
        }

        /// advance time by currentStepLength
        virtual void advance()
        {
            ++substep_;
        }

        /// Return true if timer indicates that simulation of timer interval is finished
        virtual bool done() const
        {
            return substep_ >= substep_lengths_.size();
        }
    private:
        std::vector<double> substep_lengths_;
        size_t substep_;
    };






    template <class Solver>
    class ImplicitEuler
    {
    public:
        ImplicitEuler(Solver& solver)
            : solver_(solver)
        {
        }

        template <class ReservoirState, class WellState, class Output>
        SimulatorReport
        step(SimulatorTimerInterface& timer,
             ReservoirState& reservoir_state,
             WellState& well_state)
        {
            return solver_.step(timer,
                                reservoir_state,
                                well_state);
        }

    private:
        Solver& solver_;
    };






    template <class Solver>
    class TwoSteps
    {
    public:
        TwoSteps(Solver& solver)
            : solver_(solver)
        {
        }

        template <class ReservoirState, class WellState, class Output>
        SimulatorReport
        step(SimulatorTimerInterface& timer,
             ReservoirState& reservoir_state,
             WellState& well_state)
        {
            SubStepSimulatorTimer sub_timer(timer);
            const double dt = timer.currentStepLength();
            sub_timer.setSteps({dt/2.0, dt/2.0});
            SimulatorReport rep = solver_.step(sub_timer,
                                               reservoir_state,
                                               well_state);
            sub_timer.advance();
            rep += solver_.step(sub_timer,
                                reservoir_state,
                                well_state);
            return rep;
        }

    private:
        Solver& solver_;
    };


} // namespace Opm

#endif // OPM_TIMEDISCRETIZATION_HEADER_INCLUDED
