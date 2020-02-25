/*
  Copyright 2019 Norce.

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


#ifndef OPM_WELLGROUPHELPERS_HEADER_INCLUDED
#define OPM_WELLGROUPHELPERS_HEADER_INCLUDED

#include <vector>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleTypes.hpp>

namespace Opm {


    namespace wellGroupHelpers
    {

    inline void setGroupControl(const Group& group, const Schedule& schedule, const Phase& groupInjectionPhase, const int reportStepIdx, const bool injector, WellStateFullyImplicitBlackoil& wellState, std::ostringstream& ss) {

        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            setGroupControl(groupTmp, schedule, groupInjectionPhase, reportStepIdx, injector, wellState, ss);
            if (injector) {
                wellState.setCurrentInjectionGroupControl(groupInjectionPhase, groupName, Group::InjectionCMode::FLD);
            } else {
                wellState.setCurrentProductionGroupControl(groupName, Group::ProductionCMode::FLD);
            }
        }

        /*
        const auto& end = wellState.wellMap().end();
        for (const std::string& wellName : group.wells()) {
            const auto& it = wellState.wellMap().find( wellName );
            if (it == end)  // the well is not found
                continue;

            int well_index = it->second[0];
            const auto& wellEcl = schedule.getWell(wellName, reportStepIdx);

            if (wellEcl.getStatus() == Well::Status::SHUT)
                continue;

            if (!wellEcl.isAvailableForGroupControl())
                continue;

            if (wellEcl.isProducer() && !injector) {
                if (wellState.currentProductionControls()[well_index] != Well::ProducerCMode::GRUP) {
                    wellState.currentProductionControls()[well_index] = Well::ProducerCMode::GRUP;
                    ss <<"\n Producer " << wellName << " switches to GRUP control limit";
                }
            }

            if (wellEcl.isInjector() && injector) {
                // only switch if the well phase is the same as the group phase
                // Get the current controls.
                const InjectorType& injectorType = wellEcl.getInjectionProperties().injectorType;

                if (injectorType == InjectorType::WATER && groupInjectionPhase != Phase::WATER)
                    continue;

                if (injectorType == InjectorType::OIL && groupInjectionPhase != Phase::OIL)
                    continue;

                if (injectorType == InjectorType::GAS && groupInjectionPhase != Phase::GAS)
                    continue;

                if (injectorType == InjectorType::MULTI)
                    throw("Expected WATER, OIL or GAS as type for injectors " + wellEcl.name());

                if (wellState.currentInjectionControls()[well_index] != Well::InjectorCMode::GRUP) {
                    wellState.currentInjectionControls()[well_index] = Well::InjectorCMode::GRUP;
                    ss <<"\n Injector " << wellName << " switches to GRUP control limit";
                }
            }
        }
        */
    }

    inline void setCmodeGroup(const Group& group, const Schedule& schedule, const SummaryState& summaryState, const int reportStepIdx, WellStateFullyImplicitBlackoil& wellState) {

        for (const std::string& groupName : group.groups()) {
            setCmodeGroup( schedule.getGroup(groupName, reportStepIdx), schedule, summaryState, reportStepIdx, wellState);
        }

        // use NONE as default control
        const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
        for (Phase phase : all) {
            if (!wellState.hasInjectionGroupControl(phase, group.name())) {
                wellState.setCurrentInjectionGroupControl(phase, group.name(), Group::InjectionCMode::NONE);
            }
        }
        if (!wellState.hasProductionGroupControl(group.name())) {
            wellState.setCurrentProductionGroupControl(group.name(), Group::ProductionCMode::NONE);
        }

        if (group.isInjectionGroup() && schedule.hasWellGroupEvent(group.name(),  ScheduleEvents::GROUP_INJECTION_UPDATE, reportStepIdx)) {           

            for (Phase phase : all) {
                if (!group.hasInjectionControl(phase))
                    continue;

                const auto& controls = group.injectionControls(phase, summaryState);
                wellState.setCurrentInjectionGroupControl(phase, group.name(), controls.cmode);
            }
        }

        if (group.isProductionGroup() && schedule.hasWellGroupEvent(group.name(),  ScheduleEvents::GROUP_PRODUCTION_UPDATE, reportStepIdx)) {
            const auto controls = group.productionControls(summaryState);
            wellState.setCurrentProductionGroupControl(group.name(), controls.cmode);
        }

        if (schedule.gConSale(reportStepIdx).has(group.name())) {
            wellState.setCurrentInjectionGroupControl(Phase::GAS, group.name(), Group::InjectionCMode::SALE);
            std::ostringstream ss;
            setGroupControl(group, schedule, Phase::GAS, reportStepIdx, /*injector*/true, wellState, ss);
        }
    }


    inline void accumulateGroupEfficiencyFactor(const Group& group, const Schedule& schedule, const int reportStepIdx, double& factor) {
        factor *= group.getGroupEfficiencyFactor();
        if (group.parent() != "FIELD")
            accumulateGroupEfficiencyFactor(schedule.getGroup(group.parent(), reportStepIdx), schedule, reportStepIdx, factor);
    }

/*
    inline void computeGroupTargetReduction(const Group& group, const WellStateFullyImplicitBlackoil& wellState, const Schedule& schedule, const int reportStepIdx, const int phasePos, const bool isInjector, double& groupTargetReduction )
    {

        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            computeGroupTargetReduction(groupTmp, wellState, schedule, reportStepIdx, phasePos, isInjector, groupTargetReduction);
        }
        for (const std::string& wellName : group.wells()) {
            const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);

            if (wellTmp.isProducer() && isInjector)
                 continue;

            if (wellTmp.isInjector() && !isInjector)
                 continue;

            if (wellTmp.getStatus() == Well::Status::SHUT)
                continue;

            const auto& end = wellState.wellMap().end();
                const auto& it = wellState.wellMap().find( wellName );
                if (it == end)  // the well is not found
                    continue;

            int well_index = it->second[0];
            const auto wellrate_index = well_index * wellState.numPhases();
            // add contributino from wells not under group control
            if (isInjector) {
                if (wellState.currentInjectionControls()[well_index] != Well::InjectorCMode::GRUP)
                    groupTargetReduction += wellState.wellRates()[wellrate_index + phasePos];
            } else {
                if (wellState.currentProductionControls()[well_index] !=  Well::ProducerCMode::GRUP)
                    groupTargetReduction -= wellState.wellRates()[wellrate_index + phasePos];
            }
        }
    }
*/
    inline double sumWellPhaseRates(const std::vector<double>& rates, const Group& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const int phasePos,
                                    const bool injector) {

        double rate = 0.0;
        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            rate += groupTmp.getGroupEfficiencyFactor()*sumWellPhaseRates(rates, groupTmp, schedule, wellState, reportStepIdx, phasePos, injector);
        }
        const auto& end = wellState.wellMap().end();
        for (const std::string& wellName : group.wells()) {
            const auto& it = wellState.wellMap().find( wellName );
            if (it == end)  // the well is not found
                continue;

            int well_index = it->second[0];

            const auto& wellEcl = schedule.getWell(wellName, reportStepIdx);
            //only count producers or injectors
            if ( (wellEcl.isProducer() && injector) ||  (wellEcl.isInjector() && !injector))
                continue;

            if (wellEcl.getStatus() == Well::Status::SHUT)
                continue;

            double factor = wellEcl.getEfficiencyFactor();
            const auto wellrate_index = well_index * wellState.numPhases();
            if (injector)
                rate += factor * rates[ wellrate_index + phasePos];
            else
                rate -= factor * rates[ wellrate_index + phasePos];
        }
        return rate;
    }

    inline double sumWellRates(const Group& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const int phasePos, const bool injector) {
        return sumWellPhaseRates(wellState.wellRates(), group, schedule, wellState, reportStepIdx, phasePos, injector);
    }

    inline double sumWellResRates(const Group& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const int phasePos, const bool injector) {
        return sumWellPhaseRates(wellState.wellReservoirRates(), group, schedule, wellState, reportStepIdx, phasePos, injector);
    }

    inline double sumSolventRates(const Group& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const bool injector) {

        double rate = 0.0;
        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            rate += groupTmp.getGroupEfficiencyFactor()*sumSolventRates(groupTmp, schedule, wellState, reportStepIdx, injector);
        }
        const auto& end = wellState.wellMap().end();
        for (const std::string& wellName : group.wells()) {
            const auto& it = wellState.wellMap().find( wellName );
            if (it == end)  // the well is not found
                continue;

            int well_index = it->second[0];

            const auto& wellEcl = schedule.getWell(wellName, reportStepIdx);
            //only count producers or injectors
            if ( (wellEcl.isProducer() && injector) ||  (wellEcl.isInjector() && !injector))
                continue;

            if (wellEcl.getStatus() == Well::Status::SHUT)
                continue;

            double factor = wellEcl.getEfficiencyFactor();
            if (injector)
                rate += factor * wellState.solventWellRate(well_index);
            else
                rate -= factor * wellState.solventWellRate(well_index);
        }
        return rate;
    }

    inline void updateGroupTargetReduction(const Group& group, const Schedule& schedule, const int reportStepIdx, const bool isInjector, const PhaseUsage& pu, const WellStateFullyImplicitBlackoil& wellStateNupcol, WellStateFullyImplicitBlackoil& wellState, std::vector<double>& groupTargetReduction)
    {
        const int np = wellState.numPhases();
        for (const std::string& subGroupName : group.groups()) {
            std::vector<double> subGroupTargetReduction(np, 0.0);
            const Group& subGroup = schedule.getGroup(subGroupName, reportStepIdx);
            updateGroupTargetReduction(subGroup, schedule, reportStepIdx, isInjector, pu, wellStateNupcol, wellState, subGroupTargetReduction);

            // accumulate group contribution from sub group
            if (isInjector) {
                const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
                for (Phase phase : all) {
                    const Group::InjectionCMode& currentGroupControl = wellState.currentInjectionGroupControl(phase, subGroupName);
                    int phasePos;
                    if (phase == Phase::GAS && pu.phase_used[BlackoilPhases::Vapour] )
                        phasePos = pu.phase_pos[BlackoilPhases::Vapour];
                    else if (phase == Phase::OIL && pu.phase_used[BlackoilPhases::Liquid])
                        phasePos = pu.phase_pos[BlackoilPhases::Liquid];
                    else if (phase == Phase::WATER && pu.phase_used[BlackoilPhases::Aqua] )
                        phasePos = pu.phase_pos[BlackoilPhases::Aqua];
                    else
                        continue;

                    if (currentGroupControl != Group::InjectionCMode::FLD) {
                        groupTargetReduction[phasePos] += sumWellRates(subGroup, schedule, wellStateNupcol, reportStepIdx, phasePos, isInjector);
                    } else {
                        groupTargetReduction[phasePos] += subGroupTargetReduction[phasePos];
                    }
                }
            } else {
                const Group::ProductionCMode& currentGroupControl = wellState.currentProductionGroupControl(subGroupName);
                if (currentGroupControl != Group::ProductionCMode::FLD) {
                    for (int phase = 0; phase < np; phase++) {
                        groupTargetReduction[phase] += sumWellRates(subGroup, schedule, wellStateNupcol, reportStepIdx, phase, isInjector);
                    }
                } else {
                    // or accumulate directly from the wells if controled from its parents
                    for (int phase = 0; phase < np; phase++) {
                        groupTargetReduction[phase] += subGroupTargetReduction[phase];
                    }
                }
            }
        }
        for (const std::string& wellName : group.wells()) {
            const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);

            if (wellTmp.isProducer() && isInjector)
                continue;

            if (wellTmp.isInjector() && !isInjector)
                continue;

            if (wellTmp.getStatus() == Well::Status::SHUT)
                continue;

            const auto& end = wellState.wellMap().end();
            const auto& it = wellState.wellMap().find( wellName );
            if (it == end)  // the well is not found
                continue;

            int well_index = it->second[0];
            const auto wellrate_index = well_index * wellState.numPhases();
            const double efficiency = wellTmp.getEfficiencyFactor();
            // add contributino from wells not under group control
            if (isInjector) {
                if (wellState.currentInjectionControls()[well_index] != Well::InjectorCMode::GRUP)
                    for (int phase = 0; phase < np; phase++) {
                        groupTargetReduction[phase] += wellStateNupcol.wellRates()[wellrate_index + phase] * efficiency;
                    }
            } else {
                if (wellState.currentProductionControls()[well_index] !=  Well::ProducerCMode::GRUP)
                    for (int phase = 0; phase < np; phase++) {
                        groupTargetReduction[phase] -= wellStateNupcol.wellRates()[wellrate_index + phase] * efficiency;
                    }
            }
        }
        const double groupEfficiency = group.getGroupEfficiencyFactor();
        for (double& elem : groupTargetReduction) {
            elem *= groupEfficiency;
        }
        if (isInjector)
            wellState.setCurrentInjectionGroupReductionRates(group.name(), groupTargetReduction);
        else
            wellState.setCurrentProductionGroupReductionRates(group.name(), groupTargetReduction);
    }

    template <class Comm>
    inline void updateGuideRateForGroups(const Group& group, const Schedule& schedule, const PhaseUsage& pu, const int reportStepIdx, const double& simTime, const bool isInjector, WellStateFullyImplicitBlackoil& wellState, const Comm& comm, GuideRate* guideRate, std::vector<double>& pot)
    {
        const int np = pu.num_phases;
        for (const std::string& groupName : group.groups()) {
            std::vector<double> thisPot(np, 0.0);
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            updateGuideRateForGroups(groupTmp, schedule, pu, reportStepIdx, simTime, isInjector, wellState, comm, guideRate, thisPot);

            // accumulate group contribution from sub group if FLD
            if (isInjector) {
                const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
                for (Phase phase : all) {
                    const Group::InjectionCMode& currentGroupControl = wellState.currentInjectionGroupControl(phase, groupName);
                    if (currentGroupControl != Group::InjectionCMode::FLD) {
                        continue;
                    }
                    int phasePos;
                    if (phase == Phase::GAS && pu.phase_used[BlackoilPhases::Vapour] )
                        phasePos = pu.phase_pos[BlackoilPhases::Vapour];
                    else if (phase == Phase::OIL && pu.phase_used[BlackoilPhases::Liquid])
                        phasePos = pu.phase_pos[BlackoilPhases::Liquid];
                    else if (phase == Phase::WATER && pu.phase_used[BlackoilPhases::Aqua] )
                        phasePos = pu.phase_pos[BlackoilPhases::Aqua];
                    else
                        continue;

                    pot[phasePos] += thisPot[phasePos];
                }
            } else {
                const Group::ProductionCMode& currentGroupControl = wellState.currentProductionGroupControl(groupName);
                if (currentGroupControl != Group::ProductionCMode::FLD) {
                    continue;
                }
                for (int phase = 0; phase < np; phase++) {
                    pot[phase] += thisPot[phase];
                }
            }

        }
        for (const std::string& wellName : group.wells()) {
            const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);

            if (wellTmp.isProducer() && isInjector)
                continue;

            if (wellTmp.isInjector() && !isInjector)
                continue;

            if (wellTmp.getStatus() == Well::Status::SHUT)
                continue;
            const auto& end = wellState.wellMap().end();
            const auto& it = wellState.wellMap().find( wellName );
            if (it == end)  // the well is not found
                continue;

            int well_index = it->second[0];
            const auto wellrate_index = well_index * wellState.numPhases();
            // add contribution from wells not under group control
            if (isInjector) {
                if (wellState.currentInjectionControls()[well_index] == Well::InjectorCMode::GRUP)
                    for (int phase = 0; phase < np; phase++) {
                        pot[phase] += wellState.wellPotentials()[wellrate_index + phase];
                    }
            } else {
                if (wellState.currentProductionControls()[well_index] ==  Well::ProducerCMode::GRUP)
                    for (int phase = 0; phase < np; phase++) {
                        pot[phase] -= wellState.wellPotentials()[wellrate_index + phase];
                    }
            }
        }

        double oilPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Liquid])
            oilPot = pot [ pu.phase_pos[BlackoilPhases::Liquid]];

        double gasPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Vapour])
            gasPot = pot [ pu.phase_pos[BlackoilPhases::Vapour]];

        double waterPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Aqua])
            waterPot = pot [pu.phase_pos[BlackoilPhases::Aqua]];

        const double gefac = group.getGroupEfficiencyFactor();

        oilPot = comm.sum(oilPot) * gefac;
        gasPot = comm.sum(gasPot) * gefac;
        waterPot = comm.sum(waterPot) * gefac;

        if (isInjector) {
            wellState.setCurrentGroupInjectionPotentials(group.name(), pot);
        } else {
            guideRate->compute(group.name(), reportStepIdx, simTime, oilPot, gasPot, waterPot);
        }
    }

    template <class Comm>
    inline void updateGuideRatesForWells(const Schedule& schedule, const PhaseUsage& pu, const int reportStepIdx, const double& simTime, const WellStateFullyImplicitBlackoil& wellState, const Comm& comm, GuideRate* guideRate) {

        const auto& end = wellState.wellMap().end();
        for (const auto& well : schedule.getWells(reportStepIdx)) {
            double oilpot = 0.0;
            double gaspot = 0.0;
            double waterpot = 0.0;

            const auto& it = wellState.wellMap().find( well.name());
            if (it != end) {  // the well is found

                int well_index = it->second[0];

                const auto wpot = wellState.wellPotentials().data() + well_index*wellState.numPhases();
                if (pu.phase_used[BlackoilPhases::Liquid] > 0)
                    oilpot = wpot[pu.phase_pos[BlackoilPhases::Liquid]];

                if (pu.phase_used[BlackoilPhases::Vapour] > 0)
                    gaspot = wpot[pu.phase_pos[BlackoilPhases::Vapour]];

                if (pu.phase_used[BlackoilPhases::Aqua] > 0)
                    waterpot = wpot[pu.phase_pos[BlackoilPhases::Aqua]];
            }
            const double wefac = well.getEfficiencyFactor();
            oilpot = comm.sum(oilpot) * wefac;
            gaspot = comm.sum(gaspot) * wefac;
            waterpot = comm.sum(waterpot) * wefac;
            guideRate->compute(well.name(), reportStepIdx, simTime, oilpot, gaspot, waterpot);
        }

    }


    inline void updateVREPForGroups(const Group& group, const Schedule& schedule, const int reportStepIdx, const WellStateFullyImplicitBlackoil& wellStateNupcol, WellStateFullyImplicitBlackoil& wellState) {
        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            updateVREPForGroups(groupTmp, schedule, reportStepIdx, wellStateNupcol, wellState);
        }
        const int np = wellState.numPhases();
        double resv = 0.0;
        for (int phase = 0; phase < np; ++phase) {
            resv += sumWellPhaseRates(wellStateNupcol.wellReservoirRates(), group, schedule, wellState, reportStepIdx, phase, /*isInjector*/ false);
        }
        wellState.setCurrentInjectionVREPRates(group.name(), resv);
    }

    inline void updateReservoirRatesInjectionGroups(const Group& group, const Schedule& schedule, const int reportStepIdx, const WellStateFullyImplicitBlackoil& wellStateNupcol, WellStateFullyImplicitBlackoil& wellState) {
        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            updateReservoirRatesInjectionGroups(groupTmp, schedule, reportStepIdx, wellStateNupcol, wellState);
        }
        const int np = wellState.numPhases();
        std::vector<double> resv(np, 0.0);
        for (int phase = 0; phase < np; ++phase) {
            resv[phase] = sumWellPhaseRates(wellStateNupcol.wellReservoirRates(), group, schedule, wellState, reportStepIdx, phase, /*isInjector*/ true);
        }
        wellState.setCurrentInjectionGroupReservoirRates(group.name(), resv);
    }

    inline void updateREINForGroups(const Group& group, const Schedule& schedule, const int reportStepIdx, const PhaseUsage& pu, const SummaryState& st, const WellStateFullyImplicitBlackoil& wellStateNupcol, WellStateFullyImplicitBlackoil& wellState) {
        const int np = wellState.numPhases();
        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            updateREINForGroups(groupTmp, schedule, reportStepIdx, pu, st, wellStateNupcol, wellState);
        }

        std::vector<double> rein(np, 0.0);
        for (int phase = 0; phase < np; ++phase) {
            rein[phase] = sumWellPhaseRates(wellStateNupcol.wellRates(), group, schedule, wellState, reportStepIdx, phase, /*isInjector*/ false);
        }

        // add import rate and substract consumption rate for group for gas
        if (schedule.gConSump(reportStepIdx).has(group.name())) {
            const auto& gconsump = schedule.gConSump(reportStepIdx).get(group.name(), st);
            if (pu.phase_used[BlackoilPhases::Vapour]) {
                rein[pu.phase_pos[BlackoilPhases::Vapour]] += gconsump.import_rate;
                rein[pu.phase_pos[BlackoilPhases::Vapour]] -= gconsump.consumption_rate;
            }
        }

        wellState.setCurrentInjectionREINRates(group.name(), rein);
    }


    template <typename GuideRateEnumType>
    GuideRateModel::Target convertTarget(const GuideRateEnumType target)
    {
        return GuideRateModel::convert_target(target);
    }



    inline double fractionFromGuideRates(const std::string& name,
                                         const std::string& parent,
                                         const Schedule& schedule,
                                         const WellStateFullyImplicitBlackoil& wellState,
                                         const int reportStepIdx,
                                         const GuideRate* guideRate,
                                         const GuideRateModel::Target target,
                                         const bool isInjector,
                                         const bool alwaysIncludeThis = false)
    {
        double groupTotalGuideRate = 0.0;
        const Group& parentGroup = schedule.getGroup(parent, reportStepIdx);
        // Since a group may contain either wells or group, one of these
        // for-loop ranges will be empty.
        for (const std::string& wellName : parentGroup.wells()) {
            const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);
            if (wellTmp.isProducer() && isInjector)
                 continue;
            if (wellTmp.isInjector() && !isInjector)
                 continue;
            if (wellTmp.getStatus() == Well::Status::SHUT)
                continue;
            // Only count wells under group control
            const bool checkthis = (alwaysIncludeThis && wellName == name);
            if (isInjector) {
                if (!wellState.isInjectionGrup(wellName) && !checkthis)
                    continue;
            } else {
                if (!wellState.isProductionGrup(wellName) && !checkthis)
                    continue;
            }
            groupTotalGuideRate += guideRate->get(wellName, target);
        }
        for (const std::string& groupName : parentGroup.groups()) {
            const Group::ProductionCMode& currentGroupControl = wellState.currentProductionGroupControl(groupName);
            const bool checkthis = (alwaysIncludeThis && groupName == name);
            if (currentGroupControl == Group::ProductionCMode::FLD || checkthis) {
                groupTotalGuideRate += guideRate->get(groupName, target);
            }
        }

        if (groupTotalGuideRate == 0.0) {
            return 0.0;
        }

        const double myGuideRate = guideRate->get(name, target);
        return myGuideRate / groupTotalGuideRate;
    }



    inline void accumulateGroupFractionsFromGuideRates(const std::string& groupName,
                                                       const std::string& controlGroupName,
                                                       const Schedule& schedule,
                                                       const WellStateFullyImplicitBlackoil& wellState,
                                                       const int reportStepIdx,
                                                       const GuideRate* guideRate,
                                                       const GuideRateModel::Target target,
                                                       const bool isInjector,
                                                       const bool alwaysIncludeThis,
                                                       double& fraction)
    {
        const Group& group = schedule.getGroup(groupName, reportStepIdx);
        const std::string& parent = group.parent();
        if (groupName != controlGroupName) {
            fraction *= fractionFromGuideRates(groupName,
                                               parent,
                                               schedule,
                                               wellState,
                                               reportStepIdx,
                                               guideRate,
                                               target,
                                               isInjector,
                                               alwaysIncludeThis);
            accumulateGroupFractionsFromGuideRates(parent,
                                                   controlGroupName,
                                                   schedule,
                                                   wellState,
                                                   reportStepIdx,
                                                   guideRate,
                                                   target,
                                                   isInjector,
                                                   alwaysIncludeThis,
                                                   fraction);
        }
        return;
    }

    inline double groupFractionFromInjectionPotentials(const Group& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const PhaseUsage& pu, const int reportStepIdx, const Phase& injectionPhase) {
        double groupTotalGuideRate = 0.0;
        const Group& groupParent = schedule.getGroup(group.parent(), reportStepIdx);
        int phasePos;
        if (injectionPhase == Phase::GAS && pu.phase_used[BlackoilPhases::Vapour] )
            phasePos = pu.phase_pos[ pu.phase_pos[BlackoilPhases::Vapour] ];
        else if (injectionPhase == Phase::OIL && pu.phase_used[BlackoilPhases::Liquid])
            phasePos = pu.phase_pos[ pu.phase_pos[BlackoilPhases::Liquid] ];
        else if (injectionPhase == Phase::WATER && pu.phase_used[BlackoilPhases::Aqua] )
            phasePos = pu.phase_pos[ pu.phase_pos[BlackoilPhases::Aqua] ];
        else
            throw("this should not happen");

        for (const std::string& groupName : groupParent.groups()) {
            // only count group under group control from its parent
            const Group::InjectionCMode& currentGroupControl = wellState.currentInjectionGroupControl(injectionPhase, groupName);
            if (currentGroupControl != Group::InjectionCMode::FLD)
                continue;

            groupTotalGuideRate += wellState.currentGroupInjectionPotentials(groupName)[phasePos];
        }
        if (groupTotalGuideRate == 0.0)
            return 1.0;

        double groupGuideRate = wellState.currentGroupInjectionPotentials(group.name())[phasePos];
        return groupGuideRate / groupTotalGuideRate;
    }

    inline void accumulateGroupInjectionPotentialFractions(const std::string& groupName, const std::string& controlGroupName, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const PhaseUsage& pu, const int reportStepIdx, const Phase& injectionPhase, double& fraction) {
        const Group& group = schedule.getGroup(groupName, reportStepIdx);
        if (groupName != controlGroupName) {
            fraction *= groupFractionFromInjectionPotentials(group, schedule, wellState, pu, reportStepIdx, injectionPhase);
            accumulateGroupInjectionPotentialFractions(group.parent(), controlGroupName, schedule, wellState, pu, reportStepIdx, injectionPhase, fraction);
        }
        return;
    }



    class WellOrGroup
    {
    public:
        virtual const std::string& name() const = 0;
        virtual double efficiency() const = 0;
        virtual double guideRateFraction(const std::string& controlGroup,
                                         const Schedule& schedule,
                                         const WellStateFullyImplicitBlackoil& wellState,
                                         const int reportStepIdx,
                                         const GuideRate* guideRate,
                                         const GuideRateModel::Target& target,
                                         const bool isInjector,
                                         const bool alwaysIncludeThisObject) const = 0;
        virtual const std::string& parent() const = 0;
    };

    class W : public WellOrGroup
    {
        const Well& w_;
    public:
        explicit W(const Well& w)
            : w_(w)
        {
        }

        virtual const std::string& name() const override
        {
            return w_.name();
        }
        virtual double efficiency() const override
        {
            return w_.getEfficiencyFactor();
        }
        virtual double guideRateFraction(const std::string& controlGroup,
                                         const Schedule& schedule,
                                         const WellStateFullyImplicitBlackoil& wellState,
                                         const int reportStepIdx,
                                         const GuideRate* guideRate,
                                         const GuideRateModel::Target& target,
                                         const bool isInjector,
                                         const bool alwaysIncludeThisObject) const override
        {
            double fraction = fractionFromGuideRates(name(),
                                                     parent(),
                                                     schedule,
                                                     wellState,
                                                     reportStepIdx,
                                                     guideRate,
                                                     target,
                                                     isInjector,
                                                     alwaysIncludeThisObject);
            accumulateGroupFractionsFromGuideRates(parent(),
                                                   controlGroup,
                                                   schedule,
                                                   wellState,
                                                   reportStepIdx,
                                                   guideRate,
                                                   target,
                                                   isInjector,
                                                   alwaysIncludeThisObject,
                                                   fraction);
            return fraction;
        }
        virtual const std::string& parent() const override
        {
            return w_.groupName();
        }
    };

    class G : public WellOrGroup
    {
        const Group& g_;
    public:
        explicit G(const Group& g)
            : g_(g)
        {
        }

        virtual const std::string& name() const override
        {
            return g_.name();
        }
        virtual double efficiency() const override
        {
            return g_.getGroupEfficiencyFactor();
        }
        virtual double guideRateFraction(const std::string& controlGroup,
                                         const Schedule& schedule,
                                         const WellStateFullyImplicitBlackoil& wellState,
                                         const int reportStepIdx,
                                         const GuideRate* guideRate,
                                         const GuideRateModel::Target& target,
                                         const bool isInjector,
                                         const bool alwaysIncludeThisObject) const override
        {
            double fraction = fractionFromGuideRates(name(),
                                                     parent(),
                                                     schedule,
                                                     wellState,
                                                     reportStepIdx,
                                                     guideRate,
                                                     target,
                                                     isInjector,
                                                     alwaysIncludeThisObject);
            accumulateGroupFractionsFromGuideRates(parent(),
                                                   controlGroup,
                                                   schedule,
                                                   wellState,
                                                   reportStepIdx,
                                                   guideRate,
                                                   target,
                                                   isInjector,
                                                   alwaysIncludeThisObject,
                                                   fraction);
            return fraction;
        }
        virtual const std::string& parent() const override
        {
            return g_.parent();
        }
    };




    } // namespace wellGroupHelpers

}

#endif
