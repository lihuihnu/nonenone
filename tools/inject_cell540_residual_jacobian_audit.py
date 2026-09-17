#!/usr/bin/env python3
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one anchor, found {count}")
    p.write_text(text.replace(old, new, 1))


runtime = "models/include/natural/petsc/natural_petsc_runtime.hpp"
replace_once(runtime, "#include <functional>\n", "#include <functional>\n#include <iomanip>\n#include <iostream>\n")

residual = "models/include/natural/petsc/detail/natural_petsc_runtime_residual.inc"
residual_anchor = r'''            addCellConservationTransportAndWell<Indices>(
                block,
                outwardFlux,
                cellWell.component,
                cellWell.water);

            writeOwnedResidualBlock_(
'''
residual_probe = r'''            addCellConservationTransportAndWell<Indices>(
                block,
                outwardFlux,
                cellWell.component,
                cellWell.water);

            if constexpr (Indices::fullyCompositionalThreePhase)
            {
                static int cell540ResidualAuditCount = 0;
                const double auditTimeDay = currentTime_ / 86400.0;
                if (static_cast<long long>(cell) == 540 &&
                    auditTimeDay >= 0.00550 && auditTimeDay <= 0.0055620 &&
                    cell540ResidualAuditCount < 80)
                {
                    ++cell540ResidualAuditCount;
                    std::cerr << std::setprecision(17)
                              << "[CELL540-RJ-R] n=" << cell540ResidualAuditCount
                              << " t_day=" << auditTimeDay
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " p_bar=" << scalarValue(entry.state.pressure) / 1.0e5
                              << " So=" << scalarValue(entry.state.liquidSaturation)
                              << " Sg=" << scalarValue(entry.state.vaporSaturation)
                              << " Sw=" << scalarValue(entry.state.waterSaturation)
                              << " bits=" << static_cast<int>(entry.state.phasePresence.bits())
                              << " rows_mass=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << Indices::Equation::massConservation[static_cast<std::size_t>(c)];
                    std::cerr << " rows_og=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << Indices::Equation::fugacity[static_cast<std::size_t>(c)];
                    std::cerr << " rows_ow=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << Indices::Equation::waterFugacity[static_cast<std::size_t>(c)];
                    std::cerr << " row_closure=" << Indices::Equation::volumeClosure;
                    if constexpr (Indices::hasWellUnknown)
                        std::cerr << " row_well=" << Indices::Equation::wellControl;
                    std::cerr << " cols_p=" << Indices::Primary::pressure
                              << " cols_xo=";
                    for (int c = 0; c < Indices::numIndependentCompositionsPerPhase; ++c)
                        std::cerr << (c ? "," : "") << Indices::Primary::liquidComposition[static_cast<std::size_t>(c)];
                    std::cerr << " cols_xg=";
                    for (int c = 0; c < Indices::numIndependentCompositionsPerPhase; ++c)
                        std::cerr << (c ? "," : "") << Indices::Primary::vaporComposition[static_cast<std::size_t>(c)];
                    std::cerr << " cols_xw=";
                    for (int c = 0; c < Indices::numIndependentCompositionsPerPhase; ++c)
                        std::cerr << (c ? "," : "") << Indices::Primary::waterComposition[static_cast<std::size_t>(c)];
                    std::cerr << " cols_S=" << Indices::Primary::liquidSaturation
                              << "," << Indices::Primary::vaporSaturation
                              << "," << Indices::Primary::waterSaturation;
                    if constexpr (Indices::hasWellUnknown)
                        std::cerr << " col_well=" << Indices::Primary::wellPressure;
                    std::cerr << " Rscaled=";
                    for (int row = 0; row < Indices::numEquations; ++row)
                    {
                        const double scaled = scalarValue(block.value[static_cast<std::size_t>(row)]) *
                            equationScaleFactor_(row, localBlock);
                        std::cerr << (row ? "," : "") << scaled;
                    }
                    std::cerr << " xO=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << scalarValue(entry.state.liquidMoleFraction[static_cast<std::size_t>(c)]);
                    std::cerr << " xG=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << scalarValue(entry.state.vaporMoleFraction[static_cast<std::size_t>(c)]);
                    std::cerr << " xW=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << scalarValue(entry.state.aqueousMoleFraction[static_cast<std::size_t>(c)]);
                    std::cerr << " fO=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << scalarValue(entry.properties.fugacity[0][static_cast<std::size_t>(c)]);
                    std::cerr << " fG=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << scalarValue(entry.properties.fugacity[1][static_cast<std::size_t>(c)]);
                    std::cerr << " fW=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c ? "," : "") << scalarValue(entry.properties.fugacity[2][static_cast<std::size_t>(c)]);
                    std::cerr << '\n';
                }
            }

            writeOwnedResidualBlock_(
'''
replace_once(residual, residual_anchor, residual_probe)

jacobian = "models/include/natural/petsc/detail/natural_petsc_runtime_jacobian.inc"
jacobian_anchor = r'''            prepareFaceJacobianRowsAndAccumulateSelf_(
                cell,
                localBlock,
                cache);

            insertPreparedCellJacobianColumn_(
'''
jacobian_probe = r'''            prepareFaceJacobianRowsAndAccumulateSelf_(
                cell,
                localBlock,
                cache);

            if constexpr (Indices::fullyCompositionalThreePhase)
            {
                static int cell540JacobianAuditCount = 0;
                const double auditTimeDay = currentTime_ / 86400.0;
                if (static_cast<long long>(cell) == 540 &&
                    auditTimeDay >= 0.00550 && auditTimeDay <= 0.0055620 &&
                    cell540JacobianAuditCount < 40)
                {
                    ++cell540JacobianAuditCount;
                    constexpr int auditN = Indices::numPrimaryVariables;
                    const PetscScalar *selfValues = cellJacobianBlockValueScratch_.data();
                    std::cerr << std::setprecision(17)
                              << "[CELL540-RJ-J] n=" << cell540JacobianAuditCount
                              << " t_day=" << auditTimeDay
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " rows=" << auditN << " J=";
                    for (int row = 0; row < auditN; ++row)
                    {
                        if (row) std::cerr << ";";
                        double rowMax = 0.0;
                        double rowL1 = 0.0;
                        for (int col = 0; col < auditN; ++col)
                        {
                            const double v = PetscRealPart(selfValues[static_cast<std::size_t>(row * auditN + col)]);
                            rowMax = std::max(rowMax, std::abs(v));
                            rowL1 += std::abs(v);
                        }
                        std::cerr << "r" << row << "[max=" << rowMax << ",l1=" << rowL1 << "]=";
                        for (int col = 0; col < auditN; ++col)
                        {
                            if (col) std::cerr << ",";
                            std::cerr << PetscRealPart(selfValues[static_cast<std::size_t>(row * auditN + col)]);
                        }
                    }
                    std::cerr << '\n';
                }
            }

            insertPreparedCellJacobianColumn_(
'''
replace_once(jacobian, jacobian_anchor, jacobian_probe)

print("cell 540 residual/Jacobian audit injected")
