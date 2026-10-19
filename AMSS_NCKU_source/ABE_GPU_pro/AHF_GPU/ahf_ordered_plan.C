#include "macrodef.h"
#include "AHF/driver.h"
#include "ahf_ordered_csr.h"
#include <climits>
#include <map>
#include <stdexcept>
namespace AHFinderDirect { namespace gpu {
namespace {
int checked_size(std::size_t n) {
    if (n > INT_MAX) throw std::overflow_error("AHF ordered plan exceeds int indexing");
    return static_cast<int>(n);
}
}
OrderedPlan make_ordered_plan(patch_system &ps) {
    if (ps.N_grid_points() <= 0) throw std::invalid_argument("AHF empty ordered plan");
    ps.compute_synchronize_Jacobian();
    int global_min, global_max;
    ps.synchronize_Jacobian_global_minmax_ym(global_min, global_max);
    jtutil::array1d<fp> buffer(global_min, global_max);
    OrderedPlan plan;
    for (int pn = 0; pn < ps.N_patches(); ++pn) {
        patch &p = ps.ith_patch(pn);
        for (int r = p.min_irho(); r <= p.max_irho(); ++r) {
            for (int t = p.min_isigma(); t <= p.max_isigma(); ++t) {
                const int row = ps.gpn_of_patch_irho_isigma(p, r, t);
                if (row != checked_size(plan.rows.size()) - 1)
                    throw std::runtime_error("AHF ordered plan requires original row order");
                const int begin = checked_size(plan.instructions.size());
                std::map<int, int> slots;
                auto append = [&](int col, OrderedInstruction instruction) {
                    if (col < 0 || col >= ps.N_grid_points())
                        throw std::runtime_error("AHF ordered plan column out of range");
                    auto inserted = slots.emplace(col, -1);
                    instruction.slot = col; // Resolve to sorted absolute CSR slot below.
                    instruction.first = inserted.second;
                    plan.instructions.push_back(instruction);
                    checked_size(plan.instructions.size());
                };
                auto emit = [&](int rr, int tt, OrderedInstruction instruction,
                                const ghost_zone *gz) {
                    if (!gz) {
                        append(ps.gpn_of_patch_irho_isigma(p, rr, tt), instruction);
                        return;
                    }
                    const patch_edge &edge = gz->my_edge();
                    int yperp, position, minym, maxym;
                    const patch_edge &ye = ps.synchronize_Jacobian(
                        *gz, edge.iperp_of_irho_isigma(rr, tt),
                        edge.ipar_of_irho_isigma(rr, tt), yperp, position,
                        minym, maxym, buffer);
                    // The API may combine corner dependencies internally. Preserve
                    // its returned ascending ym order verbatim; no composed dep map.
                    for (int ym = minym; ym <= maxym; ++ym) {
                        instruction.ghost = 1;
                        instruction.ghost_factor = buffer(ym);
                        append(ps.gpn_of_patch_irho_isigma(ye.my_patch(),
                            ye.irho_of_iperp_ipar(yperp, position + ym),
                            ye.isigma_of_iperp_ipar(yperp, position + ym)), instruction);
                    }
                };
                for (int m = p.molecule_min_m(); m <= p.molecule_max_m(); ++m) {
                    const int rr = r + m;
                    emit(rr, t, {-1, OrderedTermKind::rho, 0, 2,
                        p.partial_rho_coeff(m), p.partial_rho_rho_coeff(m), 1, 0, 0},
                        p.is_in_nominal_grid(rr, t) ? nullptr : &p.minmax_rho_ghost_zone(m < 0));
                }
                for (int m = p.molecule_min_m(); m <= p.molecule_max_m(); ++m) {
                    const int tt = t + m;
                    emit(r, tt, {-1, OrderedTermKind::sigma, 1, 4,
                        p.partial_sigma_coeff(m), p.partial_sigma_sigma_coeff(m), 1, 0, 0},
                        p.is_in_nominal_grid(r, tt) ? nullptr : &p.minmax_sigma_ghost_zone(m < 0));
                }
                for (int mr = p.molecule_min_m(); mr <= p.molecule_max_m(); ++mr) {
                    for (int mt = p.molecule_min_m(); mt <= p.molecule_max_m(); ++mt) {
                        const int rr = r + mr, tt = t + mt;
                        emit(rr, tt, {-1, OrderedTermKind::mixed, 3, -1,
                            p.partial_rho_sigma_coeff(mr, mt), 0, 1, 0, 0},
                            p.is_in_nominal_grid(rr, tt) ? nullptr :
                            &p.corner_ghost_zone_containing_point(mr < 0, mt < 0, rr, tt));
                    }
                }
                for (auto &entry : slots) {
                    entry.second = checked_size(plan.columns.size());
                    plan.columns.push_back(entry.first);
                }
                for (int j = begin; j < checked_size(plan.instructions.size()); ++j)
                    plan.instructions[j].slot = slots.at(plan.instructions[j].slot);
                plan.rows.push_back(checked_size(plan.columns.size()));
                plan.instruction_rows.push_back(checked_size(plan.instructions.size()));
            }
        }
    }
    return plan;
}
}}
