/**
 * @file well_source_matrix.cpp
 * @brief 编译矩阵测试：验证 `well_source_matrix` 所覆盖模板/后端组合能够实例化。
 */
#include <indices/indices.hpp>
#include <natural/physics/well_source.hpp>
#include <array>

template<class Config>
void forceWell()
{
    using I = MPMC::ScalarIndices<Config>;
    std::array<double, I::numPhases> p{}, rho{}, mob{};
    std::array<double, I::numPhases> rhoS{}, frac{};
    std::array<double, I::numComponents> inj{};
    std::array<std::array<double, I::numComponents>, I::numPhases> mf{};
    rhoS.fill(1.0);
    (void)MPMC::computePerforationWellSource<I>(
        MPMC::WellType::Producer, 1.0, 0.0,
        p, rho, mob, rhoS, frac, inj, mf,
        0.0, I::hasAqueousCO2Dissolution ? 0 : -1);
}

using Base = MPMC::CompositionalModelConfig<6,true,true,false,false,false>;
using Diss = MPMC::CompositionalModelConfig<6,true,true,true,false,false>;
using NoWater = MPMC::CompositionalModelConfig<6,false,true,false,false,false>;
using FullyCompositional = MPMC::CompositionalModelConfig<
    4,true,true,false,false,false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;

template void forceWell<Base>();
template void forceWell<Diss>();
template void forceWell<NoWater>();
template void forceWell<FullyCompositional>();
int main(){return 0;}
