#include <case/petsc_custom_hooks.hpp>
#include <case/structured_single_eos_reservoir_runner.hpp>
#include "case_config.hpp"
#include "well_config.hpp"
#include "viscosity.hpp"
#include <fstream>
#include <filesystem>
namespace Case {
struct Definition {using Config=H02::Config;static auto wells(){return H02::makeWellDefinitions(60,20);}};
using Runner=MPMC::cases::StructuredSingleEosReservoirRunner<Definition>;
using Indices=Runner::Indices;using Runtime=Runner::Runtime;using Grid=Runner::Grid;
int run(){
 PetscMPIInt rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
 char mode[8]="C";PetscInt nx=60,ny=20;
 PetscCallAbort(PETSC_COMM_WORLD,PetscOptionsGetString(nullptr,nullptr,"-h02_mode",mode,sizeof(mode),nullptr));
 PetscCallAbort(PETSC_COMM_WORLD,PetscOptionsGetInt(nullptr,nullptr,"-h02_nx",&nx,nullptr));
 PetscCallAbort(PETSC_COMM_WORLD,PetscOptionsGetInt(nullptr,nullptr,"-h02_ny",&ny,nullptr));
 if(std::string(mode)!="B"&&std::string(mode)!="C")throw std::invalid_argument("Only B/C are wired to this compositional runner; A needs a genuine no-transfer closure and is NOT approximated using kij.");
 if(nx<2||ny<2)throw std::invalid_argument("H02 requires nx,ny >= 2");
 auto run=MPMC::cases::readRunOptions<H02::Config>();
 if(run.dtDays*86400.>2.0000001)throw std::invalid_argument("Output dt also caps internal dt: use <=2 seconds");
 MPMC::cases::validateCaseConfig<Indices,H02::Config>();
 Grid grid(nx,ny,1,MPMC::GridExtent{.30,.10,.010});grid.setup();Runner::initializeRock(grid);
 auto fluid=MPMC::cases::makeFluidSystem<Indices,H02::Config>();
 using Eval=Indices::ValueType;char selected=mode[0];
 fluid.flowViscosityOverride=[selected](Eval,const MPMC::FluidSystem<Indices>::Composition& w){return H02::viscosity(w[0],selected);};
 auto opt=MPMC::cases::makeRuntimeOptions<Indices,Runtime,H02::Config>(run);
 opt.scaling.pressureScale=1e7;opt.scaling.compositionScale=1.;opt.scaling.saturationScale=1.;opt.scaling.massResidualScale=2e-5;opt.scaling.fugacityResidualScale=1.;opt.scaling.closureResidualScale=1.;opt.scaling.rateWellResidualFloor=2.5e-8;
 MPMC::cases::applyNaturalScalingPetscOptions(opt,rank);
 Runtime runtime(grid,fluid,opt);
 auto definitions=H02::makeWellDefinitions(nx,ny);
 auto wells=MPMC::cases::makeStructuredWells<Indices,PetscInt>(grid,definitions,[](Grid& g,const auto &d,int i,int j,int k){return Runner::wellIndex(g,i,j,k,d.radius,d.skin);});
 runtime.setWells(std::move(wells));
 Vec x=grid.createGlobalVector(Indices::numPrimaryVariables);
 runtime.initializeUniformFromPTZ(x,28e6,653.15,H02::InitialState::overallComposition);
 auto a=grid.vecGetArray<Indices::numPrimaryVariables>(x);auto own=grid.ownedRegion();
 for(const auto &d:definitions){int i=d.completion.i,j=d.completion.j;if(i>=own.xStart&&i<own.xStart+own.xCount&&j>=own.yStart&&j<own.yStart+own.yCount&&own.zStart==0)a[0][j][i][Indices::Primary::wellPressure]=d.initialBhp;}
 grid.vecRestoreArray<Indices::numPrimaryVariables>(x,a);runtime.initializeHistory(x);
 if(rank==0){std::filesystem::create_directories(run.resultDirectory);std::ofstream f(run.resultDirectory+"/h02_manifest.json");f<<"{\"mode\":\""<<mode<<"\",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"PV_m3\":7.5e-5,\"Q_target_m3_s\":2.5e-8,\"physics\":\"ASSUMED_PR76_HEAVY\",\"run_scope\":\"TIME_LIMITED_NATIVE_RUN_NOT_AUTOMATIC_2PVI_CERTIFICATION\"}\n";}
 MPMC::cases::NaturalSolver<Runtime> solver(runtime,grid.dm(Indices::numPrimaryVariables));
 MPMC::cases::runTimeLoop<Indices,Runtime,H02::Config>(runtime,solver.snes(),x,run,rank);
 PetscCallAbort(PETSC_COMM_WORLD,VecDestroy(&x));return 0;
}}
MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)
int main(int argc,char **argv){return MPMC::cases::runPetscCaseMain(argc,argv,[]{return Case::run();});}
