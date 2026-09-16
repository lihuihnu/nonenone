/**
 * @file dmdacontainer.cpp
 * @brief StructuredGrid 的 PETSc DMDA 对象所有权与访问封装的实现。
 */
#include <structuredgrid/dmdacontainer.hpp>

#include <stdexcept>

namespace MPMC
{
namespace
{

void validateConfiguration(int nx,
                           int ny,
                           int nz,
                           int stencilWidth,
                           MPI_Comm comm)
{
    if (nx <= 0 || ny <= 0 || nz <= 0)
    {
        throw std::invalid_argument(
            "DMDAContainer requires positive grid dimensions.");
    }
    if (stencilWidth < 0)
    {
        throw std::invalid_argument(
            "DMDAContainer requires a non-negative stencil width.");
    }
    if (comm == MPI_COMM_NULL)
    {
        throw std::invalid_argument(
            "DMDAContainer requires a valid MPI communicator.");
    }
}

void validateDof(int dof)
{
    if (dof <= 0)
    {
        throw std::invalid_argument(
            "DMDAContainer requires a positive dof.");
    }
}

} // namespace

DMDAContainer::DMDAContainer(int nx,
                             int ny,
                             int nz,
                             int stencilWidth,
                             DMBoundaryType boundaryType,
                             DMDAStencilType stencilType,
                             MPI_Comm comm)
    : nx_(nx),
      ny_(ny),
      nz_(nz),
      stencilWidth_(stencilWidth),
      boundaryType_(boundaryType),
      stencilType_(stencilType),
      comm_(comm)
{
    validateConfiguration(nx_, ny_, nz_, stencilWidth_, comm_);
}

void DMDAContainer::registerLayout(Dof dof)
{
    validateDof(dof);
    if (contains(dof))
    {
        return;
    }

    if (layouts_.empty())
    {
        auto [iterator, inserted] = layouts_.try_emplace(
            dof,
            nx_,
            ny_,
            nz_,
            dof,
            stencilWidth_,
            boundaryType_,
            stencilType_,
            comm_);
        (void)inserted;

        try
        {
            validatePrimaryLayout_(iterator->second);
            canonicalDof_ = dof;
        }
        catch (...)
        {
            layouts_.erase(iterator);
            throw;
        }
        return;
    }

    const DMDAElem &canonical = layouts_.at(canonicalDof_);
    layouts_.try_emplace(dof, canonical, dof);
}

void DMDAContainer::validatePrimaryLayout_(const DMDAElem &layout) const
{
    if (layout.dimensions() != requestedDimensions())
    {
        throw std::runtime_error(
            "PETSc options changed the DMDA global dimensions. "
            "MPMC StructuredGrid requires topology dimensions to remain equal "
            "to the dimensions supplied by the grid model.");
    }

    PetscInt stencilWidth = 0;
    DMBoundaryType bx = DM_BOUNDARY_NONE;
    DMBoundaryType by = DM_BOUNDARY_NONE;
    DMBoundaryType bz = DM_BOUNDARY_NONE;
    DMDAStencilType stencilType = DMDA_STENCIL_STAR;

    PetscCallAbort(
        comm_,
        DMDAGetInfo(layout.dm(),
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    &stencilWidth,
                    &bx,
                    &by,
                    &bz,
                    &stencilType));

    if (stencilWidth != static_cast<PetscInt>(stencilWidth_) ||
        bx != boundaryType_ || by != boundaryType_ || bz != boundaryType_ ||
        stencilType != stencilType_)
    {
        throw std::runtime_error(
            "PETSc options changed the DMDA boundary/stencil topology. "
            "MPMC StructuredGrid permits processor-layout options but requires "
            "the physical grid topology to remain unchanged.");
    }
}

} // namespace MPMC
