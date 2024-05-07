#pragma once

#include <dolfin.h>

#include "ufl/TentativeVelocity.h"
#include "ufl/PressureUpdate.h"
#include "ufl/VelocityUpdate.h"

#include "Constants.hpp"
#include "BoundaryConditions.hpp"
#include "Domains.hpp"

using namespace dolfin;

class AtmosphereModel {

    private:
    TentativeVelocity::BilinearForm* a1;
    TentativeVelocity::LinearForm* L1;
    PressureUpdate::BilinearForm* a2;
    PressureUpdate::LinearForm* L2;
    VelocityUpdate::BilinearForm* a3;
    VelocityUpdate::LinearForm* L3;

    DirichletBC* bc_vel;
    DirichletBC* bc_pres_sur;
    DirichletBC* bc_pres_space;

    Matrix A1, A2, A3;
    Vector b1, b2, b3;

    std::shared_ptr<Function> vel0;
    std::shared_ptr<Function> vel;
    std::shared_ptr<Function> pres;

    // Звучит страшно, я это не трогал
    const std::string prec{has_krylov_solver_preconditioner("amg") ? "amg" : "default"};

    public:
    AtmosphereModel(std::shared_ptr<dolfin::Mesh> mesh, 
                    std::shared_ptr<Function>& atmo_vel,
                    std::shared_ptr<Function>& atmo_pres);
    void calculate();
};
