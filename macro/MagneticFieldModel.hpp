#pragma once

#include <dolfin.h>

#include "ufl/velocity.h"
#include "ufl/concentration.h"

#include "Constants.hpp"
#include "BoundaryConditions.hpp"
#include "Domains.hpp"

using namespace dolfin;

class MagneticFieldModel {

    private:
    velocity::BilinearForm* av;
    velocity::LinearForm* Lv;

    concentration::BilinearForm* ac;
    concentration::LinearForm* Lc;

    DirichletBC* vel_bc;
    DirichletBC* conc_bc;    

    Matrix Av, Ac;
    Vector bv, bc;

    std::shared_ptr<Function> vel;
    std::shared_ptr<Function> conc;
    std::shared_ptr<Function> vel0;
    std::shared_ptr<Function> conc0;

    public:
    MagneticFieldModel(std::shared_ptr<dolfin::Mesh> mesh, 
                       std::shared_ptr<Function>& magn_vel,
                       std::shared_ptr<Function>& magn_conc);
    void calculate();
};
