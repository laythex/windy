#pragma once

#include <dolfin.h>

#include "ufl/velocity.h"
#include "ufl/concentration.h"

#include "Constants.hpp"
#include "BoundaryConditions.hpp"
#include "Domains.hpp"

using namespace dolfin;

class MagneticFieldModel{

    private:
    velocity::BilinearForm* av;
    velocity::LinearForm* Lv;

    concentration::BilinearForm* ac;
    concentration::LinearForm* Lc;

    DirichletBC* vel_bc;
    DirichletBC* conc_bc;    

    Function* vel;
    std::shared_ptr<Function> vel0;
    Function* conc;
    std::shared_ptr<Function> conc0;

    File vfile{"results/velocity.pvd"};
    File cfile{"results/concentration.pvd"};

    public:
    MagneticFieldModel(std::shared_ptr<dolfin::Mesh> mesh);
    void calculate();
};
