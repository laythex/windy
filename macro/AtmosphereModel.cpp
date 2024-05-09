#include "AtmosphereModel.hpp"

AtmosphereModel::AtmosphereModel(std::shared_ptr<dolfin::Mesh> mesh, 
                                 std::shared_ptr<Function>& atmo_vel,
                                 std::shared_ptr<Function>& atmo_pres) {
    
    auto V = std::make_shared<VelocityUpdate::FunctionSpace>(mesh);
    auto Q = std::make_shared<PressureUpdate::FunctionSpace>(mesh);

    // Давление на поверхности
    auto vel_sur = std::make_shared<SurfaceVelocity>();
    auto pres_sur = std::make_shared<Constant>(Constants::PRESSURE_ASL);
    auto pres_space = std::make_shared<Constant>(Constants::PRESSURE_ASL * exp(-Constants::DISTRIB_COEFF * Constants::ATMO_HEIGHT));
    
    auto surface = std::make_shared<SurfaceBoundary>();
    auto space = std::make_shared<SpaceBoundary>();

    bc_vel = new DirichletBC(V, vel_sur, surface);
    bc_pres_sur = new DirichletBC(Q, pres_sur, surface);
    bc_pres_space = new DirichletBC(Q, pres_space, space);

    a1 = new TentativeVelocity::BilinearForm(V, V);
    L1 = new TentativeVelocity::LinearForm(V);
    a2 = new PressureUpdate::BilinearForm (Q, Q);
    L2 = new PressureUpdate::LinearForm(Q);
    a3 = new VelocityUpdate::BilinearForm(V, V);
    L3 = new VelocityUpdate::LinearForm(V);

    // Создаем все из ufl файлов. Там три таких, на три вариационные задачи (расчет скоростей, расчет давлений, корректировка скоростей)
    vel = std::make_shared<Function>(V);
    vel0 = std::make_shared<Function>(V);
    atmo_vel = vel;

    pres = std::make_shared<Function>(Q);
    atmo_pres = pres;

    auto k = std::make_shared<Constant>(Constants::DELTA_TIME);
    auto f = std::make_shared<GravityForces>();
    auto rho = std::make_shared<Density>();

    auto rho_fun = std::make_shared<Function>(Q);
    *rho_fun = *rho;

    // Веселуха
    a1->k = k;      a1->rho = rho;
    L1->k = k;      L1->u0 = vel0;    
    L1->f = f;      L1->rho = rho;
    L2->k = k;      L2->u1 = vel;
    L2->rho = rho;  a3->rho = rho;
    L3->k = k;      L3->u1 = vel;
    L3->p1 = pres;  L3->rho = rho;

    // Решаем вариационные линейные задачи
    assemble(A1, *a1);
    assemble(A2, *a2);
    assemble(A3, *a3);

    File dfile("results/density.pvd");
    dfile << *rho_fun;
}

void AtmosphereModel::calculate() {
    // Compute tentative velocity step
    begin("Computing tentative velocity");
    assemble(b1, *L1);
    bc_vel->apply(A1, b1);
    solve(A1, *vel->vector(), b1, "gmres", "default");
    end();

    // Pressure correction
    begin("Computing pressure correction");
    assemble(b2, *L2);
    bc_pres_sur->apply(A2, b2);
    bc_pres_sur->apply(*pres->vector());
    bc_pres_space->apply(A2, b2);
    bc_pres_space->apply(*pres->vector());
    solve(A2, *pres->vector(), b2, "bicgstab", prec);
    end();

    // Velocity correction
    begin("Computing velocity correction");
    assemble(b3, *L3);
    bc_vel->apply(A3, b3);
    solve(A3, *vel->vector(), b3, "gmres", "default");
    end();
    
    // Move to next time step
    *vel0 = *vel;
}
