#include "AtmosphereModel.hpp"

AtmosphereModel::AtmosphereModel(std::shared_ptr<dolfin::Mesh> mesh) {
    auto V = std::make_shared<VelocityUpdate::FunctionSpace>(mesh);
    auto Q = std::make_shared<PressureUpdate::FunctionSpace>(mesh);

    // Давление на поверхности atm
    auto sur_p = std::make_shared<Constant>(Constants::PRESSURE_ASL);
    auto sur_v = std::make_shared<SurfaceVelocity>();
    auto space_p = std::make_shared<Constant>(Constants::PRESSURE_ASL * exp(-Constants::DISTRIB_COEFF * Constants::ATMO_HEIGHT));
    
    auto surface = std::make_shared<SurfaceBoundary>();
    auto space = std::make_shared<SpaceBoundary>();

    bcu1 = new DirichletBC(V, sur_v, surface);
    bcp1 = new DirichletBC(Q, sur_p, surface);
    bcp2 = new DirichletBC(Q, space_p, space);

    a1 = new TentativeVelocity::BilinearForm(V, V);
    L1 = new TentativeVelocity::LinearForm(V);
    a2 = new PressureUpdate::BilinearForm (Q, Q);
    L2 = new PressureUpdate::LinearForm(Q);
    a3 = new VelocityUpdate::BilinearForm(V, V);
    L3 = new VelocityUpdate::LinearForm(V);

    // Создаем все из ufl файлов. Там три таких, на три вариационные задачи (расчет скоростей, расчет давлений, корректировка скоростей)
    u0 = std::make_shared<Function>(V);
    u1 = std::make_shared<Function>(V);
    p1 = std::make_shared<Function>(Q);

    auto k = std::make_shared<Constant>(Constants::DELTA_TIME);
    auto f = std::make_shared<GravityForces>();
    auto rho = std::make_shared<Density>();

    auto rho_fun = std::make_shared<Function>(Q);
    *rho_fun = *rho;

    // Веселуха
    a1->k = k;      a1->rho = rho;
    L1->k = k;      L1->u0 = u0;    
    L1->f = f;      L1->rho = rho;
    L2->k = k;      L2->u1 = u1;
    L2->rho = rho;  a3->rho = rho;
    L3->k = k;      L3->u1 = u1;
    L3->p1 = p1;    L3->rho = rho;

    // Решаем вариационные линейные задачи
    assemble(A1, *a1);
    assemble(A2, *a2);
    assemble(A3, *a3);

    // std::filesystem::remove_all("results");
    File dfile("results/density.pvd");
    dfile << *rho_fun;
}

std::pair<std::shared_ptr<Function>, std::shared_ptr<Function>> AtmosphereModel::calculate() {
    // Compute tentative velocity step
    begin("Computing tentative velocity");
    assemble(b1, *L1);
    bcu1->apply(A1, b1);
    solve(A1, *u1->vector(), b1, "gmres", "default");
    end();

    // Pressure correction
    begin("Computing pressure correction");
    assemble(b2, *L2);
    bcp1->apply(A2, b2);
    bcp1->apply(*p1->vector());
    bcp2->apply(A2, b2);
    bcp2->apply(*p1->vector());
    solve(A2, *p1->vector(), b2, "bicgstab", prec);
    end();

    // Velocity correction
    begin("Computing velocity correction");
    assemble(b3, *L3);
    bcu1->apply(A3, b3);
    solve(A3, *u1->vector(), b3, "gmres", "default");
    end();
    

    // Move to next time step
    *u0 = *u1;

    return std::make_pair(p1, u0);
}
