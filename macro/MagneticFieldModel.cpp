#include "MagneticFieldModel.hpp"

MagneticFieldModel::MagneticFieldModel(std::shared_ptr<dolfin::Mesh> mesh) {
    // Пространства функций
    auto V = std::make_shared<velocity::FunctionSpace>(mesh);
    auto C = std::make_shared<concentration::FunctionSpace>(mesh);

    // Домен, который отвечает за вход
    auto inflow_domain = std::make_shared<InflowBoundary>();

    // Инициализация скорости
    auto inflow_vel = std::make_shared<InflowVelocity>();    // Граничное условие на скорость
    vel_bc = new DirichletBC(V, inflow_vel, inflow_domain);

    vel = new Function(V);    // Искомая функция

    InflowVelocity vel_ic;    // Начальные условия (совпадают с граничными)
    *vel = vel_ic;

    auto k = std::make_shared<Constant>(Constants::DELTA_TIME);    // Инициализация вещей из хедера
    auto B = std::make_shared<MagneticField>();
    vel0 = std::make_shared<Function>(V);
    *vel0 = vel_ic;

    av = new velocity::BilinearForm(V, V);
    Lv = new velocity::LinearForm(V);

    av->k = k;
    av->B = B;
    Lv->k = k;
    Lv->vel0 = vel0;

    // Инициализация концентрации
    auto inflow_conc = std::make_shared<InflowConcentration>();
    conc_bc = new DirichletBC(C, inflow_conc, inflow_domain);

    conc = new Function(C);

    InflowConcentration conc_ic;
    *conc = conc_ic;

    conc0 = std::make_shared<Function>(C);
    *conc0 = conc_ic;

    ac = new concentration::BilinearForm(C, C);
    Lc = new concentration::LinearForm(C);

    ac->k = k;
    ac->vel = vel0;
    Lc->k = k;
    Lc->conc0 = conc0;

    // уэээ
    assemble(Av, *av);
    assemble(Ac, *ac);
}

void MagneticFieldModel::calculate() {
    begin("Computing wind particles velocity");
    assemble(bv, *Lv);
    vel_bc->apply(Av, bv);
    solve(Av, *vel->vector(), bv, "gmres", "default");
    end();

    // begin("Computing wind particles concentration");
    // assemble(bc, *Lc);
    // conc_bc->apply(Ac, bc);
    // solve(Ac, *conc->vector(), bc, "gmres", "ilu");
    // end();

    vfile << *vel;
    cfile << *conc;

    *vel0 = *vel;
    *conc0 = *conc;
}
