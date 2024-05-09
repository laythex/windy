#include "MagneticFieldModel.hpp"

MagneticFieldModel::MagneticFieldModel(std::shared_ptr<dolfin::Mesh> mesh, 
                                       std::shared_ptr<Function>& magn_vel,
                                       std::shared_ptr<Function>& magn_conc) {

    // Пространства функций
    auto V = std::make_shared<velocity::FunctionSpace>(mesh);
    auto C = std::make_shared<concentration::FunctionSpace>(mesh);

    // Домен, который отвечает за вход
    auto inflow_domain = std::make_shared<InflowBoundary>();

    // Инициализация скорости
    auto inflow_vel = std::make_shared<InflowVelocity>();    // Граничное условие на скорость
    vel_bc = std::make_shared<DirichletBC>(V, inflow_vel, inflow_domain);

    vel = std::make_shared<Function>(V);   // Указатель на искомую функцию
    vel0 = std::make_shared<Function>(V);
    magn_vel = vel;   // Инициализиурем указатель пришедший извне

    *vel = *inflow_vel;    // Начальные условия (совпадают с граничными)
    *vel0 = *inflow_vel;

    auto k = std::make_shared<Constant>(Constants::DELTA_TIME);    // Инициализация вещей из хедера
    auto B = std::make_shared<MagneticField>();

    av = std::make_shared<velocity::BilinearForm>(V, V);
    Lv = std::make_shared<velocity::LinearForm>(V);

    av->k = k;
    av->B = B;
    Lv->k = k;
    Lv->vel0 = vel0;

    // Инициализация концентрации
    auto inflow_conc = std::make_shared<InflowConcentration>();
    conc_bc =  std::make_shared<DirichletBC>(C, inflow_conc, inflow_domain);

    conc = std::make_shared<Function>(C);
    conc0 = std::make_shared<Function>(C);
    magn_conc = conc;

    *conc = *inflow_conc;
    *conc0 = *inflow_conc;

    ac = std::make_shared<concentration::BilinearForm>(C, C);
    Lc = std::make_shared<concentration::LinearForm>(C);

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

    begin("Computing wind particles concentration");
    assemble(bc, *Lc);
    conc_bc->apply(Ac, bc);
    solve(Ac, *conc->vector(), bc, "gmres", "default");
    end();

    *vel0 = *vel;
    *conc0 = *conc;
}
