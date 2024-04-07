#include <dolfin.h>
#include <cmath>
#include <vector>
#include "TentativeVelocity.h"
#include "PressureUpdate.h"
#include "VelocityUpdate.h"
#include "mshr.h"
#include <filesystem>

using namespace dolfin;

// Радиус планеты
const double rad = 0.5;
// Высота атмосферы
const double height = 1.5;
// Коэффициент силы притяжения
const double G = 0.01;
// Коэффициент плотности в экспоненте
const double K = 5;
// Плотность на поверхности планеты
const double rho_0 = 1;
// Атмосферное давление
const double atm = 0.1;

// Вектор вращения планеты
std::vector<double> w = {0, 0, 1};
 
// Скорость у поверхности
class SurfaceVelocity : public Expression
{
public:
    SurfaceVelocity() : Expression(3) {}

    void eval(Array<double> &values, const Array<double> &coord) const
    {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        values[0] = w[1] * z - y * w[2];
        values[1] = -w[0] * z + w[2] * x;
        values[2] = w[0] * y - x * w[1];
    }
};

// Объемные силы притяжения

class GravityForces : public Expression
{
public:
    GravityForces() : Expression(3) {}

    void eval(Array<double> &values, const Array<double> &coord) const
    {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        double r = sqrt(x * x + y * y + z * z);
        values[0] = -G * x / (r * r * r);
        values[1] = -G * y / (r * r * r);
        values[2] = -G * z / (r * r * r);
    }
};

class Density : public Expression
{
public:
    void eval(Array<double> &values, const Array<double> &coord) const
    {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        double r = sqrt(x * x + y * y + z * z);
        if (r >= rad - DOLFIN_EPS)
        {
            values[0] = rho_0 * exp(-K * (r - rad));
        }
        else
        {
            values[0] = rho_0;
        }
    }
};

// Граница поверхности

class SurfaceBoundary : public SubDomain
{
    bool inside(const Array<double> &coord, bool on_boundary) const
    {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        return sqrt(x * x + y * y + z * z) < rad + DOLFIN_EPS;
    }
};

// Граница всей сетки

class SpaceBoundary : public SubDomain
{
    bool inside(const Array<double> &coord, bool on_boundary) const
    {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        return sqrt(x * x + y * y + z * z) >= rad + height - rad / 5;
    }
};

int main()
{
    unsigned resolution = 40;
    auto universe = mshr::Sphere(Point(0, 0, 0), rad + height);
    auto planet = mshr::Sphere(Point(0, 0, 0), rad);
    auto atmosphere = universe - planet;
    auto mesh = mshr::generate_mesh(atmosphere, resolution);

    auto V = std::make_shared<VelocityUpdate::FunctionSpace>(mesh);
    auto Q = std::make_shared<PressureUpdate::FunctionSpace>(mesh);

    double dt = 0.01;
    double T = 5;

    // Давление на поверхности atm
    auto sur_p = std::make_shared<Constant>(atm);
    auto sur_v = std::make_shared<SurfaceVelocity>();

    auto space_p = std::make_shared<Constant>(atm*exp(-K * (height)));

    auto surface = std::make_shared<SurfaceBoundary>();
    auto space = std::make_shared<SpaceBoundary>();

    DirichletBC bcu1(V, sur_v, surface);
    DirichletBC bcp1(Q, sur_p, surface);
    DirichletBC bcp2(Q, space_p, space);

    // Создаем все из ufl файлов. Там три таких, на три вариационные задачи (расчет скоростей, расчет давлений, корректировка скоростей)

    auto u0 = std::make_shared<Function>(V);
    auto u1 = std::make_shared<Function>(V);
    auto p1 = std::make_shared<Function>(Q);

    auto k = std::make_shared<Constant>(dt);
    auto f = std::make_shared<GravityForces>();
    auto rho = std::make_shared<Density>();

    auto rho_fun = std::make_shared<Function>(Q);
    *rho_fun = *rho;

    TentativeVelocity::BilinearForm a1(V, V);
    TentativeVelocity::LinearForm L1(V);
    PressureUpdate::BilinearForm a2(Q, Q);
    PressureUpdate::LinearForm L2(Q);
    VelocityUpdate::BilinearForm a3(V, V);
    VelocityUpdate::LinearForm L3(V);

    a1.k = k;
    a1.rho = rho;
    L1.k = k;
    L1.u0 = u0;
    L1.f = f;
    L1.rho = rho;
    L2.k = k;
    L2.u1 = u1;
    L2.rho = rho;
    a3.rho = rho;
    L3.k = k;
    L3.u1 = u1;
    L3.p1 = p1;
    L3.rho = rho;

    // Решаем вариационные линейные задачи
    Matrix A1, A2, A3;
    assemble(A1, a1);
    assemble(A2, a2);
    assemble(A3, a3);

    Vector b1, b2, b3;

    // Звучит страшно, я это не трогал
    const std::string prec(has_krylov_solver_preconditioner("amg") ? "amg" : "default");

    // std::filesystem::remove_all("results");

    File ufile("results/velocity.pvd");
    File pfile("results/pressure.pvd");
    File dfile("results/density.pvd");

    double t = dt;
    dfile << *rho_fun;

    // dfile << *rho;

    while (t < T + DOLFIN_EPS)
    {
        // Compute tentative velocity step
        begin("Computing tentative velocity");
        assemble(b1, L1);
        bcu1.apply(A1, b1);
        solve(A1, *u1->vector(), b1, "gmres", "default");
        end();

        // Pressure correction
        begin("Computing pressure correction");
        assemble(b2, L2);
        bcp1.apply(A2, b2);
        bcp1.apply(*p1->vector());
        bcp2.apply(A2, b2);
        bcp2.apply(*p1->vector());
        solve(A2, *p1->vector(), b2, "bicgstab", prec);
        end();

        // Velocity correction
        begin("Computing velocity correction");
        assemble(b3, L3);
        bcu1.apply(A3, b3);
        solve(A3, *u1->vector(), b3, "gmres", "default");
        end();

        // Save to file
        ufile << *u1;
        pfile << *p1;

        // Move to next time step
        *u0 = *u1;
        t += dt;
        cout << "t = " << t << endl;
    }

    return 0;
}