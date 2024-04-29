#include <dolfin.h>
#include <cmath>
#include <vector>

using namespace dolfin;

namespace atm{

    double rad;
    double height;
    double G;
    double K;
    double rho_0;
    double p_0;

    // Вектор вращения планеты
    std::vector<double> w = {0, 0, 1};

    //Линейные формы для задач
    TentativeVelocity::BilinearForm* a1;
    TentativeVelocity::LinearForm* L1;
    PressureUpdate::BilinearForm* a2;
    PressureUpdate::LinearForm* L2;
    VelocityUpdate::BilinearForm* a3;
    VelocityUpdate::LinearForm* L3;

    DirichletBC* bcu1;
    DirichletBC* bcp1;
    DirichletBC* bcp2;

    Matrix A1, A2, A3;
    Vector b1, b2, b3;

    std::shared_ptr<Function> u0;
    std::shared_ptr<Function> u1;
    std::shared_ptr<Function> p1;

    // Звучит страшно, я это не трогал
    const std::string prec(has_krylov_solver_preconditioner("amg") ? "amg" : "default");

    File ufile("results/velocity.pvd");
    File pfile("results/pressure.pvd");

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


    void init(std::shared_ptr<dolfin::Mesh> mesh, double dt, double rad_, double height_, double G_, double K_, double rho_0_, double p_0_)
    {
        rad = rad_;
        height = height_;
        G = G_;
        K = K_;
        rho_0 = rho_0_;
        p_0 = p_0_;

        auto V = std::make_shared<VelocityUpdate::FunctionSpace>(mesh);
        auto Q = std::make_shared<PressureUpdate::FunctionSpace>(mesh);

        // Давление на поверхности atm
        auto sur_p = std::make_shared<Constant>(p_0);
        auto sur_v = std::make_shared<SurfaceVelocity>();

        auto space_p = std::make_shared<Constant>(p_0 * exp(-K * (height)));
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

        auto k = std::make_shared<Constant>(dt);
        auto f = std::make_shared<GravityForces>();
        auto rho = std::make_shared<Density>();

        auto rho_fun = std::make_shared<Function>(Q);
        *rho_fun = *rho;

        a1->k = k;
        a1->rho = rho;
        L1->k = k;
        L1->u0 = u0;
        L1->f = f;
        L1->rho = rho;
        L2->k = k;
        L2->u1 = u1;
        L2->rho = rho;
        a3->rho = rho;
        L3->k = k;
        L3->u1 = u1;
        L3->p1 = p1;
        L3->rho = rho;

        // Решаем вариационные линейные задачи
        
        assemble(A1, *a1);
        assemble(A2, *a2);
        assemble(A3, *a3);

        // std::filesystem::remove_all("results");
        File dfile("results/density.pvd");
        dfile << *rho_fun;
    }


    void calc(){
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

        

        // Save to file
        ufile << *u1;
        pfile << *p1;

        // Move to next time step
        *u0 = *u1;
    }
}

