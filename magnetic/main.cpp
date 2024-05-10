#include <dolfin.h>
#include <mshr.h>
#include <boost/filesystem.hpp>

#include "velocity.h"
#include "concentration.h"

using namespace dolfin;

// Средний радиус Марса (Volumetric mean radius (km)) в метрах
// https://web.archive.org/web/20200317184127/https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
const double radius = 3389.5e3;

// Размеры расчетной области
const double boxsize = radius * 2;

// Отсюда приходит поток частиц
class InflowDomain : public SubDomain
{
	bool inside(const Array<double> &x, bool on_boundary) const
	{
		return (fabs(x[0] + boxsize) < DOLFIN_EPS);
	}
};

// Скорость частиц на входе
class InflowVelocity : public Expression
{
public:
	// Конструктор чтобы Expression был vector-valued
	InflowVelocity() : Expression(3) {}

	void eval(Array<double> &values, const Array<double> &x) const
	{
		// https://solarscience.msfc.nasa.gov/SolarWind.shtml
		values[0] = 800 * 1e3;
		values[1] = 0;
		values[2] = 0;
	}
};

// Концентрация частиц на входе
class InflowConcentration : public Expression
{
public:
	void eval(Array<double> &values, const Array<double> &x) const
	{
		// Относительная величина
		values[0] = 1;
	}
};

// Магнитное поле Земли
class MagneticField : public Expression
{
public:
	// Конструктор чтобы Expression был vector-valued
	MagneticField() : Expression(3) {}

	void eval(Array<double> &values, const Array<double> &coords) const
	{
		// Дипольный момент Марса
		// https://www.britannica.com/science/geomagnetic-field/Dipolar-field
		// https://www.sciencedirect.com/science/article/pii/S0032063397001864
		double m = 8.22e15;

		// Угол между экваториальной плоскостью Марса и плоскостью его орбиты (Obliquity to orbit (deg)) в радианах
		// https://web.archive.org/web/20200317184127/https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
		// Надо еще учесть угол между магнитным и обычным полюсами
		double e = 25.19 * M_PI / 180;

		// День летнего солнцестояния
		double mx = m * sin(e);
		double my = m * cos(e);
		double mz = 0;

		double x = coords[0], y = coords[1], z = coords[2];

		// Штуки для вычисления поля диполя
		double rm = x * mx + y * my + z * mz;
		double r2 = x * x + y * y + z * z;
		double r3 = r2 * sqrt(r2);
		double rm2 = 3 * rm / r2;

		// Магнитное поле диполя
		values[0] = (x * rm2 - mx) / r3;
		values[1] = (y * rm2 - my) / r3;
		values[2] = (z * rm2 - mz) / r3;
	}
};

int main()
{
	// Шаг и длительность симуляции
	double dt = 0.01;
	double T = 1;

	// Нагло позаимствовал
	unsigned resolution = 30;
	auto universe = mshr::Box(Point(-boxsize, -boxsize, -boxsize), Point(boxsize, boxsize, boxsize));
	auto planet = mshr::Sphere(Point(0, 0, 0), radius);
	auto atmosphere = universe - planet;
	auto mesh = mshr::generate_mesh(atmosphere, resolution);

	// Пространства функций
	auto V = std::make_shared<velocity::FunctionSpace>(mesh);
	auto C = std::make_shared<concentration::FunctionSpace>(mesh);

	// Домен, который отвечает за вход
	auto inflow_domain = std::make_shared<InflowDomain>();

	// Граничное условие на скорость
	auto inflow_vel = std::make_shared<InflowVelocity>();
	DirichletBC vel_bc(V, inflow_vel, inflow_domain);

	// Искомая функция
	Function vel(V);
	vel = *inflow_vel;

	// Инициализация вещей из .hpp
	auto k = std::make_shared<Constant>(dt);
	auto B = std::make_shared<MagneticField>();
	auto vel0 = std::make_shared<Function>(V);
	*vel0 = *inflow_vel;

	velocity::BilinearForm av(V, V);
	velocity::LinearForm Lv(V);

	av.k = k;
	av.B = B;
	Lv.k = k;
	Lv.vel0 = vel0;

	// Практически то же самое с концентрацией
	auto inflow_conc = std::make_shared<InflowConcentration>();
	DirichletBC conc_bc(C, inflow_conc, inflow_domain);

	Function conc(C);
	conc = *inflow_conc;

	auto conc0 = std::make_shared<Function>(C);
	*conc0 = *inflow_conc;

	concentration::BilinearForm ac(C, C);
	concentration::LinearForm Lc(C);

	ac.k = k;
	ac.vel = vel0;
	Lc.k = k;
	Lc.conc0 = conc0;

    Matrix Av, Ac;
    Vector bv, bc;

    // уэээ
    assemble(Av, av);
    assemble(Ac, ac);

    boost::filesystem::remove_all("results/");
	File vfile("results/velocity.pvd");
	File cfile("results/concentration.pvd");

	// Основной цикл
	double t = 0;
	while (t <= T + DOLFIN_EPS) {
		begin("Computing wind particles velocity");
		assemble(bv, Lv);
		vel_bc.apply(Av, bv);
		solve(Av, *vel.vector(), bv, "gmres", "default");
		end();

		begin("Computing wind particles concentration");
		assemble(bc, Lc);
		conc_bc.apply(Ac, bc);
		solve(Ac, *conc.vector(), bc, "gmres", "default");
		end();

		*vel0 = vel;
		*conc0 = conc;

		vfile << vel;
		cfile << conc;
		
		std::cout << "t = " << t << " s\t" << t / T * 1e2 << "%" << std::endl;
		t += dt;
	}

	// Костя говножоп
	return 0;
}
