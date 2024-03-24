#include <dolfin.h>
#include "velocity.h"
#include "concentration.h"
#include "mshr.h"

using namespace dolfin;

// Размеры расчетной области
const double radius = 1;
const double boxsize = 2;

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
		values[0] = 1;
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
		// Дипольный момент Земли
		double mx = 20;
		double my = 20;
		double mz = 0;

		double x = coords[0], y = coords[1], z = coords[2];

		// Штуки для вычисления поля диполя
		double rm = x * mx + y * my + z * mz;
		double r2 = x * x + y * y + z * z;
		double r = sqrt(r2);

		// Магнитное поле диполя
		values[0] = (3 * x * rm / r2 - mx) / (r2 * r);
		values[1] = (3 * y * rm / r2 - my) / (r2 * r);
		values[2] = (3 * z * rm / r2 - mz) / (r2 * r);
	}
};

int main()
{
	// Шаг и длительность симуляции
	double dt = 0.01;
	double T = 1;

	// Нагло позаимствовал
	unsigned resolution = 24;
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

	// Начальные условия (совпадают с граничными)
	InflowVelocity vel_ic;
	vel = vel_ic;

	// Инициализация вещей из .hpp
	auto k = std::make_shared<Constant>(dt);
	auto B = std::make_shared<MagneticField>();
	auto vel0 = std::make_shared<Function>(V);
	*vel0 = vel_ic;

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

	InflowConcentration conc_ic;
	conc = conc_ic;

	auto conc0 = std::make_shared<Function>(C);
	*conc0 = conc_ic;

	concentration::BilinearForm ac(C, C);
	concentration::LinearForm Lc(C);

	ac.k = k;
	ac.vel = vel0;
	Lc.k = k;
	Lc.conc0 = conc0;


	File vfile("results/velocity.pvd");
	File cfile("results/concentration.pvd");

	// Основной цикл
	double t = 0;
	while (t < T)
	{
		solve(av == Lv, vel, vel_bc);
		*vel0 = vel;
		vfile << vel;

		solve(ac == Lc, conc, conc_bc);
		*conc0 = conc;
		cfile << conc;

		cout << "t = " << t << " s\t" << t / T * 1e2 << "%" << endl;
		t += dt;
	}

	// Костя говножоп
	return 0;
}
