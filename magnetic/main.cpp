#include <dolfin.h>
#include "magnetic.h"
#include "mshr.h"

using namespace dolfin;

const double radius = 1;
const double boxsize = 2;

// Отсюда приходит поток частиц
class InflowDomain : public SubDomain
{
	bool inside(const Array<double> &x, bool on_boundary) const
	{
		// return (x[0] < DOLFIN_EPS);
		return (fabs(x[0] + boxsize) < DOLFIN_EPS);
	}
};

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

class MagneticField : public Expression
{
public:
	// Конструктор чтобы Expression был vector-valued
	MagneticField() : Expression(3) {}

	void eval(Array<double> &values, const Array<double> &coords) const
	{	
		// Дипольный момент Земли
		double mx = 5;
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
	// Для начала просто куб
	// long unsigned size = 8;
	// auto mesh = std::make_shared<Mesh>(UnitCubeMesh::create({size, size, size}, CellType::Type::tetrahedron));

	// Нагло позаимствовал
    unsigned resolution = 24;
    auto universe = mshr::Box(Point(-boxsize, -boxsize, -boxsize), Point(boxsize, boxsize, boxsize));
    auto planet = mshr::Sphere(Point(0, 0, 0), radius);
	auto atmo = universe - planet;
    auto mesh = mshr::generate_mesh(atmo, resolution);

	// Пространства функций
	auto V = std::make_shared<magnetic::FunctionSpace>(mesh);

	// Домен, который отвечает за вход
	auto inflow_domain = std::make_shared<InflowDomain>();

	// Граничное условие на скорость
	auto inflow_velocity = std::make_shared<InflowVelocity>();
	DirichletBC bc(V, inflow_velocity, inflow_domain);

	// Шаг и длительность симуляции
	double dt = 0.01;
	double T = 1;

	// Искомая функция
	Function velocity(V);

	// Начальные условия (совпадают с граничными)
	InflowVelocity ic;
	velocity = ic;

	// Инициализация вещей из .hpp
	auto k = std::make_shared<Constant>(dt);
	auto B = std::make_shared<MagneticField>();
	auto velocity0 = std::make_shared<Function>(V);

	*velocity0 = ic;

	magnetic::BilinearForm a(V, V);
	magnetic::LinearForm L(V);

	a.k = k;
	a.B = B;
	L.k = k;
	L.velocity0 = velocity0;

	File vfile("results/velocity.pvd");

	// Основной цикл
	double t = 0;
	while (t < T)
	{
		solve(a == L, velocity, bc);

		*velocity0 = velocity;

		vfile << velocity;

		cout << "t = " << t << " s\t" << t / T * 1e2 << "%" << endl;
		t += dt;
	}

	// Костя говножоп
	return 0;
}
