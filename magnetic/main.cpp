#include <dolfin.h>
#include "magnetic.h"

using namespace dolfin;

// Отсюда приходит поток частиц
class InflowDomain : public SubDomain
{
	bool inside(const Array<double> &x, bool on_boundary) const
	{
		return x[0] < DOLFIN_EPS;
	}
};

class WholeDomain : public SubDomain
{
	bool inside(const Array<double> &x, bool on_boundary) const
	{
		return true;
	}
};

class MagneticField : public Expression
{
public:

	// Конструктор чтобы Expression был vector-valued
	MagneticField() : Expression(3) {}

	void eval(Array<double> &values, const Array<double> &x) const
	{
		// Здесь можно поле диполя для начала сделать
		values[0] = 0;
		values[1] = 0;
		values[2] = 0;
	}
};

class InflowVelocity : public Expression
{
public:

	// Конструктор чтобы Expression был vector-valued
	InflowVelocity() : Expression(3) {}

	void eval(Array<double> &values, const Array<double> &x) const
	{
		// Здесь можно поле диполя для начала сделать
		values[0] = 10;
		values[1] = 0;
		values[2] = 0;
	}
};

int main()
{
	// Для начала просто куб
	long unsigned size = 4;
	auto mesh = std::make_shared<Mesh>(UnitCubeMesh::create({size, size, size}, CellType::Type::tetrahedron));

	// auto mesh = std::make_shared<Mesh>(Mesh("../thor.xml"));

	// Пространства функций
	auto M = std::make_shared<magnetic::FunctionSpace>(mesh);

	// Домен, который отвечает за вход
	auto inflow_domain = std::make_shared<InflowDomain>();
	auto whole_domain = std::make_shared<WholeDomain>();

	// Граничное условие на скорость
	auto v_in = std::make_shared<InflowVelocity>();
	DirichletBC inflow_velocity_bc(M->sub(0), v_in, inflow_domain);

	// Граничное условие на плотность
	auto d_in = std::make_shared<Constant>(1);
	DirichletBC inflow_density_bc(M->sub(1), d_in, whole_domain);

	// Вектор из граничных условий
	std::vector<const DirichletBC*> bcs = {&inflow_velocity_bc, 
										   &inflow_density_bc};

	// Шаг и длительность симуляции
	double dt = 0.01;
	double T = 60;

	// Инициализация вещей из .hpp
	auto k = std::make_shared<Constant>(dt);
	auto B = std::make_shared<MagneticField>();
	auto velocity0 = std::make_shared<Function>(M->sub(0)->collapse());

	magnetic::BilinearForm a(M, M);
	magnetic::LinearForm L(M);

	a.k = k;
	// a.B = B;
	a.velocity0 = velocity0;
	L.k = k;
	L.velocity0 = velocity0;

	// Искомые функции
	Function m(M);

	// Основной цикл
	double t = 0;
	while (t < T)
	{
		solve(a == L, m, bcs);
		velocity0 = std::make_shared<Function>(m[0]);

		t += dt;
		cout << "t = " << t << endl;
	}

	File file("magnetic.pvd");
	file << m[0];

	// Костя говножоп
	return 0;
}
