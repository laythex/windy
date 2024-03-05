#include <dolfin.h>
#include "magnetic.h"

using namespace dolfin;

// Отсюда приходит поток частиц
class InflowDomain : public SubDomain
{
	bool inside(const Array<double> &x, bool on_boundary) const
	{
		return on_boundary && (x[0] < DOLFIN_EPS);
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
		values[0] = 100;
		values[1] = 100;
		values[2] = 100;
	}
};

int main()
{
	// Для начала просто куб
	long unsigned size = 4;
	auto mesh = std::make_shared<Mesh>(UnitCubeMesh::create({size, size, size}, CellType::Type::tetrahedron));

	// auto mesh = std::make_shared<Mesh>(Mesh("../thor.xml"));

	// Пространство функций
	auto V = std::make_shared<magnetic::FunctionSpace>(mesh);

	// Единственное граничное условие
	auto u_in = std::make_shared<Constant>(1, 0, 0);
	auto inflow_domain = std::make_shared<InflowDomain>();
	DirichletBC inflow_bc(V, u_in, inflow_domain);

	// Шаг и время симуляции
	double dt = 0.01;
	double T = 1;

	// Инициализация вещей из .hpp
	auto k = std::make_shared<Constant>(dt);
	auto B = std::make_shared<MagneticField>();
	auto u0 = std::make_shared<Function>(V);

	magnetic::BilinearForm a(V, V);
	magnetic::LinearForm L(V);

	// Хз почему только эти поля
	a.k = k;
	L.k = k;
	a.B = B;
	L.u0 = u0;

	File file("magnetic.pvd");

	// Основной цикл
	auto u1 = std::make_shared<Function>(V);
	double t = 0;
	while (t < T)
	{
		solve(a == L, *u1, inflow_bc);
		file << *u1;
		u0 = u1;
		t += dt;
		cout << "t = " << t << endl;
	}

	// Костя говножоп
	return 0;
}
