#include <dolfin.h>
#include <cmath>
#include <vector>
#include "TentativeVelocity.h"
#include "PressureUpdate.h"
#include "VelocityUpdate.h"
#include "mshr.h"
#include <filesystem>
#include "atmosphere.h"
#include "magnetic.h"

using namespace dolfin;


int main()
{   
    // Радиус планеты
    const double rad = 0.5;
    // Высота атмосферы
    const double height = 3;
    // Коэффициент силы притяжения
    const double G = 0.05;
    // Коэффициент плотности в экспоненте
    const double K = 5;
    // Плотность на поверхности планеты
    const double rho_0 = 1;
    // Атмосферное давление
    const double p_0 = 0.1;

    auto planet  = mshr::Sphere(Point(0, 0, 0), rad);

    //Сетка для атмосферы
    unsigned resolution_atm = 30;
    auto universe_atm  = mshr::Sphere(Point(0, 0, 0), rad + height);
    auto atmosphere_atm  = universe_atm - planet;
    auto mesh_atm = mshr::generate_mesh(atmosphere_atm, resolution_atm);

    //Сетка для ветра
    double boxsize = 3.5;
    auto resolution_mag = 20;
	auto universe_mag = mshr::Box(Point(-boxsize, -boxsize, -boxsize), Point(boxsize, boxsize, boxsize));
	auto atmosphere_mag = universe_mag - planet;
	auto mesh_mag = mshr::generate_mesh(atmosphere_mag, resolution_mag);

    double dt = 0.01;
    double T = 5;

    //Инициализация модуля атмосферы
    atm::init(mesh_atm, dt, rad, height, G, K, rho_0, p_0);
    //Инициализация солнечного ветра и магнитного поля
    mag::init(mesh_mag, dt, boxsize);

    double t = dt;

    while (t < T + DOLFIN_EPS)
    {   
        atm::calc();
        mag::calc();
        t += dt;
        cout << "t = " << t << endl;
    }

    return 0;
}