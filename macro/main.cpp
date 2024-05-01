#include <dolfin.h>
#include <filesystem>

#include "mshr.h"

#include "MagneticFieldModel.hpp"
#include "AtmosphereModel.hpp"
#include "Constants.hpp"

using namespace dolfin;

int main()
{   
    // Создание сетки
    auto planet_mesh = mshr::Sphere(Point(0, 0, 0), Constants::PLANET_RADIUS);
    auto universe_mesh = mshr::Sphere(Point(0, 0, 0), Constants::PLANET_RADIUS + Constants::ATMO_HEIGHT);
    auto atmosphere_mesh = universe_mesh - planet_mesh;
    auto mesh = mshr::generate_mesh(atmosphere_mesh, Constants::MESH_RESOLUTION);

    MagneticFieldModel magn(mesh);
    AtmosphereModel atmo(mesh);

    double t = Constants::DELTA_TIME;
    while (t < Constants::SIM_DURATION + DOLFIN_EPS)
    {   
        magn.calculate();
        // atmo.calculate();

        t += Constants::DELTA_TIME;
        std::cout << "t = " << t << '\n';
    }

    return 0;
}
