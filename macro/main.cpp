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

    //Файлы для записи результатов
    File magn_vfile{"results/magn_velocity.pvd"};
    File magn_cfile{"results/magn_concentration.pvd"};
    File atmo_vfile{"results/atmo_velocity.pvd"};
    File atmo_pfile{"results/atmo_pressure.pvd"};

    double t = Constants::DELTA_TIME;
    while (t < Constants::SIM_DURATION + DOLFIN_EPS)
    {   
        auto magn_pair = magn.calculate();
        auto atmo_pair = atmo.calculate();
        
        atmo_pair.second->vector = (atmo_pair.second)->vector + (magn_pair.second)->vector;

        magn_cfile << *magn_pair.first;
        magn_vfile << *magn_pair.second;
        atmo_pfile << *atmo_pair.first;
        atmo_vfile << *atmo_pair.second;

        t += Constants::DELTA_TIME;
        std::cout << "t = " << t << '\n';
    }

    return 0;
}
