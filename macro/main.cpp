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

        Vector atmo_c = *((atmo_pair.first)->vector());
        Vector magn_c = *((magn_pair.first)->vector());
        Vector atmo_v = *((atmo_pair.second)->vector());
        Vector magn_v = *((magn_pair.second)->vector());

        //Для нормировки
        atmo_c/=Constants::PRESSURE_ASL;

        // info(std::to_string(atmo_c.size()));
        // info(std::to_string(magn_c.size()));
        // info(std::to_string(atmo_v.size()));
        // info(std::to_string(magn_v.size()));

        Vector result_v = Vector(atmo_v);


        for (int i = 0; i < atmo_v.size()/3; i++){
            double data[] = {atmo_v[i*3] * atmo_c[i] + magn_v[i*3] * magn_c[i], atmo_v[i*3 + 1] * atmo_c[i] + magn_v[i*3 + 1] * magn_c[i + 1],
            atmo_v[i*3 + 2] * atmo_c[i] + magn_v[i*3 + 2] * magn_c[i + 2]};
            int ind[] = {i*3, i*3+1, i*3+2};
            result_v.set(data, 3, ind);
        }
        
        //*(atmo_pair.second)->vector() = result_v;
        //*(atmo_pair.second)->vector() = *(atmo_v + magn_v);
        
        
        magn_cfile << *magn_pair.first;
        magn_vfile << *magn_pair.second;
        atmo_pfile << *atmo_pair.first;
        atmo_vfile << *atmo_pair.second;

        t += Constants::DELTA_TIME;
        std::cout << "t = " << t << '\n';
    }

    return 0;
}
