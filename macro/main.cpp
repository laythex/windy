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

    // Указатели на функции, в которые будут записаны результаты вычислений
    std::shared_ptr<Function> magn_vel;
    std::shared_ptr<Function> magn_conc;
    std::shared_ptr<Function> atmo_vel;
    std::shared_ptr<Function> atmo_pres;

    MagneticFieldModel magn(mesh, magn_vel, magn_conc);
    AtmosphereModel atmo(mesh, atmo_vel, atmo_pres);

    // Файлы для записи результатов
    File magn_vel_file("results/magn_velocity.pvd");
    File magn_conc_file("results/magn_concentration.pvd");
    File atmo_vel_file("results/atmo_velocity.pvd");
    File atmo_pres_file("results/atmo_pressure.pvd");

    double t = Constants::DELTA_TIME;
    while (t < Constants::SIM_DURATION + DOLFIN_EPS)
    {   
        // Проводим одну итерацию вычислений
        magn.calculate();
        atmo.calculate();

        Vector magn_vel_vec = *(magn_vel->vector());
        Vector magn_conc_vec = *(magn_conc->vector());
        Vector atmo_vel_vec = *(atmo_vel->vector());
        Vector atmo_pres_vec = *(atmo_pres->vector());

        // Для нормисов
        atmo_pres_vec /= Constants::PRESSURE_ASL;

        // info(std::to_string(atmo_vel_vec[0]));
        // info(std::to_string(atmo_pres_vec[0]));

        std::shared_ptr<GenericVector> result_vec = atmo_vel->vector();

        for (int i = 0; i < result_vec->size() / 3; i++){
            std::array<double, 3> total = { atmo_vel_vec[i * 3 + 0] * atmo_pres_vec[i] + magn_vel_vec[i * 3 + 0] * magn_conc_vec[i], 
                                            atmo_vel_vec[i * 3 + 1] * atmo_pres_vec[i] + magn_vel_vec[i * 3 + 1] * magn_conc_vec[i],
                                            atmo_vel_vec[i * 3 + 2] * atmo_pres_vec[i] + magn_vel_vec[i * 3 + 2] * magn_conc_vec[i] };
            std::array<int, 3> indicies = { i * 3 + 0, 
                                            i * 3 + 1, 
                                            i * 3 + 2 };

            result_vec->set(total.data(), 3, indicies.data());
        }
        
        magn_vel_file << *magn_vel;
        magn_conc_file << *magn_conc;
        atmo_vel_file << *atmo_vel;
        atmo_pres_file << *atmo_pres;

        t += Constants::DELTA_TIME;
        std::cout << "t = " << t << std::endl;
    }

    return 0;
}
