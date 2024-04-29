#pragma once

#include <array>

namespace Constants {

    // Временные параметры симуляции
    constexpr double DELTA_TIME                  = 0.01;        // Шаг по времени
    constexpr double SIM_DURATION                = 5;           // Длительность
    
    // Параметры физической модели                       
    constexpr double GRAV_PARAMETER              = 0.05;        // Гравитационный параметр (GM)
    constexpr double PLANET_RADIUS               = 0.5;         // Радиус планеты
    constexpr double ATMO_HEIGHT                 = 3;           // Высота атмосферы
    constexpr double DENSITY_ASL                 = 1;           // Плотность на поверхности планеты
    constexpr double PRESSURE_ASL                = 0.1;         // Атмосферное давление
    constexpr double DISTRIB_COEFF               = 5;           // Коэффициент в распределении плотности и давления
    constexpr std::array<double, 3> ANGULAR_VEL  = {0, 0, 1};   // Вектор угловой скорости планеты

    // Параметры вычислительной сетки   
    constexpr unsigned MESH_RESOLUTION           = 15;          // Разрешение сетки

}
 