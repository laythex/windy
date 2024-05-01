#pragma once

#include <array>

namespace Constants { // Все размерные величины заданы в СИ

    // Математика
    constexpr double deg2rad = M_PI / 180;

    // Временные параметры симуляции
    constexpr double DELTA_TIME                     = 1e-2;                 // Шаг по времени
    constexpr double SIM_DURATION                   = 1;                    // Длительность симуляции 

    // Параметры физической модели                                      
    constexpr double GRAV_PARAMETER                 = 0.05;                 // Гравитационный параметр (GM)                         [1]
    constexpr double PLANET_RADIUS                  = 1;                  // Радиус планеты                                       [1]
    constexpr double ATMO_HEIGHT                    = 3;                    // Высота атмосферы
    constexpr double DENSITY_ASL                    = 1;                    // Плотность на поверхности планеты
    constexpr double PRESSURE_ASL                   = 0.1;                  // Атмосферное давление
    constexpr double DISTRIB_COEFF                  = 5;                    // Коэффициент в распределениях плотности и давления
    constexpr std::array<double, 3> ANGULAR_VEL     = {0, 0, 1};            // Вектор угловой скорости планеты                      [1]
    constexpr double INC_TO_ORBITAL_PLANE           = 25.19 * deg2rad;      // Наклон оси вращения Марса к плоскости его орбиты     [1]
    constexpr double MAGNETIC_MOMENT                = 8.22e0;              // Магнитный дипольный момент Марса                     [2]
    constexpr double SOLAR_WIND_VELOCITY            = 8e1;                  // Скорость частиц солнечного ветра                     [3]

    // Параметры вычислительной сетки   
    constexpr unsigned MESH_RESOLUTION              = 20;                   // Разрешение сетки

}

/*
Ссылки
    [1]:
        https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
    [2]:
        https://www.britannica.com/science/geomagnetic-field/Dipolar-field
        https://www.sciencedirect.com/science/article/pii/S0032063397001864
    [3]:
        https://solarscience.msfc.nasa.gov/SolarWind.shtml

*/