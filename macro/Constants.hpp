#pragma once

#include <array>

namespace Constants { // Рабочие параметры

    // Математика
    constexpr double deg2rad = M_PI / 180;
    
    // Параметры вычислительной сетки                   
    constexpr unsigned MESH_RESOLUTION              = 20;                                   // Разрешение сетки

    // Временные параметры симуляции
    constexpr double DELTA_TIME                     = 1e-2;                                 // Шаг по времени
    constexpr double SIM_DURATION                   = 1;                                    // Длительность симуляции 

    // Параметры физической модели                                                      
    constexpr double GRAV_PARAMETER                 = 0.05;                                 // Гравитационный параметр (GM)                         
    constexpr double PLANET_RADIUS                  = 0.5;                                  // Радиус планеты                                       
    constexpr double ATMO_HEIGHT                    = 3 * PLANET_RADIUS;                    // Высота атмосферы
    constexpr double DENSITY_ASL                    = 1;                                    // Плотность на поверхности планеты
    constexpr double PRESSURE_ASL                   = 0.1;                                  // Атмосферное давление
    constexpr double DISTRIB_COEFF                  = 5;                                    // Коэффициент в распределениях плотности и давления
    constexpr double SIDERIAL_PERIOD                = 24.6229 * 60 * 60;                    // Сидерический период обращения                        
    constexpr double INC_TO_ORBITAL_PLANE           = 25.19 * deg2rad;                      // Наклон оси вращения Марса к плоскости его орбиты     
    constexpr double MAGNETIC_MOMENT                = 1e-8 * 0;        // На самом деле подгон  // Магнитный дипольный момент Марса                  
    constexpr double SOLAR_WIND_VELOCITY            = 1;                                    // Скорость частиц солнечного ветра     
    constexpr std::array<double, 3> ANGULAR_VEL     = { sin(INC_TO_ORBITAL_PLANE),          // Вектор угловой скорости планеты  
                                                        cos(INC_TO_ORBITAL_PLANE), 
                                                        0 };
    constexpr double MASS_RATIO                     = 20;                                    // Отношение массы частицы атмосферы и массы частицы ветра
}

// namespace Constants_Miha{ // Все размерные величины заданы в СИ

//     // Математика
//     constexpr double deg2rad = M_PI / 180;

//     // Временные параметры симуляции
//     constexpr double DELTA_TIME                     = 1e-2;                                 // Шаг по времени
//     constexpr double SIM_DURATION                   = 1;                                    // Длительность симуляции 

//     // Параметры физической модели                                                      
//     constexpr double GRAV_PARAMETER                 = 4.2828e13;                            // Гравитационный параметр (GM)                         [1]
//     constexpr double PLANET_RADIUS                  = 3.3895e6;                             // Радиус планеты                                       [1]
//     constexpr double ATMO_HEIGHT                    = 3 * PLANET_RADIUS;                    // Высота атмосферы
//     constexpr double DENSITY_ASL                    = 4.5e-2;                                    // Плотность на поверхности планеты
//     constexpr double PRESSURE_ASL                   = 1;                                  // Атмосферное давление
//     constexpr double DISTRIB_COEFF                  = 5;                                    // Коэффициент в распределениях плотности и давления
//     constexpr double SIDERIAL_PERIOD                = 24.6229 * 60 * 60;                    // Сидерический период обращения                        [1]
//     constexpr std::array<double, 3> ANGULAR_VEL     = {0, 0, 2 * M_PI / SIDERIAL_PERIOD};   // Вектор угловой скорости планеты                      
//     constexpr double INC_TO_ORBITAL_PLANE           = 25.19 * deg2rad;                      // Наклон оси вращения Марса к плоскости его орбиты     [1]
//     constexpr double MAGNETIC_MOMENT                = 8.22e15;  // На самом деле подгон     // Магнитный дипольный момент Марса                     [2]
//     constexpr double SOLAR_WIND_VELOCITY            = 8e5;                                  // Скорость частиц солнечного ветра                     [3]

//     // Параметры вычислительной сетки                   
//     constexpr unsigned MESH_RESOLUTION              = 30;                                   // Разрешение сетки

// }

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