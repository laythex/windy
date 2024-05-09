#include "BoundaryConditions.hpp"

/*=============================================================================*/

SurfaceVelocity::SurfaceVelocity() : Expression(3) { }

void SurfaceVelocity::eval(Array<double> &values, const Array<double> &coord) const
{
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    values[0] = Constants::ANGULAR_VEL[1] * z - Constants::ANGULAR_VEL[2] * y;
    values[1] = Constants::ANGULAR_VEL[2] * x - Constants::ANGULAR_VEL[0] * z;
    values[2] = Constants::ANGULAR_VEL[0] * y - Constants::ANGULAR_VEL[1] * x;
}

/*=============================================================================*/

GravityForces::GravityForces() : Expression(3) { }

void GravityForces::eval(Array<double> &values, const Array<double> &coord) const
{
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    double r = sqrt(x * x + y * y + z * z);
    values[0] = -Constants::GRAV_PARAMETER * x / (r * r * r);
    values[1] = -Constants::GRAV_PARAMETER * y / (r * r * r);
    values[2] = -Constants::GRAV_PARAMETER * z / (r * r * r);
}

/*=============================================================================*/

void Density::eval(Array<double> &values, const Array<double> &coord) const
{
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    double r = sqrt(x * x + y * y + z * z);
    if (r > Constants::PLANET_RADIUS - DOLFIN_EPS)
    {
        values[0] = Constants::DENSITY_ASL * exp(-Constants::DISTRIB_COEFF * (r - Constants::PLANET_RADIUS));
    }
    else
    {
        values[0] = Constants::DENSITY_ASL;
    }
}

/*=============================================================================*/

InflowVelocity::InflowVelocity() : Expression(3) { }

void InflowVelocity::eval(Array<double> &values, const Array<double> &x) const
{
    values[0] = Constants::SOLAR_WIND_VELOCITY;
    values[1] = 0;
    values[2] = 0;
}

/*=============================================================================*/

void InflowConcentration::eval(Array<double> &values, const Array<double> &x) const
{
    // Относительная величина
    values[0] = 1;
}

/*=============================================================================*/

MagneticField::MagneticField() : Expression(3) {}

void MagneticField::eval(Array<double> &values, const Array<double> &coords) const
{
    // Проекции магнитного момента на оси координат
    double mx = Constants::MAGNETIC_MOMENT * Constants::ANGULAR_VEL[0];
    double my = Constants::MAGNETIC_MOMENT * Constants::ANGULAR_VEL[1];
    double mz = Constants::MAGNETIC_MOMENT * Constants::ANGULAR_VEL[2];

    double x = coords[0], y = coords[1], z = coords[2];

    // Штуки для вычисления поля диполя
    double rm = x * mx + y * my + z * mz;
    double r2 = x * x + y * y + z * z;
    double r3 = r2 * sqrt(r2);
    double rm2 = 3 * rm / r2;

    // Магнитное поле диполя
    values[0] = (x * rm2 - mx) / r3;
    values[1] = (y * rm2 - my) / r3;
    values[2] = (z * rm2 - mz) / r3;
}

/*=============================================================================*/
