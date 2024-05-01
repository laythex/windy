#include "Domains.hpp"

bool SurfaceBoundary::inside(const Array<double> &coord, bool on_boundary) const
{
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    return sqrt(x * x + y * y + z * z) < Constants::PLANET_RADIUS + DOLFIN_EPS;
}

bool SpaceBoundary::inside(const Array<double> &coord, bool on_boundary) const
{
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    return sqrt(x * x + y * y + z * z) > Constants::PLANET_RADIUS + Constants::ATMO_HEIGHT - Constants::PLANET_RADIUS / 5; // Вот этот костыль надо пофиксить
}                                                                                                                          // Если только 5 это не мировая константа

bool InflowBoundary::inside(const Array<double> &coord, bool on_boundary) const {
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    return x < 0 && (sqrt(x * x + y * y + z * z) > Constants::PLANET_RADIUS + Constants::ATMO_HEIGHT - Constants::PLANET_RADIUS / 5);
}
