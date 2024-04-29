#pragma once

#include <dolfin.h>

#include "Constants.hpp"

using namespace dolfin;

// Граница поверхности
class SurfaceBoundary : public SubDomain
{
    bool inside(const Array<double> &coord, bool on_boundary) const;
};

// Внешняя граница атмосферы
class SpaceBoundary : public SubDomain
{
    bool inside(const Array<double> &coord, bool on_boundary) const;
};

// Подветренная сторона атмосферы
class InflowBoundary : public SubDomain
{
    bool inside(const Array<double> &coord, bool on_boundary) const;
};
