#pragma once

#include <dolfin.h>

#include "Constants.hpp"

using namespace dolfin;


// Скорость у поверхности
class SurfaceVelocity : public Expression
{
public:
    SurfaceVelocity();
    void eval(Array<double> &values, const Array<double> &coord) const;
};


// Объемные силы притяжения
class GravityForces : public Expression
{
public:
    GravityForces();
    void eval(Array<double> &values, const Array<double> &coord) const;
};


// Плотность атмосферы
class Density : public Expression
{
public:
    void eval(Array<double> &values, const Array<double> &coord) const;
};


// Скорость частиц на входе
class InflowVelocity : public Expression
{
public:
	InflowVelocity();
	void eval(Array<double> &values, const Array<double> &x) const;
};


// Концентрация частиц на входе
class InflowConcentration : public Expression
{
public:
	void eval(Array<double> &values, const Array<double> &x) const;
};


// Магнитное поле Земли
class MagneticField : public Expression
{
public:
	MagneticField();
	void eval(Array<double> &values, const Array<double> &coords) const;
};
