#pragma once
#include "mesh.h"

#define METHOD_CODE_HEAT 1
#define METHOD_CODE_GAS  2

struct Param
{
    double r;		//!< плотность
    double p;		//!< давление
    double u;		//!< первая  компонента  вектора  скорости
    double v;		//!< вторая  компонента  вектора  скорости
    double e;		//!< внутренняя   энергия
    double T;		//!< температура
    double M;		//!< число  Маха
    double mu;   //!< коэффициент вязкости

    inline double magU() { return u * u + v * v; };
};

class Method
{
protected:
    Mesh* mesh;
public:
    virtual void convertToParam(int, Param&) = 0;
    void saveVTK(int step);
    virtual void init() = 0;
    virtual void run() = 0;

    static Method* create(int methodCode);

};
