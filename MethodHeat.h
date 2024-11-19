#pragma once

#include "method.h"


class MethodHeat : public Method
{
public:
    virtual void convertToParam(int, Param&);
    virtual void init();
    virtual void run();

    ~MethodHeat();
protected:

    double* T;
    double* intT;

    double TMAX;
    double TAU;
};


