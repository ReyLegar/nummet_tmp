#pragma once

#include "method.h"


class MethodGas : public Method
{
public:
    virtual void convertToParam(int, Param&);
    virtual void init(); //начальные условия
    virtual void run(); // запуск кода

    ~MethodGas();
protected:

    void initValues(); //вызывается в init()

    /*
    pl - параметры с лева
    pr - параметры с права 

    На выходе заполняются переменные fr, fu, fv, fe
                fr = ro * u
    F1 =        fu = ro * u^2 + p
                fv = ro * u * v
                fe = u * (ro * (eps + (u^2 + v^2)/2))


                fr = ro * v
    F2 =        fu = ro * u * v
                fv = ro * v^2 + p
                fe = v * (ro * (eps + (u^2 + v^2)/2))

    Эти векторы сразу умножены на нормаль

    */

    void flux(Param pl, Param pr, Vector n, double& fr, double& fu, double& fv, double& fe); //расчет потоков


    void bnd(Edge* e, Param p1, Param& p2); // используется когда граничное ребро (нет соседнего)

    double* ro, * ru, * rv, * re;
    double* int_ro, * int_ru, * int_rv, * int_re;

    double* grad_u_x, *grad_u_y, *grad_v_x, *grad_v_y; //-------------------------

    double TMAX; // время до которого считаем
    double TAU; // шаг по времени

    const double GAM = 1.4;
};


