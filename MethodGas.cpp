#include "MethodGas.h"
#include <cstdio>
#include <cstring>
#include <cmath>

const double mu = 0.1;

void MethodGas::convertToParam(int i, Param& p)
{
    p.r = ro[i]; // ro в U
    p.u = ru[i] / p.r; // в U ro * u / ro
    p.v = rv[i] / p.r;
    p.e = re[i] / p.r - p.magU() / 2.; // по формуле в F2 вычисляем внутреннюю энергию eps
    p.p = p.r * p.e * (GAM - 1.0); // по формуле вычисляется давление
    p.T = 0.0; // не используется 
    p.M = p.magU() / sqrt(GAM * p.p / p.r); // число Маха = скорость / скорость звука 
    p.mu = 0.1;
    p.k = 0.1;
}


void MethodGas::init()
{
    // ----- Инициализируется и считывается сетка ---
    mesh = new Mesh();
    mesh->initFromFiles((char*)"carman.1");
    // ----------------------------------------------

    // ----- Создаются массивы, размерностью равной количеству ячеек
    ro = new double[mesh->cCount];
    ru = new double[mesh->cCount];
    rv = new double[mesh->cCount];
    re = new double[mesh->cCount];
    // --------------------------------------

    initValues();
    saveVTK(0);

    // ------ Инициализируются массивы для интегралов -----
    int_ro = new double[mesh->cCount];
    int_ru = new double[mesh->cCount];
    int_rv = new double[mesh->cCount];
    int_re = new double[mesh->cCount];
    // ---------------------------------------

    grad_u_x = new double[mesh->eCount];
    grad_u_y = new double[mesh->eCount];
    grad_v_x = new double[mesh->eCount];
    grad_v_y = new double[mesh->eCount];



    TMAX = 4.0;
    TAU = 1.e-4;
}


void MethodGas::initValues()
{
    for (int i = 0; i < mesh->cCount; i++) {
        ro[i] = 1.;
        ru[i] = 3.;
        rv[i] = 0.0;
        re[i] = (1. / 1.4) / (GAM - 1.0) + 0.5 * (ru[i] * ru[i]) / ro[i];
    }
}


void MethodGas::run()
{
    double t = 0.0;
    int step = 0;
    while (t < TMAX) {
        t += TAU;
        step++;
        // ---- Зануляем интегралы в каждой ячейке -----
        memset(int_ro, 0, sizeof(double) * mesh->cCount);
        memset(int_ru, 0, sizeof(double) * mesh->cCount);
        memset(int_rv, 0, sizeof(double) * mesh->cCount);
        memset(int_re, 0, sizeof(double) * mesh->cCount);
        // ----------------------

        memset(grad_u_x, 0, sizeof(double) * mesh->cCount);
        memset(grad_u_y, 0, sizeof(double) * mesh->cCount);
        memset(grad_v_x, 0, sizeof(double) * mesh->cCount);
        memset(grad_v_y, 0, sizeof(double) * mesh->cCount);


        /*
            Надо пройтись по ребрам и найти градиенты
            1. Заводим отдельные массивы grad_u, grad_v, или grad_u_x, grad_u_y, grad_v_x, grad_v_y       (+)
            2. В таком же цикле по ребрам считаем по формуле Грина-Гаусса                                 (+-)
            3. Перед тем как начнется цикл основной у нас уже есть массив градиентов                      (+-)
            4. Надо завести переменные соответствующие вязким потокам аналогично flux(p1, p2, e.n, fr, fu, fv, fe) (здесь fr, fu, fv, fe)
            5. Подсчитать переменные аналогично
                fr = 0.5 * ((frl + frr) - alpha * (ror - rol));
                fu = 0.5 * ((ful + fur) - alpha * (rur - rul));
                fv = 0.5 * ((fvl + fvr) - alpha * (rvr - rvl));
                fe = 0.5 * ((fel + fer) - alpha * (rer - rel));
                Только здесь изменится только первая часть, до alpha
                Полусуммы градиентов скоростей
        */


        double* sum_grad_u_x = new double[mesh->cCount];
        double* sum_grad_u_y = new double[mesh->cCount];
        double* sum_grad_v_x = new double[mesh->cCount];
        double* sum_grad_v_y = new double[mesh->cCount];

        memset(sum_grad_u_x, 0, sizeof(double) * mesh->cCount);
        memset(sum_grad_u_y, 0, sizeof(double) * mesh->cCount);
        memset(sum_grad_v_x, 0, sizeof(double) * mesh->cCount);
        memset(sum_grad_v_y, 0, sizeof(double) * mesh->cCount);

        for (int ie = 0; ie < mesh->eCount; ie++) {
            Edge& e = mesh->edges[ie];
            int c1 = e.c1;
            Cell& cell1 = mesh->getCell(c1);

            Param p1, p2;
            convertToParam(c1, p1);
            if (e.c2 >= 0) {
                convertToParam(e.c2, p2);
            }
            else {
                bnd(&e, p1, p2);
            }

            sum_grad_u_x[c1] += (p1.u + p2.u) / 2 * e.l * e.n.x;
            sum_grad_u_y[c1] += (p1.u + p2.u) / 2 * e.l * e.n.y;
            sum_grad_v_x[c1] += (p1.v + p2.v) / 2 * e.l * e.n.x;
            sum_grad_v_y[c1] += (p1.v + p2.v) / 2 * e.l * e.n.y;

            if (e.c2 >= 0) {
                int c2 = e.c2;
                sum_grad_u_x[c2] -= (p1.u + p2.u) / 2 * e.l * e.n.x;
                sum_grad_u_y[c2] -= (p1.u + p2.u) / 2 * e.l * e.n.y;
                sum_grad_v_x[c2] -= (p1.v + p2.v) / 2 * e.l * e.n.x;
                sum_grad_v_y[c2] -= (p1.v + p2.v) / 2 * e.l * e.n.y;
            }
        }

        for (int i = 0; i < mesh->cCount; i++) {
            Cell& cell = mesh->cells[i];
            grad_u_x[i] = sum_grad_u_x[i] / cell.S;
            grad_u_y[i] = sum_grad_u_y[i] / cell.S;
            grad_v_x[i] = sum_grad_v_x[i] / cell.S;
            grad_v_y[i] = sum_grad_v_y[i] / cell.S;
        }

        delete[] sum_grad_u_x;
        delete[] sum_grad_u_y;
        delete[] sum_grad_v_x;
        delete[] sum_grad_v_y;



        for (int ie = 0; ie < mesh->eCount; ie++) {
            Edge& e = mesh->edges[ie];
            int c1 = e.c1;
            int c2 = e.c2;
            Cell& cell1 = mesh->getCell(c1);

            Param p1, p2;
            convertToParam(c1, p1);
            if (e.c2 >= 0) {
                convertToParam(e.c2, p2);
            }
            else {
                bnd(&e, p1, p2);
            }

            double fr, fu, fv, fe;
            flux(p1, p2, e.n, fr, fu, fv, fe);

            fr *= e.l;
            fu *= e.l;
            fv *= e.l;
            fe *= e.l;
            int_ro[c1] -= fr;
            int_ru[c1] -= fu;
            int_rv[c1] -= fv;
            int_re[c1] -= fe;
            if (e.c2 >= 0) {
                int_ro[e.c2] += fr;
                int_ru[e.c2] += fu;
                int_rv[e.c2] += fv;
                int_re[e.c2] += fe;
            }

            //printf("%lf | %lf | %lf | %lf\n", int_ro[c1], int_ru[c1], int_rv[c1], int_re[c1]);



            double tau_xx_l = 2.0 * mu * grad_u_x[c1] - (2.0 / 3.0) * mu * (grad_u_x[c1] + grad_v_y[c1]);
            double tau_yy_l = 2.0 * mu * grad_v_y[c1] - (2.0 / 3.0) * mu * (grad_u_x[c1] + grad_v_y[c1]);
            double tau_xy_l = mu * (grad_u_y[c1] + grad_v_x[c1]);

            double tau_xx_r = 0, tau_yy_r = 0, tau_xy_r = 0;

            if (e.c2 >= 0) {
                tau_xx_r = 2.0 * mu * grad_u_x[c2] - (2.0 / 3.0) * mu * (grad_u_x[c2] + grad_v_y[c2]);
                tau_yy_r = 2.0 * mu * grad_v_y[c2] - (2.0 / 3.0) * mu * (grad_u_x[c2] + grad_v_y[c2]);
                tau_xy_r = mu * (grad_u_y[c2] + grad_v_x[c2]);
            }

            double fv_fu_l = tau_xx_l * e.n.x + tau_xy_l * e.n.y;
            double fv_fv_l = tau_xy_l * e.n.x + tau_yy_l * e.n.y;
            double fv_fe_l = (tau_xx_l * p1.u + tau_xy_l * p1.v) * e.n.x + (tau_xy_l * p1.u + tau_yy_l * p1.v) * e.n.y;

            double fv_fu_r = 0;
            double fv_fv_r = 0;
            double fv_fe_r = 0;

            if (e.c2 >= 0) {
                fv_fu_r = tau_xx_r * e.n.x + tau_xy_r * e.n.y;
                fv_fv_r = tau_xy_r * e.n.x + tau_yy_r * e.n.y;
                fv_fe_r = (tau_xx_r * p2.u + tau_xy_r * p2.v) * e.n.x + (tau_xy_r * p2.u + tau_yy_r * p2.v) * e.n.y;
            }

            double fu2 = 0.5 * (fv_fu_l + fv_fu_r);
            double fv2 = 0.5 * (fv_fv_l + fv_fv_r);
            double fe2 = 0.5 * (fv_fe_l + fv_fe_r);

            fu2 *= e.l;
            fv2 *= e.l;
            fe2 *= e.l;
            int_ru[c1] += fu2;
            int_rv[c1] += fv2;
            int_re[c1] += fe2;
            if (e.c2 >= 0) {
                int_ru[e.c2] -= fu2;
                int_rv[e.c2] -= fv2;
                int_re[e.c2] -= fe2;
            }


        }

        for (int i = 0; i < mesh->cCount; i++) {
            double CFL = TAU / mesh->cells[i].S;

            ro[i] += int_ro[i] * CFL;
            ru[i] += int_ru[i] * CFL;
            rv[i] += int_rv[i] * CFL;
            re[i] += int_re[i] * CFL;
            //printf("%lf | %lf | %lf | %lf\n", ro[i], ru[i], rv[i], re[i]);

        }


        if (step % 100 == 0)
        {
            try {
                saveVTK(step);
            }
            catch (int code) {
                if (code == 1) {
                    printf("Exception have been encountered: value is not a number");
                    break;
                }
                else {
                    printf("Exception have been caught, but it's impossible to identify it.");
                }
            }
            printf("Calculation results for step %d are saved.\n", step);
        }
    }
}


void MethodGas::flux(Param pl, Param pr, Vector n, double& fr, double& fu, double& fv, double& fe)
{
    double rol, rul, rvl, rel;
    double ror, rur, rvr, rer;
    double frl, ful, fvl, fel;
    double frr, fur, fvr, fer;
    double alpha, unl, unr, q1, q2;

    unl = n.x * pl.u + n.y * pl.v;
    unr = n.x * pr.u + n.y * pr.v;

    rol = pl.r;
    rul = pl.r * pl.u;
    rvl = pl.r * pl.v;
    rel = pl.r * (pl.e + 0.5 * pl.magU());

    ror = pr.r;
    rur = pr.r * pr.u;
    rvr = pr.r * pr.v;
    rer = pr.r * (pr.e + 0.5 * pr.magU());

    frl = pl.r * unl;
    ful = frl * pl.u + pl.p * n.x;
    fvl = frl * pl.v + pl.p * n.y;
    fel = (rel + pl.p) * unl;

    frr = pr.r * unr;
    fur = frr * pr.u + pr.p * n.x;
    fvr = frr * pr.v + pr.p * n.y;
    fer = (rer + pr.p) * unr;

    q1 = sqrt(pl.p * GAM / pl.r) + fabs(pl.magU());
    q2 = sqrt(pr.p * GAM / pr.r) + fabs(pr.magU());
    alpha = (q1 > q2) ? q1 : q2;

    fr = 0.5 * ((frl + frr) - alpha * (ror - rol));
    fu = 0.5 * ((ful + fur) - alpha * (rur - rul));
    fv = 0.5 * ((fvl + fvr) - alpha * (rvr - rvl));
    fe = 0.5 * ((fel + fer) - alpha * (rer - rel));
}


void MethodGas::bnd(Edge* e, Param p1, Param& p2)
{
    switch (e->type) {
    case 1: // вытекание
        p2 = p1;
        break;
    case 2: // втекание
        p2.r = 1.;
        p2.u = 3.;
        p2.v = 0.0;
        p2.p = 1. / 1.4;
        p2.e = p2.p / p2.r / (GAM - 1.0);
        p2.T = 0.0;
        p2.M = 0.0;
        break;
    case 3: // отражение
        p2 = p1;
        double Un = p1.u * e->n.x + p1.v * e->n.y;
        Vector V;
        V.x = e->n.x * Un * 2.0;
        V.y = e->n.y * Un * 2.0;
        p2.u = p1.u - V.x;
        p2.v = p1.v - V.y;
        break;
    }
}


MethodGas::~MethodGas()
{
    delete mesh;
    delete[] ro, ru, rv, re;
    delete[] int_ro, int_ru, int_rv, int_re;
}