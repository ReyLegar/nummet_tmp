#include "MethodHeat.h"
#include <cstdio>
#include <cstring>
#include <cmath>


void MethodHeat::convertToParam(int i, Param& p)
{
    p.r = 0.0;
    p.p = 0.0;
    p.u = 0.0;
    p.v = 0.0;
    p.e = 0.0;
    p.T = T[i];
    p.M = 0.0;
}


void MethodHeat::init()
{
    mesh = new Mesh();
    mesh->initFromFiles((char*)"heat.1");

    T = new double[mesh->cCount];
    memset(T, 0, sizeof(double) * mesh->cCount);

    saveVTK(0);

    intT = new double[mesh->cCount];

    TMAX = 1.0;
    TAU = 1.e-5;
}


void MethodHeat::run()
{
    double t = 0.0;
    int step = 0;
    while (t < TMAX) {
        t += TAU;
        step++;

        memset(intT, 0, sizeof(double) * mesh->cCount);
        for (int ie = 0; ie < mesh->eCount; ie++) {
            Edge& e = mesh->edges[ie];
            int c1 = e.c1;
            Cell& cell1 = mesh->getCell(c1);

            Param p;
            convertToParam(c1, p);
            double T1 = p.T;
            double T2;
            double hij;
            if (e.c2 >= 0) {
                int c2 = e.c2;
                Cell& cell2 = mesh->getCell(c2);
                convertToParam(c2, p);
                T2 = p.T;
                Vector v2 = cell2.c;
                v2 -= cell1.c;
                hij = sqrt(v2.x * v2.x + v2.y * v2.y);
            }
            else {
                switch (e.type) {
                case 1:
                    T2 = 0.0;
                    break;
                case 2:
                    T2 = 1.0;
                    break;
                }
                Vector v2 = e.c[0];
                v2 -= cell1.c;
                hij = 2.0 * sqrt(v2.x * v2.x + v2.y * v2.y);
            }

            double flux = (T2 - T1) * e.l / hij;

            intT[c1] += flux;
            if (e.c2 > -1) {
                intT[e.c2] -= flux;
            }
        }

        for (int i = 0; i < mesh->cCount; i++) {
            T[i] += intT[i] * TAU / mesh->cells[i].S;
        }

        if (step % 10000 == 0)
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

MethodHeat::~MethodHeat()
{
    delete mesh;
    delete[] T, intT;
}