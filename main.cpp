#include <iostream>
#include "method.h"
#define _CRT_SECURE_NO_WARNINGS

using namespace std;

void printUsage(char * exeName)
{
    cout << "Usage:" << endl;
    cout << "  " << exeName << " <method_code>" << endl << endl;
    cout << "  <method_code>: " << METHOD_CODE_HEAT << " - heat eq. solver; " << endl <<
            "                 " << METHOD_CODE_GAS << " - gas dynamics eq. solver." << endl;
}

int main(int argc, char** argv)
{
    /*if (argc < 2) {
        printUsage(argv[0]);
        return -1;
    }
    int mCode = atoi(argv[1]);*/
    Method * m = Method::create(2);
    if (!m) {
        return -2;
    }
    m->init();
    m->run();

    delete m;
}
