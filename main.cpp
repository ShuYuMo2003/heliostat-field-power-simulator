#include <iostream>
#include "include/type.h"
#include "include/datetime2vec.hpp"
#include "include/heliostat.hpp"

using namespace std;

void readNaddHels(Field & field){
    FILE * fptr = fopen("data/clear_mcm_data_total.csv", "r");
    assert(fptr != NULL);

    double x, y, z, w, h;
    while(fscanf(fptr, "%lf,%lf,%lf,%lf,%lf", &x, &y, &z, &w, &h) == 5){
        Vec3d pow(x, y, z);
        field.addHel(pow, w, h);
    }
    fclose(fptr);
}

int main(){
    Field field(Vec3d(0, 0, 80), 8, 7);

    readNaddHels(field);
    cerr << "read done" << endl;

    Vec3d light = datetime2vector(1, 9); // 1 月 9 点。
    field.calcWithMontecarlo(light);


    return 0;
}