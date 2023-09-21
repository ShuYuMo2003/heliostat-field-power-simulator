#include <iostream>
#include "include/type.h"
#include "include/datetime2vec.hpp"
#include "include/heliostat.hpp"

using namespace std;

void readNaddHels(Field & field){
    FILE * fptr = fopen("data/mcmdata.csv", "r");
    double x, y, z, w, h;
    while(fscanf(fptr, "%lf,%lf,%lf,%lf,%lf", &x, &y, &z, &w, &h) != EOF){
        Vec3d pow(x, y, z);
        field.addHel(pow, w, h);
    }
    fclose(fptr);
}

int main(){
    Field field(Vec3d(0, 0, 80), 8, 7);

    readNaddHels(field);

    Vec3d light = datetime2vector(1, 9); // 1 月 9 点。
    field.calcWithMontecarlo(light);


    return 0;
}