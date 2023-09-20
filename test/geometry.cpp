#include <iostream>
#include "../include/geometry.hpp"

using namespace std;

int main(){
    Vec3d vec(1, 2, 3), vvv(4, 5, 6);
    cout << vec << endl;
    cout << vec * 10 << endl;
    cout << vec - vvv << endl;
    cout << vec + vvv << endl;
    cout << vec * 1100 << endl;
    Matrix33 mat(
        Vec3d(1, 2, 3),
        Vec3d(4, 5, 6),
        Vec3d(7, 8, 9)
    );
    cout << mat * vec << endl;
    return 0;
}
