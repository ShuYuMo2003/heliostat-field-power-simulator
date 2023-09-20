#include <iostream>
#include "../include/datetime2vec.hpp"

using namespace std;

int main(){
    for(int mm = 1; mm <= 12; mm++) {
        for(int hh = 9; hh <= 15; hh += 3) {
            cout << "(" << mm << ", 21, " << hh << ", 0, 0) d = " << endl;
            sun_angle_t res = datetime2angle(mm, hh);
            cout << res.sin_alpha_s << ", " << res.cos_alpha_s << ", ";
            cout << res.sin_gamma_s << ", " << res.cos_gamma_s << endl;
        }
    }
    return 0;
}
