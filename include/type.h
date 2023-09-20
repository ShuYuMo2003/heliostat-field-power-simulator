#ifndef typeheaderfile
#define typeheaderfile
#include <cmath>
typedef double db;

const db PI = acos(-1.0);
const db EPS = 1e-8;

struct sun_angle_t {
    db sin_alpha_s, cos_alpha_s, sin_gamma_s, cos_gamma_s;
};

#endif