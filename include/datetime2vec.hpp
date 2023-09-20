#ifndef datetimevecheaderfile
#define datetimevecheaderfile

#include <cmath>
#include <cassert>
#include "geometry.hpp"
#include "type.h"

db cs2sc(db x) {
    db c2 = 1 - x * x;
    assert(c2 > -EPS);
    return sqrt(c2 > 0 ? c2 : 0);
}

/**
 * @brief month 月 21 日 hour 时 0 分 0 秒 的太阳高度/方位角。
*/
sun_angle_t datetime2angle(int month, int hour) {
    static int month_table[13] = {-1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30};
    static db phi = 0.6876597252857658;
    int ddd = 0;
    if(month > 3) {
        for(int mm = 3; mm < month; mm++) {
            ddd += month_table[mm];
        }
    } else if(month < 3) {
        for(int mm = month; mm < 3; mm++) {
            ddd -= month_table[mm];
        }
    }


    db D = ddd;
    db ST = hour;
    db omega = PI / 12 * (ST - 12);

    sun_angle_t res;
    db sin_delta = sin(2 * PI * D / 365) * sin(2 * PI / 360 * 23.45);
    db cos_delta = cs2sc(sin_delta);

    res.sin_alpha_s = cos_delta * cos(phi) * cos(omega) + sin_delta * sin(phi);
    res.cos_alpha_s = cs2sc(res.sin_alpha_s);

    res.cos_gamma_s = (sin_delta - res.sin_alpha_s * sin(phi)) / (res.cos_alpha_s * cos(phi));
    res.sin_gamma_s = cs2sc(res.cos_gamma_s);

    if(ST > 12) res.sin_gamma_s = -res.sin_gamma_s;
    return res;
}

Vec3d datetime2vector(int month, int hour){
    sun_angle_t ang = datetime2angle(month, hour);
    Vec3d res(
        ang.sin_gamma_s,
        ang.cos_gamma_s,
        ang.sin_alpha_s / ang.cos_alpha_s
    );
    res.unitized(); res = res * -1;
    return res;
}

#endif