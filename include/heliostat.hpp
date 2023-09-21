#ifndef heliostatheaderfile
#define heliostatheaderfile

#include "type.h"
#include "geometry.hpp"
#include <cassert>
#include <vector>

struct Heliostat{
    Vec3d pos;          // 中心位置。
    Vec3d normal;       // 法向量。
    db width, height;   // 长宽。
    db eta_cos;         // 余弦效率
    db eta_at;          // 大气效率
    db eta_ref;         // 镜面反射率
    db eta_sb;          // 阴影遮挡效率
    db eta_trunc;       // 截断效率
    db eta;             // 总效率
    Matrix33 trans2ground;       // 镜面坐标系 -> 地面坐标系转换矩阵。
    Heliostat() {  }
    Heliostat(Vec3d pp, db ww, db hh) : pos(pp), width(ww), height(hh) { }
    void calc_eta_cos(Vec3d light) {
        light = light * (-1);
        light.unitized();
        assert(fabs(light.length() - 1) < EPS);
        assert(fabs(normal.length() - 1) < EPS);
        this->eta_cos = normal * light;
    }
    void calc_eta_at(Vec3d center) {
        Vec3d _pos = center - this->pos;
        db dd = _pos.length();
        assert(dd <= 1000);
        this->eta_at = 0.99321 - 0.0001176 * dd + 1.97 * (1e-8) * dd * dd;
    }
    Matrix33 init_transform_matrix() {
        Vec3d z(this->normal); z.unitized();
        Vec3d x(-1, -z.x / z.y, 0); x.unitized();
        Vec3d y = cross(z, x); y.unitized();
        Matrix33 mat(x, y, z);
        this->trans2ground = mat;
        return mat;
    }
    void init_normal_with_light(const Vec3d & light, const Vec3d & center){
        Vec3d pos_inv_vector = center - this->pos;
        pos_inv_vector.unitized();
        this->normal = pos_inv_vector + light;
        this->normal.unitized();
    }
    void calc_eta() {
        eta = eta_at * eta_cos * eta_ref * eta_sb * eta_trunc;
    }
    Vec3d pt_mirror2ground(const Vec3d & pt) const {
        return this->trans2ground * pt + this->pos;
    }
    Vec3d ve_mirror2ground(const Vec3d & vec) const {
        return this->trans2ground * vec;
    }
    Vec3d pt_ground2mirror(const Vec3d & pt) const {
        return this->trans2ground.T() * (pt - this->pos);
    }
    Vec3d ve_ground2mirror(const Vec3d & vec) const {
        return this->trans2ground.T() * vec;
    }
};

struct Field{
    Vec3d center;
    db center_height, center_diameter;
    Field(Vec3d c, db h, db dd) : center(c), center_height(h), center_diameter(dd) {}
    vector<Heliostat> Hels;

    void addHel(Vec3d p, db w, db h, db eta_ref=0.92) {
        Heliostat tmp;
        tmp.pos = p;
        tmp.width = w;
        tmp.height = h;
        tmp.calc_eta_at(center);
        tmp.eta_ref = eta_ref;
        Hels.push_back(tmp);
    }

    bool checkPtNotShade(const Heliostat & hel,
                         const Vec3d & ptInHel,
                         const Vec3d & light,
                         const Vec3d & reflect) {
        bool ok = true;
        // 检查入射光线是否被遮挡。
        ok = true;
        Vec3d ptInGround0 = hel.pt_mirror2ground(ptInHel);
        for(auto & shadingHel : Hels) {
            if(&hel == &shadingHel) continue;
            Vec3d ptInSadMirror0 = shadingHel.pt_ground2mirror(ptInGround0);
            Vec3d lightInSadMirror = shadingHel.ve_ground2mirror(light);
            if(ptInSadMirror0.z < 0) {
                db x = -(ptInSadMirror0.z / lightInSadMirror.z) * lightInSadMirror.x + ptInSadMirror0.x;
                db y = -(ptInSadMirror0.z / lightInSadMirror.z) * lightInSadMirror.y + ptInSadMirror0.y;
                if(fabs(x) > shadingHel.width / 2 || fabs(y) > shadingHel.height / 2) {
                    ok = false;
                    break;
                }
            }
        }
        if(not ok) return false;

        // 检查反射光线是否被遮挡。
        ok = false;
        ptInGround0 = hel.pt_mirror2ground(ptInHel);
        for(auto & shadingHel : Hels){
            if(&hel == &shadingHel) continue;
            Vec3d ptInSadMirror0 = shadingHel.pt_ground2mirror(ptInGround0);
            Vec3d reflectInSadMirror = shadingHel.ve_ground2mirror(reflect);
            if(ptInSadMirror0.z < 0) {
                db x = -(ptInSadMirror0.z / reflectInSadMirror.z) * reflectInSadMirror.x + ptInSadMirror0.x;
                db y = -(ptInSadMirror0.z / reflectInSadMirror.z) * reflectInSadMirror.y + ptInSadMirror0.y;
                if(fabs(x) < shadingHel.width / 2 && fabs(y) < shadingHel.height / 2) {
                    ok = true;
                    break;
                }
            }
        }

        return ok;
    }

    db checkPtOnTower(const Heliostat & hel,
                        const Vec3d & ptInHel,
                        const Vec3d & reflect) {
        Vec3d ptInGround0 = hel.pt_mirror2ground(ptInHel);
        db k = - reflect.x / reflect.y;
        db x = (reflect.y * ptInGround0.x - reflect.x * ptInGround0.y) /                    \
                           (reflect.y - reflect.x * k);
        db t = (x - ptInGround0.x) / reflect.x;
        db y = t * reflect.y + ptInGround0.y;
        db z = t * reflect.z + ptInGround0.z;
        // (x, y, z) 在集热面上的坐标。
        return int((x * x + y * y < center_diameter * center_diameter / 4)                  \
                        && z > center.z - center_height / 2                                 \
                        && z < center.z + center_height / 2);
    }

    void calcWithMontecarlo(Vec3d light, db s_delta = 1.) {
        light.unitized();
        for(auto & hel : Hels) {
            // 计算单个镜面的数据。
            hel.init_normal_with_light(light, center);
            hel.init_transform_matrix();
            hel.calc_eta_cos(light);

            Vec3d reflect = center - hel.pos; reflect.unitized();
            pair<double, double> eta_sb = std::make_pair(0, 0), eta_trunc = std::make_pair(0, 0);
            for(db x = -hel.width / 2; x < hel.width / 2; x += s_delta) {
                for(db y = -hel.height / 2; y < hel.height / 2; y += s_delta) {
                    Vec3d hel_pt(x, y, 0);
                    if(checkPtNotShade(hel, hel_pt, light, reflect)) {
                        eta_sb.first += 1;
                        eta_trunc.first += checkPtOnTower(hel, hel_pt, reflect);
                        eta_trunc.second += 1;
                    }
                    eta_sb.second  += 1;
                }
            }

            hel.eta_sb = eta_sb.first / eta_sb.second;
            hel.eta_trunc = eta_trunc.first / eta_trunc.second;

            hel.calc_eta();

        }
        FILE * fptr = fopen("data/eta.csv", "w");
        fprintf(fptr, "x,y,z,eta_at,eta_cos,eta_ref,eta_sb,eta_trunc,eta\n");
        for(auto & hel : Hels) {
            fprintf(fptr, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                hel.pos.x,
                hel.pos.y,
                hel.pos.z,
                hel.eta_at,
                hel.eta_cos,
                hel.eta_ref,
                hel.eta_sb,
                hel.eta_trunc,
                hel.eta);
        }
        fclose(fptr);
    }
};


#endif
