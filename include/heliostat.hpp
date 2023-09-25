#ifndef heliostatheaderfile
#define heliostatheaderfile

#include "type.h"
#include "geometry.hpp"
#include <cassert>
#include <set>
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
    db dhr;             // 与中心距离。

    vector<int> effectSet; // 与之有影响的镜面集合。

    // @see: 塔式光热电站光学效率建模仿真及定日镜场优化布置_刘建兴
    db sigma_tot;       // 能流分布参数。

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
        this->dhr = _pos.length();
        assert(this->dhr <= 1000);
        this->eta_at = 0.99321 - 0.0001176 * this->dhr + 1.97 * (1e-8) * this->dhr * this->dhr;
    }
    void calc_sigma_tot() {
        db cos_omega = this->eta_cos;
        db f = 1e11; // TODO: make confirmation.
        db Ht = sqrt(this->width * this->height) * fabs(this->dhr / f - cos_omega);
        db Ws = sqrt(this->width * this->height) * fabs(this->dhr / f * cos_omega - 1);
        db sigma_ast = sqrt(0.5 * (Ht * Ht + Ws * Ws)) / (4 * dhr);
        db sigma_sun = 2.51e-3;
        db sigma_s   = 1.53e-3;
        db sigma_track = 1.53e-3;
        db sigma_bq = (2 * sigma_s);
        this->sigma_tot = sqrt(dhr * dhr * (sigma_sun * sigma_sun
                                          + sigma_bq * sigma_bq
                                          + sigma_ast * sigma_ast
                                          + sigma_track * sigma_track));
    }
    Matrix33 init_transform_matrix() {
        Vec3d z(this->normal); z.unitized();
        Vec3d x(1, -z.x / z.y, 0); x.unitized();
        Vec3d y = cross(z, x); y.unitized();
        assert(x * y < EPS);
        assert(x * z < EPS);
        assert(y * z < EPS);
        Matrix33 mat(x, y, z);
        this->trans2ground = mat;
        return mat;
    }
    void init_normal_with_light(const Vec3d & light, const Vec3d & center){
        Vec3d pos_inv_vector = center - this->pos;
        pos_inv_vector.unitized();
        this->normal = pos_inv_vector + (light * (-1));
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
        for(auto idx : hel.effectSet){
            Heliostat & shadingHel = Hels[idx];
            Vec3d ptInSadMirror0 = shadingHel.pt_ground2mirror(ptInGround0);
            Vec3d lightInSadMirror = shadingHel.ve_ground2mirror(light);
            if(ptInSadMirror0.z < 0) {
                db x = -(ptInSadMirror0.z / lightInSadMirror.z) * lightInSadMirror.x + ptInSadMirror0.x;
                db y = -(ptInSadMirror0.z / lightInSadMirror.z) * lightInSadMirror.y + ptInSadMirror0.y;
                if(fabs(x) <= shadingHel.width / 2 && fabs(y) <= shadingHel.height / 2) {
                    ok = false;
                    break;
                }
            }
        }
        if(not ok) {
            return false;
        }

        // 检查反射光线是否被遮挡。
        ok = true;
        ptInGround0 = hel.pt_mirror2ground(ptInHel);
        for(auto idx : hel.effectSet){
            Heliostat & shadingHel = Hels[idx];
            Vec3d ptInSadMirror0 = shadingHel.pt_ground2mirror(ptInGround0);
            Vec3d reflectInSadMirror = shadingHel.ve_ground2mirror(reflect);
            if(ptInSadMirror0.z < 0) {
                db x = -(ptInSadMirror0.z / reflectInSadMirror.z) * reflectInSadMirror.x + ptInSadMirror0.x;
                db y = -(ptInSadMirror0.z / reflectInSadMirror.z) * reflectInSadMirror.y + ptInSadMirror0.y;
                if(fabs(x) <= shadingHel.width / 2 && fabs(y) <= shadingHel.height / 2) {
                    // cerr << "Found" << x <<  " " << y << endl;
                    ok = false;
                    break;
                }
            }
        }
        return ok;
    }

    db checkPtOnTower(const Heliostat & hel,
                        const Vec3d & ptInHel,
                        const Vec3d & reflect,
                        const db & delta_with) {
        Vec3d ptInGround0 = hel.pt_mirror2ground(ptInHel);
        db k = - reflect.x / reflect.y;
        db x = (reflect.y * ptInGround0.x - reflect.x * ptInGround0.y) /                    \
                           (reflect.y - reflect.x * k);
        db t = (x - ptInGround0.x) / reflect.x;
        db y = t * reflect.y + ptInGround0.y;
        db z = t * reflect.z + ptInGround0.z;
        if((x * x + y * y < center_diameter * center_diameter / 4)                          \
                        && z > center.z - center_height / 2                                 \
                        && z < center.z + center_height / 2) {
            return 1;

            Vec3d ptOnTower(x, y, z); // 在塔集热面上的坐标。

            db distance2 = ptInHel.length(); distance2 = distance2 * distance2; // TODO: check change to (center - ptOnTower).length()

            db sigma_tot2 = hel.sigma_tot * hel.sigma_tot;
            db eta_int = (1.0 / (2 * PI * sigma_tot2)) * exp(-distance2 / (2 * sigma_tot2));
            // cerr << distance2 << " " << eta_int << endl;
            // cerr << hel.sigma_tot << " ";
            // cerr << ptOnTower << "(" << (ptOnTower - this->center).length() << ") --> " << eta_int << " area = " << (delta_with * delta_with) << endl;
            eta_int = eta_int * (delta_with * delta_with);
            return eta_int;
        } else return 0;
    }

    void initEffectSet() {
        for(int i = 0; i < Hels.size(); i++) {
            Hels[i].effectSet.clear();
            for(int j = 0; j < Hels.size(); j++) {
                if(i == j) continue;
                if((Hels[i].pos - Hels[j].pos).length() <= 70) {
                    Hels[i].effectSet.push_back(j);
                }
            }
        }
    }

    void calcWithMontecarlo(Vec3d light, db s_delta = 0.2) {
        cerr << "start montecarlo" << endl;
        initEffectSet();
        light.unitized();
        for(auto & hel : Hels) {
            cerr << "now checking " << hel.pos << endl;

            // 计算单个镜面的数据。
            hel.init_normal_with_light(light, center);
            hel.init_transform_matrix();
            hel.calc_eta_cos(light);
            hel.calc_sigma_tot();

            Vec3d reflect = center - hel.pos; reflect.unitized();
            pair<double, double> eta_sb = std::make_pair(0, 0);
            db eta_trunc = 0;
            db tot = 0;
            for(db x = -hel.width / 2; x <= hel.width / 2; x += s_delta) {
                for(db y = -hel.height / 2; y <= hel.height / 2; y += s_delta) {
                    Vec3d hel_pt(x, y, 0);
                    if(checkPtNotShade(hel, hel_pt, light, reflect)) {
                        eta_sb.first += 1;
                        eta_trunc += checkPtOnTower(hel, hel_pt, reflect, s_delta);
                        tot++;
                    }
                    eta_sb.second  += 1;
                }
            }

            hel.eta_sb = eta_sb.first / eta_sb.second;
            hel.eta_trunc = eta_trunc / tot;
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
