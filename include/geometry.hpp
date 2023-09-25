#ifndef geometryheaderfile
#define geometryheaderfile

#include "type.h"
#include <tuple>
#include <iostream>

using namespace std;

struct Vec3d{
    db x, y, z;
    Vec3d() : x(0), y(0), z(0)
        {  }
    Vec3d(db _x, db _y, db _z) : x(_x), y(_y), z(_z)
        {  }

    Vec3d operator + (const Vec3d & rhs) const {
        return Vec3d(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    Vec3d operator - (const Vec3d & rhs) const {
        return Vec3d(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    db operator * (const Vec3d & rhs) const {
        return (x * rhs.x + y * rhs.y + z * rhs.z);
    }

    void unitized() {
        db len = sqrt(x * x + y * y + z * z);
        x /= len;
        y /= len;
        z /= len;
    }

    db length() const { return sqrt(x * x + y * y + z * z); }
};

Vec3d cross(const Vec3d & lhs, const Vec3d & rhs) {
    return Vec3d(lhs.y * rhs.z - lhs.z * rhs.y,
                 lhs.z * rhs.x - lhs.x * rhs.z,
                 lhs.x * rhs.y - lhs.y * rhs.x);
}

Vec3d operator * (db val, const Vec3d & rhs) {
    return Vec3d(rhs.x * val, rhs.y * val, rhs.z * val);
}
Vec3d operator * (const Vec3d & rhs, db val) { return val * rhs; }

ostream & operator << (ostream & os, const Vec3d & now) {
    return os << "(" << now.x << ", " << now.y << ", " << now.z << ")";
}

istream & operator >> (istream & is, Vec3d & now) {
    return is >> now.x >> now.y >> now.z;
}

struct Matrix33{
    Vec3d m[3];
    Matrix33() {}
    Matrix33(Vec3d a, Vec3d b, Vec3d c) {
        m[0] = a;
        m[1] = b;
        m[2] = c;
    }
    Vec3d & operator[](const int & idx) { return m[idx]; }
    Matrix33 T() const {
        return Matrix33(Vec3d(m[0].x, m[1].x, m[2].x),
                        Vec3d(m[0].y, m[1].y, m[2].y),
                        Vec3d(m[0].z, m[1].z, m[2].z));
    }
};

Vec3d operator * (Matrix33 mat, Vec3d vec) {
    return (vec.x * mat[0] + vec.y * mat[1] + vec.z * mat[2]);
}

#endif
