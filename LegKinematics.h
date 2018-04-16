#ifndef LEG_KINEMATICS_H 
#define LEG_KINEMATICS_H 

#include <cmath>

namespace time_optimal
{
    namespace kinematics
    {
        static const double L_AB = 0.3;
        static const double L_BE = 0.3;
        static const double THETA0[2] {0, 0};
        static const double PI = 3.14159265358979;

        class Leg
        {
        public:
            static void LegIK(double *tip_pos_in, double *joint_angle_out, double leg_orient = 1);
            static void LegFK(double *joint_angle_in, double *tip_pos_out, double leg_orient = 1);
            static void LegIJ(double *tip_pos_in, double *jacobi_out, double leg_orient = 1);
            static void LegFJ(double *tip_pos_in, double *jacobi_out, double leg_orient = 1);
            static void LegIdJ(double *tip_pos_in, double *d_jacobi_out_x, double *d_jacobi_out_y, double leg_orient = 1);
        };
    }
}

#endif
