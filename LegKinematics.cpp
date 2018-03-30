#include "LegKinematics.h"

namespace robot_app
{
    namespace kinematics
    {
        void Leg::LegIK(double *tip_pos_in, double *joint_angle_out, double leg_orient)
        {
            using namespace std;

            double len   = sqrt(tip_pos_in[0] * tip_pos_in[0] + tip_pos_in[1] * tip_pos_in[1]);
            double alpha = atan2(tip_pos_in[1], leg_orient * tip_pos_in[0]);
            double beta  = acos((L_AB * L_AB + len * len - L_BE * L_BE) / (2 * L_AB * len));
            double gamma = acos((L_AB * L_AB + L_BE * L_BE - len * len) / (2 * L_AB * L_BE));

            joint_angle_out[0] = alpha - beta + PI/2.0 - THETA0[0];
            joint_angle_out[1] = alpha - beta - gamma + PI + PI/2.0 - THETA0[1];
        }

        void Leg::LegFK(double *joint_angle_in, double *tip_pos_out, double leg_orient)
        {
            using namespace std;
            double theta0 = joint_angle_in[0] + THETA0[0] - PI/2.0;
            double theta1 = joint_angle_in[1] + THETA0[1] - PI/2.0;

            tip_pos_out[0] = leg_orient * (L_AB * cos(theta0) + L_BE * cos(theta1));
            tip_pos_out[1] = L_AB * sin(theta0) + L_BE * sin(theta1);
        }

        void Leg::LegIJ(double *tip_pos_in, double *jacobi_out, double leg_orient)
        {
            //inverse of forward jacobi
            double joint_angle[2];
            Leg::LegIK(tip_pos_in, joint_angle, leg_orient);

            jacobi_out[0] = -sin(joint_angle[0]) / (L_BE * sin(joint_angle[0] - joint_angle[1]));
            jacobi_out[1] = -cos(joint_angle[0]) / (L_BE * sin(joint_angle[0] - joint_angle[1]));
            jacobi_out[2] = sin(joint_angle[1]) / (L_AB * sin(joint_angle[0] - joint_angle[1]));
            jacobi_out[3] = cos(joint_angle[1]) / (L_AB * sin(joint_angle[0] - joint_angle[1]));
        }

        void Leg::LegIdJ(double *tip_pos_in, double *d_jacobi_out_x, double *d_jacobi_out_y, double leg_orient)
        {
            double joint_angle[2];
            Leg::LegIK(tip_pos_in, joint_angle, leg_orient);
            double jacobi[4];
            Leg::LegIJ(tip_pos_in, jacobi, leg_orient);

            double d_jacobi_angle0 [4];
            double d_jacobi_angle1 [4];
            d_jacobi_angle0[0] = sin(joint_angle[1]) / (L_BE * sin(joint_angle[0] - joint_angle[1]) * sin(joint_angle[0] - joint_angle[1]));
            d_jacobi_angle0[1] = cos(joint_angle[1]) / (L_BE * sin(joint_angle[0] - joint_angle[1]) * sin(joint_angle[0] - joint_angle[1]));
            d_jacobi_angle0[2] = -sin(joint_angle[1]) * cos(joint_angle[0] - joint_angle[1]) / (L_AB * sin(joint_angle[0] - joint_angle[1]) * sin(joint_angle[0] - joint_angle[1]));
            d_jacobi_angle0[3] = -cos(joint_angle[1]) * cos(joint_angle[0] - joint_angle[1]) / (L_AB * sin(joint_angle[0] - joint_angle[1]) * sin(joint_angle[0] - joint_angle[1]));

            d_jacobi_angle1[0] = -sin(joint_angle[0]) * cos(joint_angle[0] - joint_angle[1]) / (L_BE * sin(joint_angle[0] - joint_angle[1]) * sin(joint_angle[0] - joint_angle[1]));
            d_jacobi_angle1[1] = -cos(joint_angle[0]) * cos(joint_angle[0] - joint_angle[1]) / (L_BE * sin(joint_angle[0] - joint_angle[1]) * sin(joint_angle[0] - joint_angle[1]));
            d_jacobi_angle1[2] = sin(joint_angle[0]) / (L_AB * sin(joint_angle[0] - joint_angle[1]) * sin(joint_angle[0] - joint_angle[1]));
            d_jacobi_angle1[3] = cos(joint_angle[0]) / (L_AB * sin(joint_angle[0] - joint_angle[1]) * sin(joint_angle[0] - joint_angle[1]));

            for(int i=0;i<4;i++)
            {
                d_jacobi_out_x[i] = d_jacobi_angle0[i] * jacobi[0] + d_jacobi_angle1[i] * jacobi[2];
                d_jacobi_out_y[i] = d_jacobi_angle0[i] * jacobi[1] + d_jacobi_angle1[i] * jacobi[3];
            }
        }
    }
}
