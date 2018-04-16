#ifndef TIMEOPTIMALMOTION_H
#define TIMEOPTIMALMOTION_H
#include <stdio.h>
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <algorithm>

namespace time_optimal
{
    void matrix_dot_matrix(double *matrix1, int matrix1_row, int matrix1_col, double *matrix2, int matrix2_col, double *matrix_out);
    void dlmwrite(const char *filename, const double *mtx, const int m, const int n);

    class TimeOptimalMotionSingleEffector
    {
    public:
        void GetTimeOptimalGait(double step_length, double step_height, double acc_limit, double vel_limit, double y_of_tip, double *out_tippos, double &out_period);

    private:
        void Initialize(double step_length, double step_height, double acc_limit, double vel_limit, double y_of_tip);
        void GetParam();
        void GetDsBound(int count);
        void GetSwitchPoint();
        void GetOptimalDsBySwitchPoint();
        void GetConstVelocityGait();
        void GetOptimalGait2t(double *out_tippos, double &out_period);
        void outputData();

        double GetMaxDec(int count, double ds);
        double GetMinAcc(int count, double ds);
        void GetTwoPointAtSwitch(double *lowPoint, double *upPoint);
        void ApplyExtraItegrationToBoundaryPoint(int s_count, double ds);

        double s[901];
        double stepD;
        double stepH;
        double aLmt;
        double vLmt;
        double initTipPos[2];
        double TipPos[901][2];

        double param_dds[901][2];
        double abs_param_dds[901][2];
        double param_dsds[901][2];
        double param_a2[901][2];
        double param_a0L[901][2];
        double param_a0H[901][2];
        int isParamddsExact0[901][2];

        double ds_upBound[901];
        double ds_upBound_aLmt[901];
        double ds_upBound_vLmt[901];
        double dds_upBound[901];
        double dds_lowBound[901];

        int dis_count1;
        int dis_count2;
        int switchCount;
        double switchPoint[901];
        char switchType[901];
        double slopeDelta[901];
        int switchScrewID[901];

        double real_ds[901];
        double real_dds[901];
        double real_ddsMax[901];
        double real_ddsMin[901];
        double ds_forward[901];
        double ds_backward[901];
        double dds_forward[901];
        double dds_backward[901];

        int totalCount;
        double v0;
        double vt;
    };

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
            static void LegIJ(double *tip_pos_in, double *jacobi_out, double leg_orient = 1);
            static void LegIdJ(double *tip_pos_in, double *d_jacobi_out_x, double *d_jacobi_out_y, double leg_orient = 1);
        };
    }
}


#endif
