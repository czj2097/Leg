#include "TimeOptimalMotion.h"

int main()
{
    double stepD {0.6};
    double stepH {0.05};
    double aLmt {1e6/4096*2*time_optimal::kinematics::PI / 81 * 1.0};
    double vLmt {100*time_optimal::kinematics::PI / 81 * 1.0};
    double initTipY {-0.45};
    double out_TipPos[2000][2] {0};
    double out_period {0};
    int period_count {0};

    time_optimal::TimeOptimalMotionSingleEffector planner;
    planner.GetTimeOptimalGait(stepD,stepH,aLmt,vLmt,initTipY,*out_TipPos,out_period);

    printf("out_period=%.3f\n",out_period);
    period_count=(int)(out_period*1000);
    time_optimal::dlmwrite("./log/out_TipPos.txt",*out_TipPos,period_count,2);

    return 0;
}
