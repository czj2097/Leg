#include "LegKinematics.h"
#include "TimeOptimalMotion.h"

void testMatrixDot()
{
    double matrix1[3][2] = {2,3,
                            1,5,
                            4,7,};
    double matrix2[2][4] = {1,2,3,4,
                            9,8,7,6};
    double matrix[3][4] = {0};
    matrix_dot_matrix(*matrix1,3,2,*matrix2,4,*matrix);
//    printf("matrix = \n{%d,%d,%d,%d,\n%d,%d,%d,%d,\n%d,%d,%d,%d}",matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3]
//            ,matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3]
//            ,matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3]);
    std::cout<<matrix[0][0]<<",  "<<matrix[0][1]<<",  "<<matrix[0][2]<<",  "<<matrix[0][3]<<std::endl
                <<matrix[1][0]<<",  "<<matrix[1][1]<<",  "<<matrix[1][2]<<",  "<<matrix[1][3]<<std::endl
                <<matrix[2][0]<<",  "<<matrix[2][1]<<",  "<<matrix[2][2]<<",  "<<matrix[2][3]<<std::endl;
}

void testTimeOptimal()
{
    TimeOptimalMotionSingleEffector planner;
    planner.Initialize();
    planner.GetParam();
    for(int i=0;i<901;i++)
    {
        planner.GetDsBound(i);
    }
    planner.GetSwitchPoint();
    planner.GetOptimalDsBySwitchPoint();
    planner.ApplyExtraItegration();
    planner.outputData();
    planner.GetOptimalGait2t();
    planner.GetNormalGait();
    planner.GetEntireGait();
}

int main()
{
    testTimeOptimal();

    using namespace std;
    robot_app::kinematics::Leg leg;
    
    double tip_pos_in[2] = {-0.45, 0.03};
    double joint_pos_out[2];

    double joint_pos_in[2];
    double tip_pos_out[2];

    // inverse 
    leg.LegIK(tip_pos_in, joint_pos_out);
    cout << "Joint position: (" << joint_pos_out[0] << ", " << joint_pos_out[1] << ")" << endl;
    
    // forward
    for(int i = 0; i < 2; i++)
    {
        joint_pos_in[i] = joint_pos_out[i];
    }
    leg.LegFK(joint_pos_in, tip_pos_out);
	cout << "Tip position: (" << tip_pos_out[0] << ", " << tip_pos_out[1] << ")" << endl;
	cin.get();
    return 0;
}
