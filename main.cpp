#include "LegKinematics.h"
#include <iostream>

int main()
{
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