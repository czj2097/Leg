#include "TimeOptimalMotion.h"
#include "LegKinematics.h"
using namespace robot_app::kinematics;

void matrix_dot_matrix(double *matrix1, int matrix1_row, int matrix1_col, double *matrix2, int matrix2_col, double *matrix_out)
{
    for(int i=0;i<matrix1_row;i++)
    {
        for(int k=0;k<matrix2_col;k++)
        {
            *(matrix_out+matrix2_col*i+k)=0;
            for(int j=0;j<matrix1_col;j++)
            {
                *(matrix_out+matrix2_col*i+k) += *(matrix1+matrix1_col*i+j) * *(matrix2+matrix2_col*j+k);
            }
        }
    }
}

void dlmwrite(const char *filename, const double *mtx, const int m, const int n)
{
    std::ofstream file;

    file.open(filename);

    file << std::setprecision(15);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            file << mtx[n*i + j] << "   ";
        }
        file << std::endl;
    }
}

void TimeOptimalMotionSingleEffector::Initialize()
{
    //init & set tippos here
    std::fill_n(ds_backward,901,0);
    std::fill_n(ds_forward,901,0);
    switchCount=0;
    std::fill_n(switchPoint,901,-1);
    std::fill_n(switchScrewID,901,-1);
    std::fill_n(switchType,901,'0');
    std::fill_n(*isParamddsExact0,901*2,-1);

    for(int i = 0; i < 901; i++)
    {
        s[i] = PI * i/900;
        TipPos[i][0] = initTipPos[0] + stepD/2 * cos(PI/2 * (1 - cos(s[i]))) - (stepD/4 - stepD/2 * s[i]/PI);//D/4 --> -D/4
        TipPos[i][1] = initTipPos[1] + stepH * sin(PI/2 * (1 - cos(s[i])));
    }
}

void TimeOptimalMotionSingleEffector::GetParam()
{
    double d_TipPos_s[2];
    double dd_TipPos_s[2];
    double jacobi[4];
    double d_jacobi_x[4];
    double d_jacobi_y[4];
    double d_jacobi_s[4];

    double param_dsds1[2];
    double param_dsds2[2];

    for(int i = 0; i < 901; i++)
    {
        d_TipPos_s[0] = -stepD/2 * sin(PI/2 * (1 - cos(s[i]))) * sin(s[i]) + stepD/2 / PI;
        d_TipPos_s[1] = stepH * cos(PI/2 * (1 - cos(s[i]))) * sin(s[i]);

        dd_TipPos_s[0] = -stepD/2 * (cos(PI/2 * (1 - cos(s[i]))) * sin(s[i]) * sin(s[i]) + sin(PI/2 * (1 - cos(s[i]))) * cos(s[i]));
        dd_TipPos_s[1] = stepH * (-sin(PI/2 * (1 - cos(s[i]))) * sin(s[i]) *sin(s[i]) + cos(PI/2 * (1 - cos(s[i]))) * cos(s[i]));

        Leg::LegIJ(*TipPos+2*i,jacobi,1);
        Leg::LegIdJ(*TipPos+2*i,d_jacobi_x,d_jacobi_y,1);

        for(int j = 0; j < 4; j++)
        {
            d_jacobi_s[j] = d_jacobi_x[j] * d_TipPos_s[0] + d_jacobi_y[j] * d_TipPos_s[1];
        }

        std::fill_n(param_dsds1,2,0);
        matrix_dot_matrix(jacobi,2,2,dd_TipPos_s,1,param_dsds1);
        std::fill_n(param_dsds2,2,0);
        matrix_dot_matrix(d_jacobi_s,2,2,d_TipPos_s,1,param_dsds2);
        std::fill_n(*param_dds+2*i,2,0);
        matrix_dot_matrix(jacobi,2,2,d_TipPos_s,1,*param_dds+2*i);
        std::fill_n(*param_dsds+2*i,2,0);
        for (int j = 0; j < 2; j++)
        {
            param_dsds[i][j] = param_dsds1[j] + param_dsds2[j];
            abs_param_dds[i][j]=fabs(param_dds[i][j]);

            if(param_dds[i][j]>0)
            {
                param_a2[i][j]=-param_dsds[i][j]/param_dds[i][j];
                param_a0L[i][j]=-aLmt/param_dds[i][j];
                param_a0H[i][j]=aLmt/param_dds[i][j];
            }
            else if(param_dds[i][j]<0)
            {
                param_a2[i][j]=-param_dsds[i][j]/param_dds[i][j];
                param_a0L[i][j]=aLmt/param_dds[i][j];
                param_a0H[i][j]=-aLmt/param_dds[i][j];
            }
            else
            {
                isParamddsExact0[i][j]=1;
                printf("WARNING!!! param_dds equals zero!!!");
            }
        }
    }
}

double TimeOptimalMotionSingleEffector::GetMaxDec(int count, double ds)
{
    double dec[2] {0};
    std::fill_n(dec,2,-1e6);

    for (int j=0;j<2;j++)
    {
        if(isParamddsExact0[count][j]==-1)
        {
            dec[j]=param_a2[count][j]*ds*ds+param_a0L[count][j];
        }
    }

    return *std::max_element(dec,dec+2);
}

double TimeOptimalMotionSingleEffector::GetMinAcc(int count, double ds)
{
    double acc[2] {0};
    std::fill_n(acc,2,1e6);

    for (int j=0;j<2;j++)
    {
        if(isParamddsExact0[count][j]==-1)
        {
            acc[j]=param_a2[count][j]*ds*ds+param_a0H[count][j];
        }
    }

    return *std::min_element(acc,acc+2);
}

void TimeOptimalMotionSingleEffector::GetDsBound(int count)
{
    int k_st {0};
    bool dsBoundFlag_st {false};
    const int kstCount {15000};
    while (dsBoundFlag_st==false)
    {
        double ds=0.001*k_st;
        double max_dec=GetMaxDec(count,ds);
        double min_acc=GetMinAcc(count,ds);

        k_st++;
        if(k_st==kstCount)
        {
            dsBoundFlag_st=true;
            printf("WARNING!!! kstCount=%d is too small!!!\n",kstCount);
        }
        else
        {
            if(min_acc>max_dec)
            {
                ds_upBound_aLmt[count]=ds;
            }
            else
            {
                for(int k_st2=0;k_st2<1000;k_st2++)
                {
                    ds=ds_upBound_aLmt[count]+0.001*0.001*k_st2;
                    max_dec=GetMaxDec(count,ds);
                    min_acc=GetMinAcc(count,ds);
                    if(min_acc>=max_dec)
                    {
                        dds_lowBound[count]=max_dec;
                        dds_upBound[count]=min_acc;
                        ds_upBound_aLmt[count]=ds;
                    }
                    else
                    {
                        dsBoundFlag_st=true;
                        break;
                    }
                }
            }
        }
    }

    ds_upBound_vLmt[count]=vLmt/(*std::max_element(*abs_param_dds+2*count,*abs_param_dds+2*count+2));
    ds_upBound[count]=std::min(ds_upBound_aLmt[count],ds_upBound_vLmt[count]);
}

void TimeOptimalMotionSingleEffector::GetSwitchPoint()
{
    double slopedsBound[901] {0};
    double paramdds0Point[901] {0};
    int paramdds0Count {0};

    double tangentPoint[901] {0};
    int tangentCount {0};
    int switchScrewID_tmp[901] {0};

    //initialize
    for(int i=0;i<901;i++)
    {
        tangentPoint[i]=-1;
        paramdds0Point[i]=-1;
        switchScrewID_tmp[i]=-1;
    }
    tangentCount=0;
    paramdds0Count=0;

    //known switch point
    switchPoint[switchCount]=0.0;
    switchType[switchCount]='b';
    switchCount++;
    switchPoint[switchCount]=900.0;
    switchType[switchCount]='b';
    switchCount++;

    //calculate the slope of ds_upBound
    slopedsBound[0]=(ds_upBound[1]-ds_upBound[0])/(s[1]-s[0]);
    for(int i=1;i<900;i++)
    {
        slopedsBound[i]=( (ds_upBound[i+1]-ds_upBound[i])/(s[i+1]-s[i])
                              +(ds_upBound[i]-ds_upBound[i-1])/(s[i]-s[i-1]) )/2;
    }
    slopedsBound[900]=(ds_upBound[900]-ds_upBound[900-1])/(s[900]-s[900-1]);
    for(int i=0;i<900+1;i++)
    {
        slopeDelta[i]=slopedsBound[i]-(dds_upBound[i]+dds_lowBound[i])/2;
    }

    for(int i=0;i<900;i++)
    {
        bool isParamdds0 {false};
        if(ds_upBound_vLmt[i]>=ds_upBound_aLmt[i])
        {
            for(int j=0;j<2;j++)
            {
                    if(param_dds[i+1][j]*param_dds[i][j]<0 || param_dds[i][j]==0)
                    {
                        paramdds0Point[paramdds0Count]=i
                                +fabs(param_dds[i][j])/(fabs(param_dds[i][j])+fabs(param_dds[i+1][j]));
                        switchScrewID_tmp[paramdds0Count]=j;
                        paramdds0Count++;
                        isParamdds0=true;
                    }
            }
        }

        if(slopeDelta[i]<=0 && slopeDelta[i+1]>0 && isParamdds0==false)
        {
            tangentPoint[tangentCount]=i
                    +fabs(slopeDelta[i])/(fabs(slopeDelta[i])+fabs(slopeDelta[i+1]));
            tangentCount++;
        }
    }

    printf("StanceLeg Tangent Switch Point:");
    for(int i=0;i<tangentCount+1;i++)
    {
        printf("%.2f,",tangentPoint[i]);
    }
    printf("\n");
    printf("StanceLeg ZeroInertia Switch Point:");
    for(int i=0;i<paramdds0Count+1;i++)
    {
        printf("%.2f,",paramdds0Point[i]);
    }
    printf("\n");

    //merge tangentPoint & paramdds0Point into switchPoint
    for(int i=0;i<tangentCount;i++)
    {
        switchPoint[i+switchCount]=tangentPoint[i];
        switchType[i+switchCount]='t';
    }
    switchCount+=tangentCount;
    for(int i=0;i<paramdds0Count;i++)
    {
        switchPoint[i+switchCount]=paramdds0Point[i];
        switchType[i+switchCount]='z';
        switchScrewID[i+switchCount]=switchScrewID_tmp[i];
    }
    switchCount+=paramdds0Count;

    //filtering the same point & sorting by the value
    for(int i=0;i<switchCount;i++)
    {
        for(int j=i+1;j<switchCount;j++)
        {
            if(switchPoint[j]<switchPoint[i])
            {
                double tmp1=switchPoint[i];
                switchPoint[i]=switchPoint[j];
                switchPoint[j]=tmp1;

                int tmp2=switchScrewID[i];
                switchScrewID[i]=switchScrewID[j];
                switchScrewID[j]=tmp2;

                char tmp3=switchType[i];
                switchType[i]=switchType[j];
                switchType[j]=tmp3;
            }
        }
    }

    printf("Switch Point:");
    for(int i=0;i<switchCount+1;i++)
    {
        printf("%.4f,",s[(int)switchPoint[i]]+(switchPoint[i]-(int)switchPoint[i])*(s[(int)switchPoint[i]+1]-s[(int)switchPoint[i]]));
    }
    printf("\n");
    printf("Switch Type:");
    for(int i=0;i<switchCount+1;i++)
    {
        printf("%c,",switchType[i]);
    }
    printf("\n");
}

void TimeOptimalMotionSingleEffector::GetTwoPointAtSwitch(double *lowPoint, double *upPoint)
{
    lowPoint[0]=-1;
    lowPoint[switchCount-1]=ds_upBound[900];
    upPoint[0]=ds_upBound[0];
    upPoint[switchCount-1]=-1;

    for(int i=1;i<switchCount-1;i++)
    {
        if(switchType[i]=='t')
        {
            int num=(int)switchPoint[i];
            upPoint[i]=ds_upBound[num+1];
            lowPoint[i]=ds_upBound[num];
        }
        else if(switchType[i]=='d')
        {
            int num=(int)switchPoint[i];
            upPoint[i]=lowPoint[i]=std::min(ds_upBound[num],ds_upBound[num-1]);
        }
        else if(switchType[i]=='z')
        {
            int num=round(switchPoint[i]);
            double tmp=std::min(ds_upBound[num-1],ds_upBound[num]);
            upPoint[i]=lowPoint[i]=std::min(tmp,ds_upBound[num+1]);
        }
    }
}

void TimeOptimalMotionSingleEffector::GetOptimalDsBySwitchPoint()
{
    bool stopFlag {false};
    int forwardEnd_s {-1};
    double forwardEnd_ds {0};
    double *lowPoint=new double [switchCount];
    double *upPoint=new double [switchCount];
    GetTwoPointAtSwitch(lowPoint,upPoint);

    printf("lowPoint:");
    for(int i=0;i<switchCount;i++)
    {
        printf("%.2f,",lowPoint[i]);
    }
    printf("\n");
    printf("upPoint:");
    for(int i=0;i<switchCount;i++)
    {
        printf("%.2f,",upPoint[i]);
    }
    printf("\n");

    for(int i=0;i<900+1;i++)
    {
        real_ds[i]=ds_upBound[i];
        real_dds[i]=dds_upBound[i];
    }
    for(int m=0;m<switchCount;m++)
    {
        //start of backward
        bool ignoreBackward {false};
        int k_st=(int)switchPoint[m];
        int k_st_start=k_st;
        if(switchType[m]=='z')
        {
            k_st_start=k_st=round(switchPoint[m])-1;
        }
        else if(switchType[m]=='d')
        {
            if(ds_upBound[k_st_start]>ds_upBound[k_st_start-1])
            {
                k_st_start=k_st=(int)switchPoint[m]-1;
            }
        }

        if(k_st_start==0)
        {
            ignoreBackward=true;
        }
        else if(k_st_start<forwardEnd_s)
        {
            ignoreBackward=true;
            printf("backward start at a passed point, quit switchPoint %.1f\n",switchPoint[m]);
        }
        else if(k_st_start==forwardEnd_s && lowPoint[m]>forwardEnd_ds)
        {
            ignoreBackward=true;
            if(switchType[m]=='z')
            {
                real_dds[k_st_start+1]=0;
                real_ds[k_st_start+1]=lowPoint[m];
            }
        }

        if(ignoreBackward==false)
        {
            if(switchType[m]=='z')
            {
                real_dds[k_st_start+1]=0;
                real_ds[k_st_start+1]=lowPoint[m];
            }
            stopFlag=false;
            ds_backward[k_st]=lowPoint[m];
            while(stopFlag==false)
            {
                dds_backward[k_st]=GetMaxDec(k_st,ds_backward[k_st]);
                ds_backward[k_st-1]=sqrt(ds_backward[k_st]*ds_backward[k_st]-2*dds_backward[k_st]*(s[k_st]-s[k_st-1]));

                if(ds_backward[k_st-1]>ds_upBound[k_st-1])
                {
                    real_dds[k_st-1]=(real_ds[k_st-1]-ds_backward[k_st])*(real_ds[k_st-1]+ds_backward[k_st])/2/(s[k_st]-s[k_st-1]);
                    for(int i=k_st;i<k_st_start+1;i++)
                    {
                        real_ds[i]=ds_backward[i];
                        real_dds[i]=dds_backward[i];
                    }
                    stopFlag=true;
                    printf("backward touching upBound at %d, from switchPoint %.1f\n",k_st-1,switchPoint[m]);
                }
                else if(k_st==1)
                {
                    dds_backward[k_st-1]=GetMaxDec(k_st-1,ds_backward[k_st-1]);
                    for(int i=k_st-1;i<k_st_start+1;i++)
                    {
                        real_ds[i]=ds_backward[i];
                        real_dds[i]=dds_backward[i];
                    }
                    stopFlag=true;
                    printf("StanceLeg backward touching 0, from switchPoint %.1f\n",switchPoint[m]);
                }
                else if(ds_backward[k_st-1]>=real_ds[k_st-1])
                {
                    real_dds[k_st-1]=(real_ds[k_st-1]-ds_backward[k_st])*(real_ds[k_st-1]+ds_backward[k_st])/2/(s[k_st]-s[k_st-1]);
                    for(int i=k_st;i<k_st_start+1;i++)
                    {
                        real_ds[i]=ds_backward[i];
                        real_dds[i]=dds_backward[i];
                    }
                    stopFlag=true;
                    printf("backward touching last curve at %d, from switchPoint %.1f\n",k_st-1,switchPoint[m]);
                }
                else
                {
                    k_st--;
                }
            }
        }

        //start of forward
        bool ignoreForward {false};
        k_st_start=k_st=(int)switchPoint[m];
        if(switchType[m]=='t')
        {
            k_st_start=k_st=(int)switchPoint[m]+1;
        }
        else if(switchType[m]=='z')
        {
            k_st_start=k_st=round(switchPoint[m])+1;
        }
        else if(switchType[m]=='d')
        {
            if(ds_upBound[k_st_start]>ds_upBound[k_st_start-1])
            {
                k_st_start=k_st=(int)switchPoint[m]-1;
            }
        }

        if(k_st_start==900+1 || k_st_start==900)
        {
            ignoreForward=true;
        }
        else if(k_st_start<forwardEnd_s)
        {
            ignoreForward=true;
            printf("forward start at a passed point, quit switchPoint %.1f\n",switchPoint[m]);
        }
        else if(k_st_start==forwardEnd_s)
        {
            if(upPoint[m]>forwardEnd_ds)
            {
                printf("How possible! forward curve should not stop here!\n");
            }
        }

        if(ignoreForward==false)
        {
            stopFlag=false;
            ds_forward[k_st]=upPoint[m];
            while(stopFlag==false)
            {
                dds_forward[k_st]=GetMinAcc(k_st,ds_forward[k_st]);
                ds_forward[k_st+1]=sqrt(ds_forward[k_st]*ds_forward[k_st]+2*dds_forward[k_st]*(s[k_st+1]-s[k_st]));

                if(ds_forward[k_st+1]>ds_upBound[k_st+1] || k_st==900-1)
                {
                    forwardEnd_s=k_st;
                    forwardEnd_ds=ds_forward[k_st];
                    for(int i=k_st_start;i<k_st+1;i++)
                    {
                        real_ds[i]=ds_forward[i];
                        real_dds[i]=dds_forward[i];
                    }
                    if(k_st==900-1)
                    {
                        if(ds_forward[k_st+1]>ds_upBound[k_st+1])
                        {
                            real_ds[900]=ds_upBound[900];
                            real_dds[900]=dds_upBound[900];
                        }
                        else
                        {
                            real_ds[900]=ds_forward[900];
                            real_dds[900]=GetMinAcc(900,ds_forward[900]);
                            forwardEnd_s=900;
                            forwardEnd_ds=ds_forward[900];
                        }
                    }
                    stopFlag=true;
                    printf("forward touching upBound at %d, from switchPoint %.4f\n",k_st,switchPoint[m]);
                }
                else
                {
                    k_st++;
                }
            }
        }
    }

    for(int i=0;i<900+1;i++)
    {
        real_ddsMax[i]=GetMinAcc(i,real_ds[i]);
        real_ddsMin[i]=GetMaxDec(i,real_ds[i]);
    }

    delete [] lowPoint;
    delete [] upPoint;
}

void TimeOptimalMotionSingleEffector::GetOptimalGait2t()
{
    //fot t
    double totalTime {0};
    int totalCount {0};
    double v0 {0};
    double vm {0};
    double vt {stepD/2/PI*real_ds[0]};
    double stance_begin_s;

    for (int i=0;i<901;i++)
    {
        totalTime+=2*(s[i]-s[i-1])/(real_ds[i-1]+real_ds[i]);
    }
    totalCount=(int)(totalTime*1000)+1;
    printf("totalTime is %.4f, totalCount is %d\n",totalTime,totalCount);

    double *real_s = new double [totalCount];
    double *real_Pee = new double [2*totalCount];
    double *real_Pin = new double [2*totalCount];
    real_s[0]=0;
    for (int i=1;i<totalCount;i++)
    {
        double ds=0.5*(real_ds[(int)(real_s[i-1]/PI*900)]+real_ds[(int)(real_s[i-1]/PI*900)+1]);
        real_s[i]=real_s[i-1]+ds*0.001;
        if (i==totalCount-1)
        {
            double dds=0.5*(real_dds[(int)(real_s[i-1]/PI*900)]+real_dds[(int)(real_s[i-1]/PI*900)+1]);
            v0=stepD/2/PI*(ds+dds*0.001);
            stance_begin_s=real_s[totalCount-1]+ds*0.001;
        }
    }
    vm=stepD/(0.001*totalCount)-(v0+vt)/2;

    for (int i=0;i<2*totalCount;i++)
    {
        //swing phase
        if(i<totalCount)
        {
            *(real_Pee+2*i) = initTipPos[0] + stepD/2 * cos(PI/2 * (1 - cos(real_s[i]))) - (stepD/4 - stepD/2 * real_s[i]/PI);//D/4 --> -D/4
            *(real_Pee+2*i+1) = initTipPos[1] + stepH * sin(PI/2 * (1 - cos(real_s[i])));
        }
        //stance phase
        else
        {
            *(real_Pee+2*i+1) = initTipPos[1];
            if((i-totalCount)<(double)totalCount/2)
            {
                *(real_Pee+2*i) = initTipPos[0] + stepD/2 * cos(PI/2 * (1 - cos(stance_begin_s))) - (stepD/4-stepD/2*(stance_begin_s/PI))
                        + v0*(0.001*(i-totalCount)) + 0.5*(vm-v0)/(0.001*totalCount/2) * 0.001*(i-totalCount) * 0.001*(i-totalCount);
            }
            else
            {
                *(real_Pee+2*i) = initTipPos[0] + stepD/2 * cos(PI/2 * (1-cos(stance_begin_s))) - (stepD/4 - stepD/2*(stance_begin_s/PI))
                        + v0*(0.001*totalCount/2) + 0.5*(vm-v0)/(0.001*totalCount/2) * 0.001*totalCount/2 * 0.001*totalCount/2
                        + vm*(0.001*(i-1.5*totalCount)) + 0.5*(vt-vm)/(0.001*totalCount/2) * 0.001*(i-1.5*totalCount) * 0.001*(i-1.5*totalCount);
            }
        }
    }

    dlmwrite("./real_Pee.txt",real_Pee,totalCount,2);
    dlmwrite("./real_Pin.txt",real_Pin,totalCount,2);

    delete [] real_s;
    delete [] real_Pee;
    delete [] real_Pin;

}

void TimeOptimalMotionSingleEffector::outputData()
{
    printf("Start output data...\n");
    dlmwrite("./ds_upBound_aLmt.txt",ds_upBound_aLmt,901,1);
    dlmwrite("./ds_upBound_vLmt.txt",ds_upBound_vLmt,901,1);
    dlmwrite("./dds_upBound.txt",dds_upBound,901,1);
    dlmwrite("./dds_lowBound.txt",dds_lowBound,901,1);
    dlmwrite("./ds_forward.txt",ds_forward,901,1);
    dlmwrite("./ds_backward.txt",ds_backward,901,1);
    dlmwrite("./dds_forward.txt",dds_forward,901,1);
    dlmwrite("./dds_backward.txt",dds_backward,901,1);
    printf("Finish output data.\n");
}
