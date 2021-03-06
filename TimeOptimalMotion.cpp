#include "TimeOptimalMotion.h"

namespace time_optimal
{
    using namespace kinematics;

    void TimeOptimalMotionSingleEffector::GetTimeOptimalGait(double duty_cycle, double step_length, double step_height, double acc_limit, double vel_limit, double y_of_tip, double *out_tippos, double &out_period)
    {
        Initialize(duty_cycle,step_length, step_height, acc_limit, vel_limit, y_of_tip);
        GetParam();
        for(int i=0;i<901;i++)
        {
            GetDsBound(i);
        }
        GetSwitchPoint();
        GetOptimalDsBySwitchPoint();
        ApplyExtraItegration();
        GetConstVelocityGait();
        GetOptimalGait2t(out_tippos,out_period);

        outputData();
        GetNormalGait();
        GetEntireGait();
    }

    void TimeOptimalMotionSingleEffector::Initialize(double duty_cycle, double step_length, double step_height, double acc_limit, double vel_limit, double y_of_tip)
    {
        //init & set tippos here
        dutyCycle=duty_cycle;
        stepD = step_length;
        stepH = step_height;
        aLmt = acc_limit;
        vLmt = vel_limit;
        initTipPos[0] = 0;
        initTipPos[1] = y_of_tip;

        std::fill_n(ds_backward,901,0);
        std::fill_n(ds_forward,901,0);
        dis_count1=100;
        dis_count2=800;
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
            d_TipPos_s[0] = -stepD/2 * sin(PI/2 * (1 - cos(s[i]))) * PI/2*sin(s[i]) + stepD/2 / PI;
            d_TipPos_s[1] = stepH * cos(PI/2 * (1 - cos(s[i]))) * PI/2*sin(s[i]);

            dd_TipPos_s[0] = -stepD/2 * (cos(PI/2 * (1 - cos(s[i]))) * PI/2*sin(s[i]) * PI/2*sin(s[i]) + sin(PI/2 * (1 - cos(s[i]))) * PI/2*cos(s[i]));
            dd_TipPos_s[1] = stepH * (-sin(PI/2 * (1 - cos(s[i]))) * PI/2*sin(s[i]) *PI/2*sin(s[i]) + cos(PI/2 * (1 - cos(s[i]))) * PI/2*cos(s[i]));

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
                double aLmt_tmp;
                if(i<dis_count1 || i>=dis_count2)
                {
                    aLmt_tmp=aLmt/2;
                }
                else
                {
                    aLmt_tmp=aLmt;
                }

                if(param_dds[i][j]>0)
                {
                    param_a2[i][j]=-param_dsds[i][j]/param_dds[i][j];
                    param_a0L[i][j]=-aLmt_tmp/param_dds[i][j];
                    param_a0H[i][j]=aLmt_tmp/param_dds[i][j];
                }
                else if(param_dds[i][j]<0)
                {
                    param_a2[i][j]=-param_dsds[i][j]/param_dds[i][j];
                    param_a0L[i][j]=aLmt_tmp/param_dds[i][j];
                    param_a0H[i][j]=-aLmt_tmp/param_dds[i][j];
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
        const int kstCount {50000};
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

        double vLmt_tmp;
        if(count<dis_count1 || count>=dis_count2)
        {
            vLmt_tmp=vLmt/2;
        }
        else
        {
            vLmt_tmp=vLmt;
        }
        ds_upBound_vLmt[count]=vLmt_tmp/(*std::max_element(*abs_param_dds+2*count,*abs_param_dds+2*count+2));
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
        switchPoint[switchCount] = ds_upBound[dis_count1] > ds_upBound[dis_count1-1] ? (dis_count1-1) : dis_count1;
        switchType[switchCount]='d';
        switchCount++;
        switchPoint[switchCount] = ds_upBound[dis_count2] > ds_upBound[dis_count2-1] ? (dis_count2-1) : dis_count2;
        switchType[switchCount]='d';
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

        printf("Tangent Switch Point:");
        for(int i=0;i<tangentCount+1;i++)
        {
            printf("%.2f,",tangentPoint[i]);
        }
        printf("\n");
        printf("ZeroInertia Switch Point:");
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

    void TimeOptimalMotionSingleEffector::ApplyExtraItegration()
    {
        bool stopFlag {false};
        int k {0};
        int k_start {0};
        double real_ds_tmp {0};
        if(real_ds[0]>real_ds[900])
        {
            //forward
            real_ds[0]=real_ds[900];
            k_start=k=0;
            stopFlag=false;
            while(stopFlag==false)
            {
                real_dds[k]=GetMinAcc(k,real_ds[k]);
                real_ds_tmp=sqrt(real_ds[k]*real_ds[k]+2*real_dds[k]*(s[k+1]-s[k]));

                if(real_ds_tmp>real_ds[k+1])
                {
                    stopFlag=true;
                }
                else
                {
                    real_ds[k+1]=real_ds_tmp;
                    k++;
                }
            }
        }
        else if(real_ds[0]<real_ds[900])
        {
            //backward
            real_ds[0]=real_ds[900];
            k_start=k=900;
            stopFlag=false;
            while(stopFlag==false)
            {
                real_dds[k]=GetMaxDec(k,real_ds[k]);
                real_ds_tmp=sqrt(real_ds[k]*real_ds[k]-2*real_dds[k]*(s[k]-s[k-1]));

                if(real_ds_tmp>real_ds[k-1])
                {
                    stopFlag=true;
                }
                else
                {
                    real_ds[k-1]=real_ds_tmp;
                    k--;
                }
            }
        }
        else
        {
            printf("Amazing!!!Ds[start] equal Ds[end].No need to apply extra itegration.\n");
        }

    }

    void TimeOptimalMotionSingleEffector::ApplyExtraItegrationToBoundaryPoint(int s_count, double ds)
    {
        bool stopFlag {false};
        int k {0};
        int k_start {0};
        double real_ds_tmp {0};
        if(s_count==0)
        {
            //forward
            real_ds[0]=ds;
            k_start=k=0;
            stopFlag=false;
            while(stopFlag==false)
            {
                real_dds[k]=GetMinAcc(k,real_ds[k]);
                real_ds_tmp=sqrt(real_ds[k]*real_ds[k]+2*real_dds[k]*(s[k+1]-s[k]));

                if(real_ds_tmp>real_ds[k+1])
                {
                    stopFlag=true;
                }
                else
                {
                    real_ds[k+1]=real_ds_tmp;
                    k++;
                }
            }
        }
        else if(s_count==900)
        {
            //backward
            real_ds[900]=ds;
            k_start=k=900;
            stopFlag=false;
            while(stopFlag==false)
            {
                real_dds[k]=GetMaxDec(k,real_ds[k]);
                real_ds_tmp=sqrt(real_ds[k]*real_ds[k]-2*real_dds[k]*(s[k]-s[k-1]));

                if(real_ds_tmp>real_ds[k-1])
                {
                    stopFlag=true;
                }
                else
                {
                    real_ds[k-1]=real_ds_tmp;
                    k--;
                }
            }
        }
        else
        {
            printf("unknown count for boundary point, please check!");
        }
    }

    void TimeOptimalMotionSingleEffector::GetConstVelocityGait()
    {
        double totalTime=0;
        double avgTime = 1;

        int k=0;

        while(fabs(totalTime - avgTime) > 0.0001 && k<20)
        {
            printf("real_ds=%.4f\n",real_ds[0]);

            totalTime=0;
            for (int i=1;i<901;i++)
            {
                totalTime += 2*(s[i]-s[i-1])/(real_ds[i-1]+real_ds[i]);
            }

            double min_ds = real_ds[0] < real_ds[900] ? real_ds[0] : real_ds[900];
            double avgVel = (-stepD/2 * sin(PI/2 * (1 - cos(s[0]))) * sin(s[0]) + stepD/2 / PI) * min_ds;
            avgTime = stepD/2 / avgVel /dutyCycle*(1-dutyCycle);

            if(totalTime<=avgTime)
            {
                printf("totalTime=%.4f is smaller than avgTime=%.4f. apply time-scaling to real_ds.\n",totalTime,avgTime);
                for(int i=0;i<901;i++)
                {
                    real_ds[i]=real_ds[i] * totalTime / avgTime;
                }
                break;
            }
            else
            {
                printf("totalTime=%.4f is larger than avgTime=%.4f. apply extra integration to boundary points.\n",totalTime,avgTime);
                double new_ds = min_ds * avgTime / totalTime;
                ApplyExtraItegrationToBoundaryPoint(0,new_ds);
                ApplyExtraItegrationToBoundaryPoint(900,new_ds);
                k++;
                printf("min_ds=%.4f, new_ds=%.4f, time error=%.5f\n",min_ds,new_ds,fabs(totalTime - avgTime));
            }
        }

        printf("Finish GetConstVelocityGait, iteration count is %d\n",k);

    }

    void TimeOptimalMotionSingleEffector::GetOptimalGait2t(double *out_tippos, double &out_period)
    {
        //fot t
        double timeArray[901] {0};
        double totalTime {0};

        for (int i=1;i<901;i++)
        {
            timeArray[i] = timeArray[i-1] + 2*(s[i]-s[i-1])/(real_ds[i-1]+real_ds[i]);
        }
        totalCount = (int)(timeArray[900]*1000)+1;
        totalTime = totalCount*0.001;
        out_period = totalTime;
        printf("totalTime is %.4f, totalCount is %d\n",timeArray[900],totalCount);

        double timeArray_scale[901] {0};
        double real_ds_scale[901] {0};
        double real_dds_scale[901] {0};
        for(int i=0;i<901;i++)
        {
            timeArray_scale[i] = timeArray[i] / timeArray[900] * totalTime;
            real_ds_scale[i] = real_ds[i] / totalTime * timeArray[900];
            real_dds_scale[i] = real_dds[i] / totalTime * timeArray[900];
        }

        double *s_t = new double [totalCount+1];
        double *ds_t = new double [totalCount+1];
        double *dds_t = new double [totalCount+1];
        double *real_Pee = new double [4*totalCount];
        double *real_Pin = new double [4*totalCount];
        double *real_Vee = new double [2*totalCount];
        double *real_Vin = new double [2*totalCount];
    //    double *real_Aee = new double [2*totalCount];
    //    double *real_Ain = new double [2*totalCount];

        int k_start {0};
        for (int i=0;i<totalCount;i++)
        {
            for(int k=k_start; k<900; k++)
            {
                if(0.001*i>=timeArray_scale[k] && 0.001*i<timeArray_scale[k+1])
                {
                    k_start=k;
                    s_t[i] = s[k] + (s[k+1]-s[k]) * (0.001*i-timeArray_scale[k]) / (timeArray_scale[k+1]-timeArray_scale[k]);
                    ds_t[i] = real_ds_scale[k] + (real_ds_scale[k+1]-real_ds_scale[k]) * (0.001*i-timeArray_scale[k]) / (timeArray_scale[k+1]-timeArray_scale[k]);
                    dds_t[i] = real_dds_scale[k] + (real_dds_scale[k+1]-real_dds_scale[k]) * (0.001*i-timeArray_scale[k]) / (timeArray_scale[k+1]-timeArray_scale[k]);
                    break;
                }
            }
        }
        s_t[totalCount]=s[900];
        ds_t[totalCount]=real_ds_scale[900];
        dds_t[totalCount]=real_dds_scale[900];

        v0 = (-stepD/2 * sin(PI/2 * (1 - cos(s[0]))) * sin(s[0]) + stepD/2 / PI) * real_ds_scale[0];
        vt = (-stepD/2 * sin(PI/2 * (1 - cos(s[900]))) * sin(s[900]) + stepD/2 / PI) * real_ds_scale[900];
        double vm = stepD/totalTime - (v0+vt)/2;
        double stance_begin_s {PI};

        //swing phase
        for (int i=0;i<totalCount;i++)
        {
            *(real_Pee+2*i) = initTipPos[0] + stepD/2 * cos(PI/2 * (1 - cos(s_t[i]))) - (stepD/4 - stepD/2 * s_t[i]/PI);//D/4 --> -D/4
            *(real_Pee+2*i+1) = initTipPos[1] + stepH * sin(PI/2 * (1 - cos(s_t[i])));

            *(real_Vee+2*i) = (-stepD/2 * sin(PI/2 * (1 - cos(s_t[i]))) * PI/2*sin(s_t[i]) + stepD/2 / PI) * ds_t[i];
            *(real_Vee+2*i+1) = (stepH * cos(PI/2 * (1 - cos(s_t[i]))) * PI/2*sin(s_t[i])) * ds_t[i];
        }
        memcpy(out_tippos,real_Pee,2*totalCount*sizeof(double));

        //stance phase
        for (int i=totalCount; i<2*totalCount; i++)
        {
            *(real_Pee+2*i+1) = initTipPos[1];
            if((i-totalCount)<(double)totalCount/2)
            {
                *(real_Pee+2*i) = initTipPos[0] + stepD/2 * cos(PI/2 * (1 - cos(stance_begin_s))) - (stepD/4-stepD/2*(stance_begin_s/PI))
                        + v0 * 0.001*(i-totalCount) + 0.5 * (vm-v0)/(totalTime/2) * 0.001*(i-totalCount) * 0.001*(i-totalCount);
            }
            else
            {
                *(real_Pee+2*i) = initTipPos[0] + stepD/2 * cos(PI/2 * (1-cos(stance_begin_s))) - (stepD/4 - stepD/2*(stance_begin_s/PI))
                        + v0 * 0.001*totalCount/2 + 0.5 * (vm-v0)/(totalTime/2) * totalTime/2 * totalTime/2
                        + vm * (0.001*(i-1.5*totalCount)) + 0.5 * (vt-vm)/(totalTime/2) * 0.001*(i-1.5*totalCount) * 0.001*(i-1.5*totalCount);
            }
        }

        for(int i=0;i<2*totalCount;i++)
        {
            Leg::LegIK(real_Pee+2*i,real_Pin+2*i,1);
        }
        for(int i=0;i<totalCount;i++)
        {
            double jacobi[4];
            Leg::LegIJ(real_Pee+2*i,jacobi,1);
            matrix_dot_matrix(jacobi,2,2,real_Vee+2*i,1,real_Vin+2*i);
        }

        dlmwrite("./log/timeArray.txt",timeArray_scale,totalCount,1);
        dlmwrite("./log/s_t.txt",s_t,totalCount+1,1);
        dlmwrite("./log/real_Pee.txt",real_Pee,2*totalCount,2);
        dlmwrite("./log/real_Pin.txt",real_Pin,2*totalCount,2);
        dlmwrite("./log/real_ds.txt",real_ds_scale,901,1);
        dlmwrite("./log/real_Vee.txt",real_Vee,totalCount,2);
        dlmwrite("./log/real_Vin.txt",real_Vin,totalCount,2);

        delete [] s_t;
        delete [] ds_t;
        delete [] dds_t;
        delete [] real_Pee;
        delete [] real_Pin;
        delete [] real_Vee;
        delete [] real_Vin;
    }

    void TimeOptimalMotionSingleEffector::outputData()
    {
        printf("Start output data...\n");
        dlmwrite("./log/ds_upBound_aLmt.txt",ds_upBound_aLmt,901,1);
        dlmwrite("./log/ds_upBound_vLmt.txt",ds_upBound_vLmt,901,1);
        dlmwrite("./log/dds_upBound.txt",dds_upBound,901,1);
        dlmwrite("./log/dds_lowBound.txt",dds_lowBound,901,1);
        dlmwrite("./log/ds_forward.txt",ds_forward,901,1);
        dlmwrite("./log/ds_backward.txt",ds_backward,901,1);
        dlmwrite("./log/dds_forward.txt",dds_forward,901,1);
        dlmwrite("./log/dds_backward.txt",dds_backward,901,1);
        printf("Finish output data.\n");
    }

    void TimeOptimalMotionSingleEffector::GetNormalGait()
    {
        printf("Start GetNormalGait\n");
        int normalTotalCount=totalCount;
        double jacobi[4] {0};
        double d_jacobi[4] {0};
        double d_jacobi_x[4] {0};
        double d_jacobi_y[4] {0};
        double d_TipPos_s[2] {0};
        double dd_TipPos_s[2] {0};
        double maxAin {0};
        double maxVin {0};
        int k {0};
        bool stopFlag {false};

        while(stopFlag==false)
        {
            double * normalPee = new double [normalTotalCount*2];
            double * normalVee = new double [normalTotalCount*2];
            double * normalAee = new double [normalTotalCount*2];
            double * normalPin = new double [normalTotalCount*2];
            double * normalVin = new double [normalTotalCount*2];
            double * normalAin = new double [normalTotalCount*2];
            double * abs_Ain = new double [normalTotalCount*2];
            double * abs_Vin = new double [normalTotalCount*2];

            //swing phase
            for (int i=0;i<normalTotalCount;i++)
            {
                double s_n = PI*i/normalTotalCount;
                double ds_t = PI*1000/normalTotalCount;
                *(normalPee+2*i) = initTipPos[0] + stepD/2 * cos(PI/2 * (1 - cos(s_n))) - (stepD/4 - stepD/2 * s_n/PI);//D/4 --> -D/4
                *(normalPee+2*i+1) = initTipPos[1] + stepH * sin(PI/2 * (1 - cos(s_n)));

                d_TipPos_s[0] = -stepD/2 * sin(PI/2 * (1 - cos(s_n))) * PI/2*sin(s_n) + stepD/2 / PI;
                d_TipPos_s[1] = stepH * cos(PI/2 * (1 - cos(s_n))) * PI/2*sin(s_n);

                dd_TipPos_s[0] = -stepD/2 * (cos(PI/2 * (1 - cos(s_n))) * PI/2*sin(s_n) * PI/2*sin(s_n) + sin(PI/2 * (1 - cos(s_n))) * PI/2*cos(s_n));
                dd_TipPos_s[1] = stepH * (-sin(PI/2 * (1 - cos(s_n))) * PI/2*sin(s_n) *PI/2*sin(s_n) + cos(PI/2 * (1 - cos(s_n))) * PI/2*cos(s_n));

                *(normalVee+2*i) = d_TipPos_s[0] * ds_t;
                *(normalVee+2*i+1) = d_TipPos_s[1] * ds_t;

                *(normalAee+2*i) = dd_TipPos_s[0] * ds_t * ds_t;
                *(normalAee+2*i+1) = dd_TipPos_s[1] * ds_t * ds_t;

                Leg::LegIK(normalPee+2*i,normalPin+2*i,1);

                Leg::LegIJ(normalPee+2*i,jacobi,1);
                matrix_dot_matrix(jacobi,2,2,normalVee+2*i,1,normalVin+2*i);

                Leg::LegIdJ(normalPee+2*i,d_jacobi_x,d_jacobi_y,1);
                for(int j = 0; j < 4; j++)
                {
                    d_jacobi[j] = (d_jacobi_x[j] * d_TipPos_s[0] + d_jacobi_y[j] * d_TipPos_s[1]) * ds_t;
                }
                double tmp1[2];
                double tmp2[2];
                matrix_dot_matrix(d_jacobi,2,2,normalVee+2*i,1,tmp1);
                matrix_dot_matrix(jacobi,2,2,normalAee+2*i,1,tmp2);
                for (int j=0;j<2;j++)
                {
                    *(normalAin+2*i+j)=tmp1[j]+tmp2[j];
                }

                for(int j=0;j<2;j++)
                {
                    *(abs_Ain+2*i+j)=fabs(*(normalAin+2*i+j));
                    *(abs_Vin+2*i+j)=fabs(*(normalVin+2*i+j));
                }
            }

            maxAin=*std::max_element(abs_Ain,abs_Ain+2*normalTotalCount);
            maxVin=*std::max_element(abs_Vin,abs_Vin+2*normalTotalCount);

            if(maxAin<aLmt && maxVin<vLmt)
            {
                if(k==0)
                {
                    printf("How impossible! NormalGait is faster than OptimalGait! maxAin = %.4f < aLmt = %.4f\n",maxAin,aLmt);
//                    dlmwrite("./log/normalPee.txt",normalPee,normalTotalCount,2);
//                    dlmwrite("./log/normalVee.txt",normalVee,normalTotalCount,2);
//                    dlmwrite("./log/normalAee.txt",normalAee,normalTotalCount,2);
//                    dlmwrite("./log/normalPin.txt",normalPin,normalTotalCount,2);
//                    dlmwrite("./log/normalVin.txt",normalVin,normalTotalCount,2);
//                    dlmwrite("./log/normalAin.txt",normalAin,normalTotalCount,2);
                    break;
                }
                else
                {
                    stopFlag=true;
                    printf("Ain=%f or Vin=%f reach the maximum at %d-th iteration\n",maxAin,maxVin,k);
                    dlmwrite("./log/normalPee.txt",normalPee,normalTotalCount,2);
                    dlmwrite("./log/normalVee.txt",normalVee,normalTotalCount,2);
                    dlmwrite("./log/normalPin.txt",normalPin,normalTotalCount,2);
                    dlmwrite("./log/normalVin.txt",normalVin,normalTotalCount,2);
                    dlmwrite("./log/normalAin.txt",normalAin,normalTotalCount,2);
                }
            }
            else
            {
                normalTotalCount++;
                k++;
            }

            delete [] normalPee;
            delete [] normalVee;
            delete [] normalAee;
            delete [] normalPin;
            delete [] normalVin;
            delete [] normalAin;
            delete [] abs_Ain;
            delete [] abs_Vin;
        }
        printf("Finish GetNormalGait\n");
    }

    void TimeOptimalMotionSingleEffector::GetEntireGait()
    {
        printf("Start GetEntireGait\n");
        double c3;
        double c2;
        double c1;
        double pEB;

        int gait_num=2;
        double * entirePee=new double [2 * 2*gait_num * totalCount];
        double * entirePin=new double [2 * 2*gait_num * totalCount];
        double * Pin_const=new double [2 * 2 * totalCount];

        //acc for totalCount
        c3=(-v0 + stepD/2/(totalCount*0.001)) / (totalCount*0.001) / (totalCount*0.001);
        c2=(-v0 - 3*(-v0+stepD/2/(totalCount*0.001))) / (2*totalCount*0.001);
        for(int i=0;i<totalCount;i++)
        {
            pEB = c3*1e-9*i*i*i + c2*1e-6*i*i;//0 --> -D/4

            *(entirePee+2*i) = initTipPos[0] + stepD/4 * cos(PI/2 * (1 - cos(PI*i/totalCount))) - stepD/4 - pEB;
            *(entirePee+2*i+1) = initTipPos[1] + stepH * sin(PI/2 * (1 - cos(PI*i/totalCount)));

            Leg::LegIK(entirePee+2*i,entirePin+2*i,1);
        }

        //const for 2*gait_num*totalCount
        dlmread("./log/real_Pin.txt",Pin_const);
        for(int i=0;i<gait_num-1;i++)
        {
            for(int j=0;j<2*totalCount;j++)
            {
                int k = 0;
                if(j<totalCount)
                {
                    k = totalCount + i*2*totalCount + j+totalCount;
                }
                else
                {
                    k = totalCount + i*2*totalCount + j-totalCount;
                }
                *(entirePin+2*k) = *(Pin_const+2*j);
                *(entirePin+2*k+1) = *(Pin_const+2*j+1);

                Leg::LegFK(entirePin+2*k,entirePee+2*k,1);
            }
        }

        //dec for totalCount
        c1=-vt;
        c2=(-3*stepD/4-2*c1*(totalCount*0.001))/(totalCount*0.001)/(totalCount*0.001);
        c3=(stepD/2+c1*(totalCount*0.001))/(totalCount*0.001)/(totalCount*0.001)/(totalCount*0.001);
        for(int i=0;i<totalCount;i++)
        {
            pEB = c3*1e-9*i*i*i + c2*1e-6*i*i + c1*1e-3*i;// 0 --> -D/4

            int j = i + (2*gait_num-1)*totalCount;
            *(entirePee+2*j) = initTipPos[0] - stepD/4 - pEB;
            *(entirePee+2*j+1) = initTipPos[1];

            Leg::LegIK(entirePee+2*j,entirePin+2*j,1);
        }

        dlmwrite("./log/entirePee.txt",entirePee,2*gait_num * totalCount,2);
        dlmwrite("./log/entirePin.txt",entirePin,2*gait_num * totalCount,2);

        delete [] entirePee;
        delete [] entirePin;
        delete [] Pin_const;
        printf("Finish GetEntireGait\n");
    }

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

    void dlmread(const char *FileName, double *pMatrix)
    {
        std::ifstream file;

        file.open(FileName);

        if (!file) throw std::logic_error("file not exist");

        int i = 0;
        while (!file.eof())
        {
            file >> *(pMatrix + i);
            ++i;
        }
    }
}
