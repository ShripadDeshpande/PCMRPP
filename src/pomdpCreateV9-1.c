//#define MDP 1
//#define VERBOSE_STATE_NAMES 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "pomdpV9-1.h"

enum {action_MF, action_TL45,action_TL90,action_TR45,action_TR90,action_ST}action_d;
const char *action_str[ACTION_MAX]={"ac_MF","ac_TL45","ac_TL90","ac_TR45","ac_TR90","ac_ST"};

const char *alpha_str[ALPHA_MAX]={"al_0","al_45","al_90","al_135","al_180","al_225","al_270","al_315"};

const char *theta_str[THETA_MAX]={"th_0","th_45","th_90","th_135","th_180","th_225","th_270","th_315"};

enum {r1_S, r1_N, r1_F}r1speed_d;
const char *r1speed_str[]={"r1_S", "r1_N", "r1_F"};

enum { r2_S, r2_N, r2_F}r2speed_d;
const char *r2speed_str[]={"r2_S", "r2_N", "r2_F"};

enum {prox_N, prox_P, prox_C}prox_d;
const char *prox_str[]={"pr_N","pr_P", "px_C"};

int max_ob;
int findalpha(void),findtheta(void),findprox(void);
int descretise(int ang);
int alpha, theta, r1speed, r2speed, prox, action, dist;
int alpha_descr, theta_descr, r1speed_descr, r2speed_descr, prox_descr, action, dist;
int newalpha, newtheta, newr1speed, newr2speed, newprox, action, dist;
int newalpha_descr, newtheta_descr, newr1speed_descr, newr2speed_descr, newprox_descr;
float r;
double r1x=50,r1y=50,r2x,r2y;
double targetx=50,targety=100;
int target_dist,newtarget_dist, delta_target_dist;
float target_angle, newtarget_angle;
double newr1x=50,newr1y=50,newr2x,newr2y;
int c;
int states[PROXIMITY_MAX][ALPHA_MAX][THETA_MAX][R1SPEED_MAX][R2SPEED_MAX];
int observations[PROXIMITY_MAX][ALPHA_MAX][THETA_MAX][R1SPEED_MAX][R2SPEED_MAX];
int T_matrix[ACTION_MAX][STATES_MAX][STATES_MAX];
int O_matrix[ACTION_MAX][STATES_MAX];
int R_matrix[ACTION_MAX][STATES_MAX][STATES_MAX];
int ss,es,ob,tc,oc,rc;
int T_mat_cnt[ACTION_MAX][STATES_MAX][MAX_ES_CNT]; // For each action, for each start state, which all end states are reached
int state_cnt[STATES_MAX][5];  //List of all the states in the system. The dimension 5 is for the 5 components of the state formulation - prox, alpha,theta,r1speed and r2speed
int tmp_ObsProb[ACTION_MAX][STATES_MAX];
int tmp_TransProb[ACTION_MAX][STATES_MAX];

FILE *fout,*ftmp0,*ftmp1,*ftmp2;
void write_preamble(void);
void read_TMatrix(void);
void sortTMatCnt(int a, int s);
void writeObsProb(int act,int es,int ss);

char start_str[50],end_str[50],rwrd_str[50],obs_str[50];
int ii,a,s,it,tot_es_cnt;
int prob;
int tot_prob_till_now ;
float tmp_flt;
int x[20];
double tot_iterations = ((50/PROX_INC)+1)*(360/AL_THE_INC)*(360/AL_THE_INC)*((12/SPEED_INC)+1)*((12/SPEED_INC)+1)*ACTION_MAX ;
int itPerbar = (((50/PROX_INC)+1)*(360/AL_THE_INC)*(360/AL_THE_INC)*((12/SPEED_INC)+1)*((12/SPEED_INC)+1)*ACTION_MAX )/50;
long iterations_done=0;
void parse_arg(int argc, char *argv[]);
int verbosity=5;
char filestr[50];



int main(int argc, char **argv)
{
    double time_spent = 0.0;
    int tmp_state,max_ii=0;
    int positive_alpha,positive_theta;
    clock_t begin,end;
    begin=clock();
///////////////////////////////////  init module - start /////////////////////////////////////////////////////////////////
    memset(states,0,sizeof states);
    memset(observations,0,sizeof observations);
    memset(T_matrix,0,sizeof T_matrix);
    memset(O_matrix,0,sizeof O_matrix);
    memset(R_matrix,0,sizeof R_matrix);
    memset(T_mat_cnt,0,sizeof T_mat_cnt);
    memset(state_cnt,0,sizeof state_cnt);
    memset(tmp_ObsProb,0,sizeof tmp_ObsProb);
    memset(tmp_TransProb,0,sizeof tmp_TransProb);
    max_ob=MAX_OBS; // default
    sprintf(filestr,"../../%sOBS%02d.pomdp","mrs",MAX_OBS);
    verbosity=DEFAULT_VERBOSITY ;
    parse_arg(argc,argv);
    fout=fopen(filestr,"w+");
    if(fout == NULL)
    {
        fprintf(stderr,"\nCan not open file %s\n",filestr);
        exit(-1);
    }
    ftmp1=fopen("mrs_Obs.pomdp","w+");
    if(ftmp1 == NULL)
    {
        fprintf(stderr,"\nCan not open file mrs_Obs.pomdp\n");
        fclose(fout);
        exit(-1);
    }

    ftmp2=fopen("mrs_Rwrd.pomdp","w+");
    if(ftmp2 == NULL)
    {
        fprintf(stderr,"\nCan not open file mrs_Rwrd.pomdp\n");
        fclose(fout);
        fclose(ftmp1);
        exit(-1);
    }
    write_preamble();
    fflush(fout);
//////////////////////////////////////////// Init module -- end //////////////////////////////////////////////
    printf("\nPreamble written. Writing T, O and R matrices ...\n");
    target_dist=sqrt(((targety-r1y)*(targety-r1y))+ ((targetx-r1x)*(targetx-r1x)));
    target_angle=atan2((targety-r1y),(targetx-r1x));
    ss= states[prox_N][0][0][r1_N][r2_N] ;
#ifdef VERBOSE_STATE_NAMES
    fprintf(fout,"start: st%04d-%s-%s-%s-%s-%s \n",ss,prox_str[prox_N],alpha_str[0],theta_str[0],r1speed_str[r1_N],r2speed_str[r2_N]);
#else
    fprintf(fout,"start: st%04d \n",ss);
#endif
    it=0;
    fflush(stdout);
//////////////////////////////////////// Main Iteration Module -- Start ////////////////////////////////////////////////
    for(prox=0;prox<=50;prox+=5){
        for(alpha=0;alpha<360;alpha+=1){
            // decide the position of r2
            dist=2*R_SIZE+prox;
            r2x=r1x+dist*cos((alpha+90)*PI/180);
            r2y=r1y+dist*sin((alpha+90)*PI/180);
            for(theta=0;theta<360;theta+=1){
                for(r1speed=0;r1speed<=12;r1speed+=4){
                    for(r2speed=0;r2speed<=12;r2speed+=4){
                        newr2x = r2x+r2speed*cos(theta*PI/180);
                        newr2y = r2y+r2speed*sin(theta*PI/180);
                        for(action=0;action<ACTION_MAX;action++){
                            switch (action) {
                            case action_MF :
                                newr1y = r1y+r1speed ; // r1x does not change
                                newr1x = r1x ;
                                newalpha = findalpha();
                                newtheta = findtheta();
                                newprox = findprox();
                                // no change in r1speed and r2speed
                                break ;
                            case action_TL45 :
                                // No change in r1x, r1y
                                newr1y = r1y;
                                newr1x = r1x ;
                                newalpha = findalpha();
                                newalpha -= 45 ; // since r1 turned left by 45 degrees
                                newtheta = findtheta();
                                newtheta -= 45 ;
                                newprox = findprox();
                                // no change in r1speed and r2speed
                                break ;
                            case action_TL90 :
                                // No change in r1x, r1y
                                newr1y = r1y;
                                newr1x = r1x ;
                                newalpha = findalpha();
                                newalpha -= 90 ; // since r1 turned left by 90 degrees
                                newtheta = findtheta();
                                newtheta -= 90 ;
                                newprox = findprox();
                                // no change in r1speed and r2speed
                                break ;
                            case action_TR45 :
                                // No change in r1x, r1y
                                newr1y = r1y;
                                newr1x = r1x ;
                                newalpha = findalpha();
                                newalpha += 45 ; // since r1 turned right by 45 degrees
                                newtheta = findtheta();
                                newtheta += 45 ;
                                newprox = findprox();
                                // no change in r1speed and r2speed
                                break ;
                            case action_TR90 :
                                // No change in r1x, r1y
                                newr1y = r1y;
                                newr1x = r1x ;
                                newalpha = findalpha();
                                newalpha += 90 ; // since r1 turned right by 90 degrees
                                newtheta = findtheta();
                                newtheta += 90 ;
                                newprox = findprox();
                                // no change in r1speed and r2speed
                                break ;
                            case action_ST :
                                // No change in r1x, r1y
                                newr1y = r1y;
                                newr1x = r1x ;
                                newalpha = findalpha();
                                newtheta = findtheta() ;
                                newprox = findprox();
                                // no change in r1speed and r2speed
                                break ;
                            }
                                // no change in r1speed and r2speed
                            newr1speed=r1speed;
                            newr2speed=r2speed;
                            newtarget_dist=sqrt(((targety-newr1y)*(targety-newr1y))+ ((targetx-newr1x)*(targetx-newr1x)));
                            newtarget_angle=atan2((targety-newr1y),(targetx-newr1x));
                            if(newtarget_angle < 0)newtarget_angle += 2*PI ;
                            delta_target_dist=newtarget_dist-target_dist;
                            positive_alpha=alpha;
                            if(positive_alpha<0)positive_alpha+=360;
                            alpha_descr=descretise(positive_alpha);
                            positive_theta=theta;
                            if(positive_theta<0)positive_theta+=360;
                            theta_descr=descretise(positive_theta);
                            if(r1speed<=2)r1speed_descr=r1_S;
                            else if(r1speed <= 7)r1speed_descr=r1_N;
                            else r1speed_descr=r1_F ;
                            //discretise r2speed
                            if(r2speed<=2)r2speed_descr=r2_S;
                            else if(r2speed <= 7)r2speed_descr=r2_N;
                            else r2speed_descr=r2_F ;
                            //discretise proximity
                            if(prox<0)prox_descr=prox_C;
                            else if(prox < 5)prox_descr=prox_P;
                            else prox_descr=prox_N ;

                            //discretize all state components of new state

                            //discretise newalpha
                            if(newalpha<0)newalpha+=360;
                            newalpha_descr=descretise(newalpha);
                            //discretise newtheta
                            if(newtheta<0)newtheta+=360;
                            newtheta_descr=descretise(newtheta);
                            //discretise newr1speed
                            if(newr1speed<=R_SPEED_S)newr1speed_descr=r1_S;
                            else if(newr1speed < R_SPEED_M)newr1speed_descr=r1_N;
                            else newr1speed_descr=r1_F ;
                            //discretise newr2speed
                            if(newr2speed<=R_SPEED_S)newr2speed_descr=r2_S;
                            else if(newr2speed <= R_SPEED_M)newr2speed_descr=r2_N;
                            else newr2speed_descr=r2_F ;
                            //discretise proximity
                            if(newprox<0)newprox_descr=prox_C;
                            else if(newprox < R_SIZE)newprox_descr=prox_P;
                            else newprox_descr=prox_N ;

                            if(prox_descr==prox_C){
                                    ss=STATES_MAX-1 ; // Collision state
#ifdef VERBOSE_STATE_NAMES
                                    sprintf(start_str,"st%04d-pr_C ",ss);
#else
                                    sprintf(start_str,"st%04d ",ss);
#endif
                            }
                            else {
                                    ss=states[prox_descr][alpha_descr][theta_descr][r1speed_descr][r2speed_descr]; //decide start state sequence no
#ifdef VERBOSE_STATE_NAMES
                                    sprintf(start_str,"st%04d-%s-%s-%s-%s-%s ",ss,prox_str[prox_descr],alpha_str[alpha_descr],theta_str[theta_descr],r1speed_str[r1speed_descr],r2speed_str[r2speed_descr]);
#else
                                    sprintf(start_str,"st%04d ",ss);
#endif
                            }
                            if(newprox_descr==prox_C){
                                    es=STATES_MAX-1 ; // Collision state
#ifdef VERBOSE_STATE_NAMES
                                    sprintf(end_str,"st%04d-pr_C ",STATES_MAX-1);
#else
                                    sprintf(end_str,"st%04d ",STATES_MAX-1);
#endif
                                    strcpy(obs_str,end_str);obs_str[0]='o';obs_str[1]='b';
                            }
                            else {
                                    es=states[newprox_descr][newalpha_descr][newtheta_descr][newr1speed_descr][newr2speed_descr]; // decide end state sequence no
#ifdef VERBOSE_STATE_NAMES
                                    sprintf(end_str,"st%04d-%s-%s-%s-%s-%s ",es,prox_str[newprox_descr],alpha_str[newalpha_descr],theta_str[newtheta_descr],r1speed_str[newr1speed_descr],r2speed_str[newr2speed_descr]);
#else
                                    sprintf(end_str,"st%04d ",es);
#endif
                                    strcpy(obs_str,end_str);obs_str[0]='o';obs_str[1]='b';
                            }
                            ob=es;
                            //output State transition probabilities -
                            T_matrix[action][ss][es]+=1; // Count how many such transitions. That will decide probability
                            if(T_matrix[action][ss][es]==1){
                                //this T entry is first time for this es so record the es
                                for(ii=0;ii<MAX_ES_CNT;ii++){
                                    // find empty slot at the end
                                    if(T_mat_cnt[action][ss][ii])continue;
                                    // record the es. es is stored with an offset because es can genuinely be 0
                                    T_mat_cnt[action][ss][ii]=es+STATE_STORAGE_OFFSET;
                                    if(ii>max_ii)
                                         max_ii=ii;
                                    break;
                                }
                                if(ii >= MAX_ES_CNT){
                                    fprintf(stderr,"\nEnd State count %d for start state %d, Exceeded limit. Increase MAX_ES_CNT. Aborting ...",ii,ss);
                                    exit(-2);
                                }
                                tc++;
                            }
////////////////////////////////////////  Rewards Computation Module -- Start ///////////////////////////////////
                            //output rewards-
                            r=0;
                            if(!R_matrix[action][ss][es]){
                                ////////  Rewards ////////////////////////////////
                                if(newprox_descr == prox_C){r=-5;x[0]++;}
                                else if(prox_descr==prox_N && newprox_descr==prox_P){r=-2;x[1]++;}
                                else if(prox_descr==prox_C && newprox_descr==prox_P){r=2;x[2]++;}   //never happens (?)
                                else if(prox_descr==prox_C && newprox_descr==prox_N){r=3;x[3]++;}   //never happens (?)
                                else if(prox_descr==prox_P && newprox_descr==prox_N){r=2;x[4]++;}   //never happens (?)
                                else{
                                    if(newtarget_dist < target_dist){r=5;x[5]++;} // Compare with the original target distance.
                                    if(delta_target_dist<0){r=5;x[6]++;}  // has it reduced in this iteration
                                }
                                R_matrix[action][ss][es]=r ;
                                rc++;
                            }
////////////////////////////////////////  Rewards Computation Module -- End ///////////////////////////////////
                            iterations_done++;
                        }
//                        if(verbosity==5){
//                            if(iterations_done > progressCnt*itPerbar){
//                                printf("O");
//                                progressCnt++;
//                            }
//                        }
                        if(verbosity==4)
                            printf("\rIterations %08ld Proximity=%02d alpha=%03d theta=%03d r1Speed=%03d r2Speed=%03d",iterations_done,prox,alpha,theta,r1speed,r2speed);
                    }
                    if(verbosity==3)
                        printf("\rIterations %08ld Proximity=%02d alpha=%03d theta=%03d r1Speed=%03d",iterations_done,prox,alpha,theta,r1speed);
                }
                if(verbosity==2)
                    printf("\rIterations %08ld Proximity=%02d alpha=%03d theta=%03d",iterations_done,prox,alpha,theta);
            }
            if(verbosity==1)
                printf("\rProx=%02d alpha=%03d",prox,alpha);
        }
        if(verbosity==0)
           printf("\rProx %d",prox);
//        printf("\n");
    }
//////////////////////////////////////// Main Iteration Module -- End ////////////////////////////////////////////////

//////////////////////////////////////// State Transition Probability Computation Module -- Start ////////////////////////////////////////////////

    memset(tmp_ObsProb,0,sizeof(tmp_ObsProb));
    for(a=0;a<ACTION_MAX;a++){
        for(s=0;s<STATES_MAX;s++){
            tot_es_cnt=0;
            if(s==STATES_MAX-1)
#ifdef VERBOSE_STATE_NAMES
                sprintf(start_str,"st%04d-pr_C ",STATES_MAX-1);
#else
                sprintf(start_str,"st%04d ",STATES_MAX-1);
#endif
            else
#ifdef VERBOSE_STATE_NAMES
                sprintf(start_str,"st%04d-%s-%s-%s-%s-%s ",s,prox_str[state_cnt[s][0]],alpha_str[state_cnt[s][1]],theta_str[state_cnt[s][2]],r1speed_str[state_cnt[s][3]],r2speed_str[state_cnt[s][4]]);
#else
                sprintf(start_str,"st%04d ",s);
#endif
            tot_prob_till_now = 0; //for this ss
            for(ii=0;ii<MAX_ES_CNT;ii++){
                    // find all es for which this ss has transition
                if(T_mat_cnt[a][s][ii]){
                    tot_es_cnt += T_matrix[a][s][T_mat_cnt[a][s][ii]-STATE_STORAGE_OFFSET]; // es was stored with offset 1
                    continue;
                }
                else break;
            }
            sortTMatCnt(a,s);
            for(ii=0;;ii++){
                if(tot_prob_till_now >= 100) break;
                if(T_mat_cnt[a][s][ii]){
                    prob = floor((T_matrix[a][s][T_mat_cnt[a][s][ii]-STATE_STORAGE_OFFSET]*100.0/tot_es_cnt)+0.5);
                    if(!T_mat_cnt[a][s][ii+1]) // if this is last end state for current start state then make probability 100-whatever is left
                        prob = 100-tot_prob_till_now ;
                    tmp_state = T_mat_cnt[a][s][ii]-STATE_STORAGE_OFFSET ;
                    if(tmp_state==(STATES_MAX-1)){
#ifdef VERBOSE_STATE_NAMES
                        sprintf(end_str,"st%04d-pr_C ",STATES_MAX-1);
#else
                        sprintf(end_str,"st%04d ",STATES_MAX-1);
#endif
                    }
                    else{
#ifdef VERBOSE_STATE_NAMES
                        sprintf(end_str,"st%04d-%s-%s-%s-%s-%s ",tmp_state,prox_str[state_cnt[tmp_state][0]],alpha_str[state_cnt[tmp_state][1]],theta_str[state_cnt[tmp_state][2]],r1speed_str[state_cnt[tmp_state][3]],r2speed_str[state_cnt[tmp_state][4]]);
#else
                        sprintf(end_str,"st%04d ",tmp_state);
#endif
                    }
                    strcpy(obs_str,end_str);obs_str[0]='o';obs_str[1]='b';
                    if(tot_prob_till_now+prob >= 100)
                        prob = 100-tot_prob_till_now ;
                    tot_prob_till_now += prob ;
                }
                else {
                    if(tot_prob_till_now<100)
                    {
                        prob+=(100-tot_prob_till_now);
                        if(!tot_prob_till_now){ //if probability is zero, it means there is no transition from this start state so the state remains same
                            strcpy(end_str,start_str);
                            strcpy(obs_str,end_str);obs_str[0]='o';obs_str[1]='b';
                        }
                    }
                    else
                    {
                        //Total prob more than or equal to 100% !!
                        prob-=(tot_prob_till_now-100); // assuming this number does not turn out to be negative.
                    }
                }
                if(!T_mat_cnt[a][s][ii])
                    break;
                if(prob>=100){
                    tmp_TransProb[a][s]=100;
                    fprintf(fout,"\nT: %s : %s : %s 1.00",action_str[a],start_str,end_str) ;
                }else{
                    tmp_TransProb[a][s]+=prob;
                    fprintf(fout,"\nT: %s : %s : %s 0.%02d",action_str[a],start_str,end_str,prob) ;
                }
//////////////////////////////////////// State Transition Probability Computation Module -- End ////////////////////////////////////////////////

//////////////////////////////////////////////// Observation Probability Computation module -- Start //////////////////////////////////////////
#ifndef MDP
                if(!O_matrix[a][tmp_state]){
                    writeObsProb(a,tmp_state,s);
                    O_matrix[a][tmp_state]=1;
                    oc++;
                }
#endif
                fflush(fout);
                fflush(ftmp1);
            }
#ifdef MDP
            fprintf(ftmp1,"\nO: %d : %d : %d 1.00",a,s,s) ;
#endif
        }   // state loop
    } // action loop
#ifndef MDP
    for(a=0;a<ACTION_MAX;a++){
        for(s=0;s<STATES_MAX;s++){
            if(!O_matrix[a][s]){
                O_matrix[a][s]=1;
                fprintf(ftmp1,"\nO: %d : %d : %d 1.00",a,s,s) ;
                tmp_ObsProb[a][s]=100;
                oc++;
            }
            if(tmp_ObsProb[a][s]!=100)
                if(verbosity>3)
                    printf("\nObservation Probability O for action %d and state %04d is %d. Does not add to 100",a,s,tmp_ObsProb[a][s]);
            if(tmp_TransProb[a][s]!=100)
                if(verbosity>3)
                if(verbosity>3)
                    printf("\nTransition Probability T for action %d and state %04d is %d. Does not add to 100",a,s,tmp_TransProb[a][s]);
        }
    }
    fflush(ftmp1);
#endif
    //special case for collision state
#ifdef VERBOSE_STATE_NAMES
    sprintf(start_str,"st%04d-pr_C ",STATES_MAX-1);
#else
    sprintf(start_str,"st%04d ",STATES_MAX-1);
#endif
    strcpy(end_str,start_str);
    fprintf(fout,"\nT: * : %s : %s 1.00",start_str,end_str) ;
//////////////////////////////////////////////// Observation Probability Computation module -- End //////////////////////////////////////////

    //Copy O matrix to output file
    fseek(ftmp1,0l,SEEK_SET);
    do {
      c = fgetc (ftmp1);
      if (c == EOF)break;
      fputc(c,fout);
    } while (c != EOF);
    fclose(ftmp1);
    //Copy R matrix to output file
    fseek(ftmp2,0l,SEEK_SET);
    do {
      c = fgetc (ftmp2);
      if (c == EOF)break;
      fputc(c,fout);
    } while (c != EOF);
    fclose(ftmp2);

    printf("\nPOMDP File %s ready.\n",filestr);
    if(verbosity==4){
        printf("\nTransition to collision state %d times",x[0]);
        printf("\nTransition from collision state to proximity %d times",x[1]);
        printf("\nTransition from collision state to no proximity %d times",x[2]);
        printf("\nTransition from proximity to no proximity state %d times",x[3]);
        printf("\nDistance to target reduced %d times",x[4]);
        printf("\nDelta distance to target reduced %d times",x[5]);
        printf("\nAngle to target reduced %d times",x[6]);
        printf("\nMaximum no. end states reached from one single start state %d times",max_ii);
    }
    end=clock();
    time_spent += (double)(end-begin)/(CLOCKS_PER_SEC) ;
    printf("\nTime taken : %f seconds\n\n",time_spent);
    return 0;
}

///////////////////////  Function writeObsProb() -- part of Observation Probability Computation module ///////////////////
void writeObsProb(int act,int es,int ss)
{
    int subtot[3*3*3*3*3],tmp_sortIndex[3*3*3*3*3];
    int ob_index[3*3*3*3*3];
    int x,tot_n,total;
    int prob_Current,prob_sofar;
    int prob[(PROXIMITY_MAX+1)+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+R2SPEED_MAX];
    int p,a,t,r1,r2,aplus1,aminus1,tplus1,tminus1;
    int i0,i1,i2,i3,i4;
    char tmpStr[20];
    p=state_cnt[ss][0] ;
    a=state_cnt[ss][1];
    t=state_cnt[ss][2];
    r1=state_cnt[ss][3];
    r2=state_cnt[ss][4];
#ifdef VERBOSE_STATE_NAMES
    sprintf(start_str,"st%04d-%s-%s-%s-%s-%s ",ss,prox_str[p],alpha_str[a],theta_str[t],r1speed_str[r1],r2speed_str[r2]);
#else
    sprintf(start_str,"st%04d ",ss);
#endif
    p=state_cnt[es][0] ;
    a=state_cnt[es][1];
    aplus1= (a+1)%ALPHA_MAX;
    aminus1= (a-1+ALPHA_MAX)%ALPHA_MAX;
    t=state_cnt[es][2];
    tplus1= (t+1)%THETA_MAX;
    tminus1= (t-1+THETA_MAX)%THETA_MAX;
    r1=state_cnt[es][3];
    r2=state_cnt[es][4];
    memset(prob,0,sizeof(prob));
    memset(subtot,0,sizeof(subtot));
    memset(tmp_sortIndex,0,sizeof(tmp_sortIndex));
    memset(ob_index,0,sizeof(ob_index));
    switch(p)
    {
    case 0: prob[0]=19;prob[1]=1;prob[2]=0;break;
    case 1: prob[0]=1;prob[1]=18;prob[2]=1;break;
    case 2: prob[0]=0;prob[1]=1;prob[2]=19;break;
    }
    //alpha
    prob[PROXIMITY_MAX+1+a]=16;
    prob[PROXIMITY_MAX+1+aplus1]=2;
    prob[PROXIMITY_MAX+1+aminus1]=2;
    //theta
    prob[PROXIMITY_MAX+1+ALPHA_MAX+t]=16;
    prob[PROXIMITY_MAX+1+ALPHA_MAX+tplus1]=2;
    prob[PROXIMITY_MAX+1+ALPHA_MAX+tminus1]=2;
    switch(r1)
    {
    case 0:
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX]=19;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+1]=1;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+2]=0;
        break;
    case 1:
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX]=1;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+1]=18;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+2]=1;
        break;
    case 2:
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX]=0;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+1]=1;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+2]=19;
        break;
    }
    switch(r2)
    {
    case 0:
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX]=19;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+1]=1;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+2]=0;
        break;

    case 1:
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX]=1;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+1]=18;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+2]=1;
        break;
    case 2:
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX]=0;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+1]=1;
        prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+2]=19;
        break;
    }
    x=0;
    total=0;
    for(i0=0;i0<=PROXIMITY_MAX;i0++){
        if(!prob[i0])continue;
        for(i1=0;i1<ALPHA_MAX;i1++){
            if(!prob[PROXIMITY_MAX+1+i1])continue;
            for(i2=0;i2<THETA_MAX;i2++){
                if(!prob[PROXIMITY_MAX+1+ALPHA_MAX+i2])continue;
                for(i3=0;i3<R1SPEED_MAX;i3++){
                    if(!prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+i3])continue;
                    for(i4=0;i4<R2SPEED_MAX;i4++){
                        if(!prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+i4])continue;
                        subtot[x]=prob[i0]+prob[PROXIMITY_MAX+1+i1]+prob[PROXIMITY_MAX+1+ALPHA_MAX+i2]+prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+i3]+prob[PROXIMITY_MAX+1+ALPHA_MAX+THETA_MAX+R1SPEED_MAX+i4];
                        if(i0==PROXIMITY_MAX)
                            ob_index[x]=STATES_MAX-1;
                        else
                            ob_index[x]=states[i0][i1][i2][i3][i4]; //decide state sequence no
                        tmp_sortIndex[x]=x; // using for sorting subtot below
                        x++;
                    }
                }
            }
        }
    }
    tot_n = x;
    //collate all entries of STATES_MAX into one
    i1=-1;
    for(i0=0;i0<tot_n;i0++){
        if(ob_index[i0]==STATES_MAX-1){
            if(i1==-1)
                i1=i0; //This is first entry
            else{
                // this is not the first entry so collate with first entry
                subtot[i1] += subtot[i0] ;
                subtot[i0]=0;
                ob_index[i0]=0;
            }
        }
    }

    //sort subtot array, store sorted indexes in the tmp_sortIndex array.
    i1=1; // for swap
    while(i1){
        i1=0;
        for(i0=0;i0<tot_n-1;i0++){
            if(subtot[tmp_sortIndex[i0]]<subtot[tmp_sortIndex[i0+1]]){
                i1=1;
                i2=tmp_sortIndex[i0];
                tmp_sortIndex[i0] = tmp_sortIndex[i0+1] ;
                tmp_sortIndex[i0+1]=i2;
            }
        }
    }

    total=0;
    for(i0=0;i0<max_ob && i0<tot_n;i0++){
        total += subtot[tmp_sortIndex[i0]] ;
    }
    while(i0<tot_n){
        subtot[tmp_sortIndex[i0++]] =0;
    }
    x=0;
    prob_sofar=0.0;
    if(p==2)
#ifdef VERBOSE_STATE_NAMES
        sprintf(end_str,"st%04d-pr_C ",STATES_MAX-1);
#else
        sprintf(end_str,"st%04d ",STATES_MAX-1);
#endif
    else
#ifdef VERBOSE_STATE_NAMES
        sprintf(end_str,"st%04d-%s-%s-%s-%s-%s ",es,prox_str[p],alpha_str[a],theta_str[t],r1speed_str[r1],r2speed_str[r2]);
#else
        sprintf(end_str,"st%04d ",es);
#endif
    for(i0=0;i0<max_ob;i0++){
#ifdef VERBOSE_STATE_NAMES
        printf("\nConditional Compilation will cause trouble at line %d in file %s\n",__LINE__.__FILE__);
        //Here the individual components of state formulation are required for verbose name construction of obs string
#else
        sprintf(obs_str,"ob%04d ",ob_index[tmp_sortIndex[i0]]);
#endif
        if(!subtot[tmp_sortIndex[i0+1]] || i0 == max_ob-1)
            prob_Current=100-prob_sofar ;
        else
            prob_Current=(int)(((float)subtot[tmp_sortIndex[i0]]/total*100)+0.5);
        prob_sofar += prob_Current;
        if(prob_Current>0){
            if(prob_Current==100)
                sprintf(tmpStr,"1.00");
            else
                sprintf(tmpStr,"0.%02d",prob_Current);
            tmp_ObsProb[act][es]+=prob_Current;
            fprintf(ftmp1,"\nO: %s : %s : %s %s",action_str[act],end_str,obs_str,tmpStr);
            fprintf(ftmp2,"\nR: %s : %s : %s : %s %4.2f ",action_str[act],start_str,end_str,obs_str,(float)R_matrix[act][ss][es]);
        }
    }
    if(prob[2] && prob_sofar<100){
#ifdef VERBOSE_STATE_NAMES
        sprintf(obs_str,"ob%04d-pr_C ",STATES_MAX-1);
#else
        sprintf(obs_str,"ob%04d ",STATES_MAX-1);
#endif
        tmp_ObsProb[act][es]+=100-prob_sofar;
        fprintf(ftmp1,"\nO: %s : %s : %s 0.%02d",action_str[act],end_str,obs_str,(100-prob_sofar));
//        printf("\nO: %s : %s : %s 0.%02d",action_str[act],end_str,obs_str,(100-prob_sofar));
        fprintf(ftmp2,"\nR: %s : %s : %s : %s %4.2f",action_str[act],start_str,end_str,obs_str,(float)(((100-prob_sofar)*R_matrix[act][ss][es])/100));
    }
}
/////////////// Functions sortTMatCnt(), descretise(), findalpha(), findtheta(), findprox() are part of utility module ///////////////////////
void sortTMatCnt(int a, int s)
{
    int *ptr = &T_mat_cnt[a][s][0];
    int max,maxIndex,tmp,i,j;
    int cnt=0;
    while(*ptr++)cnt++;
    if(!cnt)return;
    ptr=&T_mat_cnt[a][s][0];
    i=0;
    for(i=0;i<cnt-1;i++){
        max=T_matrix[a][s][ptr[i]-STATE_STORAGE_OFFSET];
        maxIndex=i;
        for (j=i+1;j<cnt;j++){
            if(max<T_matrix[a][s][ptr[j]-STATE_STORAGE_OFFSET]){
                max=T_matrix[a][s][ptr[j]-STATE_STORAGE_OFFSET] ;
                maxIndex = j;
            }
        }
        tmp=T_mat_cnt[a][s][i] ;
        T_mat_cnt[a][s][i]=T_mat_cnt[a][s][maxIndex];
        T_mat_cnt[a][s][maxIndex]=tmp;
    }
}

int descretise(int ang)
{
    int i,z, d_val=0;
    z=360/(2*N_DIV);
    for(i=0;i<N_DIV;i++)
    {
        if(ang<z)
        {
            d_val=i;
            break;
        }
        z+=(360/N_DIV);
    }
    return d_val;
}

int findalpha(void)
{
        double x,x_diff, y_diff ;
        x_diff = newr2x-newr1x ;
        y_diff = newr2y-newr1y ;
        if(x_diff < 0.0001 && x_diff > -0.0001)
        {
            x = 90 ;
            if(y_diff < 0) x=-x ;
        }
        else
        {
            x = atan2(y_diff,x_diff);
            x = x*180/PI ;
        }
        x-=90; // this is because alpha is measured w.r.t.y axis
        return((int)x);
}
int findtheta(void)
{
    return (theta) ; // presently R2 assumed to be going with constant speed and direction.
}
int findprox(void)
{
    double dist,d1,d2 ;
    d1 = newr2y-newr1y ;
    d2 = newr2x-newr1x ;
    dist = sqrt(d1*d1+d2*d2);
    return ((int)(dist-2*R_SIZE));
}
///////////////////////////// Function write_preamble() is part of init module
void write_preamble(void)
{
    int prox,alpha,theta,r1s,r2s,act,cnt=0;
    printf("\nwriting discount and values ...");
    fprintf(fout,"discount: 0.95\n");
    fprintf(fout,"values: reward\n");

    printf("\nwriting states ...\n");
    fprintf(fout,"\nstates: ");
//    fprintf(fout,"\n");
    for(prox=0;prox<PROXIMITY_MAX;prox++){
        for(alpha=0;alpha<ALPHA_MAX;alpha++){
            for(theta=0;theta<THETA_MAX;theta++){
                for(r1s=0;r1s<R1SPEED_MAX;r1s++){
                    for(r2s=0;r2s<R2SPEED_MAX;r2s++){
#ifdef VERBOSE_STATE_NAMES
                        fprintf(fout,"st%04d-%s-%s-%s-%s-%s ",cnt,prox_str[prox],alpha_str[alpha],theta_str[theta],r1speed_str[r1s],r2speed_str[r2s]);
#else
                        fprintf(fout,"st%04d ",cnt);
#endif
                        state_cnt[cnt][0]=prox;state_cnt[cnt][1]=alpha;state_cnt[cnt][2]=theta;state_cnt[cnt][3]=r1s;state_cnt[cnt][4]=r2s;
                        states[prox][alpha][theta][r1s][r2s]=cnt++;
                        printf("\r%04d",cnt);
//                        if((cnt%4)==0)fprintf(fout,"\n");
                    }
                }
            }
        }
    }
    // Special for collision case
    state_cnt[cnt][0]=PROXIMITY_MAX;
    states[prox][0][0][0][0]=cnt;

#ifdef VERBOSE_STATE_NAMES
    fprintf(fout,"st%04d-pr_C ",cnt);
#else
    fprintf(fout,"st%04d ",cnt);
#endif
    printf("\nwriting actions ...\n");
    fprintf(fout,"\nactions: ");
    for(act=0;act<ACTION_MAX;act++){
#ifdef VERBOSE_STATE_NAMES
        fprintf(fout,"%s ",action_str[act]);
#else
        fprintf(fout,"%s ",action_str[act]);
#endif
//        actions[act]=act;
    }
    printf("\nwriting observations ...\n");
    fprintf(fout,"\n\nobservations: ");
    cnt=0;
    for(prox=0;prox<PROXIMITY_MAX;prox++){
        for(alpha=0;alpha<ALPHA_MAX;alpha++){
            for(theta=0;theta<THETA_MAX;theta++){
                for(r1s=0;r1s<R1SPEED_MAX;r1s++){
                    for(r2s=0;r2s<R2SPEED_MAX;r2s++){
#ifdef VERBOSE_STATE_NAMES
                        fprintf(fout,"ob%04d-%s-%s-%s-%s-%s ",cnt,prox_str[prox],alpha_str[alpha],theta_str[theta],r1speed_str[r1s],r2speed_str[r2s]);
#else
                        fprintf(fout,"ob%04d ",cnt);
#endif
                        observations[prox][alpha][theta][r1s][r2s]=cnt++;
                        printf("\r%d",cnt);
//                        if((cnt%4)==0)fprintf(fout,"\n");
                    }
                }
            }
        }
    }
#ifdef VERBOSE_STATE_NAMES
    fprintf(fout,"ob%04d-pr_C ",cnt);
#else
    fprintf(fout,"ob%04d ",cnt);
#endif
    fprintf(fout,"\n");
}

