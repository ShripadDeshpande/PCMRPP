#ifndef POMDPV9_1_H_INCLUDED
#define POMDPV9_1_H_INCLUDED

#define N_DIV 8
#define ALPHA_MAX N_DIV
#define THETA_MAX N_DIV
#define R1SPEED_MAX 3
#define R2SPEED_MAX 3
#define PROXIMITY_MAX 2
#define ACTION_MAX 6
#define STATES_MAX ((PROXIMITY_MAX*ALPHA_MAX*THETA_MAX*R1SPEED_MAX*R2SPEED_MAX)+1)
#define OBSERVATIONS_MAX STATES_MAX

#define PROX_INC 5
#define AL_THE_INC 1
#define SPEED_INC 4

#define MAXNODES 1340
#define MAXBUFF 2400
#define MAX_OBS 10

#define MAX_X 620
#define MIN_X 20
#define MAX_Y 460
#define MIN_Y 20

#define rad2deg(x) ((x)*180/(PI))
#define deg2rad(x) ((x)*(PI)/180)

#define DEFAULT_VERBOSITY 0
#define R_SPEED_S 3
#define R_SPEED_M 7
#define R_SPEED_F 12
#define R_SIZE 15

#define MAX_ES_CNT 320    // This is max number of end states that a start state can transform to. Keep it at least one more that max
#define MAX_ITERATIONS 200
#define PI 3.141592654
#define TOLERANCE 0.00001
#define COA 0.2 // 0.0000100000001   // coefficient of attraction
#define COR 0.1   // coefficient of repulsion
#define STEPSIZE 1
#define STATE_STORAGE_OFFSET 1000


#endif // POMDPV9_1_H_INCLUDED


