#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "pomdpV9-1.h"

extern int verbosity;
extern char filestr[20];
char *help_str = " [OPTION...]\nCommand Line options :\n\t-v verbosity Level : 0...4 setting verbosity level of messages\n\t-o FILE : output pomdp file name without the extension .pomdp\n\t-h : this help";
void wrong_option(char *opt, char *cmd)
{
   printf("\nInvalid option %s.\nEnter %s -h for help\n\n",opt,cmd);
   exit(-1);
}

void parse_arg(int argc, char *argv[])
{
    int i;
    if(argc < 2){
        // No options specified. defaults take over
        return;
    }
    for(i=1;i<argc;i++){
        if(argv[i][0]=='-'){
            switch(argv[i][1]){
            case 'v':
                if(strlen(argv[i])==2){
                    if(argv[i+1] != NULL){
                        if(strlen(argv[i+1])==1 && argv[i+1][0]>='0' && argv[i+1][0]<='4'){
                            verbosity=argv[i+1][0]-'0';
                            printf("\nVerbosity=%d\n",verbosity);
                            ++i;
                            continue;
                        }
                        else{
                            printf("\nInvalid parameter %s for option -v\nEnter %s -h for help\n\n",argv[i+1],argv[0]);
                            exit(-1);
                        }
                    }
                    else{
                        printf("\nNo parameter specified for option %s\nEnter %s -h for help\n\n",argv[i],argv[0]);
                        exit(-1);
                    }
                }
                else
                    wrong_option(argv[i],argv[0]);
                break;
            case 'o':
                if(strlen(argv[i])==2)
                   if(argv[i+1] != NULL){
                        strcpy(filestr,argv[i+1]);
                        printf("\nOutput file name %s\n",filestr);
                        ++i;
                        continue;
                    }
                    else
                        printf("\nNo output file specified. Taking default (mrs.pomdp)\n");
                else
                    wrong_option(argv[i],argv[0]);
                break;
            case 'h':
                if(strlen(argv[i])==2){
                    printf("%s %s\n",argv[0],help_str);
                    exit(0);
                }
                else
                    wrong_option(argv[i],argv[0]);
                break;
            default:
                wrong_option(argv[i],argv[0]);
            }
        }
        else{
            printf("\nInvalid argument %s. Enter %s -h for help\n\n",argv[i], argv[0]);
            exit(-1);
        }

    }// for loop
}
