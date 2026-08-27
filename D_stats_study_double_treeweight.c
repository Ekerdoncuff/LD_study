//
//  D_stats_study_double_treeweight.c
//  
//
//  Created by elise.kerdoncuff on 09/10/2024. (First version)
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <malloc/malloc.h> //Need to modify to include <malloc.h> for Linux
#include <unistd.h>

#define MAX(a,b) (a>=b?a:b)
#define MIN(a,b) (a<=b?a:b)

#include "D_etude_lib.h"


int main(int argc, char* argv[])
{
 
    //options
    int optc;
    extern char * optarg;
    extern int optind;
    
    //default value
    int N=10;
    double TCHANG = 0.0;
    float CONTRAC = 1.0;
    int ITE = 1;
    int n=1;
    int report_type_recombi=0;

    int i;
    int tree=0; //to not print newick tree
    int SEED=0;
    int stats_D=0;
    
    
    char *outfile;
    
    while( (optc=getopt(argc,argv, "N:T:C:A:S:D:n:R:")) != -1 ){
        
        switch(optc)
        {
                
            case 'N': //number of leaves
                 N=atoi(optarg);
                 printf("N = %d\n",N);
                 break;
                
            case 'T': //TCHANG time of contraction
                TCHANG=atof(optarg);
                printf("TCHANG = %1.9lf\n",TCHANG);
                break;
                
            case 'C': // strenght of contraction
                CONTRAC=atof(optarg);
                printf("CONTRAC = %f\n",CONTRAC);
                break;
                
            case 'A': //  newick tree (0) not print, (1) print
                tree=atoi(optarg);
                printf("Tree = %d\n",tree);
                break;
                
            case 'S': //SEED
                SEED=atoi(optarg);
                printf("SEED = %d\n",SEED);
                break;
                
            case 'D': //To know which stat to print D0:0 D1:1 Dinf:8
                stats_D=atoi(optarg);
                printf("stats_D = %d\n",stats_D);
                break;
                
            case 'n': // number of iterations
                n=atoi(optarg);
                printf("n = %d\n",n);
                break;
            
            case 'R': // If 1 compute per recombination type
                report_type_recombi=atoi(optarg);
                printf("report_type_recombi = %d\n",report_type_recombi);
                break;
                
            default:
                printf("unknown option\n");
                exit(1);
        }
    }
    
    if (SEED==0)
    {
        srand(time(NULL));

    }
    else
    {
        srand(SEED);
    }
    
    //allocation global file
    struct Linkage * tab_D_global= (struct Linkage *)malloc(n*((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    //only for th estimation of var D0
    float th_varD0, var_th_D0_global=0;
    
    ////To stock information for different type of recombination
    int nb_recombi_par_type[4]={0,0,0,0};
    double weight_tree_recombi_par_type[4]={0,0,0,0};
    int type_recombi;

   /// To store configuration probs
    
    double * tab_CP0_global = (double *)malloc(n*((N*2-2)*(N*2-2))*sizeof(double));
    double * tab_CP1_global = (double *)malloc(n*((N*2-2)*(N*2-2))*sizeof(double));
    double * tab_CP2_global = (double *)malloc(n*((N*2-2)*(N*2-2))*sizeof(double));
    double * tab_CP3_global = (double *)malloc(n*((N*2-2)*(N*2-2))*sizeof(double));
    
    
       // Allocate CP arrays for each tree type (tr0, tr1, tr2, tr3)
    double *tab_CP0_global_tr0 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP1_global_tr0 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP2_global_tr0 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP3_global_tr0 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));

    double *tab_CP0_global_tr1 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP1_global_tr1 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP2_global_tr1 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP3_global_tr1 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    
    double *tab_CP0_global_tr2 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP1_global_tr2 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP2_global_tr2 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP3_global_tr2 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));

    double *tab_CP0_global_tr3 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP1_global_tr3 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP2_global_tr3 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));
    double *tab_CP3_global_tr3 = (double *)malloc(n * ((N*2-2) * (N*2-2)) * sizeof(double));

    
    struct Linkage * tab_D_global_tr0= (struct Linkage *)malloc(n*((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    struct Linkage * tab_D_global_tr1= (struct Linkage *)malloc(n*((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    struct Linkage * tab_D_global_tr2= (struct Linkage *)malloc(n*((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    struct Linkage * tab_D_global_tr3= (struct Linkage *)malloc(n*((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    
    
    if(argc-optind<1)fprintf(stderr,"usage is %s [opt] outfile, bye\n",argv[0]),exit(3);
    outfile=argv[optind];
    
    printf("seed %.0f\n", (double) time(NULL));
    
    printf("N= %d\n", N);

    //How many tree are computed:
    if (stats_D==0)
        ITE=1;
    else
        ITE=2;
    
    
    //iterations
    int ite_measure;
    for (ite_measure=0;ite_measure<n;ite_measure++)
    {
        ////
        // Initialise values
        ////
        int ite=0; //iteration compte for number of blocks
        //To compare trees
        float Ltot1,Ltot2; //total lenght
        int root1,root2; //root ID

        
        //Trees info
        struct InfosIndiv * story=NULL; //all informations about the tree, see
        story=(struct InfosIndiv *)malloc((N*2-1)*sizeof(struct InfosIndiv)); //sizeof 16
        if (story==NULL) printf("Memory issue\n");
        struct InfosIndiv * story_init;
        story_init=(struct InfosIndiv *)malloc((N*2-1)*sizeof(struct InfosIndiv));
        
        
        
        // Internal node descendant tables for partition comparison
        int * desc_node_inter_tree1 = NULL;
        desc_node_inter_tree1=(int *)malloc(N*sizeof(int));
        int ** tab_desc_node_inter_tree1;
        tab_desc_node_inter_tree1=(int **)malloc((N-1)*sizeof(desc_node_inter_tree1));
        
        int * desc_node_inter_tree2 = NULL;
        desc_node_inter_tree2=(int *)malloc(N*sizeof(int));


        int ** tab_desc_node_inter_tree2;
        tab_desc_node_inter_tree2=(int **)malloc((N-1)*sizeof(desc_node_inter_tree2));
        

        
        
        while (ite<ITE){
            
            ////
            // Tree layout
            ////
            if (ite ==0) // Tree for chromosome start
            {
            
                story=Init_ARG(story,N, TCHANG,CONTRAC);
                Ltot1=sumlbs_To(story,N);
                root1=find_root(story,N);
                
                //For D calculation between two trees
                for (i=0;i<N*2-1;i++)
                {
                    story_init[i]=story[i];
                }
            
                ////
                // Writing the tree in NEWICK format
                ////
                if (tree!=0)
                {
                    char * strtree; //to write the string
                    strtree=(char*)malloc(1000*(sizeof(char)));
                    if (strtree==NULL)
                    {
                        fprintf(stdout, "I cannot printtree"), exit(35);
                    }
                    strtree[999]=0;
                    int pos=0;
                    strtree=PrintTree(0,story,N,strtree,&pos);
                    strtree[pos]=0;
                    printf("%s;\n", strtree);
                    
                    free(strtree);
                }
                
            }
            else{
                
                
                // Internal node descendant tables for partition comparison
                int parcours,writing=0;
                for (parcours=N;parcours<N*2-1;parcours++)
                {

                        int * table = NULL;
                        table=(int *)calloc(N,sizeof(int));
                        table=leavesof(table,parcours,story);
                        tab_desc_node_inter_tree1[writing]=table;
                        writing++;
                    
                }
                if (stats_D==1)
                    story=SMC_tree(story,N,TCHANG,CONTRAC);
                else if (stats_D==8)
                    story=Init_ARG(story,N, TCHANG,CONTRAC);
                else
                    printf("Problem!! not possible\n"), exit(180);
                Ltot2=sumlbs_To(story,N);
                root2=find_root(story,N);
                
                
                // Internal node descendant tables for partition comparison
                writing=0;
                for (parcours=N;parcours<N*2-1;parcours++)
                {
                        int * table = NULL;
                        table=(int *)calloc(N,sizeof(int));
                        table=leavesof(table,parcours,story);
                        tab_desc_node_inter_tree2[writing]=table;
                        writing++;

                }
                ////
                // Recombination type
                ////
                //0:Invisible 1:Silent 2:Discrete 3:Incompatible
                type_recombi=-1;
                // Tree comparison (bi-partiton) between before and after recombination
                int similarite_arbres;
                if (ite!=0)
                {
                    similarite_arbres=cmp_tabdetabe(tab_desc_node_inter_tree1,tab_desc_node_inter_tree2,N-1,N);
                    if (similarite_arbres==0)
                    {
                        if (root1==root2)
                        {
                            if (Ltot1==Ltot2)
                                type_recombi=0;
                            else
                                type_recombi=1;
                        }
                    }
                    else if (similarite_arbres==10)
                        type_recombi=2;
                    else if (similarite_arbres==11)
                        type_recombi=3;
                
                    nb_recombi_par_type[type_recombi]++;
                    if (stats_D==0 ||stats_D==1)
                    {
                        weight_tree_recombi_par_type[type_recombi]+=Ltot1;
                    }else if (stats_D==8)
                    {
                        weight_tree_recombi_par_type[type_recombi]+=Ltot1*Ltot2;
                    }
                }
            
                
                //To continue comparison
                Ltot1=Ltot2;
                root1=root2;
                
             
            }
            
            
            
            
            ////
            // Writing tree in NEWICK
            ////
            if (tree!=0)
            {
                char * strtree; //where we are writing the caracters
                strtree=(char*)malloc(1000*(sizeof(char)));
                if (strtree==NULL)
                {
                    fprintf(stdout, "I cannot print"), exit(35);
                }
                strtree[999]=0;
                int pos=0;
                strtree=PrintTree(0,story,N,strtree,&pos);
                strtree[pos]=0;
                printf("%s;\n", strtree);
                
                free(strtree);
            }
        
        
            
            if (stats_D==0 || ite==1)
            {

                ////
                // Linkage Desequilibrium
                /////
                struct Linkage * tab_D= NULL;
                //For one tree, values and frequencies
                if (stats_D==0)
                    tab_D=compute_D_ontree(story,N);
                else
                    tab_D=compute_D_twotrees(story_init,story,N);
                
                //Theoretical value of var(D) from tree (if one tree)
                if (stats_D==0)
                {
                    th_varD0=calcul_varD0(story,N);
                    var_th_D0_global+=th_varD0/n;
                }
                
                
                //Study per tree
                //Description
                for (i=0;i<(N*2-2)*(N*2-2);i++)
                {

                    //For global study: storage
                    tab_D_global[i+ite_measure*(N*2-2)*(N*2-2)]=tab_D[i];
                    
                    
                    // Fill in config probs
                    if (tab_D[i].r==1){ // perfect pos
                        tab_CP0_global[i+ite_measure*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                        tab_CP1_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP2_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP3_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                    }else if (tab_D[i].r>=0){ // partial pos
                        tab_CP0_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP1_global[i+ite_measure*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                        tab_CP2_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP3_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                    }else if (tab_D[i].r==-1){ // perfect neg
                        tab_CP0_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP1_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP2_global[i+ite_measure*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                        tab_CP3_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                    }else if (tab_D[i].r<0){ // partial neg
                        tab_CP0_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP1_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP2_global[i+ite_measure*(N*2-2)*(N*2-2)]=0;
                        tab_CP3_global[i+ite_measure*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                    }
                    else{
                        printf("PB with assigning CP, r value %f\n",tab_D[i].r);
                    }
                        
 
                    
                    if (report_type_recombi==1){
                        //for each recombination type
                        switch(type_recombi)
                        {
                                
                            case 0:
                                tab_D_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=tab_D[i];
                                // Fill in config probs
                                
                                //Flag Dp!=1
                                if (tab_D[i].d_p<0.99999)
                                {
                                    printf("\n!!!!!!!!!!!\nPB Dp : %f\n!!!!!!!!!!!!!\n",tab_D[i].d_p);
                                }
                                
                                
                                if (tab_D[i].r==1){ // perfect pos
                                    tab_CP0_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP1_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r>=0){ // partial pos
                                    tab_CP0_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP2_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r==-1){ // perfect neg
                                    tab_CP0_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP3_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r<0){ // partial neg
                                    tab_CP0_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr0[i+(nb_recombi_par_type[0]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                }

                                break;
                                
                            case 1:
                                tab_D_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=tab_D[i];
                                                                // Fill in config probs
                                
                                //Flag Dp!=1
                                if (tab_D[i].d_p<0.99999)
                                {
                                    printf("\n!!!!!!!!!!!\nPB Dp : %f\n!!!!!!!!!!!!!\n",tab_D[i].d_p);
                                }
                                
                                if (tab_D[i].r==1){ // perfect pos
                                    tab_CP0_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP1_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r>=0){ // partial pos
                                    tab_CP0_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP2_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r==-1){ // perfect neg
                                    tab_CP0_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP3_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r<0){ // partial neg
                                    tab_CP0_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr1[i+(nb_recombi_par_type[1]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                }
                                break;
                                
                            case 2:
                                tab_D_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=tab_D[i];
                                // Fill in config probs
                                
                                //Flag Dp!=1
                                if (tab_D[i].d_p<0.99999)
                                {
                                    printf("\n!!!!!!!!!!!\nPB Dp : %f\n!!!!!!!!!!!!!\n",tab_D[i].d_p);
                                }
                                
                                if (tab_D[i].r==1){ // perfect pos
                                    tab_CP0_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP1_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r>=0){ // partial pos
                                    tab_CP0_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP2_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r==-1){ // perfect neg
                                    tab_CP0_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP3_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r<0){ // partial neg
                                    tab_CP0_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr2[i+(nb_recombi_par_type[2]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                }
                                
                                
                                break;
                                
                            case 3:
                                tab_D_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=tab_D[i];
                                
                                if (tab_D[i].r==1){ // perfect pos
                                    tab_CP0_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP1_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r>=0){ // partial pos
                                    tab_CP0_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP2_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r==-1){ // perfect neg
                                    tab_CP0_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                    tab_CP3_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                }else if (tab_D[i].r<0){ // partial neg
                                    tab_CP0_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP1_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP2_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=0;
                                    tab_CP3_global_tr3[i+(nb_recombi_par_type[3]-1)*(N*2-2)*(N*2-2)]=tab_D[i].frequency;
                                }

                                break;
                            
                                
                            default:
                                printf("unknown recombi %d\n",type_recombi);
                                exit(1);
                        }
                    }
                }
                
                //printf("moy D %f\tmoy Dp %f\tmoy r %f\n",mean_D,mean_D_p,mean_r2);
                
                
                free(tab_D);
            }
            
            ite++;
        }
        
    

        free(story);
    
    }
    
    //Global study
    //Description
    //mean
    double mean_D_global_treeweight=0,mean_D_p_global_treeweight=0,mean_r_global_treeweight=0,mean_r2_global_treeweight=0,mean_D2_global_treeweight=0,mean_pqA_pqB_global_treeweight=0,mean_CP0_global_treeweight=0,mean_CP1_global_treeweight=0,mean_CP2_global_treeweight=0,mean_CP3_global_treeweight=0, mean_treeweight=0;
    for (i=0;i<n*(N*2-2)*(N*2-2);i++)
    {
        
        if (stats_D==0||stats_D==1)
        {
            mean_D_global_treeweight+=tab_D_global[i].d*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A;
            mean_D_p_global_treeweight+=tab_D_global[i].d_p*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A;
            mean_r_global_treeweight+=tab_D_global[i].r*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A;
            mean_r2_global_treeweight+=tab_D_global[i].r2*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A;
            mean_D2_global_treeweight+=(tab_D_global[i].d*tab_D_global[i].d*tab_D_global[i].frequency)*tab_D_global[i].tot_lbranch_A;
            mean_pqA_pqB_global_treeweight+=tab_D_global[i].pqA_pqB*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A;
            
            mean_CP0_global_treeweight+=tab_CP0_global[i]*tab_D_global[i].tot_lbranch_A;
            mean_CP1_global_treeweight+=tab_CP1_global[i]*tab_D_global[i].tot_lbranch_A;
            mean_CP2_global_treeweight+=tab_CP2_global[i]*tab_D_global[i].tot_lbranch_A;
            mean_CP3_global_treeweight+=tab_CP3_global[i]*tab_D_global[i].tot_lbranch_A;
            
            mean_treeweight+=tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A;
            
        }else if (stats_D==8)
        {
            mean_D_global_treeweight+=tab_D_global[i].d*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            mean_D_p_global_treeweight+=tab_D_global[i].d_p*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            mean_r_global_treeweight+=tab_D_global[i].r*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            mean_r2_global_treeweight+=tab_D_global[i].r2*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            mean_D2_global_treeweight+=(tab_D_global[i].d*tab_D_global[i].d*tab_D_global[i].frequency)*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            mean_pqA_pqB_global_treeweight+=tab_D_global[i].pqA_pqB*tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            
            mean_CP0_global_treeweight+=tab_CP0_global[i]*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            mean_CP1_global_treeweight+=tab_CP1_global[i]*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            mean_CP2_global_treeweight+=tab_CP2_global[i]*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            mean_CP3_global_treeweight+=tab_CP3_global[i]*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
            
            mean_treeweight+=tab_D_global[i].frequency*tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
        }
        
        

         
    }
    
    mean_D_global_treeweight=mean_D_global_treeweight/mean_treeweight;
    mean_D_p_global_treeweight=mean_D_p_global_treeweight/mean_treeweight;
    mean_r_global_treeweight=mean_r_global_treeweight/mean_treeweight;
    mean_r2_global_treeweight=mean_r2_global_treeweight/mean_treeweight;
    mean_D2_global_treeweight=mean_D2_global_treeweight/mean_treeweight;
    mean_pqA_pqB_global_treeweight=mean_pqA_pqB_global_treeweight/mean_treeweight;
    
    mean_CP0_global_treeweight=mean_CP0_global_treeweight/mean_treeweight;
    mean_CP1_global_treeweight=mean_CP1_global_treeweight/mean_treeweight;
    mean_CP2_global_treeweight=mean_CP2_global_treeweight/mean_treeweight;
    mean_CP3_global_treeweight=mean_CP3_global_treeweight/mean_treeweight;
    
    
    double mean_sD_global_treeweight=mean_D2_global_treeweight/mean_pqA_pqB_global_treeweight;
    
    
    //variance
    double var_D_global_treeweight=0,var_D_p_global_treeweight=0,var_r_global_treeweight=0,var_r2_global_treeweight=0,var_sD_global_treeweight=0,var_CP0_global_treeweight=0,var_CP1_global_treeweight=0,var_CP2_global_treeweight=0,var_CP3_global_treeweight=0;
    for (i=0;i<n*(N*2-2)*(N*2-2);i++)
    {
        double tree_i_weight=0;
        
        if (stats_D==0||stats_D==1)
        {
            tree_i_weight=tab_D_global[i].tot_lbranch_A;
        }
        else if (stats_D==8)
        {
            tree_i_weight=tab_D_global[i].tot_lbranch_A*tab_D_global[i].tot_lbranch_B;
        }
        
        var_D_global_treeweight+=pow(tab_D_global[i].d-mean_D_global_treeweight,2)*tab_D_global[i].frequency*tree_i_weight;
        var_D_p_global_treeweight+=pow(tab_D_global[i].d_p-mean_D_p_global_treeweight,2)*tab_D_global[i].frequency*tree_i_weight;
        var_r_global_treeweight+=pow(tab_D_global[i].r-mean_r_global_treeweight,2)*tab_D_global[i].frequency*tree_i_weight;
        var_r2_global_treeweight+=pow(tab_D_global[i].r2-mean_r2_global_treeweight,2)*tab_D_global[i].frequency*tree_i_weight;
        var_sD_global_treeweight+=pow(tab_D_global[i].d/tab_D_global[i].pqA_pqB -mean_sD_global_treeweight,2)*tab_D_global[i].frequency*tree_i_weight;
        
        var_CP0_global_treeweight+=pow(tab_CP0_global[i]-mean_CP0_global_treeweight,2)*tab_CP0_global[i]*tree_i_weight;
        var_CP1_global_treeweight+=pow(tab_CP1_global[i]-mean_CP1_global_treeweight,2)*tab_CP1_global[i]*tree_i_weight;
        var_CP2_global_treeweight+=pow(tab_CP2_global[i]-mean_CP2_global_treeweight,2)*tab_CP2_global[i]*tree_i_weight;
        var_CP3_global_treeweight+=pow(tab_CP3_global[i]-mean_CP3_global_treeweight,2)*tab_CP3_global[i]*tree_i_weight;
        
        
    }
    
    
    var_D_global_treeweight=var_D_global_treeweight/mean_treeweight;
    var_D_p_global_treeweight=var_D_p_global_treeweight/mean_treeweight;
    var_r_global_treeweight=var_r_global_treeweight/mean_treeweight;
    var_r2_global_treeweight=var_r2_global_treeweight/mean_treeweight;
    var_sD_global_treeweight=var_sD_global_treeweight/mean_treeweight;
    
    var_CP0_global_treeweight= var_CP0_global_treeweight/mean_treeweight;
    var_CP1_global_treeweight=var_CP1_global_treeweight/mean_treeweight;
    var_CP2_global_treeweight=var_CP2_global_treeweight/mean_treeweight;
    var_CP3_global_treeweight=var_CP3_global_treeweight/mean_treeweight;

    
    if (stats_D==0)
    {
        printf("var th D0 %f\n",var_th_D0_global);
    }

    ////
    // Writing in output FILE
    ////
    
    FILE * fichier_mean_var=NULL;
    fichier_mean_var=fopen(outfile,"a");
    if (fichier_mean_var!=NULL)
    {
        fprintf(fichier_mean_var,"%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",stats_D,N,n,TCHANG,CONTRAC,mean_D_global_treeweight,var_D_global_treeweight,mean_D_p_global_treeweight,var_D_p_global_treeweight,mean_r_global_treeweight,var_r_global_treeweight,mean_r2_global_treeweight,var_r2_global_treeweight,mean_sD_global_treeweight,var_sD_global_treeweight,mean_CP0_global_treeweight,var_CP0_global_treeweight,mean_CP1_global_treeweight,var_CP1_global_treeweight,mean_CP2_global_treeweight,var_CP2_global_treeweight,mean_CP3_global_treeweight,var_CP3_global_treeweight);
        
        fclose(fichier_mean_var);
    }
    else
        printf("Problem writing/openning outfile\n");
    
    
    if (report_type_recombi==1)
    {
        //Type recombination
        FILE * fichier_tr=NULL;
        char* type_recombi_char= "type_recombi_";
        char* type_recombi_outfile;
        type_recombi_outfile = malloc(strlen(type_recombi_char)+1+strlen(outfile));
        strcpy(type_recombi_outfile, type_recombi_char);
        strcat(type_recombi_outfile,outfile);
        fichier_tr=fopen(type_recombi_outfile,"a");
        if (fichier_tr!=NULL)
        {
            fprintf(fichier_tr,"%d\t%f\t%f\t%f\t%f\t%f\t%f\n",N,TCHANG,CONTRAC,(float)nb_recombi_par_type[0]/(float)n,(float)nb_recombi_par_type[1]/(float)n,(float)nb_recombi_par_type[2]/(float)n,(float)nb_recombi_par_type[3]/(float)n);
            
            fclose(fichier_tr);
        }
        else
            printf("Problem writing/openning recombination type outfile\n");
        
        free(type_recombi_outfile);
        
        //Compute mean and variance per recombination type
        //TR0
        //mean
        double mean_D_global_tr0_treeweight=0,mean_D_p_global_tr0_treeweight=0,mean_r_global_tr0_treeweight=0,mean_r2_global_tr0_treeweight=0,mean_CP0_global_tr0_treeweight=0,mean_CP1_global_tr0_treeweight=0,mean_CP2_global_tr0_treeweight=0,mean_CP3_global_tr0_treeweight=0;
        for (i=0;i<nb_recombi_par_type[0]*(N*2-2)*(N*2-2);i++)
        {
            float tree_i_weight=0;
            
            if (stats_D==0||stats_D==1)
            {
                tree_i_weight=tab_D_global_tr0[i].tot_lbranch_A;
            }
            else if (stats_D==8)
            {
                tree_i_weight=tab_D_global_tr0[i].tot_lbranch_A*tab_D_global_tr0[i].tot_lbranch_B;
            }
            
            
            mean_D_global_tr0_treeweight+=tab_D_global_tr0[i].d*tab_D_global_tr0[i].frequency*tree_i_weight;
            mean_D_p_global_tr0_treeweight+=tab_D_global_tr0[i].d_p*tab_D_global_tr0[i].frequency*tree_i_weight;
            mean_r_global_tr0_treeweight+=tab_D_global_tr0[i].r*tab_D_global_tr0[i].frequency*tree_i_weight;
            mean_r2_global_tr0_treeweight+=tab_D_global_tr0[i].r2*tab_D_global_tr0[i].frequency*tree_i_weight;
            mean_CP0_global_tr0_treeweight+=tab_CP0_global_tr0[i]*tree_i_weight;
            mean_CP1_global_tr0_treeweight+=tab_CP1_global_tr0[i]*tree_i_weight;
            mean_CP2_global_tr0_treeweight+=tab_CP2_global_tr0[i]*tree_i_weight;
            mean_CP3_global_tr0_treeweight+=tab_CP3_global_tr0[i]*tree_i_weight;
                        
        }
        
        mean_D_global_tr0_treeweight/=weight_tree_recombi_par_type[0];
        mean_D_p_global_tr0_treeweight/=weight_tree_recombi_par_type[0];
        mean_r_global_tr0_treeweight/=weight_tree_recombi_par_type[0];
        mean_r2_global_tr0_treeweight/=weight_tree_recombi_par_type[0];
        mean_CP0_global_tr0_treeweight/=weight_tree_recombi_par_type[0];
        mean_CP1_global_tr0_treeweight/=weight_tree_recombi_par_type[0];
        mean_CP2_global_tr0_treeweight/=weight_tree_recombi_par_type[0];
        mean_CP3_global_tr0_treeweight/=weight_tree_recombi_par_type[0];
        
        
        //variance
        double var_D_global_tr0_treeweight=0,var_D_p_global_tr0_treeweight=0,var_r_global_tr0_treeweight=0,var_r2_global_tr0_treeweight=0,var_CP0_global_tr0_treeweight=0,var_CP1_global_tr0_treeweight=0,var_CP2_global_tr0_treeweight=0,var_CP3_global_tr0_treeweight=0;
        for (i=0;i<nb_recombi_par_type[0]*(N*2-2)*(N*2-2);i++)
        {
            float tree_i_weight = 0;

                if (stats_D == 0 || stats_D == 1)
                {
                    tree_i_weight = tab_D_global_tr0[i].tot_lbranch_A;
                }
                else if (stats_D == 8)
                {
                    tree_i_weight = tab_D_global_tr0[i].tot_lbranch_A * tab_D_global_tr0[i].tot_lbranch_B;
                }
            
            
            var_D_global_tr0_treeweight += pow(tab_D_global_tr0[i].d - mean_D_global_tr0_treeweight, 2) * tab_D_global_tr0[i].frequency * tree_i_weight;
            var_D_p_global_tr0_treeweight += pow(tab_D_global_tr0[i].d_p - mean_D_p_global_tr0_treeweight, 2) * tab_D_global_tr0[i].frequency * tree_i_weight;
            var_r_global_tr0_treeweight += pow(tab_D_global_tr0[i].r - mean_r_global_tr0_treeweight, 2) * tab_D_global_tr0[i].frequency * tree_i_weight;
            var_r2_global_tr0_treeweight += pow(tab_D_global_tr0[i].r2 - mean_r2_global_tr0_treeweight, 2) * tab_D_global_tr0[i].frequency * tree_i_weight;

            var_CP0_global_tr0_treeweight += pow(tab_CP0_global_tr0[i] - mean_CP0_global_tr0_treeweight, 2) * tree_i_weight;
            var_CP1_global_tr0_treeweight += pow(tab_CP1_global_tr0[i] - mean_CP1_global_tr0_treeweight, 2) * tree_i_weight;
            var_CP2_global_tr0_treeweight += pow(tab_CP2_global_tr0[i] - mean_CP2_global_tr0_treeweight, 2) * tree_i_weight;
            var_CP3_global_tr0_treeweight += pow(tab_CP3_global_tr0[i] - mean_CP3_global_tr0_treeweight, 2) * tree_i_weight;
            
        
        }
        
        var_D_global_tr0_treeweight /= weight_tree_recombi_par_type[0];
        var_D_p_global_tr0_treeweight /= weight_tree_recombi_par_type[0];
        var_r_global_tr0_treeweight /= weight_tree_recombi_par_type[0];
        var_r2_global_tr0_treeweight /= weight_tree_recombi_par_type[0];

        var_CP0_global_tr0_treeweight /= weight_tree_recombi_par_type[0];
        var_CP1_global_tr0_treeweight /= weight_tree_recombi_par_type[0];
        var_CP2_global_tr0_treeweight /= weight_tree_recombi_par_type[0];
        var_CP3_global_tr0_treeweight /= weight_tree_recombi_par_type[0];
        
        
        FILE * fichier_mean_var_tr0=NULL;
        char* tr0_char= "recombi_type0_";
        char* tr0_outfile;
        tr0_outfile = malloc(strlen(tr0_char)+1+strlen(outfile));
        strcpy(tr0_outfile, tr0_char);
        strcat(tr0_outfile,outfile);
        fichier_mean_var_tr0=fopen(tr0_outfile,"a");
        if (fichier_mean_var_tr0!=NULL)
        {
            //printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",stats_D,0,N,nb_recombi_par_type[0],TCHANG,CONTRAC,mean_D_global_tr0_treeweight,var_D_global_tr0_treeweight,mean_D_p_global_tr0_treeweight,var_D_p_global_tr0_treeweight,mean_r_global_tr0_treeweight,var_r_global_tr0_treeweight,mean_r2_global_tr0_treeweight,var_r2_global_tr0_treeweight,mean_CP0_global_tr0_treeweight,var_CP0_global_tr0_treeweight,mean_CP1_global_tr0_treeweight,var_CP1_global_tr0_treeweight,mean_CP2_global_tr0_treeweight,var_CP2_global_tr0_treeweight,mean_CP3_global_tr0_treeweight,var_CP3_global_tr0_treeweight);
            fprintf(fichier_mean_var_tr0,"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",stats_D,0,N,nb_recombi_par_type[0],TCHANG,CONTRAC,mean_D_global_tr0_treeweight,var_D_global_tr0_treeweight,mean_D_p_global_tr0_treeweight,var_D_p_global_tr0_treeweight,mean_r_global_tr0_treeweight,var_r_global_tr0_treeweight,mean_r2_global_tr0_treeweight,var_r2_global_tr0_treeweight,mean_CP0_global_tr0_treeweight,var_CP0_global_tr0_treeweight,mean_CP1_global_tr0_treeweight,var_CP1_global_tr0_treeweight,mean_CP2_global_tr0_treeweight,var_CP2_global_tr0_treeweight,mean_CP3_global_tr0_treeweight,var_CP3_global_tr0_treeweight);
            
            fclose(fichier_mean_var_tr0);
        }
        else
            printf("Problem writing/openning per recombination type outfile\n");
        
        free(tr0_outfile);
        
        //TR1
        //mean
        double mean_D_global_tr1_treeweight=0,mean_D_p_global_tr1_treeweight=0,mean_r_global_tr1_treeweight=0,mean_r2_global_tr1_treeweight=0,mean_CP0_global_tr1_treeweight=0,mean_CP1_global_tr1_treeweight=0,mean_CP2_global_tr1_treeweight=0,mean_CP3_global_tr1_treeweight=0;

        for (i=0;i<nb_recombi_par_type[1]*(N*2-2)*(N*2-2);i++)
        {
            float tree_i_weight=0;
            
            if (stats_D==0||stats_D==1)
            {
                tree_i_weight=tab_D_global_tr1[i].tot_lbranch_A;
            }
            else if (stats_D==8)
            {
                tree_i_weight=tab_D_global_tr1[i].tot_lbranch_A*tab_D_global_tr1[i].tot_lbranch_B;
            }
            
            mean_D_global_tr1_treeweight+=tab_D_global_tr1[i].d*tab_D_global_tr1[i].frequency*tree_i_weight;
            mean_D_p_global_tr1_treeweight+=tab_D_global_tr1[i].d_p*tab_D_global_tr1[i].frequency*tree_i_weight;
            mean_r_global_tr1_treeweight+=tab_D_global_tr1[i].r*tab_D_global_tr1[i].frequency*tree_i_weight;
            mean_r2_global_tr1_treeweight+=tab_D_global_tr1[i].r2*tab_D_global_tr1[i].frequency*tree_i_weight;
            mean_CP0_global_tr1_treeweight+=tab_CP0_global_tr1[i]*tree_i_weight;
            mean_CP1_global_tr1_treeweight+=tab_CP1_global_tr1[i]*tree_i_weight;
            mean_CP2_global_tr1_treeweight+=tab_CP2_global_tr1[i]*tree_i_weight;
            mean_CP3_global_tr1_treeweight+=tab_CP3_global_tr1[i]*tree_i_weight;
            
        }
        
        
        mean_D_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        mean_D_p_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        mean_r_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        mean_r2_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        mean_CP0_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        mean_CP1_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        mean_CP2_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        mean_CP3_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        
        
        //variance
        double var_D_global_tr1_treeweight=0,var_D_p_global_tr1_treeweight=0,var_r_global_tr1_treeweight=0,var_r2_global_tr1_treeweight=0,var_CP0_global_tr1_treeweight=0,var_CP1_global_tr1_treeweight=0,var_CP2_global_tr1_treeweight=0,var_CP3_global_tr1_treeweight=0;
        for (i=0;i<nb_recombi_par_type[1]*(N*2-2)*(N*2-2);i++)
        {
            float tree_i_weight=0;
            
            if (stats_D==0||stats_D==1)
            {
                tree_i_weight=tab_D_global_tr1[i].tot_lbranch_A;
            }
            else if (stats_D==8)
            {
                tree_i_weight=tab_D_global_tr1[i].tot_lbranch_A*tab_D_global_tr1[i].tot_lbranch_B;
            }
            
            
            var_D_global_tr1_treeweight+=pow(tab_D_global_tr1[i].d-mean_D_global_tr1_treeweight,2)*tab_D_global_tr1[i].frequency*tree_i_weight;
            var_D_p_global_tr1_treeweight+=pow(tab_D_global_tr1[i].d_p-mean_D_p_global_tr1_treeweight,2)*tab_D_global_tr1[i].frequency*tree_i_weight;
            var_r_global_tr1_treeweight+=pow(tab_D_global_tr1[i].r-mean_r_global_tr1_treeweight,2)*tab_D_global_tr1[i].frequency*tree_i_weight;
            var_r2_global_tr1_treeweight+=pow(tab_D_global_tr1[i].r2-mean_r2_global_tr1_treeweight,2)*tab_D_global_tr1[i].frequency*tree_i_weight;
            var_CP0_global_tr1_treeweight+=pow(tab_CP0_global_tr1[i]-mean_CP0_global_tr1_treeweight,2)*tab_CP0_global_tr1[i]*tree_i_weight;
            var_CP1_global_tr1_treeweight+=pow(tab_CP1_global_tr1[i]-mean_CP1_global_tr1_treeweight,2)*tab_CP1_global_tr1[i]*tree_i_weight;
            var_CP2_global_tr1_treeweight+=pow(tab_CP2_global_tr1[i]-mean_CP2_global_tr1_treeweight,2)*tab_CP2_global_tr1[i]*tree_i_weight;
            var_CP3_global_tr1_treeweight+=pow(tab_CP3_global_tr1[i]-mean_CP3_global_tr1_treeweight,2)*tab_CP3_global_tr1[i]*tree_i_weight;
            
        }
        
        var_D_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        var_D_p_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        var_r_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        var_r2_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        var_CP0_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        var_CP1_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        var_CP2_global_tr1_treeweight/=weight_tree_recombi_par_type[1];
        var_CP3_global_tr1_treeweight/=weight_tree_recombi_par_type[1];

        FILE * fichier_mean_var_tr1=NULL;
        char* tr1_char= "recombi_type1_";
        char* tr1_outfile;
        tr1_outfile = malloc(strlen(tr1_char)+1+strlen(outfile));
        strcpy(tr1_outfile, tr1_char);
        strcat(tr1_outfile,outfile);
        fichier_mean_var_tr1=fopen(tr1_outfile,"a");
        if (fichier_mean_var_tr1!=NULL)
        {
            fprintf(fichier_mean_var_tr1,"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",stats_D,1,N,nb_recombi_par_type[1],TCHANG,CONTRAC,mean_D_global_tr1_treeweight,var_D_global_tr1_treeweight,mean_D_p_global_tr1_treeweight,var_D_p_global_tr1_treeweight,mean_r_global_tr1_treeweight,var_r_global_tr1_treeweight,mean_r2_global_tr1_treeweight,var_r2_global_tr1_treeweight,mean_CP0_global_tr1_treeweight,var_CP0_global_tr1_treeweight,mean_CP1_global_tr1_treeweight,var_CP1_global_tr1_treeweight,mean_CP2_global_tr1_treeweight,var_CP2_global_tr1_treeweight,mean_CP3_global_tr1_treeweight,var_CP3_global_tr1_treeweight);
            
            fclose(fichier_mean_var_tr1);
        }
        else
            printf("Problem writing/openning per recombination type outfile\n");
        
        free(tr1_outfile);
        
        //TR2
        //mean
        double mean_D_global_tr2_treeweight=0,mean_D_p_global_tr2_treeweight=0,mean_r_global_tr2_treeweight=0,mean_r2_global_tr2_treeweight=0,mean_CP0_global_tr2_treeweight=0,mean_CP1_global_tr2_treeweight=0,mean_CP2_global_tr2_treeweight=0,mean_CP3_global_tr2_treeweight=0;
        for (i=0;i<nb_recombi_par_type[2]*(N*2-2)*(N*2-2);i++)
        {
            float tree_i_weight=0;
            
            if (stats_D==0||stats_D==1)
            {
                tree_i_weight=tab_D_global_tr2[i].tot_lbranch_A;
            }
            else if (stats_D==8)
            {
                tree_i_weight=tab_D_global_tr2[i].tot_lbranch_A*tab_D_global_tr2[i].tot_lbranch_B;
            }
            
            
            mean_D_global_tr2_treeweight+=tab_D_global_tr2[i].d*tab_D_global_tr2[i].frequency*tree_i_weight;
            mean_D_p_global_tr2_treeweight+=tab_D_global_tr2[i].d_p*tab_D_global_tr2[i].frequency*tree_i_weight;
            mean_r_global_tr2_treeweight+=tab_D_global_tr2[i].r*tab_D_global_tr2[i].frequency*tree_i_weight;
            mean_r2_global_tr2_treeweight+=tab_D_global_tr2[i].r2*tab_D_global_tr2[i].frequency*tree_i_weight;
            mean_CP0_global_tr2_treeweight+=tab_CP0_global_tr2[i]*tree_i_weight;
            mean_CP1_global_tr2_treeweight+=tab_CP1_global_tr2[i]*tree_i_weight;
            mean_CP2_global_tr2_treeweight+=tab_CP2_global_tr2[i]*tree_i_weight;
            mean_CP3_global_tr2_treeweight+=tab_CP3_global_tr2[i]*tree_i_weight;
            
        }
                
        mean_D_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        mean_D_p_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        mean_r_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        mean_r2_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        mean_CP0_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        mean_CP1_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        mean_CP2_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        mean_CP3_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        
        //variance
        double var_D_global_tr2_treeweight=0,var_D_p_global_tr2_treeweight=0,var_r_global_tr2_treeweight=0,var_r2_global_tr2_treeweight=0,var_CP0_global_tr2_treeweight=0,var_CP1_global_tr2_treeweight=0,var_CP2_global_tr2_treeweight=0,var_CP3_global_tr2_treeweight=0;
        for (i=0;i<nb_recombi_par_type[2]*(N*2-2)*(N*2-2);i++)
        {
            float tree_i_weight=0;
            
            if (stats_D==0||stats_D==1)
            {
                tree_i_weight=tab_D_global_tr2[i].tot_lbranch_A;
            }
            else if (stats_D==8)
            {
                tree_i_weight=tab_D_global_tr2[i].tot_lbranch_A*tab_D_global_tr2[i].tot_lbranch_B;
            }
            
            var_D_global_tr2_treeweight+=pow(tab_D_global_tr2[i].d-mean_D_global_tr2_treeweight,2)*tab_D_global_tr2[i].frequency*tree_i_weight;
            var_D_p_global_tr2_treeweight+=pow(tab_D_global_tr2[i].d_p-mean_D_p_global_tr2_treeweight,2)*tab_D_global_tr2[i].frequency*tree_i_weight;
            var_r_global_tr2_treeweight+=pow(tab_D_global_tr2[i].r-mean_r_global_tr2_treeweight,2)*tab_D_global_tr2[i].frequency*tree_i_weight;
            var_r2_global_tr2_treeweight+=pow(tab_D_global_tr2[i].r2-mean_r2_global_tr2_treeweight,2)*tab_D_global_tr2[i].frequency*tree_i_weight;
            var_CP0_global_tr2_treeweight+=pow(tab_CP0_global[i]-mean_CP0_global_tr2_treeweight,2)*tab_CP0_global_tr2[i]*tree_i_weight;
            var_CP1_global_tr2_treeweight+=pow(tab_CP1_global[i]-mean_CP1_global_tr2_treeweight,2)*tab_CP1_global_tr2[i]*tree_i_weight;
            var_CP2_global_tr2_treeweight+=pow(tab_CP2_global[i]-mean_CP2_global_tr2_treeweight,2)*tab_CP2_global_tr2[i]*tree_i_weight;
            var_CP3_global_tr2_treeweight+=pow(tab_CP3_global[i]-mean_CP3_global_tr2_treeweight,2)*tab_CP3_global_tr2[i]*tree_i_weight;
            
            
        }
        
        var_D_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        var_D_p_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        var_r_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        var_r2_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        var_CP0_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        var_CP1_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        var_CP2_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        var_CP3_global_tr2_treeweight/=weight_tree_recombi_par_type[2];
        

        FILE * fichier_mean_var_tr2=NULL;
        char* tr2_char= "recombi_type2_";
        char* tr2_outfile;
        tr2_outfile = malloc(strlen(tr2_char)+1+strlen(outfile));
        strcpy(tr2_outfile, tr2_char);
        strcat(tr2_outfile,outfile);
        fichier_mean_var_tr2=fopen(tr2_outfile,"a");
        if (fichier_mean_var_tr2!=NULL)
        {
            fprintf(fichier_mean_var_tr2,"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",stats_D,1,N,nb_recombi_par_type[1],TCHANG,CONTRAC,mean_D_global_tr2_treeweight,var_D_global_tr2_treeweight,mean_D_p_global_tr2_treeweight,var_D_p_global_tr2_treeweight,mean_r_global_tr2_treeweight,var_r_global_tr2_treeweight,mean_r2_global_tr2_treeweight,var_r2_global_tr2_treeweight,mean_CP0_global_tr2_treeweight,var_CP0_global_tr2_treeweight,mean_CP1_global_tr2_treeweight,var_CP1_global_tr2_treeweight,mean_CP2_global_tr2_treeweight,var_CP2_global_tr2_treeweight,mean_CP3_global_tr2_treeweight,var_CP3_global_tr2_treeweight);
            
            fclose(fichier_mean_var_tr2);
        }
        else
            printf("Problem writing/openning per recombination type outfile\n");
        
        free(tr2_outfile);
        
        
        //TR3
        //mean
        double mean_D_global_tr3_treeweight=0,mean_D_p_global_tr3_treeweight=0,mean_r_global_tr3_treeweight=0,mean_r2_global_tr3_treeweight=0,mean_CP0_global_tr3_treeweight=0,mean_CP1_global_tr3_treeweight=0,mean_CP2_global_tr3_treeweight=0,mean_CP3_global_tr3_treeweight=0;
        for (i=0;i<nb_recombi_par_type[3]*(N*2-2)*(N*2-2);i++)
        {
            float tree_i_weight=0;
            
            if (stats_D==0||stats_D==1)
            {
                tree_i_weight=tab_D_global_tr3[i].tot_lbranch_A;
            }
            else if (stats_D==8)
            {
                tree_i_weight=tab_D_global_tr3[i].tot_lbranch_A*tab_D_global_tr3[i].tot_lbranch_B;
            }
            
            mean_D_global_tr3_treeweight+=tab_D_global_tr3[i].d*tab_D_global_tr3[i].frequency*tree_i_weight;
            mean_D_p_global_tr3_treeweight+=tab_D_global_tr3[i].d_p*tab_D_global_tr3[i].frequency*tree_i_weight;
            mean_r_global_tr3_treeweight+=tab_D_global_tr3[i].r*tab_D_global_tr3[i].frequency*tree_i_weight;
            mean_r2_global_tr3_treeweight+=tab_D_global_tr3[i].r2*tab_D_global_tr3[i].frequency*tree_i_weight;
            mean_CP0_global_tr3_treeweight+=tab_CP0_global_tr3[i]*tree_i_weight;
            mean_CP1_global_tr3_treeweight+=tab_CP1_global_tr3[i]*tree_i_weight;
            mean_CP2_global_tr3_treeweight+=tab_CP2_global_tr3[i]*tree_i_weight;
            mean_CP3_global_tr3_treeweight+=tab_CP3_global_tr3[i]*tree_i_weight;
            
        }
        
        mean_D_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        mean_D_p_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        mean_r_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        mean_r2_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        mean_CP0_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        mean_CP1_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        mean_CP2_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        mean_CP3_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        
        
        //variance
        double var_D_global_tr3_treeweight=0,var_D_p_global_tr3_treeweight=0,var_r_global_tr3_treeweight=0,var_r2_global_tr3_treeweight=0,var_CP0_global_tr3_treeweight=0,var_CP1_global_tr3_treeweight=0,var_CP2_global_tr3_treeweight=0,var_CP3_global_tr3_treeweight=0;
        for (i=0;i<nb_recombi_par_type[3]*(N*2-2)*(N*2-2);i++)
        {
            float tree_i_weight=0;
            
            if (stats_D==0||stats_D==1)
            {
                tree_i_weight=tab_D_global_tr3[i].tot_lbranch_A;
            }
            else if (stats_D==8)
            {
                tree_i_weight=tab_D_global_tr3[i].tot_lbranch_A*tab_D_global_tr3[i].tot_lbranch_B;
            }
            
            var_D_global_tr3_treeweight+=pow(tab_D_global_tr3[i].d-mean_D_global_tr3_treeweight,2)*tab_D_global_tr3[i].frequency*tree_i_weight;
            var_D_p_global_tr3_treeweight+=pow(tab_D_global_tr3[i].d_p-mean_D_p_global_tr3_treeweight,2)*tab_D_global_tr3[i].frequency*tree_i_weight;
            var_r_global_tr3_treeweight+=pow(tab_D_global_tr3[i].r-mean_r_global_tr3_treeweight,2)*tab_D_global_tr3[i].frequency*tree_i_weight;
            var_r2_global_tr3_treeweight+=pow(tab_D_global_tr3[i].r2-mean_r2_global_tr3_treeweight,2)*tab_D_global_tr3[i].frequency*tree_i_weight;
            var_CP0_global_tr3_treeweight+=pow(tab_CP0_global[i]-mean_CP0_global_tr3_treeweight,2)*tab_CP0_global_tr3[i]*tree_i_weight;
            var_CP1_global_tr3_treeweight+=pow(tab_CP1_global[i]-mean_CP1_global_tr3_treeweight,2)*tab_CP1_global_tr3[i]*tree_i_weight;
            var_CP2_global_tr3_treeweight+=pow(tab_CP2_global[i]-mean_CP2_global_tr3_treeweight,2)*tab_CP2_global_tr3[i]*tree_i_weight;
            var_CP3_global_tr3_treeweight+=pow(tab_CP3_global[i]-mean_CP3_global_tr3_treeweight,2)*tab_CP3_global_tr3[i]*tree_i_weight;
            
        }
        
        var_D_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        var_D_p_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        var_r_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        var_r2_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        var_CP0_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        var_CP1_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        var_CP2_global_tr3_treeweight/=weight_tree_recombi_par_type[3];
        var_CP3_global_tr3_treeweight/=weight_tree_recombi_par_type[3];

        FILE * fichier_mean_var_tr3=NULL;
        char* tr3_char= "recombi_type3_";
        char* tr3_outfile;
        tr3_outfile = malloc(strlen(tr3_char)+1+strlen(outfile));
        strcpy(tr3_outfile, tr3_char);
        strcat(tr3_outfile,outfile);
        fichier_mean_var_tr3=fopen(tr3_outfile,"a");
        if (fichier_mean_var_tr3!=NULL)
        {
            fprintf(fichier_mean_var_tr3,"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",stats_D,1,N,nb_recombi_par_type[1],TCHANG,CONTRAC,mean_D_global_tr3_treeweight,var_D_global_tr3_treeweight,mean_D_p_global_tr3_treeweight,var_D_p_global_tr3_treeweight,mean_r_global_tr3_treeweight,var_r_global_tr3_treeweight,mean_r2_global_tr3_treeweight,var_r2_global_tr3_treeweight,mean_CP0_global_tr3_treeweight,var_CP0_global_tr3_treeweight,mean_CP1_global_tr3_treeweight,var_CP1_global_tr3_treeweight,mean_CP2_global_tr3_treeweight,var_CP2_global_tr3_treeweight,mean_CP3_global_tr3_treeweight,var_CP3_global_tr3_treeweight);
            
            fclose(fichier_mean_var_tr3);
        }
        else
            printf("Problem writing/openning per recombination type outfile\n");
        
        free(tr3_outfile);
        
    }
    
    
    
    
        
        
    /// If you need full distribution
    /////Uncomment if needed
    /*
     float val_min,  val_bin=0.01;
     int cmpt_bin, nb_bin=100, ite_val, init_bin;
     //D
     //sort
     qsort(tab_D_global,n*((N*2-2)*(N*2-2)),sizeof(struct Linkage), &linkD_D_comparator);
     
     
     val_min=-0.5,cmpt_bin=0,ite_val=0;
     struct Distri * distri_D_global=(struct Distri *)malloc(nb_bin*sizeof(struct Distri));
     for (init_bin=0;init_bin<nb_bin;init_bin++)
     {
         distri_D_global[init_bin].frequency=0;
         distri_D_global[init_bin].bin_inf=val_min+init_bin*val_bin;
     }
     
     while (ite_val<n*(N*2-2)*(N*2-2))
     {
         if (tab_D_global[ite_val].d<val_min+val_bin)
         {
             float tree_i_weight=0;
             
             if (stats_D==0||stats_D==1)
             {
                 tree_i_weight=tab_D_global[ite_val].tot_lbranch_A;
             }
             else if (stats_D==8)
             {
                 tree_i_weight=tab_D_global[ite_val].tot_lbranch_A*tab_D_global[ite_val].tot_lbranch_B;
             }
             
             
             distri_D_global[cmpt_bin].frequency+=tab_D_global[ite_val].frequency*tree_i_weight/mean_treeweight;
             ite_val++;
         }
         else
         {
             cmpt_bin++;
             val_min+=val_bin;
         }
     }
     
     FILE * fichier_distrib_d=NULL;
     char* distrib_d_char= "Distribution_D_";
     char* distrib_d_outfile;
     distrib_d_outfile = malloc(strlen(distrib_d_char)+1+strlen(outfile));
     strcpy(distrib_d_outfile, distrib_d_char);
     strcat(distrib_d_outfile,outfile);
     fichier_distrib_d=fopen(distrib_d_outfile,"w");
     if (fichier_distrib_d!=NULL)
     {
         for (int i = 0; i < nb_bin; i++) {
             fprintf(fichier_distrib_d, "%.6f\t%.6f\n", distri_D_global[i].bin_inf, distri_D_global[i].frequency);
         }
         
         fclose(fichier_distrib_d);
     }
     else
         printf("Problem writing/openning outfile\n");
    */
    /*
    
    //Dp
    //sort
    qsort(tab_D_global,((N*2-2)*(N*2-2)),sizeof(struct Linkage), &linkD_D_p_comparator);
    val_min=0,val_bin=0.1,cmpt_bin=0;
    struct Distri * distri_D_p_global=(struct Distri*)malloc(nb_bin*sizeof(struct Distri));
    for (init_bin=0;init_bin<nb_bin;init_bin++)
    {
        distri_D_p_global[init_bin].frequency=0;
        distri_D_p_global[init_bin].bin_inf=val_min+init_bin*val_bin;
    }
    
    while (ite_val<n*(N*2-2)*(N*2-2))
    {
        if (tab_D_global[ite_val].d_p<val_min+val_bin)
        {
            distri_D_p_global[cmpt_bin].frequency+=tab_D_global[ite_val].frequency/(float)n;
            ite_val++;
        }
        else
        {
            cmpt_bin++;
            val_min+=val_bin;
        }
    }
    
    //r
    //sort
    qsort(tab_D_global,((N*2-2)*(N*2-2)),sizeof(struct Linkage), &linkD_r_comparator);
    val_min=0,val_bin=0.1,cmpt_bin=0;
    struct Distri * distri_r_global=(struct Distri*)malloc(nb_bin*sizeof(struct Distri));
    for (init_bin=0;init_bin<nb_bin;init_bin++)
    {
        distri_r_global[init_bin].frequency=0;
        distri_r_global[init_bin].bin_inf=val_min+init_bin*val_bin;
    }
    
    while (ite_val<n*(N*2-2)*(N*2-2))
    {
        if (tab_D_global[ite_val].r<val_min+val_bin)
        {
            distri_r_global[cmpt_bin].frequency+=tab_D_global[ite_val].frequency/(float)n;
            ite_val++;
        }
        else
        {
            cmpt_bin++;
            val_min+=val_bin;
        }
    }
    
    //r2
    //sort
    qsort(tab_D_global,((N*2-2)*(N*2-2)),sizeof(struct Linkage), &linkD_r2_comparator);
    val_min=0,val_bin=0.1,cmpt_bin=0;
    struct Distri * distri_r2_global=(struct Distri*)malloc(nb_bin*sizeof(struct Distri));
    for (init_bin=0;init_bin<nb_bin;init_bin++)
    {
        distri_r2_global[init_bin].frequency=0;
        distri_r2_global[init_bin].bin_inf=val_min+init_bin*val_bin;
    }
    
    while (ite_val<n*(N*2-2)*(N*2-2))
    {
        if (tab_D_global[ite_val].r2<val_min+val_bin)
        {
            distri_r2_global[cmpt_bin].frequency+=tab_D_global[ite_val].frequency/(float)n;
            ite_val++;
        }
        else
        {
            cmpt_bin++;
            val_min+=val_bin;
        }
    }
    
     */
     //free(distri_D_global);
     //free(distri_D_p_global);
     //free(distri_r_global);
     //free(distri_r2_global);
     
    
    
    free(tab_D_global);
    free(tab_D_global_tr0);
    free(tab_D_global_tr1);
    free(tab_D_global_tr2);
    free(tab_D_global_tr3);
    
    free(tab_CP0_global);
    free(tab_CP1_global);
    free(tab_CP2_global);
    free(tab_CP3_global);
     
        
    free(tab_CP0_global_tr0);
    free(tab_CP1_global_tr0);
    free(tab_CP2_global_tr0);
    free(tab_CP3_global_tr0);
    free(tab_CP0_global_tr1);
    free(tab_CP1_global_tr1);
    free(tab_CP2_global_tr1);
    free(tab_CP3_global_tr1);
    free(tab_CP0_global_tr2);
    free(tab_CP1_global_tr2);
    free(tab_CP2_global_tr2);
    free(tab_CP3_global_tr2);
    free(tab_CP0_global_tr3);
    free(tab_CP1_global_tr3);
    free(tab_CP2_global_tr3);
    free(tab_CP3_global_tr3);

    return 0;
}
