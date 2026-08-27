//
//  D_etude_lib.c
//  
//
//  Created by elise.kerdoncuff on 10/12/2023 (First version).
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <malloc/malloc.h> //Need to modify to include <malloc.h> for Linux

#include "D_etude_lib.h"

#define MAX(a,b) (a>=b?a:b)
#define MIN(a,b) (a<=b?a:b)


//function to compute the sum of a int table
int sumTable(int table[], int tailleTable)
{
    int i;
    int tot=0;
    for (i=0;i<tailleTable;i++)
    {
        tot+=table[i];
    }
    return tot;
}

//function which fills an array of size nb leaf with a "1" when the box corresponds to a leaf of the point entered in the function
int * leavesof(int * tableau, int point, struct InfosIndiv * story)
{
    if (story[point].desc1==-1)
    {
        tableau[point]=1;
        return tableau;
    }
    else{
        tableau=leavesof(tableau,story[point].desc1, story);
        tableau=leavesof(tableau,story[point].desc2, story);
        return tableau;
    }
}

//allows you to sort an array of floats
void  tricroiss( float * tab, int tab_size)
{
    int i=0;
    float tmp=0;
    int j=0;
    
    for(i = 0; i < tab_size; i++)          //to fill table[i]
    {
        for(j = i+1; j < tab_size; j++)    //We check if there is no lower number in the following boxes
        {
            if(tab[j] < tab[i])
            {
                tmp = tab[i];              //if it's the case we invert them
                tab[i] = tab[j];
                tab[j] = tmp;
            }
        }
    }
}

//// Functions for fast sorting tri_fusion
void fusion(float tableau[],int deb1,int fin1,int fin2)
{
    float *table1;
    int deb2=fin1+1;
    int compt1=deb1;
    int compt2=deb2;
    int i;
    
    table1=malloc((fin1-deb1+1)*sizeof(float));
    
    //we copy the elements from the beginning of the table
    for(i=deb1;i<=fin1;i++)
    {
        table1[i-deb1]=tableau[i];
    }
    
    for(i=deb1;i<=fin2;i++)
    {
        if (compt1==deb2) //all the elements of the first table have been used
        {
            break;
        }
        else if (compt2==(fin2+1)) //all the elements of the second table have been used
        {
            tableau[i]=table1[compt1-deb1]; //we add the remaining elements of the first table
            compt1++;
        }
        else if (table1[compt1-deb1]<tableau[compt2])
        {
            tableau[i]=table1[compt1-deb1]; //add en element in the first table
            compt1++;
        }
        else
        {
            tableau[i]=tableau[compt2]; //add en element in the second table
            compt2++;
        }
    }
    free(table1);
}


void tri_fusion_bis(float tableau[],int deb,int fin)
{
    if (deb!=fin)
    {
        int milieu=(fin+deb)/2;
        tri_fusion_bis(tableau,deb,milieu);
        tri_fusion_bis(tableau,milieu+1,fin);
        fusion(tableau,deb,milieu,fin);
    }
}

void tri_fusion(float tableau[],int longueur)
{
    if (longueur>0)
    {
        tri_fusion_bis(tableau,0,longueur-1);
    }
}

// function which sums all the branch lengths
float sumlbs_To(struct InfosIndiv InfosIndiv[], int N)
{
    int i;
    float ltot=0.0;
    for (i=0;i<N*2-1;i++)
        ltot+=InfosIndiv[i].lbs;
    return ltot;
}

// function which returns an array of size N*2-1 with 0s if the node is not a descendant of the point entered in the function and a 1 if it is
int* descde(int * tableau, int point, struct InfosIndiv InfosIndiv[]) //a leave has one descendant : itself
{
    if (InfosIndiv[point].desc1==-1)
    {
        tableau[point]=1;

    }
    else{
        tableau[InfosIndiv[point].desc1]=1;
        tableau[InfosIndiv[point].desc2]=1;
        tableau=descde(tableau,InfosIndiv[point].desc1, InfosIndiv);
        tableau=descde(tableau,InfosIndiv[point].desc2, InfosIndiv);

    }
    return tableau;
}

// function that returns the most recent common ancestor of two individuals (its number in the table)
int MRCA(int ind1, int ind2, struct InfosIndiv * story, int N)
{
    int MRCA;
    //Find root:
    int root=find_root(story,N);
    
    while (ind1!=root && ind1 !=ind2 && ind2!=root)
    {
        if (story[ind1].tps>=story[ind2].tps)
        {
            ind2=story[ind2].anc;
        }
         else {
             ind1=story[ind1].anc;
        }
    }
    if (ind1==root || ind2==root)
    {
        MRCA=root;
    }
    else if (ind1==ind2)
    {
        MRCA=ind1;
    }
    
    return MRCA;
}

    //function which allows you to write a tree in the newick format, which takes as input an individual (initially the root), the structure containing the information of the tree, a pointer to the character string and a pointer to the position where we write in the string
char * PrintTree(int indiv, struct InfosIndiv InfosIndiv[], int N, char *strtree, int * pos)
{
    if ((*pos)==0)
    {
        int i;
        for (i=N;i<N*2-1;i++)
        {
            if (InfosIndiv[i].anc==-1)
                indiv=i;
        }
    }
    if (malloc_size(strtree)-(*pos) < 150) //Need to modify to malloc_usable_size() for Linux
    {
        size_t sizestrtree = strlen(strtree);
        strtree=realloc(strtree, 2*sizestrtree );
        if (!strtree)
        {
            fprintf(stdout, "PB allocation");
            exit(254);
        }
    }
    
    if (InfosIndiv[indiv].desc1==-1)
    {
        sprintf(strtree+(*pos),"%d",indiv);
        float ajout=ceil(log10(indiv+1));
        int ajout_pos= (int) ajout;
        if (ajout_pos==0) ajout_pos=1;
        (*pos)+=ajout_pos;
        sprintf(strtree+(*pos),":%7.6f",InfosIndiv[indiv].lbs);
        (*pos)+=9;
        
        
    }
    else
    {
        sprintf(strtree+(*pos)++,"(");
        strtree=PrintTree(InfosIndiv[indiv].desc1,InfosIndiv,N,strtree,pos);
        sprintf(strtree+(*pos)++,",");
        strtree=PrintTree(InfosIndiv[indiv].desc2,InfosIndiv,N,strtree,pos);
        sprintf(strtree+(*pos)++,")");
        sprintf(strtree+(*pos),":%7.6f",InfosIndiv[indiv].lbs);
        (*pos)+=9;
    
    }
    return strtree;
}

//function which generates an ARG from a struct InfosIndiv table of the size of the nodes, the time at which the pop size changes and the force of contraction
struct InfosIndiv * Init_ARG(struct InfosIndiv * story,int N,double tchang,float contrac)
{
    int i;
    int comp=N; //N coalescent event
    int * indiv=malloc((N*2-1)*sizeof(int)); //individual pool, leaves of the tree
    
    float * tps = calloc(N,sizeof(float));// coalescent times
    float * tint=calloc(N-1,sizeof(float)); // interval between times
    
    
    for (i=0;i<(N*2-1);i++) //to keep track of the indivs
    {
        indiv[i] = -1;
    }
    
    struct InfosIndiv newInd = {-1,-1,-1,0.0,0.0,1,0,0}; //new ind to fill later
    
    for (i=0;i<N;i++) //init leaves
    {
        indiv[i] = i;
        story[i]=newInd;
    }
    int node=N;
    
    while (comp !=1)
    {
        // compute coalescent time
        float nombre;
        do
        {
           nombre  = (float)rand() / (float)RAND_MAX;
        } while (nombre == 0 || nombre == 1);
        
        
        float tpscoal= ((-log(nombre))/((comp*(comp-1))/2)); // interval between coalescent events, exponential size 1
        if (tpscoal + tps[N-comp] < tchang) // depending of time of contraction
        {
            tps[N-comp+1]= tpscoal + tps[N-comp];
            tint[N-comp]=tpscoal;
        }
        else if ( tpscoal + tps[N-comp]>tchang && tps[N-comp]<tchang) //when the change of size happen between two nodes
        {
            tps[N-comp+1]= (tpscoal + tps[N-comp]- tchang)*contrac + tchang;
            tint[N-comp]= ((tpscoal + tps[N-comp]- tchang)*contrac + tchang - tps[N-comp]);
        }
        else if (tps[N-comp]>=tchang)
        {
            tps[N-comp+1]=tps[N-comp]+tpscoal*contrac;
            tint[N-comp]=tpscoal*contrac;
        }
        //choice of coalescing individuals
        int coalA;
        int coalB;
        do
        {
            coalA=rand()%(N*2-1);
        } while (indiv[coalA]==-1);
        
        do
        {
            coalB=rand()%(N*2-1);
        } while (coalA==coalB || indiv[coalB]==-1); // need two different individuals
        
        //update new node
        story[coalA].anc=story[coalB].anc=node;
        
        story[node].desc1=coalA;
        story[node].desc2=coalB;
        story[node].anc=-1;
        story[node].tps=tps[N-comp+1];
        
        int * table = NULL;
        table=(int *)calloc(N,sizeof(int));
        table=leavesof(table,node,story);
        story[node].nb_leaves=sumTable(table,N);
        free(table);
        
        //add node to list of possible coalescing nodes
        indiv[node]=node;
        //delete the ones that just coalesce
        indiv[coalA]=-1;
        indiv[coalB]=-1;
        
        
        node++;
        comp--;
        
    }
    
    //compute branch lenght of the nodes
    for (i=0;i<N*2-1;i++)
    {
        if (story[i].anc!=-1)
            story[i].lbs=story[story[i].anc].tps - story[i].tps;
        else
            story[i].lbs=0;
    }
    
    //compute lb1 et lb2 sum of branch length below nodes
    for (i=0; i< N*2-1;i++)
    {
        if (story[i].desc1==-1)
        {
            story[i].lb1=0;
            story[i].lb2=0;
        }
        else
        {
            int * desctest1 = NULL;
            desctest1=(int *)calloc(N*2-1,sizeof(int));
            int * desctest2 = NULL;
            desctest2=(int *)calloc(N*2-1,sizeof(int));
            desctest1=descde(desctest1,story[i].desc1,story); //all descendants marqued by 1
            desctest2=descde(desctest2,story[i].desc2,story);
            
            desctest1[story[i].desc1]=1; //add descendants
            desctest2[story[i].desc2]=1;
            float comp_lb1=0;
            float comp_lb2=0;
            int j;
            for (j=0;j<N*2-1;j++)
            {
                if (desctest1[j]==1)
                {
                    comp_lb1+=story[j].lbs;
                }
                if (desctest2[j]==1)
                {
                    comp_lb2+=story[j].lbs;
                }
            }
            story[i].lb1=comp_lb1;
            story[i].lb2=comp_lb2;
            
            free(desctest1);
            free(desctest2);
        }
    }
    
    free(indiv);
    free(tps);
    free(tint);
    
    
    return story ;
}


//function to perform a SMC algorithm from a tree structure
struct InfosIndiv * SMC_tree(struct InfosIndiv * story,int N,double tchang,float contrac)
{
    int i=0;
    int root=0;// root of the tree
    float ltot=0; //total branch length sum
    float trec=0.0; // time of recombination
    int whererec=0; // branch of recombination
    int j=0;
    
    
    float * tps = (float *)calloc(N,sizeof(float));// time of coalescent event
    float * tint= (float *)calloc(N-1,sizeof(float)); // interval between coalescent events
    
    
    for (j=N;j<N*2-1;j++)
    {
        tps[j-N]=story[j].tps;
    }
    
    for (j=0;j<(N*2-1);j++)
    {
        ltot += story[j].lbs;
        if (story[j].anc==-1)
            root=j;
    }

    tps[N-1]=0;
    tricroiss(tps, N);//sort depending of the time of coalescence
    for (i=0;i<N-1;i++)
    {
        tint[i]=(tps[i+1]-tps[i]);
    }
    ////
    //Recombination event position
    ////
    float lunhook; //unhook position - chosen uniformly on total branch length
    do
    {
        lunhook = ((float)rand() / (float)RAND_MAX)*ltot;
    } while (lunhook == 0 || lunhook == ltot);
    
    
    int x=0; // compt
    
    while (x < N*2-1)
    {
        if (lunhook <= story[x].lbs)
        {
            whererec=x; //the branch
            x=N*2; //break
        }
        else
        {
            lunhook -= story[x].lbs;
            x++;
        }
        
    }
    
    float tunhook;
    tunhook= lunhook + story[whererec].tps; //unhook time
    float newl = lunhook + story[whererec].tps;
    
    float newrec=0;
    float t=0;
    int k=0; //compt for interbal table
    float f=1.0; //factor of time change
    int flag=0; // to know if the change in size already happen
    if (tunhook>=tchang) //if time of change lower than recombination time
    {
        flag=10;
        f=contrac;
    }
    
    while (k<(N-1)) // N-1, nb interval in tint
    {
        t+=tint[k];
        if (newl>=tchang && flag==0)
        {
            f=contrac;
            newl=(newl-tchang)*f+tchang;
            flag=10;
        }
        
        if (newl<=t) //to know in which tint is newl //(newl<=t)
        {
            int i=N-k; // compute i (nb individual present at the time) to know of to draw kingamn exp
            do
            {
                newrec= (-log((float)rand() / (float)RAND_MAX)/((i*(i-1))/2))*f+newl;
            } while (newrec == 0 || newrec == t);
            
            if (newrec>=tchang && flag==0)
            {
                f=contrac;
                newrec=(newrec - tchang)*f+tchang;
                flag=10;
            }
            
            if (newrec<t)
            {
                trec=newrec;
                k=N; //break while
            }
            else if ((k==(N-2) && newrec>=t && flag!=0)||(k==(N-2) && newrec>=t && newrec<tchang))
            {

                i=N-k;
                trec= (-log((float)rand() / (float)RAND_MAX)/((i*(i-1))/2))*f+t; // i change, re-draw the exp
                
                k=N+1; //break while
            }
            
            else if (k==(N-2) && newrec>=t && newrec>=tchang)
            {
    
                i=N-k;
                f=1;
                f=contrac;
                trec=(-log((float)rand() / (float)RAND_MAX)/((i*(i-1))/2))*f+tchang; // i change, t change, re-draw the exp
                k=N+2; //break while

            }
            else
            {
                newl=t;
                
                k++;
            }
            
        }
        else
        {
           
            k++;
        }
    }
    
    if (trec==0)
    {
        if (newrec!=0)
        {
            newrec=(newrec-tchang)*f + tchang;
            trec=newrec;
        }
        else
        {
            i=2;
            trec=  (-log((float)rand() / (float)RAND_MAX)/((i*(i-1))/2))*f+newl;
        }

    }
        
    
   
    /////////
    // Erecombination event- hook bakc
    /////////
    //branche present at trec
    int * tab_anc= (int *)calloc(N*2-1,sizeof(int));
    for (i=0;i<N*2-1;i++)
    {
        if (story[i].tps>trec)
            tab_anc[i]=1;
    }
    
    int * tab_desc= (int *)calloc(N*2-1,sizeof(int));;
    for (i=0;i<N*2-1;i++)
    {
        if (tab_anc[i]==1)
        {
            tab_desc[story[i].desc1]=1;
            tab_desc[story[i].desc2]=1;
        }
    }
    if (sumTable(tab_anc, N*2-1)==0)//if no node older than trec
    {
        tab_desc[root]=1;
    }
    
    
    
    for (i=0;i<N*2-1;i++)
    {
        if (tab_anc[i]==1)
            tab_desc[i]=0; //delete ancestors and descendand to update
    }

    int nbdesc=sumTable(tab_desc,N*2-1); //nb of branches present at trec to draw the exp
    int choice;
    int wherehook = -1;
    choice=rand()%nbdesc;
    j=0;
    while (j<N*2-1 && choice!=-1)
    {
        if (tab_desc[j]==1)
            choice--; //stupid but worked
        if (choice==-1)
        {
            wherehook=j;
        }
        j++;

    }


    ////////////
    //Update of "story"
    ////////////
    if (wherehook==whererec) //no change
    {
         story[whererec].anc=story[wherehook].anc;
    }
    else if (story[whererec].anc==story[wherehook].anc || story[whererec].anc==wherehook)
    {
        story[story[whererec].anc].tps=trec;
    }
    else
    {
        int notmovedwhererec=-1; //other descendant of the ancestor of whererec
        int notmovedwherehook=-1; //aother descendant of the ancestor of  wherehook
        if (story[whererec].anc!=-1 )
        {
            if (story[story[whererec].anc].desc1==whererec)
            {
                notmovedwhererec=story[story[whererec].anc].desc2;
            }
            else
            {
                notmovedwhererec=story[story[whererec].anc].desc1;
            }
        }
        
        if (story[wherehook].anc!=-1)
        {
            if (story[story[wherehook].anc].desc1==wherehook)
            {
                notmovedwherehook=story[story[wherehook].anc].desc2;
                story[story[wherehook].anc].desc1=notmovedwherehook;
                story[story[wherehook].anc].desc2=story[whererec].anc;
            }
            else
            {
                notmovedwherehook=story[story[wherehook].anc].desc1;
                story[story[wherehook].anc].desc1=notmovedwherehook;
                story[story[wherehook].anc].desc2=story[whererec].anc;
            }
        }

        if (story[story[whererec].anc].anc!=-1)
        {
            if (story[story[story[whererec].anc].anc].desc1==story[whererec].anc)
            {
                story[story[story[whererec].anc].anc].desc1=notmovedwhererec;
            }
            else
            {
                story[story[story[whererec].anc].anc].desc2=notmovedwhererec;
            }
        }
        
        story[notmovedwhererec].anc=story[story[whererec].anc].anc;
        story[story[whererec].anc].anc=story[wherehook].anc;
        story[wherehook].anc=story[whererec].anc;
        story[story[whererec].anc].desc1=whererec;
        story[story[whererec].anc].desc2=wherehook;
        story[story[whererec].anc].tps=trec;



        for (i=N;i<(N*2-1);i++)
        {
            int * table = NULL;
            table=(int *)calloc(N,sizeof(int));
            table=leavesof(table,i,story);
            story[i].nb_leaves=sumTable(table,N);
            free(table);
        }

    }
    // compute branch length
    for (i=0;i<(N*2-1);i++)
    {
        if (story[i].anc!=-1){
            story[i].lbs= story[story[i].anc].tps - story[i].tps;
        }
        else
        {
            story[i].lbs=0;
        }
    }

    
    //compute sum lb1 lb2 below nodes
    for (i=0; i< N*2-1;i++)
    {
        if (story[i].desc1==-1)
        {
            story[i].lb1=0;
            story[i].lb2=0;
        }
        else
        {

            int * desctest1 = NULL;
            desctest1=(int *)calloc(N*2-1,sizeof(int));
            int * desctest2 = NULL;
            desctest2=(int *)calloc(N*2-1,sizeof(int));
            desctest1=descde(desctest1,story[i].desc1,story); //all descendants marqued by 1
            desctest2=descde(desctest2,story[i].desc2,story);
            
            desctest1[story[i].desc1]=1; //add this descendants
            desctest2[story[i].desc2]=1;
            float comp_lb1=0;
            float comp_lb2=0;
            int j;
            for (j=0;j<N*2-1;j++)
            {
                if (desctest1[j]==1)
                {
                    comp_lb1+=story[j].lbs;
                }
                if (desctest2[j]==1)
                {
                    comp_lb2+=story[j].lbs;
                }
            }
            story[i].lb1=comp_lb1;
            story[i].lb2=comp_lb2;
            
            free(desctest1);
            free(desctest2);
        }
    }
    return story;
}




//Compares two positions in a struct table Two_Floats (by nb2) to be able to subsequently sort them with qsort
int Two_Floats_nb2_comparator(const void *v1, const void *v2)
{
    const struct Two_Floats *p1 = (struct Two_Floats *)v1;
    const struct Two_Floats *p2 = (struct Two_Floats *)v2;
    if (p1->nb2 < p2->nb2)
        return -1;
    else if (p1->nb2 > p2->nb2)
        return +1;
    else
        return 0;
}

//function compares table of int one to one, return 1 if different otherwise return 0
int cmp_tab_egal(int *tab1,int *tab2,int taille)
{
    int i;
    for(i=0;i<taille;i++)
    {
        if (tab1[i]!=tab2[i])
        {
            return 1;
        }
    }
    return 0;
}

//function that compares two tables of tables, return 1 if different, 0 if similar
int cmp_tabdetabe(int **tab1, int **tab2, int taille, int N)
{
    int i,j, tabeg, result=0; //0: compatible (1: incompatible) 10: compatible non rooted topology 11: incomptible
    for (i=0;i<taille;i++)
    {
        j=0;
        while (j<taille)
        {
            tabeg=cmp_tab_egal(tab1[i], tab2[j],N); // table identical or not
            j++;
            if (tabeg==0) //identical tables
            {
                j=taille+2; // next table
            }
        }
        if (j==taille) //no table tab2 identical to tab1 (compare to every tab2)
        {
            if (sumTable(tab1[i],N)!=N-1) //if ancestor of every oen exept one, similar as a singleton -> same non rooted topology
            {
                result=11;
                return result;// leave function tables are different
            }
            else
            {
                result=10;
            }
            
        }
        
    }
    return result; // we always found similar bi-partition, all partitions are egal
}

//Function to know the root of the tree, from a tree matrix and the number of leaves
int find_root(struct InfosIndiv * story, int N)
{
    int i, root;
    for (i=0;i<2*N-1;i++)
    {
        if (story[i].anc==-1)
            root=i;
    }
    return root;
}

///
//////
////////
// Linkage Desiquilibrium
///////
/////
///

//Function that returns the frequency of a branch of a tree (frequency if a mutation falls on it)
float frequency_mutations(struct InfosIndiv * story,int indiv, int N)
{
    int * table=(int *)calloc(N,sizeof(int));
    table=leavesof(table,indiv,story);
    return  (float)sumTable(table,N)/(float)N;
}

//Function that returns the joint frequency of two branches of a tree (frequency if mutations fall on it)
float frequency_mutationAB(struct InfosIndiv * story,int indiv1, int indiv2, int N)
{
    int * table_1=(int *)calloc(N,sizeof(int));
    table_1=leavesof(table_1,indiv1,story);
    
    int * table_2=(int *)calloc(N,sizeof(int));
    table_2=leavesof(table_2,indiv2,story);
    
    int comp_12=0;
    int i;
    for (i=0;i<N;i++)
        if (table_1[i]==1 && table_2[i]==1)
            comp_12++;
    
    return (float)comp_12/(float)N;
    
}

//Function that returns the joined frequency of two branches TWO DIFFERENCE TREES (frequency if mutations fall on them)
float frequency_mutationAB_twotrees(struct InfosIndiv * story1,int indiv1,struct InfosIndiv * story2, int indiv2, int N)
{
    int * table_1=(int *)calloc(N,sizeof(int));
    table_1=leavesof(table_1,indiv1,story1);
    
    int * table_2=(int *)calloc(N,sizeof(int));
    table_2=leavesof(table_2,indiv2,story2);
    
    int comp_12=0;
    int i;
    for (i=0;i<N;i++)
        if (table_1[i]==1 && table_2[i]==1)
            comp_12++;
    
    return (float)comp_12/(float)N;
    
}
////
// Similar with Minor allele frequency MAF and not derived allele frequency
////
//Function that returns the joint frequency of two branches of a tree (frequency if mutations fall on it)MAF
float frequency_mutationAB_MAF(struct InfosIndiv * story,int indiv1, int indiv2, int N)
{
    int * real_table_1=(int *)calloc(N,sizeof(int));
    int * table_1=(int *)calloc(N,sizeof(int));
    real_table_1=leavesof(real_table_1,indiv1,story);
    if ((float)sumTable(real_table_1, N) <= (float)N/2) //if already MAF
    {
        table_1=real_table_1;
    }
    else
    {

        int i;
        for (i=0;i<N;i++)
        {
            if (real_table_1[i]==1) //reverse 1/0 in table
                table_1[i]=0;
            else
                table_1[i]=1;
        }
    }

    
    int * real_table_2=(int *)calloc(N,sizeof(int));
    int * table_2 = (int *)calloc(N,sizeof(int));
    real_table_2=leavesof(real_table_2,indiv2,story);
    
    if ((float)sumTable(real_table_2, N) <= (float)N/2) //if already MAF
    {

        table_2=real_table_2;

    }
    else
    {

        int i;
        for (i=0;i<N;i++)
        {
            if (real_table_2[i]==1) //reverse 1/0 in table
                table_2[i]=0;
            else
                table_2[i]=1;
        }
    }
    
    
    int comp_12=0;
    int i;
    for (i=0;i<N;i++)
        if (table_1[i]==1 && table_2[i]==1)
            comp_12++;
    
    free(table_1);
    free(table_2);
    return (float)comp_12/(float)N;
    
}

//Function that returns the joined frequency of two branches TWO DIFFERENCE TREES (frequency if mutations fall on them) MAF
float frequency_mutationAB_twotrees_MAF(struct InfosIndiv * story1,int indiv1,struct InfosIndiv * story2, int indiv2, int N)
{
    
    int * real_table_1=(int *)calloc(N,sizeof(int));
    int * table_1=(int *)calloc(N,sizeof(int));
    real_table_1=leavesof(real_table_1,indiv1,story1);
    if ((float)sumTable(real_table_1, N) <= (float)N/2) //if already MAF
        table_1=real_table_1;
    else
    {
        int i;
        for (i=0;i<N;i++)
        {
            if (real_table_1[i]==1) //reverse 1/0 in table
                table_1[i]=0;
            else
                table_1[i]=1;
        }
    }
   
    
    int * real_table_2=(int *)calloc(N,sizeof(int));
    int * table_2=(int *)calloc(N,sizeof(int));
    real_table_2=leavesof(real_table_2,indiv2,story2);
    
    if ((float)sumTable(real_table_2, N) <= (float)N/2) //if already MAF
        table_2=real_table_2;
    else
    {
        int i;
        for (i=0;i<N;i++)
        {
            if (real_table_2[i]==1) //reverse 1/0 in table
                table_2[i]=0;
            else
                table_2[i]=1;

        }
    }
    
    
    int comp_12=0;
    int i;
    for (i=0;i<N;i++)
        if (table_1[i]==1 && table_2[i]==1)
            comp_12++;
    
    
    free(table_1);
    free(table_2);
    return (float)comp_12/(float)N;
    
}



//Function to calculate the values ​​of D, Dp and r2 on a single tree, returns a Linkage struct table containing the values ​​of D and its proportion
//Input, the tree matrix and the number of leaves
struct Linkage * compute_D_ontree(struct InfosIndiv * story, int N)
{
    //init
    int nb_indiv=N*2-1;
    int ind_A,ind_B;
    float Ltot=sumlbs_To(story,N);
    struct Linkage * tab_D= (struct Linkage *)malloc(((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    int compteur_D=0;
    
    for (ind_A=0;ind_A<nb_indiv-1;ind_A++)
    {
        //Ind A
        float freq_A = frequency_mutations(story,ind_A,N); //Frequency
        float long_A = story[ind_A].lbs/Ltot; //branch length, divided by Tot Length
        
        for (ind_B=0;ind_B<nb_indiv-1;ind_B++)
        {
            //Frequency ind B
            float freq_B = frequency_mutations(story,ind_B,N); //Frequency
            float long_B = story[ind_B].lbs/Ltot;//branch length, divided by Tot Length
            
            //Frequency AB
            float freq_AB=frequency_mutationAB(story,ind_A,ind_B,N);
            
            //Final values
            //D
            float val_D = freq_AB - freq_A*freq_B;
            
            //D'
            float val_Dp;
            if (val_D>0)
                val_Dp=val_D/MIN(freq_A*(1-freq_B),(1-freq_A)*freq_B);
            else
                val_Dp=val_D/MAX(-freq_A*freq_B,-(1-freq_A)*(1-freq_B));
            
            //r
            float val_r=(val_D)/sqrt(freq_A*(1-freq_A)*freq_B*(1-freq_B));
            //r2
            float val_r2=(val_D*val_D)/(freq_A*(1-freq_A)*freq_B*(1-freq_B));
            
            float freq_D = long_A*long_B; //frequency of the D value (based of branch length)

            
            
            tab_D[compteur_D].d=val_D;
            tab_D[compteur_D].d_p=val_Dp;
            tab_D[compteur_D].r=val_r;
            tab_D[compteur_D].r2=val_r2;
            tab_D[compteur_D].pqA_pqB=(freq_A*(1-freq_A)*freq_B*(1-freq_B));
            tab_D[compteur_D].frequency=freq_D;
            tab_D[compteur_D].lbranch_A=story[ind_A].lbs;
            tab_D[compteur_D].lbranch_B=story[ind_B].lbs;
            tab_D[compteur_D].tot_lbranch_A=Ltot;
            tab_D[compteur_D].tot_lbranch_B=Ltot;
            compteur_D++;
            
        }
    }
    
    return tab_D;
}


//Function to calculate the values ​​of D, Dp and r2 on TWO DIFFERENT TREES, returns a Linkage struct table including the values ​​of D and its proportion
//Input, the tree matrices and the number of leaves
struct Linkage * compute_D_twotrees(struct InfosIndiv * story1,struct InfosIndiv * story2 , int N)
{
    //init
    int nb_indiv=N*2-1;
    int ind_A=0,ind_B=0;
    float Ltot1=sumlbs_To(story1,N);
    float Ltot2=sumlbs_To(story2,N);
    struct Linkage * tab_D= (struct Linkage *)malloc(((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    int compteur_D=0;
    int root1=find_root(story1,N);
    int root2=find_root(story2,N);
    
    while (ind_A<nb_indiv)
    {
        ind_B=0;
        if( ind_A!=root1)
        {
        
            //Ind A
            float freq_A = frequency_mutations(story1,ind_A,N); //Frequency
            float long_A = story1[ind_A].lbs/Ltot1; //branch length, divided by Tot Length
            
            while (ind_B<nb_indiv)
            {
                if (ind_B!=root2)
                {
                    //Frequency ind B
                    float freq_B = frequency_mutations(story2,ind_B,N); //Frequency
                    float long_B = story2[ind_B].lbs/Ltot2;//branch length, divided by Tot Length
                    
                    //Frequency AB
                    float freq_AB=frequency_mutationAB_twotrees(story1,ind_A,story2,ind_B,N);
                    
                    //Final values
                    //D
                    float val_D = freq_AB - freq_A*freq_B;
                    
                    //D'
                    float val_Dp;
                    if (val_D>0)
                        val_Dp=val_D/MIN(freq_A*(1-freq_B),(1-freq_A)*freq_B);
                    else
                        val_Dp=val_D/MAX(-freq_A*freq_B,-(1-freq_A)*(1-freq_B));
                    
                    //r
                    float val_r=(val_D)/sqrt(freq_A*(1-freq_A)*freq_B*(1-freq_B));
                    //r2
                    float val_r2=(val_D*val_D)/(freq_A*(1-freq_A)*freq_B*(1-freq_B));
                    
                    float freq_D = long_A*long_B; //frequency of the D value (based of branch length)
 
                    tab_D[compteur_D].d=val_D;
                    tab_D[compteur_D].d_p=val_Dp;
                    tab_D[compteur_D].r=val_r;
                    tab_D[compteur_D].r2=val_r2;
                    tab_D[compteur_D].pqA_pqB=(freq_A*(1-freq_A)*freq_B*(1-freq_B));
                    tab_D[compteur_D].frequency=freq_D;
                    tab_D[compteur_D].lbranch_A=story1[ind_A].lbs;
                    tab_D[compteur_D].lbranch_B=story2[ind_B].lbs;
                    tab_D[compteur_D].tot_lbranch_A=Ltot1;
                    tab_D[compteur_D].tot_lbranch_B=Ltot2;
                    compteur_D++;
                }
                
                ind_B++;
            }
        }
        ind_A++;
    }
    
    return tab_D;
}

////
// Similar but MAF
////
//Function to calculate the values ​​of D, Dp and r2 on a single tree, returns a Linkage struct table containing the values ​​of D and its proportion
//Input, the tree matrix and the number of leaves
struct Linkage * compute_D_ontree_MAF(struct InfosIndiv * story, int N)
{
    //init
    int nb_indiv=N*2-1;
    int ind_A,ind_B;
    float Ltot=sumlbs_To(story,N);
    struct Linkage * tab_D= (struct Linkage *)malloc(((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    int compteur_D=0;
    
    for (ind_A=0;ind_A<nb_indiv-1;ind_A++)
    {
        //Ind A
        float real_freq_A = frequency_mutations(story,ind_A,N); //Frequency
        float freq_A=MIN(real_freq_A,1-real_freq_A); //MAF
        float long_A = story[ind_A].lbs/Ltot; //branch length, divided by Tot Length
        
        for (ind_B=0;ind_B<nb_indiv-1;ind_B++)
        {
            //Frequency ind B
            float real_freq_B = frequency_mutations(story,ind_B,N); //Frequency
            float freq_B= MIN(real_freq_B,1-real_freq_B); //MAF
            float long_B = story[ind_B].lbs/Ltot;//branch length, divided by Tot Length
            
            //Frequency AB
            float freq_AB=frequency_mutationAB_MAF(story,ind_A,ind_B,N);
            
            //Final values
            //D
            float val_D = freq_AB - freq_A*freq_B;
            
            //D'
            float val_Dp;
            if (val_D>0)
                val_Dp=val_D/MIN(freq_A*(1-freq_B),(1-freq_A)*freq_B);
            else
                val_Dp=val_D/MAX(-freq_A*freq_B,-(1-freq_A)*(1-freq_B));
            
            //r
            float val_r=(val_D)/sqrt(freq_A*(1-freq_A)*freq_B*(1-freq_B));
            //r2
            float val_r2=(val_D*val_D)/(freq_A*(1-freq_A)*freq_B*(1-freq_B));
            
            float freq_D = long_A*long_B;  //frequency of the D value (based of branch length)
  
            tab_D[compteur_D].d=val_D;
            tab_D[compteur_D].d_p=val_Dp;
            tab_D[compteur_D].r=val_r;
            tab_D[compteur_D].r2=val_r2;
            tab_D[compteur_D].pqA_pqB=(freq_A*(1-freq_A)*freq_B*(1-freq_B));
            tab_D[compteur_D].frequency=freq_D;
            tab_D[compteur_D].lbranch_A=story[ind_A].lbs;
            tab_D[compteur_D].lbranch_B=story[ind_B].lbs;
            tab_D[compteur_D].tot_lbranch_A=Ltot;
            tab_D[compteur_D].tot_lbranch_B=Ltot;
            compteur_D++;
            
        }
    }
    
    return tab_D;
}


//Function to calculate the values ​​of D, Dp and r2 on TWO DIFFERENT TREES, returns a Linkage struct table including the values ​​of D and its proportion
//Input, the tree matrices and the number of leaves
struct Linkage * compute_D_twotrees_MAF(struct InfosIndiv * story1,struct InfosIndiv * story2 , int N)
{
    //init
    int nb_indiv=N*2-1;
    int ind_A=0,ind_B=0;
    float Ltot1=sumlbs_To(story1,N);
    float Ltot2=sumlbs_To(story2,N);
    struct Linkage * tab_D= (struct Linkage *)malloc(((N*2-2)*(N*2-2))*sizeof(struct Linkage));
    int compteur_D=0;
    int root1=find_root(story1,N);
    int root2=find_root(story2,N);

    
    while (ind_A<nb_indiv)
    {
        ind_B=0;
        if( ind_A!=root1)
        {
            
            //Ind A
            float real_freq_A = frequency_mutations(story1,ind_A,N); //Frequency
            float freq_A= MIN(real_freq_A,1-real_freq_A); //MAF
            float long_A = story1[ind_A].lbs/Ltot1; //branch length, divided by Tot Length
            
            while (ind_B<nb_indiv)
            {
                if (ind_B!=root2)
                {
                    //Frequency ind B
                    float real_freq_B = frequency_mutations(story2,ind_B,N); //Frequency
                    float freq_B= MIN(real_freq_B,1-real_freq_B); //MAF
                    float long_B = story2[ind_B].lbs/Ltot2;//branch length, divided by Tot Length
                    
                    //Frequency AB
                    float freq_AB=frequency_mutationAB_twotrees_MAF(story1,ind_A,story2,ind_B,N);
                    

                    //Final values
                    //D
                    float val_D = freq_AB - freq_A*freq_B;
                    
                    //D'
                    float val_Dp;
                    if (val_D>0)
                        val_Dp=val_D/MIN(freq_A*(1-freq_B),(1-freq_A)*freq_B);
                    else
                        val_Dp=val_D/MAX(-freq_A*freq_B,-(1-freq_A)*(1-freq_B));
                    
                    //r
                    float val_r=(val_D)/sqrt(freq_A*(1-freq_A)*freq_B*(1-freq_B));
                    
                    //r2
                    float val_r2=(val_D*val_D)/(freq_A*(1-freq_A)*freq_B*(1-freq_B));
                    
                    float freq_D = long_A*long_B; //frequency of the D value (based of branch length)
                    tab_D[compteur_D].d=val_D;
                    tab_D[compteur_D].d_p=val_Dp;
                    tab_D[compteur_D].r=val_r;
                    tab_D[compteur_D].r2=val_r2;
                    tab_D[compteur_D].pqA_pqB=(freq_A*(1-freq_A)*freq_B*(1-freq_B));
                    tab_D[compteur_D].frequency=freq_D;
                    tab_D[compteur_D].lbranch_A=story1[ind_A].lbs;
                    tab_D[compteur_D].lbranch_B=story2[ind_B].lbs;
                    tab_D[compteur_D].tot_lbranch_A=Ltot1;
                    tab_D[compteur_D].tot_lbranch_B=Ltot2;
                    compteur_D++;
                }
                
                ind_B++;
            }
        }
        ind_A++;
    }
    
    return tab_D;
}

//Compute theoretical computaiton of Var D0 for one tree
float calcul_varD0(struct InfosIndiv * story, int N)
{
    int i,j;
    int root=find_root(story,N);
    float Ltot=sumlbs_To(story, N); //for normalisation
    float varD0=0;
    //For each leaf, compute H_bar
    float * H_bar = calloc(N,sizeof(float));
    for (i=0;i<N;i++)
    {
        for (j=0;j<N;j++)
        {
            H_bar[i]+=(story[root].tps-story[MRCA(i,j,story,N)].tps)/Ltot;
        }
        H_bar[i]=H_bar[i]/(float)N;
    }
    //Compute variance
    for (i=0;i<N;i++)
    {
        for (j=0;j<N;j++)
        {
            float Hij=(story[root].tps-story[MRCA(i,j,story,N)].tps)/Ltot;
            varD0+=(Hij - H_bar[i])*(Hij - H_bar[j]);
        }
        
    }
    varD0=varD0/(float)(N*N);
    
    return varD0;
}


////
// Comparison, sorting
////

// Compare two values ​​of .D in a strcut linkage table to then be able to sort them with qsort
int linkD_D_comparator(const void *v1, const void *v2)
{
    const struct Linkage *p1 = (struct Linkage *)v1;
    const struct Linkage *p2 = (struct Linkage *)v2;
    if (p1->d < p2->d)
        return -1;
    else if (p1->d > p2->d)
        return +1;
    else
        return 0;
}

// Compare two values ​​of .D_p in a strcut linkage table to then be able to sort them with qsort
int linkD_D_p_comparator(const void *v1, const void *v2)
{
    const struct Linkage *p1 = (struct Linkage *)v1;
    const struct Linkage *p2 = (struct Linkage *)v2;
    if (p1->d_p < p2->d_p)
        return -1;
    else if (p1->d_p > p2->d_p)
        return +1;
    else
        return 0;
}

// Compare two values ​​of .r2 in a strcut linkage table to then be able to sort them with qsort
int linkD_r2_comparator(const void *v1, const void *v2)
{
    const struct Linkage *p1 = (struct Linkage *)v1;
    const struct Linkage *p2 = (struct Linkage *)v2;
    if (p1->r2 < p2->r2)
        return -1;
    else if (p1->r2 > p2->r2)
        return +1;
    else
        return 0;
}

// Compare two values ​​of .r in a strcut linkage table to then be able to sort them with qsort
int linkD_r_comparator(const void *v1, const void *v2)
{
    const struct Linkage *p1 = (struct Linkage *)v1;
    const struct Linkage *p2 = (struct Linkage *)v2;
    if (p1->r < p2->r)
        return -1;
    else if (p1->r > p2->r)
        return +1;
    else
        return 0;
}
