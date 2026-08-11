//
//  SMC_lib.h
//  
//
//  Created by elise.kerdoncuff on 15/12/2017.
//

#ifndef SMC_lib_h
#define SMC_lib_h


#include <stdio.h>

//How we stock node information to build a tree
struct InfosIndiv
{
    int desc1; // -1 if it's a leave
    int desc2;
    int anc; //-1 if it's the root
    float tps; //time of the node du point (0 for the leaves, Tmax for the root)
    float lbs; //branch length above the node
    int nb_leaves;
    float lb1; //total branch length below the node (below desc1)
    float lb2; //total branch length below the node (below desc2)
};


//To stock distribution
struct Distri
{
    float bin_inf; //inferior limite of the bin (cumulated)
    float frequency;
};

//To stock Linkage statitstics
struct Linkage
{
    double d;
    double d_p;
    double r;
    double r2;
    double CP;
    double pqA_pqB;
    double frequency;
    double lbranch_A;
    double lbranch_B;
    double tot_lbranch_A;
    double tot_lbranch_B;
    //float distance_norm;
};

//Use to compare numbers in
struct Two_Floats
{
    float nb1;
    float nb2;
} Two_Floats;


//
///
// Trees and ARG functions
///
//


struct InfosIndiv * Init_ARG(struct InfosIndiv story[],int N,double tchang,float contract);
//function which generates an ARG from a struct InfosIndiv table of the size of the nodes, the time at which the pop size changes and the force of contraction

struct InfosIndiv * SMC_tree(struct InfosIndiv story[],int N,double tchang,float contrac);
//function to perform a SMC algorithm from a tree structure


int sumTable(int table[], int tailleTable);
//function to compute the sum of a int table

int* leavesof(int table[], int point, struct InfosIndiv InfosIndiv[]);
//function which fills an array of size nb leaf with a "1" when the box corresponds to a leaf of the point entered in the function

void tricroiss( float tab[], int tab_size);
//allows you to sort an array of floats

//// Functions for fast sorting tri_fusion
void fusion(float table[],int deb1,int fin1,int fin2);
void tri_fusion_bis(float table[],int deb,int fin);
void tri_fusion(float table[],int longueur);

float sumlbs_To(struct InfosIndiv InfosIndiv[], int N);
// function which sums all the branch lengths

int * descde(int  table[], int point, struct InfosIndiv InfosIndiv[]);
// function which returns an array of size N*2-1 with 0s if the node is not a descendant of the point entered in the function and a 1 if it is


int MRCA(int ind1, int ind2, struct InfosIndiv * story, int N);
// function that returns the most recent common ancestor of two individuals (its number in the table)

char * PrintTree(int indiv, struct InfosIndiv InfosIndiv[],int N, char *strtree, int * pos);
//function which allows you to write a tree in the newick format, which takes as input an individual (initially the root), the structure containing the information of the tree, a pointer to the character string and a pointer to the position where we write in the string


int Two_Floats_nb2_comparator(const void *v1, const void *v2);
//Compares two positions in a struct table Two_Floats (by nb2) to be able to subsequently sort them with qsort


int cmp_tab_egal(int *tab1,int *tab2,int taille);
//function compares table of int one to one, return 1 if different otherwise return 0

int cmp_tabdetabe(int **tab1, int **tab2, int taille, int N);
//function that compares two tables of tables, return 1 if different, 0 if similar

int find_root(struct InfosIndiv * story, int N);
//Function to know the root of the tree, from a tree matrix and the number of leaves


//
////
// LD
////
//

////
/// Frequency functions
////
//Function that returns the frequency of a branch of a tree (frequency if a mutation falls on it)
float frequency_mutations(struct InfosIndiv * story,int indiv, int N);

//Function that returns the joint frequency of two branches of a tree (frequency if mutations fall on it)
float frequency_mutationAB(struct InfosIndiv * story,int indiv1, int indiv2, int N);

//Function that returns the joined frequency of two branches TWO DIFFERENCE TREES (frequency if mutations fall on them)
float frequency_mutationAB_twotrees(struct InfosIndiv * story1,int indiv1,struct InfosIndiv * story2, int indiv2, int N);

//SAME MAF
//Function that returns the joint frequency of two branches of a tree (frequency if mutations fall on it)
float frequency_mutationAB_MAF(struct InfosIndiv * story,int indiv1, int indiv2, int N);
//Function that returns the joined frequency of two branches TWO DIFFERENCE TREES (frequency if mutations fall on them)
float frequency_mutationAB_twotrees_MAF(struct InfosIndiv * story1,int indiv1,struct InfosIndiv * story2, int indiv2, int N);

////
///LD statistcs
////
//Function to calculate the values ​​of D, Dp and r2 on a single tree, returns a Linkage struct table containing the values ​​of D and its proportion
//Input, the tree matrix and the number of leaves
struct Linkage * compute_D_ontree(struct InfosIndiv * story, int N);
//Function to calculate the values ​​of D, Dp and r2 on TWO DIFFERENT TREES, returns a Linkage struct table including the values ​​of D and its proportion
//Input, the tree matrices and the number of leaves
struct Linkage * compute_D_twotrees(struct InfosIndiv * story1,struct InfosIndiv * story2 , int N);
//SAME MAF
//Function to calculate the values ​​of D, Dp and r2 on a single tree, returns a Linkage struct table containing the values ​​of D and its proportion
//Input, the tree matrix and the number of leaves
struct Linkage * compute_D_ontree_MAF(struct InfosIndiv * story, int N);
//Function to calculate the values ​​of D, Dp and r2 on TWO DIFFERENT TREES, returns a Linkage struct table including the values ​​of D and its proportion
//Input, the tree matrices and the number of leaves
struct Linkage * compute_D_twotrees_MAF(struct InfosIndiv * story1,struct InfosIndiv * story2 , int N);


//Compute theoretical computing of Var D0 for one tree
float calcul_varD0(struct InfosIndiv * story, int N);


////
///Comparison functions
////
// Compare two values ​​of .D in a strcut linkage table to then be able to sort them with qsort
int linkD_D_comparator(const void *v1, const void *v2);
// Compare two values ​​of .D_p in a strcut linkage table to then be able to sort them with qsort
int linkD_D_p_comparator(const void *v1, const void *v2);
// Compare two values ​​of .r in a strcut linkage table to then be able to sort them with qsort
int linkD_r_comparator(const void *v1, const void *v2);
// Compare two values ​​of .r2 in a strcut linkage table to then be able to sort them with qsort
int linkD_r2_comparator(const void *v1, const void *v2);


#endif /* D_etude_lib_h */
