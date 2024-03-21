#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "seq_utils.h"
#define MINF -1.0/0.0
#define MAX_CHAR 100000

// macro max
#define MAX(x, y) (x->score >= y->score ? x:y)


char nuc_char[] = {'A','C','T','G'};
char comp_nuc_char[] = {'T','G','A','C'};
int comp_nuc_int[] = {2, 3, 0, 1};

// Gencode & similarity matrix -------------------------------------------------
int gencode[4][4][4] = {{
                            {LYS, ASN, ASN, LYS},
                            {THR, THR, THR, THR},
                            {ILE, ILE, ILE, MET},
                            {ARG, SER, SER, ARG},
                           }, {
                            {GLN, HIS, HIS, GLN},
                            {PRO, PRO, PRO, PRO},
                            {LEU, LEU, LEU, LEU},
                            {ARG, ARG, ARG, ARG},
                           }, {
                            {STOP, TYR, TYR, STOP},
                            {SER, SER, SER, SER},
                            {LEU, PHE, PHE, LEU},
                            {STOP, CYS, CYS, TRP},
                           }, {
                            {GLU, ASP, ASP, GLU},
                            {ALA, ALA, ALA, ALA},
                            {VAL, VAL, VAL, VAL},
                            {GLY, GLY, GLY, GLY},
                           },};

float MATRIX[5][21][21] = {
                           //Identity
                           {
{1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0, MINF},
{MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, MINF, 1.0}
                           },
                           //BLOSUM62
                           {
{  4,   0,   0,  -1,  -2,  -2,  -3,  -1,   0,  -1,  -1,  -1,  -2,  -1,   1,  -1,  -2,  -2,  -1,   0, -10},
{  0,   9,  -1,  -4,  -3,  -2,  -2,  -1,  -1,  -1,  -3,  -1,  -3,  -3,  -1,  -3,  -2,  -3,  -3,  -3, -10},
{  0,  -1,   5,  -1,  -1,  -2,  -2,  -1,   0,  -1,  -1,  -1,   0,  -1,   1,  -1,  -2,  -2,  -1,  -2, -10},
{ -1,  -4,  -1,   5,   2,  -3,  -3,  -3,  -2,  -3,   1,  -2,   0,   2,   0,   0,  -2,   0,  -1,  -2, -10},
{ -2,  -3,  -1,   2,   6,  -3,  -4,  -3,  -3,  -4,  -1,  -3,   1,   0,   0,  -2,  -3,  -1,  -1,  -1, -10},
{ -2,  -2,  -2,  -3,  -3,   6,   1,   0,  -1,   0,  -3,   0,  -3,  -3,  -2,  -3,   3,  -1,  -4,  -3, -10},
{ -3,  -2,  -2,  -3,  -4,   1,  11,  -3,  -3,  -2,  -3,  -1,  -4,  -2,  -3,  -3,   2,  -2,  -4,  -2, -10},
{ -1,  -1,  -1,  -3,  -3,   0,  -3,   4,   3,   2,  -3,   1,  -3,  -3,  -2,  -3,  -1,  -3,  -3,  -4, -10},
{  0,  -1,   0,  -2,  -3,  -1,  -3,   3,   4,   1,  -2,   1,  -3,  -2,  -2,  -3,  -1,  -3,  -2,  -3, -10},
{ -1,  -1,  -1,  -3,  -4,   0,  -2,   2,   1,   4,  -2,   2,  -3,  -2,  -2,  -2,  -1,  -3,  -3,  -4, -10},
{ -1,  -3,  -1,   1,  -1,  -3,  -3,  -3,  -2,  -2,   5,  -1,   0,   1,   0,   2,  -2,  -1,  -1,  -2, -10},
{ -1,  -1,  -1,  -2,  -3,   0,  -1,   1,   1,   2,  -1,   5,  -2,   0,  -1,  -1,  -1,  -2,  -2,  -3, -10},
{ -2,  -3,   0,   0,   1,  -3,  -4,  -3,  -3,  -3,   0,  -2,   6,   0,   1,   0,  -2,   1,  -2,   0, -10},
{ -1,  -3,  -1,   2,   0,  -3,  -2,  -3,  -2,  -2,   1,   0,   0,   5,   0,   1,  -1,   0,  -1,  -2, -10},
{  1,  -1,   1,   0,   0,  -2,  -3,  -2,  -2,  -2,   0,  -1,   1,   0,   4,  -1,  -2,  -1,  -1,   0, -10},
{ -1,  -3,  -1,   0,  -2,  -3,  -3,  -3,  -3,  -2,   2,  -1,   0,   1,  -1,   5,  -2,   0,  -2,  -2, -10},
{ -2,  -2,  -2,  -2,  -3,   3,   2,  -1,  -1,  -1,  -2,  -1,  -2,  -1,  -2,  -2,   7,   2,  -3,  -3, -10},
{ -2,  -3,  -2,   0,  -1,  -1,  -2,  -3,  -3,  -3,  -1,  -2,   1,   0,  -1,   0,   2,   8,  -2,  -2, -10},
{ -1,  -3,  -1,  -1,  -1,  -4,  -4,  -3,  -2,  -3,  -1,  -2,  -2,  -1,  -1,  -2,  -3,  -2,   7,  -2, -10},
{  0,  -3,  -2,  -2,  -1,  -3,  -2,  -4,  -3,  -4,  -2,  -3,   0,  -2,   0,  -2,  -3,  -2,  -2,   6, -10},
{-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, 100}
                           },
                           //BLOSUM80
                           {
{  5,  -1,   0,  -1,  -2,  -3,  -3,  -2,   0,  -2,  -1,  -1,  -2,  -1,   1,  -2,  -2,  -2,  -1,   0,  -6},
{ -1,   9,  -1,  -5,  -4,  -3,  -3,  -2,  -1,  -2,  -4,  -2,  -3,  -4,  -2,  -4,  -3,  -4,  -4,  -4,  -6},
{  0,  -1,   5,  -1,  -1,  -2,  -4,  -1,   0,  -2,  -1,  -1,   0,  -1,   1,  -1,  -2,  -2,  -2,  -2,  -6},
{ -1,  -5,  -1,   6,   1,  -4,  -4,  -4,  -3,  -4,   1,  -2,  -1,   2,   0,  -1,  -3,   0,  -2,  -3,  -6},
{ -2,  -4,  -1,   1,   6,  -4,  -6,  -4,  -4,  -5,  -1,  -4,   1,  -1,  -1,  -2,  -4,  -2,  -2,  -2,  -6},
{ -3,  -3,  -2,  -4,  -4,   6,   0,  -1,  -1,   0,  -4,   0,  -4,  -4,  -3,  -4,   3,  -2,  -4,  -4,  -6},
{ -3,  -3,  -4,  -4,  -6,   0,  11,  -3,  -3,  -2,  -4,  -2,  -4,  -3,  -4,  -4,   2,  -3,  -5,  -4,  -6},
{ -2,  -2,  -1,  -4,  -4,  -1,  -3,   5,   3,   1,  -3,   1,  -4,  -3,  -3,  -3,  -2,  -4,  -4,  -5,  -6},
{  0,  -1,   0,  -3,  -4,  -1,  -3,   3,   4,   1,  -3,   1,  -4,  -3,  -2,  -3,  -2,  -4,  -3,  -4,  -6},
{ -2,  -2,  -2,  -4,  -5,   0,  -2,   1,   1,   4,  -3,   2,  -4,  -3,  -3,  -3,  -2,  -3,  -3,  -4,  -6},
{ -1,  -4,  -1,   1,  -1,  -4,  -4,  -3,  -3,  -3,   5,  -2,   0,   1,  -1,   2,  -3,  -1,  -1,  -2,  -6},
{ -1,  -2,  -1,  -2,  -4,   0,  -2,   1,   1,   2,  -2,   6,  -3,   0,  -2,  -2,  -2,  -2,  -3,  -4,  -6},
{ -2,  -3,   0,  -1,   1,  -4,  -4,  -4,  -4,  -4,   0,  -3,   6,   0,   0,  -1,  -3,   0,  -3,  -1,  -6},
{ -1,  -4,  -1,   2,  -1,  -4,  -3,  -3,  -3,  -3,   1,   0,   0,   6,   0,   1,  -2,   1,  -2,  -2,  -6},
{  1,  -2,   1,   0,  -1,  -3,  -4,  -3,  -2,  -3,  -1,  -2,   0,   0,   5,  -1,  -2,  -1,  -1,  -1,  -6},
{ -2,  -4,  -1,  -1,  -2,  -4,  -4,  -3,  -3,  -3,   2,  -2,  -1,   1,  -1,   6,  -3,   0,  -2,  -3,  -6},
{ -2,  -3,  -2,  -3,  -4,   3,   2,  -2,  -2,  -2,  -3,  -2,  -3,  -2,  -2,  -3,   7,   2,  -4,  -4,  -6},
{ -2,  -4,  -2,   0,  -2,  -2,  -3,  -4,  -4,  -3,  -1,  -2,   0,   1,  -1,   0,   2,   8,  -3,  -3,  -6},
{ -1,  -4,  -2,  -2,  -2,  -4,  -5,  -4,  -3,  -3,  -1,  -3,  -3,  -2,  -1,  -2,  -4,  -3,   8,  -3,  -6},
{  0,  -4,  -2,  -3,  -2,  -4,  -4,  -5,  -4,  -4,  -2,  -4,  -1,  -2,  -1,  -3,  -4,  -3,  -3,   6,  -6},
{ -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,   1}
                           },
                           // BLOSUM90
                           {
{  5,  -1,   0,  -1,  -3,  -3,  -4,  -2,  -1,  -2,  -1,  -2,  -2,  -1,   1,  -2,  -3,  -2,  -1,   0,  -6},
{ -1,   9,  -2,  -6,  -5,  -3,  -4,  -2,  -2,  -2,  -4,  -2,  -4,  -4,  -2,  -5,  -4,  -5,  -4,  -4,  -6},
{  0,  -2,   6,  -1,  -2,  -3,  -4,  -1,  -1,  -2,  -1,  -1,   0,  -1,   1,  -2,  -2,  -2,  -2,  -3,  -6},
{ -1,  -6,  -1,   6,   1,  -5,  -5,  -4,  -3,  -4,   0,  -3,  -1,   2,  -1,  -1,  -4,  -1,  -2,  -3,  -6},
{ -3,  -5,  -2,   1,   7,  -5,  -6,  -5,  -5,  -5,  -1,  -4,   1,  -1,  -1,  -3,  -4,  -2,  -3,  -2,  -6},
{ -3,  -3,  -3,  -5,  -5,   7,   0,  -1,  -2,   0,  -4,  -1,  -4,  -4,  -3,  -4,   3,  -2,  -4,  -5,  -6},
{ -4,  -4,  -4,  -5,  -6,   0,  11,  -4,  -3,  -3,  -5,  -2,  -5,  -3,  -4,  -4,   2,  -3,  -5,  -4,  -6},
{ -2,  -2,  -1,  -4,  -5,  -1,  -4,   5,   3,   1,  -4,   1,  -4,  -4,  -3,  -4,  -2,  -4,  -4,  -5,  -6},
{ -1,  -2,  -1,  -3,  -5,  -2,  -3,   3,   5,   0,  -3,   0,  -4,  -3,  -2,  -3,  -3,  -4,  -3,  -5,  -6},
{ -2,  -2,  -2,  -4,  -5,   0,  -3,   1,   0,   5,  -3,   2,  -4,  -3,  -3,  -3,  -2,  -4,  -4,  -5,  -6},
{ -1,  -4,  -1,   0,  -1,  -4,  -5,  -4,  -3,  -3,   6,  -2,   0,   1,  -1,   2,  -3,  -1,  -2,  -2,  -6},
{ -2,  -2,  -1,  -3,  -4,  -1,  -2,   1,   0,   2,  -2,   7,  -3,   0,  -2,  -2,  -2,  -3,  -3,  -4,  -6},
{ -2,  -4,   0,  -1,   1,  -4,  -5,  -4,  -4,  -4,   0,  -3,   7,   0,   0,  -1,  -3,   0,  -3,  -1,  -6},
{ -1,  -4,  -1,   2,  -1,  -4,  -3,  -4,  -3,  -3,   1,   0,   0,   7,  -1,   1,  -3,   1,  -2,  -3,  -6},
{  1,  -2,   1,  -1,  -1,  -3,  -4,  -3,  -2,  -3,  -1,  -2,   0,  -1,   5,  -1,  -3,  -2,  -2,  -1,  -6},
{ -2,  -5,  -2,  -1,  -3,  -4,  -4,  -4,  -3,  -3,   2,  -2,  -1,   1,  -1,   6,  -3,   0,  -3,  -3,  -6},
{ -3,  -4,  -2,  -4,  -4,   3,   2,  -2,  -3,  -2,  -3,  -2,  -3,  -3,  -3,  -3,   8,   1,  -4,  -5,  -6},
{ -2,  -5,  -2,  -1,  -2,  -2,  -3,  -4,  -4,  -4,  -1,  -3,   0,   1,  -2,   0,   1,   8,  -3,  -3,  -6},
{ -1,  -4,  -2,  -2,  -3,  -4,  -5,  -4,  -3,  -4,  -2,  -3,  -3,  -2,  -2,  -3,  -4,  -3,   8,  -3,  -6},
{  0,  -4,  -3,  -3,  -2,  -5,  -4,  -5,  -5,  -5,  -2,  -4,  -1,  -3,  -1,  -3,  -5,  -3,  -3,   6,  -6},
{ -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,   1}
                           },
                           //VTML200
                           {
{  4,   1,   1,  -1,  -1,  -3,  -4,  -1,   0,  -2,  -1,  -1,  -1,  -1,   1,  -2,  -3,  -2,   0,   0,  -6},
{  1,  12,   0,  -4,  -4,  -3,  -6,   0,   1,  -3,  -4,  -1,  -2,  -3,   1,  -3,   0,  -2,  -3,  -2,  -6},
{  1,   0,   4,  -1,  -1,  -3,  -5,  -1,   0,  -2,   0,  -1,   0,   0,   2,  -1,  -3,  -1,  -1,  -2,  -6},
{ -1,  -4,  -1,   5,   3,  -5,  -6,  -4,  -3,  -4,   1,  -3,   1,   2,   0,  -1,  -3,   0,  -1,  -1,  -6},
{ -1,  -4,  -1,   3,   6,  -6,  -6,  -5,  -4,  -5,   0,  -4,   3,   1,   0,  -2,  -4,   0,  -1,  -1,  -6},
{ -3,  -3,  -3,  -5,  -6,   8,   3,   0,  -1,   2,  -5,   1,  -4,  -3,  -3,  -4,   5,   0,  -4,  -5,  -6},
{ -4,  -6,  -5,  -6,  -6,   3,  15,  -2,  -4,  -1,  -4,  -3,  -5,  -6,  -4,  -3,   4,  -1,  -4,  -5,  -6},
{ -1,   0,  -1,  -4,  -5,   0,  -2,   5,   4,   3,  -3,   2,  -4,  -3,  -3,  -3,  -2,  -3,  -4,  -6,  -6},
{  0,   1,   0,  -3,  -4,  -1,  -4,   4,   4,   2,  -3,   2,  -3,  -2,  -2,  -3,  -2,  -3,  -3,  -4,  -6},
{ -2,  -3,  -2,  -4,  -5,   2,  -1,   3,   2,   5,  -3,   3,  -4,  -2,  -3,  -3,  -1,  -2,  -3,  -5,  -6},
{ -1,  -4,   0,   1,   0,  -5,  -4,  -3,  -3,  -3,   5,  -2,   1,   2,   0,   4,  -3,   0,  -1,  -2,  -6},
{ -1,  -1,  -1,  -3,  -4,   1,  -3,   2,   2,   3,  -2,   6,  -3,  -1,  -2,  -2,  -2,  -3,  -3,  -4,  -6},
{ -1,  -2,   0,   1,   3,  -4,  -5,  -4,  -3,  -4,   1,  -3,   6,   1,   1,   0,  -2,   1,  -2,   0,  -6},
{ -1,  -3,   0,   2,   1,  -3,  -6,  -3,  -2,  -2,   2,  -1,   1,   5,   0,   2,  -3,   2,  -1,  -2,  -6},
{  1,   1,   2,   0,   0,  -3,  -4,  -3,  -2,  -3,   0,  -2,   1,   0,   4,  -1,  -2,   0,   0,   0,  -6},
{ -2,  -3,  -1,  -1,  -2,  -4,  -3,  -3,  -3,  -3,   4,  -2,   0,   2,  -1,   7,  -2,   1,  -1,  -2,  -6},
{ -3,   0,  -3,  -3,  -4,   5,   4,  -2,  -2,  -1,  -3,  -2,  -2,  -3,  -2,  -2,   9,   3,  -5,  -5,  -6},
{ -2,  -2,  -1,   0,   0,   0,  -1,  -3,  -3,  -2,   0,  -3,   1,   2,   0,   1,   3,   8,  -2,  -2,  -6},
{  0,  -3,  -1,  -1,  -1,  -4,  -4,  -4,  -3,  -3,  -1,  -3,  -2,  -1,   0,  -1,  -5,  -2,   9,  -2,  -6},
{  0,  -2,  -2,  -1,  -1,  -5,  -5,  -6,  -4,  -5,  -2,  -4,   0,  -2,   0,  -2,  -5,  -2,  -2,   8,  -6},
{ -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,   1}
                           }
};

// Functions -------------------------------------------------------------------
float compute_score_sim_local(pair_p init, pair_p prime, int mat_type, float thr_scr) {
  float score_x = 0.0, score_y = 0.0;
  if (init->x != GAP)
    score_x = MATRIX[mat_type][init->x][prime->x];
  if (init->y != GAP)
    score_y = MATRIX[mat_type][init->y][prime->y];
  if (score_x < thr_scr || score_y < thr_scr)
    return MINF;
  else
    return 0.5 * (score_x + score_y);
};

float compute_score_sim(pair_p init, pair_p prime, int mat_type) {
  float score_x = 0.0, score_y = 0.0;
  if (init->x != GAP)
    score_x = MATRIX[mat_type][init->x][prime->x];
  if (init->y != GAP)
    score_y = MATRIX[mat_type][init->y][prime->y];
  return 0.5 * (score_x + score_y);
};


// convert_aa string to enum amino acids
Aminoacid convert_aa(char x) {
  char amino_acids_char[] = {'A' ,'C' ,'T' ,'E' ,'D' ,'F' ,'W' ,'I' ,'V' ,'L' ,'K' ,'M' ,'N' ,'Q' ,'S' ,'R' ,'Y' ,'H' ,'P' ,'G', '*', '-'};
  for (Aminoacid i = ALA; i <= GAP; i++) {
    if (amino_acids_char[i] == x)
      return i;
  }
  printf("Error in the AA sequence %c\n", x);
  exit(1);
};

// convert string to enum nucleotides
Nucleotide convert_nuc(char x) {
  for (Nucleotide i = A; i <= G; i++) {
    if (nuc_char[i] == x)
      return i;
  }
  printf("Error in the NUC sequence\n");
  exit(1);
};

// translation of a quartet
Aminoacid* translate(quartet_p quart, int* frames) {
  int f;
  Aminoacid* aa = (Aminoacid*)malloc(sizeof(Aminoacid)*2);
  for (f = 0; f < 2; f++) {
    switch (frames[f]) {
    case 1:
      aa[f] = (Aminoacid)gencode[quart->kmer[0]][quart->kmer[1]][quart->kmer[2]];
      break;
    case 2:
      aa[f] = (Aminoacid)gencode[quart->kmer[1]][quart->kmer[2]][quart->kmer[3]];
      break;
    case -1:
      aa[f] = (Aminoacid)gencode[quart->ckmer[2]][quart->ckmer[1]][quart->ckmer[0]];
      break;
    case -2:
      aa[f] = (Aminoacid)gencode[quart->ckmer[3]][quart->ckmer[2]][quart->ckmer[1]];
      break;
    }
  }
  return aa;
};

// get the maximum score
quartet_p get_max_quartet(quartet_p substitution, quartet_p deletion, quartet_p insertion) {
  quartet_p tmp;
  tmp = MAX(substitution, deletion);
  tmp = MAX(tmp, insertion);
  return tmp;
}

// copy quartet infos
quartet_p copy_to_max(quartet_p from_quartet, quartet_p to_quartet) {
  // copy values values
  from_quartet->init->x = to_quartet->init->x;
  from_quartet->init->y = to_quartet->init->y;
  from_quartet->prime->x = to_quartet->prime->x;
  from_quartet->prime->y = to_quartet->prime->y;
  from_quartet->score = to_quartet->score;
  from_quartet->previous = to_quartet->previous;
  from_quartet->gap = to_quartet->gap;
  from_quartet->kmer = to_quartet->kmer;
  from_quartet->last_gap = to_quartet->last_gap;
  return from_quartet;
}

// Initialize quartet for kmer list
void init_quartet(quartet_p quart, int i, int j, int k, int l, int frames[]) {
  Aminoacid* pair_aa;
  quart->kmer = (Nucleotide*)malloc(sizeof(Nucleotide)*4);
  quart->ckmer = (Nucleotide*)malloc(sizeof(Nucleotide)*4);
  quart->prime = (pair_p)malloc(sizeof(pair));
  quart->kmer[0] = i; quart->kmer[1] = j; quart->kmer[2] = k; quart->kmer[3] = l;
  quart->ckmer[0] = comp_nuc_int[i]; quart->ckmer[1] = comp_nuc_int[j];
  quart->ckmer[2] = comp_nuc_int[k]; quart->ckmer[3] = comp_nuc_int[l];

  // translate
  pair_aa = translate(quart, frames);
  quart->prime->x = pair_aa[0]; quart->prime->y = pair_aa[1];
  free(pair_aa);
};

// Initialize all the kmers, return an array of 256 quartets
quartet_p**** get_all_kmers(int* frames) {
  quartet_p**** all_kmers = (quartet_p****)malloc(sizeof(quartet_p)*256);
  int i, j, k, l, N=4;

  // initialize all the 256 quartets
  for (i = A; i < N; ++i) {
    all_kmers[i] = (quartet_p***)malloc(sizeof(quartet_p)*64);
    for (j = A; j < N; ++j) {
      all_kmers[i][j] = (quartet_p**)malloc(sizeof(quartet_p)*16);
      for (k = A; k < N; ++k) {
        all_kmers[i][j][k] = (quartet_p*)malloc(sizeof(quartet_p)*4);
        for (l = A; l < N; ++l) {
          all_kmers[i][j][k][l] = (quartet_p)malloc(sizeof(quartet));
          init_quartet(all_kmers[i][j][k][l], i, j, k, l, frames);
        }
      }
    }
  }
  return all_kmers;
};

// free the list of 256 kmers
void free_kmer_list(quartet_p**** kmer_list) {
  int i, j, k, l, N=4;

  // initialize all the 256 quartets
  for (i = A; i < N; ++i) {
    for (j = A; j < N; ++j) {
      for (k = A; k < N; ++k) {
        for (l = A; l < N; ++l) {
          free(kmer_list[i][j][k][l]->kmer);
          free(kmer_list[i][j][k][l]->ckmer);
          free(kmer_list[i][j][k][l]->prime);
          free(kmer_list[i][j][k][l]);
        }
        free(kmer_list[i][j][k]);
      }
      free(kmer_list[i][j]);
    }
    free(kmer_list[i]);
  }
  free(kmer_list);
};

// Initialize an empty quartet
// takes a pair of amino acids and return a quartet
quartet_p init_empty_quartet(Aminoacid x, Aminoacid y) {
  quartet_p new_quartet = (quartet_p)malloc(sizeof(quartet));
  new_quartet->init = malloc(sizeof(pair));
  new_quartet->prime = malloc(sizeof(pair));
  new_quartet->init->x = x;
  new_quartet->init->y = y;
  return new_quartet;
};

// Initialze the matrix
// for each entry, we have a list of 4 quartet_p, max for each kind of
// Nucleotide
matrix_qp init_matrix(char* seq_x, char* seq_y) {
  matrix_qp matrix = malloc(sizeof(matrix_q));
  int len_x = strlen(seq_x), len_y = strlen(seq_y);
  int i, j, nuc;

  // save lengths
  matrix->len_x = len_x;
  matrix->len_y = len_y;

  // allocate memory for the matrix now
  matrix->matrix_entries = malloc(sizeof(matrix_entry_p)*(len_x+1)*(len_y+1));

  // initialize all the quartet entries
  for (i = 0; i <= len_x; ++i) {
    matrix->matrix_entries[i] = malloc(sizeof(matrix_entry_p)*(len_y+1));
    for (j = 0; j <= len_y; ++j) {
      // initialize an array of 4 quartets
      matrix->matrix_entries[i][j] = malloc(sizeof(matrix_entry_p));
      matrix->matrix_entries[i][j]->posx = i-1; matrix->matrix_entries[i][j]->posy = j-1;
      matrix->matrix_entries[i][j]->quartet_entries = (quartet_p*)malloc(sizeof(quartet_p)*4);
      if (i > 0 && j > 0) {
        // for each element, initialize an quartet
        for (nuc = A; nuc <= G; nuc++)
          matrix->matrix_entries[i][j]->quartet_entries[nuc] = init_empty_quartet(convert_aa(seq_x[i-1]), convert_aa(seq_y[j-1]));

      } else {
        // for each element, initialize an quartet
        for (nuc = A; nuc <= G; nuc++) {
          if (i == 0 && j == 0) {
            matrix->matrix_entries[i][j] = NULL;
          } else if (i == 0)
            matrix->matrix_entries[i][j]->quartet_entries[nuc] = init_empty_quartet(GAP, convert_aa(seq_y[j-1]));
          else if (j == 0)
            matrix->matrix_entries[i][j]->quartet_entries[nuc] = init_empty_quartet(convert_aa(seq_x[i-1]), GAP);
        }
      }
    }
  }
  return matrix;
};

// free the dynamic programming matrix
void free_matrix(matrix_qp matrix) {
  int i, j, nuc;
  for (i = 0; i <= matrix->len_x; ++i) {
    for (j = 0; j <= matrix->len_y; ++j) {
      if (i > 0 || j > 0) {
        for (nuc = A; nuc <= G; nuc++) {
          free(matrix->matrix_entries[i][j]->quartet_entries[nuc]->init);
          free(matrix->matrix_entries[i][j]->quartet_entries[nuc]->prime);
          free(matrix->matrix_entries[i][j]->quartet_entries[nuc]);
        }
        free(matrix->matrix_entries[i][j]->quartet_entries);
      }
      free(matrix->matrix_entries[i][j]);
    }
    free(matrix->matrix_entries[i]);
  }
  free(matrix->matrix_entries);
  free(matrix);
};

// dca scores, lower is better -> takes the opposite of the dca scores
float compute_score_dca(pair_p init, pair_p prime, float**** mat_dca_x, float**** mat_dca_y, int pos_x, int pos_y) {
  float score_x = 0.0, score_y = 0.0;
  if (init->x != GAP)
    score_x = mat_dca_x[pos_x][pos_x][prime->x][prime->x];
  if (init->y != GAP)
    score_y = mat_dca_y[pos_y][pos_y][prime->y][prime->y];
  return -0.5 * (score_x + score_y);
};

// dca scores, lower is better -> takes the opposite of the dca scores
float compute_score_dca_pair(pair_p cur_init, pair_p cur_prime, pair_p prev_init, pair_p prev_prime,
                             float**** mat_dca_x, float**** mat_dca_y, matrix_entry_p cur_mat,
                             matrix_entry_p prev_mat) {
  float score_x = 0.0, score_y = 0.0;
  int cposx = cur_mat->posx, cposy = cur_mat->posy;
  int pposx = prev_mat->posx, pposy = prev_mat->posy;
  if (cur_init->x != GAP && prev_init->x != GAP)
    score_x = mat_dca_x[pposx][cposx][prev_prime->x][cur_prime->x];
  if (cur_init->y != GAP && prev_init->y != GAP)
    score_y = mat_dca_y[pposy][cposy][prev_prime->y][cur_prime->y];
  return -0.5 * (score_x + score_y);
};

// Initialize all the matrix score
float**** read_score(int len_seq, char* dca_score_file, int frame) {
  float**** all_dca = (float****)malloc(sizeof(float***) * len_seq);
  FILE *dca_file;
  char line[MAX_CHAR];
  char typeiS, typejS;
  float scoreij;
  int posi, posj, typei, typej, i, j, ai, aj;
  if ((dca_file = fopen(dca_score_file, "r")) == NULL) {
    printf("Error! opening file");
    exit(1);
  }

  for (i = 0; i < len_seq; i++) {
    all_dca[i] = (float***)malloc(sizeof(float**) * len_seq);
    for (j = 0; j < len_seq; j++) {
      all_dca[i][j] = (float**)malloc(sizeof(float*) * 22);
      for (ai = ALA; ai <= GAP; ai++) {
        all_dca[i][j][ai] = (float*)malloc(sizeof(float) * 22);
        for (aj = ALA; aj <= GAP; aj++) {
          if (ai == STOP || aj == STOP)
            all_dca[i][j][ai][aj] = INFINITY;
          else
            all_dca[i][j][ai][aj] = 0.0;
        }
      }
    }
  }

  // loop in the bias file
  while(fgets(line, MAX_CHAR, dca_file) != NULL) {
    if (sscanf(line,"%d %d %c %c %f", &posi, &posj, &typeiS, &typejS, &scoreij) == 5) {
      typei = convert_aa(typeiS); typej = convert_aa(typejS);
      if (typei != GAP && typej != GAP) {
        if (frame < 0) {
          posi = len_seq - posi - 1;
          posj = len_seq - posj - 1;
        }
        all_dca[posi][posj][typei][typej] = scoreij;
        all_dca[posj][posi][typej][typei] = scoreij;
      }
    }
  }

  fclose(dca_file);
  return all_dca;
};

// free the dca scores
void free_score(float**** dca_score, int len_seq) {
  int i, j, ai;
  for (i = 0; i < len_seq; i++) {
    for (j = 0; j < len_seq; j++) {
      for (ai = ALA; ai <= GAP; ai++) {
        free(dca_score[i][j][ai]);
      }
      free(dca_score[i][j]);
    }
    free(dca_score[i]);
  }
  free(dca_score);
}

// reverse a string
void reverse(char *x, int begin, int end) {
  char c;
  if (begin >= end)
    return;

  c          = *(x+begin);
  *(x+begin) = *(x+end);
  *(x+end)   = c;

  reverse(x, ++begin, --end);
}

// read a fasta file
fasta_seq_p read_fasta(char* fasta_file_str) {
  FILE *fasta_file;
  /* char line[MAX_CHAR]; */
  char* line = (char*)malloc(sizeof(char)*MAX_CHAR);
  char* seq = (char*)malloc(sizeof(char)*MAX_CHAR);
  char* name = (char*)malloc(sizeof(char)*MAX_CHAR);
  fasta_seq_p cur_seq = (fasta_seq_p)malloc(sizeof(fasta_seq));

  if ((fasta_file = fopen(fasta_file_str, "r")) == NULL) {
    printf("Error! opening file");
    exit(1);
  }
  while(fgets(line, MAX_CHAR, fasta_file) != NULL) {
    if (line[0] != '>') {
      if (strcmp(seq, "") != 0) {
        seq[strcspn(seq, "\n")] = 0;
        strcat(seq, &line[0]);
      } else {
        strcpy(seq, line);
      }
    } else {
      line++;
      strcpy(name, line);
      name[strcspn(name, "\n")] = 0;
    }
  }
  seq[strcspn(seq, "\n")] = 0;
  fclose(fasta_file);
  cur_seq->name = name;
  cur_seq->seq = seq;
  return cur_seq;
}

// return a random integer in [from_int; to_int] interval
int rand_int(int from_int, int to_int) {
  if (from_int < to_int)
    return rand() % (to_int - from_int) + from_int;
  else
    fprintf(stderr,"ERROR");
  exit(1);
  return 0;
}

// return a random float in [min; max] interval
float rand_float(int min, int max) {
  if (max == min) return min;
  else if (min < max) return (float)rand() / RAND_MAX;;
  // return 0 if min > max
  return 0.0;
}

// get the mutated aa
Aminoacid get_mut_codon(Nucleotide* dna, int frame, int posi, int icod, Nucleotide nnuc) {
  Nucleotide codon[3];
  codon[0] = dna[abs(frame) - 1 + posi*3];
  codon[1] = dna[abs(frame) - 1 + posi*3 + 1];
  codon[2] = dna[abs(frame) - 1 + posi*3 + 2];
  codon[icod] = nnuc;
  if (frame > 0) {
    return (Aminoacid)gencode[codon[0]][codon[1]][codon[2]];
  } else {
    return (Aminoacid)gencode[comp_nuc_int[codon[2]]][comp_nuc_int[codon[1]]][comp_nuc_int[codon[0]]];
  }
}

// return an array of id
int* ungapped_seq_id(Aminoacid* seq, int len_seq) {
  int i, ri = 0;
  int* ugap = (int*)malloc(sizeof(int)*len_seq);
  for (i = 0; i < len_seq; i++) {
    if (seq[i] != GAP) {
      ugap[i] = ri;
      ri++;
    }
  }
  return ugap;
}

void printProgress (int step, int traj_len) {
  int prcent = (int)(((float)step / traj_len) * 100.0);
  printf("\r\e[?25l Progress: %4d %%", prcent);
  if (prcent == 100) {
    printf("\n\e[?25h");
  }
}

int count_gap(char* seq) {
  int count = 0, i;
  for (i = 0; i < strlen(seq); i++) {
    if (seq[i] != '-')
      count++;
  }
  return count;
}

// Print the sequences
// 0 | 1 = initial X | Y
// 2 | 3 = prime X | Y
// 4 = dna sequence
// frame indicates if the sequence is in the opposite sens
// zero indicates the phase = 0
void backtrace_seq(char* tmp, int tmp_id, quartet_p best_quartet, int mode, int frame, Bool zero, Bool glob) {
  char amino_acids_char[] = {'A' ,'C' ,'T' ,'E' ,'D' ,'F' ,'W' ,'I' ,'V' ,'L' ,'K' ,'M' ,'N' ,'Q' ,'S' ,'R' ,'Y' ,'H' ,'P' ,'G', '*', '-'};
  char nuc_2_char[] = {'A','C','T','G'};
  if (best_quartet != NULL && (best_quartet->score > 0 || glob == true)) {
    quartet_p prev_quartet = (quartet_p)best_quartet->previous;
    if (frame > 0) {
      if (mode == 4)
        backtrace_seq(tmp, tmp_id+3, prev_quartet, mode, frame, zero, glob);
      else
        backtrace_seq(tmp, tmp_id+1, prev_quartet, mode, frame, zero, glob);
    }

    if (mode == 0) {
      tmp[tmp_id] = amino_acids_char[best_quartet->init->x];
    }
    if (mode == 1)
      tmp[tmp_id] = amino_acids_char[best_quartet->init->y];
    if (mode == 2) {
      if (best_quartet->init->x != GAP)
        tmp[tmp_id] = amino_acids_char[best_quartet->prime->x];
      else
        tmp[tmp_id] = amino_acids_char[GAP];
    }
    if (mode == 3) {
      if (best_quartet->init->y != GAP)
        tmp[tmp_id] = amino_acids_char[best_quartet->prime->y];
      else
        tmp[tmp_id] = amino_acids_char[GAP];
    }
    if (mode == 4) {
      if ((prev_quartet == NULL || prev_quartet->score <= 0) && zero != true) {
        tmp[tmp_id] = nuc_2_char[best_quartet->kmer[3]];
        tmp[tmp_id+1] = nuc_2_char[best_quartet->kmer[2]];
        tmp[tmp_id+2] = nuc_2_char[best_quartet->kmer[1]];
        tmp[tmp_id+3] = nuc_2_char[best_quartet->kmer[0]];
      } else {
        tmp[tmp_id] = nuc_2_char[best_quartet->kmer[3]];
        tmp[tmp_id+1] = nuc_2_char[best_quartet->kmer[2]];
        tmp[tmp_id+2] = nuc_2_char[best_quartet->kmer[1]];
      }
    }
    if (frame < 0)
      backtrace_seq(tmp, tmp_id+1, prev_quartet, mode, frame, zero, glob);
  }
}


// trace back, reverse read in tmp
void backtrace(result_quartet_p solution, char* seq_x_full, char* seq_y_full, int* frames, Bool glob) {
  char* seq_x_i = (char*)malloc(sizeof(char)*MAX_CHAR*2);
  char* seq_y_i = (char*)malloc(sizeof(char)*MAX_CHAR*2);
  char* seq_x_p = (char*)malloc(sizeof(char)*MAX_CHAR*2);
  char* seq_y_p = (char*)malloc(sizeof(char)*MAX_CHAR*2);
  char* dna = (char*)malloc(sizeof(char)*MAX_CHAR*2*4);
  quartet_p best_quartet = solution->result;
  int len_sol_x, len_sol_y;
  int len_x = strlen(seq_x_full), len_y = strlen(seq_y_full);
  int startx, endx, starty, endy;
  backtrace_seq(seq_x_i, 0, best_quartet, 0, frames[0], false, glob);
  backtrace_seq(seq_y_i, 0, best_quartet, 1, frames[1], false, glob);
  backtrace_seq(seq_x_p, 0, best_quartet, 2, frames[0], false, glob);
  backtrace_seq(seq_y_p, 0,best_quartet, 3, frames[1], false, glob);
  if ((frames[0] == 1 && frames[1] == -1) || (frames[0] == 1 && frames[1] == 1))
    backtrace_seq(dna, 0, best_quartet, 4, 1, true, glob);
  else
    backtrace_seq(dna, 0, best_quartet, 4, 1, false, glob);

  // reverse dna
  reverse(dna, 0, strlen(dna)-1);
  reverse(seq_x_i, 0, strlen(seq_x_i)-1);
  reverse(seq_x_p, 0, strlen(seq_x_p)-1);
  if (frames[1] > 0) {
    reverse(seq_y_i, 0, strlen(seq_y_i)-1);
    reverse(seq_y_p, 0, strlen(seq_y_p)-1);
  }

  len_sol_x = count_gap(seq_x_i);
  len_sol_y = count_gap(seq_y_i);

  startx = solution->posx - len_sol_x;
  endx = solution->posx;

  if (frames[1] > 0) {
    starty = solution->posy - len_sol_y;
    endy = solution->posy;
  } else {
    starty = len_y - solution->posy;
    endy = starty + len_sol_y;
  }

  printf("# STARTX = %d\n# STARTY = %d\n# ENDX = %d\n# ENDY = %d\n", startx, starty, endx, endy);
  printf(">SEQ_X_FULL\n%s\n" ,seq_x_full);
  printf(">SEQ_Y_FULL\n%s\n" ,seq_y_full);
  printf(">SEQ_X_INIT\n");
  printf("%s\n", seq_x_i);
  printf(">SEQ_Y_INIT\n");
  printf("%s\n", seq_y_i);
  printf(">SEQ_X_PRIME\n");
  printf("%s\n", seq_x_p);
  printf(">SEQ_Y_PRIME\n");
  printf("%s\n", seq_y_p);
  printf(">DNA\n");
  printf("%s\n", dna);
  free(seq_x_i); free(seq_y_i); free(seq_x_p); free(seq_y_p); free(dna);
};

// convert string to enum amino acids
int* convert_frame(char* frame) {
  static int frame_id[2];
  if (strcmp(frame, "-2") == 0) {
    frame_id[0] = 1; frame_id[1] = -2;
  }
  if (strcmp(frame, "-1") == 0) {
    frame_id[0] = 2; frame_id[1] = -1;
  }
  if (strcmp(frame, "0") == 0) {
    frame_id[0] = 1; frame_id[1] = -1;
  }
  if (strcmp(frame, "1") == 0) {
    frame_id[0] = 1; frame_id[1] = 2;
  }
  if (strcmp(frame, "2") == 0) {
    frame_id[0] = 2; frame_id[1] = 1;
  }
  if (strcmp(frame, "3") == 0) {
    frame_id[0] = 1; frame_id[1] = 1;
  }
  return frame_id;
};


int convert_matrix(char* matrix_str) {
  int matrix=0;
  if (strcmp(matrix_str, "ident") == 0)
    matrix = 0;
  if (strcmp(matrix_str, "blosum62") == 0)
    matrix = 1;
  if (strcmp(matrix_str, "blosum80") == 0)
    matrix = 2;
  if (strcmp(matrix_str, "blosum90") == 0)
    matrix = 3;
  if (strcmp(matrix_str, "vtml200") == 0)
    matrix = 4;
  return matrix;
}
