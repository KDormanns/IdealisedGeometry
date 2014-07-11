#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream>
#include <algorithm>

using namespace std;

#define RIGHT	0
#define LEFT	1
#define x_coord 0
#define y_coord 1
#define z_coord 2
#define P	0
#define L	1
#define R	2
#define endcaps 3

#define default_NUM_ECs					4
#define default_NUM_SMCs				3
#define default_max_corrections 	5
extern "C" {
void dbihar_(double* a, double* b, int* m, double bda[], double bdb[], double bdc[], double bdd[], double* c, double* d, int* n, double *f, int* idf,
		double* alpha, double* beta, int* iflag, double* tol, int* itcg, double w[], int* lw);
}

#define P2(x) ((x)*(x))
#define pi		3.141592653589793
#define true  1
#define false 0
typedef struct {
	int m, n, nx, ny;
	double a, b, c, d, alpha, beta, rb, rt, apex, neck, length, h1, h2, h3, xshift, yshift, zshift, rotation_angle, apex_scaling;

	double c1, c2, b_bnw, alfa, angle, ss0, ss1, ra, xc, h, r,			// radius of intersecting cylinder where 0.636 <= r <=0.9886
			s0, s1, A1, B1, C1, D1, A2, B2, C2, D2, v1, v2;

	int max_corrections;
	int ECs, SMCs;

} parms;

typedef struct {
	double **x, **y, **z, **cx, **cy, **cz;
} mesh_store;

double** allocate_dirichlet(int rows, int cols);
double* allocate_neuman(int);
void initialize_bd(int indx, double *bd);
void initialize_f(parms* p, double **fx, double **fy, double **fz);
void print_stdout(parms p, double **fx, double **fy, double **fz);

void format_vtk(parms p, double** fx, double** fy, double** fz, char* prefix);
void format_stl(parms p, double*** f, char* prefix, int, int);
void facet1(int i, int j, double*** f, FILE* FileWrite);
void facet2(int i, int j, double*** f, FILE* FileWrite);
double** solve(parms* p, double *bda, double *bdb, double *bdc, double *bdd, double **f, int idf, int iflag, double tol, int itcg, double *w, int lw);

void dirichlet_boundary_parent_segment(parms p, double **fx, double **fy, double **fz);
double* dirichlet_boundary_outer_wall(parms p, double **fx, double **fy, double **fz);
double* dirichlet_boundary_inner_wall(parms p, double **fx, double **fy, double **fz);
void dirichlet_boundary_end_cap_parent(parms p, double **fx, double **fy, double **fz);
void dirichlet_boundary_end_cap_daughter(parms p, double **fx, double **fy, double **fz);
void record_data(int k, parms p, double **fx, double **fy, double **fz, double ****f);
void boundary_correction(parms p, double ****f);

void print_func(int nx, int ny, double **f, char* str);
int ec_mesh(parms p, double** fx, double** fy, double** fz, char*, char*, int, int);
int smc_mesh(parms p, double** fx, double** fy, double** fz, char*, char*, int, int);
void store_data(int m, int n, double **x, double **y, double **z, char *prefix, int indx);
int format_vtk_small_mesh(char *input_file, char *input_file2, char *centeroids, char *output_file, int nx, int ny, int total_points, int, int);
void deallocate(int m, int n, double **f);
bool isEmpty(FILE *file);

void format_vtk_centeroid(char* filename, int count);
void format_vtk_unstructuredGrid(parms p, double** fx, double** fy, double** fz, char* prefix, int, int);
void store_arrays(parms* p, double** fx, double** fy, double** fz, double**** storage, int half);
double**** allocate_storage_array(int nx, int ny);
int* format_primitive(parms p, double**** storage, char *prefix, int downstream_offset, int upstream_offset);
int* smc_mesh_ver2(parms, double****, double***, int, int, int, int, char*);
int* ec_mesh_ver2(parms p, double**** storage, double*** buffer, int Total_pannels_axially, int Total_pannels_circumferentially,
		int downstream_offset, int upstream_offset, char* prefix);

void bound_v_correction_ver2(parms p, double**** storage, int downstream_offset, int upstream_offset);
void bound_v_correction_ver3(parms p, double**** storage, double**, int downstream_offset, int upstream_offset);
double* dirichlet_boundary_outer_wall_corrected_v_boundary(parms* p, double**** storage, double **fx, double **fy, double **fz, int half);

double* dirichlet_boundary_inner_wall_corrected_v_boundary(parms* p, double**** storage, double **fx, double **fy, double **fz, int half);
void bound_v_correction(parms* p, double**** storage, double**);
void tmp_fun(parms* p, double** fx, double** fy, double** fz, double** fx_new, double** fy_new, double** fz_new);
void itereate_to_correct(int, parms*, double**, double**, double**, double****, double*, double*, double*, double*, int, int, double, int, double*,
		int, double**);
void set_neumann_conditions(parms*, int, double*, double*, double*, double*, int, int);
