#include "functionality.h"

double** allocate_dirichlet(int rows, int cols) {

	double **f;

	f = (double**) malloc(rows * sizeof(double*));
	if (f == NULL) {
		printf("allocation failed for first dimension of f\n");
	}
	for (int i = 0; i < rows; i++) {
		f[i] = (double*) malloc(cols * sizeof(double));
		if (f[i] == NULL) {
			printf("allocation failed for second dimension of f\n");
		}
	}
	return (f);
}

double* allocate_neuman(int indx) {
	double *buffer = (double*) malloc(indx * sizeof(double));
	return (buffer);
}

void initialize_f(parms* p, double **fx, double **fy, double **fz) {
	for (int i = 0; i < p->nx; i++) {
		for (int j = 0; j < p->ny; j++) {
			fx[i][j] = 0.0;
			fy[i][j] = 0.0;
			fz[i][j] = 0.0;
		}
	}
}

void initialize_bd(int indx, double *bd) {
	for (int i = 0; i < indx; i++) {
		bd[i] = 0.0;
	}
}

void print_stdout(parms p, double **fx, double **fy, double **fz) {
	printf("C code output\nx coordinates\n");
	for (int i = 0; i < p.nx; i++) {
		for (int j = 0; j < p.ny; j++) {
			printf("%2.12lf\t", fx[i][j] * 1e-3);
		}
		printf("\n");
	}
	printf("y coordinates\n");
	for (int i = 0; i < p.nx; i++) {
		for (int j = 0; j < p.ny; j++) {
			printf("%2.12lf\t", fy[i][j] * 1e-3);
		}
		printf("\n");
	}
	printf("z coordinates\n");
	for (int i = 0; i < p.nx; i++) {
		for (int j = 0; j < p.ny; j++) {
			printf("%2.12lf\t", fz[i][j] * 1e-3);
		}
		printf("\n");
	}
}

void format_vtk(parms p, double** fx, double** fy, double** fz, char* prefix) {
	FILE *f;
	char filename[50];
	int err = sprintf(filename, "%s.vtk", prefix);
	f = fopen(filename, "w+");

	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "Default Set with interval=1s\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "\n");

	fprintf(f, "DATASET STRUCTURED_GRID\n");
	int d_x = p.nx, d_y = p.ny, d_z = 1;
	fprintf(f, "DIMENSIONS %d %d %d\n", d_x, d_y, d_z);
	fprintf(f, "POINTS %d double\n", d_x * d_y * d_z);

	double ***a;
	a = (double***) malloc(3 * sizeof(double**));
	for (int i = 0; i < 3; i++)
		a[i] = (double**) malloc((p.nx) * sizeof(double*));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < p.nx; j++) {
			a[i][j] = (double*) malloc((p.ny) * sizeof(double));
		}
	}
	for (int j = 0; j < p.ny; j++) {
		for (int i = 0; i < p.nx; i++) {
			fprintf(f, "%2.12lf %2.12lf %2.12lf\n", fx[i][j] * 1e-3, fy[i][j] * 1e-3, fz[i][j] * 1e-3);
			a[0][i][j] = fx[i][j] * 1e-3;
			a[1][i][j] = fy[i][j] * 1e-3;
			a[2][i][j] = fz[i][j] * 1e-3;
		}
	}
//	format_stl(p, a, prefix);
	fclose(f);
}

void format_stl(parms p, double*** f, char* prefix, int Total_pannels_axially, int Total_pannels_circumferentially) {
	FILE *FileWrite;
	char filename[50];
	int offset = p.n + 2;
	int err = sprintf(filename, "%s.stl", prefix);
	FileWrite = fopen(filename, "w+");
	int circ = 2 * p.ny;
	fprintf(FileWrite, "solid ascii\n");
	for (int j = 0; j < Total_pannels_circumferentially; j++) {
		for (int i = 0; i < Total_pannels_axially + 1; i++) {

			if ((j >= 0) && (j < Total_pannels_circumferentially)) {
				if (i == 0) {
					facet1(i, j, f, FileWrite);
				} else if (i == (Total_pannels_axially)) {
					facet2(i, j, f, FileWrite);
				} else {
					facet1(i, j, f, FileWrite);
					facet2(i, j, f, FileWrite);
				}
			} else if (j == (Total_pannels_circumferentially - 1)) {
				//do nothing
			}
		}
	}
	fprintf(FileWrite, "endsolid");
	fclose(FileWrite);
}
void facet1(int i, int j, double ***f, FILE* FileWrite) {

	fprintf(FileWrite, "facet normal 0 0 0\n");
	fprintf(FileWrite, "outer loop\n");
	fprintf(FileWrite, "vertex %2.8lf %2.8lf %2.8lf\n", f[0][i][j], f[1][i][j], f[2][i][j]);
	fprintf(FileWrite, "vertex %2.8lf %2.8lf %2.8lf\n", f[0][i][j + 1], f[1][i][j + 1], f[2][i][j + 1]);
	fprintf(FileWrite, "vertex %2.8lf %2.8lf %2.8lf\n", f[0][i + 1][j], f[1][i + 1][j], f[2][i + 1][j]);
	fprintf(FileWrite, "endloop\n");
	fprintf(FileWrite, "endfacet\n");
}
void facet2(int i, int j, double ***f, FILE* FileWrite) {

	fprintf(FileWrite, "facet normal 0 0 0\n");
	fprintf(FileWrite, "outer loop\n");
	fprintf(FileWrite, "vertex %2.8lf %2.8lf %2.8lf\n", f[0][i][j], f[1][i][j], f[2][i][j]);
	fprintf(FileWrite, "vertex %2.8lf %2.8lf %2.8lf\n", f[0][i][j + 1], f[1][i][j + 1], f[2][i][j + 1]);
	fprintf(FileWrite, "vertex %2.8lf %2.8lf %2.8lf\n", f[0][i - 1][j + 1], f[1][i - 1][j + 1], f[2][i - 1][j + 1]);
	fprintf(FileWrite, "endloop\n");
	fprintf(FileWrite, "endfacet\n");
}

double** solve(parms* p, double *bda, double *bdb, double *bdc, double *bdd, double **f, int idf, int iflag, double tol, int itcg, double *w,
		int lw) {

	double t[p->ny][p->nx];			//Transpose matrix
	double a[p->nx * p->ny];		//a vector to pass the data to dbihar function.
	int offset = p->nx;
	for (int i = 0; i < p->nx; i++) {
		for (int j = 0; j < p->ny; j++) {
			t[j][i] = f[i][j];
		}
	}
	for (int i = 0; i < p->ny; i++) {
		for (int j = 0; j < p->nx; j++) {
			a[i * offset + j] = t[i][j];
		}
	}

	dbihar(&p->a, &p->b, &p->m, bda, bdb, bdc, bdd, &p->c, &p->d, &p->n, a, &idf, &p->alpha, &p->beta, &iflag, &tol, &itcg, w, &lw);

	for (int i = 0; i < p->ny; i++) {
		for (int j = 0; j < p->nx; j++) {
			t[i][j] = a[i * offset + j];
		}
	}
	for (int i = 0; i < p->nx; i++) {
		for (int j = 0; j < p->ny; j++) {
			f[i][j] = t[j][i];
		}
	}
	return (f);
}

void dirichlet_boundary_parent_segment(parms p, double **fx, double **fy, double **fz) {
	for (int i = 1; i <= p.nx; i++) {
		double x = p.a + (i - 1) * (p.b - p.a) / (p.m + 1);
		double r = p.neck - (p.neck - p.rb) * (x);
		int j;
		///v=c
		j = 1;
		fx[i - 1][j - 1] = -r * cos(p.c);
		fy[i - 1][j - 1] = p.h2 - (p.h2 - p.h1) * x;
		fz[i - 1][j - 1] = r * sin(p.c);
		///v=d
		j = p.ny;
		fx[i - 1][j - 1] = -r * cos(p.d);
		fy[i - 1][j - 1] = p.h2 - (p.h2 - p.h1) * x;
		fz[i - 1][j - 1] = r * sin(p.d);
	}

	for (int j = 1; j <= p.ny; j++) {
		double y = p.c + (j - 1) * (p.d - p.c) / (p.n + 1);
		int i;
		///x=a
		i = 1;
		fx[i - 1][j - 1] = -p.neck * cos(y);
		fy[i - 1][j - 1] = p.h2;
		fz[i - 1][j - 1] = p.neck * sin(y);
		///x=b
		i = p.nx;
		fx[i - 1][j - 1] = -p.rb * cos(y);
		fy[i - 1][j - 1] = p.h1;
		fz[i - 1][j - 1] = p.rb * sin(y);
	}
}
double* dirichlet_boundary_outer_wall(parms p, double **fx, double **fy, double **fz) {
	for (int i = 1; i <= p.nx; i++) {
		double x = p.a + (i - 1) * (p.b - p.a) / (p.m + 1);
		int j;
		double h = p.h3 - (p.h3 - p.h2) * x;
		double r = p.neck - (p.neck - p.rt) * (1 - x);
		///v=c
		j = 1;
		fx[i - 1][j - 1] = -(h * cos(p.alfa / 2) + /*p.ra*/r * sin(p.alfa / 2) * cos(p.c));
		fy[i - 1][j - 1] = h * sin(p.alfa / 2) - /*p.ra*/r * cos(p.alfa / 2) * cos(p.c);
		fz[i - 1][j - 1] = /*p.ra*/r * sin(p.c);
		///v=d
		j = p.ny;
		fx[i - 1][j - 1] = -(h * cos(p.alfa / 2) + /*p.ra*/r * sin(p.alfa / 2) * cos(p.d));
		fy[i - 1][j - 1] = h * sin(p.alfa / 2) - /*p.ra*/r * cos(p.alfa / 2) * cos(p.d);
		fz[i - 1][j - 1] = /*p.ra*/r * sin(p.d);
	}
	double* theta_val = (double*) malloc(p.ny * sizeof(double));
	for (int j = 1; j <= p.ny; j++) {
		double y = p.c + (j - 1) * (p.d - p.c) / (p.n + 1);
		theta_val[j - 1] = y;
		int i;
		///x=a
		i = 1;
		fx[i - 1][j - 1] = -(p.h3 * cos(p.alfa / 2) + p.ra * sin(p.alfa / 2) * cos(y));
		fy[i - 1][j - 1] = p.h3 * sin(p.alfa / 2) - p.ra * cos(p.alfa / 2) * cos(y);
		fz[i - 1][j - 1] = p.ra * sin(y);

		///x=b
		i = p.nx;
		fx[i - 1][j - 1] = -p.neck * cos(y);
		fy[i - 1][j - 1] = p.h2;
		fz[i - 1][j - 1] = p.neck * sin(y);
	}
	return (theta_val);
}

double* dirichlet_boundary_inner_wall(parms p, double **fx, double **fy, double **fz) {
	for (int i = 1; i <= p.nx; i++) {
		double x = p.a + (i - 1) * (p.b - p.a) / (p.m + 1);
		int j;
		///v=c
		double h = p.h3 - (p.h3 - p.h2) * x;
		double r = p.neck - (p.neck - p.rt) * (1 - x);
		j = 1;
		fx[i - 1][j - 1] = -(h * cos(p.alfa / 2) + /*p.ra*/r * sin(p.alfa / 2) * cos(p.c));
		fy[i - 1][j - 1] = h * sin(p.alfa / 2) - /*p.ra*/r * cos(p.alfa / 2) * cos(p.c);
		fz[i - 1][j - 1] = /*p.ra*/r * sin(p.c);
		///v=d
		j = p.ny;
		fx[i - 1][j - 1] = -(h * cos(p.alfa / 2) + /*p.ra*/r * sin(p.alfa / 2) * cos(p.d));
		fy[i - 1][j - 1] = h * sin(p.alfa / 2) - /*p.ra*/r * cos(p.alfa / 2) * cos(p.d);
		fz[i - 1][j - 1] = /*p.ra*/r * sin(p.d);
	}
	double* theta_val = (double*) malloc(p.ny * sizeof(double));
	for (int j = 1; j <= p.ny; j++) {
		double y = p.c + (j - 1) * (p.d - p.c) / (p.n + 1);
		double direction = 0.0;
		theta_val[j - 1] = y;
		int i;
		///x=a
		i = 1;
		fx[i - 1][j - 1] = -(p.h3 * cos(p.alfa / 2) + p.ra * sin(p.alfa / 2) * cos(y));
		fy[i - 1][j - 1] = p.h3 * sin(p.alfa / 2) - p.ra * cos(p.alfa / 2) * cos(y);
		fz[i - 1][j - 1] = p.ra * sin(y);
		///x=b
		i = p.nx;

		fx[i - 1][j - 1] = 0.0;
		if ((y >= pi / 2) && (y <= 3 * pi / 2)) {
			fy[i - 1][j - 1] = -p.apex * cos(y);
		} else if ((y >= 3 * pi / 2) && (y <= 2 * pi + pi / 2)) {
			fy[i - 1][j - 1] = p.apex * cos(y);
		}

		fz[i - 1][j - 1] = p.neck * sin(y);
	}
	return (theta_val);
}

void dirichlet_boundary_end_cap_parent(parms p, double **fx, double **fy, double **fz) {

	double x, y;
	double* theta_val = (double*) malloc(p.ny * sizeof(double));
	for (int j = 1; j <= p.ny; j++) {
		y = p.c + (j - 1) * (p.d - p.c) / (p.n + 1);
		theta_val[j - 1] = y;
		for (int i = 1; i <= p.nx; i++) {
			x = p.a + (i - 1) * (p.b - p.a) / (p.m + 1);
			fx[i - 1][j - 1] = p.rb * (1 - x) * cos(y);
			fy[i - 1][j - 1] = p.h1;
			fz[i - 1][j - 1] = p.rb * (1 - x) * sin(y);
		}
	}
}

void dirichlet_boundary_end_cap_daughter(parms p, double **fx, double **fy, double **fz) {

	double x, y;
	double* theta_val = (double*) malloc(p.ny * sizeof(double));
	for (int j = 1; j <= p.ny; j++) {
		y = p.c + (j - 1) * (p.d - p.c) / (p.n + 1);
		theta_val[j - 1] = y;
		for (int i = 1; i <= p.nx; i++) {
			x = p.a + (i - 1) * (p.b - p.a) / (p.m + 1);
			fx[i - 1][j - 1] = -(p.h3 * cos(p.alfa / 2) + p.rt * (1 - x) * sin(p.alfa / 2) * cos(y));
			fy[i - 1][j - 1] = p.h3 * sin(p.alfa / 2) - p.rt * (1 - x) * cos(p.alfa / 2) * cos(y);
			fz[i - 1][j - 1] = p.rt * (1 - x) * sin(y);
		}
	}
}

///The following functions are used to create cellular mesh inside every quadilateral 
///constituting the subdomain. The decisions consider an EC of dimension 65um x 10um and a SMC is 
///50um x 5um.
/****************************************************************/
int ec_mesh(parms p, double** fx, double** fy, double** fz, char *prefix, char* suffix, int downstream_offset, int upstream_offset) {
	/****************************************************************/
	int indx = 0;
	int total_points = 0;
	double hx = 65e-6, hy = 50e-6;
	FILE *fw, *f_tmp;
	char tmp_file[50], master_file[50], centeroid_file[50], output_file[50];

	int err;
	err = sprintf(master_file, "%s_master_%s.txt", prefix, suffix);
	err = sprintf(tmp_file, "%s_mesh_tmp_%s.txt", prefix, suffix);
	err = sprintf(centeroid_file, "%s_centeroids_%s", prefix, suffix);
	err = sprintf(output_file, "%s_mesh_%s.vtk", prefix, suffix);

	int buffer[2];
	fw = fopen(master_file, "w");
	if (fw == NULL) {
		printf("Unable to open ec_master.txt for writing\n");
	}
	f_tmp = fopen(tmp_file, "w");
	if (f_tmp == NULL) {
		printf("Unable to open ec_master.txt for writing\n");
	}
	fseek(fw, 0, SEEK_SET);
	for (int i = 0 + downstream_offset; i < (p.nx - 1 - upstream_offset); i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			double ax = fx[i][j] * 1e-3, bx = fx[i + 1][j] * 1e-3, cx = fx[i][j + 1] * 1e-3, dx = fx[i + 1][j + 1] * 1e-3, ay = fy[i][j] * 1e-3, by =
					fy[i + 1][j] * 1e-3, cy = fy[i][j + 1] * 1e-3, dy = fy[i + 1][j + 1] * 1e-3, az = fz[i][j] * 1e-3, bz = fz[i + 1][j] * 1e-3, cz =
					fz[i][j + 1] * 1e-3, dz = fz[i + 1][j + 1] * 1e-3;
			int m, n;
			int num_ec_axially, num_smc_circumferentially;
			double tmp = P2(bx - ax) + P2(by - ay) + P2(bz - az);
			num_ec_axially = 8;		//rint(sqrt(tmp) / hx);
			tmp = P2(cx - ax) + P2(cy - ay) + P2(cz - az);
			num_smc_circumferentially = 5;		//rint(sqrt(tmp) / hy);

			double **x, **y, **z;

			/// Meshing ECs pannels
			m = num_ec_axially + 1;
			n = (num_smc_circumferentially * 5) + 1;
			printf("[%d,%d] num_ec = %d \t num_smc = %d\tm=%d\tn=%d\t%2.7lf\t%2.7lf\t%2.7lf\t%2.7lf\n", i, j, num_ec_axially,
					num_smc_circumferentially, m, n, ax, bx, cx, dx);
			total_points = total_points + (m * n);
			buffer[0] = m;
			buffer[1] = n;
			fwrite(buffer, sizeof(int), 2, fw);
			x = allocate_dirichlet(m, n);
			y = allocate_dirichlet(m, n);
			z = allocate_dirichlet(m, n);

			for (int i = 0; i < m; i++) {
				int j = 0;
				x[i][j] = ax + i * (bx - ax) / (m - 1);
				y[i][j] = ay + i * (by - ay) / (m - 1);
				z[i][j] = az + i * (bz - az) / (m - 1);
				j = n - 1;
				x[i][j] = cx + i * (dx - cx) / (m - 1);
				y[i][j] = cy + i * (dy - cy) / (m - 1);
				z[i][j] = cz + i * (dz - cz) / (m - 1);
			}

			for (int j = 0; j < n; j++) {
				int i = 0;
				x[i][j] = ax + j * (cx - ax) / (n - 1);
				y[i][j] = ay + j * (cy - ay) / (n - 1);
				z[i][j] = az + j * (cz - az) / (n - 1);
				i = m - 1;
				x[i][j] = bx + j * (dx - bx) / (n - 1);
				y[i][j] = by + j * (dy - by) / (n - 1);
				z[i][j] = bz + j * (dz - bz) / (n - 1);
			}

			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					x[i][j] = x[i][0] + j * (x[i][n - 1] - x[i][0]) / (n - 1);
					y[i][j] = y[i][0] + j * (y[i][n - 1] - y[i][0]) / (n - 1);
					z[i][j] = z[i][0] + j * (z[i][n - 1] - z[i][0]) / (n - 1);
				}
			}
			for (int jj = 0; jj < n; jj++) {
				for (int ii = 0; ii < m; ii++) {
					double tmp_buf[3];
					tmp_buf[0] = x[ii][jj];
					tmp_buf[1] = y[ii][jj];
					tmp_buf[2] = z[ii][jj];
					fwrite(tmp_buf, sizeof(double), 3, f_tmp);
				}
			}
			indx++;
		}
	}
	fclose(fw);
	fclose(f_tmp);
	int number_of_ECs = format_vtk_small_mesh(tmp_file, master_file, centeroid_file, output_file, p.nx, p.ny, total_points, downstream_offset,
			upstream_offset);
	int status;
	status = remove(master_file);
	status = remove(tmp_file);
	return (number_of_ECs);
}
/*****************************************************************/
int smc_mesh(parms p, double**** storage, char* prefix, char* suffix, int downstream_offset, int upstream_offset) {
	/*****************************************************************/
	int indx = 0;
	int total_points = 0;
	double hx = 65e-6, hy = 50e-6;
	FILE *fw, *f_tmp;
	int buffer[2];
	char tmp_file[50], master_file[50], centeroid_file[50], output_file[50];

	int err;
	err = sprintf(master_file, "%s_master_%s.txt", prefix, suffix);
	err = sprintf(tmp_file, "%s_mesh_tmp_%s.txt", prefix, suffix);
	err = sprintf(centeroid_file, "%s_centeroids_%s", prefix, suffix);
	err = sprintf(output_file, "%s_mesh_%s.vtk", prefix, suffix);
	fw = fopen(master_file, "w");
	if (fw == NULL) {
		printf("Unable to open smc_master.txt for writing\n");
	}
	f_tmp = fopen(tmp_file, "w");
	if (f_tmp == NULL) {
		printf("Unable to open smc_master.txt for writing\n");
	}
	fseek(fw, 0, SEEK_SET);
	for (int i = 0 + downstream_offset; i < (p.nx - 1 - upstream_offset); i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			double ax = storage[RIGHT][x_coord][i][j] * 1e-3, bx = storage[RIGHT][x_coord][i + 1][j] * 1e-3, cx = storage[RIGHT][x_coord][i][j + 1]
					* 1e-3, dx = storage[RIGHT][x_coord][i + 1][j + 1] * 1e-3, ay = storage[RIGHT][y_coord][i][j] * 1e-3, by =
					storage[RIGHT][y_coord][i + 1][j] * 1e-3, cy = storage[RIGHT][y_coord][i][j + 1] * 1e-3, dy =
					storage[RIGHT][y_coord][i + 1][j + 1] * 1e-3, az = storage[RIGHT][z_coord][i][j] * 1e-3, bz = storage[RIGHT][z_coord][i + 1][j]
					* 1e-3, cz = storage[RIGHT][z_coord][i][j + 1] * 1e-3, dz = storage[RIGHT][z_coord][i + 1][j + 1] * 1e-3;
			int m, n;
			int num_ec_axially, num_smc_circumferentially;
			double tmp = P2(bx-ax) + P2(by-ay) + P2(bz-az);
			num_ec_axially = 8;		//rint(sqrt(tmp) / hx);
			tmp = P2(cx-ax) + P2(cy-ay) + P2(cz-az);
			num_smc_circumferentially = 5;		//rint(sqrt(tmp) / hy);

			double **x, **y, **z;

			/// Meshing SMCs pannels
			m = (num_ec_axially * 13) + 1;
			n = num_smc_circumferentially + 1;

			total_points = total_points + (m * n);
			buffer[0] = m;
			buffer[1] = n;
			fwrite(buffer, sizeof(int), 2, fw);

			x = allocate_dirichlet(m, n);
			y = allocate_dirichlet(m, n);
			z = allocate_dirichlet(m, n);
			for (int i = 0; i < m; i++) {
				int j = 0;
				x[i][j] = ax + i * (bx - ax) / (m - 1);
				y[i][j] = ay + i * (by - ay) / (m - 1);
				z[i][j] = az + i * (bz - az) / (m - 1);
				j = n - 1;
				x[i][j] = cx + i * (dx - cx) / (m - 1);
				y[i][j] = cy + i * (dy - cy) / (m - 1);
				z[i][j] = cz + i * (dz - cz) / (m - 1);
			}
			for (int j = 0; j < n; j++) {
				int i = 0;
				x[i][j] = ax + j * (cx - ax) / (n - 1);
				y[i][j] = ay + j * (cy - ay) / (n - 1);
				z[i][j] = az + j * (cz - az) / (n - 1);
				i = m - 1;
				x[i][j] = bx + j * (dx - bx) / (n - 1);
				y[i][j] = by + j * (dy - by) / (n - 1);
				z[i][j] = bz + j * (dz - bz) / (n - 1);
			}

			for (int i = 1; i < m - 1; i++) {
				for (int j = 1; j < n - 1; j++) {
					x[i][j] = x[i][0] + j * (x[i][n - 1] - x[i][0]) / (n - 1);
					y[i][j] = y[i][0] + j * (y[i][n - 1] - y[i][0]) / (n - 1);
					z[i][j] = z[i][0] + j * (z[i][n - 1] - z[i][0]) / (n - 1);
				}
			}

			for (int jj = 0; jj < n; jj++) {
				for (int ii = 0; ii < m; ii++) {
					double tmp_buf[3];
					tmp_buf[0] = x[ii][jj];
					tmp_buf[1] = y[ii][jj];
					tmp_buf[2] = z[ii][jj];
					fwrite(tmp_buf, sizeof(double), 3, f_tmp);
				}
			}
			indx++;
		}

		for (int j = 0; j < p.ny - 1; j++) {
			double ax = storage[LEFT][x_coord][i][j] * 1e-3, bx = storage[LEFT][x_coord][i + 1][j] * 1e-3, cx = storage[LEFT][x_coord][i][j + 1]
					* 1e-3, dx = storage[LEFT][x_coord][i + 1][j + 1] * 1e-3, ay = storage[LEFT][y_coord][i][j] * 1e-3, by = storage[LEFT][y_coord][i
					+ 1][j] * 1e-3, cy = storage[LEFT][y_coord][i][j + 1] * 1e-3, dy = storage[LEFT][y_coord][i + 1][j + 1] * 1e-3, az =
					storage[LEFT][z_coord][i][j] * 1e-3, bz = storage[LEFT][z_coord][i + 1][j] * 1e-3, cz = storage[LEFT][z_coord][i][j + 1] * 1e-3,
					dz = storage[LEFT][z_coord][i + 1][j + 1] * 1e-3;
			int m, n;
			int num_ec_axially, num_smc_circumferentially;
			double tmp = P2(bx-ax) + P2(by-ay) + P2(bz-az);
			num_ec_axially = 8;		//rint(sqrt(tmp) / hx);
			tmp = P2(cx-ax) + P2(cy-ay) + P2(cz-az);
			num_smc_circumferentially = 5;		//rint(sqrt(tmp) / hy);

			double **x, **y, **z;

			/// Meshing SMCs pannels
			m = (num_ec_axially * 13) + 1;
			n = num_smc_circumferentially + 1;

			total_points = total_points + (m * n);
			buffer[0] = m;
			buffer[1] = n;
			fwrite(buffer, sizeof(int), 2, fw);

			x = allocate_dirichlet(m, n);
			y = allocate_dirichlet(m, n);
			z = allocate_dirichlet(m, n);
			for (int i = 0; i < m; i++) {
				int j = 0;
				x[i][j] = ax + i * (bx - ax) / (m - 1);
				y[i][j] = ay + i * (by - ay) / (m - 1);
				z[i][j] = az + i * (bz - az) / (m - 1);
				j = n - 1;
				x[i][j] = cx + i * (dx - cx) / (m - 1);
				y[i][j] = cy + i * (dy - cy) / (m - 1);
				z[i][j] = cz + i * (dz - cz) / (m - 1);
			}
			for (int j = 0; j < n; j++) {
				int i = 0;
				x[i][j] = ax + j * (cx - ax) / (n - 1);
				y[i][j] = ay + j * (cy - ay) / (n - 1);
				z[i][j] = az + j * (cz - az) / (n - 1);
				i = m - 1;
				x[i][j] = bx + j * (dx - bx) / (n - 1);
				y[i][j] = by + j * (dy - by) / (n - 1);
				z[i][j] = bz + j * (dz - bz) / (n - 1);
			}

			for (int i = 1; i < m - 1; i++) {
				for (int j = 1; j < n - 1; j++) {
					x[i][j] = x[i][0] + j * (x[i][n - 1] - x[i][0]) / (n - 1);
					y[i][j] = y[i][0] + j * (y[i][n - 1] - y[i][0]) / (n - 1);
					z[i][j] = z[i][0] + j * (z[i][n - 1] - z[i][0]) / (n - 1);
				}
			}

			for (int jj = 0; jj < n; jj++) {
				for (int ii = 0; ii < m; ii++) {
					double tmp_buf[3];
					tmp_buf[0] = x[ii][jj];
					tmp_buf[1] = y[ii][jj];
					tmp_buf[2] = z[ii][jj];
					fwrite(tmp_buf, sizeof(double), 3, f_tmp);
				}
			}
			indx++;
		}
	}		//end of loop on fx, fy, fz
	fclose(fw);
	fclose(f_tmp);
	int number_of_SMCs = format_vtk_small_mesh(tmp_file, master_file, centeroid_file, output_file, p.nx, p.ny, total_points, downstream_offset,
			upstream_offset);

	int status;
//	status = remove(master_file);
//	status = remove(tmp_file);

	return (number_of_SMCs);
}

/**********************************************************/
void print_func(int nx, int ny, double **f, char* str) {
	/**********************************************************/
	printf("%s coordinate\n", str);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			printf("%2.12lf\t", f[i][j]);
		}
		printf("\n");
	}
}
/**********************************************************************************/
///This function takes the data formualted by the calling function in each quadilateral 
///of the subdomain and stores it into a format such that it is accessible later while
///visualizing.
void store_data(int m, int n, double **x, double **y, double **z, char *prefix, int indx) {

	/**********************************************************************************/
	FILE *f;
	char filename[50];
	int err = sprintf(filename, "%s.vtk", prefix);
	f = fopen(filename, "a+");

	fseek(f, 0, SEEK_END);
	long int size = ftell(f);
	printf("%ld\n", size / sizeof(double));
	double buffer[3];
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			buffer[0] = x[i][j];
			buffer[1] = y[i][j];
			buffer[2] = z[i][j];
			fwrite(buffer, sizeof(double), 3, f);
		}
	}
	fclose(f);
}
/******************************************************************/
///This function produces vtk file visualizing the cellular mesh as
///vtk unstructured gird followed by the definition of the connectivity
///of points into quadilateral cells.
int format_vtk_small_mesh(char *tmp_file, char *master_file, char *centeroids, char *output_file, int nx, int ny, int total_points,
		int downstream_offset, int upstream_offset) {
	/******************************************************************/
	FILE *fr, *fw, *fm, *cell_ref;
	char centeroid_file[50];
	int err = sprintf(centeroid_file, "%s.txt", centeroids);
	fr = fopen(tmp_file, "r");
	if (fr == NULL)
		printf("failed opening smc_mesh_tmp for read operation.txt\n");
	fw = fopen(output_file, "w+");
	if (fw == NULL)
		printf("failed opening smc_mesh.txt for write operation.\n");
	fm = fopen(master_file, "r+");
	if (fm == NULL)
		printf("failed opening smc_master.txt\n");
	cell_ref = fopen(centeroid_file, "w");
	if (cell_ref == NULL)
		printf("failed opening %s\n", centeroid_file);

	printf("step 1\n");

	fprintf(fw, "# vtk DataFile Version 3.0\n");
	fprintf(fw, "smc mesh\n");
	fprintf(fw, "ASCII\n");
	fprintf(fw, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fw, "POINTS %d double\n", total_points);
	double buffer[3];
	for (int i = 0; i < total_points; i++) {
		fread(buffer, sizeof(double), 3, fr);
		fprintf(fw, "%2.10lf %2.10lf %2.10lf\n", buffer[0], buffer[1], buffer[2]);
	}
	int buf[nx - 1][ny - 1][2];
	int m, n;
	int count = 0;
	for (int i = 0; i < (nx - upstream_offset - downstream_offset - 1); i++) {
		for (int j = 0; j < ny - 1; j++) {
			int tmp[2];
			fread(&buf[i][j][0], sizeof(int), 2, fm);
			count += (buf[i][j][0] - 1) * (buf[i][j][1] - 1);
		}
	}
	fprintf(fw, "CELLS %d %d\n", count, count * 5);
	int circ_offset = 0, axial_offset = 0, offset = 0, counter = 0;
	for (int p_indx = 0 + downstream_offset; p_indx < (nx - 1 - upstream_offset); p_indx++) {
		for (int q = 0; q < ny - 1; q++) {
			int p = p_indx - downstream_offset;
			m = buf[p][q][0];
			n = buf[p][q][1];

			if (p + q > 0) {
				if (p == 0) {
					offset += (buf[p][q - 1][0]) * (buf[p][q - 1][1]);
				} else if (p > 0) {
					if (q == 0) {
						offset += (buf[p - 1][ny - 2][0]) * (buf[p - 1][ny - 2][1]);
					} else {
						offset += (buf[p][q - 1][0]) * (buf[p][q - 1][1]);
					}
				}
			}

			for (int i = 0; i < m - 1; i++) {
				for (int j = 0; j < n - 1; j++) {

					fseek(fr, 0, SEEK_SET);

					fprintf(fw, "4 %d %d %d %d\n", (offset + i + j * (m)), (offset + i + (j + 1) * (m)), (offset + i + 1 + (j + 1) * (m)),
							(offset + i + 1 + j * (m)));

					double buffer[4][3];
					long int disp;
					disp = (3 * (offset + i + j * m) * sizeof(double));
					fseek(fr, disp, SEEK_SET);
					fread(&buffer[0][0], sizeof(double), 3, fr);
					disp = (3 * (offset + i + (j + 1) * m) * sizeof(double));
					fseek(fr, disp, SEEK_SET);
					fread(&buffer[1][0], sizeof(double), 3, fr);
					disp = (3 * (offset + i + 1 + (j + 1) * m) * sizeof(double));
					fseek(fr, disp, SEEK_SET);
					fread(&buffer[2][0], sizeof(double), 3, fr);
					disp = (3 * (offset + i + 1 + j * m) * sizeof(double));
					fseek(fr, disp, SEEK_SET);
					fread(&buffer[3][0], sizeof(double), 3, fr);
					double coords[3] = { 0.0, 0.0, 0.0 };
					for (int k = 0; k < 4; k++) {
						coords[0] += buffer[k][0];
						coords[1] += buffer[k][1];
						coords[2] += buffer[k][2];
					}
					coords[0] = coords[0] / 4;
					coords[1] = coords[1] / 4;
					coords[2] = coords[2] / 4;
					fwrite(coords, sizeof(double), 3, cell_ref);

				}
			}
		}
	}

	printf("total_points=%d\n", total_points);
	fprintf(fw, "CELL_TYPES %d\n", count);
	for (int i = 0; i < count; i++) {
		fprintf(fw, "9\n");
	}
	fprintf(fw, "\n");
	fclose(fr);
	fclose(fw);
	fclose(fm);
	fclose(cell_ref);
	format_vtk_centeroid(centeroids, count);
	int status = remove(centeroid_file);
	return (count);
}

/********************************************************/
void format_vtk_centeroid(char* filename, int count) {
	/********************************************************/
	FILE *fr, *fw;
	char input_file[50], output_file[50];

	int err = sprintf(input_file, "%s.txt", filename);
	err = sprintf(output_file, "%s.vtk", filename);
	printf("%s %s\n", input_file, output_file);
	fr = fopen(input_file, "r");
	if (fr == NULL) {
		printf("Unable to open %s.txt", input_file);
	}
	fw = fopen(output_file, "w");
	if (fw == NULL) {
		printf("Unable to open %s.txt", output_file);
	}
	fprintf(fw, "# vtk DataFile Version 3.0\n");
	fprintf(fw, "%s\n", filename);
	fprintf(fw, "ASCII\n");
	fprintf(fw, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fw, "POINTS %d double\n", count);
	fseek(fr, 0, SEEK_SET);
	for (int i = 0; i < count; i++) {
		double buffer[3];
		fread(buffer, sizeof(double), 3, fr);
		fprintf(fw, "%lf\t%lf\t%lf\n", buffer[0], buffer[1], buffer[2]);
	}
	fprintf(fw, "CELLS %d %d\n", count, 2 * count);
	for (int i = 0; i < count; i++) {
		fprintf(fw, "1 %d\n", i);
	}
	fprintf(fw, "CELL_TYPES %d\n", count);
	for (int i = 0; i < count; i++) {
		fprintf(fw, "%d\n", 1);
	}
	fclose(fr);
	fclose(fw);

}

void deallocate(int m, int n, double **f) {
	for (int i = 0; i < m; i++) {
		free(f[i]);
	}
	free(f);
}

///Checking if the file is empty
bool isEmpty(FILE *file) {
	long savedOffset = ftell(file);
	fseek(file, 0, SEEK_END);
	if (ftell(file) == 0) {
		return true;
	}
	fseek(file, savedOffset, SEEK_SET);
	return false;
}

/// Rudemnetary functions. Not being used but may be useful in future, if need be.
void record_data(int k, parms p, double **fx, double **fy, double **fz, double ****f) {
	for (int i = 0; i < p.nx; i++) {
		for (int j = 0; j < p.ny; j++) {
			f[k][0][i][j] = fx[i][j];
			f[k][1][i][j] = fy[i][j];
			f[k][2][i][j] = fz[i][j];
		}
	}
}

void format_vtk_unstructuredGrid(parms p, double** fx, double** fy, double** fz, char* prefix, int downstream_offset, int upstream_offset) {

	FILE *f;
	char filename[50];
	int err = sprintf(filename, "%s.vtk", prefix);
	f = fopen(filename, "w+");

	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "Default Set with interval=1s\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "\n");

	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
	int TotalNumberOfPoints = ((p.nx - upstream_offset) - downstream_offset) * p.ny, TotalNumberOfCells = (((p.nx - upstream_offset)
			- downstream_offset) - 1) * (p.ny - 1);
	fprintf(f, "POINTS %d double\n", TotalNumberOfPoints);

	double ***a;
	a = (double***) malloc(3 * sizeof(double**));
	for (int i = 0; i < 3; i++)
		a[i] = (double**) malloc((p.nx) * sizeof(double*));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < p.nx; j++) {
			a[i][j] = (double*) malloc((p.ny) * sizeof(double));
		}
	}
	for (int i = 0 + downstream_offset; i < (p.nx - upstream_offset); i++) {
		for (int j = 0; j < p.ny; j++) {

			fprintf(f, "%2.12lf %2.12lf %2.12lf\n", fx[i][j] * 1e-3, fy[i][j] * 1e-3, fz[i][j] * 1e-3);
			a[0][i][j] = fx[i][j] * 1e-3;
			a[1][i][j] = fy[i][j] * 1e-3;
			a[2][i][j] = fz[i][j] * 1e-3;
		}
	}

	fprintf(f, "CELLS %d %d\n", TotalNumberOfCells, 5 * TotalNumberOfCells);
	for (int i = 0 + downstream_offset; i < (p.nx - 1 - upstream_offset); i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			fprintf(f, "%d %d %d %d %d\n", 4, ((i - downstream_offset) * p.ny) + j, (((i - downstream_offset) + 1) * p.ny) + j,
					(((i - downstream_offset) + 1) * p.ny) + j + 1, ((i - downstream_offset) * p.ny) + j + 1);
		}
	}
	fprintf(f, "CELL_TYPES %d\n", TotalNumberOfCells);
	for (int i = 0 + downstream_offset; i < (p.nx - 1 - upstream_offset); i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			fprintf(f, "%d\n", 9);
		}
	}

	//format_stl(p, a, prefix);
	fclose(f);
}

double**** allocate_storage_array(int nx, int ny) {
	double****storage = (double****) malloc(2 * sizeof(double***));
	for (int i = 0; i < 2; i++)
		storage[i] = (double***) malloc(3 * sizeof(double**));
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			storage[i][j] = (double**) malloc(nx * sizeof(double*));
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < nx; k++)
				storage[i][j][k] = (double*) malloc(ny * sizeof(double));
		}
	}
	return (storage);
}

void store_arrays(parms* p, double** fx, double** fy, double** fz, double**** storage, int half) {
	for (int i = 0; i < p->nx; i++) {
		for (int j = 0; j < p->ny; j++) {
			storage[half][0][i][j] = fx[i][j];
			storage[half][1][i][j] = fy[i][j];
			storage[half][2][i][j] = fz[i][j];
		}
	}
}
/*************************************************************************************************************/
int* format_primitive(parms p, double**** storage, char *prefix, int downstream_offset, int upstream_offset)
/*************************************************************************************************************/
{
	int j_offset = p.ny - 1;
	double ***stl_buffer = (double***) malloc(3 * sizeof(double**));
	for (int i = 0; i < 3; i++) {
		stl_buffer[i] = (double**) malloc(p.nx * sizeof(double*));
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < p.nx; j++) {
			stl_buffer[i][j] = (double*) malloc(2*p.ny * sizeof(double));
		}
	}
	for (int i = 0; i < p.nx; i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			stl_buffer[x_coord][i][j] = storage[RIGHT][x_coord][i][j] * 1e-3;
			stl_buffer[y_coord][i][j] = storage[RIGHT][y_coord][i][j] * 1e-3;
			stl_buffer[z_coord][i][j] = storage[RIGHT][z_coord][i][j] * 1e-3;
		}
		for (int j = 0; j < p.ny; j++) {
			stl_buffer[x_coord][i][j + j_offset] = storage[LEFT][x_coord][i][j] * 1e-3;
			stl_buffer[y_coord][i][j + j_offset] = storage[LEFT][y_coord][i][j] * 1e-3;
			stl_buffer[z_coord][i][j + j_offset] = storage[LEFT][z_coord][i][j] * 1e-3;
		}
	}
	format_stl(p, stl_buffer, prefix, p.nx - 1, 2 * (p.ny - 1));


	FILE *f, *f_tmp, *f_points, *f_cells;
	char filename[50], txtfile[50];
	int err = sprintf(filename, "%s.vtk", prefix);
	err = sprintf(txtfile, "%s_points.txt", prefix);
	int *info = (int*) malloc(19 * sizeof(int));

	f = fopen(filename, "w+");
	f_points = fopen(txtfile, "w+");
	err = sprintf(txtfile, "%s_cells.txt", prefix);
	f_cells = fopen(txtfile, "w+");

	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "Default Set with interval=1s\n");
	fprintf(f, "ASCII\n");

	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

	int TotalNumberOfPoints = (((p.nx - upstream_offset) - downstream_offset) * (2 * p.ny - 1));
	info[0] = TotalNumberOfPoints;
	info[1] = ((p.nx - upstream_offset) - downstream_offset);
	info[2] = 2 * p.ny - 1;

	int Total_pannels_axially = (((p.nx - upstream_offset) - downstream_offset) - 1);
	int Total_pannels_circumferentially = 2 * (p.ny - 1);
	int TotalNumberOfCells = Total_pannels_axially * Total_pannels_circumferentially;
	info[3] = TotalNumberOfCells;
	info[4] = Total_pannels_axially;
	info[5] = Total_pannels_circumferentially;
	double ***buffer = (double***) malloc(3 * sizeof(double**));
	for (int i = 0; i < 3; i++) {
		buffer[i] = (double**) malloc(((p.nx - upstream_offset) - downstream_offset) * sizeof(double*));
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < ((p.nx - upstream_offset) - downstream_offset); j++) {
			buffer[i][j] = (double*) malloc((2 * p.ny) * sizeof(double));
		}
	}

	j_offset = p.ny - 1;
	fprintf(f, "POINTS %d double\n", TotalNumberOfPoints);
	for (int i = 0 + downstream_offset; i < (p.nx - upstream_offset); i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			fprintf(f, "%2.12lf %2.12lf %2.12lf\n", storage[RIGHT][x_coord][i][j] * 1e-3, storage[RIGHT][y_coord][i][j] * 1e-3,
					storage[RIGHT][z_coord][i][j] * 1e-3);
			fprintf(f_points, "%2.12lf %2.12lf %2.12lf\n", storage[RIGHT][x_coord][i][j] * 1e-3, storage[RIGHT][y_coord][i][j] * 1e-3,
					storage[RIGHT][z_coord][i][j] * 1e-3);
			buffer[x_coord][i - downstream_offset][j] = storage[RIGHT][x_coord][i][j] * 1e-3;
			buffer[y_coord][i - downstream_offset][j] = storage[RIGHT][y_coord][i][j] * 1e-3;
			buffer[z_coord][i - downstream_offset][j] = storage[RIGHT][z_coord][i][j] * 1e-3;
		}

		for (int j = 0; j < p.ny; j++) {
			fprintf(f, "%2.12lf %2.12lf %2.12lf\n", storage[LEFT][x_coord][i][j] * 1e-3, storage[LEFT][y_coord][i][j] * 1e-3,
					storage[LEFT][z_coord][i][j] * 1e-3);
			fprintf(f_points, "%2.12lf %2.12lf %2.12lf\n", storage[LEFT][x_coord][i][j] * 1e-3, storage[LEFT][y_coord][i][j] * 1e-3,
					storage[LEFT][z_coord][i][j] * 1e-3);
			buffer[x_coord][i - downstream_offset][j + j_offset] = storage[LEFT][x_coord][i][j] * 1e-3;
			buffer[y_coord][i - downstream_offset][j + j_offset] = storage[LEFT][y_coord][i][j] * 1e-3;
			buffer[z_coord][i - downstream_offset][j + j_offset] = storage[LEFT][z_coord][i][j] * 1e-3;
		}
	}

	fprintf(f, "CELLS %d %d\n", TotalNumberOfCells, 5 * TotalNumberOfCells);
	int row_offset = 2 * p.ny - 1;
	for (int i = 0 + downstream_offset; i < (p.nx - 1 - upstream_offset); i++) {
		for (int j = 0; j < row_offset - 1; j++) {
			fprintf(f, "%d %d %d %d %d\n", 4, ((i - downstream_offset) * row_offset) + j, (((i - downstream_offset) + 1) * row_offset) + j,
					(((i - downstream_offset) + 1) * row_offset) + j + 1, ((i - downstream_offset) * row_offset) + j + 1);
			fprintf(f_cells, "%d %d %d %d %d\n", 4, ((i - downstream_offset) * row_offset) + j, (((i - downstream_offset) + 1) * row_offset) + j,
					(((i - downstream_offset) + 1) * row_offset) + j + 1, ((i - downstream_offset) * row_offset) + j + 1);
		}
	}
	fprintf(f, "CELL_TYPES %d\n", TotalNumberOfCells);
	for (int i = 0 + downstream_offset; i < (p.nx - 1 - upstream_offset); i++) {
		for (int j = 0; j < row_offset - 1; j++) {
			fprintf(f, "%d\n", 9);
		}
	}

	int* info_smc = (int*) malloc(6 * sizeof(int));
	info_smc = smc_mesh_ver2(p, storage, buffer, Total_pannels_axially, Total_pannels_circumferentially, downstream_offset, upstream_offset, prefix);
	int* info_ec = (int*) malloc(7 * sizeof(int));
	info_ec = ec_mesh_ver2(p, storage, buffer, Total_pannels_axially, Total_pannels_circumferentially, downstream_offset, upstream_offset, prefix);

	info[6] = info_smc[0];
	info[7] = info_smc[1];
	info[8] = info_smc[2];
	info[9] = info_smc[3];
	info[10] = info_smc[4];
	info[11] = info_smc[5];

	info[12] = info_ec[0];
	info[13] = info_ec[1];
	info[14] = info_ec[2];
	info[15] = info_ec[3];
	info[16] = info_ec[4];
	info[17] = info_ec[5];
	info[18] = info_ec[6];
	fclose(f);
	fclose(f_points);
	fclose(f_cells);

	return (info);
}

int* smc_mesh_ver2(parms p, double**** storage, double*** buffer, int Total_pannels_axially, int Total_pannels_circumferentially,
		int downstream_offset, int upstream_offset, char* prefix) {
	int* info = (int*) malloc(6 * sizeof(int));
	mesh_store **mesh = (mesh_store**) malloc(Total_pannels_axially * sizeof(mesh_store*));
	for (int i = 0; i < Total_pannels_axially; i++) {
		mesh[i] = (mesh_store*) malloc(Total_pannels_circumferentially * sizeof(mesh_store));
	}
	int m, n, num_ec_axially = p.ECs, num_smc_circumferentially = p.SMCs;
	/// Meshing SMCs pannels
	m = (num_ec_axially * 13) + 1;
	n = num_smc_circumferentially + 1;

	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			mesh[i][j].x = allocate_dirichlet(m, n);
			mesh[i][j].y = allocate_dirichlet(m, n);
			mesh[i][j].z = allocate_dirichlet(m, n);
		}

	}
	int j_offset = p.ny - 1;
	int indx = 0;
	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			double ax = buffer[x_coord][i][j], bx = buffer[x_coord][i + 1][j], cx = buffer[x_coord][i][j + 1], dx = buffer[x_coord][i + 1][j + 1],
					ay = buffer[y_coord][i][j], by = buffer[y_coord][i + 1][j], cy = buffer[y_coord][i][j + 1], dy = buffer[y_coord][i + 1][j + 1],
					az = buffer[z_coord][i][j], bz = buffer[z_coord][i + 1][j], cz = buffer[z_coord][i][j + 1], dz = buffer[z_coord][i + 1][j + 1];
			for (int p = 0; p < m; p++) {
				int q = 0;
				mesh[i][j].x[p][q] = ax + p * (bx - ax) / (m - 1);
				mesh[i][j].y[p][q] = ay + p * (by - ay) / (m - 1);
				mesh[i][j].z[p][q] = az + p * (bz - az) / (m - 1);
				q = n - 1;
				mesh[i][j].x[p][q] = cx + p * (dx - cx) / (m - 1);
				mesh[i][j].y[p][q] = cy + p * (dy - cy) / (m - 1);
				mesh[i][j].z[p][q] = cz + p * (dz - cz) / (m - 1);
			}
			for (int q = 0; q < n; q++) {
				int p = 0;
				mesh[i][j].x[p][q] = ax + q * (cx - ax) / (n - 1);
				mesh[i][j].y[p][q] = ay + q * (cy - ay) / (n - 1);
				mesh[i][j].z[p][q] = az + q * (cz - az) / (n - 1);
				p = m - 1;
				mesh[i][j].x[p][q] = bx + q * (dx - bx) / (n - 1);
				mesh[i][j].y[p][q] = by + q * (dy - by) / (n - 1);
				mesh[i][j].z[p][q] = bz + q * (dz - bz) / (n - 1);
			}

			for (int p = 1; p < m - 1; p++) {
				for (int q = 1; q < n - 1; q++) {
					mesh[i][j].x[p][q] = mesh[i][j].x[p][0] + q * (mesh[i][j].x[p][n - 1] - mesh[i][j].x[p][0]) / (n - 1);
					mesh[i][j].y[p][q] = mesh[i][j].y[p][0] + q * (mesh[i][j].y[p][n - 1] - mesh[i][j].y[p][0]) / (n - 1);
					mesh[i][j].z[p][q] = mesh[i][j].z[p][0] + q * (mesh[i][j].z[p][n - 1] - mesh[i][j].z[p][0]) / (n - 1);
				}
			}

			indx++;

		}
		for (int j = 0; j < p.ny - 1; j++) {
			double ax = buffer[x_coord][i][j + j_offset], bx = buffer[x_coord][i + 1][j + j_offset], cx = buffer[x_coord][i][j + j_offset + 1], dx =
					buffer[x_coord][i + 1][j + j_offset + 1], ay = buffer[y_coord][i][j + j_offset], by = buffer[y_coord][i + 1][j + j_offset], cy =
					buffer[y_coord][i][j + j_offset + 1], dy = buffer[y_coord][i + 1][j + j_offset + 1], az = buffer[z_coord][i][j + j_offset], bz =
					buffer[z_coord][i + 1][j + j_offset], cz = buffer[z_coord][i][j + j_offset + 1], dz = buffer[z_coord][i + 1][j + j_offset + 1];
			for (int p = 0; p < m; p++) {
				int q = 0;
				mesh[i][j + j_offset].x[p][q] = ax + p * (bx - ax) / (m - 1);
				mesh[i][j + j_offset].y[p][q] = ay + p * (by - ay) / (m - 1);
				mesh[i][j + j_offset].z[p][q] = az + p * (bz - az) / (m - 1);
				q = n - 1;
				mesh[i][j + j_offset].x[p][q] = cx + p * (dx - cx) / (m - 1);
				mesh[i][j + j_offset].y[p][q] = cy + p * (dy - cy) / (m - 1);
				mesh[i][j + j_offset].z[p][q] = cz + p * (dz - cz) / (m - 1);
			}
			for (int q = 0; q < n; q++) {
				int p = 0;
				mesh[i][j + j_offset].x[p][q] = ax + q * (cx - ax) / (n - 1);
				mesh[i][j + j_offset].y[p][q] = ay + q * (cy - ay) / (n - 1);
				mesh[i][j + j_offset].z[p][q] = az + q * (cz - az) / (n - 1);
				p = m - 1;
				mesh[i][j + j_offset].x[p][q] = bx + q * (dx - bx) / (n - 1);
				mesh[i][j + j_offset].y[p][q] = by + q * (dy - by) / (n - 1);
				mesh[i][j + j_offset].z[p][q] = bz + q * (dz - bz) / (n - 1);
			}

			for (int p = 1; p < m - 1; p++) {
				for (int q = 1; q < n - 1; q++) {
					mesh[i][j + j_offset].x[p][q] = mesh[i][j + j_offset].x[p][0]
							+ q * (mesh[i][j + j_offset].x[p][n - 1] - mesh[i][j + j_offset].x[p][0]) / (n - 1);
					mesh[i][j + j_offset].y[p][q] = mesh[i][j + j_offset].y[p][0]
							+ q * (mesh[i][j + j_offset].y[p][n - 1] - mesh[i][j + j_offset].y[p][0]) / (n - 1);
					mesh[i][j + j_offset].z[p][q] = mesh[i][j + j_offset].z[p][0]
							+ q * (mesh[i][j + j_offset].z[p][n - 1] - mesh[i][j + j_offset].z[p][0]) / (n - 1);
				}
			}
			indx++;
		}
	}

	/// Now writing a VTK file
	FILE *fw, *f_points, *f_cells;

	char filename[50], fname_points[50], fname_cells[50];
	sprintf(filename, "smc_mesh_%s.vtk", prefix);
	sprintf(fname_points, "%s_smc_mesh_points.txt", prefix);
	sprintf(fname_cells, "%s_smc_mesh_cells.txt", prefix);
	fw = fopen(filename, "w+");
	f_points = fopen(fname_points, "w+");
	f_cells = fopen(fname_cells, "w+");

	fprintf(fw, "# vtk DataFile Version 3.0\n");
	fprintf(fw, "SMC mesh\n");
	fprintf(fw, "ASCII\n");
	fprintf(fw, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fw, "POINTS %d double\n", indx * m * n);

	info[0] = m * n;
	info[1] = m;
	info[2] = n;

	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			for (int p = 0; p < m; p++) {
				for (int q = 0; q < n; q++) {
					fprintf(fw, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j].x[p][q], mesh[i][j].y[p][q], mesh[i][j].z[p][q]);
					fprintf(f_points, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j].x[p][q], mesh[i][j].y[p][q], mesh[i][j].z[p][q]);
				}
			}
		}
		for (int j = 0; j < p.ny - 1; j++) {
			for (int p = 0; p < m; p++) {
				for (int q = 0; q < n; q++) {
					fprintf(fw, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j + j_offset].x[p][q], mesh[i][j + j_offset].y[p][q],
							mesh[i][j + j_offset].z[p][q]);
					fprintf(f_points, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j + j_offset].x[p][q], mesh[i][j + j_offset].y[p][q],
							mesh[i][j + j_offset].z[p][q]);

				}
			}
		}
	}
	int cell_offset = 0;
	int num_cells = Total_pannels_axially * Total_pannels_circumferentially * (m - 1) * (n - 1);

	info[3] = (m - 1) * (n - 1);
	info[4] = m - 1;
	info[5] = n - 1;

	int row_offset = n;
	fprintf(fw, "CELLS %d %d\n", num_cells, num_cells * 5);
	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			for (int p = 0; p < m - 1; p++) {
				for (int q = 0; q < n - 1; q++) {
					fprintf(fw, "%d %d %d %d %d\n", 4, cell_offset + p * row_offset + q, cell_offset + p * row_offset + (q + 1),
							cell_offset + (p + 1) * row_offset + (q + 1), cell_offset + (p + 1) * row_offset + q);
					fprintf(f_cells, "%d %d %d %d %d\n", 4, cell_offset + p * row_offset + q, cell_offset + p * row_offset + (q + 1),
							cell_offset + (p + 1) * row_offset + (q + 1), cell_offset + (p + 1) * row_offset + q);
				}
			}
			cell_offset += m * n;
		}
	}
	fprintf(fw, "CELL_TYPES %d\n", num_cells);
	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			for (int p = 0; p < m - 1; p++) {
				for (int q = 0; q < n - 1; q++) {
					fprintf(fw, "%d \n", 9);
				}
			}
		}
	}
	fclose(fw);
	fclose(f_points);
	fclose(f_cells);
	return (info);
}
int* ec_mesh_ver2(parms p, double**** storage, double*** buffer, int Total_pannels_axially, int Total_pannels_circumferentially,
		int downstream_offset, int upstream_offset, char* prefix) {
	int* info = (int*) malloc(7 * sizeof(int));
	mesh_store **mesh = (mesh_store**) malloc(Total_pannels_axially * sizeof(mesh_store*));
	for (int i = 0; i < Total_pannels_axially; i++) {
		mesh[i] = (mesh_store*) malloc(Total_pannels_circumferentially * sizeof(mesh_store));
	}
	int m, n, num_ec_axially = p.ECs, num_smc_circumferentially = p.SMCs;
	/// Meshing SMCs pannels
	m = num_ec_axially + 1;
	n = num_smc_circumferentially * 5 + 1;

	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			mesh[i][j].x = allocate_dirichlet(m, n);
			mesh[i][j].y = allocate_dirichlet(m, n);
			mesh[i][j].z = allocate_dirichlet(m, n);
			mesh[i][j].cx = allocate_dirichlet(m, n);
			mesh[i][j].cy = allocate_dirichlet(m, n);
			mesh[i][j].cz = allocate_dirichlet(m, n);

		}

	}
	int j_offset = p.ny - 1;
	int indx = 0;
	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			double ax = buffer[x_coord][i][j], bx = buffer[x_coord][i + 1][j], cx = buffer[x_coord][i][j + 1], dx = buffer[x_coord][i + 1][j + 1],
					ay = buffer[y_coord][i][j], by = buffer[y_coord][i + 1][j], cy = buffer[y_coord][i][j + 1], dy = buffer[y_coord][i + 1][j + 1],
					az = buffer[z_coord][i][j], bz = buffer[z_coord][i + 1][j], cz = buffer[z_coord][i][j + 1], dz = buffer[z_coord][i + 1][j + 1];
			for (int p = 0; p < m; p++) {
				int q = 0;
				mesh[i][j].x[p][q] = ax + p * (bx - ax) / (m - 1);
				mesh[i][j].y[p][q] = ay + p * (by - ay) / (m - 1);
				mesh[i][j].z[p][q] = az + p * (bz - az) / (m - 1);
				q = n - 1;
				mesh[i][j].x[p][q] = cx + p * (dx - cx) / (m - 1);
				mesh[i][j].y[p][q] = cy + p * (dy - cy) / (m - 1);
				mesh[i][j].z[p][q] = cz + p * (dz - cz) / (m - 1);
			}
			for (int q = 0; q < n; q++) {
				int p = 0;
				mesh[i][j].x[p][q] = ax + q * (cx - ax) / (n - 1);
				mesh[i][j].y[p][q] = ay + q * (cy - ay) / (n - 1);
				mesh[i][j].z[p][q] = az + q * (cz - az) / (n - 1);
				p = m - 1;
				mesh[i][j].x[p][q] = bx + q * (dx - bx) / (n - 1);
				mesh[i][j].y[p][q] = by + q * (dy - by) / (n - 1);
				mesh[i][j].z[p][q] = bz + q * (dz - bz) / (n - 1);
			}

			for (int p = 1; p < m - 1; p++) {
				for (int q = 1; q < n - 1; q++) {
					mesh[i][j].x[p][q] = mesh[i][j].x[p][0] + q * (mesh[i][j].x[p][n - 1] - mesh[i][j].x[p][0]) / (n - 1);
					mesh[i][j].y[p][q] = mesh[i][j].y[p][0] + q * (mesh[i][j].y[p][n - 1] - mesh[i][j].y[p][0]) / (n - 1);
					mesh[i][j].z[p][q] = mesh[i][j].z[p][0] + q * (mesh[i][j].z[p][n - 1] - mesh[i][j].z[p][0]) / (n - 1);
				}
			}

			indx++;

		}
		for (int j = 0; j < p.ny - 1; j++) {
			double ax = buffer[x_coord][i][j + j_offset], bx = buffer[x_coord][i + 1][j + j_offset], cx = buffer[x_coord][i][j + j_offset + 1], dx =
					buffer[x_coord][i + 1][j + j_offset + 1], ay = buffer[y_coord][i][j + j_offset], by = buffer[y_coord][i + 1][j + j_offset], cy =
					buffer[y_coord][i][j + j_offset + 1], dy = buffer[y_coord][i + 1][j + j_offset + 1], az = buffer[z_coord][i][j + j_offset], bz =
					buffer[z_coord][i + 1][j + j_offset], cz = buffer[z_coord][i][j + j_offset + 1], dz = buffer[z_coord][i + 1][j + j_offset + 1];
			for (int p = 0; p < m; p++) {
				int q = 0;
				mesh[i][j + j_offset].x[p][q] = ax + p * (bx - ax) / (m - 1);
				mesh[i][j + j_offset].y[p][q] = ay + p * (by - ay) / (m - 1);
				mesh[i][j + j_offset].z[p][q] = az + p * (bz - az) / (m - 1);
				q = n - 1;
				mesh[i][j + j_offset].x[p][q] = cx + p * (dx - cx) / (m - 1);
				mesh[i][j + j_offset].y[p][q] = cy + p * (dy - cy) / (m - 1);
				mesh[i][j + j_offset].z[p][q] = cz + p * (dz - cz) / (m - 1);
			}
			for (int q = 0; q < n; q++) {
				int p = 0;
				mesh[i][j + j_offset].x[p][q] = ax + q * (cx - ax) / (n - 1);
				mesh[i][j + j_offset].y[p][q] = ay + q * (cy - ay) / (n - 1);
				mesh[i][j + j_offset].z[p][q] = az + q * (cz - az) / (n - 1);
				p = m - 1;
				mesh[i][j + j_offset].x[p][q] = bx + q * (dx - bx) / (n - 1);
				mesh[i][j + j_offset].y[p][q] = by + q * (dy - by) / (n - 1);
				mesh[i][j + j_offset].z[p][q] = bz + q * (dz - bz) / (n - 1);
			}

			for (int p = 1; p < m - 1; p++) {
				for (int q = 1; q < n - 1; q++) {
					mesh[i][j + j_offset].x[p][q] = mesh[i][j + j_offset].x[p][0]
							+ q * (mesh[i][j + j_offset].x[p][n - 1] - mesh[i][j + j_offset].x[p][0]) / (n - 1);
					mesh[i][j + j_offset].y[p][q] = mesh[i][j + j_offset].y[p][0]
							+ q * (mesh[i][j + j_offset].y[p][n - 1] - mesh[i][j + j_offset].y[p][0]) / (n - 1);
					mesh[i][j + j_offset].z[p][q] = mesh[i][j + j_offset].z[p][0]
							+ q * (mesh[i][j + j_offset].z[p][n - 1] - mesh[i][j + j_offset].z[p][0]) / (n - 1);
				}
			}
			indx++;
		}
	}

	/// Now writing a VTK file
	FILE *fw, *f_points, *f_cells;

	char filename[50], fname_points[50], fname_cells[50];
	sprintf(filename, "ec_mesh_%s.vtk", prefix);
	sprintf(fname_points, "%s_ec_mesh_points.txt", prefix);
	sprintf(fname_cells, "%s_ec_mesh_cells.txt", prefix);
	fw = fopen(filename, "w+");
	f_points = fopen(fname_points, "w+");
	f_cells = fopen(fname_cells, "w+");

	fprintf(fw, "# vtk DataFile Version 3.0\n");
	fprintf(fw, "EC mesh\n");
	fprintf(fw, "ASCII\n");
	fprintf(fw, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fw, "POINTS %d double\n", indx * m * n);

	info[0] = m * n;
	info[1] = m;
	info[2] = n;

	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < p.ny - 1; j++) {
			for (int p = 0; p < m; p++) {
				for (int q = 0; q < n; q++) {
					fprintf(fw, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j].x[p][q], mesh[i][j].y[p][q], mesh[i][j].z[p][q]);
					fprintf(f_points, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j].x[p][q], mesh[i][j].y[p][q], mesh[i][j].z[p][q]);
				}
			}
		}
		for (int j = 0; j < p.ny - 1; j++) {
			for (int p = 0; p < m; p++) {
				for (int q = 0; q < n; q++) {
					fprintf(fw, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j + j_offset].x[p][q], mesh[i][j + j_offset].y[p][q],
							mesh[i][j + j_offset].z[p][q]);
					fprintf(f_points, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j + j_offset].x[p][q], mesh[i][j + j_offset].y[p][q],
							mesh[i][j + j_offset].z[p][q]);
				}
			}
		}
	}
	int cell_offset = 0;
	int num_cells = Total_pannels_axially * Total_pannels_circumferentially * (m - 1) * (n - 1);
	info[3] = (m - 1) * (n - 1);
	info[4] = m - 1;
	info[5] = n - 1;
	int row_offset = n;
	fprintf(fw, "CELLS %d %d\n", num_cells, num_cells * 5);
	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			for (int p = 0; p < m - 1; p++) {
				for (int q = 0; q < n - 1; q++) {
					fprintf(fw, "%d %d %d %d %d\n", 4, cell_offset + p * row_offset + q, cell_offset + p * row_offset + (q + 1),
							cell_offset + (p + 1) * row_offset + (q + 1), cell_offset + (p + 1) * row_offset + q);
					fprintf(f_cells, "%d %d %d %d %d\n", 4, cell_offset + p * row_offset + q, cell_offset + p * row_offset + (q + 1),
							cell_offset + (p + 1) * row_offset + (q + 1), cell_offset + (p + 1) * row_offset + q);
					mesh[i][j].cx[p][q] = (mesh[i][j].x[p][q] + mesh[i][j].x[p][q + 1] /*+ mesh[i][j].x[(p + 1)][q + 1], mesh[i][j].x[(p + 1)][q]*/)
							/ 2;
					mesh[i][j].cy[p][q] = (mesh[i][j].y[p][q] + mesh[i][j].y[p + 1][q]) / 2; /*(mesh[i][j].y[p][q] + mesh[i][j].y[p][q + 1] + mesh[i][j].y[(p + 1)][q + 1]+ mesh[i][j].y[(p + 1)][q])
					 / 2;*/
					mesh[i][j].cz[p][q] = (mesh[i][j].z[p][q] + mesh[i][j].z[p][q + 1]) / 2;//(mesh[i][j].z[p][q] + mesh[i][j].z[p][q + 1] + mesh[i][j].z[(p + 1)][q + 1], mesh[i][j].z[(p + 1)][q])/ 2;
				}
			}
			cell_offset += m * n;
		}
	}
	fprintf(fw, "CELL_TYPES %d\n", num_cells);
	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			for (int p = 0; p < m - 1; p++) {
				for (int q = 0; q < n - 1; q++) {
					fprintf(fw, "%d \n", 9);
				}
			}
		}
	}
	fclose(fw);
	fclose(f_points);
	fclose(f_cells);

	char centeroid_file[50];
	sprintf(centeroid_file, "ec_centeroid_%s.vtk", prefix);
	sprintf(fname_points, "%s_ec_centeroid_points.txt", prefix);
	sprintf(fname_cells, "%s_ec_centeroid_cells.txt", prefix);

	fw = fopen(centeroid_file, "w+");
	f_points = fopen(fname_points, "w+");
	f_cells = fopen(fname_cells, "w+");

	fprintf(fw, "# vtk DataFile Version 3.0\n");
	fprintf(fw, "EC centeroids % branch\n", prefix);
	fprintf(fw, "ASCII\n");
	fprintf(fw, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fw, "POINTS %d double\n", num_cells);

	info[6] = num_cells / (Total_pannels_axially * Total_pannels_circumferentially);

	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			for (int p = 0; p < m - 1; p++) {
				for (int q = 0; q < n - 1; q++) {
					fprintf(fw, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j].cx[p][q], mesh[i][j].cy[p][q], mesh[i][j].cz[p][q]);
					fprintf(f_points, "%2.8lf %2.8lf %2.8lf\n", mesh[i][j].cx[p][q], mesh[i][j].cy[p][q], mesh[i][j].cz[p][q]);
				}
			}
		}
	}
	int count = 0;
	fprintf(fw, "CELLS %d %d\n", num_cells, num_cells * 2);
	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			for (int p = 0; p < m - 1; p++) {
				for (int q = 0; q < n - 1; q++) {
					fprintf(fw, "%d %d\n", 1, count);
					fprintf(f_cells, "%d %d\n", 1, count);
					count++;
				}
			}
		}
	}
	fprintf(fw, "CELL_TYPES %d\n", num_cells);
	for (int i = 0; i < Total_pannels_axially; i++) {
		for (int j = 0; j < Total_pannels_circumferentially; j++) {
			for (int p = 0; p < m - 1; p++) {
				for (int q = 0; q < n - 1; q++) {
					fprintf(fw, "%d \n", 1);
				}
			}
		}
	}
	fclose(fw);
	fclose(f_points);
	fclose(f_cells);

	return (info);
}

/***********************************************************************/

double* dirichlet_boundary_outer_wall_corrected_v_boundary(parms* p, double**** storage, double **fx, double **fy, double **fz, int half) {
	for (int i = 1; i <= p->nx; i++) {
		double x = p->a + (i - 1) * (p->b - p->a) / (p->m + 1);
		int j;
		double h = p->h3 - (p->h3 - p->h2) * x;
		double r = p->neck - (p->neck - p->rt) * (1 - x);
		///v=c
		j = 1;
		fx[i - 1][j - 1] = storage[half][x_coord][i - 1][j - 1];
		fy[i - 1][j - 1] = storage[half][y_coord][i - 1][j - 1];
		fz[i - 1][j - 1] = storage[half][z_coord][i - 1][j - 1];
		///v=d
		j = p->ny;
		fx[i - 1][j - 1] = storage[half][x_coord][i - 1][j - 1];
		fy[i - 1][j - 1] = storage[half][y_coord][i - 1][j - 1];
		fz[i - 1][j - 1] = storage[half][z_coord][i - 1][j - 1];
	}

	double* theta_val = (double*) malloc(p->ny * sizeof(double));
	for (int j = 1; j <= p->ny; j++) {
		double y = p->c + (j - 1) * (p->d - p->c) / (p->n + 1);
		double direction = 0.0;
		theta_val[j - 1] = y;
		int i;
		///x=a
		i = 1;
		fx[i - 1][j - 1] = -(p->h3 * cos(p->alfa / 2) + p->ra * sin(p->alfa / 2) * cos(y));
		fy[i - 1][j - 1] = p->h3 * sin(p->alfa / 2) - p->ra * cos(p->alfa / 2) * cos(y);
		fz[i - 1][j - 1] = p->ra * sin(y);

		///x=b
		i = p->nx;
		fx[i - 1][j - 1] = -p->neck * cos(y);
		fy[i - 1][j - 1] = p->h2;
		fz[i - 1][j - 1] = p->neck * sin(y);
	}
	return (theta_val);
}

double* dirichlet_boundary_inner_wall_corrected_v_boundary(parms* p, double**** storage, double **fx, double **fy, double **fz, int half) {
	for (int i = 1; i <= p->nx; i++) {
		double x = p->a + (i - 1) * (p->b - p->a) / (p->m + 1);
		int j;
		///v=c
		double h = p->h3 - (p->h3 - p->h2) * x;
		double r = p->neck - (p->neck - p->rt) * (1 - x);
		j = 1;
		fx[i - 1][j - 1] = storage[half][x_coord][i - 1][j - 1];
		fy[i - 1][j - 1] = storage[half][y_coord][i - 1][j - 1];
		fz[i - 1][j - 1] = storage[half][z_coord][i - 1][j - 1];

		///v=d
		j = p->ny;
		fx[i - 1][j - 1] = storage[half][x_coord][i - 1][j - 1];
		fy[i - 1][j - 1] = storage[half][y_coord][i - 1][j - 1];
		fz[i - 1][j - 1] = storage[half][z_coord][i - 1][j - 1];
	}

	double* theta_val = (double*) malloc(p->ny * sizeof(double));
	for (int j = 1; j <= p->ny; j++) {
		double y = p->c + (j - 1) * (p->d - p->c) / (p->n + 1);
		double direction = 0.0;
		theta_val[j - 1] = y;

		int i;
///x=a
		i = 1;
		fx[i - 1][j - 1] = -(p->h3 * cos(p->alfa / 2) + p->ra * sin(p->alfa / 2) * cos(y));
		fy[i - 1][j - 1] = p->h3 * sin(p->alfa / 2) - p->ra * cos(p->alfa / 2) * cos(y);
		fz[i - 1][j - 1] = p->ra * sin(y);
///x=b
		i = p->nx;

		fx[i - 1][j - 1] = 0.0;
		if ((y >= pi / 2) && (y <= 3 * pi / 2)) {
			fy[i - 1][j - 1] = -p->apex * cos(y);
		} else if ((y >= 3 * pi / 2) && (y <= 2 * pi + pi / 2)) {
			fy[i - 1][j - 1] = p->apex * cos(y);
		}

		fz[i - 1][j - 1] = p->neck * sin(y);
	}
	return (theta_val);
}

void bound_v_correction_ver3(parms p, double**** storage, double** theta_val, int downstream_offset, int upstream_offset) {

	for (int i = 1; i < p.nx - 1; i++) {
		int ja0, ja1, ja2, j, jb0, jb1, jb2;

		double xa0, xa1, xa2, x, xb0, xb1, xb2;
		double ya0, ya1, ya2, y, yb0, yb1, yb2;
		double za0, za1, za2, z, zb0, zb1, zb2;
		double b0, b1, b2, b3, b4, b5;
		double c, c0, c1, c2, c3, c4, c5;

///For RIGHT Half start point
		ja0 = 3;
		ja1 = 2;
		ja2 = 1;
		j = 0;
		jb0 = p.ny - 1 - 1;
		jb1 = p.ny - 1 - 2;
		jb2 = p.ny - 1 - 3;
		c = theta_val[RIGHT][j];
		c0 = theta_val[RIGHT][ja0];
		c1 = theta_val[RIGHT][ja1];
		c2 = theta_val[RIGHT][ja2];
		c3 = theta_val[LEFT][jb0];
		c4 = theta_val[LEFT][jb1];
		c5 = theta_val[LEFT][jb2];

		x = storage[RIGHT][x_coord][i][j];
		xa0 = storage[RIGHT][x_coord][i][ja0];
		xa1 = storage[RIGHT][x_coord][i][ja1];
		xa2 = storage[RIGHT][x_coord][i][ja2];
		xb0 = storage[LEFT][x_coord][i][jb0];
		xb1 = storage[LEFT][x_coord][i][jb1];
		xb2 = storage[LEFT][x_coord][i][jb2];
		b0 = xa0;
		b1 = (xa1 - xa0) / (c1 - c0);
		b2 = (((xa2 - xa1) / (c2 - c1)) - b1) / (c2 - c0);
		b3 = (((xb0 - xa2) / (c3 - c2)) - b2 - b1) / (c3 - c0);
		b4 = (((xb1 - xb0) / (c4 - c3)) - b3 - b2 - b1) / (c4 - c0);
		b5 = (((xb2 - xb1) / (c5 - c4)) - b4 - b3 - b2 - b1) / (c5 - c0);
		x = b0 + b1 * (c - c0) + b2 * (c - c0) * (c - c1) + b3 * (c - c0) * (c - c1) * (c - c2) + b4 * (c - c0) * (c - c1) * (c - c2) * (c - c3)
				+ b5 * (c - c0) * (c - c1) * (c - c2) * (c - c3) * (c - c4);
		storage[RIGHT][x_coord][i][j] = x;
		storage[LEFT][x_coord][i][p.ny - 1] = x;

		y = storage[RIGHT][y_coord][i][j];
		ya0 = storage[RIGHT][y_coord][i][ja0];
		ya1 = storage[RIGHT][y_coord][i][ja1];
		ya2 = storage[RIGHT][y_coord][i][ja2];
		yb0 = storage[LEFT][y_coord][i][jb0];
		yb1 = storage[LEFT][y_coord][i][jb1];
		yb2 = storage[LEFT][y_coord][i][jb2];
		b0 = ya0;
		b1 = (ya1 - ya0) / (c1 - c0);
		b2 = (((ya2 - ya1) / (c2 - c1)) - b1) / (c2 - c0);
		b3 = (((yb0 - ya2) / (c3 - c2)) - b2 - b1) / (c3 - c0);
		b4 = (((yb1 - yb0) / (c4 - c3)) - b3 - b2 - b1) / (c4 - c0);
		b5 = (((yb2 - yb1) / (c5 - c4)) - b4 - b3 - b2 - b1) / (c5 - c0);
		y = b0 + b1 * (c - c0) + b2 * (c - c0) * (c - c1) + b3 * (c - c0) * (c - c1) * (c - c2) + b4 * (c - c0) * (c - c1) * (c - c2) * (c - c3)
				+ b5 * (c - c0) * (c - c1) * (c - c2) * (c - c3) * (c - c4);
		storage[RIGHT][y_coord][i][j] = y;
		storage[LEFT][y_coord][i][p.ny - 1] = y;

		z = storage[RIGHT][z_coord][i][j];
		/*xa0 = storage[RIGHT][z_coord][i][ja0];
		 xa1 = storage[RIGHT][z_coord][i][ja1];
		 xa2 = storage[RIGHT][z_coord][i][ja2];
		 xb0 = storage[LEFT][z_coord][i][jb0];
		 xb1 = storage[LEFT][z_coord][i][jb1];
		 xb2 = storage[LEFT][z_coord][i][jb2];
		 b0 = xa0;
		 b1 = (xa1 - xa0) / (c1 - c0);
		 b2 = (((xa2 - xa1) / (c2 - c1)) - b1) / (c2 - c0);
		 b3 = (((xb0 - xa2) / (c3 - c2)) - b2 - b1) / (c3 - c0);
		 b4 = (((xb1 - xb0) / (c4 - c3)) - b3 - b2 - b1) / (c4 - c0);
		 b5 = (((xb2 - xb1) / (c5 - c4)) - b4 - b3 - b2 - b1) / (c5 - c0);
		 z = b0 + b1 * (c - c0) + b2 * (c - c0) * (c - c1) + b3 * (c - c0) * (c - c1) * (c - c2) + b4 * (c - c0) * (c - c1) * (c - c2) * (c - c3)
		 + b5 * (c - c0) * (c - c1) * (c - c2) * (c - c3) * (c - c4);*/
		storage[RIGHT][z_coord][i][j] = z;
		storage[LEFT][z_coord][i][p.ny - 1] = z;
		printf("%2.7lf\t%2.7lf\t%2.7lf\n", 1e-3 * x, 1e-3 * y, 1e-3 * z);
///For LEFT Half end point
		/*	ja0 = p.ny - 1 - 2;
		 ja1 = p.ny - 1 - 1;
		 j = p.ny - 1;
		 jb0 = 1;
		 jb1 = 2;
		 c = theta_val[LEFT][j];
		 c0 = theta_val[LEFT][ja0];
		 c1 = theta_val[LEFT][ja1];
		 c2 = theta_val[RIGHT][jb0];
		 c3 = theta_val[RIGHT][jb1];

		 x = storage[LEFT][x_coord][i][j];
		 xa0 = storage[LEFT][x_coord][i][ja0];
		 xa1 = storage[LEFT][x_coord][i][ja1];
		 xb0 = storage[RIGHT][x_coord][i][jb0];
		 xb1 = storage[RIGHT][x_coord][i][jb1];
		 b0 = xa0;
		 b1 = (xa1 - xa0) / (c1 - c0);
		 b2 = (((xb0 - xa1) / (c2 - c1)) - ((xa1 - xa0) / (c1 - c0))) / (c2 - c0);
		 b3 = (((xb1 - xb0) / (c3 - c2)) - ((xb0 - xa1) / (c2 - c1)) - ((xa1 - xa0) / (c1 - c0))) / (c3 - c0);
		 x = b0 + b1 * (c - c0) + b2 * (c - c0) * (c - c1) + b3 * (c - c0) * (c - c1) * (c - c2);
		 storage[LEFT][x_coord][i][j] = x;

		 y = storage[LEFT][y_coord][i][j];
		 ya0 = storage[LEFT][y_coord][i][ja0];
		 ya1 = storage[LEFT][y_coord][i][ja1];
		 yb0 = storage[RIGHT][y_coord][i][jb0];
		 yb1 = storage[RIGHT][y_coord][i][jb1];
		 b0 = ya0;
		 b1 = (ya1 - ya0) / (c1 - c0);
		 b2 = (((yb0 - ya1) / (c2 - c1)) - ((ya1 - ya0) / (c1 - c0))) / (c2 - c0);
		 b3 = (((yb1 - yb0) / (c3 - c2)) - ((yb0 - ya1) / (c2 - c1)) - ((ya1 - ya0) / (c1 - c0))) / (c3 - c0);
		 y = b0 + b1 * (c - c0) + b2 * (c - c0) * (c - c1) + b3 * (c - c0) * (c - c1) * (c - c2);
		 storage[LEFT][y_coord][i][j] = y;

		 z = storage[LEFT][z_coord][i][j];
		 za0 = storage[LEFT][z_coord][i][ja0];
		 za1 = storage[LEFT][z_coord][i][ja1];
		 zb0 = storage[RIGHT][z_coord][i][jb0];
		 zb1 = storage[RIGHT][z_coord][i][jb1];
		 b0 = za0;
		 b1 = (za1 - za0) / (c1 - c0);
		 b2 = (((zb0 - za1) / (c2 - c1)) - ((za1 - za0) / (c1 - c0))) / (c2 - c0);
		 b3 = (((zb1 - zb0) / (c3 - c2)) - ((zb0 - za1) / (c2 - c1)) - ((za1 - za0) / (c1 - c0))) / (c3 - c0);
		 z = b0 + b1 * (c - c0) + b2 * (c - c0) * (c - c1) + b3 * (c - c0) * (c - c1) * (c - c2);
		 //storage[LEFT][z_coord][i][j] = z;*/
	}
}

void bound_v_correction(parms* p, double**** storage, double** theta_val) {

	for (int i = 1; i < p->nx - 1; i++) {
		int ja0, ja1, ja2, j, jb0, jb1, jb2;

		double xa0, x, xb0;
		double ya0, y, yb0;
		double za0, z, zb0;
		double b0, b1;
		double c, c0, c1;
		double mu = 0.5;
///For RIGHT Half start point
		ja0 = p->ny - 1 - 1;
		j = 0;
		jb0 = 1;
		c = theta_val[RIGHT][j];
		c0 = theta_val[LEFT][ja0];
		c1 = theta_val[RIGHT][jb0];

		x = storage[RIGHT][x_coord][i][j];
		xa0 = storage[LEFT][x_coord][i][ja0];
		xb0 = storage[RIGHT][x_coord][i][jb0];
		b0 = xa0;
		b1 = (xb0 - xa0) / (c1 - c0);
		x = b0 + ((xb0 - xa0) * mu);		// (c - c0) / (c1 - c0)); //b0 + b1 * (c - c0);
		storage[RIGHT][x_coord][i][j] = x;
		storage[LEFT][x_coord][i][p->ny - 1] = x;

		y = storage[RIGHT][y_coord][i][j];
		ya0 = storage[LEFT][y_coord][i][ja0];
		yb0 = storage[RIGHT][y_coord][i][jb0];

		b0 = ya0;
		b1 = (yb0 - ya0) / (c1 - c0);
		y = b0 + ((yb0 - ya0) * mu);		//(c - c0) / (c1 - c0));//b0 + b1 * (c - c0);
		storage[RIGHT][y_coord][i][j] = y;
		storage[LEFT][y_coord][i][p->ny - 1] = y;

		z = storage[RIGHT][z_coord][i][j];
		za0 = storage[LEFT][z_coord][i][ja0];
		zb0 = storage[RIGHT][z_coord][i][jb0];

		b0 = za0;
		b1 = (zb0 - za0) / (c1 - c0);
		z = b0 + ((zb0 - za0) * mu);		//(c - c0) / (c1 - c0));//b0 + b1 * (c - c0);
		storage[RIGHT][z_coord][i][j] = z;
		///For LEFT Half start point
		storage[LEFT][z_coord][i][p->ny - 1] = z;

		///For RIGHT Half end point
		ja0 = p->ny - 1 - 1;
		j = p->ny - 1;
		jb0 = 1;
		c = theta_val[RIGHT][j];
		c0 = theta_val[RIGHT][ja0];
		c1 = theta_val[LEFT][jb0];

		x = storage[RIGHT][x_coord][i][j];
		xa0 = storage[RIGHT][x_coord][i][ja0];
		xb0 = storage[LEFT][x_coord][i][jb0];
		b0 = xa0;
		b1 = (xb0 - xa0) / (c1 - c0);
		x = b0 + ((xb0 - xa0) * mu);		// (c - c0) / (c1 - c0)); //b0 + b1 * (c - c0);
		storage[RIGHT][x_coord][i][j] = x;
		storage[LEFT][x_coord][i][0] = x;

		y = storage[RIGHT][y_coord][i][j];
		ya0 = storage[RIGHT][y_coord][i][ja0];
		yb0 = storage[LEFT][y_coord][i][jb0];

		b0 = ya0;
		b1 = (yb0 - ya0) / (c1 - c0);
		y = b0 + ((yb0 - ya0) * mu);		//(c - c0) / (c1 - c0));//b0 + b1 * (c - c0);
		storage[RIGHT][y_coord][i][j] = y;
		storage[LEFT][y_coord][i][0] = y;

		z = storage[RIGHT][z_coord][i][j];
		za0 = storage[RIGHT][z_coord][i][ja0];
		zb0 = storage[LEFT][z_coord][i][jb0];

		b0 = za0;
		b1 = (zb0 - za0) / (c1 - c0);
		z = b0 + ((zb0 - za0) * mu);		//(c - c0) / (c1 - c0));//b0 + b1 * (c - c0);
		storage[RIGHT][z_coord][i][j] = z;
		///For LEFT Half end point
		storage[LEFT][z_coord][i][0] = z;

	}
}

void itereate_to_correct(int branch, parms *p, double** fx, double** fy, double** fz, double**** storage, double* bda, double* bdb, double* bdc,
		double* bdd, int idf, int iflag, double tol, int itcg, double* w, int lw, double** theta_val) {
	if (branch == P) {
		//Do nothing
	} else if (branch == L) {
		for (int corrections = 0; corrections < p->max_corrections; corrections++) {
			p->alfa = pi - p->angle;
			p->c = 3 * pi / 2;
			p->d = 2 * pi + pi / 2;
			for (int i = 0; i < 2; i++) {
				theta_val[i] = (double*) malloc(p->ny * sizeof(double));
			}
			theta_val[LEFT] = dirichlet_boundary_outer_wall_corrected_v_boundary(p, storage, fx, fy, fz, LEFT);
			set_neumann_conditions(p, L, bda, bdb, bdc, bdd, LEFT, x_coord);
			fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
			set_neumann_conditions(p, L, bda, bdb, bdc, bdd, LEFT, y_coord);
			fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
			set_neumann_conditions(p, L, bda, bdb, bdc, bdd, LEFT, z_coord);
			fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
			store_arrays(p, fx, fy, fz, storage, LEFT);
			initialize_f(p, fx, fy, fz);
			initialize_bd(p->n, bda);
			initialize_bd(p->n, bdb);
			initialize_bd(p->m, bdc);
			initialize_bd(p->m, bdd);
			//Inner wall
			p->alfa = pi - p->angle;
			p->c = pi / 2;
			p->d = 3 * pi / 2;
			theta_val[RIGHT] = dirichlet_boundary_inner_wall_corrected_v_boundary(p, storage, fx, fy, fz, RIGHT);
			set_neumann_conditions(p, L, bda, bdb, bdc, bdd, RIGHT, x_coord);
			fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
			set_neumann_conditions(p, L, bda, bdb, bdc, bdd, RIGHT, y_coord);
			fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
			set_neumann_conditions(p, L, bda, bdb, bdc, bdd, RIGHT, z_coord);
			fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
			store_arrays(p, fx, fy, fz, storage, RIGHT);
			initialize_f(p, fx, fy, fz);
			initialize_bd(p->n, bda);
			initialize_bd(p->n, bdb);
			initialize_bd(p->m, bdc);
			initialize_bd(p->m, bdd);

			bound_v_correction(p, storage, theta_val);
		}
		/************************************************************/
		p->alfa = pi - p->angle;
		p->c = 3 * pi / 2;
		p->d = 2 * pi + pi / 2;
		for (int i = 0; i < 2; i++) {
			theta_val[i] = (double*) malloc(p->ny * sizeof(double));
		}
		theta_val[LEFT] = dirichlet_boundary_outer_wall_corrected_v_boundary(p, storage, fx, fy, fz, LEFT);
		set_neumann_conditions(p, L, bda, bdb, bdc, bdd, LEFT, x_coord);
		fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

		set_neumann_conditions(p, L, bda, bdb, bdc, bdd, LEFT, y_coord);
		fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

		set_neumann_conditions(p, L, bda, bdb, bdc, bdd, LEFT, z_coord);
		fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
		store_arrays(p, fx, fy, fz, storage, LEFT);
		initialize_f(p, fx, fy, fz);
		initialize_bd(p->n, bda);
		initialize_bd(p->n, bdb);
		initialize_bd(p->m, bdc);
		initialize_bd(p->m, bdd);

		//Inner wall
		p->alfa = pi - p->angle;
		p->c = pi / 2;
		p->d = 3 * pi / 2;
		theta_val[RIGHT] = dirichlet_boundary_inner_wall_corrected_v_boundary(p, storage, fx, fy, fz, RIGHT);
		set_neumann_conditions(p, L, bda, bdb, bdc, bdd, RIGHT, x_coord);
		fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

		set_neumann_conditions(p, L, bda, bdb, bdc, bdd, RIGHT, y_coord);
		fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
		set_neumann_conditions(p, L, bda, bdb, bdc, bdd, RIGHT, z_coord);
		fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
		store_arrays(p, fx, fy, fz, storage, RIGHT);
		initialize_f(p, fx, fy, fz);
		initialize_bd(p->n, bda);
		initialize_bd(p->n, bdb);
		initialize_bd(p->m, bdc);
		initialize_bd(p->m, bdd);
	} else if (branch == R) {
		for (int corrections = 0; corrections < p->max_corrections; corrections++) {
			p->alfa = pi + p->angle;
			p->c = pi / 2;
			p->d = 3 * pi / 2;

			theta_val[LEFT] = dirichlet_boundary_outer_wall_corrected_v_boundary(p, storage, fx, fy, fz, LEFT);
			set_neumann_conditions(p, R, bda, bdb, bdc, bdd, LEFT, x_coord);
			fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
			set_neumann_conditions(p, R, bda, bdb, bdc, bdd, LEFT, y_coord);
			fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
			set_neumann_conditions(p, R, bda, bdb, bdc, bdd, LEFT, z_coord);
			fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

			store_arrays(p, fx, fy, fz, storage, LEFT);
			initialize_f(p, fx, fy, fz);
			initialize_bd(p->n, bda);
			initialize_bd(p->n, bdb);
			initialize_bd(p->m, bdc);
			initialize_bd(p->m, bdd);

			//Inner wall
			p->alfa = pi + p->angle;
			p->c = 3 * pi / 2;
			p->d = 2 * pi + pi / 2;
			theta_val[RIGHT] = dirichlet_boundary_inner_wall_corrected_v_boundary(p, storage, fx, fy, fz, RIGHT);
			set_neumann_conditions(p, R, bda, bdb, bdc, bdd, RIGHT, x_coord);
			fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
			set_neumann_conditions(p, R, bda, bdb, bdc, bdd, RIGHT, y_coord);
			fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
			set_neumann_conditions(p, R, bda, bdb, bdc, bdd, RIGHT, z_coord);
			fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

			store_arrays(p, fx, fy, fz, storage, RIGHT);
			initialize_f(p, fx, fy, fz);
			initialize_bd(p->n, bda);
			initialize_bd(p->n, bdb);
			initialize_bd(p->m, bdc);
			initialize_bd(p->m, bdd);

			bound_v_correction(p, storage, theta_val);
		}
		p->alfa = pi + p->angle;
		p->c = pi / 2;
		p->d = 3 * pi / 2;

		theta_val[LEFT] = dirichlet_boundary_outer_wall_corrected_v_boundary(p, storage, fx, fy, fz, LEFT);
		set_neumann_conditions(p, R, bda, bdb, bdc, bdd, LEFT, x_coord);
		fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
		set_neumann_conditions(p, R, bda, bdb, bdc, bdd, LEFT, y_coord);
		fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
		set_neumann_conditions(p, R, bda, bdb, bdc, bdd, LEFT, z_coord);
		fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

		store_arrays(p, fx, fy, fz, storage, LEFT);
		initialize_f(p, fx, fy, fz);
		initialize_bd(p->n, bda);
		initialize_bd(p->n, bdb);
		initialize_bd(p->m, bdc);
		initialize_bd(p->m, bdd);

		//Inner wall
		p->alfa = pi + p->angle;
		p->c = 3 * pi / 2;
		p->d = 2 * pi + pi / 2;
		theta_val[RIGHT] = dirichlet_boundary_inner_wall_corrected_v_boundary(p, storage, fx, fy, fz, RIGHT);
		set_neumann_conditions(p, R, bda, bdb, bdc, bdd, RIGHT, x_coord);
		fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
		set_neumann_conditions(p, R, bda, bdb, bdc, bdd, RIGHT, y_coord);
		fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
		set_neumann_conditions(p, R, bda, bdb, bdc, bdd, RIGHT, z_coord);
		fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

		store_arrays(p, fx, fy, fz, storage, RIGHT);
		initialize_f(p, fx, fy, fz);
		initialize_bd(p->n, bda);
		initialize_bd(p->n, bdb);
		initialize_bd(p->m, bdc);
		initialize_bd(p->m, bdd);
	}

}

void set_neumann_conditions(parms* p, int branch, double* bda, double* bdb, double* bdc, double* bdd, int half, int coordinate) {
	int scale = 1;
	/// If branch is parent
	if (branch == P) {
		if (half == LEFT) {
			if (coordinate == x_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = p->s0;
					bdb[k] = 0.0;
				}
			} else if (coordinate == y_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = -p->h1;
					bdb[k] = p->h1;
				}
			} else if (coordinate == z_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;
					bdb[k] = 0;
				}

			}
		} else if (half == RIGHT) {
			if (coordinate == x_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = -p->s0;
					bdb[k] = 0.0;
				}
			} else if (coordinate == y_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = -p->h1;
					bdb[k] = p->h1;
				}
			} else if (coordinate == z_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;
					bdb[k] = 0;
				}
			}
		}
	}
	/// If branch is Left daughter
	else if (branch == L) {
		if (half == LEFT) {
			if (coordinate == x_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = p->s0 * cos(p->angle);
					bdb[k] = 0;
				}
			} else if (coordinate == y_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = p->s0 * sin(p->angle);
					bdb[k] = -3.5 * p->s1;	// * cos(p->angle);
				}
			} else if (coordinate == z_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;
					bdb[k] = 0;
				}
			}
		} else if (half == RIGHT) {
			if (coordinate == x_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;
					bdb[k] = p->apex_scaling * p->s1;	// * sin(p->angle);
				}
			} else if (coordinate == y_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = p->s0 * cos(p->alfa);
					bdb[k] = p->s1;
				}
			} else if (coordinate == z_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;
					bdb[k] = 0;
				}
			}
		}
	}
	/// If branch is Right daughter
	else if (branch == R) {
		if (half == LEFT) {
			if (coordinate == x_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;	//* cos(p->angle);
					bdb[k] = 0;	//p->s0;
				}
			} else if (coordinate == y_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;	//p->s1;	//* sin(p->angle);
					bdb[k] = -3 * p->s1;	//p->s0;
				}
			} else if (coordinate == z_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;
					bdb[k] = 0;
				}
			}
		} else if (half == RIGHT) {
			if (coordinate == x_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = -p->s0 * sin(p->alfa);	//-p->s1 * cos(p->angle);
					bdb[k] = -p->apex_scaling * p->s1;	// * sin(p->angle);
				}
			} else if (coordinate == y_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = p->s0 * cos(p->alfa);
					bdb[k] = p->s1;	//p->s0;	//-p->s0 * cos(p->angle);
				}
			} else if (coordinate == z_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0;
					bdb[k] = 0;
				}
			}
		}
	}
	/// If branch is Endcaps
	else if (branch == endcaps) {
		if (half == LEFT) {
			if (coordinate == x_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0.0;
					bdb[k] = 0.0;
				}
			} else if (coordinate == y_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0.0;
					bdb[k] = 0.0;
				}
			} else if (coordinate == z_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0.0;
					bdb[k] = 0.0;
				}
			}
		} else if (half == RIGHT) {
			if (coordinate == x_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0.0;
					bdb[k] = 0.0;
				}
			} else if (coordinate == y_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0.0;
					bdb[k] = 0.0;
				}
			} else if (coordinate == z_coord) {
				for (int k = 0; k < p->n; k++) {
					bda[k] = 0.0;
					bdb[k] = 0.0;
				}
			}
		}
	}
}
