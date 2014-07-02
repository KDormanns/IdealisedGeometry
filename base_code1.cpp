#include "functionality.h"
parms p;

int main(int argc, char* argv[]) {

	///All entities are in mm dimensions. (the scaling is done while writing the results in vtk and stl files)
	p.alpha = 0.0;
	p.beta = 0.0;

	p.a = 0.0;
	p.b = 1.0;
	p.c = 0.0;
	p.d = 0.0;

	p.rb = 2.5;
	p.length = 10 * p.rb;
	p.h1 = -p.length;
	p.h2 = 0;
	p.h3 = p.length;
	p.max_corrections = default_max_corrections;
	p.ECs = default_NUM_ECs;
	p.SMCs = default_NUM_SMCs;
	p.m = 11;
	p.n = 5;
	int parent_grid_point_skip = 0, daughter_grid_point_skip = 0;
	int downstream_offset, upstream_offset;

	if (argc > 1) {
		for (int i = 0; i < argc; i++) {
			if (argv[i][0] == '-') {
				if (argv[i][1] == 'c') {
					p.max_corrections = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'E') { // Number of ECs in axial direction. E.g. 4.
					p.ECs = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'S') { // Number of SMCs in circumferential direction. E.g. 4.
					p.SMCs = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'm') { // Number of grid points in axial direction. Must be odd. E.g. 33.
					p.m = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'n') { // Number of grid points in circumferential direction for a half-circle. Must be odd. E.g. 15.
					p.n = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'D') {      // Downstream skip. E.g. 0.
					parent_grid_point_skip = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'U') {        // Upstream skip. E.g. 0.
					daughter_grid_point_skip = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'l') { // Length of a branch in mm. Must be langer than radius r. E.g. 25.
					p.h1 = atof(argv[i + 1]) * (-1);
					p.h3 = atof(argv[i + 1]);
				} else if (argv[i][1] == 'r') {           // Raduis r. E.g. 2.5.
					p.rb = atof(argv[i + 1]);
				}
			}
		}
	}

	p.rt = p.rb * 0.6;

	p.neck = p.rb * 1.2;
	p.apex = p.neck * 1.5;

	p.ra = p.rt;

	p.xshift = 0.0;
	p.c1 = 1.0;
	p.c2 = 1.5;
	p.b_bnw = p.rb;
	p.angle = pi / 1.8;
	p.ss0 = 1;
	p.ss1 = 1;

	p.xc = 0.8;

	p.r = 0.95; // radius of intersecting cylinder where 0.636 <= r <=0.9886

	p.s0 = 1;
	p.s1 = 1;
	p.apex_scaling = 7.5;

	p.A1 = p.c1;
	p.B1 = 0;
	p.C1 = 3 * (p.c2 * cos(p.alfa) - p.c1) - (p.ss1 * sin(p.alfa));
	p.D1 = 2 * (p.c2 * cos(p.alfa) - p.c1) + (p.ss1 * sin(p.alfa));

	p.A2 = 0;
	p.B2 = 3;
	p.C2 = (3 * p.c2 * sin(p.alfa)) - (2 * p.ss0) - (p.ss1 * cos(p.alfa));
	p.D2 = (p.ss1 * cos(p.alfa)) + (p.ss0) - (2 * p.c2 * sin(p.alfa));

	p.v1 = acos(p.xc / p.c1);
	p.v2 = acos(-(p.xc - p.r * sin(p.alfa)) / (p.c2 * cos(p.alfa)));

	p.nx = p.m + 2;
	p.ny = p.n + 2;
	cout << p.max_corrections << "\t" << p.ECs << "\t" << p.SMCs << "\t" << p.m << "\t" << p.n << "\t" << parent_grid_point_skip << "\t"
			<< daughter_grid_point_skip << endl;
	cout << "length per branch = " << p.h1 << "\t" << "radius=" << p.rb << endl;
	double tol = 1e-7;
	int iflag = 4;
	int itcg = 50;
	int mm, nn, lw;
	if (iflag == 2) {
		mm = 3 * p.m;
		nn = 7 * p.n;
		lw = (int) (fmax(mm, nn)) + 2 * (p.n + p.m) + 19; //for iflag = 2
	} else if (iflag == 4) {
		mm = 3 * p.m;
		nn = 4 * p.n;
		lw = (int) (fmax(mm, nn)) + 4 * p.n + 2 * p.m + 0.5 * ((p.n + 1) * (p.n + 1)) + 19;	//for iflag = 4
	}
	lw = 10000;
	double *w = (double*) malloc(lw * sizeof(double));
	int idf = p.nx;

	double **fx, **fy, **fz;
	double *bda, *bdb, *bdc, *bdd;

	fx = allocate_dirichlet(p.nx, p.ny);
	fy = allocate_dirichlet(p.nx, p.ny);
	fz = allocate_dirichlet(p.nx, p.ny);
	bda = allocate_neuman(p.n);
	bdb = allocate_neuman(p.n);
	bdc = allocate_neuman(p.m);
	bdd = allocate_neuman(p.m);

	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);
	double ****storage = allocate_storage_array(p.nx, p.ny);
	double** theta_val = (double**) malloc(2 * sizeof(double*));
	for (int i = 0; i < 2; i++) {
		theta_val[i] = (double*) malloc(p.ny * sizeof(double));
	}
	int* info = (int*) malloc(2 * sizeof(int));

	double* total_points = (double*) malloc(3 * info[0] * sizeof(double));
	int* total_cells = (int*) malloc(3 * info[1] * sizeof(int));

	/***** Solving for parent segment ******/
	printf("Solving for parent segment\n");

	///Setting up conditions for left half of the geometry
	printf("Processing left half\n");

	p.alfa = pi - p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;

	dirichlet_boundary_parent_segment(p, fx, fy, fz);

	set_neumann_conditions(&p, P, bda, bdb, bdc, bdd, LEFT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	set_neumann_conditions(&p, P, bda, bdb, bdc, bdd, LEFT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	set_neumann_conditions(&p, P, bda, bdb, bdc, bdd, LEFT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(&p, fx, fy, fz, storage, LEFT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	//Right of parent segment
	printf("Processing right half\n");
	p.alfa = pi + p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;

	dirichlet_boundary_parent_segment(p, fx, fy, fz);

	set_neumann_conditions(&p, P, bda, bdb, bdc, bdd, RIGHT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	set_neumann_conditions(&p, P, bda, bdb, bdc, bdd, RIGHT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	set_neumann_conditions(&p, P, bda, bdb, bdc, bdd, RIGHT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(&p, fx, fy, fz, storage, RIGHT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = 0;
	upstream_offset = parent_grid_point_skip;
	info = format_primitive(p, storage, "parent", downstream_offset, upstream_offset);
	/******************************************************************************************************/
	/***** Solving for Left daughter segment ******/
	printf("\nSolving for Left daughter segment\n");
	//Outer wall
	printf("Processing Outer wall of Left daughter segment\n");
	p.alfa = pi - p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;

	theta_val[LEFT] = dirichlet_boundary_outer_wall(p, fx, fy, fz);

	set_neumann_conditions(&p, L, bda, bdb, bdc, bdd, LEFT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, L, bda, bdb, bdc, bdd, LEFT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, L, bda, bdb, bdc, bdd, LEFT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(&p, fx, fy, fz, storage, LEFT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	//Inner wall
	printf("Inner wall of Left daughter segment\n");
	p.alfa = pi - p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	theta_val[RIGHT] = dirichlet_boundary_inner_wall(p, fx, fy, fz);

	set_neumann_conditions(&p, L, bda, bdb, bdc, bdd, RIGHT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, L, bda, bdb, bdc, bdd, RIGHT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, L, bda, bdb, bdc, bdd, RIGHT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(&p, fx, fy, fz, storage, RIGHT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = daughter_grid_point_skip;
	upstream_offset = 0;
	itereate_to_correct(L, &p, fx, fy, fz, storage, bda, bdb, bdc, bdd, idf, iflag, tol, itcg, w, lw, theta_val);
	format_primitive(p, storage, "left_daughter", downstream_offset, upstream_offset);

	/***** Solving for Right daughter segment ******/
	printf("Solving for Right daughter segment\n");
	//Outer wall
	printf("Outer wall of Right daughter segment\n");
	p.alfa = pi + p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;

	dirichlet_boundary_outer_wall(p, fx, fy, fz);

	set_neumann_conditions(&p, R, bda, bdb, bdc, bdd, LEFT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, R, bda, bdb, bdc, bdd, LEFT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, R, bda, bdb, bdc, bdd, LEFT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(&p, fx, fy, fz, storage, LEFT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	//Inner wall
	printf("Inner wall of Left daughter segment\n");
	p.alfa = pi + p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;
	dirichlet_boundary_inner_wall(p, fx, fy, fz);
	set_neumann_conditions(&p, R, bda, bdb, bdc, bdd, RIGHT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, R, bda, bdb, bdc, bdd, RIGHT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, R, bda, bdb, bdc, bdd, RIGHT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(&p, fx, fy, fz, storage, RIGHT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = daughter_grid_point_skip;
	upstream_offset = 0;
	itereate_to_correct(R, &p, fx, fy, fz, storage, bda, bdb, bdc, bdd, idf, iflag, tol, itcg, w, lw, theta_val);
	format_primitive(p, storage, "right_daughter", downstream_offset, upstream_offset);
	/******************************************************************************************************/
	/******* Forming end_caps on inlet and outlets ********/
	//Parent segment inlet endcap
	printf("Solving for Parent segment inlet endcap\n");

	printf("Processing left half of Parent segment inlet endcap\n");
	p.alfa = pi - p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;
	dirichlet_boundary_end_cap_parent(p, fx, fy, fz);

	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, z_coord);

	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(&p, fx, fy, fz, storage, LEFT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	printf("Processing right half of Parent segment inlet endcap\n");
	p.alfa = pi - p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	dirichlet_boundary_end_cap_parent(p, fx, fy, fz);

	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(&p, fx, fy, fz, storage, RIGHT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = 0;
	upstream_offset = 0;
	format_primitive(p, storage, "endcap_parent", downstream_offset, upstream_offset);

	//Left daughter segment inlet endcap
	printf("Processing left half of Left daughter segment inlet endcap\n");
	p.alfa = pi - p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(&p, fx, fy, fz, storage, LEFT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	printf("Processing right half of Left daughter segment inlet endcap\n");
	p.alfa = pi - p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(&p, fx, fy, fz, storage, RIGHT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = 0;
	upstream_offset = 0;
	format_primitive(p, storage, "endcap_left_daughter", downstream_offset, upstream_offset);

	//Right daughter segment inlet endcap
	printf("Processing left half of Right daughter segment inlet endcap\n");
	p.alfa = pi + p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, LEFT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(&p, fx, fy, fz, storage, LEFT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	printf("Processing right half of Right daughter segment inlet endcap\n");
	p.alfa = pi + p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, x_coord);
	fx = solve(&p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, y_coord);
	fy = solve(&p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	set_neumann_conditions(&p, endcaps, bda, bdb, bdc, bdd, RIGHT, z_coord);
	fz = solve(&p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(&p, fx, fy, fz, storage, RIGHT);
	initialize_f(&p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = 0;
	upstream_offset = 0;
	format_primitive(p, storage, "endcap_right_daughter", downstream_offset, upstream_offset);

	FILE * fw;
	char filename[50];
	sprintf(filename, "configuration_info.txt");
	fw = fopen(filename, "w+");

	fprintf(fw, "Processors information\n");
	fprintf(fw, "Total number of points per branch (vtk points)= %d\t m = %d n= %d\n", info[0], info[1], info[2]);
	fprintf(fw, "Total number of cells per branch (vtk cells)= %d\t m = %d n= %d\n", info[3], info[4], info[5]);

	fprintf(fw, "Total number of SMC mesh points per processor mesh (vtk points)= %d\t m = %d n= %d\n", info[6], info[7], info[8]);
	fprintf(fw, "Total number of SMC mesh cells per processor mesh (vtk points)= %d\t m = %d n= %d\n", info[9], info[10], info[11]);

	fprintf(fw, "Total number of EC mesh points per processor mesh (vtk points)= %d\t m = %d n= %d\n", info[12], info[13], info[14]);
	fprintf(fw, "Total number of EC mesh cells per processor mesh (vtk points)= %d\t m = %d n= %d\n", info[15], info[16], info[17]);

	fprintf(fw, "Total number of EC mesh centeroid points per processor mesh (vtk points)= %d\t m = %d n= %d\n", info[18], info[16], info[17]);
	fprintf(fw, "Total number of EC mesh centeroid cells per processor mesh (vtk points)= %d\t m = %d n= %d\n", info[18], info[16], info[17]);

	fclose(fw);
}

