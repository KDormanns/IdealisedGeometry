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

	p.rb = 1.65;
	p.rt = 1.65;
	p.h1 = -16.575;
	p.h2 = 0;
	p.h3 = 16.575;
	p.ra = p.rt;

	p.xshift = 0.0;
	p.c1 = 1.0;
	p.c2 = 1.5;
	p.b_bnw = p.rb;
	p.angle = pi / 1.8;
	p.ss0 = 3;
	p.ss1 = 3;

	p.xc = 0.8;

	p.r = 0.95; // radius of intersecting cylinder where 0.636 <= r <=0.9886

	p.s0 = 1.5;
	p.s1 = 1.0;

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

	p.m = 13;
	p.n = 5;
	int parent_grid_point_skip = 3, daughter_grid_point_skip = 3;
	int downstream_offset, upstream_offset;

	p.nx = p.m + 2;
	p.ny = p.n + 2;
	int num_SMCs = 0, num_ECs = 0;
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
	double w[lw];
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

	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);
	double ****storage = allocate_storage_array(p.nx, p.ny);
	int* info = (int*) malloc(2 * sizeof(int));

	double* total_points = (double*) malloc(3 * info[0] * sizeof(double));
	int* total_cells = (int*) malloc(3 * info[1] * sizeof(int));

	/***** Solving for parent segment ******/

	///Setting up conditions for left half of the geometry
	p.alfa = pi - p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;

	dirichlet_boundary_parent_segment(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.h1;
		bdb[k] = p.h1;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(p, fx, fy, fz, storage, LEFT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	//Right of parent segment

	p.alfa = pi + p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;

	dirichlet_boundary_parent_segment(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.h1;
		bdb[k] = p.h1;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(p, fx, fy, fz, storage, RIGHT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = 0;
	upstream_offset = parent_grid_point_skip;
	bound_v_correction(p, storage, downstream_offset, upstream_offset);
	info = format_primitive(p, storage, "parent", downstream_offset, upstream_offset);

	/***** Solving for Left daughter segment ******/

	//Outer wall
	p.alfa = pi - p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;

	dirichlet_boundary_outer_wall(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = p.s1 * cos(p.angle);
		bdb[k] = 0;
	}

	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = p.s1 * sin(p.angle);
		bdb[k] = -p.s0;
	}

	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(p, fx, fy, fz, storage, LEFT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	//Inner wall
	p.alfa = pi - p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	dirichlet_boundary_inner_wall(p, fx, fy, fz);
	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s1 * cos(p.angle);
		bdb[k] = p.s0 * sin(p.angle);
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s1 * sin(p.angle);
		bdb[k] = -p.s0 * cos(p.angle);
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(p, fx, fy, fz, storage, RIGHT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = daughter_grid_point_skip;
	upstream_offset = 0;
	bound_v_correction(p, storage, downstream_offset, upstream_offset);
	format_primitive(p, storage, "left_daughter", downstream_offset, upstream_offset);

	/***** Solving for Right daughter segment ******/

	//Outer wall
	p.alfa = pi + p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;

	dirichlet_boundary_outer_wall(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = p.s1 * cos(p.angle);
		bdb[k] = 0;
	}

	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = p.s1 * sin(p.angle);
		bdb[k] = -p.s0;
	}

	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(p, fx, fy, fz, storage, LEFT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	//Inner wall
	p.alfa = pi + p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;
	dirichlet_boundary_inner_wall(p, fx, fy, fz);
	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s1 * cos(p.angle);
		bdb[k] = -p.s0 * sin(p.angle);
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s1 * sin(p.angle);
		bdb[k] = -p.s0 * cos(p.angle);
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	store_arrays(p, fx, fy, fz, storage, RIGHT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = daughter_grid_point_skip;
	upstream_offset = 0;
	bound_v_correction(p, storage, downstream_offset, upstream_offset);
	format_primitive(p, storage, "right_daughter", downstream_offset, upstream_offset);

	/******* Forming end_caps on inlet and outlets ********/
	//Parent segment inlet endcap
	p.alfa = pi - p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;
	dirichlet_boundary_end_cap_parent(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(p, fx, fy, fz, storage, LEFT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	p.alfa = pi - p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	dirichlet_boundary_end_cap_parent(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(p, fx, fy, fz, storage, RIGHT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = 0;
	upstream_offset = 0;
	format_primitive(p, storage, "endcap_parent", downstream_offset, upstream_offset);

	//Left daughter segment inlet endcap
	p.alfa = pi - p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(p, fx, fy, fz, storage, LEFT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	p.alfa = pi - p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(p, fx, fy, fz, storage, RIGHT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = 0;
	upstream_offset = 0;
	format_primitive(p, storage, "endcap_left_daughter", downstream_offset, upstream_offset);

	//Right daughter segment inlet endcap
	p.alfa = pi + p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(p, fx, fy, fz, storage, LEFT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	p.alfa = pi + p.angle;
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);
	store_arrays(p, fx, fy, fz, storage, RIGHT);
	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);

	downstream_offset = 0;
	upstream_offset = 0;
	format_primitive(p, storage, "endcap_right_daughter", downstream_offset, upstream_offset);

}

