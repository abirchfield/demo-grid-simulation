#include "example1.h"

/* TINY EXAMPLE 2 - SINGLE GENROU MACHINE, INFINITE BUS in pure IDA */

#define CHECK_OUTPUT(F) retval=F; if (retval < 0) { \
     printf("ERROR: Flag %d in function " #F, retval); return 1; }

#define CHECK_ALLOCATION(A, F) A = F; if (A == NULL) { \
    printf("ERROR: Bad Alloc from function " #F); return 1; }

typedef struct 
{
    double omega0, H, D, Ra, Tdop, Tdopp, Tqopp, Tqop, Xd, Xdp, Xdpp, Xq, Xqp,
        Xqpp, Xl, S10, S12, SA, SB, Xd1, Xd2, Xd3, Xd4, Xd5, Xq1, Xq2, Xq3,
        Xq4, Xq5, Xqd, g, b, branch_x, fault_b, pmech_set, efd_set,
        vr_inf_set, vi_inf_set;
} ex1_user_data;

int ex1_residual(sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector rr, 
    void *user_data)
{
    sunrealtype *yy_ = NV_DATA_S(yy);
    sunrealtype *yp_ = NV_DATA_S(yp);
    sunrealtype *rr_ = NV_DATA_S(rr);
    ex1_user_data* u = (ex1_user_data*) user_data;
    double delta, omega, Eqp, psidp, psiqp, Edp, psiqpp, psidpp, psipp, ksat,
        vd, vq, telec, id, iq, ir, ii, pmech, efd, inr, ini, vr, vi, vr_inf, 
        vi_inf, delta_dot, omega_dot, Eqp_dot, psidp_dot, psiqp_dot, Edp_dot,
        g, b;

    /* Read variables */
    delta = yy_[0];
    omega = yy_[1];
    Eqp = yy_[2];
    psidp = yy_[3];
    psiqp = yy_[4];
    Edp = yy_[5];
    psiqpp = yy_[6]; 
    psidpp = yy_[7]; 
    psipp = yy_[8]; 
    ksat = yy_[9];
    vd = yy_[10];
    vq = yy_[11];
    telec = yy_[12];
    id = yy_[13];
    iq = yy_[14];
    ir = yy_[15];
    ii = yy_[16];
    pmech = yy_[17];
    efd = yy_[18];
    inr = yy_[19];
    ini = yy_[20];
    vr = yy_[21]; 
    vi = yy_[22]; 
    vr_inf = yy_[23];
    vi_inf = yy_[24];

    /* Read derivatives */
    delta_dot = yp_[0];
    omega_dot = yp_[1];
    Eqp_dot = yp_[2];
    psidp_dot = yp_[3];
    psiqp_dot = yp_[4];
    Edp_dot = yp_[5];

    /* 6 GENROU differential equations */
    rr_[0] = delta_dot - omega*u->omega0;
    rr_[1] = omega_dot - (1/(2*u->H)) * ((pmech - u->D*omega) / (1 + omega) 
        - telec);
    rr_[2] = Eqp_dot - (1/u->Tdop) * (efd - (Eqp + u->Xd1*(id + u->Xd3*(Eqp 
        - psidp - u->Xd2*id)) + psidpp*ksat ));
    rr_[3] = psidp_dot - (1/u->Tdopp) * (Eqp - psidp - u->Xd2*id);
    rr_[4] = psiqp_dot - (1/u->Tqopp) * (Edp - psiqp + u->Xq2*iq);
    rr_[5] = Edp_dot - (1/u->Tqop) * (-Edp + u->Xqd*psiqpp*ksat 
        + u->Xq1*(iq - u->Xq3*(Edp + iq*u->Xq2 - psiqp)));
    
    /* 11 GENROU algebraic equations */
    rr_[6] = psiqpp - (-psiqp * u->Xq4 - Edp * u->Xq5);
    rr_[7] = psidpp - (psidp * u->Xd4 + Eqp * u->Xd5);
    rr_[8] = psipp - sqrt(pow(psidpp, 2.0) + pow(psiqpp, 2.0));
    rr_[9] = ksat - u->SB*pow(psipp - u->SA, 2.0);
    rr_[10] = vd + psiqpp * (1 + omega);
    rr_[11] = vq - psidpp * (1 + omega);
    rr_[12] = telec - ((psidpp - id*u->Xdpp)*iq - (psiqpp - iq*u->Xdpp)*id);
    rr_[13] = id - (ir*sin(delta) - ii*cos(delta));
    rr_[14] = iq - (ir*cos(delta) + ii*sin(delta));
    rr_[15] = ir + u->g*vr - u->b*vi - inr;
    rr_[16] = ii + u->b*vr + u->g*vi - ini;
    
    /* 2 GENROU control inputs are set to constant for this example */
    rr_[17] = pmech - u->pmech_set;
    rr_[18] = efd - u->efd_set;

    /* 2 GENROU current source definitions */
    rr_[19] = inr - (u->g*(sin(delta)*vd + cos(delta)*vq) 
        - u->b*(-cos(delta)*vd + sin(delta)*vq));
    rr_[20] = ini - (u->b*(sin(delta)*vd + cos(delta)*vq) 
        + u->g*(-cos(delta)*vd + sin(delta)*vq));

    /* Bus 1 network constraints (current balance) */
    g = u->g;
    double bbr = - 1/u->branch_x;
    b = u->b + bbr + u->fault_b;
    rr_[21] = -inr + vr*g - vi*b + vi_inf*bbr;
    rr_[22] = -ini + vr*b + vi*g - vr_inf*bbr;

    /* Bus 2 network constraints (fixed infinite bus voltage) */
    rr_[23] = vr_inf - u->vr_inf_set;
    rr_[24] = vi_inf - u->vi_inf_set;

    return 0;
}

#define J_ADD(i, j, a) Ti[nz] = i; Tj[nz] = j; Tx[nz] = a; nz++;

int ex1_jacobian(sunrealtype t, sunrealtype cj, N_Vector yy, N_Vector yp,
    N_Vector resvec, SUNMatrix J, void *user_data, N_Vector tmp1,
    N_Vector tmp2, N_Vector tmp3)
{
    sunrealtype *yy_ = NV_DATA_S(yy);
    sunrealtype *yp_ = NV_DATA_S(yp);
    ex1_user_data* u = (ex1_user_data*) user_data;
    double delta, omega, Eqp, psidp, psiqp, Edp, psiqpp, psidpp, psipp, ksat,
        vd, vq, telec, id, iq, ir, ii, pmech, efd, inr, ini, vr, vi, vr_inf, 
        vi_inf, delta_dot, omega_dot, Eqp_dot, psidp_dot, psiqp_dot, Edp_dot,
        g, b;

    /* Read variables */
    delta = yy_[0];
    omega = yy_[1];
    Eqp = yy_[2];
    psidp = yy_[3];
    psiqp = yy_[4];
    Edp = yy_[5];
    psiqpp = yy_[6]; 
    psidpp = yy_[7]; 
    psipp = yy_[8]; 
    ksat = yy_[9];
    vd = yy_[10];
    vq = yy_[11];
    telec = yy_[12];
    id = yy_[13];
    iq = yy_[14];
    ir = yy_[15];
    ii = yy_[16];
    pmech = yy_[17];
    efd = yy_[18];
    inr = yy_[19];
    ini = yy_[20];
    vr = yy_[21]; 
    vi = yy_[22]; 
    vr_inf = yy_[23];
    vi_inf = yy_[24];

    /* Read derivatives */
    delta_dot = yp_[0];
    omega_dot = yp_[1];
    Eqp_dot = yp_[2];
    psidp_dot = yp_[3];
    psiqp_dot = yp_[4];
    Edp_dot = yp_[5];
    
    g = u->g;
    double bbr = - 1/u->branch_x;
    b = u->b + bbr + u->fault_b;

    /* Set up triplet structure for entry of matrix values */
    const int nzmax = 200;
    int n = SUNSparseMatrix_Columns(J), nz = 0;
    int Ti[nzmax];
    int Tj[nzmax];
    double Tx[nzmax];


    /* Equation 1: delta_dot */
    J_ADD(0, 0, cj);
    J_ADD(0, 1, -u->omega0);

    /* Equation 2: omega_dot */
    J_ADD(1, 1, cj - (1/(2*u->H)) * (1 + omega + u->D*(-1 + omega) 
        - pmech) / pow(omega + 1, 2) );
    J_ADD(1, 12, -1/(2*u->H));
    J_ADD(1, 17, - (1/(2*u->H)) / (1 + omega));

    /* Equation 3: Eqp_dot */
    J_ADD(2, 2, cj - (1/u->Tdop) * (-1 -u->Xd1*u->Xd3));
    J_ADD(2, 3,  - (1/u->Tdop) * u->Xd1*u->Xd3*psidp);
    J_ADD(2, 7, - (1/u->Tdop) * ksat);
    J_ADD(2, 9, - (1/u->Tdop) * psidpp);
    J_ADD(2, 13, - (1/u->Tdop) * (-u->Xd1 + u->Xd1*u->Xd3*u->Xd2));
    J_ADD(2, 18, - (1/u->Tdop));

    /* Equation 4: psidp_dot */
    J_ADD(3, 2, - 1/u->Tdopp);
    J_ADD(3, 3, cj + 1/u->Tdopp);
    J_ADD(3, 13, - 1/u->Tdopp * -u->Xd2);

    /* Equation 5: psiqp_dot */
    J_ADD(4, 5, - 1/u->Tqopp);
    J_ADD(4, 4, cj + 1/u->Tqopp);
    J_ADD(4, 14, - 1/u->Tqopp*u->Xq2);

    /* Equation 6: Edp_dot */
    J_ADD(5, 4, - 1/u->Tqop * u->Xq1*u->Xq3*psiqp);
    J_ADD(5, 5, cj - 1/u->Tqop * (-1 -u->Xq1*u->Xq3));
    J_ADD(5, 6, - 1/u->Tqop * u->Xqd*ksat);
    J_ADD(5, 9, - 1/u->Tqop * u->Xqd*psiqpp);
    J_ADD(5, 14, - 1/u->Tqop * (u->Xq1 - u->Xq1*u->Xq3*u->Xq2));
    
    /* Equation 7: psiqpp */
    J_ADD(6, 4, u->Xq4);
    J_ADD(6, 5, u->Xq5);
    J_ADD(6, 6, 1);

    /* Equation 8: psidpp */
    J_ADD(7, 2, -u->Xd5);
    J_ADD(7, 3, -u->Xd4);
    J_ADD(7, 7, 1);

    /* Equation 9: psipp */
    J_ADD(8, 6, psiqpp / sqrt(pow(psidpp, 2.0) + pow(psiqpp, 2.0)));
    J_ADD(8, 7, psidpp / sqrt(pow(psidpp, 2.0) + pow(psiqpp, 2.0)));
    J_ADD(8, 8, 1);

    /* Equation 10: ksat */
    J_ADD(9, 9, 1);
    J_ADD(9, 8, -u->SB*2*(psipp -u->SA));

    /* Equation 11: vd */
    J_ADD(10, 1, psiqpp);
    J_ADD(10, 6, 1+omega);
    J_ADD(10, 10, 1);

    /* Equation 12: vq */
    J_ADD(11, 11, 1);
    J_ADD(11, 7, -(1 + omega));
    J_ADD(11, 1, -psidpp);

    /* Equation 13: telec */
    J_ADD(12, 6, id);
    J_ADD(12, 7, -iq);
    J_ADD(12, 12, 1);
    J_ADD(12, 13, psiqpp);
    J_ADD(12, 14, -psidpp);

    /* Equation 14: id */
    J_ADD(13, 0, -ir*cos(delta) -ii*sin(delta));
    J_ADD(13, 13, 1);
    J_ADD(13, 15, -sin(delta));
    J_ADD(13, 16, cos(delta));

    /* Equation 15: iq */
    J_ADD(14, 0, ir*sin(delta) - ii*cos(delta));
    J_ADD(14, 14, 1);
    J_ADD(14, 15, -cos(delta));
    J_ADD(14, 16, -sin(delta));

    /* Equation 16: ir */
    J_ADD(15, 15, 1);
    J_ADD(15, 21, u->g);
    J_ADD(15, 22, -u->b);
    J_ADD(15, 19, -1);

    /* Equation 17: ii */
    J_ADD(16, 16, 1);
    J_ADD(16, 21, u->b);
    J_ADD(16, 22, u->g);
    J_ADD(16, 20, -1);

    /* Equation 18: pmech */
    J_ADD(17, 17, 1);

    /* Equation 19: efd */
    J_ADD(18, 18, 1);

    /* Equation 20: inr */
    J_ADD(19, 19, 1);
    J_ADD(19, 10, -u->g*sin(delta) - u->b*cos(delta));
    J_ADD(19, 11, -u->g*cos(delta) + u->b*sin(delta));
    J_ADD(19, 0, - (u->g*(cos(delta)*vd + -sin(delta)*vq) 
        - u->b*(sin(delta)*vd + cos(delta)*vq)));

    /* Equation 21: ini */
    J_ADD(20, 20, 1);
    J_ADD(20, 10, -u->b*sin(delta) - u->g*cos(delta));
    J_ADD(20, 11, -u->b*cos(delta) + u->g*sin(delta));
    J_ADD(20, 0, -(u->b*(cos(delta)*vd - sin(delta)*vq)
        + u->g*(sin(delta)*vd + cos(delta)*vq)));

    /* Equation 22: real current balance at bus 1 */
    J_ADD(21, 19, -1);
    J_ADD(21, 21, g);
    J_ADD(21, 22, -b);
    J_ADD(21, 24, bbr);

    /* Equation 23: imag current balance at bus 1 */
    J_ADD(22, 20, -1);
    J_ADD(22, 21, b);
    J_ADD(22, 22, g);
    J_ADD(22, 23, -bbr);

    /* Equation 24: bus 2 network constaint (infinite bus) real */
    J_ADD(23, 23, 1);

    /* Equation 25: bus 2 network constaint (infinite bus) imag */
    J_ADD(24, 24, 1);

    /* Compress matrix into CSC format */
    compress_matrix(J, n, nz, Ti, Tj, Tx);
    if (0)
    {
        FILE * f = fopen("out.txt", "w");
        SUNSparseMatrix_Print(J, f);
        fclose(f);
    }
    return 0;
}

int example1()
{
    printf("Example 1\n");

    SUNContext context;
    N_Vector yy, yp, yy0, yp0, yid, rr;
    void *ida_mem;
    SUNMatrix J;
    SUNLinearSolver linsol;

    int retval;
    int ny = 25, nnz = 200;
    ex1_user_data u;
    double rtol = 1e-7, atol = 1e-9;
    long maxsteps = 2000;

    CHECK_OUTPUT( SUNContext_Create(SUN_COMM_NULL, &context) );
    CHECK_ALLOCATION( ida_mem, IDACreate(context) );
    CHECK_ALLOCATION( yy, N_VNew_Serial(ny, context) );
    CHECK_ALLOCATION( yp, N_VClone(yy) );
    CHECK_ALLOCATION( rr, N_VClone(yy) );


    /* Set given parameters */
    u.omega0 = 2*M_PI*60;
    u.H = 3;
    u.D = 0;
    u.Ra = 0;
    u.Tdop = 7;
    u.Tdopp = 0.04;
    u.Tqopp = 0.05;
    u.Tqop = 0.75;
    u.Xd = 2.1;
    u.Xdp = 0.2;
    u.Xdpp = 0.18;
    u.Xq = 0.5;
    u.Xqp = 0.5;
    u.Xqpp = 0.18;
    u.Xl = 0.15;
    u.S10 = 0;
    u.S12 = 0;

    /* Set calculated constants from given parameters */
    u.SA = 0;
    u.SB = 0;
    if (u.S12 != 0)
    {
		double s112 = sqrt(u.S10 / u.S12);
		u.SA = (1.2*s112 + 1) / (s112 + 1);
		u.SB = (1.2*s112 - 1) / (s112 - 1);
		if (u.SB < u.SA) u.SA = u.SB;
		u.SB = u.S12 / pow(u.SA - 1.2, 2);
    }
    u.Xd1 = u.Xd - u.Xdp;
    u.Xd2 = u.Xdp - u.Xl;
    u.Xd3 = (u.Xdp - u.Xdpp) / (u.Xd2 * u.Xd2);
    u.Xd4 = (u.Xdp - u.Xdpp) / u.Xd2;
    u.Xd5 = (u.Xdpp - u.Xl) / u.Xd2;
    u.Xq1 = u.Xq - u.Xqp;
    u.Xq2 = u.Xqp - u.Xl;
    u.Xq3 = (u.Xqp - u.Xqpp) / (u.Xq2 * u.Xq2);
    u.Xq4 = (u.Xqp - u.Xqpp) / u.Xq2;
    u.Xq5 = (u.Xqpp - u.Xl) / u.Xq2;
    u.Xqd = (u.Xq - u.Xl) / (u.Xd - u.Xl);
    u.g = u.Ra / (u.Ra*u.Ra + u.Xqpp*u.Xqpp);
    u.b = -u.Xqpp / (u.Ra*u.Ra + u.Xqpp*u.Xqpp);
    
    /* Initial scenario data */
    u.branch_x = 0.1;
    u.fault_b = 0;
    u.vr_inf_set = 1.0;
    u.vi_inf_set = 0;

    /* Set initial conditions */
    double delta, omega, Eqp, psidp, psiqp, Edp, psiqpp, psidpp, psipp, vd, vq,
        id, iq, vr, vi, ir, ii, p, q, Te, g, b, ksat, *yy_;

    yy_ = NV_DATA_S(yy);

    /* Start by assuming state variables */
    yy_[0] = delta = 0.55399038;
    yy_[1] = omega = 0;
    yy_[2] = Eqp = 0.995472581; 
    yy_[3] = psidp = 0.971299567;
    yy_[4] = psiqp = 0.306880069; 
    yy_[5] = Edp = 0;

    /* Assumed terminal conditions for system */
    vr = 0.9949877346411762;
    vi = 0.09999703952427966;
    p = 1.0;
    q = 0.05013;

    /* Calculate the rest of the initial conditions */
    yy_[6] = psiqpp = -psiqp*u.Xq4 - Edp*u.Xq5;
    yy_[7] = psidpp = psidp*u.Xd4 + Eqp*u.Xd5;
    yy_[8] = psipp = sqrt(psiqpp*psiqpp + psidpp*psidpp);
    yy_[9] = ksat = u.SB*pow(psipp - u.SA, 2);
    yy_[10] = vd = -psiqpp*(1 + omega); 
    yy_[11] = vq = psidpp*(1 + omega); 
    ir = (p*vr + q*vi) / (vr*vr + vi*vi);
    ii = (p*vi - q*vr) / (vr*vr + vi*vi);
    id = ir*sin(delta) -ii*cos(delta);
    iq = ir*cos(delta) +ii*sin(delta);
    yy_[12] = Te = (psidpp - id*u.Xdpp)*iq - (psiqpp - iq*u.Xdpp)*id; 
    yy_[13] = id;
    yy_[14] = iq;
    yy_[15] = ir;
    yy_[16] = ii;
    yy_[17] = u.pmech_set = Te; 
    yy_[18] = u.efd_set = Eqp + u.Xd1*(id 
        + u.Xd3*(Eqp - psidp - u.Xd2*id)) + psidpp*ksat;
    yy_[19] = u.g*(vd*sin(delta)+vq*cos(delta)) 
        - u.b*(vd*-cos(delta) + vq*sin(delta)); /* machine inort, real */
    yy_[20] = u.b*(vd*sin(delta)+vq*cos(delta)) 
        + u.g*(vd*-cos(delta) + vq*sin(delta)); /* machine inort, imag */
    yy_[21] = vr;   /* v1r - bus1 has GENROU */
    yy_[22] = vi;   /* v1i - bus1 has GENROU  */
    yy_[23] = 1.0;   /* v2r - bus2 has infinite bus */
    yy_[24] = 0;   /* v2i - bus2 has infinite bus  */

    CHECK_ALLOCATION( yy0, N_VClone(yy) );
    CHECK_ALLOCATION( yp0, N_VClone(yp) );
    CHECK_OUTPUT( IDAInit(ida_mem, ex1_residual, 0, yy, yp) );
    CHECK_OUTPUT( IDASetUserData(ida_mem, &u) );
    CHECK_OUTPUT( IDASStolerances(ida_mem, rtol, atol) );
    CHECK_OUTPUT( IDASetMaxNumSteps(ida_mem, maxsteps) );
    CHECK_ALLOCATION( yid, N_VClone(yy) );
    
    /* First six variables are differential */
    for (int i = 0; i < 6; ++i) NV_DATA_S(yid)[i] = 1;
    for (int i = 6; i < ny; ++i) NV_DATA_S(yid)[i] = 0;

    /* Initialize IDA solver */
    CHECK_OUTPUT( IDASetId(ida_mem, yid) );
    CHECK_OUTPUT( IDASetSuppressAlg(ida_mem, SUNTRUE) );
    CHECK_ALLOCATION( J, SUNSparseMatrix(ny, ny, nnz, CSC_MAT, context) );
    CHECK_ALLOCATION(linsol, SUNLinSol_KLU(yy, J, context) );
    CHECK_OUTPUT( IDASetLinearSolver(ida_mem, linsol, J) );
    CHECK_OUTPUT( IDASetJacFn(ida_mem, ex1_jacobian) );
    CHECK_OUTPUT( IDAReInit(ida_mem, 0, yy, yp) );

    /* Check for any residual issues starting out */
    ex1_residual(0, yy, yp, rr, &u);
    double *r = NV_DATA_S(rr);
    for (int i = 0 ; i < ny; ++i)
    {
        if (fabs(r[i]) > 1e-4) printf("R[%d] = %.16g\n", i, r[i]);
    }

    /* Output file with headers */
    FILE *f = fopen("example1_results.csv", "w");
    if (!f) printf("ERROR writing to output file!\n");
    fprintf(f, "%s,%s,%s", "i", "IDA Return Value", "t");
    for (int j = 0; j < ny; ++j) fprintf(f, ",y[%d]", j);
    fprintf(f, "\n");
    double dt = 1/4.0/60.0;
    fprintf(f, "%d,%d,%.16g", 0, 0, 0.0);
    for (int j = 0; j < ny; ++j) fprintf(f, ",%.16g", yy_[j]);
    fprintf(f, "\n");

    /* Main solution loop */

    double start = (double) clock();
    for (int i = 0; i < 7200; ++i)
    {
        if (i == 240) /* Fault occurs at t=1.0 seconds */
        {
            u.fault_b = -1000;
            CHECK_OUTPUT( IDAReInit(ida_mem, (i)*dt, yy, yp) );
        }
        if (i == 264) /* Fault is cleared at t=1.1 seconds */
        {
            u.fault_b = 0;
            CHECK_OUTPUT( IDAReInit(ida_mem, (i)*dt, yy, yp) );
        }
        sunrealtype tout = (i + 1)*dt, tret = tout;
        retval = IDASolve(ida_mem, tout, &tret, yy, yp, IDA_NORMAL);
        fprintf(f, "%d,%d,%.16g", i+1, retval, tret);
        for (int j = 0; j < ny; ++j) fprintf(f, ",%.16g", yy_[j]);
        fprintf(f, "\n");
    }
    printf("Complete in %.4g seconds\n", (clock() - start) / CLOCKS_PER_SEC);
    fclose(f);

    /* Clean up */
    N_VDestroy(yy);
    N_VDestroy(yp);
    N_VDestroy(yy0);
    N_VDestroy(yp0);
    SUNLinSolFree_KLU(linsol);
    SUNMatDestroy_Sparse(J);
    IDAFree(&ida_mem);
    CHECK_OUTPUT( SUNContext_Free(&context) );

    return 0;
}