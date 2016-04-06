#include <omp.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <stdio.h>

#include "p3dBlob.h"
#include "p3dTime.h"
#include "p3dAuth.h"


#define MAXROT 512  /*  The maximum # of MIL rotations	*/
#define MAXPOINTS 1024

#define FOO_I(i,j,N) ( (j)*(N) + (i) ) 

#define BAD   1
#define GOOD  0
#define ABEND 1


#define EIGVALS(i,j) *(eigvals+(i)*dim+(j))
#define EIGVECS(i,j) *(eigvecs+(i)*dim+(j))
#define A(X,Y)  *(a + (X)*dim + (Y))
#define A_IN(X,Y)  *(in + (X)*dim + (Y))
#define BIGNUM  1.0E15

/*  The next four subroutines are used with permission of Triaxis, Inc,
Los Alamos, NM 87544.  They are from  Version 2.0 of their Math++ library */

/**********************************************************************
 *		Copyright (c) Triakis, Inc. 1989		      *
 *								      *
 *  Function: linSystemSolve                                          *
 *                                                                    *
 *  Version:  2.0                                                     *
 *                                                                    *
 *  Purpose:  Solution of linear system: Ax = b.                      *
 *            Used in conjunction with function linLuDecomp.          *
 *            Not used if linLuDecomp has detected a singularity.     *
 *                                                                    *
 *            linSystemSolve does the forward elimination             *
 *            and back substitution.  Splitting the process into      *
 *            two parts allows solutions for multiple b vectors       *
 *            with a minimum of recomputation.                        *
 *                                                                    *
 *  Returns:  Status ( 0 = success and 1 = failure ).                 *
 *                                                                    *
 *  Call:     linSystemSolve( dim, n, a, b, pivot );                  *
 *                                                                    *
 *      input:                                                        *
 *                                                                    *
 *          dim = Maximum dimension of A in main.                     *
 *            n = order of matrix.                                    *
 *            a = Triangularized matrix obtained from linLuDecomp.    *
 *            b = Right hand side vector.                             *
 *        pivot = Pivot vector obtained from linLuDecomp.             *
 *                                                                    *
 *      output:                                                       *
 *                                                                    *
 *            b = solution vector, x.                                 *
 *                                                                    *
 *  External functions:                                               *
 *            none.                                                   *
 *                                                                    *
 * Algorithm: Forsythe, Malcolm, and Moler from the text              *
 *            "Computer Methods for Mathematical Computations",       *
 *            1977, Prentice Hall.                                    *
 *                                                                    *
 *            Coded in structured 'C' by Tim Julian and Eric Blommer. *
 *                                                                    *
 **********************************************************************/


int linSystemSolve(int dim, int order, double* a, double* b, int* pivot) {

    int i,
            k,
            m;
    double temp;

    /* Forward elimination, trivial case. */
    if (order == 1) {
        if (A(0, 0) == 0.0) {
            return BAD;
        } else {
            b[0] /= A(0, 0);
            return GOOD;
        }
    }
    /* Forward elimination, non-trivial case. */
    for (k = 0; k < order - 1; ++k) {
        m = pivot[k];
        temp = b[m];
        b[m] = b[k];
        b[k] = temp;
        for (i = k + 1; i < order; ++i) {
            b[i] += A(i, k) * temp;
        }
    }
    /* Back substitution */
    for (k = order - 1; k >= 0; --k) {
        if (A(k, k) == 0) {
            return BAD;
        }
        b[k] /= A(k, k);
        temp = -b[k];
        for (i = 0; i < k; ++i) {
            b[i] += A(i, k) * temp;
        }
    }

    return GOOD;

} /* End linSystemSolve */

/**********************************************************************
 *		Copyright (c) Triakis, Inc. 1989		      *
 *								      *
 *  Function: linLuDecomp                                             *
 *                                                                    *
 *  Version:  2.0                                                     *
 *                                                                    *
 *  Purpose:  Decomposes a real general matrix by LU decomposition    *
 *            and estimates the condition of the matrix.              *
 *                                                                    *
 *            Function linSystemSolve computes solution to linear     *
 *	      systems after decomposition by this function.           *
 *                                                                    *
 *  Returns:  Condnum = An estimate of the condition number of "A".   *
 *            A good rule of thumb in using condnum is that by taking *
 *            the log10 of the condnum and subtracting that value     *
 *            from the number of double precision significant digits  *
 *            availible, the remainder is the approximate number of   *
 *            significant digits in your answer.                      *
 *            Condnum is set to 3.33e+33 if exact singularity         *
 *            is detected.                                            *
 *                                                                    *
 *  Call:     status = linLuDecomp( dim, order, a, pivot, &condnum ); *
 *                                                                    *
 *       input:                                                       *
 *          dim  = Maximum dimensioned columns of "A".                *
 *         order = Order of the matrix  (actual size of "A").         *
 *            a  = Square matrix to be decomposed.                    *
 *                                                                    *
 *      output:                                                       *
 *            a = Contains an upper triangular matrix u and a         *
 *                lower triangular matrix of multipliers.             *
 *      condnum = An estimate of the condition number of A.           *
 *                For the linear system A*x = b, changes in A and b   *
 *                may cause changes condnum times as large in x.      *
 *                Condnum is set to 3.33e+33 if exact singularity     *
 *                is detected.                                        *
 *        pivot = The pivot vector.                                   *
 *                pivot[i] is the index of the ith pivot row.         *
 *                pivot[order-1] is (-1) ** (number of interchanges)  *
 *                                                                    *
 *  External functions:                                               *
 *            linSystemSolve(),                                       *
 *            fabs(),                                                 *
 *            calloc(),                                               *
 *            free().                                                 *
 *                                                                    *
 * Algorithm: Forsythe, Malcolm, and Moler from the text              *
 *            "Computer Methods for Mathematical Computations",       *
 *            1977, Prentice Hall( pages 30 - 62 ).                   *
 *            Also see the text "Linpack Users' Guide" by Dongarra,   *
 *            Moler, Bunch, and Stewart, 1979 by Siam.                *
 *            ( pages 1.1 - 1.34 and source listings in back )        *
 *            Note the routines dgeco, dgefa, dgedi, dgesl,           *
 *            ddot, idamax, dscal, daxpy, and dasum.                  *
 *            Corresponding single precision routines are sgeco,      *
 *            sgefa, sgedi, sgesl, sdot, isamax, sscal, saxpy,        *
 *            and sasum.                                              *
 *                                                                    *
 *  Misc:    The determinant of A can be obtained on output by        *
 *           det_a = pivot[order-1] * A(0,0) * A(1,1) *               *
 *                                ... * A(order-1,order-1).           *
 *           See the book Computer Methods for Mathematical           *
 *           Computations by Forsythe, Malcolm and Moler              *
 *           pages 30 - 62 for further explanation.                   *
 *                                                                    *
 *           Coded in structured 'C' by Tim Julian and Eric Blommer.  *
 *                                                                    *
 **********************************************************************/


int linLuDecomp(int dim, int order, double* a, int* pivot, double* condnum) {
    /* Local variable declarations. */
    double *work,
            ek,
            temp,
            anorm,
            ynorm,
            znorm;
    int i,
            j,
            k,
            kplus1,
            m,
            flag;
    /* Function declarations. */
    //char           *calloc ();
    //int		   linSystemSolve ();

    /* Check that "a" matrix order not larger than maximum size */
    if (order > dim) {
        (void) fprintf(stderr, " In linLuDecomp: order > Max. dim.\n");
        (void) exit(ABEND);
    }
    /* Check that "a" matrix order not smaller than one. */
    if (order < 1) {
        (void) fprintf(stderr, " In linLuDecomp: order < 1\n");
        (void) exit(ABEND);
    }

    /* Dynamically allocate work vector appropriate for problem. */
    if ((work = (double *) calloc(order, sizeof (double))) == NULL) {
        (void) fprintf(stderr, " In linLuDecomp: can't allocate work vector.\n");
        (void) exit(ABEND);
    }

    pivot[order - 1] = 1;
    /* Treat the trivial case. */
    if (order == 1) {
        if (A(0, 0) != 0.0) {
            *condnum = 1.0;
            free(work);
            return GOOD;
        } else {
            /* Exact singularity detected for trivial case */
            *condnum = 3.33e33;
            free(work);
            return BAD;
        }
    }
    /* Compute the 1-norm of matrix 'a'.                    */
    /* The 1-norm of matrix 'a' is defined as the largest   */
    /* sum of the absolute values of each of the elements   */
    /* of a single column vector across all column vectors. */
    anorm = 0.0;
    for (j = 0; j < order; ++j) {
        temp = 0.0;
        for (i = 0; i < order; ++i) {
            temp += fabs(A(i, j));
        }
        if (temp > anorm) {
            anorm = temp;
        }
    }
    /* Gaussian elimination with partial pivoting. */
    /* See Linpack function dgefa or sgefa         */
    for (k = 0; k < order - 1; ++k) {
        kplus1 = k + 1;
        /* Find pivot. */
        m = k;
        for (i = kplus1; i < order; ++i) {
            if (fabs(A(i, k)) > fabs(A(m, k))) {
                m = i;
            }
        }
        pivot[k] = m;
        /* See Linpack function dgedi or sgedi */
        if (m != k) {
            pivot[order - 1] = -pivot[order - 1];
        }
        /* Skip step if pivot is zero. */
        if (A(m, k) != 0.0) {

            if (m != k) {
                temp = A(m, k);
                A(m, k) = A(k, k);
                A(k, k) = temp;
            }
            /* Compute multipliers. */
            temp = -1.0 / A(k, k);

            for (i = kplus1; i < order; ++i) {
                A(i, k) *= temp;
            }
            /* Row elimination with column indexing */
            for (j = kplus1; j < order; ++j) {
                temp = A(m, j);
                A(m, j) = A(k, j);
                A(k, j) = temp;
                if (temp != 0.0) {
                    for (i = kplus1; i < order; ++i) {
                        A(i, j) += A(i, k) * temp;
                    }
                }
            }
        }
    } /* End k loop. */
    /* Calculate condition number of matrix.                         */
    /* condnum = (1-norm of a)*(an estimate of 1-norm of a-inverse)  */
    /* Estimate obtained by one step of inverse iteration for the    */
    /* small singular vector.  This involves solving two systems     */
    /* of equations, (a-transpose)*y = e and a*z = y where e         */
    /* is a vector of +1 or -1 chosen to cause growth in y.          */
    /* Estimate = (1-norm of z)/(1-norm of y)                        */
    /* First, solve (a-transpose)*y = e                              */
    for (k = 0; k < order; ++k) {
        temp = 0.0;
        if (k != 0) {
            for (i = 0; i < k - 1; ++i) {
                temp += A(i, k) * work[i];
            }
        }
        if (temp < 0.0) {
            ek = -1.0;
        } else {
            ek = 1.0;
        }
        if (A(k, k) == 0.0) {
            /* Exact singularity detected for non-trivial case */
            *condnum = 3.33e33;
            free(work);
            return BAD;
        }
        work[k] = -(ek + temp) / A(k, k);
    }
    for (k = order - 2; k >= 0; --k) {
        temp = 0.0;
        for (i = k + 1; i < order; ++i) {
            temp += A(i, k) * work[k];
        }
        work[k] = temp;
        m = pivot[k];
        if (m != k) {
            temp = work[m];
            work[m] = work[k];
            work[k] = temp;
        }
    }

    ynorm = 0.0;
    for (i = 0; i < order; ++i) {
        ynorm += fabs(work[i]);
    }
    /* Solve a*z = y. */
    flag = linSystemSolve(dim, order, a, work, pivot);
    if (flag == 1) {
        (void) fprintf(stderr, "In linLuDecomp: error detected in function linSystemSolve.\n");
        (void) exit(ABEND);
    }

    znorm = 0.0;
    for (i = 0; i < order; ++i) {
        znorm += fabs(work[i]);
    }
    /* Estimate condition number. */
    *condnum = anorm * znorm / ynorm;
    if (*condnum < 1.0) {
        *condnum = 1.0;
    }

    free(work);

    return GOOD;

} /* End linLuDecomp */

/**********************************************************************
 *		Copyright (c) Triakis, Inc. 1989		      *
 *								      *
 *  Function: linInverseMatrix                                        *
 *                                                                    *
 *  Version:  2.0                                                     *
 *                                                                    *
 *  Purpose:  Determines the inverse of a matrix via LU decomposition.*
 *                                                                    *
 *  Returns:  Nothing of value.                                       *
 *                                                                    *
 *  Call:     linInverseMatrix( dim, order, a, in );                  *
 *                                                                    *
 *       input:                                                       *
 *          dim = maximum dimensioned columns of "A".                 *
 *         order = order of the matrix  (actual size of "A").         *
 *            a = square matrix to be decomposed.                     *
 *                                                                    *
 *         output:                                                    *
 *           in = Inverse of matrix "a".                              *
 *                                                                    *
 *  External functions:                                               *
 *            linSystemSolve(),					      *
 *            linLuDecomp(),                                          *
 *            calloc(),                                               *
 *            free().                                                 *
 *                                                                    *
 * Algorithm: Tim Julian and Eric Blommer.                            *
 *                                                                    *
 **********************************************************************/


void linInverseMatrix(int dim, int order, double* a, double* in) {
    int i,
            j; /* work variables */
    int status; /* return flag to check functions */
    int *pivot; /* pivoting from function linLuDecomp */
    //int             linLuDecomp ();		/* decompose matrix a into LU matrix */
    //int             linSystemSolve ();		/* solve system of simultaneous eqns. */
    double *b; /* B vector */
    double condnum; /* matrix condition number */
    //char           *calloc ();

    if (order > dim) {
        (void) fprintf(stderr, "Inverse: order is greater than dim!\n");
        (void) exit(ABEND);
    }

    /* allocate work vectors */

    if ((b = (double *) calloc(order, sizeof (double))) == NULL) {
        (void) fprintf(stderr, " In linInverseMatrix, can't allocate b vector.\n");
        (void) exit(ABEND);
    }

    if ((pivot = (int *) calloc(order, sizeof (int))) == NULL) {
        (void) fprintf(stderr, " In linInverseMatrix, can't allocate pivot vector.\n");
        (void) exit(ABEND);
    }

    /* Initialize data values */
    condnum = 0.0;

    /* decompose matrix */
    status = linLuDecomp(dim, order, a, pivot, &condnum);

    if (status == 1) {
        (void) fprintf(stderr, "In Inverse, error detected in function linLuDecomp.\n");
        (void) exit(ABEND);
    }

    if (condnum > BIGNUM) {
        (void) fprintf(stderr, "Inverse: Matrix singular to double precision.\n");
        (void) fprintf(stderr, "Condition number : %e\n", condnum);
        (void) exit(ABEND);
    }
    /* backsolve matrix */
    for (i = 0; i < order; ++i) {
        for (j = 0; j < order; ++j) {
            if (j == i) {
                b[j] = 1.0;
            } else {
                b[j] = 0.0;
            }
        }
        status = linSystemSolve(dim, order, a, b, pivot);

        if (status == 1) {
            (void) fprintf(stderr, "In Inverse, error detected in function linSystemSolve.\n");
            (void) exit(ABEND);
        }
        for (j = 0; j < order; ++j) {
            A_IN(j, i) = b[j];
        }
    }

    free(b);
    free(pivot);

} /* end linInverseMatrix */

/**********************************************************************
 *		Copyright (c) Triakis, Inc. 1989		      *
 *                                                                    *
 *  Function: eigJacobi                                               *
 *                                                                    *
 *  Version:  2.0                                                     *
 *                                                                    *
 *  Purpose:  Computes the eigenvalues and eigenvectors of a real     *
 *            symmetric matrix via the Jacobi method.                 *
 *                                                                    *
 *  Returns:  status = (0 = success, 1 = failure)                     *
 *                                                                    *
 *  Call:     status =  eigJacobi(  d, e1, e2, o, e, it, mit ), where:*
 *                                                                    *
 *       input:                                                       *
 *              d = maximum dimension of matrix eigenvalue matrix     *
 *             e1 = eigenvalue matrix.                                *
 *             e2 = eigenvector matrix.                               *
 *              o = order of the matrices                             *
 *              e = convergence criterion                             *
 *            mit = maximum number of iterations allowed              *
 *                                                                    *
 *      output:                                                       *
 *             e1 = eigenvalue matrix.                                *
 *             e2 = eigenvector matrix.                               *
 *             it = actual number of iterations                       *
 *                                                                    *
 *  External functions:                                               *
 *            fabs(),                                                 *
 *            calloc(),                                               *
 *            free(),                                                 *
 *            sqrt().                                                 *
 *                                                                    *
 *  Author(s): Tim Julian and Eric Blommer                            *
 *                                                                    *
 *  Misc:    None.                                                    *
 *                                                                    *
 *           Coded in structured 'C' by Tim Julian and Eric Blommer.  *
 *                                                                    *
 **********************************************************************/



int eigJacobi(int dim,
        double* eigvals,
        double* eigvecs,
        int order,
        double eps,
        int* numiters,
        int maxiters) {

    double vo,
            f,
            u,
            uf,
            alpha,
            beta,
            c,
            s;
    double *temp1,
            *temp2,
            *temp3;
    int i,
            j,
            l,
            m,
            p,
            q;
    //char           *calloc ();

    /* dynamically allocate temporary vectors */
    if ((temp1 = (double *) calloc(dim, sizeof (double))) == NULL) {
        (void) fprintf(stderr, "Can't allocate vector temp1\n");
        (void) exit(ABEND);
    }

    if ((temp2 = (double *) calloc(dim, sizeof (double))) == NULL) {
        (void) fprintf(stderr, "Can't allocate vector temp2\n");
        (void) exit(ABEND);
    }

    if ((temp3 = (double *) calloc(dim, sizeof (double))) == NULL) {
        (void) fprintf(stderr, "Can't allocate vector temp3\n");
        (void) exit(ABEND);
    }

    /* make identity matrix */
    for (i = 0; i < order; ++i)
        for (j = 0; j <= i; ++j) {
            if (i == j) {
                EIGVECS(i, j) = 1.0;
            } else {
                EIGVECS(i, j) = EIGVECS(j, i) = 0.0;
            }
        }


    vo = 0.0;

    for (i = 0; i < order; ++i)
        for (j = 0; j < order; ++j) {
            if (i != j) {
                vo += EIGVALS(i, j) * EIGVALS(i, j);
            }
        }


    u = sqrt(vo) / (double) order;
    uf = eps * u;

    for (uf = eps * u; uf < u; u /= order) {
        for (l = 0; l < order - 1; ++l) {
            for (m = l + 1; m < order; ++m) {
                if (fabs(EIGVALS(l, m)) >= u) {
                    p = l;
                    q = m;
                    for (i = 0; i < order; ++i) {
                        temp1[i] = EIGVALS(p, i);
                        temp2[i] = EIGVALS(i, p);
                        temp3[i] = EIGVECS(i, p);
                    }

                    f = EIGVALS(p, p);
                    alpha = 0.5 * (EIGVALS(p, p) - EIGVALS(q, q));
                    beta = sqrt(EIGVALS(p, q) * EIGVALS(p, q) + alpha * alpha);
                    c = sqrt(0.5 + fabs(alpha) / (2.0 * beta));
                    s = alpha * (-EIGVALS(p, q)) / (2.0 * beta * fabs(alpha) * c);

                    for (j = 0; j < order; ++j) {
                        if (j != p) {
                            if (j != q) {
                                EIGVALS(p, j) = c * EIGVALS(p, j) - s * EIGVALS(q, j);
                                EIGVALS(q, j) = s * temp1[j] + c * EIGVALS(q, j);
                                EIGVALS(j, p) = c * EIGVALS(j, p) - s * EIGVALS(j, q);
                                EIGVALS(j, q) = s * temp2[j] + c * EIGVALS(j, q);
                            }
                        }
                        EIGVECS(j, p) = c * EIGVECS(j, p) - s * EIGVECS(j, q);
                        EIGVECS(j, q) = s * temp3[j] + c * EIGVECS(j, q);
                    }

                    EIGVALS(p, p) = c * c * EIGVALS(p, p) + s * s *
                            EIGVALS(q, q) - 2.0 * c * s * EIGVALS(p, q);
                    EIGVALS(q, q) = s * s * f + c * c * EIGVALS(q, q) +
                            2.0 * c * s * EIGVALS(p, q);
                    EIGVALS(p, q) = 0.0;
                    EIGVALS(q, p) = 0.0;

                    if (++*numiters > maxiters) {
                        (void) fprintf(stderr, "number of iterations exceeded %d", maxiters);
                        return BAD;
                    }
                }
            }
        }
    }

    free(temp1);
    free(temp2);
    free(temp3);

    return GOOD;

} /* end eigJacobi */

/*int customPrint(const char* str, ...) {
    printf(str);
    return printf("\n");
}*/

/************************************************************************/
/*                                                                      */
/* 	Software to perform three-dimensional stereology.  The methods 	    */
/*	are based on a three-dimensional version of the directed secant	    */
/*	algorithm (Saltykov, 1958). 3-D mean intercept length vectors,      */
/*  MIL(theta,phi), are calculated and fit to an ellipsoid from         */
/*  which the three principal MIL vectors are determined.			    */
/*									                                    */
/************************************************************************/


/************************************************************************/
/*  Function:  _getMILs                                                 */
/*									                                    */
/*	Function to scan the image v_array using a three-dimensional 	    */
/*	version of the directed secant method.  The image is scanned by     */
/*	a 3-D test grid at randomly determined orientations.  The 	        */
/*	orientations are defined by spherical angles (theta, phi).  	    */
/*	Theta is the rotation about the x-axis and phi is the rotation 	    */
/*	about the z-axis.  Based on the threshold value, the image is 	    */
/*	converted to a binary format and the bone volume fraction is 	    */
/*	determined.  The image is then systematically scanned and 	        */
/*	intersections are recorded when the binary value of the current	    */
/*	voxel differs from the binary value of the previous voxel.	        */
/*	Standard morphology parameters are calculated based on the 	        */
/*	parallel plate model:						                        */
/*		BV/TV = # bone voxels/total # voxels in test sphere	            */
/*		Tb.N = total # intersections/total length of test lines	        */
/*		Tb.Th = BV/TV / Tb.N					                        */
/*		Tb.Sp = (1-BV/TV) / Tb.N				                        */
/*		BS/BV = 2 * Tb.N / BV/TV				                        */
/*	Three-dimensional mean intercept length vectors are calculated 	    */
/*	for each rotation as:						                        */
/*		MIL(theta, phi) = (BV/TV) / ((1/2) * #intersections	            */
/*					per unit test line length)	                        */
/*									                                    */

/************************************************************************/

void _getMILs(
        unsigned char* in_im, // A_IN: Input segmented (binary) volume
        unsigned char* msk_im,
        double* rot_theta, // OUT: rotations in THETA angle
        double* rot_phi, // OUT: rotations in PHI angle
        double* mil, // OUT: Mean Intercept Lengths
        const unsigned int dimx,
        const unsigned int dimy,
        const unsigned int dimz,
        const double voxelsize
        ) {
    int i;
    double s;

    // Indexes for volume scanning:
    double x, y, z;
    double prec_x, prec_y, prec_z;

    double end_x1, end_y1, end_z1;
    double end_x2, end_y2, end_z2;
    int start_x, start_y, start_z;


    // Test line direction vectors:
    double r_0, r_1, r_2;
    double normr_0, normr_1, normr_2;
    int ct;

    // Variables for length management:
    double length, tot_length;
    double bvf; // Bone Volume Fraction (i.e. BV/TV)    
    double intersect_ct;

    /*double	t_totlength, t_intersect_ct;
    double	pl, sv, tb, tpd, tps;*/

    unsigned char prec_status, curr_status;
    unsigned int iseed;


    // Variables for rotation cycle:
    int ct_rot;

    // Init random seeds:
    iseed = (unsigned int) time(NULL);
    srand(iseed);


    // Get BV/TV:
    s = 0.0;
#pragma omp parallel for reduction (+ : s)
    for (i = 0; i < (dimx * dimy * dimz); i++)
        if (in_im[i] == OBJECT) s++;

    bvf = s / ((double) (dimx * dimy * dimz));

    /*t_totlength = 0.0;
    t_intersect_ct = 0.0;*/

    // For each rotations:
#pragma omp parallel for private( ct, tot_length, intersect_ct, start_x, start_y, start_z, x, y, z, prec_x, prec_y, prec_z, end_x1, end_y1, end_z1, end_x2, end_y2, end_z2, prec_status, curr_status, length, r_0, r_1, r_2, normr_0, normr_1, normr_2 ) 
    for (ct_rot = 0; ct_rot < MAXROT; ct_rot++) {
        rot_theta[ct_rot] = (rand() / ((double) RAND_MAX + 1)) * M_PI;
        rot_phi[ct_rot] = (rand() / ((double) RAND_MAX + 1)) * M_PI;

        /*fprintf(stderr, "rot_theta[ct_rot]: %0.5f\n", rot_theta[ct_rot]);
        fprintf(stderr, "rot_phi[ct_rot]: %0.5f\n", rot_phi[ct_rot]);*/

        //   Define direction vector (versor):
        r_0 = cos(rot_phi[ct_rot]) * sin(rot_theta[ct_rot]);
        r_1 = sin(rot_phi[ct_rot]) * sin(rot_theta[ct_rot]);
        r_2 = cos(rot_theta[ct_rot]);

        // Normalize direction vector to +1 for the maximum component,
        // doing so we can save iteration on next cycle:
        if ((fabs(r_0) >= fabs(r_1)) && (fabs(r_0) >= fabs(r_2))) {
            normr_0 = 1.0;
            normr_1 = r_1 / r_0;
            normr_2 = r_2 / r_0;
        } else if ((fabs(r_1) >= fabs(r_0)) && (fabs(r_1) >= fabs(r_2))) {
            normr_1 = 1.0;
            normr_0 = r_0 / r_1;
            normr_2 = r_2 / r_1;
        } else {
            normr_2 = 1.0;
            normr_0 = r_0 / r_2;
            normr_1 = r_1 / r_2;
        }

        // Make sure signs are correct:		
        normr_0 = (r_0 > 0) ? fabs(normr_0) : (-fabs(normr_0));
        normr_1 = (r_1 > 0) ? fabs(normr_1) : (-fabs(normr_1));
        normr_2 = (r_2 > 0) ? fabs(normr_2) : (-fabs(normr_2));


        tot_length = 0.0;
        intersect_ct = 0.0;

        for (ct = 0; ct < MAXPOINTS; ct++) {
            start_x = (int) (rand() / (((double) RAND_MAX + 1) / dimx));
            start_y = (int) (rand() / (((double) RAND_MAX + 1) / dimy));
            start_z = (int) (rand() / (((double) RAND_MAX + 1) / dimz));
            
            // Check if inside MSK (if mask):
            if (msk_im != NULL)
            {
                if (msk_im[ I(start_x, start_y, start_z, dimx, dimy)] == BACKGROUND)
                {
                    continue;
                }
                    
            }

            x = start_x;
            y = start_y;
            z = start_z;

            // Reset status:
            prec_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];
            curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

            // Init variables:
            prec_x = x;
            prec_y = y;
            prec_z = z;

            // Explore the incremental and decremental sides while edges of 
            // VOI are reached:
            if (msk_im == NULL) {
                // Explore the incremental versus while edges of VOI are reached:
                while ((x >= 0) && (x < dimx) &&
                        (y >= 0) && (y < dimy) &&
                        (z >= 0) && (z < dimz)) {
                    curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

                    // If we reach the object save coords:
                    if (curr_status != prec_status) {
                        intersect_ct++;
                        prec_status = curr_status;
                    }

                    prec_x = x;
                    prec_y = y;
                    prec_z = z;

                    x = x + normr_0;
                    y = y + normr_1;
                    z = z + normr_2;
                }

                // Get end point of the test line:
                end_x1 = prec_x;
                end_y1 = prec_y;
                end_z1 = prec_z;

                // Reset "center" of the line:
                x = start_x;
                y = start_y;
                z = start_z;

                // Reset status:
                prec_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];
                curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

                // Init variables:
                prec_x = x;
                prec_y = y;
                prec_z = z;

                // Explore the decremental versus while edges of VOI are reached:
                while ((x >= 0) && (x < dimx) &&
                        (y >= 0) && (y < dimy) &&
                        (z >= 0) && (z < dimz)) {
                    curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

                    // If we reach the object save coords:
                    if (curr_status != prec_status) {
                        intersect_ct++;
                        prec_status = curr_status;
                    }

                    prec_x = x;
                    prec_y = y;
                    prec_z = z;

                    x = x - normr_0;
                    y = y - normr_1;
                    z = z - normr_2;
                }

                // Get the other end point of the test line:
                end_x2 = prec_x;
                end_y2 = prec_y;
                end_z2 = prec_z;
            } else {
                // Explore the incremental and decremental sides while edges of 
                // mask are reached of image edges are reached:

                // Explore the incremental versus while edges of VOI are reached:
                while ((x >= 0) && (x < dimx) &&
                        (y >= 0) && (y < dimy) &&
                        (z >= 0) && (z < dimz) &&
                        (msk_im[ I((int) x, (int) y, (int) z, dimx, dimy)] == OBJECT)
                        ) {
                    curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

                    // If we reach the object save coords:
                    if (curr_status != prec_status) {
                        intersect_ct++;
                        prec_status = curr_status;
                    }

                    prec_x = x;
                    prec_y = y;
                    prec_z = z;

                    x = x + normr_0;
                    y = y + normr_1;
                    z = z + normr_2;
                }

                // Get end point of the test line:
                end_x1 = prec_x;
                end_y1 = prec_y;
                end_z1 = prec_z;

                // Reset "center" of the line:
                x = start_x;
                y = start_y;
                z = start_z;

                // Reset status:
                prec_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];
                curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

                // Init variables:
                prec_x = x;
                prec_y = y;
                prec_z = z;

                // Explore the decremental versus while edges of VOI are reached:
                while ( (x >= 0) && (x < dimx) &&
                        (y >= 0) && (y < dimy) &&
                        (z >= 0) && (z < dimz) &&
                        (msk_im[ I((int) x, (int) y, (int) z, dimx, dimy)] == OBJECT)
                        ) {
                    curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

                    // If we reach the object save coords:
                    if (curr_status != prec_status) {
                        intersect_ct++;
                        prec_status = curr_status;
                    }

                    prec_x = x;
                    prec_y = y;
                    prec_z = z;

                    x = x - normr_0;
                    y = y - normr_1;
                    z = z - normr_2;
                }

                // Get the other end point of the test line:
                end_x2 = prec_x;
                end_y2 = prec_y;
                end_z2 = prec_z;

            }

            // Get distance:
            length = (end_x1 - end_x2)*(end_x1 - end_x2) +
                    (end_y1 - end_y2)*(end_y1 - end_y2) +
                    (end_z1 - end_z2)*(end_z1 - end_z2);
            length = sqrt(length);

            // Add current length to array for this grid:
            tot_length += length*voxelsize;
        }

        // Save current MIL:
        mil[ ct_rot ] = (2.0 * bvf) * (tot_length / intersect_ct);

        /*t_totlength += tot_length;
        t_intersect_ct += intersect_ct;*/

    } // end for each rotation

    /*pl = t_intersect_ct / t_totlength;

    sv = 2 * pl / bvf;		// BS/BV, [mm^2 mm^-3] 
    tb = bvf / pl;			// Tb.Th, [mm] 
    tpd = pl;				// Tb.N, [mm^-1] 
    tps = (1.0-bvf)/pl;		// Tb.Sp, [mm] */

}

/************************************************************************/
/* Function:  _MILs2fit()						*/
/*									*/
/*	This function fits the 3-D MIL data to an ellipsoid, and 	*/
/*	determines the goodness of the fit.  These methods are based on */
/*	2-D methods of Whitehouse (JMicroscopy 101:153-168, 1974) and	*/
/*	Harrigan and Mann (JMatSci 19:761-767, 1984).  The approach is 	*/
/*	to plot the locus of the end points of the MIL vectors issuing 	*/
/*	from a common center and fit them to an ellipsoid of general 	*/
/*	formula:							*/
/*									*/
/*	A*n1^2 + B*n2^2 + C*n3^2 + D*n1*n2 + E*n1*n3 + F*n2*n3 = 1/L^2	*/
/*									*/
/*	where L is the length of the MIL, ni are the direction cosines	*/
/*	between L and the base vectors in an arbitrary coordinate 	*/
/*	system and A...F are the ellipsoid coefficients.		*/
/*									*/
/*	For this code, a multivariable linear least squares fitting 	*/
/*	technique is used to fit the data by solving the linear system 	*/
/*	A * x = b, where A contains the projection data for each 	*/
/*	rotation, x is a column vector of the ellipsoid coefficients 	*/
/*	and b is a column vector of 1/L^2, where L is the magnitude of 	*/
/*	the MIL vector.  The solution is formulated as:			*/
/*									*/
/*			x = (At * A)^-1 * At * b			*/
/*									*/

/************************************************************************/
void _MILs2fit(
        double* rot_theta, // A_IN: rotations in THETA angle
        double* rot_phi, // A_IN: rotations in PHI angle
        double* mil, // A_IN: Mean Intercept Lengths
        double* xvector // OUT: vector of six ellipsoid coefficients
        ) {
    double x, y, z; /* direction vector projections	*/
    double a[MAXROT][6]; /* A matrix 			*/
    double at[6][MAXROT]; /* A transpose 	(At)		*/
    //double	ata[6][6];				/* A transpose * A  (AtA)	*/
    //double	atai[6][6];				/* inverse of AtA  (AtA)^-1	*/
    double* ata;
    double* atai;
    double c[6][MAXROT]; /* (AtA)^-1 * At		*/
    double b[MAXROT]; // vector of 1/(MIL)^2		

    /*double	bmn, bfitmn;		// mean of measured and fit b	
    double	bfit[MAXROT];			// b vector based on fit	
    double	expdif, expfit, expsqr;	// stats for measured data	
    double	fitsqr, fitdif;			// stats for fit data		*/

    int i, j, k; // loop variables 		
    int degree = 6; // degree of vectors (6)	


    // Allocate memory:
    ata = (double*) calloc(6 * 6, sizeof (double));
    atai = (double*) calloc(6 * 6, sizeof (double));

    // Initialize the matrices:
    /*#pragma omp parallel for private( j ) 
    for ( i = 0; i < degree; i++ ) 
            for ( j = 0; j < degree; j++ ) 
            {
                    ata[i][j]  = 0.0;
                    atai[i][j] = 0.0;
            }*/

#pragma omp parallel for private( i )
    for (k = 0; k < MAXROT; k++)
        for (i = 0; i < degree; i++)
            c[i][k] = 0.0;

    for (i = 0; i < degree; i++)
        xvector[i] = 0.0;


    // Determine A and b matrices of Ax = b:
#pragma omp parallel for private( x, y, z ) 
    for (i = 0; i < MAXROT; i++) {
        x = cos(rot_phi[i]) * sin(rot_theta[i]);
        y = sin(rot_phi[i]) * sin(rot_theta[i]);
        z = cos(rot_theta[i]);

        a[i][0] = x * x;
        a[i][1] = y * y;
        a[i][2] = z * z;
        a[i][3] = x * y;
        a[i][4] = y * z;
        a[i][5] = x * z;

        b[i] = 1.0 / (*(mil + i) * *(mil + i));
    }

    // Determine At of AtAx = Atb:
#pragma omp parallel for
    for (i = 0; i < MAXROT; i++) {
        at[0][i] = a[i][0];
        at[1][i] = a[i][1];
        at[2][i] = a[i][2];
        at[3][i] = a[i][3];
        at[4][i] = a[i][4];
        at[5][i] = a[i][5];
    }

    // Multiply At * A:
#pragma omp parallel for private( i, j ) 
    for (k = 0; k < MAXROT; k++)
        for (j = 0; j < degree; j++)
            for (i = 0; i < degree; i++)
                //ata[i][j] += at[i][k]*a[k][j];
                ata[ FOO_I(i, j, 6) ] += at[i][k] * a[k][j];

    // Find inverse of AtA:
    linInverseMatrix(degree, degree, ata, atai);

    // Multiply (AtA)^-1 * At :
#pragma omp parallel for private( i, j )
    for (k = 0; k < MAXROT; k++)
        for (j = 0; j < degree; j++)
            for (i = 0; i < degree; i++)
                //c[i][k] += atai[i][j]*at[j][k];
                c[i][k] += atai[ FOO_I(i, j, 6) ] * at[j][k];

    // Multiply (AtA)^-1 * At * b:
#pragma omp parallel for private( i )
    for (k = 0; k < MAXROT; k++)
        for (i = 0; i < degree; i++)
            xvector[i] += c[i][k] * b[k];

    // Release resources:
    if (ata != NULL) free(ata);
    if (atai != NULL) free(atai);

    /*
    // Determine the goodness of fit.  Begin by calculating the sum of the
    // 	square of the residuals. 				
    sqres = 0.0;

    for(i=0;i<MAXROT;i++) {
            bfit[i] = xvector[0]*a[i][0] + xvector[1]*a[i][1] + 
                      xvector[2]*a[i][2] + xvector[3]*a[i][3] + 
                      xvector[4]*a[i][4] + xvector[5]*a[i][5];

            sqres += pow((b[i] - bfit[i]), 2.0);
    }

    //  Find the mean of b[] and bfit[] 
    bmn = 0.0;
    bfitmn = 0.0;

    for(i=0; i<MAXROT; i++) {
            bmn += b[i];
            bfitmn += bfit[i];
    }

    bmn /= (double)MAXROT;
    bfitmn /= (double)MAXROT;

    // Find the variance of the error in the ellipsoid fit: 
    variance = sqres / (double)(MAXROT - 6 - 1);
    stdev	 = sqrt(variance);

    //  Calculate the correlation coefficients: 
    expfit=0.0;
    expsqr=0.0;
    fitsqr=0.0;

    for ( i = 0; i < MAXROT; i++ ) 
    {
            expdif = b[i] - bmn;
            fitdif = bfit[i] - bfitmn;
            expfit += expdif * fitdif;
            expsqr += expdif * expdif;
            fitsqr += fitdif * fitdif;
    }

    // Avoid division by zero: 
    if(( expsqr < 0.0) || (fitsqr < 0.0) ) 
            corr = 1.0;  
    else	
            corr = expfit / sqrt(expsqr*fitsqr);*/
}


/************************************************************************/
/* Function: _getDegreesOfAnisotropy                                    */
/*                                                                      */
/*	This function defines the MIL tensor based on the ellipsoid fit	    */
/*	and determines the principal eigenvectors and eigenvalues, 	        */
/*	which are used to define the principal MIL vector orientations 	    */
/*	and magnitudes, respectively.  The MIL tensor assumes material	    */
/*	orthotropy and is defined as (Harrigan and Mann, JMatSci 	        */
/*	19:761-767, 1984):						                            */
/*									                                    */
/*				    | ae	  de/2	  ee/2 |		                    */
/*			    M = | de/2	  be	  fe/2 |		                    */
/*				    | ee/2	  fe/2	  ce   |		                    */
/*		where ae, be, ..,fe are the ellipsoid coefficients	            */
/*									                                    */
/*	The orientations of the principal MILs are defined by spherical     */
/*	angles,	a and b, where a and b correspond to theta and phi,	        */
/*	respectively.  The magnitudes of the principal MILs are 	        */
/*	calculated as the inverse of the square root of the absolute 	    */
/*	values of the eigenvalues.					                        */
/*                                                                      */
/*	The three degrees of anisotropy are defined according to Benn,      */
/*  Journ Sed Res., 64(4): 910 (1994) as:   	     			        */
/*      I = t3/t1                                                       */
/*	    E = 1 - t2/t1   			                                    */
/*  where the isotropy index I measures the similarity of a fabric to   */
/*  a uniform distribution and varies between 0 (all observations       */
/*  confined to a single plane or axis) and 1 (perfect isotropy). The   */
/*  elongation index E measures the preferred orientation of a fabric   */
/*  in the V1/V2 plane and varies between 0 (no preferred orientation   */
/*  and 1 (a perfect preferred orientation with all observations        */
/*  parallel.                                                           */
/*									                                    */

/************************************************************************/
int _getDegreesOfAnisotropy(
        double* coeffs, // A_IN: vector of six ellipsoid coefficients	
        double* degI, // OUT: isotropy index I
        double* degE, // OUT: elongation index E
        int (*wr_log)(const char*, ...)
        ) {
    double eps = 1.0E-6; /* convergence criterion 		*/

    double prin1, prin2, prin3; /* principal MIL vectors		*/
    double min_prin, med_prin, max_prin;
    double orient1a, orient1b; /* orientation of 1st principal		*/
    double orient2a, orient2b; /* orientation of 2nd principal		*/
    double orient3a, orient3b; /* orientation of 3rd principal		*/
    double prinmn; /* mean principal value of MIL		*/

    double eigvals[9];
    double eigvecs[9];

    int dim = 3; /* dimension of evalue and evector v_array*/
    int order = 3; /* order of matrices 			*/
    int numiters = 0; /* number of iterations required	*/
    int maxiters = 50; /* maximum number of iterations allowed	*/
    int status; /* status returned by eigJacobi		*/

    // Return ellipse coefficients through external variables: 
    double ae = coeffs[0];
    double be = coeffs[1];
    double ce = coeffs[2];
    double de = coeffs[3];
    double ee = coeffs[4];
    double fe = coeffs[5];

    /* Create the MIL tensor */
    *(eigvals + 0) = ae;
    *(eigvals + 1) = de / 2.0;
    *(eigvals + 2) = fe / 2.0;
    *(eigvals + 3) = de / 2.0;
    *(eigvals + 4) = be;
    *(eigvals + 5) = ee / 2.0;
    *(eigvals + 6) = fe / 2.0;
    *(eigvals + 7) = ee / 2.0;
    *(eigvals + 8) = ce;

    /* Calculate the eigenvalues and eigenvectors of the matrix/tensor */
    status = eigJacobi(dim, eigvals, eigvecs, order, eps, &numiters, maxiters);
    if (status == 1) return P3D_MEM_ERROR;

    /* Determine the angles of orientation for the structure */
    orient1a = atan2(*(eigvecs + 6), sqrt(pow(*(eigvecs + 0), 2.0) + pow(*(eigvecs + 3), 2.0)));
    orient1b = atan2(*(eigvecs + 3), *(eigvecs + 0));

    orient2a = atan2(*(eigvecs + 7), sqrt(pow(*(eigvecs + 1), 2.0) + pow(*(eigvecs + 4), 2.0)));
    orient2b = atan2(*(eigvecs + 4), *(eigvecs + 1));

    orient3a = atan2(*(eigvecs + 8), sqrt(pow(*(eigvecs + 2), 2.0) + pow(*(eigvecs + 5), 2.0)));
    orient3b = atan2(*(eigvecs + 5), *(eigvecs + 2));

    /* Convert the angles from radians -> degrees */
    orient1a *= (180.0 / M_PI);
    if (orient1a < 0.0) orient1a += 180.0;
    orient1b *= (180.0 / M_PI);
    if (orient1b < 0.0) orient1b += 180.0;
    orient2a *= (180.0 / M_PI);
    if (orient2a < 0.0) orient2a += 180.0;
    orient2b *= (180.0 / M_PI);
    if (orient2b < 0.0) orient2b += 180.0;
    orient3a *= (180.0 / M_PI);
    if (orient3a < 0.0) orient3a += 180.0;
    orient3b *= (180.0 / M_PI);
    if (orient3b < 0.0) orient3b += 180.0;

    /* Determine the principal values */
    prin1 = 1.000 / sqrt(fabs(*(eigvals + 0)));
    prin2 = 1.000 / sqrt(fabs(*(eigvals + 4)));
    prin3 = 1.000 / sqrt(fabs(*(eigvals + 8)));

    /* Normalize principal values: */
    prinmn = prin1 + prin2 + prin3;

    prin1 = prin1 / prinmn;
    prin2 = prin2 / prinmn;
    prin3 = prin3 / prinmn;

    // Find minimum (t1), medium (t2) and maximum (t3) normalized principal value:
    min_prin = MIN(MIN(prin1, prin2), prin3);
    max_prin = MAX(MAX(prin1, prin2), prin3);

    if ((prin1 != min_prin) && (prin1 != max_prin)) med_prin = prin1;
    if ((prin2 != min_prin) && (prin2 != max_prin)) med_prin = prin2;
    if ((prin3 != min_prin) && (prin3 != max_prin)) med_prin = prin3;

    // Compute degrees of anisotropy:
    *degI = (min_prin / max_prin); // I = t3 / t1
    *degE = 1 - (med_prin / max_prin); // E = 1 - t2 / t1	
    // DA = t1 / t3

    if (wr_log != NULL) {
        wr_log("\tPore3D - Anisotropy analysis verbose mode: ");
        wr_log("\t\t Principal MILs: %.4lf\t%.4lf\t%.4lf", prin1, prin2, prin3);
        wr_log("\t\t Orientation 1 (theta,phi): (%.2lf, %.2lf)", orient1a, orient1b);
        wr_log("\t\t Orientation 2 (theta,phi): (%.2lf, %.2lf)", orient2a, orient2b);
        wr_log("\t\t Orientation 3 (theta,phi): (%.2lf, %.2lf)", orient3a, orient3b);
        wr_log("\t");
        wr_log("\t\t Ellipsoid fit parameters:");
        wr_log("\t\t\t ae = %.4lf\tbe = %.4lf\tce = %.4lf", ae, be, ce);
        wr_log("\t\t\t de = %.4lf\tee = %.4lf\tfe = %.4lf", de, ee, fe);
        wr_log("\t\t Eigenvalues:");
        wr_log("\t\t \t%.4lf\t%.4lf\t%.4lf", *(eigvals + 0), *(eigvals + 4), *(eigvals + 8));
        wr_log("\t\t Eigenvectors:");
        wr_log("\t\t \tE-vector 1: (%.4lf,%.4lf,%.4lf)", *(eigvecs + 0), *(eigvecs + 3), *(eigvecs + 6));
        wr_log("\t\t \tE-vector 2: (%.4lf,%.4lf,%.4lf)", *(eigvecs + 1), *(eigvecs + 4), *(eigvecs + 7));
        wr_log("\t\t \tE-vector 3: (%.4lf,%.4lf,%.4lf)", *(eigvecs + 2), *(eigvecs + 5), *(eigvecs + 8));
    }

    // Return success:
    return P3D_SUCCESS;
}

int p3dAnisotropyAnalysis(
        unsigned char* in_im,
        unsigned char* msk_im,
        struct AnisotropyStats* out_stats,
        const int dimx,
        const int dimy,
        const int dimz,
        const double voxelsize,
        const int verbose,
        int (*wr_log)(const char*, ...)
        ) {
    double* rot_theta = (double *) malloc(MAXROT * sizeof (double));
    double* rot_phi = (double *) malloc(MAXROT * sizeof (double));
    double* mil = (double *) malloc(MAXROT * sizeof (double));
    double* coeffs = (double *) malloc(6 * sizeof (double));

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dAnisotropyAnalysis");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Performing anisotropy analysis...");
        wr_log("\tAdopted voxelsize: %0.6f mm.", voxelsize);
    }



    // Call routine to determine MILs based on random rotations: 
    _getMILs(in_im, msk_im, rot_theta, rot_phi, mil, dimx, dimy, dimz, voxelsize);

    /* Call routine to determine the best fit ellipsoid equation for data	*/
    _MILs2fit(rot_theta, rot_phi, mil, coeffs);

    /* Call routine to calculate the properties of the ellipsoid tensor  */
    if (verbose == P3D_TRUE)
        _getDegreesOfAnisotropy(coeffs, &(out_stats->I), &(out_stats->E), wr_log);
    else
        _getDegreesOfAnisotropy(coeffs, &(out_stats->I), &(out_stats->E), NULL);

    if (wr_log != NULL) {
        wr_log("\t----");
        wr_log("\tIsotropy index (I): %0.3f [-].", out_stats->I);
        wr_log("\tElongation index (E): %0.3f [-].", out_stats->E);
    }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Anisotropy analysis computed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }


    // Release resources:
    free(rot_theta);
    free(rot_phi);
    free(mil);
    free(coeffs);

    // Return OK:
    return P3D_SUCCESS;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

