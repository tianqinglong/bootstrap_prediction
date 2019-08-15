#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "likelihood_fn.h"
#include "likelihood_fn.c"

SEXP findMLE(SEXP fTime, SEXP cTime)
{
	int i;
	double out[2];
	double *fTime_p = REAL(fTime);
	double *cTime_p = REAL(cTime);
	SEXP Rout;

	Rout = PROTECT(allocVector(REALSXP, 2));

	int fail_len = length(fTime);
	int censor_len = length(cTime);

	double fail_times[fail_len], censor_times[censor_len];

	for(i=0; i<fail_len; i++)
	{
		fail_times[i] = fTime_p[i];
	}

	for(i=0; i<censor_len; i++)
	{
		censor_times[i] = cTime_p[i];
	}

	int status;
	int iter = 0, max_iter = 100;

	gsl_set_error_handler_off();

	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 2;
	double x_lo = 0.01, x_hi = 10;
	gsl_function F;
	struct Fail_Censor_Time time_data = {fail_len, censor_len, fail_times, censor_times};

	F.function = &derivative_equation;
	F.params = &time_data;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	// printf("using %s method\n",
	// 		gsl_root_fsolver_name (s));

	// printf("%5s [%9s, %9s] %9s\n", 
	// 		"iter", "lower", "upper", "root");

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);

		if (status) {break;}

		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

		// if (status == GSL_SUCCESS)
		// 	printf("Converged:\n");

		// printf("%5d [%.7f, %.7f] %.7f\n", 
		// 		iter, x_lo, x_hi,
		// 		r);
	}
	while ( status == GSL_CONTINUE && iter < max_iter );

	gsl_root_fsolver_free (s);

	if (status == GSL_SUCCESS) {
		int i;
		double lambda = 0;

		for (i=0; i<fail_len; i++)
		{
			lambda += pow (fail_times[i], r);
		}

		for (i=0; i<censor_len; i++)
		{
			lambda += pow (censor_times[i], r);
		}

		lambda = lambda/fail_len;
		lambda = pow (lambda, 1/r);

		out[0] = r;
		out[1] = lambda;
	}
	else {
		out[0] = -1;
		out[1] = -1;
	}

	REAL(Rout)[0] = out[0];
	REAL(Rout)[1] = out[1];

	UNPROTECT(1);
	return Rout;
}