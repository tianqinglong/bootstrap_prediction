double
derivative_equation (double beta, void *params)
{
	struct Fail_Censor_Time *p
		= (struct Fail_Censor_Time *) params;

	int num_f = p->num_of_failures;
	int num_c = p->num_of_censor;

	double *f_times = p->failure_times;
	double *c_times = p->censor_times;

	// Given beta, compute lambda
	int i;
	double lambda = 0, first_derivative;

	for (i=0; i<num_f; i++)
	{
		lambda += pow (f_times[i], beta);
	}

	for (i=0; i<num_c; i++)
	{
		lambda += pow (c_times[i], beta);
	}

	lambda = lambda/num_f;
	lambda = pow (lambda, 1/beta);

	first_derivative = num_f/beta - num_f*log(lambda);
	// compute the first derivative
	for (i=0; i<num_f; i++)
	{
		first_derivative += log ( f_times[i] );
		first_derivative -= pow ( f_times[i]/lambda, beta ) * ( log(f_times[i]) - log(lambda) );
	}

	for (i=0; i<num_c; i++)
	{
		first_derivative -= pow(c_times[i]/lambda, beta) * ( log(c_times[i]) - log(lambda) );
	}

	return first_derivative;
}