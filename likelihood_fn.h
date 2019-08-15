struct Fail_Censor_Time
{
	int num_of_failures;
	int num_of_censor;

	double *failure_times;
	double *censor_times;
};

double derivative_equation (double beta, void *params);