#pragma once


void ympr_to_lame(double ym, double pr, double& mu, double& lam)
{
	lam = (ym *pr) / ((1.0 + pr)*(1.0 - 2.0*pr));
	mu = ym / (2.0 *(1.0 + pr));
}