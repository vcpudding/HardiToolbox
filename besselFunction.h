#ifndef BESSELFUNCTION_H
#define BESSELFUNCTION_H
#include<iostream>
#include<cmath>
/*
 * Returns the modified Bessel Function
 * I_0(x) of zeroth order first kind
 * for any real x.
 */

long double bessi0(long double x)
{
	long double ax,ans;
	long double y; //Accumulate polynomials in long doubleprecision.
		if ((ax=fabs(x)) <3.75) { //Polynomial fit.
			y=x/3.75;
			y*=y;
			ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
							+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
		}else {
			y=3.75/ax;
			ans=(expl(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
						+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
									+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
												+y*0.392377e-2))))))));
		}
	return ans;
}
/*
 * Returns the modified Bessel Function
 * I_1(x) of zeroth order first kind
 * for any real x.
 */
long double bessi1(long double x)
{
	long double ax,ans;
	long double y; //Accumulate polynomials in long double precision.
		if ((ax=fabs(x)) <3.75) { //Polynomial fit.
			y=x/3.75;
			y*=y;
			ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
								+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
		}else {
			y=3.75/ax;
			ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
						-y*0.420059e-2));
			ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
						+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
			ans *= (expl(ax)/sqrt(ax));
		}
	return x< 0.0 ? -ans : ans;
}
/*
 * Returns the modified Bessel Function
 * I_1(x)/I_0(x) of zeroth order first kind
 * for any real x.
 */
 long double bessiRatio(long double x)
 {
   return x>700?1:bessi1(x)/bessi0(x);
 }
#endif
