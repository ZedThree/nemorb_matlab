function [dy] = get_deriv(y,x)

fun = spline(x,y);
dfun= fnder(fun);
dy  = ppval(x,dfun);

