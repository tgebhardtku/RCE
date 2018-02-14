function J = fxn(x,a)
J = fdjac( @fn, x, a );
return
