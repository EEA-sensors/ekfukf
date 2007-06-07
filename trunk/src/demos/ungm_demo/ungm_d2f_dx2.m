function d = ungm_d2f_dx2(x,param)
    d = 25*(-6*x+2*x.^5-4*x.^3)/((1+x.^2).^4);