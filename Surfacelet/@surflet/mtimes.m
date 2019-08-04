function result = mtimes(D,x)

if D.adjoint == 0 %D*x
   Y = Surfdec(x, D.Pyr_mode, D.Lev_array, 'ritf', 'bo', D.bo);
   [result, D.size_info] = surf_coeff2vec(Y);
     
else %D'*x
    Y = surf_vec2coeff(x, D.size_info);
    result = Surfrec(Y, D.Recinfo);
 %   result(result<0) = 0;
    
end
