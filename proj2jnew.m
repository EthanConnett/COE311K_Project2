%inputs:
%delta_vec=  ; a=0.0005; b=0.00005; c=0.001; T=50;
%g_vec and f_vec computed in main


function j_new=costFunctionNew(delta_vec, a, b, c, g_vec, f_vec, T);

toxic_effects=sum(delta_vec.^2);                %first term of Jnew
cellPopulation1=a*g_vec(end);                     %second term
cellPopulation2=b*trap_int(dt, g_vec);          %third term
cellPopulation3=c*trap_int(dt,f_vec);           %fourth term

j_new=toxic_effects+cellPopulation1+cellPopulation2+cellPopulation3;

end
