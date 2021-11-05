function j = j_gen(sigma_vec,alpha, beta, g_vector)
%J_GEN defines the cost function

left = sum(sigma_vec^2);
mid = alpha*g_vector;
right = beta*trap_int(dt, g_vector);

j = left + mid + right;
end

