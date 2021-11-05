%% J_GEN
% Inputs
% sigma - number of drug inputs
% 

function j = j_gen(sigma_vec,alpha, beta, g_vector)
%J_GEN defines the cost function

toxic = sum(sigma_vec^2);
tumor1 = alpha*g_vector;
tumor2 = beta*trap_int(dt, g_vector);

j = toxic + tumor1 + tumor2;
end

