%% J_GEN
% Inputs
% sigma - number of drug inputs
% alpha - 
% beta - 
% g_vector - 

function j = j_gen(delta_vec, alpha, beta, g_vec, dt)
%J_GEN defines the cost function
%
toxic = sum(delta_vec.^2);
tumor1 = alpha*g_vec(end);
tumor2 = beta*trap_int(dt, g_vec);

j = toxic + tumor1 + tumor2;
end

