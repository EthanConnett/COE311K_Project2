%% J_GEN
% Inputs
% sigma - vector containing durg volume fraction inputs
% alpha - 
% beta - 
% g_vector - tum growth vector

function j = j_gen(delta_vec, params, tum)
%J_GEN defines the cost function
[f, g, p, t] = tum(delta_vec);

toxic = sum(delta_vec.^2);
tumor1 = params.a*g(end);
tumor2 = params.b*trap_int(params.dt, g);

j = toxic + tumor1 + tumor2;
end

