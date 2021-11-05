%% trap_int
% The trap_int function can be used to integrate a function using the
% trapezoidal rule

%Inputs
% dx - Spacing between points
% y - Vector to be integrated

function val = trap_int(dx,y)
    
    val = sum(y(2:end-1)) + (y(1) + y(end))/2;
    val = val*dx;
end

