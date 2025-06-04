%% 3.piemçrs. M-fails - st_fun_prob3
function dy_dx = st_fun_prob3(x,y)
   dy_dx = zeros(3,1);
   dy_dx(1) = y(2).*y(3);
   dy_dx(2) = -y(1).*y(3);
   dy_dx(3) = -0.5*y(1).*y(2);
end
