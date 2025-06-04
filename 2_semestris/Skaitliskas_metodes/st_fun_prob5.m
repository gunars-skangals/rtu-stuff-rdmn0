%% 5.piemçrs. M-fails - st_fun_prob5
function dy_dx = st_fun_prob5(x,y)
   dy_dx = zeros(3,1);
   dy_dx(1) = y(2);
   dy_dx(2) = y(3);
   dy_dx(3) = y(3)./x-2*y(2)./x.^2+2*y(1)./x.^3+1;
end
