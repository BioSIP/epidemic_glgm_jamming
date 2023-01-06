function dC = logisticGrowth(t,C,r,p,K)
   
    dC = zeros(1,1);
    dC(1,1) = (r*(C(1,1).^p)*(1-(C(1,1)/K)));

end