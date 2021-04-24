function burst_size = tradeoff_value(k)
    sigma = 50;
    if k>0
        burst_size = 100*exp(-k^2/sigma^2);
    elseif k == 0 
        burst_size = 120;
    else
        burst_size =0;
    end    
end

