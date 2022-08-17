function result = proj(x, theta_hat, theta_min, theta_max)
    %Parameter projection for updatting law
    if(x < 0) && ( theta_hat <= theta_min)
        result = 0;
    elseif (x > 0) && ( theta_hat >= theta_max)
            result = 0;
    else
        result = x;
    end
end