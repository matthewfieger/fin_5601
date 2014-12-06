% Hull 331, 335, 350

function [Price] = BSPrice(S,K,T,r,vol,q,IsCall)
	% Returns European Black Scholes Call & Put Price
	d_1 = (log(S/K)+(r-q+0.5*vol^2)*T)/(vol*sqrt(T));
	d_2 = d_1-vol*sqrt(T);
	if IsCall
		Delta = exp(-q*T)*normcdf(d_1);
		Price = S*Delta - K*exp(-r*T)*normcdf(d_2);
	else
		Delta = exp(-q*T)*normcdf(-d_1);
		Price = -S*Delta + K*exp(-r*T)*normcdf(-d_2);
	end