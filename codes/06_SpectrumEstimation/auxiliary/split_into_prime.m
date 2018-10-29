function [y_matrix] = split_into_prime(y)

N = length(y);
primes = load('prime_numbers');
prime_numbers=primes.prime_numbers;

kk=1;
while and(N>=prime_numbers(kk),N<prime_numbers(end))
    kk=kk+1;
end

if kk==1
    N_prime=1;
elseif N<prime_numbers(end)
    N_prime=prime_numbers(kk-1);
elseif N==prime_numbers(end)
    N_prime=prime_numbers(end);
elseif N>prime_numbers(end)
    N_prime=prime_numbers(end);
    display 'make a longer list of prime numbers'
end  

if N==N_prime
y_matrix = y;
else
   y_matrix = [y(1:N_prime) y(end-N_prime+1:end)];
end
       
