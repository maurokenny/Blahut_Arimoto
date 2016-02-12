function [C0 r] = BlahutArimotoConstraint(p,x,P,rho)
if ~isempty(find(p < 0))
 disp('Error: some entry in the input matrix is negative')
 C0 = 0; return;
end

% Check that the input matrix p does not have zero column
column_sum = sum(p);
if ~isempty(find(column_sum == 0))
 disp('Error: there is a zero column in the input matrix');
 C0 = 0; return;
end
% Check that the input matrix p does not have zero row
row_sum = sum(p,2);
if ~isempty(find(row_sum == 0))
 disp('Error: there is a zero row in the input matrix');
 C0= 0; return;
else
 p = diag(sum(p,2))^(-1) * p; % Make sure that the row sums are 1
end
[m n] = size(p);


% (1) Choose any r(x) such that 0 < r(x) < 1 and
% sum_x r(x) = 1.
r = ones(1,m)/m; % initial distribution for channel input
% (2) Initialize capacity C0, C?1.
C0 = 1; C1 = 0;
error_tolerance = 1e-9/m;
%error_tolerance = 1e-9;
error_tolerance2 = 1e-9/m;
%error_tolerance2 = 1e-9;
q = zeros(m,n);
while ((C0-C1) > error_tolerance)
    C1 = C0;
    for j = 1:n
        q(:,j) = r'.*p(:,j);
        q(:,j) = q(:,j)/sum(q(:,j));
    end
    %(3)
    C0 = 0;
    for i = 1:m
        for j = 1:n
            if r(i) > 0 && q(i,j) > 0
                C0 = C0 + r(i)*p(i,j)*log(q(i,j)/r(i));
            end
        end
    end
    C0 = C0/log(2); % Capacity in bits
    
    %Step(3)
    pro = zeros(m,1);
    for i = 1:m
        pro(i) = prod(q(i,:).^p(i,:));
    end
    break_cond = 1000;
    B0 = 0; B1 = 0; cont = 0;
    while (abs(B1-B0) > error_tolerance2 || cont == 0)
        if (cont > break_cond)
            break;
        end
        cont = cont + 1;
        B1 = B0;
        
        temp = exp(B1.*abs(x).^rho).*(1-abs(x).^rho./P).*pro';
        B0 = B1 - sum(temp)/sum(abs(x).^rho.*temp);
    end
    r = exp(B0.*abs(x).^rho).*pro'/sum(exp(B0.*abs(x).^rho).*pro');
end

