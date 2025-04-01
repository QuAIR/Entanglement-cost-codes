function [cvx_optval, X] = LRbound(rho, varargin)

if nargin<2
    omega = rho;
    d = sqrt(max(size(rho)))*[1,1];
elseif nargin==2
    if min(size(varargin{1})) == 1
        omega = rho;
        d = varargin{1} * ones(1,3-max(size(varargin{1})));
    else
        omega = varargin{1};
        d = sqrt(max(size(rho)))*[1,1];
    end
elseif nargin==3
    omega = varargin{1};
    d = varargin{2} * ones(1,3-max(size(varargin{2})));;
end

dd = prod(d);

rho = (rho+rho')/2; % to avoid numerical issues
omega = (omega+omega')/2;

cvx_begin sdp quiet
    cvx_solver sedumi
    cvx_precision best
    variable X(dd,dd) hermitian

    maximize trace(X*rho)

    -eye(dd) <= PartialTranspose(X,2,d) <= eye(dd)

    -trace(X*omega)*eye(dd) <= X <= trace(X*omega)*eye(dd)

cvx_end

cvx_optval = log2(cvx_optval);

end