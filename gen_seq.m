function [x, n] = gen_seq(type, n1, n2, varargin)
n = n1:n2;
switch lower(type)
    case 'impulse'
        no=varargin{1};
        A=varargin{2};
        x = A.*double(n == no);
    case 'step'
        no=varargin{1};
        A=varargin{2};
        x = A.*double(n >= no);
    case 'ramp'
        n0 = varargin{1};
        A = varargin{2};
        x = A.*(n - n0) .* (n >= n0);
    case 'rect'
        n0 = varargin{1};
        W = varargin{2};
        A = varargin{3};
        x = A.*double(n >= n0 & n <= n0 + W);
    case 'exp'
        n0 = varargin{1};
        A = varargin{2};
        N = varargin{3};
        wc=2*pi/N;
        x = A.*exp(1j*wc*(n-n0));
    case 'sin'
        N = varargin{3};
        wc=2*pi/N;
        x = varargin{2}.* sin(wc*(n - varargin{1}));
    case 'tri'
        A = varargin{2}; W = varargin{3};
        x = A * max(1 - abs(n) / W, 0);
    case 'sign'
        x = double(n>=0) - double(n<0);
    case 'sinc'
        N = varargin{3};
        wc=(2*pi)/N;
        if wc <= 0 || wc > pi
            error('gen_seq: sinc cutoff wc must be in (0, pi]');
        end
        n0 = varargin{1};
        A = varargin{2};
        m  = n - n0;
        x  = (wc/pi) * ones(size(n));        % limit value at m = 0
        nz = (m ~= 0);
        x(nz) = A.*sin(wc * m(nz)) ./ (pi * m(nz));
    otherwise
        error('gen_seq: unknown sequence type "%s"', type);
end
end