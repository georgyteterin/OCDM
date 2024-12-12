function dFnT_output = DFnTmtrx(len, varargin)
%Muneeb Ahmad
%Supported by: BK21FOUR MERIT
%WENS, KIT
    
    % Determine the value of N
    if nargin == 1
        N = len;
    else
        N = varargin{1};
    end
    
    % Precompute constants
    sqrtN = sqrt(N);
    exp_neg_j_pi_over_4 = exp(-1j * pi / 4);
    pi_over_N = pi / N;

    % Constructing the DFnT matrix
    dFnT_matrix = zeros(N, len);
    for m = 1:N
        for n = 1:len
            if mod(N, 2) == 0  % For even N
                phase = pi_over_N * (m - n)^2;
            else  % For odd N
                phase = pi_over_N * (m + 0.5 - n)^2;
            end
            dFnT_matrix(m, n) = (1 / sqrtN) * exp_neg_j_pi_over_4 * exp(1j * phase);
        end
    end

    dFnT_output = dFnT_matrix;
    clear tmp;
end