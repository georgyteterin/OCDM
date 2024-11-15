function dFnT_output = DFnT(input_signal, varargin)
%Muneeb Ahmad
%Supported by: BK21FOUR MERIT
%WENS, KIT
    sz = size(input_signal);
    % Ensure the input_signal is a column vector
    input_signal = input_signal(:);
    
    % Determine the value of N
    if nargin == 1
        N = length(input_signal);
    else
        N = varargin{1};
    end
    
    % Precompute constants
    sqrtN = sqrt(N);
    exp_neg_j_pi_over_4 = exp(-1j * pi / 4);
    pi_over_N = pi / N;

    % Constructing the DFnT matrix
    dFnT_matrix = zeros(N, length(input_signal));
    for m = 1:N
        for n = 1:length(input_signal)
            if mod(N, 2) == 0  % For even N
                phase = pi_over_N * (m - n)^2;
            else  % For odd N
                phase = pi_over_N * (m + 0.5 - n)^2;
            end
            dFnT_matrix(m, n) = (1 / sqrtN) * exp_neg_j_pi_over_4 * exp(1j * phase);
        end
    end

    % Applying the DFnT
    tmp = dFnT_matrix * input_signal;
    dFnT_output = reshape(tmp, sz);
    clear tmp;
end