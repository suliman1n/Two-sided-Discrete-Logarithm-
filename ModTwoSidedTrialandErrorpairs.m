function [t_1, t_2, alpha] = ModTwoSidedTrialandErrorpairs(U, A, B, max_t)
    t_1 = -1;
    t_2 = -1;
    
    % Set the size of T based on the size of U
    n = size(U, 1);
    T = TropId(n);
    
    for i = 0:max_t
        for j = 0:max_t
            % Check if all entries of T are shifted by the same amount from U
            alpha = U - T;
            if all(alpha == alpha(1))
                t_1 = i;
                t_2 = j;
                alpha=alpha(1);
                return;
            end
            T = TropMulti(T, B);
        end
        T = TropMulti(fastpowermaxplus(A, i + 1), TropId(n));
    end
end

%TropMulti(TropMatPower(A,100),TropMatPower(B,100))