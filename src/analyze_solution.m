function analyze_solution(pars)
    fprintf('The original matrix:\n')
    disp(pars.A_original)
    fprintf('The reduced matrix:\n')
    disp(pars.A)
    already_handled = false(length(pars.analysis),1);
    for i = 1:length(pars.analysis)
        if ~already_handled(i)
            f1 = pars.analysis{i};
            fprintf('#####################\n');
            fprintf('I0 = [');
            fprintf(' %s ', string(int8(f1.I0')))
            fprintf(']\n');
            for j = i:length(pars.analysis)
                f2 = pars.analysis{j};
                if isequal(f1.I0, f2.I0)
                    already_handled(j) = true;
                    fprintf('Mode = %s\n', f2.text);
                    fprintf('Count = %d\n', f2.count);
                    fprintf('Solving system Dx=d with D with rank %d\n', rank(f2.D));
                    if isfield(f2, 'D')
                        disp(f2.D)
                    end
                end
            end
        end
    end
end