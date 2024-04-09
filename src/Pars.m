classdef Pars < handle
    properties
        A
        A_original
        U
        zero_columns
        multiples
        analysis
        analysis_counter
    end

    methods
        function self = Pars(A)
            s = get_u(A);
            s_fields = fieldnames(s);
            for i = 1:length(s_fields)
                self.(s_fields{i}) = s.(s_fields{i});
            end
            self.analysis = {};
            self.analysis_counter = 0;
        end

        function x_new = expand_solution(self, x)
            % Expand the solution to the original space
            idx = find(~self.zero_columns);
            idx = idx(setdiff(1:length(idx), self.multiples(:,2)));
            x_new = zeros(size(self.A_original,2),1);
            x_new(idx) = x;
            % Distribute the values into the multiples columns
            for k = 1:size(self.multiples,1)
                x_new(self.multiples(k,2)) = x_new(self.multiples(k,1)) * sign(self.multiples(k,3));
            end
        end

        function x_new = shrink_solution(self, x)
            % Shrink the solution to the reduced space
            idx = find(~self.zero_columns);
            idx = idx(setdiff(1:length(idx), self.multiples(:,2)));
            x_new = x(idx);
        end

        function i_found = analysis_get_i(self, I0, text)
            i_found = -1;
            for i = 1:length(self.analysis)
                if isequal(self.analysis{i}.I0, I0) && isequal(self.analysis{i}.text, text)
                    i_found = i;
                    break
                end
            end
        end

        function increase_counter(self, I0, text, varargin)
            i = self.analysis_get_i(I0, text);
            if i < 0
                self.analysis = [self.analysis; self.create_new_struct(I0, text)];
                i = length(self.analysis);
            end
            self.analysis_counter = self.analysis_counter + 1;
            self.update_analysis_i(i, varargin{:});
        end

        function output = create_new_struct(~, I0, text)
            output = struct('I0', I0, 'text', text, 'count', 0, 'i', []);
        end

        function update_analysis_i(self, i, varargin)
            self.analysis{i}.count = self.analysis{i}.count + 1;
            self.analysis{i}.i = [self.analysis{i}.i; self.analysis_counter];
            if length(varargin) >= 1
                n_constant = varargin{1};
                for j = 1:n_constant
                    self.analysis{i}.(varargin{2*j}) = varargin{2*j+1};
                end
                for j = n_constant+1:(length(varargin)-1)/2
                    if isfield(self.analysis{i}, varargin{2*j})
                        self.analysis{i}.(varargin{2*j}) = [self.analysis{i}.(varargin{2*j}); varargin{2*j+1}'];
                    else
                        self.analysis{i}.(varargin{2*j}) = varargin{2*j+1}';
                    end
                end
            end
        end

        function plot_s_min_s_max(self, x_scale, varargin)
            if nargin < 2
                x_scale = 1;
            end
            figure();
            hold on;
            idx = [];
            for i = 1:length(self.analysis)
                analysis_i = self.analysis{i};
                if isfield(analysis_i, 's_min')
                    idx = [idx, i];
                end
            end
            lines_color = lines(length(idx));
            for i = 1:length(lines_color)
                analysis_i = self.analysis{idx(i)};
                if isfield(analysis_i, 's_min')
                    j_start = 1;
                    j_end = 1;
                    while j_start <= length(analysis_i.i)
                        if j_end < length(analysis_i.i) && analysis_i.i(j_end + 1) == analysis_i.i(j_end) + 1
                            j_end = j_end + 1;
                        else
                            jdx = j_start:j_end;
                            if length(jdx) > 1
                                plot(x_scale*analysis_i.i(jdx), [analysis_i.s_min(jdx), analysis_i.s_max(jdx)], 'color', lines_color(i,:), varargin{:});
                                plot(x_scale*analysis_i.i(jdx), 0.5*(analysis_i.s_min(jdx)+analysis_i.s_max(jdx)), '--', 'color', lines_color(i,:));
                            else
                                scatter(x_scale*analysis_i.i(jdx), [analysis_i.s_min(jdx), analysis_i.s_max(jdx)], [], lines_color(i,:), '.');
                                scatter(x_scale*analysis_i.i(jdx), 0.5*(analysis_i.s_min(jdx)+analysis_i.s_max(jdx)), [], lines_color(i,:), '.');
                            end
                            j_start = j_end + 1;
                            j_end = j_start;
                        end
                    end
                end
            end
        end

        function analyze_solution(self, full_print)
            if nargin < 2
                full_print = false;
            end
            fprintf('The original matrix:\n')
            disp(self.A_original)
            fprintf('The reduced matrix:\n')
            disp(self.A)
            already_handled = false(length(self.analysis),1);
            for i = 1:length(self.analysis)
                if ~already_handled(i)
                    f1 = self.analysis{i};
                    fprintf('#####################\n');
                    fprintf('I0 = [');
                    fprintf(' %s ', string(int8(f1.I0')))
                    fprintf(']\n');
                    for j = i:length(self.analysis)
                        f2 = self.analysis{j};
                        if isequal(f1.I0, f2.I0)
                            already_handled(j) = true;
                            fprintf('Mode = %s\n', f2.text);
                            fprintf('Count = %d\n', f2.count);
                            if isfield(f2, 'D')
                                fprintf('Solving system Dx=d with D with rank %d\n', rank(f2.D));
                                disp(f2.D)
                            end
                            if full_print
                                if isfield(f2, 'D_pse')
                                    fprintf('Pseudoinverse of D\n');
                                    disp(f2.D_pse)
                                end
                                if isfield(f2, 'D_idx')
                                    fprintf('Kept indices of D\n');
                                    disp(f2.D_idx)
                                end
                                if isfield(f2, 'D_v')
                                    fprintf('Null space of D\n');
                                    disp(f2.D_v)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

