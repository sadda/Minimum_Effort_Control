classdef Pars < handle
    properties
        A
        A_original
        m
        n
        U
        zero_columns
        multiples
        analysis
        analysis_update
        analysis_i
    end
    
    methods
        function self = Pars(A)
            s = get_u(A);
            s_fields = fieldnames(s);
            for i = 1:length(s_fields)
                self.(s_fields{i}) = s.(s_fields{i});
            end
        end

        function solution_part(self, I0, text, varargin)
            if self.analysis_update
                if isequal(self.analysis, [])
                    self.analysis_i = self.analysis_i + 1;
                    self.analysis = {self.create_new_struct(I0, text)};
                    self.update_self(1, varargin{:});
                else
                    found = false;
                    for i = 1:length(self.analysis)
                        if isequal(self.analysis{i}.I0, I0) && isequal(self.analysis{i}.text, text)
                            found = true;
                            break
                        end
                    end
                    if found
                        self.analysis_i = self.analysis_i + 1;
                        self.update_self(i, varargin{:});
                    else
                        self.analysis_i = self.analysis_i + 1;
                        self.analysis = [self.analysis; self.create_new_struct(I0, text)];
                        self.update_self(length(self.analysis), varargin{:});
                    end
                end
            end
        end

        function output = create_new_struct(~, I0, text)
            output = struct('I0', I0, 'text', text, 'count', 0, 'i', []);
        end

        function update_self(self, i, varargin)
            self.analysis{i}.count = self.analysis{i}.count + 1;
            self.analysis{i}.i = [self.analysis{i}.i; self.analysis_i];
            if length(varargin) >= 1
                n_constant = varargin{1};
                for j = 1:n_constant
                    self.analysis{i}.(varargin{2*j}) = varargin{2*j+1};
                end
                for j = n_constant+1:(length(varargin)-1)/2
                    if isfield(self.analysis{i}, varargin{2*j})
                        self.analysis{i}.(varargin{2*j}) = [self.analysis{i}.(varargin{2*j}); varargin{2*j+1}'];
                    else
                        self.analysis{i}.(varargin{2*j}) = [varargin{2*j+1}'];
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
                                if isfield(f2, 'idx')
                                    fprintf('Kept indices of D\n');
                                    disp(f2.idx)
                                end
                                if isfield(f2, 'idx')
                                    fprintf('Null space of D\n');
                                    disp(f2.direction)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

