function [errorAbs, errorRel] = test_sensitivity_value_function(name, sens, sensFD)
    % Compute the error and relative error of the sensitivity.
    % Inputs:
    %   name: name of the variable.
    %   sens: analytical sensitivity.
    %   sensFD: numerical sensitivity.
    % Outputs:
    %   sensErr: absolute error.
    %   sensErrRel: relative error.

    % Compute absolute and relative errors
    errorAbs = abs(sens - sensFD);
    errorRel = errorAbs ./ abs(sens);

    % Check if the sensitivity is a vector
    if any(size(sens) > 1)
        sensDisp = norm(sens);
        sensFDDisp = norm(sensFD);
        errorAbsDisp = norm(errorAbs);
        errorRelDisp = norm(errorRel);
    else
        sensDisp = sens;
        sensFDDisp = sensFD;
        errorAbsDisp = errorAbs;
        errorRelDisp = errorRel;
    end

    % Display
    fprintf('\t\t%7s, %12.4e, %12.4e, %12.4e, %12.4e\n', name, sensDisp, sensFDDisp, errorAbsDisp, errorRelDisp)
end
