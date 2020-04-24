

# compute the chi^2 or the reduced chi^2

function chi2(iexp::AMat, fitdata::AMat; option::AS="log")

    qmax = min(iexp[end, 1], fitdata[end, 1]);
    fitidx = findall(fitdata[:, 1] .<= qmax);
    expidx = findall(iexp[:, 1] .<= qmax);
    ifit = Spline1D(fitdata[fitidx, 1], fitdata[fitidx, 2])(iexp[expidx,1]);

    χ2 = -1.0;
    if option == "log"
        σ = iexp[expidx, 3] ./ iexp[expidx, 2] / log(10);
        χ2 = mean(((log10.(iexp[expidx, 2]) - log10.(ifit)) ./ σ) .^ 2);
    elseif option == "lin"
        χ2 = mean(((iexp[expidx, 2] - ifit) ./ iexp[expidx, 3]) .^ 2);
    elseif option == "abs"
        χ2 = mean(abs.(iexp[expidx, 2] - ifit));
    else
        error("--- χ2: Unknown option = $option ... ");
    end

    return χ2;
end
