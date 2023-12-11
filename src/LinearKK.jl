"""
    cut_inductance(measurements,frequencies)

Removes inductive loops in the beginning or end of impedance spectra. 

Inputs : a set of complex-valued impedance values and their corresponding frequencies.
Outputs: the non-inductive parts of the impedance measurements and frequencies.
 """
function cut_inductance(measurements,frequencies)
    positive_index = findfirst(x -> x < 0, imag(measurements))
    if !isnothing(positive_index)
        measurements = measurements[positive_index:end]
        frequencies = frequencies[positive_index:end]
    end
    negative_index = findfirst(x -> x > 0, imag(measurements))
    if !isnothing(negative_index)
        measurements = measurements[1:negative_index]
        frequencies = frequencies[1:negative_index]
    end
    return measurements,frequencies
end

function voigt_model(parameters,ω)
    @assert !iseven(length(parameters))
    R_ohm = parameters[1]
    Vep = parameters[2:end]
    total_impedance = R_ohm + 0*im
    for i in 1:2:length(Vep)
        total_impedance += (Vep[i])/(1 + im*ω*Vep[i+1])
    end
    return total_impedance
end

function voigt_model_lin(resistances,ω,taus)
    R_ohm = resistances[1]
    voigt_resistances = resistances[2:end]
    M = length(voigt_resistances)
    total_impedance = R_ohm + 0*im
    for (R,tau) in zip(voigt_resistances,taus)
        total_impedance += (R)/(1 + im*ω*tau)
    end
    return total_impedance
end

voigt_model_lin(resistances,ωs::Vector,taus) = [voigt_model_lin(resistances,ω,taus) for ω in ωs]

function generate_time_constants(ω,M)
    @assert M > 3
    τ_min = 1/(maximum(ω))
    τ_max = 1/(minimum(ω))
    τs = [10^(log10(τ_min) + ((k - 1)/(M - 1))*log10(τ_max/τ_min)) for k in 2:M-1]
    return vcat(τ_min,τs,τ_max)
end

function objective_real(measurements,model_estimates)
    return sum(((real(measurements) .- real(model_estimates)) ./ abs.(measurements)).^2)
end

function objective_imag(measurements,model_estimates)
    return sum(((imag(measurements) .- imag(model_estimates)) ./ abs.(measurements)).^2)
end

function objective_all(measurements,model_estimates)
    return objective_real(measurements,model_estimates) + objective_imag(measurements,model_estimates)
end

function ECM_residuals(measurements,model_estimates)
    ΔRe = ((real(measurements) .- real(model_estimates)) ./ abs.(measurements))
    ΔIm = ((imag(measurements) .- imag(model_estimates)) ./ abs.(measurements))
    return ΔRe, ΔIm
end

function χ2_im(measurements,model_estimates)
    return (1/length(measurements))*sum(((imag(measurements) .- imag(model_estimates)) ./(abs.(measurements))).^2)
end

function fit_ECM(measurements,ωs,M)
    taus = generate_time_constants(ωs,M)
    A_real = zeros(length(ωs),M+1)
    A_real[:,1] .= 1
    A_imag = zeros(length(ωs),M+1)
    A_imag[:,1] .= 0
    for i in eachindex(ωs)
        for j in 1:M
            A_real[i,j+1] = real((1)/(1 + im*ωs[i]*taus[j]))
            A_imag[i,j+1] = imag((1)/(1 + im*ωs[i]*taus[j]))
        end
    end
    Z_real  = real(measurements)
    Z_imag  = imag(measurements)
    resistances_real = pinv(A_real'*A_real)*A_real'*Z_real
    resistances_imag = pinv(A_imag'*A_imag)*A_imag'*Z_imag
    return A_real * resistances_real .+ (A_imag * resistances_imag) .* im
end

function get_resistances(measurements,ωs,M)
    taus = generate_time_constants(ωs,M)
    A_real = zeros(length(ωs),M+1)
    A_real[:,1] .= 1
    for i in eachindex(ωs)
        for j in 1:M
            A_real[i,j+1] = real((1)/(1 + im*ωs[i]*taus[j]))
        end
    end  
    Z_real  = real(measurements)
    resistances_real = pinv(A_real'*A_real)*A_real'*Z_real
    return resistances_real
end

function compute_μ(resistances)
   R_negatives = count(x->sign(x) == -1, resistances)
   R_positives = length(resistances) - R_negatives + eps(Float64)
   μ = 1 - (R_negatives/R_positives)
   return μ
end

function find_optimal_M(measurements,ωs,c = 0.50, max_M = 50)
    @assert 200 > max_M > 4
    @assert 0 < c < 1
    μs = [compute_μ(get_resistances(measurements,ωs,M_)) for M_ in 4:max_M]
    findlast(x-> x > 0.5, μs)
    M_opt = findlast(x-> x > 0.5, mus)
    return M_opt
end

"""
    linearKK(measurements,frequencies;c = 0.5,max_M = 80)

Conduct a linear Kramers-Kronig evaluation for the validation of impedance spectroscopy measurements.
The essential inputs are a set of frequencies and the impedance measurements conducted at those frequencies.
There are also two of keyword arguments to fine-tune the calculation. 

## Keyword arguments
- `c::Float64=0.5`: Hyperparameter for automatic number of RC-elements determination. Lower values of c promote a higher number of RC elements. 
- `max_M::Int64=80`: Hyperparameter limiting the possible number of RC-elements of the Voigt circuit model.
 """
function linearKK(measurements,frequencies;c = 0.5,max_M = 80) #Alternatively: path to csv file.
    #output: real and imaginary residuals
    measurements,frequencies = cut_inductance(measurements,frequencies) #remove inductive parts.
    ωs = 2π*frequencies
    M = find_optimal_M(measurements,ωs,c,max_M)
    @info("$(M) RC-units were used to fit the measurements")
    fitted_spectrum = fit_ECM(measurements,ωs,M)
    realres, imres = ECM_residuals(measurements,fitted_spectrum)
    return  measurements, frequencies, fitted_spectrum, realres, imres
end

"""
    linearKK(path;c = 0.5,max_M = 80)

Conduct a linear Kramers-Kronig evaluation for the validation of impedance spectroscopy measurements.
The essential inputs are the path to a set of frequencies and the impedance measurements conducted at those frequencies.
There are also two of keyword arguments to fine-tune the calculation. 

## Keyword arguments
- `c::Float64=0.5`: Hyperparameter for automatic number of RC-elements determination. Lower values of c promote a higher number of RC elements. 
- `max_M::Int64=80`: Hyperparameter limiting the possible number of RC-elements of the Voigt circuit model.
 """
function linearKK(path,c = 0.5,max_M = 80) #Alternatively: path to csv file.  
    df = CSV.read(path,DataFrame)
    frequencies = df[:,1]
    measurements = df[:,2] .+ im .* df[:,3]
    measurements,frequencies = cut_inductance(measurements,frequencies) #remove inductive parts.
    ωs = 2π*frequencies
    M = find_optimal_M(measurements,ωs,c,max_M)
    @info("$(M) RC-units were used to fit the measurements")
    fitted_spectrum = fit_ECM(measurements,ωs,M)
    realres, imres = ECM_residuals(measurements,fitted_spectrum)
    return  measurements, frequencies, fitted_spectrum, realres, imres
end

"""
    summary_plot(measurements, frequencies, fitted_spectrum, realres, imres)

Visually evaluate the results of a linear Kramers-Kronig data validation.

The inputs are the outputs of the linearKK function, being the processed measurements, frequencies, Voigt-circuit fitted impedance spectrum,
the real residuals and the imaginary residuals. 

The output is a composite figure displaying the results of the linear Kramers-Kronig analysis.
 """
function summary_plot(measurements, frequencies, fitted_spectrum, realres, imres)
    # Spectrum fit
    specfit = scatter(real(measurements),-imag(measurements), label = "measured",markershape=:circle,markercolor=:white)
    scatter!(real(fitted_spectrum),-imag(fitted_spectrum), label = "fitted",markershape=:cross,markersize=6)
    xlabel!("Re(Z) [Ω]")
    ylabel!("Im(Z) [Ω]")
    title!("Measured and ECM fitted spectra")
    # Residuals fit
    res_plot = scatter(frequencies,realres .*100,label = "real residual",xaxis = :log, color = "blue", title="Residuals")
    plot!(frequencies,realres .*100,label = nothing, color = "blue")
    scatter!(frequencies,imres .*100,label = "imag residual", color = "red")
    plot!(frequencies,imres .*100,label = nothing, color = "red")
    ylims!(-2,2)
    ylabel!("Residual [%]")
    xlabel!("Frequency [Hz]")
    title!("Residuals")
    real_specfit = scatter(frequencies,real(measurements), label = "measured",markershape=:circle,markercolor=:white, xaxis = :log)
    scatter!(frequencies,real(fitted_spectrum), label = "fitted",markershape=:cross,markersize=6,legend=false)
    xlabel!("Frequency [Hz]")
    ylabel!("Re(Z) [Ω]")
    imag_specfit = scatter(frequencies,-imag(measurements), label = "measured",markershape=:circle,markercolor=:white, xaxis = :log)
    scatter!(frequencies,-imag(fitted_spectrum), label = "fitted",markershape=:cross,markersize=6,legend=false)
    xlabel!("Frequency [Hz]")
    ylabel!("-Im(Z) [Ω]")
    l = @layout [a{0.5w} grid(2,1); b{0.3h}]
    summary_fig = plot(specfit,real_specfit,imag_specfit,res_plot, layout = l,size=(1000,800))
    return summary_fig
end

"""
    summary_plot(path,c = 0.5,max_M = 80)

Visually evaluate the results of a linear Kramers-Kronig data validation.

The essential inputs are the path to a set of frequencies and the impedance measurements conducted at those frequencies.
There are also two of keyword arguments to fine-tune the calculation.

The output is a composite figure displaying the results of the linear Kramers-Kronig analysis.

## Keyword arguments
- `c::Float64=0.5`: Hyperparameter for automatic number of RC-elements determination. Lower values of c promote a higher number of RC elements. 
- `max_M::Int64=80`: Hyperparameter limiting the possible number of RC-elements of the Voigt circuit model.
 """
function summary_plot(path,c = 0.5,max_M = 80)
    measurements, frequencies, fitted_spectrum, realres, imres = linearKK(path,c,max_M)
    return summary_plot(measurements, frequencies, fitted_spectrum, realres, imres)
end


"""
    save_valid_measurements(path,measurements,frequencies,threshold = 1, c = 0.5,max_M = 80)

    remove part of EIS measurements that doesn't comply with EIS validity criteria as evaluated by the linear Kramers-Kronig relations and
    save the resulting impedance spectrum.

    The essential inputs are the path to where the valid part of the measurements should be saved and a set of frequencies and the impedance measurements conducted at those frequencies.
    There are also two of keyword arguments to fine-tune the calculation.

    ## Keyword arguments
    - `c::Float64=0.5`: Hyperparameter for automatic number of RC-elements determination. Lower values of c promote a higher number of RC elements. 
    - `max_M::Int64=80`: Hyperparameter limiting the possible number of RC-elements of the Voigt circuit model.
    - `threshold::Float64=1`: Threshold percentage under which residuals are considered to be Kramers-Kronig compliant.
 """
function save_valid_measurements(path,measurements,frequencies,threshold = 1, c = 0.5,max_M = 80) #remove part of EIS measurements that doesn't comply with linKK
    measurements, frequencies, fitted_spectrum, realres, imres = linearKK(measurements,frequencies,c,max_M)
    thres = threshold*0.01
    invalid_inds = findall((realres .> thres) .| (imres .> thres))
    meas = measurements[setdiff(1:end, invalid_inds)]
    freq = frequencies[setdiff(1:end, invalid_inds)]
    df_valid = DataFrame(frequencies = freq, Real = real(meas), Imag = imag(meas))
    CSV.write(path,df_valid)
end

"""
    save_lin_kk_residuals(path,freqs,real_res,imag_res)

    Save the residuals of the linear Kramers-Kronig test to a CSV file.

    The essential inputs are the path to where the residuals ought to be saved, the frequency values, and the real and imaginary residuals
    obtained by the linearKK function.

 """
function save_lin_kk_residuals(path,freqs,real_res,imag_res)
    df = DataFrame(frequecies=freqs,real_residuals = real_res, imag_residuals=imag_res)
    CSV.write(path,df)
end