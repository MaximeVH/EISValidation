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

function find_optimal_M(measurements,ωs,c = 0.85, max_M = 50)
    @assert max_M > 4
    @assert 0 < c < 1
    μ = 1
    M = 4
    while (μ > c) && (M < 50)
        M += 1
        μ =  compute_μ(get_resistances(measurements,ωs,M))
    end
    return M, μ
end


    