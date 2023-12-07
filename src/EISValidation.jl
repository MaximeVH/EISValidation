module EISValidation

    export linearKK, summary_plot, save_valid_measurements

    using Distributions,LinearAlgebra 
    using Plots, CSV, DataFrames

    include("LinearKK.jl")
    include("THD.jl")
end