using BiologicalOscillations, Statistics, DataFrames, JLD2, CSV

"""
This script extracts frequency and amplitude data for oscillating solutions obtained via the find oscillations algorithms.
The data is stored in a .csv file for each network. The input to the script is a directory with folders for each network.
It is assumed that networks are run in repetitions, so each network folder contains a folder for each repetition named 
with integers starting from 1.

Author: Franco Tavella
Date: 09/11/2023
"""

# User-defined parameters
folder_location = "E:/Project 2 - Tunability/ForkBiologicalOscillations.jl/results_analyses/T0_1_add/find_oscillations_result/T0/simulations"
save_path = "E:/Project 2 - Tunability/ForkBiologicalOscillations.jl/results_analyses/T0_1_add/find_oscillations_result/T0/new_result"

# Get list of network names and their repetitions
network_names = readdir(folder_location)
network_repetitions = [readdir(joinpath(folder_location, network)) for network in network_names]

for (idx, network_name) in enumerate(network_names)
    print("Processing $(network_name)...\n")
    batch_size = size(network_repetitions[idx], 1)
    total_samples = 0
    oscillatory_samples = 0
    features_df = DataFrame([Int64[], Int64[], Float64[],Float64[]], 
                            ["Batch", "Index", "Frequency", "Amplitude"])
    for batch in 1:batch_size
        fpath = "$(folder_location)/$(network_name)/$(batch)/$(network_name)_$(batch).jld2"
        data = load(fpath)

        samples = size(data["pin_result"]["parameter_sets"], 1)
        oscillatory_df = filter(row -> row["is_oscillatory"] == true, data["pin_result"]["simulation_result"])
        oscillatory_solutions = size(oscillatory_df, 1)

        # Calculate average frequency and amplitude across nodes
        amplitude_df = select(oscillatory_df, r"amplitude_*")
        frequency_df = select(oscillatory_df, r"frequency_*")
        avg_frequencies = mean.(eachrow(frequency_df))
        avg_amplitudes = mean.(eachrow(amplitude_df))

        parameter_index = oscillatory_df[!, "parameter_index"]
        batch = ones(Int64, oscillatory_solutions) .* batch
        partial_df = DataFrame([batch, parameter_index, avg_frequencies, avg_amplitudes], 
                               ["Batch", "Index", "Frequency", "Amplitude"])

        append!(features_df, partial_df)
        total_samples += samples
        oscillatory_samples += oscillatory_solutions
    end
    # Save feature dataframe as a csv file
    features_df[!, "Batch"] = convert(Vector{Int64}, features_df[!, "Batch"])
    features_df[!, "Index"] = convert(Vector{Int64}, features_df[!, "Index"])
    output_path = joinpath(save_path, "$(network_name)_features.csv")
    CSV.write(output_path, features_df)
end