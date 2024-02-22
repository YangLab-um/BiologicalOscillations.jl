using BiologicalOscillations, Statistics, StatsBase, DataFrames, JLD2, CSV

"""
This script provides the following summary data for each network in a single-addition simulation:
    - Network: Informal name of the network
    - Connectivity: Connectivity matrix of the network
    - Addition coherence: Coherence created by the addition of a single edge (coherent or incoherent)
    - Addition type: Type of feedback loop created by the addition of a single edge (positive or negative)
    - Addition length: Length of the feedback loop created by the addition of a single edge. A self-loop has length 1.
    - Oscillation probability: Ratio of oscillatory solutions to total solutions
    - Mean logFreq: Mean log(Frequency) value across oscillatory solutions
    - Mean amplitude: Mean amplitude value across oscillatory solutions
    - Stdev logFreq: Standard deviation of log(Frequency) values across oscillatory solutions
    - Stdev amplitude: Standard deviation of amplitude values across oscillatory solutions
    - Batches: Number of batches in which the network was simulated
    - Total samples: Total number of samples in which the network was simulated
    - Oscillatory solutions: Total number of oscillatory solutions
    - IQR logFreq: Interquartile range of log(Frequency) values across oscillatory solutions
    - IQR amplitude: Interquartile range of amplitude values across oscillatory solutions

The data is stored in a .csv file and contains the information for all networks within the folder provided.

Inputs:
    - Path to folder containing the results of the single-addition simulations
    - Path to file with network names and their connectivity matrices
    - Path where the summary file will be saved
    - Name for the summary file
    - Connectivity matrix of the reference network
    - Whether the data in the folders is the full simulation data or the reduced version with only `simulation_results`

Note: A single-addition simulation is a simulation where the network is simulated with a single-edge-addition 
      with respect to the negative-feedback-only network.

Author: Franco Tavella
Date: 09/11/2023
"""

# User-defined parameters
network = "P0"
simulation_data_location = "E:/Project 2 - Tunability/ForkBiologicalOscillations.jl/results_analyses/$(network)_1_add/find_oscillations_result/$(network)/simulations"
# simulation_data_location = "H:/Franco/Project 2 - Tunability - Data/S3_simulation_result/simulations"
connectivity_file_location = "E:/Project 2 - Tunability/ForkBiologicalOscillations.jl/results_analyses/$(network)_1_add/$(network)_1_add_connectivity_list.csv"
save_path = "E:/Project 2 - Tunability/ForkBiologicalOscillations.jl/results_analyses/$(network)_1_add/"
save_name = "$(network)_one_add_simulation_summary.csv"
reference_connectivity = [0 0 0 0 -1; -1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0]
is_full_simulation_data = true

# Get list of network names and their repetitions
network_names = readdir(simulation_data_location)
network_repetitions = [readdir(joinpath(simulation_data_location, network)) for network in network_names]

# Initialize summary dataframe
summary_column_names = [
    "Network", "Connectivity", "Addition coherence", "Addition type", "Addition length",
    "Oscillation probability", "Mean logFreq", "Mean amplitude", "Stdev logFreq", "Stdev amplitude",
    "Batches", "Total samples", "Oscillatory solutions", "IQR logFreq", "IQR amplitude"
    ]
summary_column_types = [
    String[], String[], String[], String[], Int64[],
    Float64[], Float64[], Float64[], Float64[], Float64[],
    Int64[], Int64[], Int64[], Float64[], Float64[]
    ]
summary_df = DataFrame(summary_column_types, summary_column_names)

# Load connectivity data
connectivity_df = CSV.read(connectivity_file_location, DataFrame, header=false)
rename!(connectivity_df, Symbol("Column1") => "Network")
rename!(connectivity_df, Symbol("Column2") => "Connectivity")

for (idx, network_name) in enumerate(network_names)
    print("Processing $(network_name)...\n")
    # Calculate addition coherence, type and length
    connectivity_str = connectivity_df[connectivity_df.Network .== network_name, "Connectivity"][1]
    connectivity_array = connectivity_string_to_matrix(connectivity_str)
    loop_properties = classify_single_addition(reference_connectivity, connectivity_array)
    # Merge all batches into a single dataframe
    batch_size = size(network_repetitions[idx], 1)
    total_samples = 0
    oscillatory_samples = 0
    features_df = DataFrame([Int64[], Int64[], Float64[],Float64[]], 
                            ["Batch", "Index", "Frequency", "Amplitude"])
    for batch in 1:batch_size
        if is_full_simulation_data
            fpath = "$(simulation_data_location)/$(network_name)/$(batch)/$(network_name)_$(batch).jld2"
            data = load(fpath)
            samples = size(data["pin_result"]["parameter_sets"], 1)
            simulation_result = data["pin_result"]["simulation_result"]
        else
            fpath = "$(simulation_data_location)/$(network_name)/$(batch)/$(network_name)_$(batch)_simulation_result.jld2"
            data = load(fpath)
            samples = data["number_of_samples"]
            simulation_result = data["simulation_result"]
        end


        oscillatory_df = filter(row -> row["is_oscillatory"] == true, simulation_result)
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
    # Process merged features
    mean_logfreq = mean(log10.(features_df[!, "Frequency"]))
    mean_amplitude = mean(features_df[!, "Amplitude"])
    stdev_logfreq = std(log10.(features_df[!, "Frequency"]))
    stdev_amplitude = std(features_df[!, "Amplitude"])
    iqr_logfreq = iqr(log10.(features_df[!, "Frequency"]))
    iqr_amplitude = iqr(features_df[!, "Amplitude"])
    # Oscillation probability
    oscillation_probability = oscillatory_samples / total_samples
    # Append data to summary dataframe
    network_summary_data = [
        [network_name], [connectivity_str], [loop_properties.loop_coherence[1]], 
        [loop_properties.loop_type[1]], [loop_properties.loop_length[1]], 
        [oscillation_probability], [mean_logfreq], [mean_amplitude], [stdev_logfreq],
        [stdev_amplitude], [batch_size], [total_samples], [oscillatory_samples], 
        [iqr_logfreq], [iqr_amplitude]
        ]
    partial_summary_df = DataFrame(network_summary_data, summary_column_names)
    append!(summary_df, partial_summary_df)
end

# Save summary dataframe as a csv file
summary_df[!, "Batches"] = convert(Vector{Int64}, summary_df[!, "Batches"])
summary_df[!, "Total samples"] = convert(Vector{Int64}, summary_df[!, "Total samples"])
summary_df[!, "Oscillatory solutions"] = convert(Vector{Int64}, summary_df[!, "Oscillatory solutions"])
output_path = joinpath(save_path, save_name)
CSV.write(output_path, summary_df)