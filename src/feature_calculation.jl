"""
        calculate_main_frequency(ode_solution::ODESolution, sampling::Int, fft_points::Int)

Obtains the dominant frequencies and their powers for each signal on `ode_solution` via the fast fourier transform

# Arguments
- `ode_solution::ODESolution`: The full output of the `solve()` function from `DifferentialEquations.jl`
- `sampling::Int`: Number of points to be used for uniform re-sampling of `ode_solution`. The same sampling is used for all signals
- `fft_points::Int`: Number of points in the frequency spectrum

# Returns
- `frequency_data::Dict`: Dictionary containing the dominant frequencies and their powers for each signal. Here power means the amplitude of that frequency on the FFT spectrum. The output is encoded as:
```julia
frequency_data = Dict('frequency' => [freq_signal_1, ..., freq_signal_N], 
                      'power' => [power_signal_1, ..., power_signal_N])
```

# Note
- The provided `ode_solution` has to be solved with the `dense = true` setting to be suitable for this function since uniform re-sampling is performed.
"""
function calculate_main_frequency(ode_solution::ODESolution, sampling::Int, fft_points::Int)
    number_of_signals = length(ode_solution.u[1])
    main_frequency = fill(NaN, number_of_signals)
    main_power = fill(NaN, number_of_signals)
    signal_sampling = LinRange(ode_solution.t[1], ode_solution.t[end], sampling)
    sampling_frequency = sampling / (ode_solution.t[end] - ode_solution.t[1])
    for i in 1:number_of_signals
        # Remove mean -- Avoids a dominant low frequency caused by mean ≠ 0
        signal = ode_solution(signal_sampling)[i,:] .- mean(ode_solution(signal_sampling)[i,:])
        spectrum = periodogram(signal; fs=sampling_frequency, nfft=fft_points)
        frequency_peaks, ~ = findmaxima(spectrum.power)
        # Remove small peaks -- Discard any peak with height less than 1/3 of the largest spectral value
        frequency_peaks, ~ = peakproms(frequency_peaks, spectrum.power, minprom=maximum(spectrum.power)/3.0)
        if length(frequency_peaks) >= 1
            # From the top peaks save the one with the lowest frequency
            freq_power_sorted = sort(collect(zip(spectrum.freq[frequency_peaks], spectrum.power[frequency_peaks])); by=first)
            main_frequency[i] = freq_power_sorted[1][1]
            main_power[i] = freq_power_sorted[1][2]
        end
    end

    return Dict("frequency" => main_frequency, "power" => main_power)
end


"""
    calculate_amplitude(ode_solution::ODESolution)

Calculate the amplitude of each signal in `ode_solution` and the relative peak/trough variation.

# Arguments
- `ode_solution::ODESolution`: The full output of the `solve()` function from `DifferentialEquations.jl`

# Returns
- `amplitude_data::Dict`: Dictionary containing the amplitude for each signal and the peak/trough variation. The output is encoded as: 
```julia
amplitude_data = Dict('amplitude' => [amp_signal_1, ..., amp_signal_N], 
                      'peak_variation' => [peak_var_signal_1, ..., peak_var_signal_N], 
                      'trough_variation' => [trough_var_signal_1, ..., trough_var_signal_N])
```

# Notes
- Variation of peaks and troughs is calculated as the standard deviation of all the peaks/troughs in the signal divided by the amplitude. This quantity allows to distinguish stable oscillatory solution from dampened ones.
- To find the most relevant peaks/troughs, we discard any peak/trough that has a prominence less than 1/3 of the total amplitude of the signal (calculated as the maximum - minimum of the signal). Therefore, it is better to use this function with a signal that has little to none transient behavior.
"""
function calculate_amplitude(ode_solution::ODESolution)
    number_of_signals = length(ode_solution.u[1])
    amplitude = fill(NaN, number_of_signals)
    peak_variation = fill(NaN, number_of_signals)
    trough_variation = fill(NaN, number_of_signals)

    for i in 1:number_of_signals
        rough_amp = maximum(ode_solution[i,:]) - minimum(ode_solution[i,:])
        # Peaks
        peak_positions, ~ = findmaxima(ode_solution[i,:])
        peak_positions, ~ = peakproms(peak_positions, ode_solution[i,:], minprom=rough_amp/3.0)
        peak_values = ode_solution[i, peak_positions]
        # Troughs
        trough_positions, ~ = findminima(ode_solution[i,:])
        trough_positions, ~ = peakproms(trough_positions, ode_solution[i,:], minprom=rough_amp/3.0)
        trough_values = ode_solution[i, trough_positions]
        if length(peak_values) > 0 && length(trough_values) > 0
            amplitude[i] = mean(peak_values) - mean(trough_values)
            peak_variation[i] = std(peak_values)/amplitude[i]
            trough_variation[i] = std(trough_values)/amplitude[i]
        end
    end
    return Dict("amplitude" => amplitude, "peak_variation" => peak_variation, "trough_variation" => trough_variation)
end


"""
    is_ODE_oscillatory(frequency_data::Dict, amplitude_data::Dict, freq_variation_threshold=0.05, power_threshold=1e-7, amp_variation_threshold=0.01)

Returns true if an `ODESolution` is oscillatory based on the calculated frequency and amplitude. False otherwise.

# Arguments (Required)
- `frequency_data::Dict`: The output from [`calculate_main_frequency`](@ref)
- `amplitude_data::Dict`: The output from [`calculate_amplitude`](@ref)

# Arguments (Optional)
- `freq_variation_threhsold::Real`: Maximum tolerated variation in frequency between species in the solution to be declared oscillatory. Default value 0.05.
- `power_threshold::Real`: Minimum spectral power that the main peak has to have to be declared oscillatory. Default value 1e-7.
- `amp_variation_threhsold::Real`: Maximum tolerated value for peak/trough variation to be declared oscillatory. Default value 0.01.
"""
function is_ODE_oscillatory(frequency_data::Dict, amplitude_data::Dict; freq_variation_threshold=0.05, power_threshold=1e-7, amp_variation_threshold=0.01)
    # If any quantity is NaN, the system isn't considered oscillatory
    freq_is_nan = any(isnan.(frequency_data["frequency"]))
    amp_is_nan = any(isnan.(amplitude_data["amplitude"]))
    if freq_is_nan || amp_is_nan
        return false
    end
    # Frequency stability
    mean_freq = mean(frequency_data["frequency"])
    freq_variation = mean(frequency_data["frequency"] .- mean_freq) / mean_freq
    freq_is_stable = freq_variation < freq_variation_threshold
    # Power amplitude
    power_is_significant = all(frequency_data["power"] .> power_threshold)
    # Peak variation
    peak_is_stable = all(amplitude_data["peak_variation"] .< amp_variation_threshold)
    # Trough variation
    trough_is_stable = all(amplitude_data["trough_variation"] .< amp_variation_threshold)

    if freq_is_stable && power_is_significant && peak_is_stable && trough_is_stable
        return true
    else
        return false
    end
end