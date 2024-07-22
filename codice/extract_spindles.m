
function ss = extract_spindles(eeg_data, instants, samples, sampling_rate)

    ss = zeros(samples,  size(eeg_data, 2) , length(instants));  

    for i = 1 : length(instants)
        start_spindles = (instants(i) * sampling_rate);
        end_spindles = start_spindles + samples - 1;
        ss(:, :, i) = eeg_data(start_spindles : end_spindles, :);
    end
end

