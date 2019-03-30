struct FFTData <: JAXTAMData
    mission        :: Mission
    instrument     :: Symbol
    obsid          :: String
    e_range        :: Tuple{Float64,Float64}
    bin_time       :: Real
    bin_size       :: Int
    bin_count      :: Int
    gti_index      :: Int
    gti_start_time :: Real
    group          :: Int
    src_ctrate     :: Array
    bkg_ctrate     :: Array
    power          :: Array
    avg_power      :: Array
    freq           :: Array
end

function _FFTData(gti, freq, power, fspec_bin_size, bin_count, src_ctrate)

    return FFTData(gti.mission, gti.instrument, gti.obsid, gti.e_range, gti.bin_time, fspec_bin_size, bin_count,
        gti.gti_index, gti.gti_start_time, gti.group, src_ctrate, zeros(size(src_ctrate)), power, mean(power, dims=2)[:], freq)
end

function _fft(counts::Array, times::StepRangeLen, bin_time::Real, fspec_bin_size; leahy=true)
    FFTW.set_num_threads(4)

    spec_no = Int(floor(length(counts)/fspec_bin_size))

    counts = counts[1:Int(spec_no*fspec_bin_size)]
    counts = reshape(counts, fspec_bin_size, spec_no)

    src_ctrate = mean(counts, dims=1)./bin_time

    power = abs.(rfft(counts))
    freq = Array(rfftfreq(fspec_bin_size, 1/bin_time))

    if leahy
        power = (2 .*(power.^2)) ./ sum(counts)
    end

    return freq, power, spec_no, src_ctrate
end

function _best_gti_pow2(gtis::Dict{Int64,JAXTAM.GTIData})
    gti_pow2s = [prevpow(2, length(gti[2].counts)) for gti in gtis]
    min_pow2  = Int(median(gti_pow2s))
    exclude   = gti_pow2s .< min_pow2
    
    instrument = gtis[collect(keys(gtis))[1]].instrument

    @info "`$instrument` median `pow2` length is $min_pow2, excluded $(count(exclude)) gtis,
        $(length(gti_pow2s) - count(exclude)) remain"

    good_gtis = [gti for gti in gtis][.!exclude]

    return Dict(good_gtis), min_pow2
end

function gtis_pow2(mission_name::Symbol, gtis::Dict{Symbol,Dict{Int64,JAXTAM.GTIData}})
    instruments = config(mission_name).instruments

    gtis_pow2 = Dict{Symbol,Dict{Int64,JAXTAM.GTIData}}()

    for instrument in instruments
        gtis_pow2[Symbol(instrument)] = _best_gti_pow2(gtis[Symbol(instrument)])
    end

    return gtis_pow2
end

function _fspec_save(fspec_data::Dict{Int64,JAXTAM.FFTData}, fspec_dir::String)
    mkpath(fspec_dir)
    log_entries = Dict{Union{Symbol,Int},String}()

    fspec_indecies = collect(keys(fspec_data))
    fspec_starts   = [t.gti_start_time for t in values(fspec_data)]
    fspec_b_counts = [t.bin_count for t in values(fspec_data)]
    fspec_groups   = [o.group for o in values(fspec_data)]
    fspec_example  = fspec_data[fspec_indecies[1]]

    fspec_basename = string("$(fspec_example.instrument)_lc_$(fspec_example.bin_time)_fspec")
    fspec_meta     = DataFrame(mission=_mission_name(fspec_example.mission),
    instrument=String(fspec_example.instrument), obsid=fspec_example.obsid,
    bin_time=fspec_example.bin_time, bin_size=fspec_example.bin_size, groups=fspec_groups,
    indecies=fspec_indecies, starts=fspec_starts, b_counts=fspec_b_counts)
    
    Feather.write(joinpath(fspec_dir, "$(fspec_basename)_meta.feather"), fspec_meta)
    
    for index in fspec_indecies
        fspec_savepath  = joinpath(fspec_dir, "$(fspec_basename)_$(index).feather")

        fspec_data_df = DataFrame(fspec_data[index].power)
        fspec_data_df[:freq]   = fspec_data[index].freq
        fspec_data_df[:avg_power] = fspec_data[index].avg_power
        
        Feather.write(fspec_savepath, fspec_data_df)
    end
    @warn "Logging not implemented yet as saving is off by default" "Contact dev to request fix"
    # TODO: Fix this... eventually
end

function _fspec_load(fspec_dir::String, instrument::Symbol, e_range::Tuple{Float64,Float64}, bin_time::Number, fspec_bin::Number)
    bin_time = float(bin_time)

    fspec_basename  = string("$(instrument)_lc_$(bin_time)_fspec")
    fspec_meta_path = joinpath(fspec_dir, "$(fspec_basename)_meta.feather")

    fspec_meta = Feather.read(fspec_meta_path)

    fspec_data = Dict{Int64,JAXTAM.FFTData}()

    for fspec_row in eachrow(fspec_meta)
        fspec_mission = _mission_symbol_to_type(Symbol(fspec_row[:mission]))
        fspec_inst    = Symbol(fspec_row[1, :instrument])
        fspec_obsid   = fspec_row[1, :obsid]
        fspec_bin_t   = fspec_row[1, :bin_time]
        fspec_bin_sze = fspec_row[1, :bin_size]
        fspec_idx     = fspec_row[1, :indecies]
        fspec_starts  = fspec_row[1, :starts]
        fspec_b_count = fspec_row[1, :b_counts]
        fspec_group   = fspec_row[1, :groups]
        fspec_file    = Feather.read(joinpath(fspec_dir, "$(fspec_basename)_$(fspec_idx).feather"))
        fspec_freq   = Array(fspec_file[:freq]); delete!(fspec_file, :freq)
        fspec_ampavg  = Array(fspec_file[:avg_power]); delete!(fspec_file, :avg_power)
        if fspec_idx != -1
            fspec_power = convert(Array, fspec_file[:, names(fspec_file)])
        else # If avg_power is loaded, fspec_power = [], throws Argument Error on convert
            fspec_power = []
        end

        fspec_data[fspec_idx] = FFTData(fspec_mission, fspec_inst, fspec_obsid, e_range, fspec_bin_t,
            fspec_bin_sze, fspec_b_count, fspec_idx, fspec_starts, fspec_group, fspec_power, fspec_ampavg,
            fspec_freq)
    end

    return fspec_data
end

function _fspec(gtis::Dict{Int64,JAXTAM.GTIData}, fspec_bin::Real; pow2=true, fspec_bin_type=:time)
    powspec = Dict{Int64,FFTData}()

    first_gti = gtis[collect(keys(gtis))[1]]

    if fspec_bin_type == :time
        fspec_bin_size = Int(fspec_bin/first_gti.bin_time)
    elseif fspec_bin_type == :length
        fspec_bin_size = Int(fspec_bin)
    end

    @info "                       -> fspec_bin_size = $fspec_bin_size"

    if !ispow2(fspec_bin_size)
        @warn "`fspec_bin_size` not `pow2`: $fspec_bin_size, `fspec_bin` = $fspec_bin ($(1/fspec_bin))"
    end

    for gti in values(gtis)
        if length(gti.counts) < fspec_bin_size
            @info "                       -> skipped gti $(lpad(gti.gti_index, 3)) (under by: $(length(gti.counts) - fspec_bin_size))"
            continue
        end

        freq, power, bin_count, src_ctrate = _fft(gti.counts, gti.times, gti.bin_time, fspec_bin_size)

        powspec[gti.gti_index] = _FFTData(gti, freq, power, fspec_bin_size, bin_count, src_ctrate)
    end

    return powspec
end

function _scrunch_sections(instrument_data::Dict{Int64,JAXTAM.FFTData}; append_mean=true)
    gtis = keys(instrument_data)
    gti_example = instrument_data[collect(gtis)[1]]

    joined_powspecs = Dict{Symbol,Array}()
    joined_together = Array{Float64,2}(undef, size(gti_example.power, 1), 0)

    bin_count_sum = 0

    mean_src_ctrate = Array{Float64,1}(undef, 1)
    mean_bkc_ctrate = Array{Float64,1}(undef, 1)
    for gti in values(instrument_data)
        if gti.gti_index <= 0 # Don't try and scrunch special, non-gti values
            continue
        end

        bin_count_sum += gti.bin_count
        joined_together = hcat(joined_together, gti.power)
        mean_src_ctrate = hcat(mean_src_ctrate, gti.src_ctrate)
        mean_bkc_ctrate = hcat(mean_bkc_ctrate, gti.bkg_ctrate)
    end

    mean_src_ctrate = [mean(mean_src_ctrate[2:end])]
    mean_bkc_ctrate = [mean(mean_bkc_ctrate[2:end])]

    if append_mean
        instrument_data[-1] = FFTData(
            gti_example.mission, gti_example.instrument, gti_example.obsid, gti_example.e_range,
            gti_example.bin_time, gti_example.bin_size, bin_count_sum, -1, -1, -1,
            mean_src_ctrate, mean_bkc_ctrate,
            [], mean(joined_together, dims=2)[:], gti_example.freq
        )
    end

    return instrument_data
end

function fspec(mission_name::Symbol, gtis::Dict{Symbol,Dict{Int64,JAXTAM.GTIData}}, fspec_bin::Real; pow2=true, fspec_bin_type=:time, scrunched=true)
    instruments = config(mission_name).instruments

    instrument_fspecs = Dict{Symbol,Dict{Int,JAXTAM.FFTData}}()

    for instrument in instruments
        instrument_fspecs[Symbol(instrument)] = _fspec(gtis[Symbol(instrument)], fspec_bin, pow2=pow2, fspec_bin_type=fspec_bin_type)

        if scrunched
            instrument_fspecs = _scrunch_sections(instrument_fspecs[Symbol(instrument)])
        end
    end

    return instrument_fspecs
end

function fspec(mission::Mission, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame},
        bin_time::Real, fspec_bin::Real; overwrite=true, overwrite_gtis=false, save_fspec=false,
        pow2=true, fspec_bin_type=:time, scrunched=true, e_range::Tuple{Float64,Float64}=_mission_good_e_range(mission),
        gtis_data::Dict{Symbol,Dict{Int64,JAXTAM.GTIData}}=Dict{Symbol,Dict{Int64,JAXTAM.GTIData}}()
    )
    obsid       = obs_row[:obsid]
    instruments = _mission_instruments(mission)
    fspec_files = _log_query(mission, obs_row, "data", :fspec, e_range)

    fspec_path  = abspath(
        _obs_path_local(mission, obs_row; kind=:jaxtam),
        "JAXTAM",
        _log_query_path(; category=:data, kind=:fspec, e_range=e_range, bin_time=bin_time)
    )

    fspecs = Dict{Symbol,Dict{Int64,JAXTAM.FFTData}}()
    for instrument in instruments
        if ismissing(fspec_files) || !haskey(fspec_files, instrument) || overwrite
            @info "Missing fspec files for $instrument"

            if !haskey(gtis_data, instrument)
                gtis_data = JAXTAM.gtis(mission, obs_row, bin_time, overwrite=overwrite_gtis, e_range=e_range)
            end

            @assert _recursive_first(gtis_data).e_range == e_range "gti energy range does not match input: $(_recursive_first(gtis_data).e_range) != $e_range"

            @info "Generating fspec for $instrument"
            fspec_data = _fspec(gtis_data[instrument], fspec_bin; pow2=pow2, fspec_bin_type=fspec_bin_type)

            if length(fspec_data) == 0
                throw(JAXTAMError("No GITs left in fspec", :fspec))
            end

            if scrunched
                @info "                       -> scrunching gtis"
                fspec_data = _scrunch_sections(fspec_data)
            end

            if save_fspec
                @info "                       -> saving $instrument -> $fspec_path"
                _fspec_save(fspec_data, fspec_path)
            else
                @info "                       -> not saving fspec file, `save_fspec` is false"
            end

            fspecs[instrument] = fspec_data
        else
            @info "Loading $instrument fspecs"
            fspecs[instrument] = _fspec_load(fspec_path, instrument, bin_time, fspec_bin)
        end
    end

    return fspecs
end

function fspec(mission::Mission, obsid::String, bin_time::Number, fspec_bin::Real;
    overwrite_fs=false, overwrite_gtis=false, pow2=true, fspec_bin_type=:time, scrunched=true, e_range=_mission_good_e_range(mission))

    obs_row = master_query(mission, :obsid, obsid)

    fs = fspec(mission, obs_row, bin_time, fspec_bin; overwrite=overwrite_fs, overwrite_gtis=overwrite_gtis,
        pow2=pow2, fspec_bin_type=fspec_bin_type, scrunched=scrunched, e_range=e_range)
    
    return fs
end

function _fspec_rebin(power, freq, bin_count::Int, bin_size, bin_time, rebin=(:log10, 0.01))
    rebin_type     = rebin[1]
    rebin_factor   = rebin[2]
    
    freq_intervals = missing
    if ismissing(bin_time)
        freq_intervals = mean(diff(freq))
    else
        freq_intervals = 1/(bin_size*bin_time) #freq[3] - freq[2]
    end
    
    # Frequency rebinning (e.g. binning to 1Hz, 2Hz, etc... intervals) is 
    # just linear rebinning, with the factor determined by the desired frequency
    binary_thresh  = missing
    if rebin_type == :freq
        rebin_type   = :linear
        rebin_factor = round(Int, bin_size*bin_time*rebin_factor)
    elseif rebin_type == :freq_binary
        rebin_factor  = round(Int, bin_size*bin_time*rebin_factor)
        binary_thresh = rebin[3]
    end

    if rebin_type == :log10
        freq_max = freq[end]

        scale_start  = log10(rebin_factor)

        scale = exp10.(range(scale_start, step=rebin_factor, stop=log10(freq_max)))
        scale = [0; scale[1:end-1]]
        scale = scale./freq_intervals
        scale = ceil.(Int, scale[:, 1])
        scale[1, 1] = 1
        scale = unique(scale)

        final_scale_value = ceil(Int, exp10(log10(scale[end])+rebin_factor))
        final_scale_value > length(power) ? final_scale_value=length(power) : ""
        scale = [scale [scale[2:end]; final_scale_value].+1]

        errors = power
        errors = [mean(errors[scale[i, 1]:scale[i, 2]]) for i = 1:size(scale, 1)]
        errors = errors ./ sqrt.(diff(scale, dims=2).*bin_count)

        rebinned_fspec = [mean(power[scale[i, 1]:scale[i, 2]]) for i = 1:size(scale, 1)]
        
        freq_scale = freq[round.(Int, (scale[:, 2]+scale[:, 1])./2)]
        
        rebinned_fspec[rebinned_fspec.==0] .= NaN
    elseif rebin_type == :linear
        if rebin_factor < 1 || typeof(rebin_factor) != Int
            error("Rebin factor must be a >1 integer for linear rebinning")
        end

        in_rebin = floor(Int, size(power, 1)/rebin_factor).*rebin_factor

        errors = power ./ sqrt(bin_count*rebin_factor)
        errors = mean(reshape(errors[1:in_rebin], rebin_factor, :), dims=1)'

        rebinned_fspec = mean(reshape(power[1:in_rebin], rebin_factor, :), dims=1)'
        rebinned_fspec = rebinned_fspec[:]

        freq_scale = mean(reshape(freq[1:in_rebin], rebin_factor, :), dims=1)'
        freq_scale = freq_scale[:]

        rebinned_fspec[freq_scale .== 0] .= NaN
    elseif rebin_type == :freq_binary
        if rebin_factor < 1 || typeof(rebin_factor) != Int
            error("Rebin factor must be a >1 integer for linear rebinning")
        end

        in_rebin = floor(Int, size(power, 1)/rebin_factor).*rebin_factor

        errors = power ./ sqrt(bin_count*rebin_factor)
        errors = zeros(in_rebin) # Errors don't really have a meaning for this

        binned_power    = reshape(power[1:in_rebin], rebin_factor, :)
        rebinned_fspec = zeros(1, size(binned_power, 2))
        rebinned_fspec[any(binned_power .>= binary_thresh, dims=1)] .= binary_thresh
        rebinned_fspec = rebinned_fspec'

        freq_scale = mean(reshape(freq[1:in_rebin], rebin_factor, :), dims=1)'
    else
        error("Invalid rebin type $rebin_type, must be :log10, :linear, :freq, or :freq_binary")
    end

    return freq_scale[:], rebinned_fspec, errors
end

function fspec_rebin(fs::JAXTAM.FFTData; rebin=(:log10, 0.01))
    rebinned_data = _fspec_rebin(fs.avg_power, fs.freq, fs.bin_count, fs.bin_size, fs.bin_time, rebin)

    return rebinned_data
end

function fspec_rebin_sgram(fs::Dict{Int64,JAXTAM.FFTData}; rebin=(:log10, 0.01))
    fs = sort(fs) # Ensure po=wer spectra is in GTI order
    fs_groups       = unique([f[2].group for f in fs if f[1]>0])
    fs_group_bounds = cumsum([sum([f[2].bin_count for f in fs if f[2].group==g]) for g in fs_groups])
    #fs_group_bounds = cumsum([f[2].bin_count for f in fs if f[1]>0])

    fs_freq   = fs[-1].freq
    fs_power  = hcat([f[2].power for f in fs if f[1] > 0]...) # Skip special < 0 keys
    fs_power[1, :] .= NaN # NaN 0 Hz

    bin_size = fs[-1].bin_size; bin_time = fs[-1].bin_time
    fs_rebinned_power  = hcat([_fspec_rebin(fs_power[:, i], fs_freq, 1, bin_size, bin_time, rebin)[2] for i in 1:size(fs_power, 2)]...)
    fs_rebinned_freq = fspec_rebin(fs[-1], rebin=rebin)[1]

    return fs_rebinned_freq, fs_rebinned_power, fs_group_bounds, fs_groups
end

function fspec_pulses(fs::Dict{Int64,JAXTAM.FFTData};
        freq_bin=1, power_limits=[10, 20, 100])

    rebinned_data = Dict{Int,Tuple}()
    for p in power_limits
        rebinned_data[p] = fspec_rebin_sgram(fs;
            rebin=(:freq_binary, freq_bin, p))
    end

    return rebinned_data
end