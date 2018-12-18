struct BinnedData <: JAXTAMData
    mission::Symbol
    instrument::Symbol
    obsid::String
    bin_time::Real
    counts::Array{Int,1}
    times::StepRangeLen
    gtis::DataFrames.DataFrame
end

"""
    _lcurve_filter_time(event_times::Arrow.Primitive{Float64}, event_energies::Arrow.Primitive{Float64},

Largely not used, as the GTI filtering is enough to deal with out-of-time-range events, and 
manually filtering those events out early is computationally intensive

Function takes in event times and energies, then filters any events outside of the `start` and `stop` times

Optionally performs early filtering to remove low count (under 1/sec) GTIs, disabled by default as this is 
performed later anyway

Returns array of filtered times, enegies, and GTIs
"""
function _lcurve_filter_time(event_times::Arrow.Primitive{<:AbstractFloat}, event_energies::Arrow.Primitive{<:AbstractFloat},
    gtis::DataFrames.DataFrame, start_time::Number, stop_time::Number)
    @info "               -> Filtering times"

    if abs(event_times[1] - start_time) < 2.0^-8 && abs(event_times[end] - stop_time) < 2.0^-8
        # Times within tolerence
    else
        start_time_idx = findfirst(event_times .> start_time)
        stop_time_idx  = findfirst(event_times .>= stop_time)
        if stop_time_idx == nothing
            stop_time_idx = length(event_times)
        else
            stop_time_idx = stop_time_idx - 1
        end
        event_times    = event_times[start_time_idx:stop_time_idx]
        event_energies = event_energies[start_time_idx:stop_time_idx]
    end

    event_times = event_times .- start_time
    gtis = hcat(gtis[:START], gtis[:STOP])
    gtis = gtis .- start_time

    # WARNING: rounding GTIs up/down to get integers, shouldn't(?) cause a problem
    # TODO: check with Diego
    # gtis[:, 1] = ceil.(gtis[:, 1])
    # gtis[:, 2] = floor.(gtis[:, 2])

    if filter_low_count_gtis # Move GTI count filtering to post-binning stages
        counts_per_gti_sec = [count(gtis[g, 1] .<= event_times .<= gtis[g, 2])/(gtis[g, 2] - gtis[g, 1]) for g in 1:size(gtis, 1)]
        mask_min_counts = counts_per_gti_sec .>= 1

        @info "Excluded $(size(gtis, 1) - count(mask_min_counts)) gtis under 1 count/sec"

        gtis = gtis[mask_min_counts, :]
    end
    
    return Array(event_times), Array(event_energies), gtis
end

"""
    _lc_filter_energy(event_times::Array{T,1}, event_energies::Array{T,1}, good_energy_min::T, good_energy_max::T) where T <: AbstractFloat

Optionally filters events by a custom energy range, not just that given in the RMF files for a mission

Returns filtered event times and energies
"""
function _lc_filter_energy(event_times::Array{<:AbstractFloat,1}, event_energies::Array{<:AbstractFloat,1}, good_energy_min::AbstractFloat, good_energy_max::AbstractFloat)
    @info "               -> Filtering energies"
    mask_good_energy = good_energy_min .<= event_energies .<= good_energy_max

    @info "               -> excluded $(count(mask_good_energy.==false)) energies out of $good_energy_min -> $good_energy_max range"
    
    event_times = event_times[mask_good_energy]
    event_energies = event_energies[mask_good_energy]

    return event_times, event_energies
end

"""
    _lc_bin(event_times::Array{T,1}, bin_time::Union{T,W}, time_start::Union{T,W}, time_stop::Union{T,W}) where {T<:AbstractFloat, W<:Integer}

Bins the event times to bins of `bin_time` [sec] lengths

Performs binning out of memory for speed via `OnlineStats.jl` Hist function

Returns a range of `times`, with associated `counts` per time
"""
function _lc_bin(event_times::Array{<:AbstractFloat,1}, bin_time::Real, time_start::Real, time_stop::Real)
    @info "               -> Running OoM binning"

    if bin_time < 1
        if !ispow2(Int(1/bin_time))
            @warn "Bin time not pow2"
        end
    elseif bin_time > 1
        if !ispow2(Int(bin_time))
            @warn "Bin time not pow2"
        end
    elseif bin_time == 0
        throw(ArgumentError("Bin time cannot be zero"))
    end

    times = 0:bin_time:(time_stop-time_start)
    online_hist_fit = OnlineStats.Hist(times)

    histogram = fit!(online_hist_fit, event_times)
    counts    = value(histogram)[2]
    times     = times[1:end-1] # One less count than times

    return times, counts
end

"""
    _group_select(data::BinnedData)

Checks the differennce in time between GTIs, if the difference is under 
a `group_period` (128 [sec] by default) then the GTIs are in the same `group`

Done as data frequently has small breaks between GTIs, even though there is no 
significant gap in the lightcurve. Groups are used during plotting, periodograms, 
and when grouping/averaging together power spectra

Returns GTIs with an extra `:group` column added in
"""
function _group_select(data::BinnedData)
    group_period = 128

    gti_gap_time = data.gtis[:start][2:end] .- data.gtis[:stop][1:end-1]
    gti_gap_time = [0; gti_gap_time]

    gti_group_index = zeros(Int, size(data.gtis, 1))

    group = 1
    for (i, gap) in enumerate(gti_gap_time)
        if gap > group_period
            group = group + 1
        end

        gti_group_index[i] = group
    end

    data.gtis[:group] = gti_group_index

    return data
end

"""
    _group_return(data::BinnedData)

Uses the `_group_select` function to find the group each GTI belongs in

Using these groups, splits sections of light curve up into an group, then 
creates a new `BinnedData` for just the lightcurve of that one group, finally 
returns a `Dict{Int64,JAXTAM.BinnedData}` where each `Int64` is for a different group
"""
function _group_return(data::BinnedData)
    gtis = data.gtis
    # gtis = gtis[(gtis[:, :stop]-gtis[:, :start]) .>= 16, :] # Ignore gtis under 16s long
    available_groups = unique(gtis[:group])

    data_group = Dict{Int64,JAXTAM.BinnedData}()
    for group in available_groups
        if group == 0
            end_gti_idx = findfirst(data.gtis[:group] .== maximum(available_groups))

            first_gti_idx = findfirst(data.gtis[:group][end_gti_idx:end] .== group) + end_gti_idx
            last_gti_idx  = findlast(data.gtis[:group][end_gti_idx:end] .== group) + end_gti_idx

            if last_gti_idx > size(data.gtis, 1)
                last_gti_idx = size(data.gtis, 1)
            end
        else
            first_gti_idx = findfirst(data.gtis[:group] .== group)
            last_gti_idx  = findlast(data.gtis[:group] .== group)
        end

        first_gti = data.gtis[first_gti_idx, :]
        last_gti  = data.gtis[last_gti_idx,  :]

        first_idx = findfirst(data.times .>= first_gti[:start])
        last_idx  = findfirst(data.times .>= last_gti[:stop])

        if last_gti[1, :stop] >= data.times[end]
            last_idx = length(data.times)
        end
        
        data_group[group] = BinnedData(
            data.mission,
            data.instrument,
            data.obsid,
            data.bin_time,
            data.counts[first_idx:last_idx],
            data.times[first_idx:last_idx],
            data.gtis[first_gti_idx:last_gti_idx, :]
        )
    end

    return data_group
end

"""
    _lcurve(instrument_data::InstrumentData, bin_time::Union{Float64,Int64})

Takes in the `InstrumentData` and desired `bin_time`

Runs functions to perform extra time (`_lcurve_filter_time`) and energy (`_lc_filter_energy`) filtering

Runs the binning (`_lc_bin`) function, then finally `_group_select` to append group numbers to each GTI

Returns a `BinnedData` lightcurve
"""
function _lcurve(instrument_data::InstrumentData, bin_time::Union{Float64,Int64})
    event_times, event_energies, gtis = _lcurve_filter_time(instrument_data.events[:TIME], instrument_data.events[:E], instrument_data.gtis, instrument_data.start, instrument_data.stop)

    mission_name = instrument_data.mission
    good_energy_min, good_energy_max = (Float64(config(mission_name).good_energy_min), Float64(config(mission_name).good_energy_max))

    event_times, event_energies = _lc_filter_energy(event_times, event_energies, good_energy_min, good_energy_max)

    binned_times, binned_counts = _lc_bin(event_times, bin_time, instrument_data.start, instrument_data.stop)

    gtis = DataFrame(start=gtis[:, 1], stop=gtis[:, 2])

    lc = BinnedData(instrument_data.mission, instrument_data.instrument, instrument_data.obsid, bin_time, binned_counts, binned_times, gtis)

    lc = _group_select(lc)

    return lc
end

"""
    _lcurve_save(lightcurve_data::BinnedData, lc_dir::String)

Takes in `BinnedData` and splits the information up into three files, `meta`, `gtis`, and `data`

Saves the files in a lightcurve directory (`/JAXTAM/lc/\$bin_time/*`) per-instrument
"""
function _lcurve_save(mission_name::Symbol, obs_row::DataFrames.DataFrame, instrument::Symbol, 
        lightcurve_data::BinnedData, lc_dir::String; log=true)

    lc_basename = string("$(lightcurve_data.instrument)_lc_$(lightcurve_data.bin_time)")
    mkpath(lc_dir)

    lc_meta = DataFrame(mission=string(lightcurve_data.mission), instrument=string(lightcurve_data.instrument),
        obsid=lightcurve_data.obsid, bin_time=lightcurve_data.bin_time, times=[lightcurve_data.times[1], 
        lightcurve_data.times[2] - lightcurve_data.times[1], lightcurve_data.times[end]])

    lc_gtis = lightcurve_data.gtis # DataFrame(start=lightcurve_data.gtis[:, 1], stop=lightcurve_data.gtis[:, 2])

    lc_data = DataFrame(counts=lightcurve_data.counts)

    path_meta = joinpath(lc_dir, "$(lc_basename)_meta.feather")
    path_gtis = joinpath(lc_dir, "$(lc_basename)_gtis.feather")
    path_data = joinpath(lc_dir, "$(lc_basename)_data.feather")

    Feather.write(path_meta, lc_meta)
    Feather.write(path_gtis, lc_gtis)
    Feather.write(path_data, lc_data)

    if log
        _log_add(mission_name, obs_row, 
            Dict("data" =>
                Dict(:lcurve =>
                    Dict(instrument =>
                        Dict(
                            :path_meta => path_meta,
                            :path_gtis => path_gtis,
                            :path_data => path_data
                        )
                    )
                )
            )
        )
    end
end

"""
    _lc_read(lc_dir::String, instrument::Symbol, bin_time)

Reads the split files saved by `_lcurve_save`, combines them to return a single `BinnedData` type
"""
function _lc_read(lc_dir::String, instrument::Symbol, bin_time)
    lc_basename = string("$(instrument)_lc_$(string(bin_time))")

    lc_meta = Feather.read(joinpath(lc_dir, "$(lc_basename)_meta.feather"))
    lc_gtis = Feather.read(joinpath(lc_dir, "$(lc_basename)_gtis.feather"))
    lc_data = Feather.read(joinpath(lc_dir, "$(lc_basename)_data.feather"))

    return BinnedData(Symbol(lc_meta[:mission][1]), Symbol(lc_meta[:instrument][1]), lc_meta[:obsid][1], 
        lc_meta[:bin_time][1], lc_data[:counts], lc_meta[:times][1]:lc_meta[:times][2]:lc_meta[:times][3], 
        lc_gtis)
end

"""
    lcurve(mission_name::Symbol, obs_row::DataFrame, bin_time::Number; overwrite=false)

Main function, handles all the lightcurve binning

Runs binning functions if no files are found, then saves the generated `BinnedData`

Loads saved files if they exist

Returns `Dict{Symbol,BinnedData}`, with the instrument as a symbol, e.g. `lc[:XTI]` for NICER, 
`lc[:FPMA]`/`lc[:FPMB]` for NuSTAR
"""
function lcurve(mission_name::Symbol, obs_row::DataFrame, bin_time::Number; overwrite=false)
    obsid          = obs_row[:obsid][1]
    JAXTAM_path    = abspath(string(obs_row[:obs_path][1], "/JAXTAM/")); mkpath(JAXTAM_path)
    JAXTAM_content = readdir(JAXTAM_path)
    JAXTAM_lc_path = joinpath(JAXTAM_path, "lc/$bin_time/")
    
    mkpath(JAXTAM_lc_path)
    
    JAXTAM_c_files  = count(contains.(JAXTAM_content, "calib"))
    JAXTAM_lc_files = count(contains.(readdir(JAXTAM_lc_path), "lc_"))/3

    instruments = Symbol.(config(mission_name).instruments)

    lightcurves = Dict{Symbol,BinnedData}()

    if JAXTAM_c_files == 0
        @warn "No calibrated files found, running calibration"
        calibrate(mission_name, obs_row)
        return lcurve(mission_name, obs_row, bin_time; overwrite=overwrite)
    end

    if (JAXTAM_c_files > 0 && JAXTAM_c_files > JAXTAM_lc_files) || (JAXTAM_c_files > 0 && overwrite)
        calibrated_data = calibrate(mission_name, obs_row)
    
        for instrument in instruments
            @info "Binning LCURVE"
            lightcurve_data = _lcurve(calibrated_data[instrument], bin_time)

            @info "               -> Saving `$instrument` $(bin_time)s lightcurve data"

            _lcurve_save(lightcurve_data, JAXTAM_lc_path)

            lightcurves[instrument] = lightcurve_data
        end
    elseif JAXTAM_c_files > 0 && JAXTAM_c_files == JAXTAM_lc_files
        @info "Loading LC $(obsid): from $JAXTAM_lc_path"

        for instrument in instruments
            lightcurves[instrument] = _lc_read(JAXTAM_lc_path, instrument, bin_time)
        end
    end
    
    return lightcurves
end

"""
    lcurve(mission_name::Symbol, obsid::String, bin_time::Number; overwrite=false)

Runs `master_query` to find the desired `obs_row` for the observation

Calls main `lcurve(mission_name::Symbol, obs_row::DataFrame, bin_time::Number; overwrite=false)` function
"""
function lcurve(mission_name::Symbol, obsid::String, bin_time::Number; overwrite=false)
    obs_row = JAXTAM.master_query(mission_name, :obsid, obsid)

    return lcurve(mission_name, obs_row, bin_time; overwrite=overwrite)
end