struct GTIData <: JAXTAMData
    mission        :: Mission
    instrument     :: Symbol
    obsid          :: String
    e_range        :: Tuple{Float64,Float64}
    bin_time       :: Real
    gti_index      :: Int
    gti_start_time :: Real
    counts         :: Array
    times          :: StepRangeLen
    group          :: Int
end

"""
    _lc_filter_gtis(binned_times::StepRangeLen, binned_counts::Array{Int,1}, gtis::DataFrames.DataFrame, mission::Symbol, instrument::Symbol, obsid::String; min_gti_sec=16)

Splits the lightcurve (count) data into GTIs

First, removes GTIs under `min_gti_sec`, then puts the lightcurve data 
into a `Dict{Int,GTIData}`, with the key as the `index` of the GTI
"""
function _lc_filter_gtis(mission::Mission, instrument::Symbol, obsid::String, e_range::Tuple{Float64,Float64},
        binned_times::StepRangeLen, binned_counts::Array{Int,1}, gtis::DataFrames.DataFrame;
        min_gti_sec=16
    )
    gti_data = Dict{Int,GTIData}()

    gti_groups = gtis[:group]
    gti_times  = hcat(gtis[:start], gtis[:stop])
    
    bin_time = binned_times[2] - binned_times[1]
    
    gti_mask = (gtis[:stop] .- gtis[:start]) .>= min_gti_sec # Exclude GTIs under `min_gti_sec`
    gti_times[gti_mask.==false, :] .= -1 # Set excluded GTIs to -1, GTIs aren't just discarded so that their absolute index can be kept track of

    excluded_gti_count = count(gti_mask.==false)

    @info "               -> prelim. excluded $excluded_gti_count GTIs under $(min_gti_sec)s"

    if excluded_gti_count == length(gti_groups)
        # All GTIs have been excluded
        throw(JAXTAMError("All GTIs under $(min_gti_sec)s, none left", :gtis, BoundsError(Any[], (1,))))
    end

    # Do this as the first GTI being zero screws with the initial start index
    # in the upcoming `for gtis` loop
    if gti_times[1, 1] == 0
        gti_times[1, 1] = eps()
    end

    # Dodgy way to convert a matrix into an array of arrays
    # so each GTI is stored as an array of [start; finish]
    # and each of those GTI arrays is an array itself
    # makes life a bit easier for the following `for gti in gtis` loop
    gti_times = [gti_times[x, :] for x in 1:size(gti_times, 1)]

    @info "               -> sorting GTIs"

    # Offset if binned_times are a subarray of the lightcurve
    if typeof(binned_times) <: SubArray
        subarray_offset = binned_times.offset1
    else
        subarray_offset = 0
    end

    for (i, gti) in enumerate(gti_times) # For each GTI, store the selected times and count rate within that GTI
        if gti[1] == -1 # Short GTIs have their start time set to -1, skipped here
            continue
        end
        
        start = ceil(Int, gti[1]/bin_time)+1  - subarray_offset
        stop  = floor(Int, gti[2]/bin_time)-1 - subarray_offset

        # Subtract GTI start time from all times, so all start from t=0
        array_times = binned_times[start:stop].-gti[1]
        range_times = array_times[1]:bin_time:array_times[end]

        gti_data[Int(i)] = GTIData(mission, instrument, obsid, e_range, bin_time, i, start, 
            binned_counts[start:stop], range_times, gti_groups[i])
    end

    total_counts = sum(binned_counts)[1]
    gti_counts   = sum([sum(gti.counts) for gti in values(gti_data)])
    count_delta  = gti_counts-total_counts
    delta_prcnt  = round(count_delta/total_counts*100, digits=2)

    @info "               -> original counts: $total_counts | remaining: $gti_counts | delta: $count_delta ($delta_prcnt %)"

    if delta_prcnt > 10
        @warn "Count delta > 10% of total counts"
    end

    return gti_data
end

"""
    _gtis(lc::BinnedData)

Calls `_lc_filter_gtis` using `BinnedData` input
"""
function _gtis(lc::BinnedData)
    gti_data = _lc_filter_gtis(lc.mission, lc.instrument, lc.obsid, lc.e_range, lc.times, lc.counts, lc.gtis)

    return gti_data
end

"""
    _gtis_save(gtis, gti_dir::String)

Splits up the GTI data into a `_meta.feather` file containing non-array variables for each GTI (index, start/stop times, etc...) 
and multiple `_gti.feather` files for each GTI containing the counts and times
"""
function _gtis_save(mission::Mission, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame}, 
        instrument::Union{String,Symbol}, e_range::Tuple,
        gtis_dir::String, gtis::Dict{Int64,JAXTAM.GTIData}; log=true
    )
    @info "Saving $instrument gti data"
    mkpath(gtis_dir)
    log_entries = Dict{Union{Symbol,Int},String}()

    gti_indecies = [k for k in keys(gtis)]
    gti_starts   = [t.gti_start_time for t in values(gtis)]
    gti_example  = gtis[gti_indecies[1]]
    gti_groups   = [o.group for o in values(gtis)]

    gti_basename = string("$(gti_example.instrument)_lc_$(gti_example.bin_time)_gti")
    gtis_meta    = DataFrame(mission=string(_mission_name(gti_example.mission)), instrument=String(gti_example.instrument),
        obsid=gti_example.obsid, bin_time=gti_example.bin_time, indecies=gti_indecies, starts=gti_starts, groups=gti_groups)

    Feather.write(joinpath(gtis_dir, "$(gti_basename)_meta.feather"), gtis_meta)
    log_entries[:path_meta] = joinpath(gtis_dir, "$(gti_basename)_meta.feather")

    for index in gti_indecies
        gti_savepath = joinpath(gtis_dir, "$(gti_basename)_$(index).feather")
        Feather.write(gti_savepath, DataFrame(counts=gtis[index].counts, times=gtis[index].times))
        log_entries[index] = gti_savepath
    end

    if log
        _log_add(mission, obs_row, 
            Dict("data" =>
                Dict(:gtis =>
                    Dict(e_range =>
                        Dict(gti_example.bin_time =>
                            Dict(instrument =>
                                log_entries
                            )
                        )
                    )
                )
            )
        )
    end
    @info "Saving $instrument gti data finished"
end

"""
    _gtis_load(gti_dir, instrument, bin_time)

Loads and parses the `_meta` and `_gti` files, puts into a `GTIData` constructor, returns `Dict{Int,GTIData}`
"""
function _gtis_load(gti_dir::String, instrument, e_range::Tuple{Float64,Float64}, bin_time)
    @info "Loading $instrument gti data"
    bin_time = float(bin_time)

    gti_basename  = string("$(instrument)_lc_$(bin_time)_gti")
    gti_meta_path = joinpath(gti_dir, "$(gti_basename)_meta.feather")

    gti_meta = Feather.read(gti_meta_path)

    gti_data = Dict{Int,GTIData}()

    for gti_row in DataFrames.eachrow(gti_meta)
        gti_mission = _mission_symbol_to_type(Symbol(gti_row[:mission]))
        gti_inst    = Symbol(gti_row[:instrument])
        gti_obsid   = gti_row[:obsid]
        gti_bin_t   = gti_row[:bin_time]
        gti_idx     = gti_row[:indecies]
        gti_starts  = gti_row[:starts]
        current_gti = Feather.read(joinpath(gti_dir, "$(gti_basename)_$(gti_idx).feather"))
        gti_counts  = current_gti[:counts]
        gti_times   = current_gti[:times]
        gti_times   = gti_times[1]:gti_bin_t:gti_times[end] # Convert Array to Step Range
        gti_group   = gti_row[:groups]
        
        gti_data[gti_idx] = GTIData(gti_mission, gti_inst, gti_obsid, e_range, gti_bin_t, gti_idx,
            gti_starts, gti_counts, gti_times, gti_group)
    end

    @info "Loading $instrument gti data finished"

    return gti_data
end

"""
    gtis(mission_name::Symbol, obs_row::DataFrames.DataFrame, bin_time::Number;  =false)

Handles file management, checks to see if GTI files exist already and loads them, if files do not 
exist then the `_gits` function is ran, then the data is saved

TODO: FIX E_RANGE PROPAGATION WHEN LCURVE_DATA IS GIVEN
"""
function gtis(mission::Mission, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame}, bin_time::Number; overwrite=false,
        e_range::Tuple{Float64,Float64}=_mission_good_e_range(mission),
        lcurve_data::Dict{Symbol,BinnedData}=Dict{Symbol,BinnedData}()
    )
    obsid       = obs_row[:obsid]
    instruments = _mission_instruments(mission)
    gti_files   = _log_query(mission, obs_row, "data", :gtis, e_range, bin_time)
    gtis_path   = abspath(_obs_path_local(mission, obs_row; kind=:jaxtam), "JAXTAM", _log_query_path(; category=:data, kind=:gtis, e_range=e_range, bin_time=bin_time))

    gtis = Dict{Symbol,Dict{Int64,JAXTAM.GTIData}}()
    for instrument in instruments
        if ismissing(gti_files) || !haskey(gti_files, instrument) || overwrite
            @info "Missing gti files for $instrument"

            if !haskey(lcurve_data, instrument)
                lcurve_data = lcurve(mission, obs_row, bin_time, e_range=e_range)
            end

            @assert _recursive_first(lcurve_data).e_range == e_range "lcurve energy range does not match input: $(_recursive_first(gtis_data).e_range) != $e_range"

            @info "Selecting GTIs for $instrument"
            gtis_data = _gtis(lcurve_data[Symbol(instrument)])

            _gtis_save(mission, obs_row, instrument, e_range, gtis_path, gtis_data)

            gtis[instrument] = gtis_data
        else
            gtis[instrument] = _gtis_load(gtis_path, instrument, e_range, bin_time)
        end
    end

    if ismissing(_log_query(mission, obs_row, "meta", :countrates, e_range))
        countrates = Dict{Symbol,Float64}()
        for instrument in instruments
            countrates[instrument] = mean([sum(gti[2].counts)/gti[2].times[end] for gti in gtis[instrument]])
        end
        _log_add(mission, obs_row, Dict("meta" => Dict(:countrates => Dict(e_range => countrates))))
    end

    return gtis
end

function gtis(mission::Mission, obsid::String, bin_time::Number; overwrite=false)
    obs_row = master_query(mission, :obsid, obsid)

    return gtis(mission, obs_row, bin_time, overwrite=overwrite)
end