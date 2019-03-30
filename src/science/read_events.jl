struct InstrumentData <: JAXTAMData
    mission    :: Mission
    instrument :: Symbol
    obsid      :: String
    events     :: DataFrames.DataFrame
    gtis       :: DataFrames.DataFrame
    start      :: Number
    stop       :: Number
    src_ctrate :: Union{Missing, Number}
    bkg_ctrate :: Union{Missing, Number}
    header     :: DataFrames.DataFrame
end

"""
    _read_fits_hdu(fits_file::FITS, hdu_id::String; cols="auto")

Reads the HDU `hdu_id` from the loaded FITS file `fits_file`, returns the HDU data

Cannot read `BitArray` type columns due to `FITSIO` limitations
"""
function _read_fits_hdu(fits_file::FITS, hdu_id::String; cols="auto")
    # Columns read individually, instead of using the FITSIO function to read
    # all the header, as if header contains non-supported column then that function
    # just fails. This avoids that problem and throws warning for non-supported columns
    if cols == "auto"
        fits_cols_events = Array{String,1}

        fits_cols_events = FITSIO.colnames(fits_file[hdu_id])
    else
        fits_cols_events = cols
    end

    fits_hdu_data = DataFrame()

    for col in fits_cols_events
        try
            fits_col_data = read(fits_file[hdu_id], col)

            if ndims(fits_col_data) == 1
                fits_hdu_data[Symbol(col)] = fits_col_data
            end
        catch error
            @warn "$col not supported by FITSIO, skipped - $error"
        end
    end

    return fits_hdu_data
end

"""
    _read_fits_event(fits_path::String)

Reads the standard columns for timing analysis ("TIME", "PI", "GTI") from a FITS file, 
returns `InstrumentData` type filled with the relevant data
"""
function _read_fits_event(mission::Mission, fits_path::String)
    @info "Loading $fits_path"
    fits_file   = FITS(fits_path)

    fits_header = read_header(fits_file["EVENTS"])

    instrument_name = Symbol(fits_header["INSTRUME"])
    fits_telescope  = lowercase(fits_header["TELESCOP"])
    @assert fits_telescope == _mission_name(mission) "Mission names do not match"
    fits_telescope  = _mission_symbol_to_type(Symbol(lowercase(fits_telescope)))
    
    fits_events_df = _read_fits_hdu(fits_file, "EVENTS"; cols=String["TIME", "PI"])
    
    fits_gtis_df = _read_fits_hdu(fits_file, "GTI")
    
    fits_obsid = fits_header["OBS_ID"]
    fits_start = fits_header["TSTART"]
    fits_stop  = fits_header["TSTOP"]

    fits_header_df = DataFrame()

    for (i, key) in enumerate(keys(fits_header))
        key = Symbol(replace(key, '-' => '_')) # Colnames with `-` don't behave well with DataFrames/Feather
        if typeof(fits_header[i]) == Nothing || fits_header[i] == ""
            fits_header_df[Symbol(key)] =  "empty" # Have to write something, or Feahter.jl errors saving ""
        else
            fits_header_df[Symbol(key)] = fits_header[i]
        end
    end

    # Hacky fix to not create DF with col type of missing
    src_rt = Array{Union{Float64,Missing},1}(undef, 1)
    bkg_rt = Array{Union{Float64,Missing},1}(undef, 1)

    total_gti_time    = sum(fits_gtis_df[:STOP] .- fits_gtis_df[:START])
    total_event_count = size(fits_events_df, 1)

    src_rt[1] = total_event_count/total_gti_time
    bkg_rt[1] = missing # TODO: Add in support for background count rates

    fits_header_df[:SRC_RT] = src_rt
    fits_header_df[:BKG_RT] = bkg_rt
    fits_header_df[:TELESCOP] = lowercase.(fits_header_df[:TELESCOP])

    close(fits_file)
    @info "Finished FITS loading"

    return InstrumentData(fits_telescope, instrument_name, fits_obsid, 
        fits_events_df, fits_gtis_df, fits_start, fits_stop,
        src_rt[1], bkg_rt[1],
        fits_header_df)
end

"""
    read_cl_fits(mission::Symbol, obs_row::DataFrames.DataFrame)

Reads in FITS data for an observation, returns a `Dict{Symbol,InstrumentData}`, with the 
symbol as the instrument name. So `instrument_data[:XTI]` works for NICER, and either 
`instrument_data[:FPMA]` or `instrument_data[:FPMB]` work for NuSTAR
"""
function read_cl_fits(mission::Mission, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame}, file_path::String)  
    if isfile(file_path)
        @info "Found: $file_path"
    elseif isfile(file_path[1:end-3]) # Check for files ending in .evt.gz too
        file_path = file_path[1:end-1]
        @info "Found: $file_path"
    else
        throw(SystemError("attempted to open:\n\t$file_path\n\t$(file_path[1:end-3])\n", 2, nothing))
    end
    
    instrument_data = _read_fits_event(mission, file_path)
    return Dict(instrument_data.instrument => instrument_data)
end

"""
    _save_cl_feather(feather_dir::String, instrument::Union{String,Symbol},
        fits_events_df::DataFrame, fits_gtis_df::DataFrame, fits_meta_df::DataFrame)

Due to Feather file restrictions, cannot save all the event and GTI data in one, so 
they are split up into three files: `events`, `gtis`, and `meta`. The `meta` file contains 
just the mission name, obsid, and observation start and stop times
"""
function _save_cl_feather(mission::Mission, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame}, instrument::Union{String,Symbol},
        feather_dir::String, fits_events_df::DataFrames.DataFrame, fits_gtis_df::DataFrames.DataFrame, fits_meta_df::DataFrames.DataFrame;
        log=true)
    @info "Saving '$instrument' cl files to $feather_dir"
    mkpath(feather_dir)

    path_events = joinpath(feather_dir, "$(instrument)_events.feather")
    path_gtis   = joinpath(feather_dir, "$(instrument)_gtis.feather")
    path_meta   = joinpath(feather_dir, "$(instrument)_meta.feather")

    Feather.write(path_events, fits_events_df)
    Feather.write(path_gtis, fits_gtis_df)
    Feather.write(path_meta, fits_meta_df)

    if log
        _log_add(mission, obs_row, 
            Dict("data" =>
                Dict(:feather_cl =>
                    Dict(instrument =>
                        Dict(
                            :path_events => path_events,
                            :path_gtis   => path_gtis,
                            :path_meta   => path_meta,
                        )
                    )
                )
            )
        )
    end

    @info "Saved $instrument feather file"

    return path_events, path_gtis, path_meta
end

function _read_cl_feather(path_events::String, path_gtis::String, path_meta::String)
    data_events = Feather.read(path_events)
    data_gtis   = Feather.read(path_gtis)
    data_meta   = Feather.read(path_meta)

    meta_missn = _mission_symbol_to_type(Symbol(data_meta[1, :TELESCOP]))
    meta_inst  = data_meta[1, :INSTRUME]
    meta_obsid = data_meta[1, :OBS_ID]
    meta_start = data_meta[1, :TSTART]
    meta_stop  = data_meta[1, :TSTOP]
    meta_srcrt = data_meta[1, :SRC_RT]
    meta_bkgrt = data_meta[1, :BKG_RT]

    return InstrumentData(meta_missn, Symbol(meta_inst), string(meta_obsid),
        data_events, data_gtis, meta_start, meta_stop,
        meta_srcrt, meta_bkgrt, data_meta)
end

"""
    read_cl(mission::Symbol, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame}; overwrite=false)

Attempts to read saved (feather) data, if none is found then the `read_cl_fits` function is ran 
and the data is saved with `_save_cl_feather` for future use
"""
function read_cl(mission::Mission, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame}; overwrite=false)
    obsid       = obs_row[:obsid]
    instruments = _mission_instruments(mission)

    cl_files_feather = _log_query(mission, obs_row, "data", :feather_cl)
    cl_files_raw     = _obs_files_cl(mission, obs_row)

    total_src_ctrate = 0.0
    instrument_data = Dict{Symbol,JAXTAM.InstrumentData}()
    for instrument in instruments
        if ismissing(cl_files_feather) || !haskey(cl_files_feather, instrument) || overwrite
            file_path = cl_files_raw[instrument]

            @info "Missing feather_cl for $instrument"
            @info "Processing $instrument cl fits"

            current_instrument = read_cl_fits(mission, obs_row, file_path)[instrument]
            
            if size(current_instrument.events, 1) == 0
                @warn "No events found in $instrument observation $obsid"
                throw(JAXTAMError("Unable to construct InstrumentData from empty DataFrame. Observation likely has no events", :read_cl))
            else
                cl_path = abspath(_obs_path_local(mission, obs_row; kind=:jaxtam), "JAXTAM", _log_query_path(; category=:data, kind=:feather_cl))
                path_events, path_gtis, path_meta = _save_cl_feather(mission, obs_row, instrument, cl_path,
                    current_instrument.events, current_instrument.gtis, current_instrument.header)
            end
            
            total_src_ctrate += current_instrument.header[1, :SRC_RT]
            instrument_data[instrument] = current_instrument
        else
            @info "Loading feather cl files for '$instrument'"
            feather_paths = cl_files_feather[instrument]
            instrument_data[instrument] = _read_cl_feather(feather_paths[:path_events], feather_paths[:path_gtis], feather_paths[:path_meta])
        end
    end

    if ismissing(_log_query(mission, obs_row, "meta", :countrates, :raw))
        total_src_ctrate = total_src_ctrate/length(instruments)
        _log_add(mission, obs_row, Dict("meta" => Dict(:countrates => Dict(:raw => total_src_ctrate))))
    end

    return instrument_data
end

function read_cl(mission::Mission, master_df::DataFrames.DataFrame, obsid::String; overwrite=false)
    obs_row = master_query(master_df, :obsid, obsid)

    return read_cl(mission, obs_row; overwrite=overwrite)
end

function read_cl(mission::Mission, obsid::String; overwrite=false)
    return read_cl(mission, master(mission), obsid; overwrite=overwrite)
end