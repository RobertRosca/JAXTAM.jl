"""
    _master_download(master_path::String, master_url::String)

Downloads (and unzips) a master table from HEASARC given its `url`
and a destination `path`
"""
function _master_download(master_path::String, master_url::String)
    if !isdir(dirname(master_path))
        mkpath(dirname(master_path))
    end
    
    @info "Downloading latest master catalog"
    Base.download(master_url, master_path)

    # Windows (used to) unzip .gz during download, unzip now if Linux
    unzip!(master_path)
end

"""
    _type_master_df!(master_df)

Slightly janky way to strongly type columns in the master table, this
needs to be done to ensure the `.feather` file is saved/read correctly

TODO: Make this less... stupid
"""
function _type_master_df!(master_df)
    nt_converter = (name=string, ra=float, dec=float, lii=float, bii=float, roll_angle=float,
        time=Dates.DateTime, end_time=Dates.DateTime, obsid=string, exposure=float, exposure_a=float,
        exposure_b=float, ontime_a=float, ontime_b=float, observation_mode=string, instrument_mode=string,
        spacecraft_mode=string, slew_mode=string, time_awarded=float, num_fpm=Meta.parse,
        processing_status=string, processing_date=Dates.DateTime,public_date=Dates.DateTime,
        processing_version=string, num_processed=Meta.parse, caldb_version=String, software_version=string,
        prnb=string, abstract=string, subject_category=string, category_code=Meta.parse,priority=string,
        country=string, data_gap=Meta.parse, nupsdout=Meta.parse, solar_activity=string, coordinated=string,
        issue_flag=Meta.parse, comments=string, satus=string,  pi_lname=string, pi_fname=string,
        cycle=Meta.parse, obs_type=string, title=string, remarks=string)

    for (name, coltype) in pairs(nt_converter)
        try
            master_df[name] = coltype.(master_df[name])
        catch e
            if typeof(e) != KeyError
                @warn e
            end
        end
    end

    return master_df
end

"""
    _master_read_tdat(master_path::String)

Reads a raw `.tdat` table from HEASARC mastertable archives,
parses the ASCII data, finds and stores column names, cleans punctuation,
converts to `DataFrame`, strongly types the columsn, and finally returns
cleaned table as `DataFrame`
"""
function _master_read_tdat(master_path::String)
    master_ascii = readdlm(master_path, '\n')

    data_start = Int(findfirst(master_ascii .== "<DATA>")[1] + 1)
    data_end   = Int(findfirst(master_ascii .== "<END>")[1] - 1)
    keys_line  = Int(findfirst(master_ascii .== "# Data Format Specification")[1] + 2)
    field_line = Int(findfirst(master_ascii .== "# Table Parameters")[1] + 2)

    # Key names are given on the keys_line, split and make into symbols for use later
    key_names = Symbol.(split(master_ascii[keys_line][11:end])) # 11:end to remove 'line[1] = '
    no_cols   = length(key_names)
    key_obsid = findfirst(key_names .== :obsid)[1]
    key_archv = findfirst(key_names .== :processing_status)

    key_types = [line[3] for line in split.(master_ascii[field_line:field_line+no_cols-1], " ")]

    master_ascii_data = master_ascii[data_start:data_end]

    master_df = DataFrame(zeros(1, no_cols), key_names)

    deleterows!(master_df, 1) # Remove row, only made to get column names

    for (row_i, row) in enumerate(master_ascii_data)
        obs_values = split(row, "|")[1:end - 1] # Split row by | delims

        if length(obs_values) != no_cols # Some rows don't have the proper no. of columns, skip them
            @warn "Skipped row $row_i due to malformed columns, ObsID: $(obs_values[key_obsid])"
            continue
        end

        df_tmp = DataFrame()

        for (itr, key) in enumerate(key_names) # Create DataFrame of key and val for row
            cleaned = replace(obs_values[itr], "," => ".. ") # Remove some punctuation, screw with CSV
            cleaned = replace(cleaned, ";" => ".. ")

            if cleaned != ""
                if key in [:time, :end_time, :processing_date, :public_date]
                    cleaned = _mjd2datetime(Meta.parse(obs_values[itr]))
                end
            else
                cleaned = missing
            end

            df_tmp[key] = cleaned
        end

        master_df = [master_df; df_tmp] # Concat
    end

    sort!(master_df, :name)

    _type_master_df!(master_df)

    return master_df
end

"""
    master(mission::Mission; update=false)

Reads in a previously created `.feather` master table for a specific `mission_name`
using a path provided by `_mission_master_url(mission))`
"""
function master_base(mission::Mission; update=false, reload_cache=false)
    path_jaxtam         = mission_paths(mission)[:jaxtam]
    path_master_tdat    = joinpath(path_jaxtam, "master.tdat")
    path_master_feather = joinpath(path_jaxtam, "master.feather")

    if (!isfile(path_master_tdat) && !isfile(path_master_tdat)) || update
        _master_download(path_master_tdat, _mission_master_url(mission))
    end
    
    if isfile(path_master_feather) && !update
        @info "Loading $path_master_feather"
        master_data = Feather.read(path_master_feather)
    elseif isfile(path_master_tdat)
        @info "Loading $(path_master_tdat)"
        master_data = _master_read_tdat(path_master_tdat)
        @info "Saving $path_master_feather"
        Feather.write(path_master_feather, master_data)
    end

    # Convert Primitive's to Dates.DateTime so they don't get returned as an Arrow timestamp 
    for (col_key, col_data) in DataFrames.eachcol(master_data, true)
        col_type = typeof(col_data) # 'col' is a Pair of (:key, col_data)
        if col_type <: Primitive{T} where T <: Timestamp
            master_data[col_key] = convert(Array{Dates.DateTime,1}, master_data[col_key])
        end
    end

    if isdefined(JAXTAM, Symbol(mission, "_master_df")) && update
        if reload_cache || update
            master(mission; cache=true, reload_cache=true, update=false) # Reload cache if master is updated
        end
    end

    return master_data
end

"""
    _add_append_publicity!(append_df::DataFrames.DataFrame, master_df::DataFrames.DataFrame)

Appends column of `Union{Bool,Missing}`, true if `public_date <=`now()`
"""
function _add_append_publicity!(mission::Mission, append_df::DataFrames.DataFrame, master_df::DataFrames.DataFrame)
    n  = Dates.now()
    pd = Array{DateTime,1}(master_df[:public_date])

    append_df[:publicity] = map(t->n>t, pd)

    return append_df
end

function _add_append_logged_vars!(mission::Mission, append_df::DataFrames.DataFrame, master_df::DataFrames.DataFrame)
    obs_count = size(append_df, 1)

    append_logged      = falses(obs_count)
    append_downloaded  = falses(obs_count)
    append_error       = falses(obs_count)
    append_error_stage = Array{String,1}(undef, obs_count)

    for (i, obs_row) in enumerate(DataFrames.eachrow(master_df))
        log_path         = _log_path(mission, obs_row)
        log_exists       = isfile(log_path)
        append_logged[i] = log_exists

        append_error_stage[i] = ""

        if log_exists
            log_contents = _log_query(mission, master_df[i, :])

            if haskey(log_contents, "meta")
                append_downloaded[i] = haskey(log_contents["meta"], :downloaded) ? log_contents["meta"][:downloaded] : false
            end

            if haskey(log_contents, "errors")
                append_error[i] = true
                append_error_stage[i] = join([string(k) for k in keys(log_contents["errors"])], ", ")
            end
        else
            continue
        end
    end

    append_df[:logged]      = append_logged
    append_df[:downloaded]  = append_downloaded
    append_df[:error]       = append_error
    append_df[:error_stage] = append_error_stage
end

function _add_append_report!(mission::Mission, append_df::DataFrames.DataFrame, master_df::DataFrames.DataFrame)
    append_report_path   = Array{String,1}(undef, size(append_df, 1))
    append_report_exists = falses(size(append_df, 1))

    full_e_range = _mission_good_e_range(mission)

    logged_indecies     = findall(append_df[:logged])
    not_logged_indecies = findall(append_df[:logged].==false)

    append_report_path[not_logged_indecies] .= ""

    for i in logged_indecies
        web_reports = _log_query(mission, master_df[i, :], "web"; surpress_warn=true)

        if ismissing(web_reports)
            append_report_path[i] = ""
        elseif haskey(web_reports, full_e_range)
            append_report_path[i]   = web_reports[full_e_range]
            append_report_exists[i] = true
        else
            append_report_path[i]   = first(web_reports)[2] # First value in dict
            append_report_exists[i] = true
        end
    end

    append_df[:report_exists] = append_report_exists
    append_df[:report_path]   = append_report_path
    return 
end

"""
    _append_gen(mission, master_df)

Runs all the `_add_append` functions, returns the full `append_df`
"""
function _append_gen(mission::Mission, master_df::DataFrames.DataFrame)
    append_df = DataFrame(obsid=master_df[:obsid])

    _add_append_publicity!(mission, append_df, master_df)
    _add_append_logged_vars!(mission, append_df, master_df)
    _add_append_report!(mission, append_df, master_df)

    return append_df
end

function master_append(mission::Mission; update=false)
    path_jaxtam         = mission_paths(mission)[:jaxtam]
    path_append_feather = joinpath(path_jaxtam, "append.feather")

    if !isfile(path_append_feather) || update
        append_df = _append_gen(mission, master_base(mission))
        Feather.write(path_append_feather, append_df)
        
        if isdefined(JAXTAM, Symbol(mission, "_master_df"))
            master(mission; cache=true, reload_cache=true, update=false) # Reload cache if append is updated
        end
        return append_df
    else
        @info "Loading $path_append_feather"
        return Feather.read(path_append_feather)
    end
end

"""
    master(::Mission; cache::Bool, reload_cache::Bool, update::Bool)

Wrapper function which deals with downloading and updating both the `master` and `append` files, as well as 
the caching functionality

Will download and set up the master and append files when first called

Subsequent calls will use the downloaded files

If you want to update the master table, simply call `master(::Mission, update=true)`
"""
function master(mission::Mission; cache::Bool=true, reload_cache::Bool=false, update::Bool=false)
    master_df_var = Symbol(_mission_name(mission), "_master_df")

    if update
        master_base(mission; update=true)
        master_append(mission; update=true)
    end

    if cache
        if reload_cache
            @info "Reloading master_append cache"
        end
        
        if isdefined(JAXTAM, master_df_var) && !reload_cache
            @assert typeof(getproperty(JAXTAM, master_df_var)) == DataFrames.DataFrame
            return getproperty(JAXTAM, master_df_var)
        else cache
            master_df = master_base(mission)
            append_df = master_append(mission)

            # If reloading cache after master update, then current append_df may not be the correct size
            # for some odd edge cases, this is a hacky fix
            if size(master_df, 1) != size(append_df, 1)
                append_df = master_append(mission; update=true)
            end

            master_df_a = join(master_df, append_df, on=:obsid)

            eval(:(global const $master_df_var = $master_df_a))

            return getproperty(JAXTAM, master_df_var)
        end
    else
        master_df = master_base(mission)
        append_df = master_append(mission)
        master_df_a = join(master_df, append_df, on=:obsid)
        return master_df_a
    end
end