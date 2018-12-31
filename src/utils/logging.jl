function _log_entry(; kwargs...)
    template = OrderedDict(
        :category     => String,
        :e_range      => Union{String,Missing},
        :bin_time     => Union{String,Missing},
        :bin_size_sec => Union{Real,Missing},
        :kind         => Symbol,
        :kind_order   => Union{Int,Missing},
        :value        => Union{Any,Missing},
        :group        => Union{Int,Missing},
        :file_name    => Union{String,Missing},
        :path         => Union{String,Missing},
        :notes        => Union{String,Missing},
    )

    kind_order = Dict{Symbol,Int}(
        :lcurve => 1,
        :pgram  => 2,
        :fspec  => 3,
        :sgram  => 4,
        :pulsec => 5,
        :espec  => 6,
        :DELETE => 0,
        :feather_cl => -1,
        :gtis   => -1
    )

    t = Array{Type,1}(collect(values(template)))
    n = collect(keys(template))

    log_entry = DataFrame(t, n, 1)
    for col in names(log_entry)
        if haskey(kwargs, col)
            if col == :e_range
                v = "$(kwargs[:e_range][1])_$(kwargs[:e_range][2])"
            elseif col == :bin_time
                v = "2pow$(Int(log2(kwargs[:bin_time])))"
            elseif col == :file_name
                if haskey(kwargs, :group)
                    v = string("groups/", kwargs[:group], "_", kwargs[:file_name])
                else
                    v = kwargs[:file_name]
                end
            else
                v = kwargs[col]
            end

            log_entry[col] = v
        elseif col == :kind_order
            log_entry[:kind_order] = kind_order[kwargs[:kind]]
        else
            log_entry[col] = missing
        end
    end

    if !haskey(kwargs, :path)
        path = "./"
        [path=joinpath(path, string(v)) for (k, v) in zip(names(log_entry), convert(Array, log_entry)) if !in(k, [:kind_order, :value, :group, :path, :notes]) && !ismissing(v)]
    end

    if !haskey(kwargs, :file_name)
        path = string(path, "/") # No file name? Make path a directory path
    end
    
    log_entry[:path] = path

    return log_entry
end

function _log_query_path(; kwargs...)
    return _log_entry(; kwargs...)[1, :path]
end

function _log_gen(mission_name::Symbol, obs_row::DataFrames.DataFrame)
    path = joinpath(obs_row[1, :obs_path], "JAXTAM/obs_log.jld2")

    if !ispath(path)
        mkpath(dirname(path))
    end

    meta_info = Dict{Any,Any}(
        :obs_row => obs_row
    )

    save(path, Dict{Any,Any}("meta"=>meta_info))
end

function _log_read(mission_name::Symbol, obs_row::DataFrames.DataFrame)
    path = joinpath(obs_row[1, :obs_path], "JAXTAM/obs_log.jld2")

    if !isfile(path)
        _log_gen(mission_name, obs_row)
    end

    Dict{Any,Any}(load(path))
end

function _merge(log::Dict, k::Union{Symbol,Tuple}, entry::Dict{T,DataFrames.DataFrame}) where T <: Any
    template_log = _log_entry(; category="DELETE", kind=:DELETE)
    deleterows!(template_log, 1)

    append!(template_log, log[k])
    append!(template_log, entry[k])
    unique!(template_log)

    log[k] = template_log

    return log
end

function _merge(log::Dict, k::Union{Symbol,Tuple}, entry::Dict)
    log   = convert(Dict{Any,Any}, log)
    entry = convert(Dict{Any,Any}, entry)
    
    merge!(log, entry)
    
    return log
end

function _log_add_recursive(log::Dict, entry::Dict)
    log = Dict{Any,Any}(log)
    for (k, v) in entry
        if typeof(v) <: Dict
            if haskey(log, k)
                log[k] = _log_add_recursive(log[k], entry[k])
            else
                log[k] = entry[k]
                return log
            end
        else
            return _merge(log, k, entry)
        end
    end

    return log
end

function _log_add(mission_name::Symbol, obs_row::DataFrames.DataFrame, entry::Dict{String,T}) where T <: Any
    log = _log_read(mission_name, obs_row)

    log_new = _log_add_recursive(log, entry)
    save(joinpath(obs_row[1, :obs_path], "JAXTAM/obs_log.jld2"), log_new)
    return log_new
end

function _log_query(mission_name::Symbol, obs_row::DataFrames.DataFrame, args...)
    log = _log_read(mission_name, obs_row)

    reply = log
    for (i, key) in enumerate(args)
        if reply isa Dict
            if haskey(reply, key)
                reply = reply[key]
            else
                @warn "Key '$key' not found in log, available: $(keys(reply))"
                return missing
            end
        else
            @warn "Entry not Dict, gone too deep? Accessing: $(args[1:i]) gives log at '$key' as $(typeof(reply)) not Dict"
            return missing
        end
    end

    return reply
end