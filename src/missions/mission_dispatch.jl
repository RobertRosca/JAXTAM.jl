"""
_clean_path_dots(dir)

FTP directories use hidden dot folders frequently, function removes dots for local use
"""
function _clean_path_dots(dir)
    return abspath(replace(dir, "." => "")) # Remove . from folders to un-hide them
end

"""
    Mission

Each mission must have some defined functions:

`_mission_name(<:Mission)` - returns the name of the mission

`_mission_master_url(<:Mission)` - returns the url to the HEASARC master table (.tdat.gz) file

`_mission_paths(<:Mission)` - returns a `NamedTuple` with entries corresponding to folder and file paths:
* download - path the observation folder will be downloaded to from the FTP
* jaxtam   - path JAXTAM-generated files will be saved to
* web      - path the web report files will be saved to
* rmf      - path to the **local** copy of the mission RMF file from caldb

`_obs_path_server(<:Mission, obs_row::Union{DataFrame,DataFrameRow{DataFrame}})` - returns the **server-side** path to an observation

`_mission_good_e_range(<:Mission)` - returns a `Tuple{Float64,Float64}` of the **low** and **high** good energy ranges, respectively

`_mission_instruments(<:Mission)` - returns a `Array{Symbol1,}` with each element as a symbol of the instrument name, e.g. [:FPMA, :FPMB] for NuSTAR

Look at /src/missions/nicer.jl for an example. Put any custom mission functions in /src/missions/custom/ folder.
"""
abstract type Mission end

if !isdir(joinpath(JAXTAM.__sourcedir__, "src/missions/custom"))
    mkdir(joinpath(JAXTAM.__sourcedir__, "src/missions/custom"))
end

[include("base/$mission")   for mission in readdir(joinpath(JAXTAM.__sourcedir__,   "src/missions/base"))]
[include("custom/$mission") for mission in readdir(joinpath(JAXTAM.__sourcedir__, "src/missions/custom"))]

_mission_log(log_path::String=joinpath(JAXTAM.__sourcedir__, "mission_paths.json")) = JSON.parsefile(log_path)

mission_paths(; log_path::String=joinpath(JAXTAM.__sourcedir__, "mission_paths.json")) = _mission_log(log_path)

function mission_paths(mission::Mission; log_path::String=joinpath(JAXTAM.__sourcedir__, "mission_paths.json"), overwrite::Bool=false)
    mission_name = _mission_name(mission)

    if !isfile(log_path)
        log_dict = Dict()
    else
        log_dict = JSON.parsefile(log_path)
    end

    if !haskey(log_dict, mission_name) || overwrite
        @info "Mission not found in $log_path, please enter paths:"
        print("Download path: ");                 path_download = readline()
        print("JAXTAM (processed data) path: ");  path_jaxtam   = readline()
        print("Web (html reports) path: ");       path_web      = readline()
        print("RMF (caldb mission file) path: "); path_rmf      = readline()
        log_dict[mission_name] = (download=path_download, jaxtam=path_jaxtam, web=path_web, rmf=path_rmf)
        write(log_path, JSON.json(log_dict, 4))
        @info "Wrote to $log_path \n Add custom keys in to JSON file manually if required"
    end

    # Return named tuple to prevent user from attempting to change the paths
    nt_keys = Tuple(Symbol.(collect(keys(log_dict[mission_name]))))
    return NamedTuple{nt_keys}((values(log_dict[mission_name])))
end

function _obs_path_local(mission::Mission, obs_row::Union{DataFrames.DataFrame,DataFrameRow{DataFrame}}; kind::Symbol)
    return joinpath(mission_paths(mission)[kind], _clean_path_dots(_obs_path_server(mission, obs_row))[2:end])
end

function mission_summary(mission::Mission)
    name = string("Name: ", _mission_name(mission))
    master_url = string("URL: ", _mission_master_url(mission))
    paths = string("Paths: \n\t\t", join(["$k: $v" for (k,v) in pairs(mission_paths(mission))], "\n\t\t"))
    e_range = string("Good energy range: ", _mission_good_e_range(mission), " keV")
    instruments = string("Instruments: ", join([i for i in _mission_instruments(mission)], ", "))

    println("$name\n\t$master_url\n\t$paths\n\t$e_range\n\t$instruments")
end

show(mission::Mission) = JAXTAM.mission_sumary(mission)

function _mission_symbol_to_type(mission_symbol::Symbol)
    if isdefined(JAXTAM, mission_symbol)
        return getproperty(JAXTAM, mission_symbol)
    end
end