function _plot_fspec_grid(fs::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}},
        obs_row)
    example = _recursive_first(fs)

    plt_1 = plot(fs; norm=:rms,   save=false, save_csv=true)
    plt_2 = plot(fs; norm=:leahy, save=false)

    plt_3 = plot(fs, norm=:leahy, save=false, freq_lims=(0, 1),     rebin=(:linear, 1),
        logx=false, logy=false, title_append=" - 0 to 1 Hz")
    plt_4 = plot(fs, norm=:leahy, save=false, freq_lims=(1, :end),  rebin=(:linear, 1),
        logx=false, logy=false, title_append=" - 1 to :end Hz")
    plt_4 = plot(fs, norm=:leahy, save=false, freq_lims=(50, :end), rebin=(:linear, 1),
        logx=false, logy=false, title_append=" - 50 to :end Hz")

    Plots.plot(plt_1, plt_2, plt_3, plt_4, layout=grid(4,1), size=(1140,600*4))

    savefig(example, obs_row, example.e_range; plot_name="grid")
end

# Sort of works but doesn't
# macro error_wrapper(mission, obs_row, func)
#     local escaped_obs_row = esc(obs_row)

#     return quote
#         local _wrapped_result = nothing
#         try
#             _wrapped_result = $func
#         catch err
#             if typeof(err) == JAXTAMError
#                 _wrapped_result = err
#                 _log_add($mission, $escaped_obs_row, Dict{String,Any}("errors"=>Dict(err.step=>err)))
#                 @warn err
#             else
#                 rethrow(err)
#             end
#         end

#         _wrapped_result
#     end
# end

function _error_wrapper(func, args...; mission=nothing, obs_row=nothing, kwargs...)
    if mission == nothing
        # Assume first arg is mission
        @assert args[1] isa Mission
        mission = args[1]
    end

    if obs_row == nothing
        # Assume second arg is obs_row or obsid
        if args[2] isa String
            obs_row = master_query(mission, :obsid, args[2])
        elseif args[2] isa DataFrameRow
            obs_row = args[2]
        end
    end

    try
        return func(args...; kwargs...)
    catch err
        if typeof(err) == JAXTAMError
            _log_add(mission, obs_row, Dict{String,Any}("errors"=>Dict(err.step=>err)))
            @warn err
            return nothing
        else
            rethrow(err)
        end
    end
end

"""
    report(::Mission, ::DataFrameRow; e_range::Tuple{Float64,Float64}, overwrite::Bool, nuke::Bool, update_masterpage::Bool)

Generates a report for the default energy range

Creates plots for:
    * Lightcurve (+ grouped lightcurves)
    * Periodigram (+ grouped periodograms)
    * Power Spectra
        * :rms, full range, log-rebinned, log-log plot
        * :leahy, full range, log-rebinned, log-log plot
        * :leahy, 0 to 1 Hz, no rebin, linear-linear plot
        * :leahy, 1 to end Hz, no rebin, linear-linear plot
        * :leahy, 50 to end Hz, no regin, linear-linear plot
    * Spectrogram
    * Pulsation search plot

Produces HTML report page

Updates the homepage
"""
function report(mission, obs_row; e_range=_mission_good_e_range(mission), overwrite=false, nuke=false, update_masterpage=true)
    path_jaxtam = abspath(_obs_path_local(mission, obs_row; kind=:jaxtam), "JAXTAM")
    path_web    = abspath(_obs_path_local(mission, obs_row; kind=:web), "JAXTAM")

    if nuke
        GC.gc() # Required due to Feather.jl loading files lazily, meaning they can't be removed from disk
                # until garbace collection runs and un-lazily-loads them
        ispath(path_jaxtam) ? rm(path_jaxtam, recursive=true) : false
        ispath(path_web)    ? rm(path_web,    recursive=true) : false
    end

    if ismissing(_log_query(mission, obs_row, "meta", :downloaded; surpress_warn=true)) || !_log_query(mission, obs_row, "meta", :downloaded)
        _error_wrapper(download, (mission, obs_row))
    end

    if !ismissing(JAXTAM._log_query(mission, obs_row, "errors", :read_cl; surpress_warn=true))
        @warn "Error logged at :read_cl stage, no files to be analysed, skipping report gen"
        return nothing
    end

    images = _log_query(mission, obs_row, "images", e_range)

    img_count_groupless = ismissing(images) ? 0 : size(filter(x->ismissing(x[:group]), images), 1)
    # Expect five 'groupless' plots: lightcurve, periodogram, powerspectra, spectrogram, pulsations
    if img_count_groupless < 5 || overwrite
        lc = _error_wrapper(JAXTAM.lcurve, mission, obs_row, 2.0^0; e_range=e_range)
        _error_wrapper(JAXTAM.plot, lc; save=true, mission=mission, obs_row=obs_row)
        _error_wrapper(JAXTAM.plot_groups, lc; save=true, size_in=(1140,400/2), mission=mission, obs_row=obs_row)
        pg = _error_wrapper(JAXTAM.pgram, lc, mission=mission, obs_row=obs_row)
        _error_wrapper(JAXTAM.plot, pg; save=true, mission=mission, obs_row=obs_row)
        pg = _error_wrapper(JAXTAM.pgram, lc; per_group=true, mission=mission, obs_row=obs_row)
        _error_wrapper(JAXTAM.plot_groups, pg; save=true, size_in=(1140,600/2), mission=mission, obs_row=obs_row)
        lc = nothing; pg = nothing; GC.gc()

        lc = _error_wrapper(JAXTAM.lcurve, mission, obs_row, 2.0^-13; e_range=e_range)
        gtis = _error_wrapper(JAXTAM.gtis, mission, obs_row, 2.0^-13; lcurve_data=lc, e_range=e_range)
        lc = nothing; 

        fs = _error_wrapper(JAXTAM.fspec, mission, obs_row, 2.0^-13, 128; gtis_data=gtis, e_range=e_range)
        if fs != nothing # Only try plotting if fspec didn't fail
            @info "Plotting fspec grid"
            _error_wrapper(JAXTAM._plot_fspec_grid, fs, obs_row, mission=mission, obs_row=obs_row)
            @info "Plotting fspec groups"
            _error_wrapper(JAXTAM.plot_groups, fs; save=true, size_in=(1140,600/2), mission=mission, obs_row=obs_row)
            @info "Plotting sgram"
            _error_wrapper(JAXTAM.plot_sgram, fs;  save=true, size_in=(1140,600), mission=mission, obs_row=obs_row)
            @info "Plotting pulses"
            _error_wrapper(JAXTAM.plot_pulses_candle, fs; save=true, size_in=(1140,600/2), mission=mission, obs_row=obs_row)
            @info "Plotting pulses groups"
            _error_wrapper(JAXTAM.plot_pulses_candle_groups, fs; save=true, size_in=(1140,600/2), mission=mission, obs_row=obs_row)
            fs = 0; GC.gc()
        end
    end

    subpage_path = _error_wrapper(_webgen_subpage, mission, obs_row; e_range=e_range)

    if update_masterpage
        webgen_mission(mission)
    end

    return subpage_paths
end

function report(mission::Mission, obsid::String; e_range=_mission_good_e_range(mission), overwrite=false, nuke=false, update_masterpage=true)
    obs_row = master_query(mission, :obsid, obsid)

    return report(mission, obs_row; e_range=e_range, overwrite=overwrite, nuke=nuke, update_masterpage=update_masterpage)
end

"""
    report_all(::Mission, ::DataFrameRow; e_ranges=[(0.2,12.0), (2.0,10.0), (0.2,2.0)], overwrite::Bool, nuke::Bool, update_masterpage::Bool)

Calls `report` with three default energy ranges
"""
function report_all(mission::Mission, obs_row::DataFrames.DataFrameRow; e_ranges=[(0.2,12.0), (2.0,10.0), (0.2,2.0)], overwrite=false, nuke=false, update_masterpage=true)
    if nuke
        path_jaxtam = abspath(_obs_path_local(mission, obs_row; kind=:jaxtam), "JAXTAM")
        path_web    = abspath(_obs_path_local(mission, obs_row; kind=:web), "JAXTAM")
    
        GC.gc() # Required due to Feather.jl loading files lazily, meaning they can't be removed from disk
                # until garbace collection runs and un-lazily-loads them
        ispath(path_jaxtam) ? rm(path_jaxtam, recursive=true) : false
        ispath(path_web)    ? rm(path_web,    recursive=true) : false
    end

    for e_range in e_ranges
        report(mission, obs_row; e_range=e_range, overwrite=overwrite, update_masterpage=update_masterpage)
    end

    # Run it again to re-generate report pages so that they all have each others links in them
    for e_range in e_ranges
        println("Report path:\n\t$(report(mission, obs_row; e_range=e_range, overwrite=false))")
    end
end

function report_all(mission::Mission, obsid::String; e_ranges=[(0.2,12.0), (2.0,10.0), (0.2,2.0)], overwrite=false, nuke=false, update_masterpage=true)
    obs_row = _master_query(mission, :obsid, obsid)

    return report_all(mission, obs_row; e_ranges=e_ranges, overwrite=overwrite, nuke=nuke, update_masterpage=update_masterpage)
end