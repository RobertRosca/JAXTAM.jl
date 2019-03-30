# Helper functions

function savefig(mission::Mission, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame},
        e_range::Tuple{Float64,Float64}, kind::Symbol, plot_name::String;
        kwargs...
    )
    kwargs = Dict(kwargs)

    obs_log = _log_entry(; category=:images, e_range=e_range, kind=kind, file_name=string(plot_name, ".png"), kwargs...)

    plot_path = abspath(_obs_path_local(mission, obs_row; kind=:web), "JAXTAM", _log_query_path(; category=:images, kind=kind, e_range=e_range, file_name=string(plot_name, ".png"), kwargs...))

    mkpath(dirname(plot_path)); Plots.savefig(plot_path)
    @info "Saved $plot_path"

    _log_add(mission, obs_row, Dict("images"=>Dict(e_range=>obs_log)))

    return plot_path
end

function savefig(plot_data::JAXTAMData, obs_row::DataFrames.DataFrameRow{DataFrames.DataFrame},
        e_range::Tuple; plot_name::String="",
        kwargs...
    )
    mission = plot_data.mission

    types_unbinned = [ESpec]
    types_binned   = [BinnedData, PgramData]
    types_binsec   = [FFTData, FFTSpec, PulseC]

    plot_type_symbols = Dict(
        ESpec      => :espec,
        BinnedData => :lcurve,
        PgramData  => :pgram,
        FFTData    => :fspec,
        FFTSpec    => :sgram,
        PulseC     => :pulsec
    )

    plot_type = plot_type_symbols[typeof(plot_data)]

    if plot_name == ""
        plot_name = string(plot_type)
    end

    plot_path = if any([isa(plot_data, tu) for tu in types_unbinned])
        savefig(mission, obs_row, e_range, plot_type, plot_name; kwargs...)
    elseif any([isa(plot_data, tb) for tb in types_binned])
        savefig(mission, obs_row, e_range, plot_type, plot_name; bin_time=plot_data.bin_time, kwargs...)
    elseif any([isa(plot_data, tbs) for tbs in types_binsec])
        bss = plot_data.bin_time * plot_data.bin_size
        savefig(mission, obs_row, e_range, plot_type, plot_name; bin_time=plot_data.bin_time, bin_size_sec=bss, kwargs...)
    else
        @error "Data type not found in savefig configuration"
    end

    return plot_path
end

"""
    _recursive_first(instrument_data::Dict)

Recursively uses the `first()` function to pull out the first actual element 
of a dictionary

Useful when dealing with nested `insturment_data` such as `Dict{Symbol,Dict{Int64,JAXTAM.FFTData}}` 
as it will pull out a value of type `FFTData`, typically when getting other fields like `bin_time`
"""
function _recursive_first(instrument_data::Dict)
    f = first(instrument_data)[2]

    if typeof(f) <: Dict
        return _recursive_first(f)
    else
        return f
    end
end

function _plot_formatter!()
    # Plots.plot!(title_location=:left, titlefontsize=12, margin=2mm, xguidefontsize=10, yguidefontsize=10)
    return Plots.plot!(title_location=:left, margin=2mm, xguidefontsize=12, yguidefontsize=12)
end

function _log10_minor_ticks!(start=0.01, stop=4096)
    start < 1 ? start = 1/nextpow(10, 1/start) : start = prevpow(10, start)
    stop  < 1 ? stop  = 1/prevpow(10, 1/stop)  : stop  = nextpow(10, stop)

    segments = Int(log10(stop./start)+1)
    ranges   = repeat(1:9, inner=(1, segments)) # Skip powers of 10 (1, 10) to avoid redrawing

    segment_mults = 10 .^(0:segments-1) .* start
    segment_mults = repeat(segment_mults', outer=(9,1))

    log10_minors = ranges .* segment_mults

    vline!(log10_minors[2:end, :][:], color=:grey, alpha=0.2, lab="")

    return log10_minors[1, :]
end

function _save_plot_data_csv(obs_row, bin_time, subfolder, file_name; kwargs...)
    @warn "FIX THIS"
    file_dir  = joinpath(obs_row[1, :obs_path], "JAXTAM/lc/$bin_time/images/", subfolder)
    file_path = joinpath(file_dir, file_name)

    mkpath(file_dir)

    data = DataFrame(kwargs)
    data = hcat([data[i] for i in 1:size(data,2)]...)

    writedlm(file_path, data)

    println("CSV: $file_path")

    return data
end

# Lightcurve plotting functions

function plot!(data::BinnedData; lab="", size_in=(1140,400), title_append="", logy=false)
    bin_time_pow2 = Int(log2(data.bin_time))

    Plots.plot()

    plot_title = "Lightcurve - 2^$(bin_time_pow2) - $(data.bin_time) bt$title_append"

    Plots.plot!(data.times, data.counts,
        xlab="Time (s)",
        lab=lab, alpha=1, title=plot_title)

    # NOTE: Disabled GTI start/stop label
    Plots.vline!(data.gtis[2], lab="",  alpha=0.75)
    Plots.vline!(data.gtis[1], lab="", alpha=0.75)

    count_min = maximum([minimum(data.counts[data.counts .!= 0]), 0.1])
    count_max = maximum(data.counts)
    if logy
        log10_min = count_min > 1 ? prevpow(10, count_min) : 1/prevpow(10, 1/count_min)
        yticks = range(log10(log10_min), stop=log10(nextpow(10, count_max)), length=5)
    
        yticks = round.(exp10.(yticks), sigdigits=3)
    
        ylim = (log10_min, nextpow(10, count_max))
        yaxis!(yticks=yticks, ylims=ylim, yscale=:log10, ylab="Counts (log10)")
    else
        width = count_max - count_min
        yaxis!(ylims=(count_min-(width*0.1),count_max+(width*0.1)), ylab="Counts")
    end

    try
        yaxis!(yformatter=yi->round(yi, sigdigits=3))
    catch
        yaxis!(yformatter=yi->round(Int, yi))
    end

    Plots.plot!(size=size_in)

    _plot_formatter!()
    return Plots.plot!()
end

function plot(instrument_data::Dict{Symbol,JAXTAM.BinnedData}; size_in=(1140,400), save=false, title_append="")
    instruments = keys(instrument_data)

    example = _recursive_first(instrument_data)
    e_range = example.e_range
    title_append = string(" - $e_range keV", title_append)

    plt = Plots.plot()
    
    # NOTE: Disable instrument label
    for instrument in instruments
        plt = plot!(instrument_data[Symbol(instrument)]; lab="", title_append=title_append)
    end

    if save
        obs_row = master_query(example.mission, :obsid, example.obsid)
        savefig(example, obs_row, e_range)
    end

    return Plots.plot!(size=size_in)
end

function plot_groups(instrument_data::Dict{Symbol,JAXTAM.BinnedData}; size_in=(1140,400), save=false)
    instruments = keys(instrument_data)

    group_plots = Dict{Symbol,Dict{Int64,Plots.Plot}}()

    example = _recursive_first(instrument_data)
    e_range = example.e_range
    obs_row = master_query(example.mission, :obsid, example.obsid)

    for instrument in instruments
        instrument_group_data  = _group_return(instrument_data[instrument])
        instrument_group_plots = Dict{Int64,Plots.Plot}()

        availabel_groups = collect(keys(instrument_group_data))

        for group in availabel_groups
            group_data = instrument_group_data[group]

            title_append = " - $e_range keV - group $group/$(maximum(availabel_groups))"

            instrument_group_plots[group] = plot(Dict(instrument=>group_data); save=false,
                    title_append=title_append, size_in=size_in)

            if save
                savefig(group_data, obs_row, e_range; group=group)
            end
        end

        group_plots[instrument] = instrument_group_plots
    end

    return group_plots
end

# Power spectra plotting functions

function plot!(data::FFTData; title_append="", norm=:rms, rebin=(:log10, 0.01), freq_lims=(:start, :end),
        lab="", logx=true, logy=true, show_errors=true,
        size_in=(1140,900)
    )
    bin_time_pow2 = Int(log2(data.bin_time))

    avg_power = data.avg_power
    freq      = data.freq

    if freq_lims[1] == :start
        # Don't plot the 0Hz amplitude
        idx_min  = 2
        freq_min = freq[idx_min]
    else
        idx_min  = findfirst(freq .>= freq_lims[1])
        idx_min <= 0 ? idx_min=1 : idx_min=idx_min
        freq_min = freq[idx_min]
    end

    if freq_lims[2] == :end
        idx_max  = length(freq)
        freq_max = freq[end]
    else
        idx_max  = findfirst(freq .> freq_lims[2]) - 1 # Find first > max freq, take the freq before that
        freq_max = freq[idx_max]
    end

    avg_power = avg_power[idx_min:idx_max]
    freq   = freq[idx_min:idx_max]

    freq, avg_power, errors = _fspec_rebin(avg_power, freq, data.bin_count, data.bin_size, data.bin_time, rebin)
    ylab = ""
    
    if norm == :rms
        src_ctrate = mean(data.src_ctrate); bkg_ctrate = mean(data.bkg_ctrate)
        rms_factor = 1
        if bkg_ctrate == 0.0
            rms_factor = 1/src_ctrate
        else
            rms_factor = (src_ctrate + bkg_ctrate) ./ src_ctrate^2
        end

        avg_power = avg_power .- 2
        avg_power = (avg_power.*rms_factor).*freq
        errors = errors.*freq*rms_factor
        power_max = maximum(avg_power[2:end]); power_min = maximum([0.0001, minimum(abs.(avg_power[2:end]))])
        avg_power[avg_power .<=0] .= NaN
        ylab = "Amplitude (RMS*freq)"
    elseif norm == :leahym2
        errors = errors.*freq
        avg_power = (avg_power.-2).*freq
        power_max = maximum(avg_power[2:end]); power_min = minimum(abs.(avg_power[2:end]))
        avg_power[avg_power .<=0] .= NaN
        ylab = "Amplitude (Leahy - 2)*freq"
    elseif norm == :leahy
        power_max = maximum(avg_power[2:end]); power_min = minimum(avg_power[2:end])
        hline!([2], line=:dash, lab="")
        ylab = "Amplitude (Leahy)"
    else
        @error "Plot norm type '$norm' not found" 
    end

    if show_errors
        Plots.plot!(freq, avg_power, color=:black,
            yerr=errors, lab=lab)
    else
        Plots.plot!(freq, avg_power, color=:black, lab=lab)
    end

    if logx
        xaxis!(xscale=:log10, xformatter=xi->xi, xlim=(freq_min, freq[end]), xlab="Freq (Hz) - log10")
        _log10_minor_ticks!(freq_min, freq[end])
    else
        xaxis!(xlab="Freq (Hz)", xlim=(freq_min,freq_max))
    end

    if logy
        # If power_min < 1, can't use prevpow10 for ylims, hacky little fix is 1/prevpow(10, 1/power_min)
        power_min > 1 ? ylim = (prevpow(10, power_min), nextpow(10, power_max)) : ylim = (1/prevpow(10, 1/power_min), nextpow(10, power_max))
        yaxis!(yscale=:log10, yformatter=yi->round(yi, sigdigits=3), ylims=ylim, ylab="$ylab - log10")
    else
        yaxis!(ylab=ylab)
    end
    
    Plots.plot!(alpha=1,
        title="FFT - $(data.obsid) - $(data.e_range) keV - 2^$(bin_time_pow2) bt - $(data.bin_size*data.bin_time) bs - $rebin rebin - $(data.bin_count) - sections averaged$title_append")

    _plot_formatter!()
    return Plots.plot!(size=size_in)
end

function plot(instrument_data::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}}; title_append="", freq_lims=(:start, :end),
        size_in=(1140,900), norm=:rms, rebin=(:log10, 0.01), logx=true, logy=true, show_errors=true, save=false, save_csv=false)
    instruments = keys(instrument_data)

    plt = Plots.plot()
    
    example = _recursive_first(instrument_data)
    e_range = example.e_range

    for instrument in instruments # NOTE: Disabled instrument label `String(instrument)`
        plt = JAXTAM.plot!(instrument_data[Symbol(instrument)][-1]; title_append=title_append, freq_lims=freq_lims,
            norm=norm, rebin=rebin, logx=logx, logy=logy, lab="", show_errors=show_errors)

        if save_csv
            obs_row = master_query(example.mission, :obsid, example.obsid)
            # _save_plot_data_csv(obs_row, example.bin_time, # TODO: RE-ENABLE THIS WHEN FIXED
            #     "fspec/$(example.bin_size*example.bin_time)", "fspec.csv";
            #     freq=plt.series_list[1].d[:x], power=plt.series_list[1].d[:y], power_error=plt.series_list[1].d[:yerror])
        end
    end

    if save
        obs_row = master_query(example.mission, :obsid, example.obsid)
        savefig(example, obs_row, e_range)
    end

    return Plots.plot!(size=size_in)
end

function plot_groups(instrument_data::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}};
    size_in=(1140,600), norm=:rms, rebin=(:log10, 0.01), freq_lims=(:start,:end), logx=true, logy=true, show_errors=true, save=false)

    instruments = keys(instrument_data)

    group_plots = Dict{Symbol,Dict{Int64,Plots.Plot}}()

    example = _recursive_first(instrument_data)
    obs_row = master_query(example.mission, :obsid, example.obsid)
    e_range = example.e_range

    for instrument in instruments
        availabel_groups = unique([gti.group for gti in values(instrument_data[instrument])])
        availabel_groups = availabel_groups[availabel_groups .>= 0] # Excluse -1, -2, etc... for scrunched/mean FFTData

        instrument_group_plots = Dict{Int64,Plots.Plot}()

        for group in availabel_groups
            instrument_data_group = Dict{Int64,JAXTAM.FFTData}()
            gtis_in_group = [gti.gti_index for gti in values(instrument_data[instrument]) if gti.group == group]

            for gti_no in gtis_in_group
                instrument_data_group[gti_no] = instrument_data[instrument][gti_no]
            end

            title_append = " - group $group/$(maximum(availabel_groups))"

            # Have to run _scrunch_sections on group data, since it lacks the -1 indexed average amplitudes
            instrument_group_plots[group] = plot(Dict(instrument=>_scrunch_sections(instrument_data_group));
                    size_in=size_in, norm=norm, rebin=rebin, freq_lims=freq_lims, logx=logx, logy=logy, save=false, show_errors=show_errors,
                    title_append=title_append)
            
            if save
                group_data = instrument_data_group[gtis_in_group[1]]
                savefig(group_data, obs_row, e_range; group=group)
            end

        end

        group_plots[instrument] = instrument_group_plots
    end

    return group_plots
end

# Energy Spectra plotting functions

struct ESpec end

function plot_energy!(data::JAXTAM.InstrumentData; e_lims=:auto, size_in=(1140,400), save=false, title_append="")    
    if e_lims == :auto
        e_lims = _mission_good_e_range(data.mission)
    end

    exposure = sum(data.gtis[:STOP] .- data.gtis[:START])

    energies = data.events[:E]
    energies = energies[e_lims[1] .<= energies .<= e_lims[2]]

    energy_bins     = e_lims[1]:1e-2:e_lims[2]
    energy_hist_fit = OnlineStats.Hist(energy_bins)
    energy_hist     = fit!(energy_hist_fit, sort(energies))

    Plots.plot(value(energy_hist)[1], value(energy_hist)[2]./exposure;
        size=size_in, lab="",
        title="Energy Spectra - $(data.obsid) - $e_lims keV")

    Plots.xaxis!(xlab="Energy [keV]")
    Plots.yaxis!(yscale=:log10, yformatter=yi->yi, ylab="Normalised Counts - log10")

    _plot_formatter!()
    return Plots.plot!()
end

function plot_energy(instrument_data::Dict{Symbol,JAXTAM.InstrumentData}; size_in=(1140,400), save=false, title_append="")
    instruments = keys(instrument_data)

    example_lc = _recursive_first(instrument_data)

    if !haskey(example_lc.events, :E)
        @error "`instrument_data` does not contain the energy column. Are you passing data from read_cl or calibrate?"
    end

    for instrument in instruments
        plt = plot_energy!(instrument_data[instrument], save=false, title_append=title_append)
    end

    return Plots.plot!(size=size_in)
end

# Periodogram plotting functions

function plot!(data::PgramData; title_append="", rebin=(:linear, 1),
        lab="", logx=false, logy=true, size_in=(1140,600))
    bin_time_pow2 = Int(log2(data.bin_time))

    # Don't plot the 0Hz amplitude
    power = data.power
    freq  = data.freq
    power[1] = NaN
    freq[1]   = NaN

    if rebin == (:linear, 1)
        # Do... nothing
    else
        freq, power, errors = _fspec_rebin(power, freq, 1, 1, missing, rebin)
    end

    Plots.plot!(xlab="Freq (Hz)", alpha=1)

    if logy      
        power_min = maximum([0.0001 minimum(power[2:end])])
        power_max = maximum(power[2:end])
        
        power = power[:]
        power[power .<= 0] .= NaN
        
        # If power_min < 1, can't use prevpow10 for ylims, hacky little fix is 1/prevpow(10, 1/power_min)
        # removed that anyway, set ylim to 1 if power_min < 1
        # TODO: Look at/fix the manual ylim settings, since it seems to... make things worse usually
        # power_min > 1 ? ylim = (prevpow(10, power_min), nextpow(10, power_max)) : ylim = (1, nextpow(10, power_max))
        yaxis!(yscale=:log10, yformatter=yi->round(yi, sigdigits=3), ylims=(power_min,power_max))
    else
        yaxis!(ylims=(maximum([0.0, minimum(power[2:end])]), maximum(power[2:end])*1.1))
    end

    xlim = (freq[2], freq[end])
    if logx
        xaxis!(xscale=:log10, xformatter=xi->xi, xlim=xlim, xlab="Freq (Hz) - log10")
        pow10_ticks = _log10_minor_ticks!(xlim[1], xlim[2])
        xticks!(pow10_ticks)
    else
        xaxis!(xlab="Freq (Hz)", xlim=xlim)
    end
    
    # NOTE: Disabled instrument and pgram type label `lab="$lab - $(data.pg_type)"`
    Plots.plot!(freq, power, color=:black, ylab="Amplitude", lab="",
        title="Periodogram - $(data.obsid) - $(data.e_range) keV - 2^$(bin_time_pow2) bt - $rebin rebin$title_append")
        
    _plot_formatter!()
    return Plots.plot!(size=size_in)
end

function plot(instrument_data::Dict{Symbol,JAXTAM.PgramData};
        rebin=(:linear, 1), logx=true, logy=false, save=false,
        size_in=(1140,600), title_append="")

    instruments = keys(instrument_data)

    plt = Plots.plot()

    for instrument in instruments # NOTE: Disabled instrument label `string(instrument)`
        plt = plot!(instrument_data[instrument], lab="",
            rebin=rebin, size_in=size_in, logx=logx, logy=logy, title_append=title_append)
    end

    if save
        example = _recursive_first(instrument_data)
        obs_row = master_query(example.mission, :obsid, example.obsid)
        e_range = example.e_range

        savefig(example, obs_row, e_range)
    end

    return Plots.plot!(size=size_in)
end

function plot_groups(instrument_data::Dict{Symbol,Dict{Int64,JAXTAM.PgramData}};
    rebin=(:linear, 1), logx=true, logy=false, save=false,
    size_in=(1140,600), title_append="")

    instruments = keys(instrument_data)

    example = _recursive_first(instrument_data)
    obs_row = master_query(example.mission, :obsid, example.obsid)
    e_range = example.e_range

    group_plots = Dict{Symbol,Dict{Int64,Plots.Plot}}()
    for instrument in instruments
        instrument_group_plots = Dict{Int64,Plots.Plot}()
        availabel_groups = unique([gti.group for gti in values(instrument_data[instrument])])
        availabel_groups = availabel_groups[availabel_groups .>= 0] # Excluse -1, -2, etc... for scrunched/mean FFTData

        for group in availabel_groups
            title_append = " - group $group/$(maximum(availabel_groups))"

            instrument_group_plots[group] = plot(Dict(instrument=>instrument_data[instrument][group]);
                title_append=title_append, size_in=size_in, save=false)

            if save
                savefig(example, obs_row, e_range; group=group)
            end
        end

        group_plots[instrument] = instrument_group_plots
    end

    return group_plots
end

# Spectrogram plotting functions

struct FFTSpec <: JAXTAMData 
    mission    :: Mission
    instrument :: Symbol
    obsid      :: String
    e_range    :: Tuple{Float64,Float64}
    bin_time   :: Real
    bin_size   :: Int
end
FFTSpec(fs::JAXTAM.FFTData) = FFTSpec(fs.mission, fs.instrument, fs.obsid, fs.e_range, fs.bin_time, fs.bin_size)

function _plot_sgram(sgram_freq, sgram_power, sgram_bounds, sgram_groups,
        e_range, obsid, bin_time_pow2, bin_size, bin_time, rebin, size_in, disable_x=true)

    # if size(sgram_power) == 

    x_tick_approx = exp10.(-1:1:4)
    x_ticks = [findfirst(f->f>=x_lim, sgram_freq) for x_lim in x_tick_approx]
    x_ticks = unique(x_ticks)
    x_ticks = x_ticks[x_ticks.!=nothing]
    
    try
        heatmap(sgram_power, size=size_in, fill=true)
    catch err
        if typeof(err) != InterruptException
            throw(JAXTAMError("Error plotting spectrogram", :plot, err))
        else
            rethrow(err)
        end
    end
    xaxis!(xticks=x_ticks, xformatter=yi->sgram_freq[Int(yi)], xlab="Freq (Hz - log10)")

    title!("Spectrogram - $(obsid) - $e_range keV - 2^$(bin_time_pow2) bt - $(bin_size*bin_time) bs - $rebin rebin")

    hline!(sgram_bounds.+0.5, alpha=0.75, line=:dot, lab="")
    if length(sgram_bounds) < 25
        yaxis_bounds_to_group = Dict(diag([(bound,group) for bound in sgram_bounds, group in sgram_groups]))
        yaxis!(yticks=sgram_bounds, yformatter=yi->yaxis_bounds_to_group[Int(yi)], ylab="Group")
    else
        yaxis!(ylab="Groups")
    end

    _plot_formatter!()
end

function plot_sgram(fs::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}}; 
        rebin=(:log10, 0.01), size_in=(1140,600), save=false)
    instruments = keys(fs)

    sgram_instrument_plots = Dict{Symbol,Plots.Plot}()
    for instrument in instruments
        example = fs[instrument][-1]
        bin_time_pow2 = Int(log2(example.bin_time))

        sgram_freq, sgram_power, sgram_bounds, sgram_groups = fspec_rebin_sgram(fs[instrument], rebin=rebin)

        if size(sgram_power, 2) == 1
            @warn "sgram_power only has one dimension, cannot plot spectrogram"
            return sgram_instrument_plots
        end

        src_ctrate = mean(fs[instrument][-1].src_ctrate); bkg_ctrate = mean(fs[instrument][-1].bkg_ctrate)
        rms_factor = 1

        if bkg_ctrate == 0.0
            rms_factor = 1/src_ctrate
        else
            rms_factor = (src_ctrate + bkg_ctrate) ./ src_ctrate^2
        end

        sgram_power = sgram_power .- 2
        sgram_power = (sgram_power.*rms_factor).*sgram_freq
        sgram_power[sgram_power .<= 0] .= 0
        sgram_power = sgram_power'

        sgram_instrument_plots[instrument] = _plot_sgram(sgram_freq, sgram_power, sgram_bounds, sgram_groups,
            example.e_range, example.obsid, bin_time_pow2, example.bin_size, example.bin_time, rebin,
            size_in)

        if save
            obs_row = master_query(example.mission, :obsid, example.obsid)
            savefig(FFTSpec(example), obs_row, example.e_range)
        end
    end

    return sgram_instrument_plots
end

# Pulsation Check Candle Plotting

struct PulseC  <: JAXTAMData
    mission    :: Mission
    instrument :: Symbol
    obsid      :: String
    e_range    :: Tuple{Float64,Float64}
    bin_time   :: Real
    bin_size   :: Int
end
PulseC(fs::JAXTAMData) = PulseC(fs.mission, fs.instrument, fs.obsid, fs.e_range, fs.bin_time, fs.bin_size)

function _plot_pulses_candle(power, freq, power_limit, f_min, f_max)
    base_freq = freq[:, 1]

    if f_min == :start
        f_min = base_freq[2]
    end

    if f_max == :end
        f_max = base_freq[end]
    end

    # Mask for power limits, and fo excluding zero-freq
    mask = (power .>= power_limit) .* (freq .!= 0)

    power = power[mask]
    freq  = freq[mask]

    Plots.plot()
    _log10_minor_ticks!(0.01, 4096)
    Plots.plot!(freq, power, line=:sticks, lab="")

    return Plots.scatter!(freq, power, lab="", alpha=0.5)
end

function plot_pulses_candle(fs::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}};
    power_limit=30, size_in=(1140,600), save=false)

    plots_candle = Dict{Symbol,Plots.Plot}()

    example = _recursive_first(fs)

    e_range = example.e_range
    bin_time_pow2 = Int(log2(example.bin_time))
    bin_time = example.bin_time
    bin_size = example.bin_size
    obsid = example.obsid
    
    instruments = keys(fs)
    for instrument in instruments
        power = hcat([f[2].power for f in fs[instrument] if f[1] > 0]...)
        freq  = repeat(fs[instrument][-1].freq, inner=(1, size(power, 2)))

        _plot_pulses_candle(power, freq, power_limit, 0.01, :end)

        Plots.title!("Pulsations - $obsid - $e_range keV - 2^$(bin_time_pow2) bt - $(bin_size*bin_time) bs")

        xaxis!(xscale=:log10, xformatter=xi->xi, xlim=(0.01, freq[end, 1]), xlab="Freq [Hz] - log10")
        yaxis!(ylab=("Amplitude (Leahy) >= $power_limit"), ylims=(power_limit, maximum(power[2:end, :])*1.1))

        _plot_formatter!()

        plots_candle[instrument] = Plots.plot!(size=size_in)

        if save
            obs_row = master_query(example.mission, :obsid, example.obsid)
            savefig(PulseC(example), obs_row, e_range)
        end
    end

    return plots_candle
end

function plot_pulses_candle_groups(fs::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}};
    power_limit=30, size_in=(1140,600), save=false, f_lims=(0.01, :end))

    instruments = keys(fs)

    example = _recursive_first(fs)
    obs_row = master_query(example.mission, :obsid, example.obsid)
    e_range = example.e_range
    bin_time_pow2 = Int(log2(example.bin_time))
    bin_time = example.bin_time
    bin_size = example.bin_size
    obsid = example.obsid

    group_plots = Dict{Symbol,Dict{Int64,Plots.Plot}}()
    for instrument in instruments
        instrument_group_plots = Dict{Int64,Plots.Plot}()
        availabel_groups = unique([f[2].group for f in fs[instrument] if f[1] > 0])

        for group in availabel_groups
            title_append = " - group $group/$(maximum(availabel_groups))"

            power = hcat([f[2].power for f in fs[instrument] if f[2].group == group]...)
            freq  = repeat(fs[instrument][-1].freq, inner=(1, size(power, 2)))

            f_lims[1] == :start ? f_lims = (freq[1, 1], f_lims[2])    : ""
            f_lims[2] == :end   ? f_lims = (f_lims[1] , freq[end, 1]) : ""

            _plot_pulses_candle(power, freq, power_limit, f_lims[1], f_lims[2])

            Plots.title!("Pulsations - $obsid - $e_range keV - 2^$(bin_time_pow2) bt - $(bin_size*bin_time) bs$title_append")

            xaxis!(xscale=:log10, xformatter=xi->xi, xlim=f_lims, xlab="Freq [Hz] - log10")
            yaxis!(ylab=("Amplitude (Leahy) >= $power_limit"), ylims=(power_limit, maximum(power[2:end, :])*1.1))
    
            _plot_formatter!()

            instrument_group_plots[group] = Plots.plot!(size=size_in)

            if save
                savefig(PulseC(example), obs_row, e_range; group=group)
            end
        end

        group_plots[instrument] = instrument_group_plots
    end

    return group_plots
end

# Pulsation Check Spectrogram Plotting

function plot_pulses(fs::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}};
        freq_bin=10, power_limits=[10, 25, 50], size_in=(1140,600), save=false)
    
    pulsation_instrument_plots = Dict{Symbol,Plots.Plot}()
    instruments = keys(fs)
    for instrument in instruments
        example_data = fs[instrument][-1]
        bin_time_pow2 = Int(log2(example_data.bin_time))
        e_range = (config(example_data.mission).good_energy_min, config(example_data.mission).good_energy_max)

        pulsation_freq, pulsation_power, pulsation_bounds, pulsation_groups = 0, 0, 0, 0
        for p in power_limits
            rebin = (:freq_binary, freq_bin, p)
            
            pulsation_freq, pulsation_power_new, pulsation_bounds, pulsation_groups = fspec_rebin_sgram(fs[instrument]; rebin=rebin)

            if pulsation_power == 0
                pulsation_power = pulsation_power_new
            else
                pulsation_power[pulsation_power_new .!= 0] .= p
            end
        end

        rebin = (:freq_binary, freq_bin, power_limits)

        pulsation_instrument_plots[instrument] = _plot_sgram(pulsation_freq, pulsation_power', pulsation_bounds, pulsation_groups,
            e_min, e_max, example_data.obsid, bin_time_pow2, example_data.bin_size, example_data.bin_time,
            rebin, size_in, false)

        if save
            obs_row = master_query(example_data.mission, :obsid, example_data.obsid)
            _savefig_obsdir(example_data.mission, obs_row, example_data.bin_time, "pulse/$(example_data.bin_size*example_data.bin_time)", "pulse.png")
        end
    end

    return pulsation_instrument_plots
end

# Covariance plotting

function plot_fspec_cov1d(fs::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}}; size_in=(1140,600), rebin=(:log10, 0.01))
    instruments = keys(fs)

    example_data = _recursive_first(fs)
    obsid = example_data.obsid
    bin_time = example_data.bin_time
    bin_time_pow2 = Int(log2(example_data.bin_time))
    bin_size = example_data.bin_size

    cov1d_plots = Dict{Symbol,Plots.Plot}()
    for instrument in instruments
        fspec_freq, fspec_power = JAXTAM.fspec_rebin_sgram(fs[instrument]; rebin=rebin)

        fspec_power[1, :] .= 0
        src_ctrate = mean(example_data.src_ctrate); bkg_ctrate = mean(example_data.bkg_ctrate)
        rms_factor = 1
        if bkg_ctrate == 0.0
            rms_factor = 1/src_ctrate
        else
            rms_factor = (src_ctrate + bkg_ctrate) ./ src_ctrate.^2
        end
        
        fspec_power = fspec_power .- 2
        fspec_power = (fspec_power.*rms_factor).*fspec_freq        

        fspec_diag = diag(cov(fspec_power, dims=2))
        yaxis_max  = nextpow(10, maximum(fspec_diag))

        bin_count = size(fspec_power, 2)

        fspec_diag[fspec_diag .<= 1e-8] .= NaN

        Plots.plot(fspec_freq, fspec_diag, lab="", #NOTE: Disabled instrument label `lab=instrument`
            color=:black, size=size_in,
            title="FFT 1D Covariance - $(obsid) - 2^$(bin_time_pow2) bt - $(bin_size*bin_time) bs - $rebin rebin - $(bin_count) sections averaged")

        xaxis!(xscale=:log10, xformatter=xi->xi, xlab="Freq (Hz - log10)", xlims=(0.01, nextpow(2, maximum(fspec_freq))))
        yaxis!(yscale=:log10, yformatter=xi->xi, ylab="Cov (diag - log10)") # , ylims=(10, yaxis_max)
        # hline!([4000], lab="4000 - Threshold")
        
        cov1d_plots[instrument] = _plot_formatter!()
    end

    return cov1d_plots
end

function _plot_cov2d(fs::Dict{Int64,JAXTAM.FFTData}, rebin::Tuple, zoom_log10=false)
    fspec_freq, fspec_power = JAXTAM.fspec_rebin_sgram(fs; rebin=rebin) 

    fspec_power[1, :] .= 0
    src_ctrate = mean(fs[-1].src_ctrate); bkg_ctrate = mean(fs[-1].bkg_ctrate)
    rms_factor = 1
    if bkg_ctrate == 0.0
        rms_factor = 1/src_ctrate
    else
        rms_factor = (src_ctrate + bkg_ctrate) ./ src_ctrate.^2
    end
    
    fspec_power = fspec_power .- 2
    fspec_power = (fspec_power.*rms_factor).*fspec_freq    

    if rebin[1] == :linear
        fspec_cov_2d = cov(reverse(reverse(fspec_power, dims=2), dims=1), dims=2)
    elseif rebin[1] == :log10
        fspec_cov_2d = cov(fspec_power, dims=2)
    else
        error("Invalid rebin type: $(rebin[1])")
    end

    fspec_cov_2d[isnan.(fspec_cov_2d)] .= 0
    
    max_ind_2d = findmax(fspec_cov_2d)[2]
    max_freq_x = fspec_freq[max_ind_2d[1]]
    max_freq_y = fspec_freq[max_ind_2d[2]]

    if zoom_log10!=false && rebin[1]==:log10
        zoom_ind_start = findfirst(fspec_freq .>= max_freq_x - max_freq_x/2)
        zoom_ind_stop  = findfirst(fspec_freq .>= max_freq_x + max_freq_x/2)
        println(zoom_log10)
        println(zoom_ind_start)
        println(zoom_ind_stop)
        fspec_cov_2d = fspec_cov_2d[zoom_ind_start:zoom_ind_stop, zoom_ind_start:zoom_ind_stop]

        max_ind_2d = findmax(fspec_cov_2d)[2]
        max_freq_x = fspec_freq[max_ind_2d[1]]
        max_freq_y = fspec_freq[max_ind_2d[2]]

        rebin = (rebin[1], rebin[2], "$(zoom_log10) zoom")
    end

    heatmap(fspec_cov_2d, legend=false, aspect_ratio=:equal, xlab="$rebin cov")

    vline!([max_ind_2d[1]], color=:cyan, alpha=0.25, line=:dot)
    hline!([max_ind_2d[2]], color=:cyan, alpha=0.25, line=:dot)
    xticks!([max_ind_2d[1]]); yticks!([0])
    xaxis!(xformatter=xi->"$(round(max_freq_x, sigdigits=3)) Hz")
    yaxis!(yformatter=xi->"$(round(max_freq_y, sigdigits=3)) Hz")
end

function plot_fspec_cov2d(fs::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}}; size_in=(1140,600*2))
    instruments = keys(fs)

    example_data = _recursive_first(fs)
    obsid = example_data.obsid
    bin_time = example_data.bin_time
    bin_time_pow2 = Int(log2(example_data.bin_time))
    bin_size = example_data.bin_size
    fs_length = length(example_data.avg_power)

    rebin_lin = (:linear, maximum([floor(Int, fs_length/1024), 1]))
    rebin_log = (:log10, 0.01)

    for instrument in instruments
        fspec_freq, fspec_power = JAXTAM.fspec_rebin_sgram(fs[instrument]; rebin=rebin_lin)

        cov2d_linear_plot = _plot_cov2d(fs[instrument], rebin_lin)

        cov2d_log_plot = _plot_cov2d(fs[instrument], rebin_log, false)

        cov2d_log_plot_x1 = _plot_cov2d(fs[instrument], rebin_log, 1)
        cov2d_log_plot_x2 = _plot_cov2d(fs[instrument], rebin_log, 0.5)

        bin_count = size(fspec_power, 2)

        l = @layout [b c; d e] # [a{.001h}; [b c; d e]]
        dual_cov_plot = Plots.plot(
            # Plots.plot(title="FFT 2D Covariances - $(obsid) - 2^$(bin_time_pow2) bt - $(bin_size*bin_time) bs - $(bin_count) sections averaged",
            #     #annotation=(0.25, 0.5, 
            #     #"FFT 2D Covariances - $(obsid) - 2^$(bin_time_pow2) bt - $(bin_size*bin_time) bs - $(bin_count) sections averaged", 12),
            #     framestyle = :none), 
            cov2d_linear_plot,
            cov2d_log_plot,
            cov2d_log_plot_x1,
            cov2d_log_plot_x2,
            layout=l,
            size=size_in)

        _plot_formatter!()

        return dual_cov_plot
    end
end