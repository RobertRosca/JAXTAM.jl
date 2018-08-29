# Helper functions

function _savefig_obsdir(obs_row, mission_name, obsid, bin_time, bin_size, fig_name)
    if bin_size == NaN
        plot_dir = joinpath(obs_row[1, :obs_path], "JAXTAM/lc/$bin_time/images/$(bin_time*bin_size)/")
    else
        plot_dir = joinpath(obs_row[1, :obs_path], "JAXTAM/lc/$bin_time/images/")
    end

    mkpath(plot_dir)
    savefig(joinpath(plot_dir, fig_name))
end

function _savefig_obsdir(mission_name, obsid, bin_time, bin_size, fig_name)
    obs_row = master_query(mission_name, :obsid, obsid)
    _savefig_obsdir(obs_row, mission_name, obsid, bin_time, bin_size, fig_name)
end

function _pull_data(instrument_data::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}})
    inst1 = collect(keys(instrument_data))[1]
    gti1  = collect(keys(instrument_data[inst1]))[1]
    row1  = instrument_data[inst1][gti1]

    return row1
end

function _pull_data(instrument_data::Dict{Symbol,JAXTAM.BinnedData})
    inst1 = collect(keys(instrument_data))[1]
    row1  = instrument_data[inst1]

    return row1
end

# Lightcurve plotting functions

function plot!(data::BinnedData; lab="", size_in=(1140,400), save_plt=true, title_append="")
    bin_time_pow2 = Int(log2(data.bin_time))

    Plots.plot()

    plot_title = "Lightcurve - 2^$(bin_time_pow2) - $(data.bin_time) bt$title_append"

    Plots.plot!(data.times, data.counts,
        xlab="Time (s)", ylab="Counts (log10)",
        lab=lab, alpha=1, title=plot_title)

    Plots.vline!(data.gtis[:, 2], lab="GTI Stop",  alpha=0.75)
    Plots.vline!(data.gtis[:, 1], lab="GTI Start", alpha=0.75)

    count_min = maximum([minimum(data.counts[data.counts != 0]), 0.1])
    count_max = maximum(data.counts)
    log10_min = count_min > 1 ? prevpow(10, count_min) : 1/prevpow(10, 1/count_min)
    yticks = range(log10(log10_min), stop=log10(nextpow(10, count_max)), length=5)

    yticks = round.(exp10.(yticks), sigdigits=3)

    ylim = (log10_min, nextpow(10, count_max))
    yaxis!(yscale=:log10, yticks=yticks, ylims=ylim)

    try
        yaxis!(yformatter=yi->round(Int, yi))
    catch
        yaxis!(yformatter=yi->yi)
    end

    if(save_plt)
        _savefig_obsdir(data.mission, data.obsid, data.bin_time, NaN, "$(data.instrument)_lcurve.png")
    end

    return Plots.plot!(size=size_in)
end

function plot(instrument_data::Dict{Symbol,JAXTAM.BinnedData}; size_in=(1140,400), save_plt=true, title_append="")
    instruments = keys(instrument_data)

    plt = Plots.plot()
    
    for instrument in instruments
        plt = plot!(instrument_data[Symbol(instrument)]; lab=String(instrument), save_plt=false, title_append=title_append)
    end

    if(save_plt)
        example_lc = _pull_data(instrument_data)
        _savefig_obsdir(example_lc.mission, example_lc.obsid, example_lc.bin_time, NaN, "lcurve.png")
    end

    return Plots.plot!(size=size_in)
end

function plot_orbits(instrument_data::Dict{Symbol,JAXTAM.BinnedData}; size_in=(1140,400), save_plt=true)
    instruments = keys(instrument_data)

    orbit_plots = Dict{Symbol,Dict{Int64,Plots.Plot}}()

    example_gti = _pull_data(instrument_data)
    obs_row     = master_query(example_gti.mission, :obsid, example_gti.obsid)

    for instrument in instruments
        gtis = instrument_data[instrument].gtis
        gtis = gtis[(gtis[:, :stop]-gtis[:, :start]) .>= 16, :]
        available_orbits = unique(gtis[:orbit])

        instrument_orbit_plots = Dict{Int64,Plots.Plot}()
        for orbit in available_orbits
            instrument_data_orbit = Dict{Int64,JAXTAM.FFTData}()

            first_gti_idx = findfirst(gtis[:orbit] .== orbit)
            last_gti_idx  = findlast(gtis[:orbit] .== orbit)

            first_gti = gtis[first_gti_idx, :]
            last_gti  = gtis[last_gti_idx,  :]

            first_idx = findfirst(instrument_data[instrument].times .>= first_gti[:start])
            last_idx  = findfirst(instrument_data[instrument].times .>= last_gti[:stop])

            if last_gti[1, :stop] >= instrument_data[instrument].times[end]
                last_idx = length(instrument_data[instrument].times)
            end

            orbit_lc = BinnedData(
                instrument_data[instrument].mission,
                instrument_data[instrument].instrument,
                instrument_data[instrument].obsid,
                instrument_data[instrument].bin_time,
                instrument_data[instrument].counts[first_idx:last_idx],
                instrument_data[instrument].times[first_idx:last_idx],
                instrument_data[instrument].gtis[first_gti_idx:last_gti_idx, :]
            )

            if orbit == 0
                title_append = " - orbit unknown (likely end of lc)"
            else
                title_append = " - orbit $orbit/$(maximum(available_orbits))"
            end

            instrument_orbit_plots[orbit] = plot(Dict(instrument=>orbit_lc); save_plt=false,
                    title_append=title_append)

            if save_plt
                _savefig_obsdir(obs_row, orbit_lc.mission, orbit_lc.obsid, orbit_lc.bin_time, NaN, "$(instrument)_$(orbit)_lcurve.png")
            end
        end

        orbit_plots[instrument] = instrument_orbit_plots
    end

    return orbit_plots
end

# Power spectra plotting functions

function plot!(data::FFTData; title_append="", norm=:rms, rebin=(:log10, 0.01),
        lab="", logx=true, logy=true, show_errors=true,
        size_in=(1140,600), save_plt=true
    )
    bin_time_pow2 = Int(log2(data.bin_time))

    # Don't plot the 0Hz amplitude
    avg_amp = data.avg_amp
    freqs   = data.freqs
    avg_amp[1] = NaN
    #freqs[2] > 1 ? freq_min = prevpow(10, freqs[2]) : freq_min =  1/prevpow(10, 1/freqs[2]) # Skip zero freq
    freq_min = 0.01 # Manually set minimum fequency, assume anything lower isn't useful in power spectra
    freqs[1]   = NaN

    freqs, avg_amp, errors = fspec_rebin(data, rebin=rebin)
    ylab = ""
    
    if norm == :rms
        freqs, avg_amp, errors = fspec_rebin(data, rebin=rebin)
        errors = errors.*freqs
        avg_amp = (avg_amp.-2).*freqs
        amp_max = maximum(avg_amp[2:end]); amp_min = minimum(abs.(avg_amp[2:end]))
        avg_amp[avg_amp .<=0] .= NaN
        ylab = "Amplitude (Leahy - 2)*freq"
    elseif norm == :leahy
        amp_max = maximum(avg_amp[2:end]); amp_min = minimum(avg_amp[2:end])
        ylab="Amplitude (Leahy)"
    else
        @error "Plot norm type '$norm' not found" 
    end

    if show_errors
        Plots.plot!(freqs, avg_amp, color=:black,
            yerr=errors, #yerr=errors, marker=stroke(0.01, :black, :none),
            ylab=ylab, lab=lab)
    else
        Plots.plot!(freqs, avg_amp, color=:black,
            ylab=ylab, lab=lab)
    end

    if logx
        xaxis!(xscale=:log10, xformatter=xi->xi, xlim=(freq_min, freqs[end]))
    end

    if logy
        # If amp_min < 1, can't use prevpow10 for ylims, hacky little fix is 1/prevpow(10, 1/amp_min)
        # removed that anyway, set ylim to 1 if amp_min < 1
        amp_min > 1 ? ylim = (prevpow(10, amp_min), nextpow(10, amp_max)) : ylim = (1, nextpow(10, amp_max))
        yaxis!(yscale=:log10, yformatter=yi->yi, ylims=ylim)
    end
    
    Plots.plot!(xlab="Freq (Hz)", alpha=1,
        title="FFT - $(data.obsid) - 2^$(bin_time_pow2) bt - $(data.bin_size*data.bin_time) bs - $rebin rebin - $(data.bin_count) sections averaged $title_append")

    if(save_plt)
        _savefig_obsdir(data.mission, data.obsid, data.bin_time, data.bin_size, "fspec.png")
    end

    return Plots.plot!(size=size_in)
end

function plot(instrument_data::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}}; title_append="",
        size_in=(1140,600), norm=:rms, rebin=(:log10, 0.01), logx=true, logy=true, save_plt=true)
    instruments = keys(instrument_data)

    plt = Plots.plot()
    
    for instrument in instruments
        plt = JAXTAM.plot!(instrument_data[Symbol(instrument)][-1]; title_append=title_append,
            norm=norm, rebin=rebin, logx=logx, logy=logy, lab=String(instrument), save_plt=false)
    end

    Plots.plot!(title_location=:left, titlefontsize=12, margin=2mm, xguidefontsize=10, yguidefontsize=10)

    if(save_plt)
        example_gti = _pull_data(instrument_data)
        _savefig_obsdir(example_gti.mission, example_gti.obsid, example_gti.bin_time, example_gti.bin_size, "$fspec.png")
    end

    return Plots.plot!(size=size_in)
end

function plot_orbits(instrument_data::Dict{Symbol,Dict{Int64,JAXTAM.FFTData}};
    size_in=(1140,600), norm=:rms, rebin=(:log10, 0.01), logx=true, logy=true, save_plt=true)

    instruments = keys(instrument_data)

    orbit_plots = Dict{Symbol,Dict{Int64,Plots.Plot}}()

    example_gti = _pull_data(instrument_data)
    obs_row     = master_query(example_gti.mission, :obsid, example_gti.obsid)

    for instrument in instruments
        available_orbits = unique([gti.orbit for gti in values(instrument_data[instrument])])
        available_orbits = available_orbits[available_orbits .>= 0] # Excluse -1, -2, etc... for scrunched/mean FFTData

        instrument_orbit_plots = Dict{Int64,Plots.Plot}()
        for orbit in available_orbits
            instrument_data_orbit = Dict{Int64,JAXTAM.FFTData}()
            gtis_in_orbit = [gti.gti_index for gti in values(instrument_data[instrument]) if gti.orbit == orbit]

            for gti_no in gtis_in_orbit
                instrument_data_orbit[gti_no] = instrument_data[instrument][gti_no]
            end

            if orbit == 0
                title_append = " - orbit unknown (likely end of lc)"
            else
                title_append = " - orbit $orbit/$(maximum(available_orbits))"
            end

            # Have to run _scrunch_sections on orbit data, since it lacks the -1 indexed average amplitudes
            instrument_orbit_plots[orbit] = plot(Dict(instrument=>_scrunch_sections(instrument_data_orbit));
                    size_in=size_in, norm=norm, rebin=rebin, logx=logx, logy=logy, save_plt=false,
                    title_append=title_append)
            
            if save_plt
                data = instrument_data_orbit[gtis_in_orbit[1]]
                _savefig_obsdir(obs_row, data.mission, data.obsid, data.bin_time, data.bin_size, "$(instrument)_$(orbit)_fspec.png")
            end
        end

        orbit_plots[instrument] = instrument_orbit_plots
    end

    return orbit_plots
end

function plot_summary(mission_name::Symbol, obsid::String)
    lc_1s = lcurve(mission_name, obsid, 1)

    fs_2m13 = fspec(mission_name, obsid, 2.0^-13, 256)
    fs_low1 = fspec(mission_name, obsid, 2.0^-3, 512)
    fs_low2 = fspec(mission_name, obsid, 2.0^-1, 512)

    plt_lc = plot(lc_1s; size_in=(1000, 400))
end