struct PgramData <: JAXTAMData
    mission    :: Mission
    instrument :: Symbol
    obsid      :: String
    e_range    :: Tuple{Float64,Float64}
    bin_time   :: Real
    pg_type    :: Symbol
    power      :: Array
    freq       :: Array
    group      :: Int
end

function _pgram(counts, times, bin_time, pg_type=:standard; maximum_frequency=:auto)
    if maximum_frequency == :auto
        maximum_frequency = 0.5/bin_time
    elseif maximum_frequency > 0.5/bin_time
        @warn "Maximum frequency cannot be greater than 0.5/bin_time ($(0.5/bin_time))"
        maximum_frequency = 0.5/bin_time
    end

    pg_plan = LombScargle.plan(
        times, float(counts), 
        normalization=pg_type,
        minimum_frequency=1/(length(counts)*bin_time), # TODO: Check if setting min freq to 1/exposure is correct
        maximum_frequency=maximum_frequency,
        # samples_per_peak=1,
    )

    pg = lombscargle(pg_plan)

    freq, power = freqpower(pg)

    return freq, power
end

function _pgram(lc::BinnedData, pg_type, group; maximum_frequency=:auto)
    freq, power = (missing, missing)
    if group == 0 # Using whole lightcurve
        lc_groups = _group_return(lc; min_length_sec=16)

        if length(lc_groups) == 0
            throw(JAXTAMError("No groups found for `min_length_sec=16`", :pgram))
        end

        lc_group_counts = vcat([lc[2].counts for lc in lc_groups]...)
        lc_group_times  = vcat([lc[2].times for lc in lc_groups]...)
        freq, power = _pgram(lc_group_counts, lc_group_times, lc.bin_time, pg_type; maximum_frequency=maximum_frequency)
    else
        freq, power = _pgram(lc.counts, lc.times, lc.bin_time, pg_type; maximum_frequency=maximum_frequency)
    end
    
    return PgramData(lc.mission, lc.instrument, lc.obsid, lc.e_range, lc.bin_time,
        pg_type, power, freq, group)
end

function _pgram_lc_group_pad(group_lc::Dict{Int64,JAXTAM.BinnedData})
    groups = keys(group_lc)
    longest_group = maximum([length(x.counts) for x in values(group_lc)])

    padded_groups = Dict{Int64,JAXTAM.BinnedData}()

    for group in groups
        mission    = group_lc[group].mission
        instrument = group_lc[group].instrument
        obsid      = group_lc[group].obsid
        bin_time   = group_lc[group].bin_time
        counts     = group_lc[group].counts
        times      = group_lc[group].times
        gtis       = group_lc[group].gtis

        if length(counts) < longest_group
            counts = [counts; zeros(longest_group-length(counts))]
        end

        padded_groups[group] = BinnedData(mission, instrument, obsid, bin_time, counts, times, gtis)
    end

    return padded_groups
end

function pgram(instrument_lc::Dict{Symbol,BinnedData}; pg_type=:standard, per_group=false, maximum_frequency=:auto)
    instruments = keys(instrument_lc)

    instrument_pgram = per_group ? Dict{Symbol,Dict{Int,PgramData}}() : Dict{Symbol,PgramData}()
    for instrument in instruments
        if per_group
            group_pgram = Dict{Int,PgramData}()

            lc_groups = _group_return(instrument_lc[instrument])
            #lc_groups = _pgram_lc_group_pad(lc_groups)
            groups = keys(lc_groups)

            for group in groups
                group_pgram[group] = _pgram(lc_groups[group], pg_type, group, maximum_frequency=maximum_frequency)
            end

            #group_pgram[-1] = PgramData(group_pgram[1].mission, group_pgram[1].instrument, group_pgram[1].obsid,
            #    group_pgram[1].bin_time, pg_type, mean([pgram.power for pgram in values(group_pgram)]), group_pgram[1].freq, -1)

            instrument_pgram[instrument] = group_pgram
        else
            instrument_pgram[instrument] = _pgram(instrument_lc[instrument], pg_type, 0, maximum_frequency=maximum_frequency)
        end
    end
    
    return instrument_pgram
end