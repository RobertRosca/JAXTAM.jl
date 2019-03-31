function auto_queue(mission::Mission)
    master_df = master_query_public(mission)

    queue = filter(x->x[:publicity]    == true,  master_df)
    queue = filter(x->x[:report_exists]== false, queue)
    queue = filter(x->x[:obs_type]     != "CAL", queue)
    queue = filter(x->x[:error]        == false, queue)

    return queue
end

"""
    auto_report(::Mission; limit::Union{Bool,Int}, update::Bool, nuke::Bool)

Calls `auto_queue` function to generate a queue of reports to make, the queue filters:

    * Public-only

    * Reportless

    * Not 'CAL' type observations

    * Error free

Leaving only suitable observations to be analysed

Calls 'report_all' using the queued observations

Will continue to generate reports until the `limit` is reached (if there is one)
"""
function auto_report(mission::Mission; limit=false, update=true, nuke=false)
    if update
        master(mission; update=true)
    end

    obs_queue = auto_queue(mission)

    i = 0
    if typeof(obs_queue) == DataFrames.DataFrameRow
        report_all(mission, obs_queue, update_masterpage=false)
    else
        for obs_row in DataFrames.eachrow(obs_queue)
            println(repeat("-", 8), "  ", obs_row[:obsid], "  ", repeat("-", 8), "\n")

            errors_dl = JAXTAM._log_query(mission, obs_row, "errors", :download; surpress_warn=true)

            if !ismissing(errors_dl)
                @warn "Download errors found, skipping"
                continue
            end
            
            report_all(mission, obs_row; nuke=nuke, update_masterpage=false)

            i+=1
            if i == limit
                break
            end
        end
    end

    master_append(mission; update=true)
    webgen_mission(mission)
end