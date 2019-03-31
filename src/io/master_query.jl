
"""
    master_query(master_df::DataFrame, key_type::Symbol, key_value::Any)

Wrapper for a query, takes in an already loaded DataFrame `master_df`, a `key_type` to
search over (e.g. `obsid`), and a `key_value` to find (e.g. `0123456789`)

Returns the full row for any observations matching the search criteria
"""
function master_query(master_df::DataFrame, key_type::Symbol, key_value::Any)
    observations = filter(row -> row[key_type] == key_value, master_df)

    if size(observations, 1) == 0
        @warn "master_query returned no results for $key_type with $key_value search"
    end

    if size(observations, 1) == 1
        # Return DataFrameRow for single observation result
        observations = observations[1, :]
    end

    return observations
end

function master_query(mission::Mission, key_type::Symbol, key_value::Any)
    return master_query(master(mission), key_type, key_value)
end

function master_query_public(master_df::DataFrame, key_type::Symbol, key_value::Any)
    query_result = master_query(master_df, key_type, key_value)

    if typeof(query_result) == DataFrames.DataFrameRow
        # Single row returned by master query
        if !query_result[:publicity]
            query_result = DataFrame() # empty dataframe returned
        end
    else
        query_result = filter(x->x[:publicity], query_result)
    end

    return query_result    
end

function master_query_public(mission::Mission, key_type::Symbol, key_value::Any)
    return master_query_public(master(mission), key_type, key_value)
end

function master_query_public(mission::Mission)
    return filter(x->x[:publicity], master(mission))
end