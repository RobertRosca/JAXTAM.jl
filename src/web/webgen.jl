@tags html head body title meta div p h1 h2 h3 h4 h5 hr intro table thead tbody tr th td img a
@tags_noescape script style
@tags intro

function _webgen_offline_css(path_web)
    sheet_urls_css = [
        "cdn/css/",
        "https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.1.1/css/bootstrap.css",
        "https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap4.min.css",
        "https://cdn.datatables.net/v/bs4/dt-1.10.18/b-1.5.2/b-html5-1.5.2/fc-3.2.5/fh-3.1.4/sc-1.5.0/datatables.min.css"
    ]

    for sheet_urls in sheet_urls_css
        dir  = sheet_urls[1]
        urls = sheet_urls[2:end]
        file_names = basename.(urls)
        file_paths = string.(path_web, dir, file_names)
        println(file_paths)
    end
end

function _webgen_offline_s(path_web)
    sheet_urls_js  = [
        "cdn/js/",
        "https://code.jquery.com/jquery-3.3.1.js",
        "https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js",
        "https://cdn.datatables.net/1.10.19/js/dataTables.bootstrap4.min.js",
        "https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/pdfmake.min.js",
        "https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/vfs_fonts.js",
        "https://cdn.datatables.net/v/bs4/dt-1.10.18/b-1.5.2/b-html5-1.5.2/fc-3.2.5/fh-3.1.4/sc-1.5.0/datatables.min.js"
    ]

    for sheet_urls in sheet_urls_js
        dir  = sheet_urls[1]
        urls = sheet_urls[2:end]
        file_names = basename.(urls)
        file_paths = string.(path_web, dir, file_names)
        println(file_paths)
    end
end

function _webgen_head(;title_in="")
    node = m("head",
        m("title", title_in),
        m("meta"; charset="utf-8"),
        m("link"; rel="stylesheet", :type=>"text/css", href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.1.1/css/bootstrap.css"),
        m("link"; rel="stylesheet", :type=>"text/css", href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap4.min.css"),
        m("link"; rel="stylesheet", :type=>"text/css", href="https://cdn.datatables.net/v/bs4/dt-1.10.18/b-1.5.2/b-html5-1.5.2/fc-3.2.5/fh-3.1.4/sc-1.5.0/datatables.min.css"),
        m("style"; :type=>"text/css", class="init"),
        m("script"; :type=>"text/javascript", language="javascript", src="https://code.jquery.com/jquery-3.3.1.js"),
        m("script"; :type=>"text/javascript", language="javascript", src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"),
        m("script"; :type=>"text/javascript", language="javascript", src="https://cdn.datatables.net/1.10.19/js/dataTables.bootstrap4.min.js"),
        m("script"; :type=>"text/javascript", src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/pdfmake.min.js"),
        m("script"; :type=>"text/javascript", src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/vfs_fonts.js"),
        m("script"; :type=>"text/javascript", src="https://cdn.datatables.net/v/bs4/dt-1.10.18/b-1.5.2/b-html5-1.5.2/fc-3.2.5/fh-3.1.4/sc-1.5.0/datatables.min.js"),
        m("script"; :type=>"text/javascript", class="init"),
        script("
        \$(document).ready(function() {
            \$('#example').DataTable();
        } );"; :type=>"text/javascript", class="init")
    )
end

function _webgen_home_intro(mission_name::String)
    node_intro = intro(
        m("div"; class="se-pre-con"),
        div(class="container",
            m("div"; class="container"),
            m("h1", "JAXTAM.jl WebView - $mission_name"),
            m("hr"),
            m("p", "JAXTAM reports summary page for $mission_name")
        )
    )
end

function _add_obsid_url(df::Union{DataFrames.DataFrame, DataFrames.DataFrameRow}, path_web)
    obsid_url = Array{Union{Hyperscript.Node{Hyperscript.HTMLSVG},String},1}(undef, size(df, 1))

    report_path_rel = replace(df[:report_path], (path_web=>"."))

    obsid_url[  df[:report_exists]] = [a(df[i, :obsid], href=report_path_rel[i]) for i in findall(df[:report_exists])]
    obsid_url[.!df[:report_exists]] = df[.!df[:report_exists], :obsid]

    return obsid_url
end

function _webgen_table(df::DataFrames.DataFrame, path_web; table_id="example")
    # Replaces plain string obsid with hyperlink to the report
    if haskey(df, :obsid)
        obsid_url = _add_obsid_url(df, path_web)
        deletecols!(df, [:obsid, :report_path])
        df[:obsid] = obsid_url
        permutecols!(df, [:obsid; names(df)[1:end-1][:]])
    end

    if haskey(df, :countrate)
        df[:countrate][df[:countrate] .== Inf] .= -2.0
        df[:countrate][isnan.(df[:countrate])] .= -3.0
        df[:countrate] = round.(Int, df[:countrate])
    end

    if haskey(df, :time) && table_id=="report_page"
        # Include MJD time in report summary
        df[:time] = string.(df[:time], " | ", round.(JAXTAM._datetime2mjd.(convert.(Dates.DateTime, df[:time])), digits=6), " (MJD 6dp)")
    end
    
    rows, cols = size(df)
    headers = names(df)

    replace!(headers, :subject_category=>:cat)
    replace!(headers, :reports_exist=>:report)
    
    node_table = div(class="container",
        table(id=table_id, class="table table-striped table-bordered", style="width:100%", 
            thead(
                tr(
                    th.(headers)
                )
            ),
            tbody(
                tr.([td.([df[r, c] for c in 1:cols]) for r in 1:rows])
            )
        )
    )
end

function _webgen_table(df_row::DataFrames.DataFrameRow, path_web; table_id="example")
    df = DataFrame(Dict(zip(keys(df_row), values(df_row)))) # Convert DataFrameRow to DataFrame for table gen

    return _webgen_table(df, path_web; table_id=table_id)
end

function webgen_mission(mission::Mission; master_cache=false)
    c = 0

    master_df = master(mission, cache=master_cache)
    
    included_cols = [
        :name,
        :obsid,
        :subject_category,
        :obs_type,
        :publicity,
        :downloaded,
        :time,
        :report_path,
        :report_exists
    ]

    if haskey(master_df, :countrate)
        append!(included_cols, [:countrate])

        countrate = master_df[:countrate]
        countrate = floor.(Int, countrate)
    end

    html_out = html(
        _webgen_head(;title_in="JAXTAM $(_mission_name(mission)) homepage"),
        body(
            _webgen_home_intro(_mission_name(mission)),
            _webgen_table(master_df[:, included_cols], mission_paths(mission).web)
        )
    )

    path_web_mission = joinpath(mission_paths(mission).web, "index.html")

    write(path_web_mission, string(Pretty(html_out)))

    return path_web_mission
end